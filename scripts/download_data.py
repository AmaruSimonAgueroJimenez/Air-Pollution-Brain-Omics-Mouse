#!/usr/bin/env python3
"""Download the manifest with curl resume, source MD5 checks, and atomic files."""
import argparse
import concurrent.futures
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
LOCK = threading.Lock()
RESULTS = {}


def digest(path, algorithm="md5"):
    result = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while block := handle.read(16 * 1024 * 1024):
            result.update(block)
    return result.hexdigest()


def emit(record):
    record["updated_at"] = datetime.now(timezone.utc).isoformat()
    with LOCK:
        RESULTS[record["path"]] = record
        state = {"updated_at": record["updated_at"], "files": list(RESULTS.values())}
        temporary = ROOT / "logs/download_status.json.tmp"
        temporary.write_text(json.dumps(state, indent=2) + "\n")
        os.replace(temporary, ROOT / "logs/download_status.json")
        with (ROOT / "logs/download_events.jsonl").open("a") as handle:
            handle.write(json.dumps(record) + "\n")
    if record["state"] in ("verified", "failed"):
        print(json.dumps(record), flush=True)


def validate(row, path):
    size = path.stat().st_size
    if row["bytes"] is not None and size != row["bytes"]:
        raise ValueError(f"Size mismatch: {size} != {row['bytes']}")
    if size == 0:
        raise ValueError("Empty file")
    actual_md5 = digest(path)
    if row["md5"] and actual_md5 != row["md5"]:
        raise ValueError(f"MD5 mismatch: {actual_md5} != {row['md5']}")
    with path.open("rb") as handle:
        head = handle.read(1024).lower().lstrip()
    if head.startswith((b"<!doctype html", b"<html")):
        raise ValueError("Received an HTML page instead of scientific data")
    if row["path"].endswith(".gz") and not head.startswith(b"\x1f\x8b"):
        raise ValueError("Invalid gzip file signature")
    return size, actual_md5


def download(row):
    destination = ROOT / row["path"]
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    record = {"dataset": row["dataset"], "kind": row["kind"], "path": row["path"], "expected_bytes": row["bytes"]}
    try:
        if destination.exists():
            try:
                size, actual_md5 = validate(row, destination)
            except ValueError:
                os.replace(destination, destination.with_name(destination.name + f".invalid.{time.time_ns()}"))
            else:
                emit({**record, "state": "verified", "bytes": size, "md5": actual_md5})
                return True
        for attempt in range(2):
            emit({**record, "state": "downloading"})
            errorlog = ROOT / "logs" / (destination.name + ".curl.log")
            command = ["curl", "--fail", "--location", "--continue-at", "-", "--retry", "8",
                       "--retry-all-errors", "--retry-delay", "5", "--connect-timeout", "30",
                       "--speed-time", "120", "--speed-limit", "1024", "--silent", "--show-error",
                       "--output", str(partial), row["url"]]
            with errorlog.open("a") as handle:
                completed = subprocess.run(command, stdout=handle, stderr=handle)
            if completed.returncode:
                raise RuntimeError(f"curl failed with code {completed.returncode}; see {errorlog.name}")
            emit({**record, "state": "checking"})
            try:
                size, actual_md5 = validate(row, partial)
            except ValueError as exc:
                os.replace(partial, partial.with_name(partial.name + f".invalid.{time.time_ns()}"))
                if attempt == 0:
                    emit({**record, "state": "retrying", "error": str(exc)})
                    continue
                raise
            os.replace(partial, destination)
            emit({**record, "state": "verified", "bytes": size, "md5": actual_md5})
            return True
    except Exception as exc:
        emit({**record, "state": "failed", "error": str(exc)})
        return False


def progress_loop(stop, rows):
    while not stop.wait(30):
        with LOCK:
            states = list(RESULTS.values())
        verified = [r for r in states if r["state"] == "verified"]
        completed_bytes = sum(r["bytes"] for r in verified)
        partial_bytes = 0
        for path in (ROOT / "data").rglob("*.part"):
            try:
                partial_bytes += path.stat().st_size
            except FileNotFoundError:
                pass  # A worker finished and atomically renamed this file.
        expected = sum(r["bytes"] or 0 for r in rows)
        print(json.dumps({"progress": {"verified_files": len(verified), "total_files": len(rows),
                                       "completed_GB": round(completed_bytes / 1e9, 3),
                                       "partial_GB": round(partial_bytes / 1e9, 3),
                                       "expected_GB": round(expected / 1e9, 3)}}), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--kind", action="append", help="Only selected kinds; default is all data")
    args = parser.parse_args()
    (ROOT / "logs").mkdir(exist_ok=True)
    process_lock = (ROOT / "logs/download.lock").open("a+")
    try:
        fcntl.flock(process_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("Another downloader already holds this repository's lock.")
    (ROOT / "logs/download.pid").write_text(str(os.getpid()) + "\n")
    running = {"complete": False, "running": True, "started_at": datetime.now(timezone.utc).isoformat()}
    completion_temporary = ROOT / "logs/completion.json.tmp"
    completion_temporary.write_text(json.dumps(running, indent=2) + "\n")
    os.replace(completion_temporary, ROOT / "logs/completion.json")
    manifest = json.loads((ROOT / "metadata/download_manifest.json").read_text())
    previous = ROOT / "logs/download_status.json"
    if previous.exists():
        RESULTS.update({row["path"]: row for row in json.loads(previous.read_text())["files"]})
    rows = manifest
    if args.kind:
        rows = [r for r in rows if r["kind"] in args.kind]
    stop = threading.Event()
    reporter = threading.Thread(target=progress_loop, args=(stop, manifest), daemon=True)
    reporter.start()
    started = time.time()
    try:
        # Processed data and metadata finish before sequencing downloads begin.
        small = [r for r in rows if r["kind"] != "fastq"]
        large = [r for r in rows if r["kind"] == "fastq"]
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            results.extend(pool.map(download, small))
            results.extend(pool.map(download, large))
        full_manifest_verified = all(RESULTS.get(row["path"], {}).get("state") == "verified"
                                     for row in manifest)
        report = {"finished_at": datetime.now(timezone.utc).isoformat(), "files": len(manifest),
                  "selected_files": len(rows), "selected_verified": sum(results),
                  "failed": len(results) - sum(results),
                  "elapsed_seconds": round(time.time() - started, 1), "complete": full_manifest_verified}
        with (ROOT / "metadata/download_receipts.csv").open("w", newline="") as handle:
            fields = ["dataset", "kind", "path", "state", "expected_bytes", "bytes", "md5", "updated_at", "error"]
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(RESULTS.values())
        if full_manifest_verified:
            audit = subprocess.run([sys.executable, str(ROOT / "scripts/verify_repository.py"), "--all"],
                                   cwd=ROOT)
            if audit.returncode:
                report["complete"] = False
                report["audit_failed"] = True
        report["running"] = False
        completion_temporary.write_text(json.dumps(report, indent=2) + "\n")
        os.replace(completion_temporary, ROOT / "logs/completion.json")
        print(json.dumps(report), flush=True)
        raise SystemExit(0 if all(results) and not report.get("audit_failed") else 1)
    finally:
        stop.set()
        process_lock.close()


if __name__ == "__main__":
    main()
