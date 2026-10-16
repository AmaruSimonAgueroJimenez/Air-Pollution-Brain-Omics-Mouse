#!/usr/bin/env python3
"""Fetch the four Zhuang MERFISH matrices and their public Allen reference data."""
import argparse
from concurrent.futures import ThreadPoolExecutor
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
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "metadata/spatial_atlas"
LOGS = ROOT / "logs/spatial_atlas"
DATA = ROOT / "data/reference_atlases/Zhuang_2023"
MANIFEST = METADATA / "download_manifest.json"
STATE = LOGS / "status.json"
RELEASE = "20261031"
SOURCE = f"https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/releases/{RELEASE}/manifest.json"
MUTEX = threading.Lock()
RESULTS = {}


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(obj, indent=2) + "\n")
    os.replace(temporary, path)


def leaves(value):
    if isinstance(value, dict):
        if "relative_path" in value and "size" in value:
            yield value
        else:
            for child in value.values():
                yield from leaves(child)


def prepare():
    content = urlopen(SOURCE, timeout=60).read()
    source = json.loads(content)
    METADATA.mkdir(parents=True, exist_ok=True)
    (METADATA / f"allen_release_{RELEASE}.json").write_bytes(content)
    selected = []
    for animal in range(1, 5):
        for suffix in ("", "-CCF"):
            selected.extend(leaves(source["file_listing"][f"Zhuang-ABCA-{animal}{suffix}"]))
    selected.extend(leaves(source["file_listing"]["WMB-taxonomy"]["metadata"]))
    selected.extend(leaves(source["file_listing"]["Allen-CCF-2020"]))
    rows = []
    for item in selected:
        relative = Path(item["relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe source path: {relative}")
        kind = relative.parts[0]
        if kind == "expression_matrices":
            kind = "measured_counts" if "-raw.h5ad" in relative.name else "normalized_expression"
        rows.append({"dataset": relative.parts[1], "kind": kind,
                     "path": str((DATA / relative).relative_to(ROOT)),
                     "url": item["url"], "bytes": item["size"],
                     "md5": item["file_hash"], "source_version": item["version"]})
    if len({r["path"] for r in rows}) != len(rows):
        raise ValueError("Duplicate source paths")
    save(MANIFEST, {"retrieved_at": now(), "source_manifest": SOURCE,
                    "release_label": RELEASE, "doi": "10.1038/s41586-023-06808-9",
                    "experimental_pollutant_exposure": False,
                    "scope": "Four measured MERFISH datasets (raw counts and log2 matrices), cell/gene metadata, CCF coordinates, taxonomy and anatomical reference volumes. Microscopy image stacks and transcriptome-wide imputation are not included.",
                    "total_bytes": sum(r["bytes"] for r in rows), "files": rows})
    print(json.dumps({"manifest_files": len(rows), "GB": sum(r["bytes"] for r in rows) / 1e9}))


def emit(row, state, **extra):
    with MUTEX:
        RESULTS[row["path"]] = {"path": row["path"], "dataset": row["dataset"],
                                "state": state, "updated_at": now(), **extra}
        save(STATE, {"updated_at": now(), "files": list(RESULTS.values())})
    if state in ("verified", "failed"):
        print(json.dumps(RESULTS[row["path"]]), flush=True)


def validate(row, path):
    if path.stat().st_size != row["bytes"]:
        raise ValueError(f"Unexpected size for {path.name}")
    digest = hashlib.md5()
    with path.open("rb") as handle:
        head = handle.read(8)
        digest.update(head)
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    if digest.hexdigest() != row["md5"]:
        raise ValueError(f"Source MD5 mismatch for {path.name}")
    if row["path"].endswith(".h5ad") and head != b"\x89HDF\r\n\x1a\n":
        raise ValueError("Invalid HDF5 signature")
    if row["path"].endswith(".gz") and head[:2] != b"\x1f\x8b":
        raise ValueError("Invalid gzip signature")
    return {"bytes": row["bytes"], "md5": digest.hexdigest(),
            "source_md5_match": True, "mtime_ns": path.stat().st_mtime_ns}


def transfer(row, verify_only=False):
    destination = ROOT / row["path"]
    if not destination.resolve().is_relative_to(DATA.resolve()):
        raise ValueError("Destination outside the reference atlas directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    try:
        if destination.exists():
            try:
                receipt = validate(row, destination)
            except ValueError:
                if verify_only:
                    raise
                os.replace(destination, destination.with_name(destination.name + f".invalid.{time.time_ns()}"))
            else:
                emit(row, "verified", **receipt)
                return True
        if verify_only:
            raise FileNotFoundError(str(destination))
        for attempt in range(2):
            if partial.exists() and partial.stat().st_size >= row["bytes"]:
                try:
                    receipt = validate(row, partial)
                except ValueError:
                    os.replace(partial, partial.with_name(partial.name + f".invalid.{time.time_ns()}"))
                else:
                    os.replace(partial, destination)
                    emit(row, "verified", **receipt)
                    return True
            emit(row, "downloading")
            command = ["curl", "--fail", "--location", "--continue-at", "-", "--retry", "8",
                       "--retry-all-errors", "--retry-delay", "5", "--connect-timeout", "30",
                       "--speed-time", "120", "--speed-limit", "1024", "--silent", "--show-error",
                       "--output", str(partial), row["url"]]
            logname = row["dataset"] + "_" + destination.name + ".curl.log"
            with (LOGS / logname).open("a") as log:
                subprocess.run(command, stdout=log, stderr=log, check=True)
            emit(row, "checking")
            try:
                receipt = validate(row, partial)
            except ValueError as exc:
                os.replace(partial, partial.with_name(partial.name + f".invalid.{time.time_ns()}"))
                if attempt == 0:
                    emit(row, "retrying", error=str(exc))
                    continue
                raise
            os.replace(partial, destination)
            emit(row, "verified", **receipt)
            return True
    except Exception as exc:
        emit(row, "failed", error=str(exc))
        return False


def status():
    manifest = json.loads(MANIFEST.read_text())
    states = json.loads(STATE.read_text())["files"] if STATE.exists() else []
    known = {r["path"]: r for r in states}
    complete = 0
    present = 0
    failed = 0
    for row in manifest["files"]:
        destination = ROOT / row["path"]
        partial = destination.with_name(destination.name + ".part")
        record = known.get(row["path"], {})
        if destination.exists():
            stat = destination.stat()
            present += stat.st_size
            complete += (record.get("state") == "verified" and stat.st_size == row["bytes"]
                         and stat.st_mtime_ns == record.get("mtime_ns"))
        elif partial.exists():
            present += partial.stat().st_size
        failed += record.get("state") == "failed"
    active = False
    pidfile = LOGS / "download.pid"
    if pidfile.exists():
        result = subprocess.run(["ps", "-p", pidfile.read_text().strip(), "-o", "args="],
                                capture_output=True, text=True)
        active = result.returncode == 0 and "download_spatial_atlas.py" in result.stdout
    print(json.dumps({"verified_files": complete, "total_files": len(manifest["files"]),
                      "present_GB": round(present / 1e9, 3),
                      "expected_GB": round(manifest["total_bytes"] / 1e9, 3),
                      "failed": failed, "process_running": active,
                      "complete": complete == len(manifest["files"])}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--background", action="store_true")
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.prepare:
        prepare()
        return
    if args.status:
        status()
        return
    LOGS.mkdir(parents=True, exist_ok=True)
    if args.background:
        with (LOGS / "download.log").open("a") as log:
            process = subprocess.Popen([sys.executable, "-u", str(Path(__file__).resolve()),
                                        "--workers", str(args.workers)], cwd=ROOT,
                                       stdin=subprocess.DEVNULL, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        print(json.dumps({"started_pid": process.pid, "log": str(LOGS / "download.log")}))
        return
    with (LOGS / "download.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("A spatial atlas downloader is already running")
        (LOGS / "download.pid").write_text(str(os.getpid()) + "\n")
        manifest = json.loads(MANIFEST.read_text())
        if STATE.exists():
            RESULTS.update({r["path"]: r for r in json.loads(STATE.read_text())["files"]})
        save(LOGS / "completion.json", {"complete": False, "running": True, "started_at": now()})
        # Finish the smaller reference files before the largest expression matrices.
        rows = sorted(manifest["files"], key=lambda r: r["bytes"])
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            outcomes = list(pool.map(lambda row: transfer(row, args.verify_only), rows))
        report = {"complete": all(outcomes), "running": False, "finished_at": now(),
                  "verified_files": sum(outcomes), "total_files": len(rows),
                  "failed": len(rows) - sum(outcomes), "total_bytes": manifest["total_bytes"]}
        save(LOGS / "completion.json", report)
        save(METADATA / "download_receipts.json", {"source_manifest": SOURCE,
                                                   "files": list(RESULTS.values()), **report})
        print(json.dumps(report), flush=True)
        if not all(outcomes):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
