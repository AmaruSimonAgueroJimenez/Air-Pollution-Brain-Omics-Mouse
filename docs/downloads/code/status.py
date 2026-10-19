#!/usr/bin/env python3
"""Summarize verified files and unfinished byte ranges, without rescanning FASTQ."""
import json
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / "metadata/download_manifest.json").read_text())
status_path = root / "logs/download_status.json"
statuses = json.loads(status_path.read_text())["files"] if status_path.exists() else []
known = {row["path"]: row for row in statuses}
summary = {}
for row in manifest:
    group = summary.setdefault(row["dataset"], {"files": 0, "verified_files": 0, "present_GB": 0,
                                                "expected_GB": 0, "failed": 0})
    group["files"] += 1
    group["expected_GB"] += (row["bytes"] or 0) / 1e9
    destination = root / row["path"]
    partial = destination.with_name(destination.name + ".part")
    if destination.exists():
        group["present_GB"] += destination.stat().st_size / 1e9
    elif partial.exists():
        group["present_GB"] += partial.stat().st_size / 1e9
    result = known.get(row["path"], {})
    group["verified_files"] += result.get("state") == "verified"
    group["failed"] += result.get("state") == "failed"
for group in summary.values():
    for key in ("present_GB", "expected_GB"):
        group[key] = round(group[key], 3)
complete_file = root / "logs/completion.json"
process_running = False
pidfile = root / "logs/download.pid"
if pidfile.exists():
    process = subprocess.run(["ps", "-p", pidfile.read_text().strip(), "-o", "args="], capture_output=True, text=True)
    process_running = process.returncode == 0 and "download_data.py" in process.stdout
print(json.dumps({"datasets": summary, "process_running": process_running,
                  "completion": json.loads(complete_file.read_text()) if complete_file.exists() else None},
                 indent=2))
