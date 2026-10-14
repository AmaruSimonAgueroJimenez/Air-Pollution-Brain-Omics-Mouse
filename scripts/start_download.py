#!/usr/bin/env python3
"""Run downloads independently of the chat session, logging all progress."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--workers", type=int, default=8)
args = parser.parse_args()
(root / "logs").mkdir(exist_ok=True)
pidfile = root / "logs/download.pid"
if pidfile.exists():
    old_pid = int(pidfile.read_text().strip())
    result = subprocess.run(["ps", "-p", str(old_pid), "-o", "args="], capture_output=True, text=True)
    if result.returncode == 0 and "download_data.py" in result.stdout:
        raise SystemExit(f"Downloader already running: PID {old_pid}")
with (root / "logs/download.log").open("a") as log:
    process = subprocess.Popen([sys.executable, "-u", str(root / "scripts/download_data.py"),
                                "--workers", str(args.workers)], cwd=root, stdin=subprocess.DEVNULL,
                               stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
print(f"Started downloader PID {process.pid}. Progress: python3 scripts/status.py", flush=True)
