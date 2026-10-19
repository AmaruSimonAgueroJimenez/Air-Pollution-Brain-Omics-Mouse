#!/usr/bin/env python3
"""Run the fixed four-study analysis and English Quarto report, in dependency order.

This entry point does not download the original study collection or install tools.
The ozone workflow resumes validated quantifications; other analysis stages rerun.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
STAGES = ("hypothalamus", "cortex", "mirna", "ozone", "enrichment", "spatial", "spatial-reference", "report")


def commands(cpus):
    python = sys.executable
    return {
        "hypothalamus": [python, "scripts/analyze_hypothalamus.py", "--cpus", str(cpus)],
        "cortex": [python, "scripts/analyze_cortex.py", "--workers", str(min(cpus, 6))],
        "mirna": [python, "scripts/analyze_mirna.py"],
        "ozone": [python, "scripts/ozone_pipeline.py", "--stage", "all", "--cpus", str(cpus)],
        "enrichment": [python, "scripts/analyze_enrichment.py", "--threads", str(min(cpus, 4))],
        "spatial": [python, "scripts/build_spatial_viewer.py"],
        "spatial-reference": [python, "scripts/analyze_spatial_reference.py"],
        "report": [python, "scripts/render_docs.py"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("all", *STAGES), default="all")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true", help="Print the execution plan without running it")
    args = parser.parse_args()
    if args.cpus < 1 or args.cpus > (os.cpu_count() or 1):
        parser.error("--cpus must be between 1 and the available logical CPU count")
    plan = commands(args.cpus)
    selected = STAGES if args.stage == "all" else (args.stage,)
    for stage in selected:
        command = plan[stage]
        print(f"[{stage}] {shlex.join(command)}", flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
