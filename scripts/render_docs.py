#!/usr/bin/env python3
"""Copy reviewed outputs and execute the Python Quarto website in docs/."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def prepare_assets():
    for source in (ROOT / "results").rglob("*"):
        if not source.is_file() or "quality_audit" in source.parts:
            continue
        relative = source.relative_to(ROOT / "results")
        # Read files and transient quantifier workspaces never belong in the report.
        if any(part in {"trimmed", "quant", "qc_subsamples"} for part in relative.parts):
            continue
        if any(part.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")) for part in relative.parts):
            continue
        if source.suffix in {".npz", ".pkl", ".h5ad"}:
            continue
        if "figures" in relative.parts and source.suffix in {".png", ".svg", ".pdf"}:
            parts = list(relative.parts)
            parts.remove("figures")
            destination = DOCS / "assets/figures" / Path(*parts)
        elif relative.parts[0] == "spatial_viewer" and source.suffix in {".png", ".svg", ".pdf"}:
            destination = DOCS / "assets/figures" / relative
        elif source.suffix in {".tsv", ".csv", ".json", ".gz", ".txt"}:
            destination = DOCS / "downloads" / relative
        else:
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    # Bundle source scripts for an inspectable, portable report copy.
    for source in (ROOT / "scripts").glob("*"):
        if source.is_file() and source.suffix in {".py", ".R"}:
            destination = DOCS / "downloads/code" / source.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    for source in (ROOT / "metadata").rglob("*"):
        if source.is_file() and source.suffix in {".json", ".csv", ".tsv", ".txt", ".soft", ".bib", ".ris"}:
            destination = DOCS / "downloads/provenance" / source.relative_to(ROOT / "metadata")
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    for source in (ROOT / "papers").glob("bibliography.*"):
        destination = DOCS / "downloads/papers" / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    packages = ["numpy", "pandas", "scipy", "scikit-learn", "statsmodels", "matplotlib", "seaborn",
                "h5py", "pydeseq2", "anndata", "gseapy", "nibabel", "scikit-image", "plotly", "ipykernel", "adjustText"]
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    environment = {"python": sys.version, "executable": str(Path(sys.executable).relative_to(ROOT))
                   if Path(sys.executable).is_relative_to(ROOT) else sys.executable,
                   "platform": platform.platform(), "packages": versions}
    destination = ROOT / "metadata/analysis_environment.json"
    destination.write_text(json.dumps(environment, indent=2) + "\n")
    shutil.copy2(destination, DOCS / "downloads/analysis_environment.json")
    requirements = "# Versions used for the Python analyses and Quarto render.\n" + "\n".join(
        f"{name}=={version}" for name, version in versions.items() if version is not None) + "\n"
    (ROOT / "requirements-analysis.txt").write_text(requirements)
    shutil.copy2(ROOT / "requirements-analysis.txt", DOCS / "downloads/requirements-analysis.txt")
    sources = list(DOCS.glob("*.qmd")) + [DOCS / name for name in (
        "_quarto.yml", "_report_utils.py", "references.bib", "assets/report.css",
        "assets/viewer/index.html", "assets/viewer/viewer.css", "assets/viewer/viewer.js")]
    with zipfile.ZipFile(DOCS / "downloads/report-source.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for source in sources:
            if source.exists():
                archive.write(source, source.relative_to(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--page")
    args = parser.parse_args()
    prepare_assets()
    if args.prepare_only:
        return
    prefix = ROOT / ".venv-analysis"
    subprocess.run([sys.executable, "-m", "ipykernel", "install", "--prefix", str(prefix),
                    "--name", "omics-analysis", "--display-name", "Omics Analysis (Python)"], check=True)
    environment = os.environ.copy()
    environment["QUARTO_PYTHON"] = sys.executable
    environment["JUPYTER_PATH"] = str(prefix / "share/jupyter") + os.pathsep + environment.get("JUPYTER_PATH", "")
    command = ["quarto", "render"] + ([args.page] if args.page else [])
    subprocess.run(command, cwd=DOCS, env=environment, check=True)


if __name__ == "__main__":
    main()
