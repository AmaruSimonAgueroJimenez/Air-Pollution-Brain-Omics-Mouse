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
import hashlib

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def promote_site():
    """Place the compiled static website in /docs for branch-based GitHub Pages."""
    built = DOCS / "_site"
    if not (built / "index.html").exists() or not (built / "es/index.html").exists():
        raise RuntimeError("Render both language versions before preparing GitHub Pages")
    published = []
    for source in sorted(built.rglob("*")):
        if not source.is_file():
            continue
        relative = source.relative_to(built)
        if relative.as_posix() == "pages-manifest.json":
            continue  # A manifest cannot include its own content hash.
        if any(part in {".quarto", "_freeze", "_site", "private_inputs"} for part in relative.parts):
            continue
        destination = DOCS / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        published.append({"path": relative.as_posix(), "bytes": source.stat().st_size,
                          "sha256": hashlib.sha256(source.read_bytes()).hexdigest()})
    (DOCS / ".nojekyll").write_text("")
    manifest = {"publishing_source": "main:/docs", "entry_point": "index.html",
                "spanish_entry_point": "es/index.html", "files": published,
                "total_bytes": sum(item["bytes"] for item in published),
                "excludes": ["original FASTQ/H5AD", "reference genomes", "fitted R objects",
                             "article PDFs/full text", "private_inputs", "local environments"]}
    manifest["manifest_self_excluded"] = True
    manifest_text = json.dumps(manifest, indent=2) + "\n"
    (DOCS / "pages-manifest.json").write_text(manifest_text)
    (built / "pages-manifest.json").write_text(manifest_text)
    print(f"GitHub Pages ready in docs/: {len(published)} files, {manifest['total_bytes']/1e6:.1f} MB")


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
            # Keep unrelated in-progress downloader maintenance out of this release.
            if source.name in {"download_data.py", "start_download.py"}:
                committed = subprocess.run(["git", "show", "HEAD:scripts/" + source.name],
                                           cwd=ROOT, check=True, capture_output=True)
                destination.write_bytes(committed.stdout)
            else:
                shutil.copy2(source, destination)
    for source in (ROOT / "metadata").rglob("*"):
        if "private_inputs" in source.relative_to(ROOT / "metadata").parts:
            continue
        if source.is_file() and source.suffix in {".json", ".csv", ".tsv", ".txt", ".soft", ".bib", ".ris"}:
            destination = DOCS / "downloads/provenance" / source.relative_to(ROOT / "metadata")
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    for source in (ROOT / "papers").glob("bibliography.*"):
        destination = DOCS / "downloads/papers" / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    # Preserve source QC locally, but publish portable paths instead of account paths.
    for destination in (DOCS / "downloads").rglob("*"):
        if destination.is_file() and destination.suffix in {".json", ".csv", ".tsv", ".txt", ".soft"}:
            original = destination.read_text()
            portable = original.replace(str(ROOT), ".").replace(str(Path.home()), "~")
            if portable != original:
                destination.write_text(portable)
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
    sources = list(DOCS.glob("*.qmd")) + list((DOCS / "es").glob("*.qmd")) + [DOCS / name for name in (
        "_quarto.yml", "_report_utils.py", "references.bib", "assets/report.css",
        "assets/viewer/index.html", "assets/viewer/index-es.html", "assets/viewer/viewer.css", "assets/viewer/viewer.js",
        "assets/viewer/annotations.js", "assets/viewer-cortex/index.html", "assets/viewer-cortex/index-es.html", "assets/viewer-cortex/viewer.css",
        "assets/spatial-explorer/index.html", "assets/spatial-explorer/index-es.html",
        "assets/spatial-explorer/explorer.js", "assets/spatial-explorer/explorer.css",
        "assets/defense-slides/defense.css", "assets/defense-slides/defense.js",
        "es/_quarto.yml", "es/_report_utils.py")]
    with zipfile.ZipFile(DOCS / "downloads/report-source.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for source in sources:
            if source.exists():
                archive.write(source, source.relative_to(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--refresh-static", action="store_true",
                        help="Refresh static assets/downloads and promote an existing bilingual render")
    parser.add_argument("--page")
    parser.add_argument("--promote-only", action="store_true",
                        help="Copy an already rendered bilingual site into /docs for GitHub Pages")
    parser.add_argument("--language", choices=("both", "en", "es"), default="both",
                        help="Render both language versions by default, or select one")
    args = parser.parse_args()
    if args.promote_only:
        promote_site()
        return
    prepare_assets()
    if args.prepare_only:
        return
    if args.refresh_static:
        for folder in ("assets", "downloads"):
            shutil.copytree(DOCS / folder, DOCS / "_site" / folder, dirs_exist_ok=True)
        promote_site()
        return
    prefix = ROOT / ".venv-analysis"
    subprocess.run([sys.executable, "-m", "ipykernel", "install", "--prefix", str(prefix),
                    "--name", "omics-analysis", "--display-name", "Omics Analysis (Python)"], check=True)
    environment = os.environ.copy()
    environment["QUARTO_PYTHON"] = sys.executable
    environment["JUPYTER_PATH"] = str(prefix / "share/jupyter") + os.pathsep + environment.get("JUPYTER_PATH", "")
    command = ["quarto", "render"] + ([args.page] if args.page else [])
    languages = ("en", "es") if args.language == "both" else (args.language,)
    for language in languages:
        source_dir = DOCS if language == "en" else DOCS / "es"
        subprocess.run(command, cwd=source_dir, env=environment, check=True)
    promote_site()


if __name__ == "__main__":
    main()
