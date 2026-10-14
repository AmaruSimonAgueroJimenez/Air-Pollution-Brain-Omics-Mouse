#!/usr/bin/env python3
"""Collect original metadata and an exact, brain-only mouse download manifest."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
ROWS = []
SOURCES = []
GEO = {
    "GSE297423": {"project": "PRJNA1263837", "gms": None, "tissue": "anterior cingulate cortex", "pollutant": "PM2.5"},
    "GSE145840": {"project": "PRJNA608484", "gms": {f"GSM{x}" for x in range(4337093, 4337104)}, "tissue": "hypothalamus", "pollutant": "PM2.5"},
    "GSE93967": {"project": None, "gms": None, "tissue": "hippocampus", "pollutant": "PM2.5"},
}


def save_remote(url, relative, dataset):
    request = urllib.request.Request(url, headers={"User-Agent": "MouseBrainOmicsRepository/1.0"})
    with urllib.request.urlopen(request, timeout=90) as response:
        data = response.read()
    destination = ROOT / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(data)
    ROWS.append({"dataset": dataset, "kind": "metadata", "path": relative, "url": url,
                 "bytes": len(data), "md5": hashlib.md5(data).hexdigest()})
    return data


def add_data(dataset, kind, url, name=None, size=None, md5=None):
    if url.startswith("ftp://"):
        url = "https://" + url[6:]
    if not url.startswith("https://"):
        raise ValueError(f"Unexpected download scheme: {url}")
    filename = name or url.rsplit("/", 1)[-1]
    relative = f"data/{dataset}/{kind}/{filename}"
    if any(row["path"] == relative for row in ROWS):
        raise ValueError(f"Duplicate destination: {relative}")
    ROWS.append({"dataset": dataset, "kind": kind, "path": relative, "url": url,
                 "bytes": int(size) if size else None, "md5": md5 or None})


def parse_soft(text):
    series, samples, current = {}, {}, None
    for line in text.splitlines():
        if line.startswith("^SERIES = "):
            current = series
        elif line.startswith("^SAMPLE = "):
            accession = line.split(" = ", 1)[1]
            current = samples.setdefault(accession, {})
        elif line.startswith("^"):
            current = None
        elif current is not None and line.startswith("!") and " = " in line:
            key, value = line.split(" = ", 1)
            current.setdefault(key, []).append(value)
    return series, samples


def ena(dataset, project, allowed_samples=None):
    fields = ["run_accession", "experiment_accession", "sample_accession", "secondary_sample_accession",
              "sample_alias", "sample_title", "scientific_name", "library_strategy", "library_layout",
              "fastq_ftp", "fastq_md5", "fastq_bytes"]
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + project
           + "&result=read_run&fields=" + ",".join(fields) + "&format=tsv&download=false")
    blob = save_remote(url, f"metadata/{dataset}/ena_all_project_runs.tsv", dataset)
    records = list(csv.DictReader(io.StringIO(blob.decode()), delimiter="\t"))
    selected = [r for r in records if allowed_samples is None or
                set([r["sample_accession"], r["secondary_sample_accession"], r["sample_alias"]]) & allowed_samples]
    if not selected:
        raise ValueError(f"No scoped ENA records for {dataset}")
    if any(r["scientific_name"] not in ("Mus musculus", "Mus musculus domesticus") for r in selected):
        raise ValueError(f"Non-mouse FASTQ included in {dataset}")
    path = ROOT / f"metadata/{dataset}/ena_selected_brain_runs.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)
    for record in selected:
        addresses = record["fastq_ftp"].split(";")
        checksums = record["fastq_md5"].split(";")
        sizes = record["fastq_bytes"].split(";")
        if not addresses[0] or not (len(addresses) == len(checksums) == len(sizes)):
            raise ValueError(f"Incomplete ENA file metadata: {record['run_accession']}")
        for address, digest, size in zip(addresses, checksums, sizes):
            add_data(dataset, "fastq", "https://" + address, size=size, md5=digest)
    return selected


def geo_dataset(accession, config):
    bucket = accession[:-3] + "nnn"
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{accession}"
    compressed = save_remote(base + f"/soft/{accession}_family.soft.gz",
                             f"metadata/{accession}/{accession}_family.soft.gz", accession)
    text = gzip.decompress(compressed).decode()
    (ROOT / f"metadata/{accession}/{accession}_family.soft").write_text(text)
    series, samples = parse_soft(text)
    selected = {gsm: item for gsm, item in samples.items() if config["gms"] is None or gsm in config["gms"]}
    if config["gms"] is not None and set(selected) != config["gms"]:
        raise ValueError(f"Missing scoped samples for {accession}")
    if any(item.get("!Sample_organism_ch1") != ["Mus musculus"] for item in selected.values()):
        raise ValueError(f"Non-mouse GEO sample in {accession}")
    export = []
    allowed = set(selected)
    for gsm, item in selected.items():
        relations = item.get("!Sample_relation", [])
        for value in relations:
            allowed.update(re.findall(r"\bSAM[NED][A-Z]?\d+", value))
        export.append({"gsm": gsm, "title": item.get("!Sample_title", [""])[0],
                       "organism": item["!Sample_organism_ch1"][0],
                       "source": item.get("!Sample_source_name_ch1", [""])[0],
                       "description": " | ".join(item.get("!Sample_description", [])),
                       "characteristics": " | ".join(item.get("!Sample_characteristics_ch1", [])),
                       "relations": " | ".join(relations)})
        if accession == "GSE145840":
            supplementary = [url for key, urls in item.items() if key.startswith("!Sample_supplementary_file") for url in urls]
            if not supplementary:
                raise ValueError(f"Missing processed counts for {gsm}")
            for url in supplementary:
                if url.lower() != "none":
                    add_data(accession, "processed", url)
    with (ROOT / f"metadata/{accession}/samples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(export[0]))
        writer.writeheader()
        writer.writerows(export)
    if accession != "GSE145840":
        for url in series.get("!Series_supplementary_file", []):
            add_data(accession, "raw_array" if url.endswith("_RAW.tar") else "processed", url)
    if accession == "GSE93967":
        add_data(accession, "processed", base + f"/matrix/{accession}_series_matrix.txt.gz")
    runs = ena(accession, config["project"], allowed) if config["project"] else []
    SOURCES.append({"dataset": accession, "title": series["!Series_title"][0], "organism": "Mus musculus",
                    "tissue": config["tissue"], "pollutant": config["pollutant"],
                    "selected_samples": len(selected), "selected_runs": len(runs),
                    "source": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"})


def main():
    (ROOT / "metadata").mkdir(parents=True, exist_ok=True)
    for accession, config in GEO.items():
        geo_dataset(accession, config)
        print(f"Prepared {accession}", flush=True)
    dataset = "E-MTAB-14917"
    for filename in ["E-MTAB-14917.idf.txt", "E-MTAB-14917.sdrf.txt"]:
        save_remote(f"https://www.ebi.ac.uk/biostudies/files/{dataset}/{filename}",
                    f"metadata/{dataset}/{filename}", dataset)
    save_remote(f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{dataset}",
                f"metadata/{dataset}/biostudies.json", dataset)
    save_remote("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11942628/fullTextXML",
                f"metadata/{dataset}/article.xml", dataset)
    add_data(dataset, "processed", f"https://www.ebi.ac.uk/biostudies/files/{dataset}/E4_vs_E3_005_total.csv")
    add_data(dataset, "supplements", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11942628/supplementaryFiles",
             name="PMC11942628_supplementary.zip", size=1996686)
    runs = ena(dataset, "PRJEB86782")
    if len(runs) != 24:
        raise ValueError(f"Unexpected O3 library count: {len(runs)}")
    SOURCES.append({"dataset": dataset, "title": "Age- and ApoE genotype-dependent transcriptomic responses to O3 in the hippocampus of mice",
                    "organism": "Mus musculus", "tissue": "hippocampus", "pollutant": "O3",
                    "selected_samples": 24, "selected_runs": len(runs),
                    "source": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14917",
                    "processed_note": "The deposited CSV is a genotype differential-expression result, not a sample-level count matrix."})
    (ROOT / "metadata/datasets.json").write_text(json.dumps(SOURCES, indent=2, ensure_ascii=False) + "\n")
    (ROOT / "metadata/download_manifest.json").write_text(json.dumps(ROWS, indent=2) + "\n")
    with (ROOT / "metadata/download_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset", "kind", "path", "url", "bytes", "md5"])
        writer.writeheader()
        writer.writerows(ROWS)
    total = sum(r["bytes"] or 0 for r in ROWS if r["kind"] == "fastq")
    print(json.dumps({"files": len(ROWS), "fastq_files": sum(r["kind"] == "fastq" for r in ROWS),
                      "fastq_bytes": total, "fastq_GB": round(total / 1e9, 2)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
