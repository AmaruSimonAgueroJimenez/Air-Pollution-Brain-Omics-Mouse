#!/usr/bin/env python3
"""Verify brain sample selection, matrices, mapping, and download receipts."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import tarfile

from prepare_manifest import parse_soft

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true", help="Require every original FASTQ to be verified")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "metadata/download_manifest.json").read_text())
    if len({r["path"] for r in manifest}) != len(manifest):
        raise ValueError("Duplicate manifest paths")
    results = {}
    with gzip.open(ROOT / "data/GSE297423/processed/GSE297423_mouse_gene_count.txt.gz", "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames[1:21]
        row_count = 0
        for row in reader:
            for column in columns:
                value = int(row[column])
                if value < 0:
                    raise ValueError("Negative RNA-seq count")
            row_count += 1
    _, samples = parse_soft((ROOT / "metadata/GSE297423/GSE297423_family.soft").read_text())
    mapping = []
    for gsm, sample in samples.items():
        descriptions = sample.get("!Sample_description", [])
        names = [text.split(":", 1)[1].strip() for text in descriptions if text.startswith("Library name:")]
        if len(names) != 1:
            raise ValueError(f"Cannot identify library name for {gsm}: {descriptions}")
        original = names[0]
        resolved = original
        inferred = False
        if original not in columns:
            if gsm == "GSM8990962" and original == "PMPFFPM2778" and "PMPFFPMA2778" in columns:
                resolved = "PMPFFPMA2778"
                inferred = True
            else:
                raise ValueError(f"Unknown matrix column: {gsm} {original}")
        mapping.append({"gsm": gsm, "title": sample["!Sample_title"][0], "source_library_name": original,
                        "matrix_column": resolved, "mapping_status": "inferred_typo_needs_confirmation" if inferred else "explicit_source_match"})
    if len(mapping) != 20 or len({r["matrix_column"] for r in mapping}) != 20:
        raise ValueError("Incomplete or duplicate GSE297423 matrix mapping")
    with (ROOT / "metadata/GSE297423/library_mapping.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mapping[0]))
        writer.writeheader()
        writer.writerows(mapping)
    results["GSE297423"] = {"genes": row_count, "sample_columns": 20, "mapping_inferred_aliases": 1}

    hypo = list((ROOT / "data/GSE145840/processed").glob("*.txt.gz"))
    if len(hypo) != 11:
        raise ValueError(f"Expected 11 hypothalamus matrices, found {len(hypo)}")
    gene_sets = []
    for path in sorted(hypo):
        genes = set()
        with gzip.open(path, "rt") as handle:
            for line in handle:
                gene, count = line.split()
                if gene in genes or int(count) < 0:
                    raise ValueError(f"Invalid gene/count entry in {path.name}")
                genes.add(gene)
        gene_sets.append(genes)
    common_genes = set.intersection(*gene_sets)
    all_genes = set.union(*gene_sets)
    selected = list(csv.DictReader((ROOT / "metadata/GSE145840/ena_selected_brain_runs.tsv").open(), delimiter="\t"))
    if len(selected) != 11 or {r["sample_alias"] for r in selected} != {f"GSM{x}" for x in range(4337093, 4337104)}:
        raise ValueError("ENA hypothalamus scope mismatch")
    results["GSE145840"] = {"gene_rows_per_file": {path.name: len(genes) for path, genes in zip(sorted(hypo), gene_sets)},
                            "common_genes": len(common_genes), "union_genes": len(all_genes),
                            "gene_sets_identical": all(genes == gene_sets[0] for genes in gene_sets),
                            "hypothalamus_samples": 11, "fastq_runs": 11,
                            "note": "Source files have different gene lists. Preserve missing entries; do not assume missing genes have zero counts."}

    with gzip.open(ROOT / "data/GSE93967/processed/GSE93967_series_matrix.txt.gz", "rt") as handle:
        text = handle.read()
    table = text.split("!series_matrix_table_begin\n", 1)[1].split("!series_matrix_table_end", 1)[0]
    records = list(csv.reader(table.splitlines(), delimiter="\t"))
    if len(records[0]) != 10 or len(records) - 1 != 1268:
        raise ValueError("Unexpected miRNA matrix shape")
    with tarfile.open(ROOT / "data/GSE93967/raw_array/GSE93967_RAW.tar") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        if len(members) != 9:
            raise ValueError("Expected nine source microarray files")
    results["GSE93967"] = {"probe_rows": 1268, "samples": 9, "raw_array_files": len(members)}

    with (ROOT / "data/E-MTAB-14917/processed/E4_vs_E3_005_total.csv").open() as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "diffexpressed"]:
            raise ValueError("Unexpected ozone result table columns")
        genes = sum(1 for _ in reader)
    ozone_runs = list(csv.DictReader((ROOT / "metadata/E-MTAB-14917/ena_selected_brain_runs.tsv").open(), delimiter="\t"))
    if len(ozone_runs) != 24:
        raise ValueError("Expected 24 hippocampal O3 libraries")
    results["E-MTAB-14917"] = {"genes_in_genotype_result_table": genes, "runs": 24,
                               "processed_table_is_counts": False}

    status = json.loads((ROOT / "logs/download_status.json").read_text())["files"]
    indexed = {r["path"]: r for r in status}
    completed = 0
    for row in manifest:
        record = indexed.get(row["path"], {})
        if record.get("state") == "verified":
            path = ROOT / row["path"]
            if not path.exists() or path.stat().st_size != record["bytes"]:
                raise ValueError(f"Verified file size changed: {row['path']}")
            if row["md5"] and record["md5"] != row["md5"]:
                raise ValueError(f"Source MD5 receipt mismatch: {row['path']}")
            completed += 1
        elif args.all or row["kind"] != "fastq":
            raise ValueError(f"Required file not verified: {row['path']}")
    report = {"processed_data_valid": True, "datasets": results, "verified_files": completed,
              "manifest_files": len(manifest), "all_downloads_complete": completed == len(manifest)}
    (ROOT / "metadata/validation_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
