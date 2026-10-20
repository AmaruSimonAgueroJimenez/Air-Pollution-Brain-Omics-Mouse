#!/usr/bin/env python3
"""Publish frozen gene products/functions and regional priorities to both viewers."""
from pathlib import Path
import hashlib
import json
import math

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main():
    annotation_path = ROOT / "results/gene_priorities/annotations.tsv"
    priority_path = ROOT / "results/gene_priorities/prioritized_genes.tsv"
    annotations = pd.read_csv(annotation_path, sep="\t").fillna("")
    priorities = pd.read_csv(priority_path, sep="\t")
    if annotations.gene.duplicated().any():
        raise ValueError("Gene annotations must have unique source symbols")
    expected = {"hypothalamus_primary_FDR": 11, "hypothalamus_measured_top10": 10,
                "cortex_primary_top20": 20, "cortex_measured_top10": 10}
    for name, n in expected.items():
        if int(priorities.selection_group.eq(name).sum()) != n:
            raise ValueError("Unexpected selected-gene population: " + name)
    mapping = {row["gene"]: row for row in annotations.to_dict("records")}
    # Preserve assay symbols while allowing a verified current symbol in the UI.
    for row in annotations.to_dict("records"):
        mapping.setdefault(row["official_symbol"], row)
    regional = {"HY": [], "ACA": []}
    for region_key, groups in (("HY", ("hypothalamus_primary_FDR", "hypothalamus_measured_top10")),
                               ("ACA", ("cortex_primary_top20", "cortex_measured_top10"))):
        seen = {}
        for group in groups:
            part = priorities.loc[priorities.selection_group.eq(group)].sort_values("selection_rank")
            for row in part.to_dict("records"):
                gene = row["gene"]
                if gene in seen:
                    seen[gene]["selection_groups"].append(group)
                    continue
                entry = {**row, "symbol": gene, "measured": bool(row["atlas_measured"]),
                         "significant": bool(row["primary_significant"]), "q": row["padj"],
                         "selection_groups": [group]}
                seen[gene] = entry
                regional[region_key].append(entry)
    for gene in annotations.to_dict("records"):
        if not str(gene["source_url"]).startswith("https://"):
            raise ValueError("A curated gene requires its official source URL")
    target = ROOT / "docs/assets/viewer/annotations.js"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("window.GENE_ANNOTATIONS=" + json.dumps(clean(mapping), ensure_ascii=False, allow_nan=False, separators=(",", ":")) + ";\n"
                      + "window.GENE_PRIORITIES=" + json.dumps(clean(regional), ensure_ascii=False, allow_nan=False, separators=(",", ":")) + ";\n")
    audit = {"annotation_rows": len(annotations), "viewer_entries": {k: len(v) for k, v in regional.items()},
             "measured_shortlist_genes": {k: sum(r["measured"] for r in v) for k, v in regional.items()},
             "source_annotation_sha256": hashlib.sha256(annotation_path.read_bytes()).hexdigest(),
             "source_priority_sha256": hashlib.sha256(priority_path.read_bytes()).hexdigest(),
             "viewer_annotation_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
             "no_expression_imputation": True, "source_symbols_preserved": True, "passed": True}
    (ROOT / "results/report_validation/gene_annotation_assets.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
