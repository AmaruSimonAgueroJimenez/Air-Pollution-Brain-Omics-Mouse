#!/usr/bin/env python3
"""Mouse-specific GO Biological Process ORA and exploratory preranked GSEA.

Run .venv-analysis/bin/python scripts/analyze_enrichment.py.
Add a later result with --dataset NAME=path/to/DE.tsv --only NAME.
References are downloaded once from GO and thereafter checked against their saved
SHA256 values. No human orthology or Enrichr symbol conversion is used.
"""
from __future__ import annotations

import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import lru_cache
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import textwrap

import gseapy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/enrichment"
REF = ROOT / "data/references/go"
META = ROOT / "metadata/functional_annotation"
SOURCE_URLS = {"mgi.gaf.gz": "https://current.geneontology.org/annotations/mgi.gaf.gz",
               "go-basic.obo": "https://current.geneontology.org/ontology/go-basic.obo"}
DEFAULT_DATASETS = {"hypothalamus": ROOT / "results/hypothalamus/differential_expression.tsv",
                    "cortex": ROOT / "results/cortex/DE_PM_vs_PBS.tsv"}
SEED = 20261018
FDR = .05
MIN_SIZE = 15
MAX_SIZE = 500
PERMUTATIONS = 1000
COLORS = {"up": "#D86F12", "down": "#008579"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_clean(value):
    if isinstance(value, dict): return {str(k): json_clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [json_clean(v) for v in value]
    if isinstance(value, np.generic): return json_clean(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    if value is pd.NA: return None
    return value


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(json_clean(data), indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def save_table(frame, path, index=False):
    frame.to_csv(path, sep="\t", index=index, na_rep="NA", float_format="%.10g")


def summary_count_table(summaries):
    """Derive the report count table from every accumulated dataset summary."""
    rows = []
    for name, summary in summaries.items():
        ora, gsea = summary["ORA"], summary["GSEA"]
        rows.append({
            "dataset": name,
            "ORA_universe": ora["universe_size"],
            "ORA_up_genes": ora["up_foreground_size"],
            "ORA_down_genes": ora["down_foreground_size"],
            "ORA_up_significant_terms": ora["significant_up_terms_BH_FDR_0_05"],
            "ORA_down_significant_terms": ora["significant_down_terms_BH_FDR_0_05"],
            "GSEA_rank_genes": gsea["ranking_size"],
            "GSEA_terms_tested": gsea["tested_terms"],
            "GSEA_terms_algorithmic_FDR_0_05": gsea["significant_terms_algorithmic_FDR_0_05"],
            "GSEA_positive_NES_significant": gsea["positive_NES_significant_terms"],
            "GSEA_negative_NES_significant": gsea["negative_NES_significant_terms"],
            "GSEA_nonfinite_original_pvalue_genes_excluded": gsea["genes_excluded_for_nonfinite_original_pvalue"],
            "GSEA_finite_p_missing_padj_genes_retained": gsea["genes_retained_with_missing_adjusted_pvalue"],
        })
    return pd.DataFrame(rows)


def load_references():
    REF.mkdir(parents=True, exist_ok=True)
    META.mkdir(parents=True, exist_ok=True)
    manifest_path = META / "source_manifest.json"
    previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
    artifacts = []
    for name, url in SOURCE_URLS.items():
        path = REF / name
        if not path.exists():
            temporary = path.with_suffix(path.suffix + ".part")
            subprocess.run(["curl", "--fail", "--location", "--retry", "3", "--connect-timeout", "30",
                            "--max-time", "300", "--silent", "--show-error", "--output", str(temporary), url], check=True)
            os.replace(temporary, path)
        fingerprint = sha(path)
        if previous:
            old = next((r for r in previous["artifacts"] if r["name"] == name), None)
            if old and old["sha256"] != fingerprint:
                raise ValueError(f"Frozen GO reference changed: {name}; restore the recorded source or explicitly create a separately versioned analysis.")
        if name.endswith("gz"):
            header = []
            with gzip.open(path, "rt") as handle:
                for line in handle:
                    if not line.startswith("!"): break
                    header.append(line.rstrip())
            if not header or not header[0].startswith("!gaf-version:"):
                raise ValueError("Invalid GAF header")
        else:
            header = path.read_text().split("[Term]", 1)[0].strip().splitlines()
            if not header[0].startswith("format-version:"): raise ValueError("Invalid OBO header")
        artifacts.append({"name": name, "url": url, "path": str(path.relative_to(ROOT)),
                          "bytes": path.stat().st_size, "sha256": fingerprint, "header": header})
    manifest = {"retrieved_at": previous.get("retrieved_at") if previous else datetime.now(timezone.utc).isoformat(),
                "artifacts": artifacts, "scope": "Official GO current distribution files, frozen locally by SHA256. Header dates describe each delivered product; annotation source headers can reference older ontology builds.",
                "documentation": ["https://geneontology.org/docs/download-go-annotations/",
                                  "https://geneontology.org/docs/download-ontology/",
                                  "https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/"]}
    write_json(manifest_path, manifest)
    return manifest


def parse_go():
    terms, alt, current = {}, {}, None
    def store(term):
        if term and "id" in term:
            terms[term["id"]] = term
            for value in term["alt_ids"]: alt[value] = term["id"]
    with (REF / "go-basic.obo").open() as handle:
        for raw in handle:
            line = raw.rstrip()
            if line == "[Term]":
                store(current)
                current = {"parents": [], "alt_ids": [], "replaced_by": [], "obsolete": False}
            elif line.startswith("["):
                store(current); current = None
            elif current is not None:
                if line.startswith("id: "): current["id"] = line[4:]
                elif line.startswith("name: "): current["name"] = line[6:]
                elif line.startswith("namespace: "): current["namespace"] = line[11:]
                elif line.startswith("alt_id: "): current["alt_ids"].append(line[8:])
                elif line.startswith("is_a: "): current["parents"].append(line[6:].split()[0])
                elif line.startswith("relationship: part_of "): current["parents"].append(line.split()[2])
                elif line == "is_obsolete: true": current["obsolete"] = True
                elif line.startswith("replaced_by: "): current["replaced_by"].append(line[13:])
    store(current)
    def resolve(term):
        term = alt.get(term, term)
        seen = set()
        while term in terms and terms[term]["obsolete"]:
            if term in seen or len(terms[term]["replaced_by"]) != 1: return None
            seen.add(term)
            term = alt.get(terms[term]["replaced_by"][0], terms[term]["replaced_by"][0])
        return term if term in terms and terms[term].get("namespace") == "biological_process" else None
    bp = {key: value for key, value in terms.items() if not value["obsolete"] and value.get("namespace") == "biological_process"}
    @lru_cache(maxsize=None)
    def ancestors(term):
        output = {term}
        for parent in bp[term]["parents"]:
            parent = resolve(parent)
            if parent in bp: output.update(ancestors(parent))
        return frozenset(output)
    canonical, aliases, direct = set(), defaultdict(set), defaultdict(set)
    annotations = Counter()
    evidence = Counter()
    with gzip.open(REF / "mgi.gaf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("!"): continue
            annotations["total_gaf_rows"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 17: raise ValueError("Unexpected GAF column count")
            if fields[12].split("|")[0] != "taxon:10090":
                annotations["non_mouse_rows_excluded"] += 1; continue
            if "NOT" in fields[3].split("|"):
                annotations["negative_NOT_rows_excluded"] += 1; continue
            symbol = fields[2]
            canonical.add(symbol)
            for synonym in fields[10].split("|"):
                if synonym: aliases[synonym].add(symbol)
            if fields[8] != "P": continue
            annotations["positive_mouse_BP_rows"] += 1
            term = resolve(fields[4])
            if term is None:
                annotations["unresolvable_or_non_BP_term_rows"] += 1; continue
            if term != fields[4]: annotations["alt_id_or_single_replaced_by_rows_resolved"] += 1
            direct[symbol].add(term)
            evidence[fields[6]] += 1
    gene_terms = {}
    term_genes = defaultdict(set)
    for symbol, ids in direct.items():
        propagated = set()
        for term in ids: propagated.update(ancestors(term))
        gene_terms[symbol] = propagated
        for term in propagated: term_genes[term].add(symbol)
    gmt = REF / "mouse_BP_is_a_part_of.gmt"
    with gmt.open("w") as handle:
        for term in sorted(term_genes):
            handle.write("\t".join([term, bp[term]["name"], *sorted(term_genes[term])]) + "\n")
    audit = {"gaf_counts": dict(annotations), "evidence_codes_retained": dict(evidence),
             "canonical_mouse_symbols": len(canonical), "unique_aliases": sum(len(v) == 1 for v in aliases.values()),
             "ambiguous_aliases": sum(len(v) > 1 for v in aliases.values()),
             "active_BP_terms_in_ontology": len(bp), "BP_terms_with_mouse_annotations": len(term_genes),
             "genes_with_BP_annotations_including_root": len(gene_terms),
             "propagation_relations": ["is_a", "part_of"],
             "ontology_ID_policy": "Resolve alternative IDs and obsolete IDs with exactly one replacement leading to an active BP term; ambiguous or missing replacements are excluded.",
             "synonym_field": "GAF column 11 (zero-based index 10); column 10 is the descriptive name and is never used as a symbol alias.",
             "derived_gmt": {"path": str(gmt.relative_to(ROOT)), "bytes": gmt.stat().st_size, "sha256": sha(gmt)}}
    write_json(META / "annotation_audit.json", audit)
    return bp, canonical, aliases, gene_terms, term_genes, audit


def map_de(path, canonical, aliases, gene_terms, output):
    data = pd.read_csv(path, sep="\t")
    if "gene" not in data.columns: raise ValueError(f"Required gene column is absent: {path}")
    required = ["baseMean", "log2FoldChange", "stat", "pvalue", "padj"]
    if any(c not in data for c in required): raise ValueError(f"DE schema must include {required}: {path}")
    data["feature_id"] = data["gene_id"].astype(str) if "gene_id" in data else data["gene"].astype(str)
    data["source_gene"] = data.gene
    mapped, status, candidates = [], [], []
    for value in data.gene:
        symbol = str(value).strip() if pd.notna(value) else ""
        if symbol in canonical:
            mapped.append(symbol); status.append("direct_canonical_symbol"); candidates.append(symbol)
        elif len(aliases.get(symbol, set())) == 1:
            target = next(iter(aliases[symbol]))
            mapped.append(target); status.append("unique_GAF_synonym"); candidates.append(target)
        else:
            mapped.append(None); status.append("ambiguous_synonym" if aliases.get(symbol) else "unmapped_symbol")
            candidates.append("|".join(sorted(aliases.get(symbol, set()))))
    data["canonical_gene"] = mapped
    data["mapping_status"] = status
    data["mapping_candidates"] = candidates
    # Selection is independent of direction, test statistic, p-value or significance.
    mapped_data = data[data.canonical_gene.notna()].sort_values(["baseMean", "feature_id"], ascending=[False, True], kind="stable")
    representative_ids = mapped_data.drop_duplicates("canonical_gene", keep="first").index
    data["selected_gene_representative"] = data.index.isin(representative_ids)
    data["in_ORA_universe"] = data.selected_gene_representative & np.isfinite(data.padj)
    data["in_GSEA_rank_before_pvalue_QC"] = data.selected_gene_representative & np.isfinite(data.stat)
    data["in_GSEA_rank"] = data.in_GSEA_rank_before_pvalue_QC & np.isfinite(data.pvalue)
    data["GSEA_exclusion_reason"] = np.select(
        [data.canonical_gene.isna(), ~data.selected_gene_representative,
         ~np.isfinite(data.stat), ~np.isfinite(data.pvalue)],
        ["unmapped_or_ambiguous_symbol", "lower_baseMean_duplicate_representative",
         "nonfinite_Wald_statistic", "nonfinite_original_Wald_pvalue_Cooks_or_unreliable_test"], default="included")
    data["has_specific_BP_annotation"] = data.canonical_gene.map(lambda x: bool(set(gene_terms.get(x, set())) - {"GO:0008150"}))
    save_table(data, output / "gene_mapping.tsv.gz")
    save_table(data.loc[data.in_GSEA_rank_before_pvalue_QC & ~data.in_GSEA_rank], output / "GSEA_excluded_for_nonfinite_original_pvalue.tsv")
    duplicate_symbols = mapped_data.canonical_gene[mapped_data.canonical_gene.duplicated(keep=False)].unique()
    save_table(data[data.canonical_gene.isin(duplicate_symbols)].sort_values(["canonical_gene", "baseMean"], ascending=[True, False]), output / "duplicate_symbol_candidates.tsv")
    selected = data[data.selected_gene_representative].set_index("canonical_gene", drop=False)
    coverage = {"input_features": len(data), "input_unique_nonmissing_symbols": int(data.gene.nunique()),
                "feature_mapping_status_counts": data.mapping_status.value_counts().to_dict(),
                "mapped_features": int(data.canonical_gene.notna().sum()), "unique_mapped_gene_representatives": len(selected),
                "mapped_duplicate_symbols": len(duplicate_symbols), "features_discarded_by_gene_deduplication": int(data.canonical_gene.notna().sum() - len(selected)),
                "unmapped_significant_feature_symbols": data.loc[data.padj.lt(FDR) & data.canonical_gene.isna(), "source_gene"].dropna().tolist(),
                "ORA_universe_genes": int(data.in_ORA_universe.sum()),
                "ORA_universe_with_specific_BP_annotation": int((data.in_ORA_universe & data.has_specific_BP_annotation).sum()),
                "GSEA_rank_genes": int(data.in_GSEA_rank.sum()),
                "GSEA_rank_before_pvalue_QC_genes": int(data.in_GSEA_rank_before_pvalue_QC.sum()),
                "GSEA_genes_removed_for_nonfinite_original_pvalue": int((data.in_GSEA_rank_before_pvalue_QC & ~data.in_GSEA_rank).sum()),
                "GSEA_genes_retained_with_missing_adjusted_pvalue": int((data.in_GSEA_rank & ~np.isfinite(data.padj)).sum()),
                "GSEA_rank_with_nonfinite_Wald_pvalue": int((data.in_GSEA_rank & ~np.isfinite(data.pvalue)).sum()),
                "source_significant_features": int(data.padj.lt(FDR).sum()),
                "ORA_significant_mapped_gene_representatives": int(selected.padj.lt(FDR).sum()),
                "Pcnx_mapping": data.loc[data.source_gene.eq("Pcnx"), ["feature_id", "source_gene", "canonical_gene", "mapping_status", "has_specific_BP_annotation"]].to_dict("records")}
    return selected, coverage


def ora_table(selected, universe, eligible_sets, terms, direction):
    foreground = set(selected.index[selected.padj.lt(FDR) & (selected.log2FoldChange.gt(0) if direction == "up" else selected.log2FoldChange.lt(0))]) & universe
    rows = []
    N, n = len(universe), len(foreground)
    for term in sorted(eligible_sets):
        members = eligible_sets[term]
        hits = foreground & members
        M, k = len(members), len(hits)
        pvalue = float(hypergeom.sf(k-1, N, M, n)) if n else 1.
        rows.append({"GO_ID": term, "term_name": terms[term]["name"], "direction": direction,
                     "universe_genes": N, "term_genes_in_universe": M, "foreground_genes": n,
                     "overlap": k, "overlap_genes": ";".join(sorted(hits)),
                     "expected_overlap": n*M/N if N else np.nan,
                     "fold_enrichment": k/(n*M/N) if n and M else np.nan, "pvalue": pvalue})
    frame = pd.DataFrame(rows)
    frame["padj_BH"] = multipletests(frame.pvalue, method="fdr_bh")[1]
    return frame.sort_values(["padj_BH", "pvalue", "GO_ID"]), foreground


def export(fig, directory, stem):
    for extension in ("png", "svg", "pdf"):
        fig.savefig(directory / f"{stem}.{extension}", dpi=300, facecolor="white")
    plt.close(fig)


def draw_ora(frames, directory, dataset):
    fig, axes = plt.subplots(1, 2, figsize=(13, 7), layout="constrained")
    for ax, direction in zip(axes, ("up", "down")):
        frame = frames[direction]
        significant = frame[frame.padj_BH.lt(FDR)].head(10).iloc[::-1]
        n = int(frame.foreground_genes.iloc[0])
        ax.set_title(f"{direction.capitalize()}regulated genes (n={n})")
        if significant.empty:
            ax.text(.5, .5, f"No GO-BP term at BH FDR <0.05\n{len(frame):,} eligible terms tested\nincluding zero-overlap terms", transform=ax.transAxes, ha="center", va="center", fontsize=12)
            ax.set_axis_off(); continue
        values = -np.log10(significant.padj_BH.clip(lower=1e-300))
        ax.scatter(values, range(len(significant)), color=COLORS[direction], s=55)
        labels = [textwrap.fill(f"{r.term_name} ({r.overlap}/{n})", 44) for r in significant.itertuples()]
        ax.set(yticks=range(len(significant)), yticklabels=labels, xlabel="-log10(BH-adjusted p-value)")
        ax.tick_params(axis="y", labelsize=8)
        ax.axvline(-np.log10(FDR), lw=.7, ls="--", color="#60676D")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(left=0)
    fig.suptitle(f"{dataset}: GO Biological Process overrepresentation", fontsize=14)
    export(fig, directory, "ORA_significant_terms")


def draw_gsea(table, output, dataset):
    significant = table[table.fdr_qvalue.lt(FDR)]
    shown = (significant if len(significant) else table).sort_values(["fdr_qvalue", "nominal_pvalue"]).head(12).sort_values("NES")
    fig, ax = plt.subplots(figsize=(10.5, 7.8), layout="constrained")
    ax.scatter(shown.NES, range(len(shown)), c=[COLORS["up"] if n > 0 else COLORS["down"] for n in shown.NES], s=55)
    ax.axvline(0, color="#555555", lw=.7)
    labels = [textwrap.fill(f"{row.term_name} [{row.GO_ID}]", 54) for row in shown.itertuples()]
    ax.set(yticks=range(len(shown)), yticklabels=labels, xlabel="Normalized enrichment score (positive: exposed-side rank)")
    ax.tick_params(axis="y", labelsize=8)
    header = f"{dataset}: exploratory preranked GO-BP GSEA"
    detail = f"{len(significant)} terms at GSEA algorithmic FDR <0.05" if len(significant) else "No term at algorithmic FDR <0.05; lowest-q terms shown descriptively"
    ax.set_title(header + "\n" + detail, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    export(fig, output / "figures", "GSEA_terms")


def enrichment_curves(pre, table, sets, output):
    rank = pre.ranking
    selected = []
    for positive in (True, False):
        options = table[table.NES.gt(0) if positive else table.NES.lt(0)].sort_values(["fdr_qvalue", "nominal_pvalue"])
        if len(options): selected.append(options.iloc[0])
    for row in selected:
        members = sets[row.GO_ID]
        hits = rank.index.isin(members)
        weights = np.abs(rank.values.astype(float)) * hits
        if weights.sum() == 0: weights = hits.astype(float)
        running = np.cumsum(weights/weights.sum() - (~hits)/(len(rank)-hits.sum()))
        source = pd.DataFrame({"rank_position": np.arange(1, len(rank)+1), "gene": rank.index,
                               "wald_statistic": rank.values, "gene_set_member": hits, "running_enrichment_score": running})
        stem = "GSEA_curve_" + row.GO_ID.replace(":", "_")
        save_table(source, output / f"{stem}.tsv.gz")
        fig, axes = plt.subplots(2, 1, figsize=(8, 5), height_ratios=[4, 1], sharex=True, layout="none")
        fig.subplots_adjust(left=.15, right=.98, bottom=.16, top=.83, hspace=.08)
        color = COLORS["up"] if row.NES > 0 else COLORS["down"]
        axes[0].plot(source.rank_position, running, color=color, lw=1.5)
        axes[0].axhline(0, color="black", lw=.6)
        axes[0].set(ylabel="Running enrichment score", title=textwrap.fill(row.term_name, 72) + f"\nNES={row.NES:.2f}; algorithmic FDR q={row.fdr_qvalue:.3g}")
        axes[1].vlines(source.rank_position[hits], 0, 1, color=color, lw=.5)
        axes[1].set(yticks=[], xlabel="Gene rank (descending exposed versus control Wald statistic)")
        axes[0].spines[["top", "right"]].set_visible(False)
        axes[1].spines[["top", "right", "left"]].set_visible(False)
        export(fig, output / "figures", stem)


def analyze_dataset(name, path, terms, canonical, aliases, gene_terms, term_genes, threads):
    output = OUT / name
    output.mkdir(parents=True, exist_ok=True)
    (output / "figures").mkdir(exist_ok=True)
    selected, coverage = map_de(path, canonical, aliases, gene_terms, output)
    universe = set(selected.index[np.isfinite(selected.padj)])
    eligible_sets = {term: genes & universe for term, genes in term_genes.items() if MIN_SIZE <= len(genes & universe) <= MAX_SIZE}
    if not universe or not eligible_sets: raise ValueError(f"No eligible GO universe / terms: {name}")
    save_table(pd.DataFrame({"gene": sorted(universe), "has_specific_BP_annotation": [bool(gene_terms.get(g, set()) - {"GO:0008150"}) for g in sorted(universe)]}), output / "ORA_universe.tsv")
    ora, foregrounds = {}, {}
    for direction in ("up", "down"):
        ora[direction], foregrounds[direction] = ora_table(selected, universe, eligible_sets, terms, direction)
        save_table(ora[direction], output / f"ORA_{direction}_all_terms.tsv")
        save_table(ora[direction][ora[direction].padj_BH.lt(FDR)], output / f"ORA_{direction}_significant.tsv")
    save_table(pd.DataFrame({"GO_ID": sorted(eligible_sets), "term_name": [terms[x]["name"] for x in sorted(eligible_sets)],
                             "genes_in_ORA_universe": [len(eligible_sets[x]) for x in sorted(eligible_sets)]}), output / "eligible_GO_terms.tsv")
    draw_ora(ora, output / "figures", name)
    rank = selected.loc[np.isfinite(selected.stat) & np.isfinite(selected.pvalue), ["stat"]].copy()
    rank["gene"] = rank.index
    rank = rank.sort_values(["stat", "gene"], ascending=[False, True], kind="stable")
    rank_table = rank[["gene", "stat"]]
    previous_rank = output / "GSEA_ranked_genes.tsv"
    if previous_rank.exists() and pd.read_csv(previous_rank, sep="\t").gene.tolist() != rank_table.gene.tolist():
        audit = output / "quality_audit"
        audit.mkdir(exist_ok=True)
        for old_name in ("GSEA_ranked_genes.tsv", "GSEA_all_terms.tsv", "GSEA_significant.tsv", "summary.json"):
            old = output / old_name
            archived = audit / ("initial_including_nonfinite_pvalue_" + old_name)
            if old.exists() and not archived.exists(): shutil.copy2(old, archived)
        old_raw = output / "gseapy_raw"
        if old_raw.exists() and not (audit / "initial_gseapy_raw").exists():
            os.replace(old_raw, audit / "initial_gseapy_raw")
        initial_figures = audit / "initial_figures"
        initial_figures.mkdir(exist_ok=True)
        for old in list((output / "figures").glob("GSEA*")):
            target = initial_figures / old.name
            if not target.exists(): os.replace(old, target)
        for old in list(output.glob("GSEA_curve_*.tsv.gz")):
            target = audit / old.name
            if not target.exists(): os.replace(old, target)
        write_json(audit / "README.json", {"status": "superseded_diagnostic_results_do_not_use_as_main_analysis",
                   "reason": "Initial exploratory GSEA ranks included finite Wald statistics with nonfinite original p-values. The final analysis excludes these unreliable tests, including Cook-filtered genes, while retaining genes with finite raw p-values and missing adjusted p-values from independent filtering.",
                   "final_rank_criterion": "finite Wald statistic AND finite original Wald p-value; no FDR or effect-size selection"})
    save_table(rank_table, output / "GSEA_ranked_genes.tsv")
    # Recompute sizes in the complete finite-statistic rank, retaining the GO family
    # selected independently by finite-adjusted-p-value testability above.
    rank_set = set(rank.gene)
    gsea_sets = {term: term_genes[term] & rank_set for term in eligible_sets
                 if MIN_SIZE <= len(term_genes[term] & rank_set) <= MAX_SIZE}
    labels = {f"{term}|{terms[term]['name']}": sorted(genes) for term, genes in sorted(gsea_sets.items())}
    with (output / "GO_BP_tested_gene_sets.gmt").open("w") as handle:
        for term, genes in sorted(gsea_sets.items()):
            handle.write("\t".join([term, terms[term]["name"], *sorted(genes)]) + "\n")
    print(f"{name}: ORA universe={len(universe)}, terms={len(eligible_sets)}, foreground up={len(foregrounds['up'])}, down={len(foregrounds['down'])}.", flush=True)
    print(f"{name}: GSEA ranking={len(rank)} genes; {len(gsea_sets)} GO-BP sets; {PERMUTATIONS} permutations; {threads} threads.", flush=True)
    pre = gseapy.prerank(rnk=rank_table, gene_sets=labels, organism="mouse", outdir=str(output / "gseapy_raw"),
                        permutation_num=PERMUTATIONS, min_size=MIN_SIZE, max_size=MAX_SIZE, weight=1,
                        ascending=None, threads=threads, seed=SEED, no_plot=True, graph_num=0,
                        verbose=True, method="permutation", pheno_pos="exposed-side", pheno_neg="control-side")
    result = pre.res2d.copy()
    result["GO_ID"] = result.Term.str.split("|", n=1).str[0]
    result["term_name"] = result.Term.str.split("|", n=1).str[1]
    result = result.rename(columns={"ES": "ES", "NES": "NES", "NOM p-val": "nominal_pvalue",
                                    "FDR q-val": "fdr_qvalue", "FWER p-val": "FWER_pvalue",
                                    "Lead_genes": "leading_edge_genes", "Tag %": "tag_percent", "Gene %": "gene_percent"})
    for col in ("ES", "NES", "nominal_pvalue", "fdr_qvalue", "FWER_pvalue"):
        result[col] = pd.to_numeric(result[col], errors="raise")
    result["rank_member_genes"] = result.GO_ID.map(lambda x: len(gsea_sets[x]))
    result["FDR_method"] = "GSEA permutation-based normalized-ES FDR estimate (not Benjamini-Hochberg)"
    result = result.sort_values(["fdr_qvalue", "nominal_pvalue", "GO_ID"])
    save_table(result, output / "GSEA_all_terms.tsv")
    significant = result[result.fdr_qvalue.lt(FDR)]
    save_table(significant, output / "GSEA_significant.tsv")
    draw_gsea(result, output, name)
    enrichment_curves(pre, result, gsea_sets, output)
    stats = {"dataset": name, "DE_source": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
             "DE_source_sha256": sha(path), "annotation_mapping": coverage,
             "ORA": {"universe_size": len(universe), "terms_per_test_family": len(eligible_sets),
                     "up_foreground_size": len(foregrounds["up"]), "down_foreground_size": len(foregrounds["down"]),
                     "significant_up_terms_BH_FDR_0_05": int(ora["up"].padj_BH.lt(FDR).sum()),
                     "significant_down_terms_BH_FDR_0_05": int(ora["down"].padj_BH.lt(FDR).sum()),
                     "up_foreground_genes": sorted(foregrounds["up"]), "down_foreground_genes": sorted(foregrounds["down"]),
                     "foreground_genes_without_specific_BP_annotation": {direction: sorted(g for g in foregrounds[direction] if not gene_terms.get(g, set()) - {"GO:0008150"}) for direction in ("up", "down")},
                     "top_up_terms": ora["up"].loc[ora["up"].overlap.gt(0)].head(10).to_dict("records"), "top_down_terms": ora["down"].loc[ora["down"].overlap.gt(0)].head(10).to_dict("records"),
                     "interpretation": ""},
             "GSEA": {"ranking_size": len(rank), "tested_terms": len(result), "permutations": PERMUTATIONS,
                      "rank_quality_filter": "finite Wald statistic AND finite original p-value; adjusted p-value may be NA; no FDR/effect-size selection",
                      "genes_excluded_for_nonfinite_original_pvalue": coverage["GSEA_genes_removed_for_nonfinite_original_pvalue"],
                      "genes_retained_with_missing_adjusted_pvalue": coverage["GSEA_genes_retained_with_missing_adjusted_pvalue"],
                      "seed": SEED, "ranking_tied_score_fraction": float(rank.stat.duplicated(keep=False).mean()),
                      "significant_terms_algorithmic_FDR_0_05": len(significant),
                      "positive_NES_significant_terms": int(significant.NES.gt(0).sum()),
                      "negative_NES_significant_terms": int(significant.NES.lt(0).sum()),
                      "zero_nominal_permutation_pvalues": int(result.nominal_pvalue.eq(0).sum()),
                      "zero_estimated_FDR_values": int(result.fdr_qvalue.eq(0).sum()),
                      "top_terms": result.head(15).to_dict("records"),
                      "interpretation": "Exploratory gene-set permutation analysis of Wald-statistic ranks; gene-gene correlations and biological-replicate uncertainty are not represented by gene-set permutations. Algorithmic q-values are not ORA BH q-values. A term association does not demonstrate pathway activation, suppression or causality."},
             "completed_at": datetime.now(timezone.utc).isoformat()}
    if stats["ORA"]["significant_up_terms_BH_FDR_0_05"] + stats["ORA"]["significant_down_terms_BH_FDR_0_05"] == 0:
        stats["ORA"]["interpretation"] = "No GO-BP term met BH FDR <0.05 in either direction. Small significant-gene foregrounds can provide very little power; absence of enrichment does not show absence of a biological effect."
    else:
        stats["ORA"]["interpretation"] = "Overrepresentation of GO annotations among FDR-selected genes relative to the mapped testable gene universe. Related GO terms overlap and do not constitute independent biological discoveries."
    write_json(output / "summary.json", stats)
    print(f"{name}: finished; ORA significant up/down={stats['ORA']['significant_up_terms_BH_FDR_0_05']}/{stats['ORA']['significant_down_terms_BH_FDR_0_05']}; exploratory GSEA significant terms={len(significant)}.", flush=True)
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", action="append", default=[], metavar="NAME=DE_TABLE", help="Add or override a dataset; tables must contain gene, baseMean, log2FoldChange, stat, pvalue and padj")
    parser.add_argument("--only", action="append", help="Analyze only the named dataset(s), including newly supplied ozone contrasts")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    threads = min(4, max(1, args.threads))
    datasets = dict(DEFAULT_DATASETS)
    for item in args.dataset:
        if "=" not in item: raise ValueError("--dataset must use NAME=PATH")
        name, path = item.split("=", 1)
        if not re.fullmatch(r"[A-Za-z0-9_-]+", name): raise ValueError("Dataset names may contain letters, numbers, underscores or hyphens")
        datasets[name] = Path(path).resolve()
    if args.only:
        datasets = {name: datasets[name] for name in args.only}
    for name, path in datasets.items():
        if not path.exists(): raise FileNotFoundError(f"DE table not yet available for {name}: {path}")
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none", "savefig.transparent": False})
    sources = load_references()
    print("Parsing official mouse GAF and propagating GO BP through is_a / part_of.", flush=True)
    terms, canonical, aliases, gene_terms, term_genes, audit = parse_go()
    summaries = {}
    combined = OUT / "summary.json"
    if combined.exists(): summaries = json.loads(combined.read_text()).get("datasets", {})
    for name, path in datasets.items():
        summaries[name] = analyze_dataset(name, path, terms, canonical, aliases, gene_terms, term_genes, threads)
        write_json(combined, {"datasets": summaries, "reference_source_manifest": str((META / "source_manifest.json").relative_to(ROOT)),
                              "annotation_audit": audit, "updated_at": datetime.now(timezone.utc).isoformat()})
        save_table(summary_count_table(summaries), OUT / "summary_counts.tsv")
    methods = {
        "scope": "Mouse GO Biological Process enrichment for primary gene-level RNA-seq contrasts. miRNA features are intentionally excluded because no validated target mapping was supplied.",
        "reference_sources": sources, "annotation_audit": audit,
        "annotation_filter": "Taxon of annotated entity must be 10090; exclude NOT relations; biological_process aspect P only for gene sets; all positive annotation evidence codes retained.",
        "propagation": "Transitive closure through active BP is_a and part_of parents only; excludes regulates, positively_regulates, negatively_regulates, has_part and cross-aspect edges.",
        "gene_identifiers": "Case-sensitive exact canonical mouse symbols or uniquely resolving GAF DB Object Synonyms (column 11 / index 10). No uppercasing, orthology transfer, free-text name matching or human pathway libraries.",
        "gene_deduplication": "Map all DE features first; retain the representative feature with highest baseMean for each canonical symbol, regardless of statistic, p-value, direction or significance. Ties resolved by feature_id lexical order. All candidate decisions are exported.",
        "ORA_universe": "Unique mapped canonical gene representatives with finite adjusted DE p-value, including those without specific GO-BP annotations. Unmapped or ambiguous symbols are excluded and their coverage reported.",
        "ORA_foreground": "DE BH-adjusted p-value <0.05; up/down from the unshrunken log2 fold change sign, separately. No exploratory effect-size threshold.",
        "term_size_filter": {"minimum": MIN_SIZE, "maximum": MAX_SIZE, "basis": "Intersect each propagated GO-BP gene set with the ORA testable gene universe before selecting the eligible family. GSEA also enforces 15-500 members in its complete finite-statistic rank."},
        "ORA_test": "One-sided hypergeometric survival P(X>=observed overlap), using the mapped finite-adjusted-p-value gene universe. Terms with zero overlap remain in the tested family, with p=1.",
        "ORA_multiplicity": "Benjamini-Hochberg across all eligible GO-BP terms, including zero overlaps, separately within each dataset and up/down family; FDR threshold 0.05.",
        "GSEA": {"implementation": "GSEApy prerank", "ranking": "All unique mapped representative genes with finite unshrunken Wald statistic AND finite original Wald p-value, descending; ties resolved by gene symbol. Genes with finite raw p-value and missing adjusted p-value from independent filtering remain eligible. No FDR or effect-size selection is used.",
                 "Cook_filter_quality_control": "Nonfinite original p-values indicate Cook-filtered or otherwise unreliable DE tests and are excluded even when their Wald statistic remains finite. Excluded genes and initial exploratory ranks are preserved in the quality audit. This is a gene-test quality filter; no animal is deleted.",
                 "method": "permutation", "permutation_type": "gene-set / gene-label permutations, not phenotype permutations", "permutations": PERMUTATIONS,
                 "weight": 1, "seed": SEED, "threads": threads,
                 "FDR": "GSEA algorithmic FDR estimate based on normalized enrichment scores and permutation null distributions, not Benjamini-Hochberg. Threshold 0.05; nominal permutation p=0 or q=0 is a finite-permutation estimate, not zero error probability.",
                 "interpretation": "Exploratory association with the ranked contrast; no pathway activation, suppression or causal claim. Gene-set permutations do not account for sample-level uncertainty or preserve gene-gene correlation."},
        "caveats": ["GO annotations are incomplete and unequal across genes; mapping and specific BP annotation coverage are exported.",
                    "Gene-length and detectability biases are not modeled by this classical ORA. The measured/testable universe limits but does not remove these biases.",
                    "GO terms share genes and ancestors; enriched terms are correlated and should not be counted as independent mechanisms.",
                    "Small biological sample sizes and sample heterogeneity in the source DE analyses carry into enrichment.",
                    "Input gene lists, mapping and tested term families are fixed before each enrichment test; no search over thresholds or removal of inconvenient samples is performed."],
        "figure_captions": {"ORA_significant_terms": "Only GO-BP terms meeting within-family BH FDR <0.05 are plotted. Parenthetical counts are overlap / foreground. Empty panels explicitly report no significant terms.",
                            "GSEA_terms": "At most twelve lowest-q significant terms; when none pass, lowest-q terms are shown and labelled descriptive. Point position is NES and colors indicate NES direction.",
                            "GSEA_curve": "Running weighted enrichment and gene-set hit positions along the QC-eligible finite-Wald-statistic rank with finite original p-value. One lowest-q term per NES direction is shown for exploration; its estimated algorithmic FDR is displayed."},
        "versions": {pkg: importlib.metadata.version(pkg) for pkg in ("gseapy", "numpy", "pandas", "scipy", "statsmodels", "matplotlib")},
        "python": platform.python_version(), "script": "scripts/analyze_enrichment.py", "script_sha256": sha(Path(__file__)),
        "documentation": ["https://geneontology.org/docs/download-go-annotations/", "https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/", "https://geneontology.org/docs/download-ontology/", "https://gseapy.readthedocs.io/en/latest/run.html#prerank"]}
    write_json(OUT / "methods.json", methods)


if __name__ == "__main__":
    main()
