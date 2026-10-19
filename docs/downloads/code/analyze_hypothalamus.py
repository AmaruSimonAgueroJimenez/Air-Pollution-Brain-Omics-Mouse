#!/usr/bin/env python3
"""Reproducible bulk RNA-seq analysis of GSE145840 mouse hypothalamus.

Primary comparison: chow + PM (n=3) versus chow + filtered air (n=4).
All source files are read-only; results are written under results/hypothalamus.
Run using .venv-analysis/bin/python scripts/analyze_hypothalamus.py.
"""
from __future__ import annotations

import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/hypothalamus"
FIGURES = OUT / "figures"
COLORS = {"FA": "#008579", "PM": "#D86F12"}
MARKERS = {"FA": "o", "PM": "^"}
ALPHA = 0.05
FIGURE_RECORDS = []


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean_json(value.tolist())
    if isinstance(value, np.generic):
        return clean_json(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if value is pd.NA:
        return None
    return value


def write_json(name, record):
    (OUT / name).write_text(json.dumps(clean_json(record), indent=2, allow_nan=False) + "\n")


def save_table(frame, name, index=True):
    frame.to_csv(OUT / name, sep="\t", index=index, na_rep="NA", float_format="%.10g")


def export_figure(fig, stem, caption, source_tables):
    for extension in ("png", "svg", "pdf"):
        fig.savefig(FIGURES / f"{stem}.{extension}", dpi=300, facecolor="white")
    FIGURE_RECORDS.append({"figure": stem, "caption": caption, "source_tables": source_tables,
                           "formats": ["png", "svg", "pdf"], "width_inches": float(fig.get_figwidth()),
                           "height_inches": float(fig.get_figheight()), "png_dpi": 300})
    plt.close(fig)


def group_legend(ax):
    handles = [Line2D([], [], marker=MARKERS[g], linestyle="", color=COLORS[g],
                      label=f"{g} (n={4 if g == 'FA' else 3})", markersize=7) for g in ("FA", "PM")]
    ax.legend(handles=handles, frameon=False)


def load_counts():
    source_meta = pd.read_csv(ROOT / "metadata/GSE145840/samples.csv")
    expected = [f"GSM{i}" for i in range(4337097, 4337104)]
    metadata = source_meta.set_index("gsm").loc[expected].copy()
    metadata.index.name = "sample"
    metadata["label"] = metadata["title"].str.replace("HYPO-", "", regex=False)
    metadata["condition"] = pd.Categorical(["FA"] * 4 + ["PM"] * 3, categories=["FA", "PM"])
    metadata["comparison_included"] = True
    files, series, sources = [], [], []
    for gsm, sample in metadata.iterrows():
        matches = list((ROOT / "data/GSE145840/processed").glob(f"{gsm}_Hypothalamus_*.txt.gz"))
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one source count file for {gsm}: {matches}")
        path = matches[0]
        frame = pd.read_csv(path, sep=r"\s+", names=["gene", "count"], header=None)
        if frame.gene.isna().any() or frame.gene.duplicated().any():
            raise ValueError(f"Missing or duplicate gene identifiers in {path.name}")
        values = pd.to_numeric(frame["count"], errors="raise")
        if values.isna().any() or (values < 0).any() or not np.equal(values, values.astype(np.int64)).all():
            raise ValueError(f"Counts must be finite nonnegative integers in {path.name}")
        values.index = frame.gene
        values.name = gsm
        series.append(values.astype(np.int64))
        files.append(path)
        sources.append({"sample": gsm, "label": sample["label"], "path": str(path.relative_to(ROOT)),
                        "bytes": path.stat().st_size, "sha256": sha256(path), "gene_rows": len(frame)})
    if any(set(s.index) != set(series[0].index) for s in series[1:]):
        raise ValueError("The seven primary samples do not have identical gene sets; missing values cannot be filled with zero.")
    counts = pd.concat(series, axis=1).loc[series[0].index]
    counts.index.name = "gene"
    if counts.shape != (50417, 7) or counts.isna().any().any():
        raise ValueError(f"Unexpected primary count matrix dimensions: {counts.shape}")
    return counts, metadata, sources


def make_figures(counts, normalized, metadata, qc, de, vst, pca_table, explained, corr, euclidean, shrink_ok):
    labels = metadata["label"].tolist()
    # Raw sample observations and zero-based bar axes preserve sample-level variation.
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    for ax, column, ylabel, scale in [(axes[0], "library_size", "Total gene-assigned counts (millions)", 1e6),
                                     (axes[1], "detected_genes", "Detected genes (count > 0)", 1)]:
        for i, (sample, row) in enumerate(metadata.iterrows()):
            ax.bar(i, qc.loc[sample, column] / scale, color=COLORS[row.condition], alpha=.78)
            ax.scatter(i, qc.loc[sample, column] / scale, color="black", marker=MARKERS[row.condition], s=18)
        ax.set(xticks=range(7), xticklabels=labels, ylabel=ylabel, ylim=(0, None))
        ax.tick_params(axis="x", rotation=45)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_title("Library size")
    axes[1].set_title("Gene detection")
    export_figure(fig, "sample_quality", "Raw gene-assigned library counts and detected genes for all seven samples. Each bar is one biological sample; FA n=4, PM n=3. No sample was excluded.", ["qc_samples.tsv"])

    fig, ax = plt.subplots(figsize=(6.6, 5.2), layout="constrained")
    pca_offsets = {"FA1": (7, 6), "FA2": (7, 6), "FA3": (24, 40), "FA4": (24, -29),
                   "PM1": (-43, 13), "PM2": (-57, 36), "PM3": (-43, -30)}
    for sample, row in metadata.iterrows():
        x, y = pca_table.loc[sample, ["PC1", "PC2"]]
        ax.scatter(x, y, color=COLORS[row.condition], marker=MARKERS[row.condition], s=75, edgecolor="white", linewidth=.7)
        ax.annotate(row.label, (x, y), xytext=pca_offsets[row.label], textcoords="offset points", fontsize=9,
                    arrowprops={"arrowstyle": "-", "color": COLORS[row.condition], "lw": .65},
                    ha="left" if pca_offsets[row.label][0] >= 0 else "right")
    ax.set(xlabel=f"PC1 ({explained[0]:.1%} variance)", ylabel=f"PC2 ({explained[1]:.1%} variance)", title="Hypothalamus: sample PCA")
    ax.margins(.18)
    ax.spines[["top", "right"]].set_visible(False)
    group_legend(ax)
    export_figure(fig, "pca", "PCA of the 500 highest-variance genes after DESeq2 variance stabilization; genes were mean-centered, without variance scaling. Each point is a biological sample. Sample separation is descriptive and does not establish a biological mechanism or batch effect.", ["pca.tsv", "pca_gene_variances.tsv.gz", "sample_metadata.tsv"])

    order = leaves_list(linkage(squareform(np.clip(1 - corr.values, 0, 2), checks=False), method="average"))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    arrays = [corr.values[np.ix_(order, order)], euclidean.values[np.ix_(order, order)]]
    for ax, values, name, cmap, barlabel in zip(axes, arrays, ["Sample Spearman correlation", "Euclidean distance after VST"], ["viridis", "magma"], ["Spearman rho", "Euclidean distance"]):
        im = ax.imshow(values, cmap=cmap, interpolation="nearest")
        ax.set(xticks=range(7), yticks=range(7), xticklabels=np.array(labels)[order], yticklabels=np.array(labels)[order], title=name)
        ax.tick_params(axis="x", rotation=45)
        for tick in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
            tick.set_color(COLORS["FA" if tick.get_text().startswith("FA") else "PM"])
        for i in range(7):
            for j in range(7):
                if name.startswith("Sample"):
                    value = values[i, j]
                    normed = im.norm(value)
                    ax.text(j, i, f"{value:.3f}", ha="center", va="center", fontsize=7, color="black" if normed > .60 else "white")
        fig.colorbar(im, ax=ax, label=barlabel, shrink=.74)
    export_figure(fig, "sample_similarity", "Left: Spearman correlation of DESeq2-normalized counts across the retained genes. Right: Euclidean distances using all variance-stabilized retained genes. Both panels use the same sample order obtained by average-linkage clustering of 1-rho. Group colors supplement explicit FA/PM labels.", ["sample_spearman_correlation.tsv", "sample_spearman_distance.tsv", "sample_vst_euclidean_distance.tsv"])

    estimate = "log2FoldChange_shrunk" if shrink_ok else "log2FoldChange"
    estimate_label = "Shrunken log2 fold change (PM / FA)" if shrink_ok else "MLE log2 fold change (PM / FA)"
    sig = de.padj.lt(ALPHA).fillna(False)
    valid = de[estimate].notna()
    fig, ax = plt.subplots(figsize=(7, 4.6), layout="constrained")
    ax.scatter(de.loc[valid & ~sig, "baseMean"], de.loc[valid & ~sig, estimate], s=5, color="#AEB7BE", alpha=.5, rasterized=True)
    for direction, color in [(1, COLORS["PM"]), (-1, COLORS["FA"])]:
        selected = valid & sig & (np.sign(de[estimate]) == direction)
        ax.scatter(de.loc[selected, "baseMean"], de.loc[selected, estimate], s=10, color=color, alpha=.8, rasterized=True)
    ax.axhline(0, color="black", lw=.7)
    ax.axhline(1, color="#767676", lw=.65, ls="--")
    ax.axhline(-1, color="#767676", lw=.65, ls="--")
    ax.set(xscale="log", xlabel="Mean normalized count (log scale)", ylabel=estimate_label,
           title=f"Hypothalamus: MA plot ({int(sig.sum())} genes at FDR < 0.05)")
    ax.spines[["top", "right"]].set_visible(False)
    export_figure(fig, "ma", "Mean normalized abundance versus effect estimate. Orange/teal points are significant increases/decreases at BH FDR <0.05; grey genes do not meet this threshold. Dashed lines indicate an exploratory absolute log2 effect of one. Shrinkage does not alter Wald p-values or BH adjusted p-values.", ["differential_expression.tsv"])

    fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
    y = -np.log10(de.padj.clip(lower=1e-300))
    testable = valid & y.notna()
    ax.scatter(de.loc[testable & ~sig, estimate], y[testable & ~sig], s=7, color="#ADB6BD", alpha=.45, rasterized=True)
    for direction, color in [(1, COLORS["PM"]), (-1, COLORS["FA"])]:
        selected = testable & sig & (np.sign(de[estimate]) == direction)
        ax.scatter(de.loc[selected, estimate], y[selected], s=11, color=color, alpha=.8, rasterized=True)
    ax.axhline(-np.log10(ALPHA), color="#555555", ls="--", lw=.8)
    ax.axvline(0, color="black", lw=.6)
    volcano_offsets = {"Hipk3": (6, 6), "Ksr2": (7, 22), "Tll1": (-22, 10), "Pcnx": (7, -12)}
    for gene in de.loc[testable].sort_values("padj").head(6).index:
        offset = volcano_offsets.get(gene, (6, 5))
        ax.annotate(gene, (de.loc[gene, estimate], y.loc[gene]), xytext=offset, textcoords="offset points", fontsize=8,
                    arrowprops={"arrowstyle": "-", "color": "#687078", "lw": .5})
    ax.set(xlabel=estimate_label, ylabel="-log10(BH-adjusted p-value)", title="Hypothalamus: differential expression")
    ax.spines[["top", "right"]].set_visible(False)
    export_figure(fig, "volcano", "Volcano plot uses the selected effect estimate on the x-axis and the unmodified two-sided Wald BH-adjusted p-value on the y-axis. Genes with missing adjusted p-values are omitted from this panel and retained as NA in the full table. Six smallest adjusted p-values are labelled; labels do not imply replication.", ["differential_expression.tsv"])

    top = de.sort_values(["padj", "pvalue"], na_position="last").head(30).index
    values = vst.loc[:, top].T
    z = values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    save_table(z, "heatmap_top30_gene_zscores.tsv")
    fig, ax = plt.subplots(figsize=(7, 8.7), layout="constrained")
    limit = max(float(np.nanmax(np.abs(z.values))), .01)
    image = ax.imshow(z.values, aspect="auto", interpolation="nearest", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
    ax.set(xticks=range(7), xticklabels=labels, yticks=range(30), yticklabels=top, title="Top 30 genes ranked by adjusted p-value")
    ax.tick_params(axis="y", labelsize=8)
    for tick, condition in zip(ax.get_xticklabels(), metadata.condition):
        tick.set_color(COLORS[condition])
    fig.colorbar(image, ax=ax, label="Within-gene VST z-score", shrink=.45)
    export_figure(fig, "heatmap_top30", "Thirty genes ranked by adjusted p-value, followed by raw p-value, shown after variance stabilization and within-gene mean centering / sample SD scaling. Fixed sample order FA1-FA4, PM1-PM3. This selection is outcome-dependent and descriptive; included genes are not necessarily all significant.", ["differential_expression.tsv", "heatmap_top30_gene_zscores.tsv"])

    targets = [g for g in ("Gfap", "Ahr") if g in normalized.index]
    if targets:
        fig, axes = plt.subplots(1, len(targets), figsize=(3.9 * len(targets), 4.2), squeeze=False, layout="constrained")
        for ax, gene in zip(axes[0], targets):
            for j, group in enumerate(("FA", "PM")):
                ids = metadata.index[metadata.condition == group]
                offsets = np.linspace(-.12, .12, len(ids))
                vals = normalized.loc[gene, ids]
                ax.scatter(j + offsets, vals, color=COLORS[group], marker=MARKERS[group], s=45, edgecolor="white", linewidth=.6)
                ax.hlines(vals.mean(), j-.20, j+.20, color="black", lw=1.5)
                for x, sid, val in zip(j + offsets, ids, vals):
                    ax.annotate(metadata.loc[sid, "label"], (x, val), xytext=(3, 3), textcoords="offset points", fontsize=7)
            text = f"BH FDR = {de.loc[gene, 'padj']:.3g}" if gene in de.index and pd.notna(de.loc[gene, "padj"]) else "BH FDR = NA (not testable / filtered)"
            ax.set(xticks=[0, 1], xticklabels=["FA (n=4)", "PM (n=3)"], ylabel="DESeq2-normalized count", title=f"{gene}\n{text}", ylim=(0, None))
            ax.spines[["top", "right"]].set_visible(False)
        export_figure(fig, "selected_gene_counts", "Preselected genes Gfap and Ahr, when available. Each labelled point is one biological sample; black lines show group means. Counts are divided by DESeq2 size factors. The displayed FDR comes from the genome-wide comparison, not a separate gene-specific test.", ["normalized_counts.tsv.gz", "sample_metadata.tsv", "differential_expression.tsv"])


def leave_one_out_sensitivity(filtered, metadata, primary, cpus):
    """Describe sensitivity to every single sample without changing the primary model."""
    records, run_summaries = [], []
    hit_genes = primary.index[primary.padj.lt(ALPHA)]
    for omitted in metadata.index:
        print(f"Sensitivity: omitting {metadata.loc[omitted, 'label']} (primary model unchanged).", flush=True)
        retained = metadata.index[metadata.index != omitted]
        inference = DefaultInference(n_cpus=cpus)
        model = DeseqDataSet(counts=filtered.loc[:, retained].T, metadata=metadata.loc[retained, ["condition"]],
                             design="~condition", refit_cooks=True, inference=inference, n_cpus=cpus, quiet=True)
        model.deseq2()
        test = DeseqStats(model, contrast=["condition", "PM", "FA"], alpha=ALPHA, cooks_filter=True,
                          independent_filter=True, inference=inference, n_cpus=cpus, quiet=True)
        test.summary()
        result = test.results_df.copy()
        result.index.name = "gene"
        result["omitted_sample"] = omitted
        result["omitted_label"] = metadata.loc[omitted, "label"]
        records.append(result.reset_index())
        common = pd.concat([primary.stat.rename("primary"), result.stat.rename("sensitivity")], axis=1).dropna()
        run_summaries.append({"omitted_sample": omitted, "omitted_label": metadata.loc[omitted, "label"],
                              "FA_samples": int((metadata.loc[retained, "condition"] == "FA").sum()),
                              "PM_samples": int((metadata.loc[retained, "condition"] == "PM").sum()),
                              "genes_testable": int(result.pvalue.notna().sum()),
                              "genes_fdr_below_0_05": int(result.padj.lt(ALPHA).sum()),
                              "primary_hits_still_fdr_below_0_05": int(result.loc[hit_genes, "padj"].lt(ALPHA).sum()),
                              "primary_hits_same_lfc_direction": int((np.sign(result.loc[hit_genes, "log2FoldChange"]) == np.sign(primary.loc[hit_genes, "log2FoldChange"])).sum()),
                              "wald_stat_spearman_with_primary": float(common.corr(method="spearman").iloc[0, 1])})
    long_table = pd.concat(records, ignore_index=True)
    save_table(long_table, "leave_one_out_all_genes.tsv.gz", index=False)
    save_table(pd.DataFrame(run_summaries), "leave_one_out_summary.tsv", index=False)
    stability = []
    for gene in hit_genes:
        values = long_table.loc[long_table.gene == gene]
        stability.append({"gene": gene, "primary_log2FoldChange": primary.loc[gene, "log2FoldChange"],
                          "primary_padj": primary.loc[gene, "padj"], "leave_one_out_runs": len(values),
                          "runs_same_lfc_direction": int((np.sign(values.log2FoldChange) == np.sign(primary.loc[gene, "log2FoldChange"])).sum()),
                          "runs_fdr_below_0_05": int(values.padj.lt(ALPHA).sum()),
                          "minimum_leave_one_out_lfc": float(values.log2FoldChange.min()),
                          "maximum_leave_one_out_lfc": float(values.log2FoldChange.max()),
                          "maximum_leave_one_out_padj": float(values.padj.max())})
    save_table(pd.DataFrame(stability), "leave_one_out_primary_hits.tsv", index=False)
    return {"status": "completed", "runs": run_summaries, "primary_hit_stability": stability,
            "method": "Seven separate condition-only models each omit one sample; the full-model prefiltered gene universe is held fixed. Each run re-estimates size factors and dispersions and applies the same Wald test and FDR settings. No sensitivity run replaces the primary comparison.",
            "interpretation": "Descriptive influence assessment, not independent validation. Reduced sample sizes and changes in dispersion estimation / Cook's filtering alter power and the number of testable genes; loss of significance alone does not establish an erroneous primary result."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpus", type=int, default=6)
    parser.add_argument("--skip-sensitivity", action="store_true", help="Skip systematic leave-one-out influence assessment")
    args = parser.parse_args()
    cpus = min(max(args.cpus, 1), 6)
    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                         "savefig.transparent": False})
    start = datetime.now(timezone.utc)
    print("Loading seven hypothalamus count files (FA=4; PM=3; HFD excluded).", flush=True)
    counts, metadata, sources = load_counts()
    keep = counts.ge(10).sum(axis=1).ge(3)
    filtered = counts.loc[keep]
    if len(filtered) != 19611:
        raise ValueError(f"Unexpected prefilter count: {len(filtered)} != 19611")
    save_table(counts, "counts.tsv.gz")
    save_table(metadata, "sample_metadata.tsv")
    save_table(pd.DataFrame({"retained": keep, "samples_count_at_least_10": counts.ge(10).sum(axis=1)}), "gene_filter.tsv.gz")
    qc = pd.DataFrame({"label": metadata.label, "condition": metadata.condition,
                       "library_size": counts.sum(axis=0), "detected_genes": counts.gt(0).sum(axis=0),
                       "retained_gene_counts": filtered.sum(axis=0), "retained_detected_genes": filtered.gt(0).sum(axis=0)})
    qc.index.name = "sample"
    infer = DefaultInference(n_cpus=cpus)
    print(f"Fitting negative-binomial model for {len(filtered)} genes with {cpus} CPUs.", flush=True)
    recorded_warnings = []
    with warnings.catch_warnings(record=True) as warning_records:
        warnings.simplefilter("always")
        dds = DeseqDataSet(counts=filtered.T, metadata=metadata[["condition"]], design="~condition",
                           refit_cooks=True, inference=infer, n_cpus=cpus, quiet=True)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=["condition", "PM", "FA"], alpha=ALPHA,
                           cooks_filter=True, independent_filter=True, inference=infer, n_cpus=cpus, quiet=True)
        stats.summary()
        mle = stats.results_df.copy(deep=True)
        mle.index.name = "gene"
        save_table(mle, "differential_expression_mle.tsv")
        print(f"Wald tests complete: {int(mle.padj.lt(ALPHA).sum())} genes at BH FDR < 0.05.", flush=True)
        shrink_status = "not_attempted"
        shrink_error = None
        shrunk = pd.DataFrame(index=mle.index, columns=["lfc_shrunk", "lfcSE_shrunk"], dtype=float)
        try:
            coef = [x for x in stats.LFC.columns if x != "Intercept"]
            if len(coef) != 1 or "PM" not in str(coef[0]):
                raise ValueError(f"Unexpected coefficient names for PM-vs-FA: {stats.LFC.columns.tolist()}")
            print(f"Estimating apeGLM shrinkage for {coef[0]}; preserving MLE tests.", flush=True)
            stats.lfc_shrink(coeff=coef[0], adapt=True)
            shrunk["lfc_shrunk"] = stats.results_df.log2FoldChange
            shrunk["lfcSE_shrunk"] = stats.results_df.lfcSE
            if hasattr(stats, "_LFC_shrink_converged"):
                shrunk["shrink_converged"] = stats._LFC_shrink_converged
                bad = shrunk.shrink_converged.fillna(False).eq(False)
                shrunk.loc[bad, ["lfc_shrunk", "lfcSE_shrunk"]] = np.nan
            if not np.isfinite(shrunk[["lfc_shrunk", "lfcSE_shrunk"]].values).any():
                raise ValueError("No finite converged shrinkage estimates")
            shrink_status = "completed" if shrunk.lfc_shrunk.notna().all() else "completed_with_nonconverged_genes_as_NA"
            if not mle[["pvalue", "padj"]].equals(stats.results_df[["pvalue", "padj"]]):
                raise ValueError("Shrinkage unexpectedly changed inferential p-values")
        except Exception as exc:
            shrink_status = "failed_mle_preserved"
            shrink_error = f"{type(exc).__name__}: {exc}"
            shrunk[["lfc_shrunk", "lfcSE_shrunk"]] = np.nan
        save_table(shrunk, "lfc_shrinkage.tsv")
        de = mle.join(shrunk[["lfc_shrunk", "lfcSE_shrunk"]].rename(columns={"lfc_shrunk": "log2FoldChange_shrunk", "lfcSE_shrunk": "lfcSE_shrunk"}))
        save_table(de, "differential_expression.tsv")
        qc["size_factor"] = dds.obs["size_factors"]
        normalized = counts.div(qc.size_factor, axis=1)
        save_table(normalized, "normalized_counts.tsv.gz")
        cooks = np.asarray(dds.layers["cooks"])
        cooks_threshold = float(f.ppf(.99, dds.obsm["design_matrix"].shape[1], len(metadata)-dds.obsm["design_matrix"].shape[1]))
        qc["max_cooks_distance"] = np.nanmax(cooks, axis=1)
        qc["genes_above_raw_cooks_cutoff"] = (cooks > cooks_threshold).sum(axis=1)
        qc["raw_cooks_cutoff"] = cooks_threshold
        qc["sample_excluded"] = False
        qc["eligible_for_cooks_count_replacement"] = False
        save_table(qc, "qc_samples.tsv")
        gene_diagnostics = dds.var.copy()
        gene_diagnostics["max_cooks_distance"] = np.nanmax(cooks, axis=0)
        gene_diagnostics["cooks_outlier_for_inference"] = np.asarray(dds.cooks_outlier())
        gene_diagnostics.index.name = "gene"
        save_table(gene_diagnostics, "gene_fit_diagnostics.tsv.gz")
        print("Variance stabilization and sample-level diagnostics.", flush=True)
        dds.vst(use_design=True)
        vst = pd.DataFrame(dds.layers["vst_counts"], index=metadata.index, columns=filtered.index)
        if not np.isfinite(vst.values).all():
            raise ValueError("Nonfinite values in the variance-stabilized matrix")
        save_table(vst.T, "vst_counts.tsv.gz")
        gene_variances = vst.var(axis=0, ddof=1).sort_values(ascending=False)
        pca_features = gene_variances.head(500).index
        pca_model = PCA(n_components=min(6, len(metadata)-1), svd_solver="full")
        scores = pca_model.fit_transform(vst.loc[:, pca_features])
        pca_table = pd.DataFrame(scores, index=metadata.index, columns=[f"PC{i+1}" for i in range(scores.shape[1])])
        pca_table.index.name = "sample"
        save_table(pca_table.join(metadata[["label", "condition"]]), "pca.tsv")
        save_table(gene_variances.rename("vst_variance").to_frame().assign(used_for_pca=lambda x: x.index.isin(pca_features)), "pca_gene_variances.tsv.gz")
        save_table(pd.DataFrame(pca_model.components_.T, index=pca_features, columns=[f"PC{i+1}" for i in range(scores.shape[1])]).rename_axis("gene"), "pca_loadings.tsv.gz")
        corr = normalized.loc[filtered.index].corr(method="spearman")
        corr.index.name = "sample"
        save_table(corr, "sample_spearman_correlation.tsv")
        save_table(1-corr, "sample_spearman_distance.tsv")
        euclidean = pd.DataFrame(squareform(pdist(vst.values, metric="euclidean")), index=metadata.index, columns=metadata.index)
        save_table(euclidean, "sample_vst_euclidean_distance.tsv")
        recorded_warnings.extend(sorted({f"{type(w.message).__name__}: {w.message}" for w in warning_records}))
    print("Writing figures and analysis summary.", flush=True)
    make_figures(counts, normalized, metadata, qc, de, vst, pca_table, pca_model.explained_variance_ratio_, corr, euclidean, shrink_status.startswith("completed"))
    sensitivity = ({"status": "not performed (explicit command-line option)"} if args.skip_sensitivity else
                   leave_one_out_sensitivity(filtered, metadata, mle, cpus))
    sig = de.padj.lt(ALPHA).fillna(False)
    large_effect = sig & de.log2FoldChange_shrunk.abs().ge(1)
    summary = {
        "dataset": "GSE145840", "tissue": "hypothalamus", "organism": "Mus musculus",
        "comparison": "PM versus FA (chow-fed samples only)", "biological_samples": 7,
        "group_sizes": {"FA": 4, "PM": 3}, "excluded_comparator": "HFD (4 samples; different diet)",
        "input_gene_rows": len(counts), "prefiltered_genes": len(filtered),
        "prefilter": "count >= 10 in at least 3 of 7 primary samples",
        "genes_with_finite_wald_pvalue": int(de.pvalue.notna().sum()),
        "genes_with_finite_adjusted_pvalue": int(de.padj.notna().sum()),
        "genes_missing_wald_pvalue": int(de.pvalue.isna().sum()),
        "genes_missing_adjusted_pvalue": int(de.padj.isna().sum()),
        "fdr_threshold": ALPHA, "genes_fdr_below_0_05": int(sig.sum()),
        "upregulated_fdr_below_0_05_mle_sign": int((sig & de.log2FoldChange.gt(0)).sum()),
        "downregulated_fdr_below_0_05_mle_sign": int((sig & de.log2FoldChange.lt(0)).sum()),
        "exploratory_fdr_and_abs_shrunken_lfc_at_least_1": int(large_effect.sum()),
        "shrinkage_status": shrink_status, "shrinkage_error": shrink_error,
        "finite_converged_shrunken_effects": int(de.log2FoldChange_shrunk.notna().sum()),
        "pca_genes": 500, "pca_variance_explained": pca_model.explained_variance_ratio_.tolist(),
        "cooks_threshold_99_percent_F": cooks_threshold,
        "genes_cooks_filtered_for_inference": int(gene_diagnostics.cooks_outlier_for_inference.sum()),
        "sample_exclusions": [], "count_replacement_min_replicates": 7,
        "count_replacement_performed": False,
        "normalization_size_factors": qc.size_factor.to_dict(),
        "top_genes_by_adjusted_pvalue": de.sort_values(["padj", "pvalue"]).head(20).reset_index().to_dict("records"),
        "selected_gene_results": {g: de.loc[g].to_dict() if g in de.index else {"status": "not retained by count prefilter"} for g in ("Gfap", "Ahr") if g in counts.index},
        "warnings": recorded_warnings,
        "leave_one_out_sensitivity": sensitivity,
        "observed_qc": "FA1 and FA2 are separated from the remaining five samples on the first two PCs; FA3/FA4 lie close to the PM samples. Thus the leading PCs do not show a clean exposure-only separation. All seven samples remain in the primary analysis.",
        "caveats": [
            "Only 3 PM and 4 FA biological samples are available; effect estimates, dispersion estimates and enrichment results have limited precision and power.",
            "The design includes condition only because the deposited sample metadata do not supply validated sample-level batch, age or sex covariates. Unmeasured confounding cannot be excluded.",
            "Samples were retained regardless of PCA position. PCA and correlation describe sample heterogeneity and do not prove a treatment mechanism or justify deleting an outlier.",
            "Cook's inference filtering is enabled. The default minimum of seven replicates per design group for count replacement is not met; no sample or extreme count was manually removed.",
            "The HFD comparator is excluded from this primary contrast. Its source files include a different gene list, and missing source entries must not be treated as measured zeros.",
            "Bulk tissue differences can reflect cell composition, cell-state changes or other factors; they do not establish a cell-specific effect or causality.",
            "Absolute shrunken log2 fold change >=1 is an exploratory effect-size summary, not an independently tested threshold.",
            "Source processed gene-count tables are analyzed as deposited; FASTQ alignment and gene quantification were not repeated."
        ],
        "started_at": start.isoformat(), "completed_at": datetime.now(timezone.utc).isoformat()
    }
    methods = {
        "analysis": "GSE145840 hypothalamus primary PM versus filtered-air contrast",
        "input_sources": sources, "source_metadata": "metadata/GSE145840/samples.csv",
        "sample_unit": "one deposited mouse hypothalamus RNA-seq library per biological replicate",
        "count_validation": "50417 unique gene symbols per sample; identical gene sets across seven primary samples; integer nonnegative counts; no missing values",
        "analysis_gene_universe": "19611 genes with count >=10 in >=3 primary samples",
        "design": "~condition", "reference_condition": "FA", "contrast": ["condition", "PM", "FA"],
        "estimation": "PyDESeq2 negative-binomial generalized linear model with gene-wise dispersion estimation, dispersion trend and MAP shrinkage",
        "size_factor_method": "DESeq2 median-of-ratios on prefiltered genes (ratio)",
        "normalization_export": "All 50417 source genes divided by the size factors estimated from prefiltered genes",
        "test": "Two-sided Wald test of log2 fold change = 0; default independent filtering; Benjamini-Hochberg FDR alpha=0.05",
        "cooks_filter": True, "refit_cooks": True, "min_replicates_for_count_replacement": 7,
        "manual_sample_or_gene_outlier_deletion": False,
        "lfc_shrinkage": "PyDESeq2 apeGLM-style adaptive Cauchy prior; failed/nonconverged estimates are NA. Original MLE log2FC, SE, statistic, p-value and adjusted p-value are saved independently before shrinkage.",
        "shrinkage_status": shrink_status,
        "pca": "DESeq2 VST with use_design=True, top 500 genes by VST variance; mean-centered genes, no variance scaling; full singular value decomposition",
        "sample_correlation": "Spearman correlation of normalized retained genes; correlation distance=1-rho; average-linkage ordering for figure only",
        "sample_distance": "Euclidean distance across all variance-stabilized retained genes",
        "heatmap": "Top 30 genes by adjusted p-value then raw p-value; within-gene VST mean centering and sample-SD scaling",
        "selected_genes": ["Gfap", "Ahr"], "cpu_process_limit": cpus, "blas_threads_per_process": 1,
        "hfd_analysis": "Not performed; excluded from the prespecified primary contrast",
        "leave_one_out_sensitivity": sensitivity["status"],
        "random_seed": "No stochastic method used (PCA full SVD, deterministic plotting offsets)",
        "versions": {pkg: importlib.metadata.version(pkg) for pkg in ("pydeseq2", "anndata", "numpy", "pandas", "scipy", "scikit-learn", "matplotlib")},
        "python": platform.python_version(), "platform": platform.platform(),
        "script": "scripts/analyze_hypothalamus.py", "script_sha256": sha256(Path(__file__)),
        "publisher_target": "Unspecified; provisional general scientific figures, not certified for a journal",
        "warnings": recorded_warnings
    }
    write_json("summary.json", summary)
    write_json("methods.json", methods)
    write_json("figure_manifest.json", {"figures": FIGURE_RECORDS})
    print(json.dumps({k: summary[k] for k in ("genes_fdr_below_0_05", "upregulated_fdr_below_0_05_mle_sign", "downregulated_fdr_below_0_05_mle_sign", "exploratory_fdr_and_abs_shrunken_lfc_at_least_1", "shrinkage_status", "genes_cooks_filtered_for_inference", "pca_variance_explained")}), flush=True)


if __name__ == "__main__":
    main()
