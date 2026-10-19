#!/usr/bin/env python3
"""Reproduce a prespecified PM versus PBS analysis of GSE297423 ACC RNA-seq."""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parents[1]
COLORS = {"PBS": "#0072B2", "PM": "#D55E00", "Control": "#0072B2", "PM1": "#009E73", "PM5": "#D55E00"}
MARKERS = {"PBS": "o", "PM": "^", "Control": "o", "PM1": "s", "PM5": "^"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def style() -> None:
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white",
                         "axes.titleweight": "bold", "savefig.facecolor": "white"})


def export(fig, out: Path, name: str, caption: str, source_tables: list[str], captions: dict) -> None:
    directory = out / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(directory / f"{name}.{suffix}", dpi=300, facecolor="white")
    captions[name] = {"caption": caption, "source_tables": source_tables,
                      "formats": ["PNG at 300 dpi", "SVG with editable text", "PDF with embedded TrueType fonts"],
                      "size_inches": [float(v) for v in fig.get_size_inches()]}
    plt.close(fig)


def sample_plots(values: pd.DataFrame, metadata: pd.DataFrame, qc: pd.DataFrame,
                 out: Path, title: str, units: str, captions: dict) -> tuple[pd.DataFrame, np.ndarray]:
    """Use all available, finite, nonconstant features without selecting by DE outcome."""
    assert values.columns.equals(metadata.index)
    finite = np.isfinite(values.to_numpy()).all(axis=1)
    variable = values.var(axis=1) > 0
    usable = values.loc[finite & variable]
    if usable.shape[0] < 2:
        raise ValueError("Fewer than two variable finite features for PCA")
    n_top = min(500, len(usable))
    top = usable.var(axis=1).nlargest(n_top).index
    model = PCA(n_components=2, svd_solver="full")
    scores = pd.DataFrame(model.fit_transform(usable.loc[top].T), index=metadata.index, columns=["PC1", "PC2"])
    scores = scores.join(metadata)
    scores.index.name = "sample"
    scores.to_csv(out / "pca_scores.tsv", sep="\t")
    pd.Series(top, name="feature_id").to_csv(out / "pca_features.tsv", sep="\t", index=False)
    variance = model.explained_variance_ratio_ * 100
    fig, ax = plt.subplots(figsize=(7.6, 5.3), layout="constrained")
    for group in metadata.condition.drop_duplicates():
        selection = scores.condition == group
        ax.scatter(scores.loc[selection, "PC1"], scores.loc[selection, "PC2"], s=65,
                   color=COLORS[group], marker=MARKERS[group], label=f"{group} (n={selection.sum()})",
                   edgecolor="white", linewidth=0.6)
    for sample, row in scores.iterrows():
        ax.annotate(str(sample), (row.PC1, row.PC2), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set(xlabel=f"PC1 ({variance[0]:.1f}% variance)", ylabel=f"PC2 ({variance[1]:.1f}% variance)",
           title=f"{title}: sample PCA")
    ax.legend(frameon=False, loc="best")
    export(fig, out, "pca", f"PCA of the {n_top} most variable features on {units}. Features are centered, not variance-scaled. Each point is one biological sample. No sample was removed by PCA.",
           ["pca_scores.tsv", "pca_features.tsv", "sample_metadata.tsv"], captions)

    corr = usable.corr(method="pearson")
    corr.to_csv(out / "sample_correlation.tsv", sep="\t", index_label="sample")
    fig, ax = plt.subplots(figsize=(7, 6), layout="constrained")
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r", interpolation="nearest")
    ax.set_xticks(range(len(corr)), corr.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(corr)), corr.index, fontsize=8)
    for i in range(len(corr)):
        for j in range(len(corr)):
            ax.text(j, i, f"{corr.iloc[i, j]:.3f}", ha="center", va="center", fontsize=7,
                    color="white" if abs(corr.iloc[i, j]) > .65 else "black")
    fig.colorbar(im, ax=ax, label="Pearson correlation", shrink=.8)
    ax.set_title(f"{title}: sample correlations")
    export(fig, out, "sample_correlation", f"Pearson correlations across all {len(usable)} finite, variable features on {units}; fixed color limits from -1 to 1.", ["sample_correlation.tsv"], captions)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
    groups = metadata.condition.tolist()
    boxes = axes[0].boxplot([values[s].to_numpy() for s in metadata.index], tick_labels=metadata.index,
                            showfliers=False, patch_artist=True)
    for patch, group in zip(boxes["boxes"], groups):
        patch.set_facecolor(COLORS[group]); patch.set_alpha(.55)
    axes[0].tick_params(axis="x", rotation=60, labelsize=8)
    axes[0].set(ylabel=units, title="Feature distributions")
    if "library_size" in qc:
        metric, label = "library_size", "Total assigned gene counts (millions)"
        bar_values = qc[metric] / 1e6
    else:
        metric, label = "present_miRNA_features", "miRNA features called Present"
        bar_values = qc[metric]
    axes[1].bar(range(len(qc)), bar_values, color=[COLORS[g] for g in groups])
    axes[1].set_xticks(range(len(qc)), qc.index, rotation=60, ha="right", fontsize=8)
    axes[1].set(ylabel=label, title="Sample quality overview")
    fig.suptitle(title)
    export(fig, out, "sample_qc", "Left: distributions of the analysis values, with extreme boxplot points omitted only from this display. Right: sample-level assay quality summaries. All biological samples are retained; individual samples are identified by name.",
           ["sample_qc.tsv", "sample_metadata.tsv"], captions)
    return scores, variance


def heatmap(values: pd.DataFrame, ranking: pd.DataFrame, feature_col: str, label_col: str,
            out: Path, name: str, title: str, captions: dict, rank_description: str) -> None:
    selected = ranking.loc[ranking[feature_col].isin(values.index)].head(30).copy()
    x = values.loc[selected[feature_col]].copy()
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=1).replace(0, np.nan), axis=0).fillna(0)
    z.to_csv(out / f"{name}_zscores.tsv", sep="\t", index_label=feature_col)
    lim = max(float(np.abs(z.to_numpy()).max()), .1)
    fig, ax = plt.subplots(figsize=(8.4, max(5, .20 * len(z) + 2)), layout="constrained")
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim), interpolation="nearest")
    labels = selected[label_col].fillna(selected[feature_col]).astype(str)
    labels = [f"{label} [{identifier}]" if (labels == label).sum() > 1 else label
              for label, identifier in zip(labels, selected[feature_col])]
    ax.set_yticks(range(len(z)), labels, fontsize=8)
    ax.set_xticks(range(len(z.columns)), z.columns, rotation=55, ha="right", fontsize=8)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Within-feature z score", shrink=.65)
    export(fig, out, name, f"Top {len(z)} features ranked by {rank_description}; rows are centered and scaled by their sample standard deviation for display only. Selection is descriptive and is not independent validation. Fixed sample order; no clustering of samples.", [f"{name}_zscores.tsv"], captions)


def run(workers: int = 4) -> None:
    np.random.seed(20261018)
    style()
    out = ROOT / "results/cortex"
    out.mkdir(parents=True, exist_ok=True)
    source = ROOT / "data/GSE297423/processed/GSE297423_mouse_gene_count.txt.gz"
    source_metadata = ROOT / "metadata/GSE297423/samples.csv"
    data = pd.read_csv(source, sep="\t").set_index("gene_id")
    if not data.index.is_unique:
        raise ValueError("Ensembl gene identifiers are not unique; no silent aggregation is permitted")
    metadata = pd.read_csv(source_metadata)
    metadata["sample"] = metadata.description.str.extract(r"Library name: (.+)")
    chosen = ["Control5", "Control6", "Control7", "Control8", "PMA2775", "PMA2776", "PMA2778", "PMA2779"]
    metadata = metadata.set_index("sample").loc[chosen].copy()
    metadata["condition"] = pd.Categorical(["PBS"] * 4 + ["PM"] * 4, categories=["PBS", "PM"])
    metadata["tissue"] = "anterior cingulate cortex"
    metadata["PFF_injection"] = False
    metadata.to_csv(out / "sample_metadata.tsv", sep="\t")
    counts = data[chosen]
    array = counts.to_numpy()
    if not np.isfinite(array).all() or (array < 0).any() or not np.equal(array, np.floor(array)).all():
        raise ValueError("Raw RNA-seq input contains missing, negative or noninteger values")
    counts = counts.astype("int64")
    keep = (counts >= 10).sum(axis=1) >= 4
    annotation = data[["gene_name", "gene_biotype", "gene_chr"]].rename(columns={"gene_name": "gene"})
    annotation.to_csv(out / "gene_annotation.tsv", sep="\t")
    filter_table = annotation.join(pd.DataFrame({"samples_with_count_ge_10": (counts >= 10).sum(axis=1), "retained": keep}))
    filter_table.to_csv(out / "gene_filter.tsv", sep="\t")
    dds = DeseqDataSet(counts=counts.loc[keep].T, metadata=metadata, design="~condition",
                       refit_cooks=True, n_cpus=workers, quiet=False)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "PM", "PBS"], alpha=.05, n_cpus=workers,
                       cooks_filter=True, independent_filter=True)
    stats.summary()
    # lfc_shrink modifies results_df; preserve the original Wald MLE statistics first.
    mle = stats.results_df.copy(deep=True)
    mle.to_csv(out / "DE_PM_vs_PBS_MLE.tsv", sep="\t", index_label="gene_id")
    coeffs = list(dds.varm["LFC"].columns)
    coeff = next(c for c in coeffs if c != "Intercept")
    stats.lfc_shrink(coeff=coeff)
    shrunk = stats.results_df[["log2FoldChange", "lfcSE"]].rename(columns={"log2FoldChange": "log2FoldChange_shrunk", "lfcSE": "lfcSE_shrunk"})
    result = annotation.join(mle, how="inner").join(shrunk)
    result["significant_FDR05"] = result.padj.lt(.05)
    result["direction"] = np.where(result.log2FoldChange > 0, "higher_in_PM", "lower_in_PM")
    result = result.sort_values(["padj", "pvalue"], na_position="last")
    result.index.name = "gene_id"
    result.to_csv(out / "DE_PM_vs_PBS.tsv", sep="\t")
    normalized = pd.DataFrame(dds.layers["normed_counts"].T, index=dds.var_names, columns=dds.obs_names)
    normalized.to_csv(out / "normalized_counts.tsv.gz", sep="\t", index_label="gene_id", compression="gzip")
    dds.vst(use_design=False)
    vst = pd.DataFrame(dds.layers["vst_counts"].T, index=dds.var_names, columns=dds.obs_names)
    vst.to_csv(out / "vst_counts.tsv.gz", sep="\t", index_label="gene_id", compression="gzip")
    qc = pd.DataFrame(index=metadata.index)
    qc["condition"] = metadata.condition
    qc["library_size"] = counts.sum()
    qc["detected_genes_count_gt_0"] = (counts > 0).sum()
    qc["retained_genes_count_ge_10"] = (counts.loc[keep] >= 10).sum()
    qc["size_factor"] = np.asarray(dds.obs["size_factors"])
    qc["mitochondrial_count_fraction"] = counts.loc[data.gene_chr.eq("MT")].sum() / qc.library_size
    qc["rRNA_count_fraction"] = counts.loc[data.gene_biotype.str.contains("rRNA", case=False, na=False)].sum() / qc.library_size
    qc.to_csv(out / "sample_qc.tsv", sep="\t")
    captions = {}
    scores, variance = sample_plots(vst, metadata, qc, out, "ACC RNA-seq: PM versus PBS", "blind variance-stabilized counts", captions)
    significant = result.significant_FDR05
    fig, ax = plt.subplots(figsize=(7, 5.2), layout="constrained")
    ax.scatter(result.loc[~significant, "baseMean"], result.loc[~significant, "log2FoldChange_shrunk"], s=7, color="#9A9A9A", alpha=.55, rasterized=True)
    ax.scatter(result.loc[significant, "baseMean"], result.loc[significant, "log2FoldChange_shrunk"], s=10, color="#D55E00", marker="^", alpha=.65, rasterized=True)
    ax.axhline(0, color="black", lw=.8)
    ax.set(xscale="log", xlabel="Mean size-factor normalized count", ylabel="Shrunken log2 fold change (PM / PBS)", title="ACC RNA-seq: MA plot")
    export(fig, out, "ma", "Mean normalized counts versus apeGLM-style shrunken log2 fold change. Orange triangles: BH-adjusted Wald p < 0.05; gray circles: all other retained genes, including genes without adjusted p values. Significance is derived from the preserved unshrunken Wald test.", ["DE_PM_vs_PBS.tsv"], captions)
    fig, ax = plt.subplots(figsize=(7, 5.2), layout="constrained")
    valid = result.pvalue.notna() & result.log2FoldChange.notna()
    plot = result.loc[valid]
    for selected, color, marker, size in [(~plot.significant_FDR05, "#9A9A9A", "o", 8), (plot.significant_FDR05, "#D55E00", "^", 12)]:
        ax.scatter(plot.loc[selected, "log2FoldChange"], -np.log10(plot.loc[selected, "pvalue"].clip(lower=np.finfo(float).tiny)), s=size,
                   color=color, marker=marker, alpha=.6, rasterized=True)
    ax.axvline(0, color="black", lw=.8)
    ax.set(xlabel="Unshrunken log2 fold change (PM / PBS)", ylabel="-log10(Wald p value)", title="ACC RNA-seq: volcano plot")
    for position, (gene_id, row) in enumerate(plot.head(3).iterrows()):
        ax.annotate(str(row.gene) if pd.notna(row.gene) else gene_id,
                    (row.log2FoldChange, -np.log10(max(row.pvalue, np.finfo(float).tiny))), fontsize=8,
                    xytext=((-40, -5) if position == 2 else (5, 5)), textcoords="offset points")
    export(fig, out, "volcano", "Each point is one retained gene with a finite Wald p value; x: unshrunken model coefficient, y: unadjusted Wald p value. Orange triangles indicate BH-adjusted p < 0.05; gray circles denote the remaining genes. Labels identify the top three adjusted-p-value ranked rows; no horizontal nominal-p cutoff is used to imply multiplicity-adjusted significance.", ["DE_PM_vs_PBS.tsv"], captions)
    heatmap(vst, result.reset_index(), "gene_id", "gene", out, "top30_heatmap", "ACC RNA-seq: top 30 ranked genes", captions,
            "adjusted Wald p value, then nominal p value (PM versus PBS)")
    valid_p = result.pvalue.notna()
    summary = {"dataset": "GSE297423", "tissue": "anterior cingulate cortex", "contrast": "PM versus PBS without PFF injection",
               "n_samples": 8, "n_PM": 4, "n_PBS": 4, "genes_input": len(data), "genes_retained": int(keep.sum()),
               "genes_with_Wald_pvalue": int(valid_p.sum()), "genes_with_adjusted_pvalue": int(result.padj.notna().sum()),
               "FDR05_genes": int(significant.sum()), "FDR05_higher_in_PM": int((significant & result.log2FoldChange.gt(0)).sum()),
               "FDR05_lower_in_PM": int((significant & result.log2FoldChange.lt(0)).sum()),
               "FDR05_abs_shrunken_log2FC_ge_1": int((significant & result.log2FoldChange_shrunk.abs().ge(1)).sum()),
               "Cook_filtered_genes": result.index[result.pvalue.isna()].tolist(),
               "library_size_range": [int(qc.library_size.min()), int(qc.library_size.max())],
               "QC_notes": ["PM samples show greater spread than PBS in the first two principal components; all samples were retained.", "Sample-level correlations are high on VST values and do not constitute proof of absence of batch effects or other confounding."],
               "PCA_PC1_percent": float(variance[0]), "PCA_PC2_percent": float(variance[1]),
               "samples_excluded_by_analysis": [], "interpretation": "Exploratory analysis of the supplied bulk gene-count matrix; n=4 biological samples per group. Effects are specific to the studied humanized alpha-synuclein mouse model and nasal PM exposure."}
    write_json(out / "summary.json", summary)
    pd.DataFrame([summary]).to_csv(out / "de_summary.tsv", sep="\t", index=False)
    methods = {"created_at_UTC": datetime.now(timezone.utc).isoformat(), "source": str(source.relative_to(ROOT)),
               "source_sha256": sha256(source), "sample_metadata_sha256": sha256(source_metadata),
               "samples": chosen, "design": "~condition", "contrast": ["condition", "PM", "PBS"],
               "gene_identifiers": "Unique supplied Ensembl IDs; supplied gene symbols retained as annotation; no collapse by symbol",
               "prefilter": "At least 10 counts in at least 4 of 8 selected samples", "normalization": "DESeq2 median-of-ratios size factors (ratio method)",
               "differential_expression": "PyDESeq2 negative-binomial GLM, two-sided Wald test; Benjamini-Hochberg FDR; alpha=0.05; independent filtering and Cook's filtering enabled",
               "outliers": "No manual sample exclusion. Default Cook's replacement needs >=7 replicates within a group; with n=4 per group replacement is unavailable and Cook's p-value filtering remains enabled.",
               "effect_size_shrinkage": "PyDESeq2 lfc_shrink with adapt=True (apeGLM-style heavy-tailed prior); original MLE coefficients, SEs, Wald statistics and p values are preserved before mutation",
               "LFC_coefficient": coeff, "PCA": "Blind DESeq2 variance-stabilizing transform, top 500 variable genes, feature centered without variance scaling",
               "heatmap": "Top 30 adjusted-p ranked retained genes; VST values standardized within gene; descriptive display, not validation",
               "multiplicity": "BH correction across retained genes after independent filtering within the single prespecified contrast",
               "units": {"counts": "nonnegative integer assigned gene counts", "log2FoldChange": "log2 PM/PBS ratio", "normalized_counts": "size-factor normalized assigned gene counts"},
               "workers": workers, "random_seed": 20261018, "software": {p: importlib.metadata.version(p) for p in ["pydeseq2", "numpy", "pandas", "scipy", "matplotlib", "scikit-learn"]},
               "references": ["https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297423", "https://doi.org/10.1126/science.adu4132", "https://doi.org/10.1186/s13059-014-0550-8", "https://doi.org/10.1093/bioinformatics/btad547"],
               "limitations": ["Only supplied sample annotations enter the design; no verified sex, batch or other unprovided covariate is invented.", "The analysis starts from the depositor's counts and does not rerun FASTQ alignment or read-level QC.", "No cell composition or spatial effects can be directly estimated from bulk RNA-seq alone."]}
    write_json(out / "methods.json", methods)
    write_json(out / "figure_captions.json", captions)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 6:
        parser.error("--workers must be between 1 and 6")
    run(args.workers)
