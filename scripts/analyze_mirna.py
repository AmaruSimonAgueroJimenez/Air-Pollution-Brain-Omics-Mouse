#!/usr/bin/env python3
"""Prepare GSE93967 probe-level signals, run limma, and export English figures."""
from __future__ import annotations
import gzip
import io
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone
import importlib.metadata

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from analyze_cortex import ROOT, COLORS, style, export, sample_plots, heatmap, sha256, write_json


def table_between(text: str, start: str, end: str) -> pd.DataFrame:
    content = text.split(start + "\n", 1)[1].split(end, 1)[0]
    return pd.read_csv(io.StringIO(content), sep="\t", dtype={"ID": str, "ID_REF": str})


def main() -> None:
    np.random.seed(20261018)
    style()
    out = ROOT / "results/mirna"
    out.mkdir(parents=True, exist_ok=True)
    source = ROOT / "data/GSE93967/processed/GSE93967_series_matrix.txt.gz"
    soft_path = ROOT / "metadata/GSE93967/GSE93967_family.soft"
    with gzip.open(source, "rt") as handle:
        raw = table_between(handle.read(), "!series_matrix_table_begin", "!series_matrix_table_end").set_index("ID_REF")
    soft = soft_path.read_text()
    annotation = table_between(soft, "!platform_table_begin", "!platform_table_end").dropna(subset=["ID"]).set_index("ID")
    if not raw.index.is_unique or not annotation.index.is_unique:
        raise ValueError("Duplicate probe IDs; no implicit collapsing is permitted")
    annotation = annotation.reindex(raw.index).rename(columns={"miRNA_ID": "miRNA", "SPOT_ID": "platform_spot_id"})
    annotation.index.name = "probe_id"
    # Exact supplied annotation identifies the biological miRNA rows. Controls are preserved separately.
    biological = annotation.miRNA.fillna("").str.startswith("mmu-")
    valid = np.isfinite(raw.to_numpy()).all(axis=1)
    variable = raw.var(axis=1).gt(0)
    selected = biological & valid & variable
    annotation["retained"] = selected
    annotation["exclusion_reason"] = np.where(~biological, "platform_control_or_unannotated", np.where(~valid, "nonfinite_signal", np.where(~variable, "zero_variance", "")))
    annotation.to_csv(out / "probe_annotation.tsv", sep="\t")
    metadata = pd.read_csv(ROOT / "metadata/GSE93967/samples.csv").set_index("gsm").loc[raw.columns].copy()
    conditions = np.where(metadata.title.str.contains("1 mg/kg", regex=False), "PM1",
                          np.where(metadata.title.str.contains("5 mg/kg", regex=False), "PM5", "Control"))
    metadata["condition"] = pd.Categorical(conditions, categories=["Control", "PM1", "PM5"])
    metadata["dose_mg_per_kg"] = [0 if c == "Control" else 1 if c == "PM1" else 5 for c in conditions]
    metadata["tissue"] = "hippocampus"
    metadata.index.name = "sample"
    if metadata.condition.value_counts().to_dict() != {"Control": 3, "PM1": 3, "PM5": 3}:
        raise ValueError("Unexpected group allocation")
    metadata.to_csv(out / "sample_metadata.tsv", sep="\t")
    raw.to_csv(out / "published_normalized_signal.tsv", sep="\t", index_label="probe_id")
    values = raw.loc[selected].copy()
    values.to_csv(out / "analysis_signal.tsv", sep="\t", index_label="probe_id")
    detection = pd.DataFrame(index=raw.index, columns=raw.columns, dtype=str)
    for part in soft.split("^SAMPLE = ")[1:]:
        sample = part.splitlines()[0].strip()
        if sample in detection.columns:
            tab = table_between(part, "!sample_table_begin", "!sample_table_end").set_index("ID_REF")
            detection[sample] = tab.Detection.reindex(detection.index)
    detection.to_csv(out / "detection_calls.tsv", sep="\t", index_label="probe_id")
    qc = pd.DataFrame(index=metadata.index)
    qc["condition"] = metadata.condition
    qc["signal_min"] = values.min()
    qc["signal_median"] = values.median()
    qc["signal_max"] = values.max()
    qc["negative_signal_features"] = values.lt(0).sum()
    qc["present_miRNA_features"] = detection.loc[selected].eq("P").sum()
    qc["marginal_miRNA_features"] = detection.loc[selected].eq("M").sum()
    qc["absent_miRNA_features"] = detection.loc[selected].eq("A").sum()
    qc["features_equal_to_sample_minimum"] = values.eq(values.min()).sum()
    qc.to_csv(out / "sample_qc.tsv", sep="\t")
    subprocess.run(["Rscript", str(ROOT / "scripts/analyze_mirna.R"), str(ROOT)], check=True)
    features = annotation.loc[selected, ["miRNA"]].copy()
    features["n_samples_Present"] = detection.loc[selected].eq("P").sum(axis=1)
    for condition in metadata.condition.cat.categories:
        group_samples = metadata.index[metadata.condition == condition]
        features[f"n_Present_{condition}"] = detection.loc[selected, group_samples].eq("P").sum(axis=1)
        features[f"mean_{condition}"] = values[group_samples].mean(axis=1)
    features["low_detection_flag"] = features[[f"n_Present_{c}" for c in metadata.condition.cat.categories]].max(axis=1).lt(3)
    features.to_csv(out / "miRNA_feature_summary.tsv", sep="\t")
    captions = {}
    scores, variance = sample_plots(values, metadata, qc, out, "Hippocampal miRNA: saline and PM doses", "published normalized signal", captions)
    summaries = []
    for contrast in ("PM1_vs_Control", "PM5_vs_Control"):
        test = pd.read_csv(out / f"limma_{contrast}.tsv", sep="\t").set_index("probe_id")
        test = features.join(test).sort_values(["padj", "pvalue"])
        test["significant_FDR05"] = test.padj.lt(.05)
        test["direction"] = np.where(test.effect_mean_difference > 0, "higher_in_PM", "lower_in_PM")
        test.to_csv(out / f"DE_{contrast}.tsv", sep="\t")
        dose = "1 mg/kg" if contrast.startswith("PM1_") else "5 mg/kg"
        fig, ax = plt.subplots(figsize=(7.8, 5.2), layout="constrained")
        sig = test.significant_FDR05
        for selection, color, marker in [(~sig, "#969696", "o"), (sig, "#D55E00", "^")]:
            ax.scatter(test.loc[selection, "effect_mean_difference"], -np.log10(test.loc[selection, "pvalue"].clip(lower=np.finfo(float).tiny)),
                       s=14, color=color, marker=marker, alpha=.65, rasterized=True)
        ax.axvline(0, color="black", lw=.8)
        ax.set(xlabel="Mean difference in published normalized signal (PM minus saline)",
               ylabel="-log10(moderated t-test p value)", title=f"miRNA: PM {dose} versus saline")
        for probe, row in test.head(1).iterrows():
            ax.annotate(str(row.miRNA), (row.effect_mean_difference, -np.log10(max(row.pvalue, np.finfo(float).tiny))),
                        fontsize=7, xytext=(3, 3), textcoords="offset points")
        export(fig, out, f"volcano_{contrast}", f"Probe-level limma moderated t tests for PM {dose} minus saline, n=3 biological samples per group. Effect units are published normalized signal units, not an asserted log2 fold change. Orange triangles: BH FDR <0.05 across tested miRNA rows for this contrast; gray circles: remaining rows. Only the top-ranked feature is labeled. No additional log transform or normalization was applied.", [f"DE_{contrast}.tsv"], captions)
        fig, ax = plt.subplots(figsize=(7.2, 5.2), layout="constrained")
        for selection, color, marker in [(~sig, "#969696", "o"), (sig, "#D55E00", "^")]:
            ax.scatter(test.loc[selection, "published_signal_mean"], test.loc[selection, "effect_mean_difference"], s=14,
                       color=color, marker=marker, alpha=.6, rasterized=True)
        ax.axhline(0, color="black", lw=.8)
        ax.set(xlabel="Mean published normalized signal", ylabel="Mean difference (PM minus saline)",
               title=f"miRNA: mean–difference plot, PM {dose}")
        export(fig, out, f"mean_difference_{contrast}", "Mean published normalized signal versus group mean difference; linear units are retained, including negative values. Orange triangles mark BH FDR <0.05; gray circles: remaining rows. This is a mean–difference display, not a conventional log-ratio MA plot.", [f"DE_{contrast}.tsv"], captions)
        heatmap(values, test.reset_index(), "probe_id", "miRNA", out, f"top30_heatmap_{contrast}",
                f"miRNA: top 30 ranked features, PM {dose}", captions, f"limma BH-adjusted p value for {contrast}")
        summaries.append({"contrast": contrast, "n_control": 3, "n_exposed": 3, "features_tested": len(test),
                          "FDR05_features": int(sig.sum()), "FDR05_higher": int((sig & test.effect_mean_difference.gt(0)).sum()),
                          "FDR05_lower": int((sig & test.effect_mean_difference.lt(0)).sum()),
                          "FDR05_low_detection_flag": int((sig & test.low_detection_flag).sum()),
                          "effect_units": "mean difference in published normalized signal"})
    omnibus = pd.read_csv(out / "limma_omnibus.tsv", sep="\t").set_index("probe_id")
    omnibus = features.join(omnibus).sort_values(["padj", "pvalue"])
    omnibus.to_csv(out / "DE_omnibus.tsv", sep="\t")
    heatmap(values, omnibus.reset_index(), "probe_id", "miRNA", out, "top30_heatmap_omnibus",
            "miRNA: top 30 features across PM doses", captions, "two-degree-of-freedom moderated F-test BH-adjusted p value")
    pd.DataFrame(summaries).to_csv(out / "de_summary.tsv", sep="\t", index=False)
    summary = {"dataset": "GSE93967", "tissue": "hippocampus", "n_samples": 9, "n_per_group": 3,
               "features_input": len(raw), "platform_controls_or_unannotated_removed": int((~biological).sum()),
               "nonfinite_or_constant_miRNA_removed": int((biological & ~selected).sum()), "miRNA_features_tested": len(values),
               "scale": "Published quantile-normalized signal; log2 preprocessing is not explicitly confirmed by GEO, so effects are mean differences in published signal units",
               "contrasts": summaries, "omnibus_FDR05_features": int(omnibus.padj.lt(.05).sum()),
               "PCA_PC1_percent": float(variance[0]), "PCA_PC2_percent": float(variance[1]),
               "low_detection_features_flagged_not_removed": int(features.low_detection_flag.sum()),
               "negative_signal_values_retained": int(values.lt(0).sum().sum()), "samples_excluded": [],
               "QC_notes": ["The deposited matrix contains a large sample-specific lower-bound pile-up near -3.3; negative values and detection flags are retained explicitly.", "PCA does not cleanly separate the three groups. No sample was removed based on its position in PCA.", "Low-detection features remain in the test family, so nominal or adjusted significance must be interpreted together with detection flags."],
               "interpretation": "Exploratory probe-level secondary analysis; low detection is explicitly flagged. The two pairwise contrast families and the omnibus family have separate BH correction; no pooled FDR guarantee is claimed."}
    write_json(out / "summary.json", summary)
    methods = {"created_at_UTC": datetime.now(timezone.utc).isoformat(), "source": str(source.relative_to(ROOT)),
               "source_sha256": sha256(source), "GEO_SOFT_sha256": sha256(soft_path), "platform": "GPL17912",
               "annotation": "Exact probe ID to supplied miRNA_ID mapping; historical miRBase v19 labels retained; no silent identifier update or collapse",
               "exclusions": "Non-mmu annotated platform controls, nonfinite and zero-variance miRNA rows only. All nine biological samples retained.",
               "normalization": "None added; retain deposited quantile-normalized GeneSpring signal unchanged",
               "scale_caution": "GEO documents quantile normalization but does not explicitly document log2 units. Signed deposited values are retained; no log transform, count model or fold-change interpretation is imposed.",
               "model": "limma lmFit with ~0+condition (Control, PM1, PM5), contrasts PM1-Control and PM5-Control, eBayes(trend=TRUE, robust=FALSE)",
               "multiplicity": "Benjamini-Hochberg adjusted p within each pairwise contrast across all retained miRNA rows; separately within the two-coefficient omnibus moderated F test. No joint correction across contrasts is claimed.",
               "uncertainty": "SE_moderated and pointwise 95% t intervals use posterior residual variance and moderated total degrees of freedom; intervals are not simultaneous or multiplicity-adjusted",
               "PCA": "Top 500 variable retained miRNA features; feature centered, not standardized; deposited signal units",
               "detection": "Agilent P/M/A detection calls retained. Low-detection flag: fewer than three Present calls in every group. Detection flags do not enter significance tests or remove miRNA rows.",
               "controls": "21 non-miRNA/control rows remain available in published_normalized_signal.tsv and probe_annotation.tsv but do not enter differential analysis",
               "references": ["https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93967", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17912", "https://doi.org/10.1186/s12989-017-0215-3", "https://doi.org/10.1093/nar/gkv007"],
               "software": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "matplotlib", "scikit-learn"]},
               "R_software": "See R_session_info.txt", "random_seed": 20261018, "limitations": ["n=3 biological samples per condition; no verified batch or pairing variables are available.", "Results do not establish miRNA-target interactions; no mRNA GO enrichment is inferred from miRNA names.", "Normalization and log-scale ambiguity cannot be resolved solely from the deposited normalized signal. Effect sizes are conservatively labeled in published signal units."]}
    write_json(out / "methods.json", methods)
    write_json(out / "figure_captions.json", captions)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
