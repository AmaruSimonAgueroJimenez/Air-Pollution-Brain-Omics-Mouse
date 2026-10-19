"""Small display helpers for the English Quarto research report."""
from __future__ import annotations
import html
import json
from pathlib import Path
import pandas as pd
from IPython.display import Markdown, display

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent

DISPLAY_LABELS = {
    "gene_id": "Gene ID", "gene": "Gene", "gene_symbol": "Gene", "baseMean": "Mean normalized count",
    "log2FoldChange": "MLE log₂FC", "lfcSE": "MLE SE", "stat": "Test statistic",
    "log2FoldChange_shrunk": "Shrunken log₂FC", "lfcSE_shrunk": "Shrunken posterior SD",
    "pvalue": "Nominal p", "padj": "BH FDR", "padj_BH": "Term BH FDR",
    "padj_across_four_strata": "BH FDR across four contrasts", "fdr_qvalue": "GSEA estimated FDR",
    "effect_mean_difference": "Mean difference (signal units)", "published_signal_mean": "Mean published signal",
    "CI95_lower": "95% CI lower", "CI95_upper": "95% CI upper", "probe_id": "Probe ID",
    "sample": "Sample", "sample_id": "Sample", "label": "Animal label", "condition": "Group",
    "gsm": "GEO sample", "title": "Source sample title", "PFF_injection": "PFF injection",
    "sample_ids": "Sample IDs", "age": "Starting age", "genotype": "Genotype", "exposure": "Exposure",
    "class": "Atlas class", "substructure": "Substructure", "n_cells": "Cells",
    "positive_cells": "Positive-count cells", "mean_imaging_count": "Mean imaging count",
    "detection_percent": "Detection (%)", "percent_of_reference": "Reference cells (%)",
    "GO_ID": "GO ID", "term_name": "GO term", "NES": "NES", "leading_edge_genes": "Leading-edge genes",
    "overlap": "Overlapping genes", "fold_enrichment": "Fold enrichment",
    "primary_padj": "Primary BH FDR", "runs_same_lfc_direction": "Same direction (of 7 fits)",
    "runs_fdr_below_0_05": "FDR < 0.05 (of 7 fits)",
    "minimum_leave_one_out_lfc": "Lowest omitted-animal log₂FC",
    "maximum_leave_one_out_lfc": "Highest omitted-animal log₂FC",
    "omitted_sample": "Omitted sample", "omitted_label": "Omitted animal",
}


def summary(dataset):
    return json.loads((ROOT / "results" / dataset / "summary.json").read_text())


def read_table(dataset, filename, **kwargs):
    path = ROOT / "results" / dataset / filename
    kwargs.setdefault("sep", "," if ".csv" in path.name else "\t")
    return pd.read_csv(path, **kwargs)


def fmt_p(value):
    if value is None or pd.isna(value):
        return "unavailable"
    return f"{value:.2e}" if value < .01 else f"{value:.3f}"


def show_table(frame, n=20):
    view = frame.head(n).copy()
    for col in view.select_dtypes(include="floating"):
        if col.lower() in {"pvalue", "padj", "p_value", "fdr", "adjusted_p_value", "p.value", "adj.p.val"}:
            view[col] = view[col].map(fmt_p)
        else:
            view[col] = view[col].map(lambda v: "—" if pd.isna(v) else f"{v:,.3g}")
    view.columns = [DISPLAY_LABELS.get(str(col), str(col).replace("_", " ")) for col in view.columns]
    display(Markdown('<div class="table-scroll">' + view.to_html(index=False, border=0,
                 classes="table table-sm table-striped evidence-table", escape=True) + '</div>'))


def figure(dataset, stem, caption, alt=""):
    image = f"assets/figures/{dataset}/{stem}.png"
    svg = f"assets/figures/{dataset}/{stem}.svg"
    pdf = f"assets/figures/{dataset}/{stem}.pdf"
    display(Markdown('<figure class="research-figure"><a href="' + html.escape(svg) + '">'
                 '<img loading="lazy" src="' + html.escape(image) + '" alt="' +
                 html.escape(alt or caption) + '"></a><figcaption>' + html.escape(caption) +
                 ' <span class="figure-links"><a href="' + html.escape(svg) + '">SVG</a> · '
                 '<a href="' + html.escape(pdf) + '">PDF</a></span></figcaption></figure>'))


def metric_cards(items):
    display(Markdown('<div class="metric-cards">' + ''.join(
        '<div><span>' + html.escape(str(label)) + '</span><strong>' + html.escape(str(value)) +
        '</strong><small>' + html.escape(str(detail)) + '</small></div>'
        for label, value, detail in items) + '</div>'))
