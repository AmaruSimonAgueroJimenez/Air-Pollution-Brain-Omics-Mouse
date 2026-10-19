#!/usr/bin/env python3
"""Describe the measured HY reference without performing cell-level inference."""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/spatial_viewer"
OUT = ROOT / "results/spatial_reference"
BLUE = "#0072B2"
TEAL = "#008A80"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def export(fig, stem: str, caption: str, sources: list[str], records: dict) -> None:
    for extension in ("png", "svg", "pdf"):
        fig.savefig(OUT / "figures" / f"{stem}.{extension}", dpi=300,
                    bbox_inches="tight", facecolor="white", pad_inches=.15)
    records[stem] = {"caption": caption, "source_tables": sources,
                     "formats": ["300-dpi PNG", "editable-text SVG", "PDF"]}
    plt.close(fig)


def composition(meta: pd.DataFrame, group: str) -> pd.DataFrame:
    result = meta.groupby(group, sort=False).size().rename("n_cells").reset_index()
    result["fraction_of_reference"] = result.n_cells / len(meta)
    result["percent_of_reference"] = 100 * result.fraction_of_reference
    return result.sort_values(["n_cells", group], ascending=[False, True]).reset_index(drop=True)


def grouped_expression(counts: np.ndarray, meta: pd.DataFrame, genes: pd.DataFrame,
                       group: str, sizes: pd.DataFrame) -> pd.DataFrame:
    blocks = []
    for label in sizes[group]:
        values = counts[meta[group].eq(label).to_numpy()]
        positive = np.count_nonzero(values, axis=0)
        totals = values.sum(axis=0, dtype=np.float64)
        block = pd.DataFrame({group: label, "gene_id": genes.id.to_numpy(),
                              "gene": genes.symbol.to_numpy(), "n_cells": len(values),
                              "sum_imaging_counts": totals,
                              "mean_imaging_count": totals / len(values),
                              "positive_cells": positive,
                              "detection_fraction": positive / len(values)})
        blocks.append(block)
    return pd.concat(blocks, ignore_index=True)


def fraction_plot(frame: pd.DataFrame, group: str, title: str, stem: str,
                  caption: str, source: str, records: dict, footnote: str) -> None:
    y = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(9.7, .31 * len(frame) + 1.9), layout="constrained")
    percentages = 100 * frame.detection_fraction.to_numpy()
    ax.barh(y, percentages, color=TEAL, height=.7)
    for yi, row in enumerate(frame.itertuples(index=False)):
        ax.text(1.02, yi, f"{row.positive_cells:,} / {row.n_cells:,}",
                transform=ax.get_yaxis_transform(), va="center", fontsize=8.5)
    ax.text(1.02, 1.015, "Positive / all cells", transform=ax.transAxes, fontsize=9)
    ax.set_yticks(y, frame[group])
    ax.invert_yaxis()
    ax.set(xlim=(0, max(1, percentages.max() * 1.09)),
           xlabel="Cells with a positive measured Tll1 imaging count (%)", title=title)
    ax.grid(axis="x", color="#E8ECEE", linewidth=.6)
    ax.set_axisbelow(True)
    fig.supxlabel(footnote, fontsize=8.5, color="#465A64")
    export(fig, stem, caption, [source], records)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9.5, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "figure.facecolor": "white"})
    paths = {"counts": INPUT / "reference_counts.npz",
             "cells": INPUT / "registered_hypothalamus_cells.tsv.gz",
             "genes": INPUT / "gene_panel_coverage.tsv",
             "viewer_summary": INPUT / "summary.json"}
    with np.load(paths["counts"], allow_pickle=False) as bundle:
        counts = bundle["counts"]
        labels = bundle["cell_labels"].astype(str)
        ids = bundle["gene_ids"].astype(str)
    meta = pd.read_csv(paths["cells"], sep="\t", dtype={"cell_label": str})
    genes = pd.read_csv(paths["genes"], sep="\t", dtype={"id": str})
    viewer = json.loads(paths["viewer_summary"].read_text())
    assert viewer["reference_dataset"] == "Zhuang-ABCA-2" and viewer["reference_sex"] == "male"
    assert len(set(labels)) == len(labels) and meta.cell_label.is_unique
    assert set(labels) == set(meta.cell_label)
    meta = meta.set_index("cell_label").loc[labels].reset_index()
    assert meta.cell_label.to_numpy().tolist() == labels.tolist()
    assert genes.id.is_unique and genes.symbol.is_unique
    assert genes.id.to_numpy().tolist() == ids.tolist()
    assert counts.shape == (len(meta), len(genes)) == (37811, 1122)
    assert np.isfinite(counts).all() and (counts >= 0).all()
    assert np.equal(counts, np.floor(counts)).all()
    assert meta.division.eq("HY").all()
    assert meta[["class", "substructure", "structure"]].notna().all().all()
    assert meta.groupby("substructure").structure.nunique().max() == 1
    records = {}
    class_comp = composition(meta, "class")
    structure_comp = composition(meta, "substructure")
    structure_map = meta[["substructure", "structure"]].drop_duplicates()
    structure_comp = structure_comp.merge(structure_map, on="substructure", validate="one_to_one")
    class_comp.to_csv(OUT / "composition_by_class.tsv", sep="\t", index=False)
    structure_comp.to_csv(OUT / "composition_by_substructure.tsv", sep="\t", index=False)
    by_class = grouped_expression(counts, meta, genes, "class", class_comp)
    by_structure = grouped_expression(counts, meta, genes, "substructure", structure_comp)
    by_structure = by_structure.merge(structure_map, on="substructure", validate="many_to_one")
    by_class.to_csv(OUT / "all_genes_by_class.tsv.gz", sep="\t", index=False, float_format="%.12g")
    by_structure.to_csv(OUT / "all_genes_by_substructure.tsv.gz", sep="\t", index=False, float_format="%.12g")
    tll_class = by_class.loc[by_class.gene.eq("Tll1")].copy()
    tll_structure = by_structure.loc[by_structure.gene.eq("Tll1")].copy()
    tll_class.to_csv(OUT / "Tll1_by_class.tsv", sep="\t", index=False, float_format="%.12g")
    tll_structure.to_csv(OUT / "Tll1_by_substructure.tsv", sep="\t", index=False, float_format="%.12g")
    positive = np.count_nonzero(counts, axis=0)
    gene_summary = pd.DataFrame({"gene_id": ids, "gene": genes.symbol,
                                 "n_cells": len(meta),
                                 "sum_imaging_counts": counts.sum(axis=0, dtype=np.float64),
                                 "mean_imaging_count": counts.mean(axis=0, dtype=np.float64),
                                 "positive_cells": positive,
                                 "detection_fraction": positive / len(meta)})
    gene_summary.to_csv(OUT / "all_genes_reference_summary.tsv", sep="\t", index=False, float_format="%.12g")
    tll = gene_summary.loc[gene_summary.gene.eq("Tll1")].iloc[0]

    fig, ax = plt.subplots(figsize=(10.1, 10.2), layout="constrained")
    y = np.arange(len(class_comp))
    ax.barh(y, class_comp.n_cells, color=BLUE, height=.7)
    for yi, row in enumerate(class_comp.itertuples(index=False)):
        pct = f"{row.percent_of_reference:.2f}" if row.percent_of_reference >= .01 else "<0.01"
        ax.text(row.n_cells + class_comp.n_cells.max() * .012, yi,
                f"{row.n_cells:,} ({pct}%)", va="center", fontsize=8.5)
    ax.set_yticks(y, class_comp["class"])
    ax.invert_yaxis()
    ax.set(xlim=(0, class_comp.n_cells.max() * 1.29), xlabel="Registered reference cells (count)",
           title="Hypothalamic reference: all 26 atlas cell classes")
    ax.grid(axis="x", color="#E8ECEE", linewidth=.6)
    ax.set_axisbelow(True)
    fig.supxlabel("One male ABCA-2 donor, one registered hemisphere; descriptive composition of selected cells.",
                  fontsize=8.5, color="#465A64")
    export(fig, "cell_class_composition", "Counts and percentages of all selected HY reference cells in each of the 26 deposited atlas classes. No classes or cells are omitted.",
           ["composition_by_class.tsv"], records)
    fraction_plot(tll_class, "class", "Tll1 in the measured reference: all atlas classes",
                  "Tll1_detection_by_class", "Percentage of cells with measured Tll1 count >0 in every atlas class; annotations give exact positive/all-cell denominators. Classes follow decreasing cell count, not detection rank. No intervals or tests are calculated.",
                  "Tll1_by_class.tsv", records,
                  "All classes shown; small denominators are descriptive. A measured zero does not establish biological absence.")
    eligible = tll_structure.loc[tll_structure.n_cells.ge(100)].copy()
    display_nuclei = eligible.sort_values(["detection_fraction", "positive_cells", "substructure"],
                                          ascending=[False, False, True]).head(15)
    display_nuclei.to_csv(OUT / "Tll1_displayed_substructures.tsv", sep="\t", index=False, float_format="%.12g")
    fraction_plot(display_nuclei, "substructure", "Tll1: highest reference detection fractions among areas with ≥100 cells",
                  "Tll1_detection_by_substructure", "The 15 highest Tll1 detection fractions among deposited HY substructures with at least 100 selected cells. The complete table retains all 49 substructures, including unassigned labels and those below the display threshold. This ranking is descriptive, without exposure data or inference.",
                  "Tll1_displayed_substructures.tsv", records,
                  f"Display: top 15 of {len(eligible)} areas with ≥100 cells; full tables include all 49 areas. Positive / all cells shown.")

    total_by_class = by_class.groupby("gene_id").agg(total=("sum_imaging_counts", "sum"), positive=("positive_cells", "sum"))
    total_by_structure = by_structure.groupby("gene_id").agg(total=("sum_imaging_counts", "sum"), positive=("positive_cells", "sum"))
    expected = gene_summary.set_index("gene_id").loc[total_by_class.index]
    checks = {"exact_unique_cell_ID_join": True, "exact_gene_ID_order": True,
              "finite_nonnegative_integer_imaging_counts": True,
              "all_cells_deposited_HY": True,
              "class_composition_sums_to_37811": int(class_comp.n_cells.sum()) == len(meta),
              "substructure_composition_sums_to_37811": int(structure_comp.n_cells.sum()) == len(meta),
              "all_1122_genes_in_each_class": bool(by_class.groupby("class").gene_id.nunique().eq(1122).all()),
              "all_1122_genes_in_each_substructure": bool(by_structure.groupby("substructure").gene_id.nunique().eq(1122).all()),
              "grouped_sums_and_detections_reconcile": bool(np.allclose(total_by_class.total, expected.sum_imaging_counts) and
                  np.array_equal(total_by_class.positive, expected.positive_cells) and
                  np.allclose(total_by_structure.loc[expected.index, "total"], expected.sum_imaging_counts) and
                  np.array_equal(total_by_structure.loc[expected.index, "positive"], expected.positive_cells)),
              "Tll1_positive_cells_reconcile_to_667": int(tll_class.positive_cells.sum()) == int(tll_structure.positive_cells.sum()) == int(tll.positive_cells) == 667,
              "no_cell_level_tests_or_exposure_effects": True}
    assert all(checks.values())
    write_json(OUT / "summary.json", {"created_at_UTC": datetime.now(timezone.utc).isoformat(),
        "reference_dataset": "Zhuang-ABCA-2", "reference_donors": 1, "sex": "male",
        "registered_hemispheres": 1, "cells": len(meta), "measured_genes": len(genes),
        "cell_classes": len(class_comp), "substructures": len(structure_comp),
        "all_gene_class_rows": len(by_class), "all_gene_substructure_rows": len(by_structure),
        "Tll1": {"gene_id": str(tll.gene_id), "positive_cells": int(tll.positive_cells),
            "detection_fraction": float(tll.detection_fraction), "sum_imaging_counts": float(tll.sum_imaging_counts),
            "mean_imaging_count": float(tll.mean_imaging_count)},
        "substructures_with_at_least_100_cells": len(eligible), "substructures_displayed": len(display_nuclei),
        "HY_unassigned_cells": int(structure_comp.loc[structure_comp.substructure.eq("HY-unassigned"), "n_cells"].item()),
        "ZI_unassigned_cells": int(structure_comp.loc[structure_comp.substructure.eq("ZI-unassigned"), "n_cells"].item()),
        "statistical_scope": "Descriptive measurements from one reference donor; no per-cell differential tests, confidence intervals, clustering or exposure-response inference."})
    write_json(OUT / "methods.json", {"created_at_UTC": datetime.now(timezone.utc).isoformat(),
        "source": "Zhang et al. 2023, DOI 10.1038/s41586-023-06808-9; exact inputs retained from the spatial viewer",
        "input_files_SHA256": {str(p.relative_to(ROOT)): digest(p) for p in paths.values()},
        "cell_join": "Exact string cell_label equality; reorder annotation to NPZ labels; reject missing or duplicate identifiers",
        "gene_join": "Exact Ensembl ID order with unique deposited gene symbols",
        "groups": ["deposited atlas class", "deposited CCF substructure"],
        "anatomy": "Substructures include named nuclei, other anatomical areas, and unassigned source labels; they are not all distinct nuclei.",
        "mean_imaging_count": "Arithmetic sum of linear imaging counts divided by all selected cells in the group, including measured zeros",
        "positive_cells": "Number of selected cells with a measured imaging count strictly greater than zero",
        "detection_fraction": "positive_cells / n_cells, not a population prevalence estimate",
        "Tll1_substructure_display": "Top 15 detection fractions among substructures with n_cells >=100; all 49 retained in exported table",
        "inference": "None; no p values, confidence intervals, cell-level DE or exposure modeling",
        "source_limits": ["One male donor and one registered hemisphere", "Selected registered HY cells are not an unbiased whole-organ composition survey",
            "Taxonomy names are inherited labels, not a replacement for deposited anatomical locations", "Measured zeros do not establish biological absence",
            "Detection differences may reflect assay sensitivity, cell size or segmentation as well as underlying biology", "No exposed single cells or cell-type vulnerability estimates are added"],
        "software": {"Python": __import__("sys").version.split()[0], "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__}})
    write_json(OUT / "figure_captions.json", records)
    write_json(OUT / "validation_checks.json", {"all_checks_passed": all(checks.values()), "checks": checks,
        "scope": "ID joins, arithmetic conservation and complete group/gene export; descriptive reference only"})
    print(json.dumps({"cells": len(meta), "genes": len(genes), "classes": len(class_comp),
                      "substructures": len(structure_comp), "Tll1_positive": int(tll.positive_cells),
                      "all_checks_passed": all(checks.values())}))


if __name__ == "__main__":
    main()
