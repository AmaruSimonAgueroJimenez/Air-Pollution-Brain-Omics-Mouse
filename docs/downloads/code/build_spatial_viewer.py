#!/usr/bin/env python3
"""Build an offline, measured-reference / regional PM2.5-scenario atlas viewer.

No single-cell exposure effect is fitted. A bulk shrunken log2 fold change is
transferred uniformly to the reference counts in registered hypothalamus cells.
"""
from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import nibabel as nib
import numpy as np
import pandas as pd
from skimage.measure import marching_cubes

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "data/reference_atlases/Zhuang_2023"
OUT = ROOT / "results/spatial_viewer"
WEB = ROOT / "docs/assets/viewer"


def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return safe(value.tolist())
    if isinstance(value, np.generic):
        return safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(safe(value), separators=(",", ":"), allow_nan=False) + "\n")


def first(pattern):
    return next(ATLAS.glob(pattern))


def surface(mask, step=6):
    reduced = mask[::step, ::step, ::step]
    vertices, faces, _, _ = marching_cubes(reduced.astype(np.uint8), .5,
                                          spacing=(.025 * step,) * 3)
    # The atlas cell coordinates use positive AP/DV/ML distances from the CCF
    # origin, not the centered RAS affine in the NIfTI header. Plot ML/AP/-DV.
    return {"x": np.round(vertices[:, 2], 4), "y": np.round(vertices[:, 0], 4),
            "z": np.round(-vertices[:, 1], 4), "i": faces[:, 0],
            "j": faces[:, 1], "k": faces[:, 2]}


def extract_reference(dataset, cells, gene_ids):
    cache = OUT / "reference_counts.npz"
    labels = cells.cell_label.to_numpy(dtype=str)
    if cache.exists():
        saved = np.load(cache)
        if np.array_equal(saved["cell_labels"], labels) and np.array_equal(saved["gene_ids"], gene_ids):
            return saved["counts"]
    source = first(f"expression_matrices/{dataset}/*/*-raw.h5ad")
    with h5py.File(source, "r") as h:
        matrix_ids = h["obs"][h["obs"].attrs["_index"]].asstr()[:]
        positions = pd.Index(matrix_ids).get_indexer(labels)
        if (positions < 0).any():
            raise ValueError("A registered atlas cell is absent from the measured matrix")
        observed_gene_ids = h["var"][h["var"].attrs["_index"]].asstr()[:]
        if not np.array_equal(gene_ids, observed_gene_ids):
            raise ValueError("Gene metadata and measured matrix orders differ")
        order = np.argsort(positions)
        ordered_positions = positions[order]
        counts = np.empty((len(labels), len(gene_ids)), dtype=np.float32)
        for start in range(0, h["X"].shape[0], 65536):
            end = min(start + 65536, h["X"].shape[0])
            lo, hi = np.searchsorted(ordered_positions, [start, end])
            if hi > lo:
                block = h["X"][start:end, :]
                counts[order[lo:hi], :] = block[ordered_positions[lo:hi] - start, :]
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("Invalid measured expression values")
    np.savez_compressed(cache, counts=counts, cell_labels=labels, gene_ids=gene_ids)
    return counts


def static_preview(cells, counts, gene_index, beta, gene):
    reference = counts[:, gene_index]
    predicted = reference * 2 ** beta
    values = [np.log2(1 + reference), np.log2(1 + predicted)]
    combined = np.concatenate(values)
    positive = combined[combined > 0]
    limit = max(float(np.quantile(positive, .99)) if len(positive) else 0, .5)
    fig = plt.figure(figsize=(11.8, 6.5))
    fig.subplots_adjust(left=.025, right=.89, bottom=.17, top=.85, wspace=.08)
    for k, (value, title) in enumerate(zip(values, ["Measured MERFISH reference", "Uniform regional PM2.5 simulation"])):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        # Preserve every measured cell while making sparse positive counts legible.
        for selected, size, opacity in [(reference == 0, 1.3, .12), (reference > 0, 4, .95)]:
            ax.scatter(cells.z[selected], cells.x[selected], -cells.y[selected], c=value[selected], s=size,
                       cmap="viridis", norm=Normalize(0, limit), rasterized=True, alpha=opacity)
        ax.set(xlabel="ML (mm)", ylabel="AP (mm)", zlabel="−DV (mm)", title=title)
        ax.tick_params(labelsize=8, pad=0)
        ax.view_init(22, -58)
        ax.set_box_aspect(np.ptp(cells[["z", "x", "y"]].to_numpy(), axis=0))
    fig.suptitle(f"{gene} · same cells and geometry · bulk PM/FA multiplier = {2 ** beta:.3f}", fontsize=12, y=.97)
    fig.text(.5, .035, "All registered HY cells shown; measured zeros are faint and positive measured counts are emphasized.",
             ha="center", fontsize=9, color="#4b6269")
    fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, limit), cmap="viridis"),
                 ax=fig.axes[:2], label="log₂(1 + imaging count); common positive-value 99th-percentile cap", shrink=.65)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"reference_vs_simulation.{ext}", dpi=220, bbox_inches="tight", pad_inches=.2)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="Zhuang-ABCA-2", choices=["Zhuang-ABCA-2"],
                        help="The validated primary reference is the registered male ABCA-2 donor")
    args = parser.parse_args()
    if args.dataset != "Zhuang-ABCA-2":
        raise SystemExit("This release fixes the male ABCA-2 donor; use a separate output for donor sensitivity")
    OUT.mkdir(parents=True, exist_ok=True)
    (WEB / "genes").mkdir(parents=True, exist_ok=True)
    lookup = pd.read_csv(first("metadata/Allen-CCF-2020/*/views/parcellation_to_parcellation_term_membership_acronym.csv"))
    ccf = pd.read_csv(first(f"metadata/{args.dataset}-CCF/*/ccf_coordinates.csv"), dtype={"cell_label": str})
    ccf = ccf.merge(lookup, on="parcellation_index", how="left", validate="many_to_one")
    cells = ccf.loc[ccf.division.eq("HY")].copy()
    annotation = pd.read_csv(first(f"metadata/{args.dataset}/*/views/cell_metadata_with_cluster_annotation.csv"),
                             usecols=["cell_label", "class", "subclass", "cluster_alias"], dtype={"cell_label": str})
    cells = cells.merge(annotation, on="cell_label", how="left", validate="one_to_one")
    if cells[["class", "subclass", "x", "y", "z"]].isna().any().any():
        raise ValueError("Incomplete coordinate / cell-type join")
    if cells.cell_label.duplicated().any():
        raise ValueError("Non-unique registered cell identifiers")
    genes = pd.read_csv(first(f"metadata/{args.dataset}/*/gene.csv"))
    counts = extract_reference(args.dataset, cells, genes.gene_identifier.to_numpy(dtype=str))
    de = pd.read_csv(ROOT / "results/hypothalamus/differential_expression.tsv", sep="\t").set_index("gene")
    normalized = pd.read_csv(ROOT / "results/hypothalamus/normalized_counts.tsv.gz", sep="\t").set_index("gene")
    samples = pd.read_csv(ROOT / "results/hypothalamus/sample_metadata.tsv", sep="\t")
    stability = pd.read_csv(ROOT / "results/hypothalamus/leave_one_out_primary_hits.tsv", sep="\t").set_index("gene")
    gene_meta = []
    for i, row in enumerate(genes.itertuples()):
        symbol = row.gene_symbol
        values = np.ascontiguousarray(counts[:, i], dtype="<f4")
        packed = base64.b64encode(gzip.compress(values.tobytes(), compresslevel=6, mtime=0)).decode()
        gene_file = WEB / "genes" / f"{row.gene_identifier}.js"
        gene_file.write_text("window.receiveAtlasGene(" + json.dumps({"id": row.gene_identifier,
                            "n": len(values), "gzip_f32le": packed}, separators=(",", ":")) + ");\n")
        effect = de.loc[symbol].to_dict() if symbol in de.index else {}
        estimable = bool(effect and all(np.isfinite(effect.get(c, np.nan))
                         for c in ("log2FoldChange_shrunk", "log2FoldChange", "lfcSE", "padj")))
        normalized_values = []
        if symbol in normalized.index:
            normalized_values = [{"sample": s.sample, "label": s.label, "condition": s.condition,
                                  "count": float(normalized.at[symbol, s.sample])} for s in samples.itertuples()]
        entry = {"id": row.gene_identifier, "symbol": symbol, "name": row.name,
                 "reference_nonzero_cells": int(np.count_nonzero(values)),
                 "reference_max": float(values.max()), "estimable": estimable,
                 "significant": bool(estimable and effect["padj"] < .05),
                 "effect": effect, "bulk_values": normalized_values,
                 "leave_one_out": stability.loc[symbol].to_dict() if symbol in stability.index else None}
        if estimable:
            entry["ci_low"] = float(effect["log2FoldChange"] - 1.96 * effect["lfcSE"])
            entry["ci_high"] = float(effect["log2FoldChange"] + 1.96 * effect["lfcSE"])
        gene_meta.append(safe(entry))
    significant = sorted([g for g in gene_meta if g["significant"]], key=lambda g: g["effect"]["padj"])
    default_gene = significant[0]["symbol"] if significant else "Gfap"
    types, type_codes = np.unique(cells["class"], return_inverse=True)
    regions, region_codes = np.unique(cells["substructure"], return_inverse=True)
    annotation_image = nib.load(first("image_volumes/Allen-CCF-2020/*/annotation_25.nii.gz"))
    volume = np.asarray(annotation_image.dataobj)
    hy_labels = lookup.loc[lookup.division.eq("HY"), "parcellation_index"].to_numpy()
    geometry = {"brain": surface(volume > 0), "hypothalamus": surface(np.isin(volume, hy_labels), step=3)}
    diagnostic = ccf.loc[ccf.parcellation_index.gt(0)].sample(min(10000, len(ccf)), random_state=20261018)
    voxels = np.floor(diagnostic[["x", "y", "z"]].to_numpy() / .025).astype(int)
    inside = ((voxels >= 0) & (voxels < volume.shape)).all(axis=1)
    v = voxels[inside]
    agreement = np.mean(volume[v[:, 0], v[:, 1], v[:, 2]] == diagnostic.loc[inside, "parcellation_index"].to_numpy())
    if agreement < .85:
        raise ValueError("Atlas coordinates and anatomical volume do not share the expected frame")
    data = {"dataset": args.dataset, "donor_sex": "male", "default_gene": default_gene,
            "cell_count": len(cells), "genes": gene_meta, "type_labels": types,
            "region_labels": regions, "type_codes": type_codes, "region_codes": region_codes,
            "cell_ids": cells.cell_label.to_numpy(), "x": np.round(cells.z.to_numpy(), 4),
            "y": np.round(cells.x.to_numpy(), 4), "z": np.round(-cells.y.to_numpy(), 4),
            "meshes": geometry}
    (WEB / "data.js").write_text("window.ATLAS_DATA=" + json.dumps(safe(data), separators=(",", ":"), allow_nan=False) + ";\n")
    import plotly
    shutil.copy(Path(plotly.__file__).parent / "package_data/plotly.min.js", WEB / "plotly.min.js")
    cells.to_csv(OUT / "registered_hypothalamus_cells.tsv.gz", sep="\t", index=False)
    pd.DataFrame([{k: g[k] for k in ["id", "symbol", "reference_nonzero_cells", "estimable", "significant"]}
                  for g in gene_meta]).to_csv(OUT / "gene_panel_coverage.tsv", sep="\t", index=False)
    default = next(g for g in gene_meta if g["symbol"] == default_gene)
    default_index = int(np.flatnonzero(genes.gene_symbol.eq(default_gene))[0])
    static_preview(cells, counts, default_index, default["effect"]["log2FoldChange_shrunk"], default_gene)
    summary = {"created_at": datetime.now(timezone.utc).isoformat(), "reference_dataset": args.dataset,
               "reference_sex": "male", "registered_reference_cells": len(ccf),
               "hypothalamus_cells": len(cells), "cell_subsampling": False,
               "cell_classes": len(types), "anatomical_substructures": len(regions),
               "measured_genes": len(genes), "genes_with_estimable_bulk_effect": sum(g["estimable"] for g in gene_meta),
               "primary_fdr_hits_available_in_measured_panel": [g["symbol"] for g in significant],
               "default_gene": default_gene, "default_multiplier": 2 ** default["effect"]["log2FoldChange_shrunk"],
               "anatomical_coordinate_check_n": int(inside.sum()), "coordinate_parcellation_agreement": float(agreement),
               "coordinate_frame": "CCF AP,DV,ML in mm; plot ML,AP,-DV; NIfTI centered RAS affine is not applied to CCF cell distances",
               "whole_brain_surface_step_mm": .15, "hypothalamus_surface_step_mm": .075,
               "model": "E_scenario(c,g)=E_reference(c,g)*2**shrunken_bulk_log2FC(g), only in HY",
               "color_transform": "log2(1+value), common scale across both panels, optional common 99th percentile of positive reference/scenario values",
               "uncertainty": "Wald MLE coefficient +/-1.96 SE defines alternative group-level scenarios, not single-cell confidence intervals",
               "limitations": ["One registered hemisphere; no mirroring or interpolation of expression", "Uniform regional multiplier; no cell-type-specific exposure effect is identified", "Atlas is not an experimental control for the seven bulk animals", "Cell positions and abundance remain unchanged", "Significance is sensitive to individual animals; systematic leave-one-out results are retained"]}
    dump(OUT / "summary.json", summary)
    dump(OUT / "validation_checks.json", {"complete_hypothalamus_coordinate_join": True,
          "unique_cell_identifiers": True, "exact_cell_to_matrix_match": True, "same_gene_order": True,
          "finite_nonnegative_reference_counts": True, "no_cell_subsampling_or_mirroring": True,
          "common_reference_and_scenario_geometry": True, "anatomical_coordinate_check_pass": bool(agreement >= .85)})
    print(json.dumps(safe(summary), indent=2))


if __name__ == "__main__":
    main()
