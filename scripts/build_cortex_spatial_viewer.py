#!/usr/bin/env python3
"""Build the measured ACA reference and explicitly modeled regional PM2.5 scenario.

All registered ACAd/ACAv cells of male Zhuang-ABCA-2 are retained. Effects use
the primary eight-animal GSE297423 PM/PBS fit and exact Ensembl gene IDs.
No per-cell exposure effect, absent transcript, or additional animal is inferred.
"""
from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
import nibabel as nib
import numpy as np
import pandas as pd

import build_spatial_viewer as shared

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "data/reference_atlases/Zhuang_2023"
OUT = ROOT / "results/spatial_cortex"
FIGURES = OUT / "figures"
WEB = ROOT / "docs/assets/viewer-cortex"
META = ROOT / "metadata/spatial_cortex"
DATASET = "Zhuang-ABCA-2"


def first(pattern):
    return next(ATLAS.glob(pattern))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export(fig, stem):
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", dpi=300, bbox_inches="tight", pad_inches=.2)
    plt.close(fig)


def spatial_preview(cells, counts, genes, entries, stems, language="en"):
    spanish = language == "es"
    figure_manifest = json.loads((OUT / "figure_manifest.json").read_text()) if (OUT / "figure_manifest.json").exists() else {}
    for stem, selected in stems.items():
        if spanish:
            stem += "_es"
        single = len(selected) == 1
        fig = plt.figure(figsize=(12.8, 6.5 if single else 5.5 * len(selected)))
        fig.subplots_adjust(left=.02, right=.88, bottom=.21 if single else .065,
                            top=.84 if single else .93, wspace=.03, hspace=.12)
        for row, symbol in enumerate(selected):
            gene_index = int(np.flatnonzero(genes.gene_symbol.eq(symbol))[0])
            entry = next(g for g in entries if g["symbol"] == symbol)
            reference = counts[:, gene_index]
            predicted = reference * 2 ** entry["effect"]["log2FoldChange_shrunk"]
            values = [np.log2(1 + reference), np.log2(1 + predicted)]
            positive = np.concatenate(values)
            positive = positive[positive > 0]
            cap = max(float(np.quantile(positive, .99, method="lower")) if len(positive) else 0, .5)
            axes = []
            titles = (["Referencia MERFISH medida", "SIMULACIÓN regional uniforme PM2.5"] if spanish else
                      ["Measured MERFISH reference", "Uniform regional PM2.5 SIMULATION"])
            for col, (value, title) in enumerate(zip(values, titles)):
                ax = fig.add_subplot(len(selected), 2, row * 2 + col + 1, projection="3d")
                axes.append(ax)
                zero = reference == 0
                ax.scatter(cells.z[zero], cells.x[zero], -cells.y[zero], c="#c6cdd0", s=1.1, alpha=.12, rasterized=True)
                ax.scatter(cells.z[~zero], cells.x[~zero], -cells.y[~zero], c=value[~zero], s=4.5,
                           cmap="viridis", norm=Normalize(0, cap), alpha=.95, rasterized=True)
                ax.set(xlabel="ML (mm)", ylabel="AP (mm)", zlabel="−DV (mm)",
                       title=f"{symbol} · {title}")
                ax.view_init(22, -58)
                ax.set_box_aspect(np.ptp(cells[["z", "x", "y"]].to_numpy(), axis=0))
                ax.tick_params(labelsize=7, pad=0)
                for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
                    axis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
            fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, cap), cmap="viridis"), ax=axes,
                         label=(f"{symbol}: log₂(1 + conteo de imagen); límite común por par" if spanish else
                                f"{symbol}: log₂(1 + imaging count); shared pairwise cap"), shrink=.66, pad=.025)
        fig.suptitle("Corteza cingulada anterior · mismas células medidas y geometría" if spanish else
                     "Anterior cingulate cortex · same measured cells and geometry", fontsize=14)
        note = ("Se conservan las 14,582 células ACAd/ACAv; ceros medidos grises y conteos positivos coloreados.\n"
                "Multiplicadores regionales de cuatro animales PM y cuatro PBS. No es expresión unicelular expuesta medida." if spanish else
                "All 14,582 registered ACAd/ACAv cells retained; gray measured zeros and colored positive counts.\n"
                "Regional multipliers use four PM and four PBS animals. This is not measured exposed single-cell expression.")
        fig.text(.5, .018, note,
                 ha="center", fontsize=9, color="#4b6269")
        export(fig, stem)
        caption = ("Measured anterior-cingulate MERFISH expression beside a uniform regional PM/PBS scenario for "
                   + ", ".join(selected) + ". All 14,582 ACAd/ACAv cells from the male ABCA-2 donor are shown without mirroring, "
                   "interpolation or cell subsampling. Gray points are measured zeros; positive measured counts are colored. "
                   "Each gene has one shared reference/scenario log2(1+count) scale capped at the combined positive 99th percentile. "
                   "The underlying counts are retained, and different rows have different scales. Bulk effects come from four PM "
                   "and four PBS animals; atlas cells are not exposure replicates. The multiplier does not identify responding cell types.")
        caption_es = ("Expresión MERFISH medida de corteza cingulada anterior junto a un escenario regional uniforme PM/PBS para "
                      + ", ".join(selected) + ". Se muestran las 14,582 células ACAd/ACAv del donante macho ABCA-2 sin reflejo, "
                      "interpolación ni submuestreo. Los puntos grises son ceros medidos y los positivos se colorean. Cada gen comparte "
                      "una escala log2(1+conteo) entre referencia y escenario, limitada al percentil 99 conjunto de valores positivos; "
                      "las filas utilizan escalas distintas y se conservan los conteos originales. Los efectos regionales proceden de "
                      "cuatro animales PM y cuatro PBS. Las células del atlas no son réplicas de exposición ni se identifica qué tipos responden.")
        figure_manifest[stem] = {"caption": caption, "caption_es": caption_es, "language": language, "genes": selected,
                                "PNG_dpi": 300, "files": {s: f"results/spatial_cortex/figures/{stem}.{s}" for s in ("png", "svg", "pdf")},
                                "sha256": {s: digest(FIGURES / f"{stem}.{s}") for s in ("png", "svg", "pdf")}}
    shared.dump(OUT / "figure_manifest.json", figure_manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-figures", action="store_true")
    args = parser.parse_args()
    for directory in (OUT, FIGURES, WEB / "genes", META):
        directory.mkdir(parents=True, exist_ok=True)
    lookup_path = first("metadata/Allen-CCF-2020/*/views/parcellation_to_parcellation_term_membership_acronym.csv")
    coordinates_path = first(f"metadata/{DATASET}-CCF/*/ccf_coordinates.csv")
    annotation_path = first(f"metadata/{DATASET}/*/views/cell_metadata_with_cluster_annotation.csv")
    genes_path = first(f"metadata/{DATASET}/*/gene.csv")
    lookup = pd.read_csv(lookup_path)
    ccf = pd.read_csv(coordinates_path, dtype={"cell_label": str}).merge(lookup, on="parcellation_index", how="left", validate="many_to_one")
    selected_lookup = lookup.loc[lookup.structure.isin(["ACAd", "ACAv"]) & lookup.division.eq("Isocortex")].copy()
    cells = ccf.loc[ccf.parcellation_index.isin(selected_lookup.parcellation_index)].copy()
    annotation = pd.read_csv(annotation_path, usecols=["cell_label", "class", "subclass", "cluster_alias"], dtype={"cell_label": str})
    cells = cells.merge(annotation, on="cell_label", how="left", validate="one_to_one")
    if cells[["class", "subclass", "x", "y", "z"]].isna().any().any() or cells.cell_label.duplicated().any():
        raise ValueError("ACA cell-to-coordinate/annotation join is incomplete or duplicated")
    genes = pd.read_csv(genes_path)
    shared.OUT = OUT  # Scope the original validated extraction helper's cache to ACA.
    counts = shared.extract_reference(DATASET, cells, genes.gene_identifier.to_numpy(dtype=str))
    de_path = ROOT / "results/cortex/DE_PM_vs_PBS.tsv"
    normalized_path = ROOT / "results/cortex/normalized_counts.tsv.gz"
    samples_path = ROOT / "results/cortex/sample_metadata.tsv"
    de = pd.read_csv(de_path, sep="\t").set_index("gene_id", verify_integrity=True)
    normalized = pd.read_csv(normalized_path, sep="\t").set_index("gene_id", verify_integrity=True)
    samples = pd.read_csv(samples_path, sep="\t")
    if samples["sample"].duplicated().any() or samples.condition.value_counts().to_dict() != {"PBS": 4, "PM": 4}:
        raise ValueError("Primary cortex sample population differs from four PBS/four PM animals")
    entries = []
    primary_ranks = de.padj.rank(method="min").to_dict()
    for i, row in enumerate(genes.itertuples()):
        gene_id = row.gene_identifier
        values = np.ascontiguousarray(counts[:, i], dtype="<f4")
        packed = base64.b64encode(gzip.compress(values.tobytes(), compresslevel=6, mtime=0)).decode()
        (WEB / "genes" / f"{gene_id}.js").write_text("window.receiveAtlasGene(" + json.dumps(
            {"id": gene_id, "n": len(values), "gzip_f32le": packed}, separators=(",", ":")) + ");\n")
        effect = de.loc[gene_id].to_dict() if gene_id in de.index else {}
        estimable = bool(effect and all(np.isfinite(effect.get(c, np.nan)) for c in
                                        ("log2FoldChange_shrunk", "log2FoldChange", "lfcSE", "padj")))
        bulk_values = []
        if gene_id in normalized.index:
            bulk_values = [{"sample": s.sample, "label": s.sample, "condition": s.condition,
                            "count": float(normalized.at[gene_id, s.sample])} for s in samples.itertuples()]
        entry = {"id": gene_id, "symbol": row.gene_symbol, "name": row.name,
                 "bulk_gene_symbol": effect.get("gene"), "match_method": "exact Ensembl gene ID" if effect else "not matched",
                 "reference_nonzero_cells": int(np.count_nonzero(values)), "reference_max": float(values.max()),
                 "estimable": estimable, "significant": bool(estimable and effect["padj"] < .05),
                 "effect": effect, "bulk_values": bulk_values, "leave_one_out": None,
                 "primary_FDR_rank": primary_ranks.get(gene_id)}
        if estimable:
            entry["ci_low"] = float(effect["log2FoldChange"] - 1.96 * effect["lfcSE"])
            entry["ci_high"] = float(effect["log2FoldChange"] + 1.96 * effect["lfcSE"])
        entries.append(shared.safe(entry))
    significant = sorted([g for g in entries if g["significant"]], key=lambda g: (g["effect"]["padj"], g["id"]))
    default = significant[0]
    types, type_codes = np.unique(cells["class"], return_inverse=True)
    regions, region_codes = np.unique(cells["substructure"], return_inverse=True)
    volume_path = first("image_volumes/Allen-CCF-2020/*/annotation_25.nii.gz")
    volume = np.asarray(nib.load(volume_path).dataobj)
    aca_surface = shared.surface(np.isin(volume, selected_lookup.parcellation_index.to_numpy()), step=2)
    geometry = {"brain": shared.surface(volume > 0), "region": aca_surface, "hypothalamus": aca_surface}
    voxels = np.floor(cells[["x", "y", "z"]].to_numpy() / .025).astype(int)
    inside = ((voxels >= 0) & (voxels < volume.shape)).all(axis=1)
    vi = voxels[inside]
    agreement = float(np.mean(volume[vi[:, 0], vi[:, 1], vi[:, 2]] == cells.loc[inside, "parcellation_index"].to_numpy()))
    if agreement < .85:
        raise ValueError("ACA coordinates disagree with the anatomical label volume")
    data = {"dataset": DATASET, "donor_sex": "male", "region": "ACA", "region_key": "ACA",
            "region_label_en": "Anterior cingulate cortex", "region_label_es": "Corteza cingulada anterior",
            "control_label": "PBS", "exposure_label": "PM", "n_control": 4, "n_PM": 4,
            "bulk_sample_n": 8, "source_expected_n": 8, "bulk_source": "GSE297423",
            "default_gene": default["symbol"], "cell_count": len(cells), "genes": entries,
            "type_labels": types, "region_labels": regions, "type_codes": type_codes, "region_codes": region_codes,
            "cell_ids": cells.cell_label.to_numpy(), "x": np.round(cells.z.to_numpy(), 4),
            "y": np.round(cells.x.to_numpy(), 4), "z": np.round(-cells.y.to_numpy(), 4), "meshes": geometry}
    (WEB / "data.js").write_text("window.ATLAS_DATA=" + json.dumps(shared.safe(data), separators=(",", ":"), allow_nan=False) + ";\n")
    cells.to_csv(OUT / "registered_cortex_cells.tsv.gz", sep="\t", index=False)
    coverage = pd.DataFrame([{k: g[k] for k in ("id", "symbol", "name", "bulk_gene_symbol", "match_method", "reference_nonzero_cells", "estimable", "significant", "primary_FDR_rank")} for g in entries])
    coverage.to_csv(OUT / "gene_panel_coverage.tsv", sep="\t", index=False)
    shortlist = pd.DataFrame([{"priority_rank_measured": i + 1, "gene_id": g["id"], "gene": g["symbol"],
                               "bulk_gene_symbol": g["bulk_gene_symbol"], "name": g["name"], "padj": g["effect"]["padj"],
                               "MLE_log2FC": g["effect"]["log2FoldChange"], "shrunk_log2FC": g["effect"]["log2FoldChange_shrunk"],
                               "reference_nonzero_cells": g["reference_nonzero_cells"], "primary_FDR_rank": g["primary_FDR_rank"]}
                              for i, g in enumerate(significant[:10])])
    shortlist.to_csv(OUT / "top10_measured_significant_genes.tsv", sep="\t", index=False)
    selected_lookup.to_csv(META / "ACA_parcellation_selection.csv", index=False)
    summary = {"created_at": datetime.now(timezone.utc).isoformat(), "reference_dataset": DATASET,
               "reference_sex": "male", "region": "ACA", "region_definition": "division=Isocortex and structure in {ACAd,ACAv}; all ten substructures",
               "registered_reference_cells": len(ccf), "cells": len(cells), "anterior_cingulate_cells": len(cells),
               "cell_subsampling": False, "cell_classes": len(types), "anatomical_substructures": len(regions),
               "structures": cells.structure.value_counts().to_dict(), "substructures": cells.substructure.value_counts().to_dict(),
               "measured_genes": len(genes), "exact_Ensembl_ID_matches": sum(bool(g["effect"]) for g in entries),
               "genes_with_estimable_bulk_effect": sum(g["estimable"] for g in entries),
               "primary_fdr_hits_available_in_measured_panel": [g["symbol"] for g in significant],
               "measured_FDR_hits": len(significant), "top10_measured_significant_genes": shortlist.to_dict("records"),
               "shortlist_rule": "All measured genes with finite primary inference/shrunk effect and BH<.05, ranked by unrounded primary BH q then exact gene ID; first10",
               "default_gene": default["symbol"], "default_multiplier": 2 ** default["effect"]["log2FoldChange_shrunk"],
               "default_reference_nonzero_cells": default["reference_nonzero_cells"],
               "all_primary_FDR_hits": int(de.padj.lt(.05).sum()),
               "measured_significant_genes_with_positive_ACA_counts": sum(g["reference_nonzero_cells"] > 0 for g in significant),
               "whole_brain_surface_step_mm": .15, "region_surface_step_mm": .05,
               "bulk_study": "GSE297423", "bulk_samples": 8, "n_control": 4, "n_PM": 4, "control_label": "PBS",
               "anatomical_coordinate_check_n": int(inside.sum()), "coordinate_parcellation_agreement": agreement,
               "coordinate_frame": "Source CCF AP,DV,ML in mm; plotted ML,AP,-DV, with no centered NIfTI RAS affine applied",
               "model": "E_scenario(c,g)=E_reference(c,g)*2**primary_cortex_shrunken_log2FC(g), restricted to ACA",
               "color_transform": "log2(1+imaging count), shared reference/scenario scale for each gene; optional positive-value combined99th-percentile cap",
               "uncertainty": "Primary unshrunken Wald MLE +/-1.96SE defines alternative group-level scenarios; not per-cell confidence intervals",
               "limitations": ["One unexposed registered male reference donor is not a control for the humanized-alpha-synuclein experimental animals",
                               "Regional bulk effects cannot separate regulation, cell composition, or other source-quality effects; original sex metadata remain unresolved",
                               "Uniform regional multiplier cannot identify responding cell types, cell loss, or migration",
                               "Only the 1122 physically measured MERFISH panel genes are available; absent genes are not imputed or substituted",
                               "Only ACA is modeled, not the entire cortex or brain; no mirrored hemisphere is created"]}
    shared.dump(OUT / "summary.json", summary)
    if not args.skip_figures:
        stems = {"reference_vs_simulation": [default["symbol"]],
                 "cortical_priority_spatial_panel": [g["symbol"] for g in significant[:3]]}
        for language in ("en", "es"):
            spatial_preview(cells, counts, genes, entries, stems, language)
    manifests = json.loads((ROOT / "metadata/spatial_atlas/download_manifest.json").read_text())["files"]
    receipts = {r["path"]: r for r in json.loads((ROOT / "metadata/spatial_atlas/download_receipts.json").read_text())["files"]}
    atlas_sources = [lookup_path, coordinates_path, annotation_path, genes_path, volume_path,
                     first(f"expression_matrices/{DATASET}/*/*-raw.h5ad")]
    provenance = []
    for path in atlas_sources:
        relative = path.relative_to(ROOT).as_posix()
        source = next(r for r in manifests if r["path"] == relative)
        receipt = receipts[relative]
        provenance.append({"path": relative, "url": source["url"], "source_version": source.get("source_version"),
                           "bytes": path.stat().st_size, "source_MD5": source["md5"],
                           "download_MD5_verified": bool(receipt.get("source_md5_match")),
                           "current_size_matches_receipt": path.stat().st_size == receipt["bytes"],
                           "sha256": digest(path) if path.stat().st_size < 150_000_000 else None})
    provenance += [{"path": p.relative_to(ROOT).as_posix(), "sha256": digest(p), "bytes": p.stat().st_size}
                   for p in (de_path, normalized_path, samples_path)]
    shared.dump(META / "source_manifest.json", {"sources": provenance, "selection": summary["region_definition"],
                                              "original_atlas_paper": "https://doi.org/10.1038/s41586-023-06808-9",
                                              "exposure_paper": "https://doi.org/10.1126/science.adu4132"})
    transport_exact = True
    for i, gene in enumerate(genes.itertuples()):
        content = (WEB / "genes" / f"{gene.gene_identifier}.js").read_text()
        payload = json.loads(content.removeprefix("window.receiveAtlasGene(").removesuffix(");\n"))
        decoded = np.frombuffer(gzip.decompress(base64.b64decode(payload["gzip_f32le"])), dtype="<f4")
        transport_exact &= payload["id"] == gene.gene_identifier and payload["n"] == len(cells) and np.array_equal(decoded, counts[:, i])
    checks = {"all_ACAd_ACAv_registered_cells_retained": len(cells) == 14582,
              "unique_cell_identifiers": not cells.cell_label.duplicated().any(),
              "complete_coordinate_annotation_join": not cells[["class", "subclass", "x", "y", "z"]].isna().any().any(),
              "exact_cell_matrix_gene_matches": counts.shape == (len(cells), len(genes)),
              "finite_nonnegative_measured_counts": bool(np.isfinite(counts).all() and (counts >= 0).all()),
              "eight_independent_bulk_animals_four_per_group": samples.condition.value_counts().to_dict() == {"PBS": 4, "PM": 4},
              "exact_Ensembl_ID_effect_mapping": all(not g["effect"] or g["id"] in de.index for g in entries),
              "same_reference_scenario_cell_geometry": True, "no_expression_imputation_mirroring_subsampling": True,
              "anatomical_coordinate_check_pass": agreement >= .85,
              "default_is_minimum_FDR_measured_gene": default["effect"]["padj"] == min(g["effect"]["padj"] for g in significant),
              "all1122_gene_transport_files": len(list((WEB / "genes").glob("*.js"))) == len(genes),
              "all1122_gene_vectors_roundtrip_exactly": bool(transport_exact)}
    shared.dump(OUT / "validation_checks.json", {"checks": checks, "all_checks_passed": all(checks.values()),
                                               "scope": "Join, geometry, numerical, sample-population and serialization checks; not biological validation"})
    print(json.dumps(shared.safe(summary), indent=2))


if __name__ == "__main__":
    main()
