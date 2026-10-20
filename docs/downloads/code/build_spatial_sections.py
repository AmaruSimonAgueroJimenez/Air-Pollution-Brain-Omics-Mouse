#!/usr/bin/env python3
"""Build measured sectional views; no interpolation, relocation or exposure inference.

The existing ACA/HY viewer cell order and raw imaging counts are immutable inputs.
Exact cell IDs add the deposited section and experimental within-slice coordinates.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/spatial_sections"
META = ROOT / "metadata/spatial_sections"
WEB = ROOT / "docs/assets/spatial-sections"
ATLAS = ROOT / "data/reference_atlases/Zhuang_2023"
BIN_MM = .1
REGIONS = {
    "ACA": {"source": "spatial_cortex", "cells": "registered_cortex_cells.tsv.gz",
            "viewer": "viewer-cortex", "group": "cortex_measured_top10",
            "genes": ["Fosb", "Egr3", "S1pr1"], "name_en": "Anterior cingulate cortex",
            "name_es": "Corteza cingulada anterior", "anatomy_kind": "layer",
            "anatomy_kind_en": "ACA layer / subdivision", "anatomy_kind_es": "Capa / subdivisión de ACA"},
    "HY": {"source": "spatial_viewer", "cells": "registered_hypothalamus_cells.tsv.gz",
           "viewer": "viewer", "group": "hypothalamus_measured_top10",
           "genes": ["Tll1", "Fosl2", "Nr4a3"], "name_en": "Hypothalamus",
           "name_es": "Hipotálamo", "anatomy_kind": "nucleus_subregion",
           "anatomy_kind_en": "HY nucleus / subregion", "anatomy_kind_es": "Núcleo / subregión de HY"},
}


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(value), ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n")


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def source(path, url=None):
    return {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
            "sha256": sha(path), "url": url}


def read_viewer(path):
    return json.loads(path.read_text().removeprefix("window.ATLAS_DATA=").strip().removesuffix(";"))


def aggregate(values, codes, denominator):
    count = np.bincount(codes, minlength=len(denominator))
    if not np.array_equal(count, denominator):
        raise ValueError("Aggregation denominator does not include every measured cell")
    positive = np.bincount(codes, weights=(values > 0).astype(int), minlength=len(count)).astype(int)
    total = np.bincount(codes, weights=values, minlength=len(count))
    return {"n": count, "positive": positive, "detection_fraction": positive / count,
            "sum_count": total, "mean_count": total / count}


def write_table(path, frame):
    frame.to_csv(path, sep="\t", index=False, compression="gzip" if path.suffix == ".gz" else None)


def build_region(region, config, metadata, genes, priorities):
    directory = OUT / region
    directory.mkdir(parents=True, exist_ok=True)
    (WEB / f"genes-{region}").mkdir(parents=True, exist_ok=True)
    inputs = ROOT / "results" / config["source"]
    cells_path = inputs / config["cells"]
    cache_path = inputs / "reference_counts.npz"
    viewer_path = ROOT / "docs/assets" / config["viewer"] / "data.js"
    cells = pd.read_csv(cells_path, sep="\t", dtype={"cell_label": str})
    data = read_viewer(viewer_path)
    cache = np.load(cache_path)
    counts = cache["counts"]
    cell_ids = cells.cell_label.to_numpy(str)
    assert np.array_equal(cell_ids, cache["cell_labels"])
    assert list(cell_ids) == data["cell_ids"]
    assert np.array_equal(cache["gene_ids"], genes.gene_identifier.to_numpy(str))
    assert counts.shape == (len(cells), len(genes))
    assert np.isfinite(counts).all() and (counts >= 0).all()
    assert np.array_equal(counts, np.floor(counts))
    np.testing.assert_allclose(cells.z, data["x"], atol=.00005001, rtol=0)
    np.testing.assert_allclose(cells.x, data["y"], atol=.00005001, rtol=0)
    np.testing.assert_allclose(-cells.y, data["z"], atol=.00005001, rtol=0)
    cells = cells.rename(columns={"x": "ccf_ap", "y": "ccf_dv", "z": "ccf_ml"})
    joined = cells.merge(metadata, on="cell_label", how="left", sort=False, validate="one_to_one")
    assert np.array_equal(joined.cell_label, cell_ids)
    assert not joined[["brain_section_label", "experimental_x", "experimental_y", "experimental_z"]].isna().any().any()
    section_rows = []
    for label, frame in joined.groupby("brain_section_label", sort=False):
        section_rows.append({"label": label, "n": len(frame), "ap_median_mm": frame.ccf_ap.median(),
                             "ccf_ap_min_mm": frame.ccf_ap.min(), "ccf_ap_max_mm": frame.ccf_ap.max(),
                             "experimental_z_mm": frame.experimental_z.median(),
                             "experimental_z_min_mm": frame.experimental_z.min(),
                             "experimental_z_max_mm": frame.experimental_z.max(),
                             "global_section_order": int(frame.global_section_order.iloc[0]),
                             "source_section_index": int(frame.source_section_index.iloc[0])})
    sections = pd.DataFrame(section_rows).sort_values(["ap_median_mm", "label"]).reset_index(drop=True)
    section_labels = sections.label.to_list()
    section_codes = joined.brain_section_label.map({label: i for i, label in enumerate(section_labels)}).to_numpy(int)
    if region == "ACA":
        preferred = [f"ACA{s}{layer}" for s in ["d", "v"] for layer in ["1", "2/3", "5", "6a", "6b"]]
        anatomy_labels = [x for x in preferred if x in set(joined.substructure)]
        assert set(anatomy_labels) == set(joined.substructure)
    else:
        anatomy_labels = sorted(joined.substructure.unique())
    anatomy_codes = joined.substructure.map({label: i for i, label in enumerate(anatomy_labels)}).to_numpy(int)
    type_labels, type_codes = np.unique(joined["class"], return_inverse=True)
    joined["section_code"] = section_codes
    joined["anatomy_code"] = anatomy_codes
    joined["type_code"] = type_codes
    section_anatomy_codes, pairs = pd.factorize(pd.MultiIndex.from_arrays([section_codes, anatomy_codes]), sort=True)
    section_anatomy_n = np.bincount(section_anatomy_codes)
    section_anatomy = pd.DataFrame(list(pairs), columns=["section_code", "anatomy_code"])
    section_anatomy["n"] = section_anatomy_n
    bin_x = np.floor(joined.experimental_x.to_numpy() / BIN_MM).astype(int)
    bin_y = np.floor(joined.experimental_y.to_numpy() / BIN_MM).astype(int)
    bin_codes, bin_keys = pd.factorize(pd.MultiIndex.from_arrays([section_codes, bin_x, bin_y]), sort=True)
    bins = pd.DataFrame(list(bin_keys), columns=["section_code", "ix", "iy"])
    bins["x0_mm"] = bins.ix * BIN_MM
    bins["y0_mm"] = bins.iy * BIN_MM
    bins["n"] = np.bincount(bin_codes)
    joined["bin_code"] = bin_codes
    joined["section_anatomy_code"] = section_anatomy_codes
    write_table(directory / "cells_with_sections.tsv.gz", joined)
    write_table(directory / "sections.tsv", sections)
    write_table(directory / "occupied_bins.tsv.gz", bins)
    write_table(directory / "section_anatomy_groups.tsv", section_anatomy)
    selected = priorities.loc[priorities.selection_group.eq(config["group"])].sort_values("selection_rank").copy()
    assert len(selected) == 10
    gene_indices = []
    gene_entries = []
    gene_payloads = {}
    detection_rows = []
    for annotation in selected.to_dict("records"):
        name = annotation["gene"]
        matches = genes.index[genes.gene_symbol.eq(name)].to_numpy()
        if region == "ACA":
            matches = genes.index[genes.gene_identifier.eq(annotation["gene_id_ensembl"])].to_numpy()
        assert len(matches) == 1
        j = int(matches[0])
        gene_indices.append(j)
        gene_id = str(genes.iloc[j].gene_identifier)
        values = counts[:, j]
        item = {"id": gene_id, "symbol": name, "name": genes.iloc[j]["name"], "priority_rank": annotation["selection_rank"],
                "protein_name_en": annotation["protein_name_en"], "protein_name_es": annotation["protein_name_es"],
                "function_en": annotation["function_en"], "function_es": annotation["function_es"],
                "protein_coding": annotation["protein_coding"], "annotation_source": annotation["source_url"],
                "uniprot_url": annotation["uniprot_url"], "evidence_label_en": annotation["evidence_label_en"],
                "evidence_label_es": annotation["evidence_label_es"], "primary_padj": annotation["padj"],
                "primary_significant": annotation["primary_significant"], "ranking_rule": annotation["ranking_rule"],
                "source_gene_index": j, "n_positive": int(np.count_nonzero(values)),
                "detection_fraction": float(np.count_nonzero(values) / len(values))}
        gene_entries.append(item)
        payload = {"id": gene_id, "symbol": name, "n": len(values), "counts": values.astype(int),
                   "anatomy": aggregate(values, anatomy_codes, np.bincount(anatomy_codes)),
                   "section": aggregate(values, section_codes, sections.n.to_numpy()),
                   "section_anatomy": aggregate(values, section_anatomy_codes, section_anatomy_n),
                   "bins": aggregate(values, bin_codes, bins.n.to_numpy())}
        gene_payloads[name] = payload
        js = "window.registerSpatialGene(" + json.dumps(region) + "," + json.dumps(gene_id) + "," + json.dumps(clean(payload), ensure_ascii=False, separators=(",", ":"), allow_nan=False) + ");\n"
        (WEB / f"genes-{region}" / f"{gene_id}.js").write_text(js)
        for kind, labels in [("anatomy", anatomy_labels), ("section", section_labels),
                             ("section_anatomy", list(range(len(section_anatomy)))), ("bins", list(range(len(bins))))]:
            stats = payload[kind]
            for k, label in enumerate(labels):
                row = {"gene_id": gene_id, "gene": name, "group_kind": kind, "group_code": k, "label": label,
                       **{key: val[k] for key, val in stats.items()}}
                if kind == "section_anatomy":
                    row.update(section_anatomy.iloc[k].to_dict())
                    row["section_label"] = section_labels[int(row["section_code"])]
                    row["anatomy_label"] = anatomy_labels[int(row["anatomy_code"])]
                if kind == "bins":
                    row.update(bins.iloc[k].to_dict())
                    row["section_label"] = section_labels[int(row["section_code"])]
                detection_rows.append(row)
    stats_frame = pd.DataFrame(detection_rows)
    for kind in ["anatomy", "section", "section_anatomy", "bins"]:
        write_table(directory / f"detection_by_{kind}.tsv.gz", stats_frame.loc[stats_frame.group_kind.eq(kind)].dropna(axis=1, how="all"))
    write_table(directory / "genes.tsv", pd.DataFrame(gene_entries))
    spacing = np.diff(sections.experimental_z_mm.sort_values().to_numpy())
    coord_note_en = ("Within-section x/y are the deposited experimentally measured coordinates after tissue rotation/alignment (mm); "
                     "they are not CCF ML/DV coordinates. Experimental z is an estimated slice AP position. "
                     "The 2D y-axis is displayed downward; numeric coordinates remain unchanged. "
                     "The separate CCF AP/DV/ML coordinates are retained for 3D. No jitter, relocation, interpolation or mirroring.")
    coord_note_es = ("Las x/y dentro de cada corte son las coordenadas experimentales depositadas tras la rotación/alineación del tejido (mm); "
                     "no son coordenadas CCF ML/DV. La z experimental estima la posición AP del corte. "
                     "El eje y 2D se muestra hacia abajo; los valores numéricos no cambian. "
                     "Se conservan por separado las coordenadas CCF AP/DV/ML para 3D. Sin jitter, reubicación, interpolación ni espejo.")
    result = {"schema_version": 1, "dataset": "Zhuang-ABCA-2", "donor_sex": "male", "region": region,
              "region_label_en": config["name_en"], "region_label_es": config["name_es"],
              "cell_count": len(joined), "cell_ids": cell_ids, "section_labels": section_labels,
              "section_codes": section_codes, "sections": sections.to_dict("records"),
              "experimental_x": joined.experimental_x.to_numpy(), "experimental_y": joined.experimental_y.to_numpy(),
              "experimental_z": joined.experimental_z.to_numpy(),
              "global_section_order": joined.global_section_order.to_numpy(),
              "source_section_index": joined.source_section_index.to_numpy(),
              "ccf": {"ap": joined.ccf_ap.to_numpy(), "dv": joined.ccf_dv.to_numpy(), "ml": joined.ccf_ml.to_numpy()},
              "type_labels": type_labels, "type_codes": type_codes,
              "anatomy_kind": config["anatomy_kind"], "anatomy_kind_en": config["anatomy_kind_en"],
              "anatomy_kind_es": config["anatomy_kind_es"], "anatomy_labels": anatomy_labels, "anatomy_codes": anatomy_codes,
              "section_anatomy_groups": section_anatomy.to_dict("records"), "section_anatomy_codes": section_anatomy_codes,
              "bin_size_mm": BIN_MM, "bins": bins.to_dict("records"), "bin_codes": bin_codes,
              "bin_semantics_en": "Occupied 100 × 100 µm slice bins only; positive measured cells / all measured cells in that bin. Not cells/mm² or tissue density; absent bins are unobserved.",
              "bin_semantics_es": "Solo cuadrículas ocupadas de 100 × 100 µm por corte; células medidas positivas / todas las células medidas de la cuadrícula. No es células/mm² ni densidad tisular; las cuadrículas ausentes no se observaron.",
              "coordinate_note_en": coord_note_en, "coordinate_note_es": coord_note_es,
              "genes": gene_entries, "default_gene": config["genes"][0],
              "source_description_url": "https://github.com/AllenInstitute/abc_atlas_access/blob/main/descriptions/Zhuang-ABCA-2.md",
              "source_paper_doi": "https://doi.org/10.1038/s41586-023-06808-9",
              "downloads_path": f"../../downloads/spatial_sections/{region}/",
              "baseline_only": True, "independent_reference_animals": 1,
              "no_cell_level_inference": True, "empty_bin_value": None}
    (WEB / f"data-{region}.js").write_text(f"window.SPATIAL_SECTIONS_{region}=" + json.dumps(clean(result), ensure_ascii=False, separators=(",", ":"), allow_nan=False) + ";\n")
    summary = {"region": region, "cells": len(joined), "sections": len(sections), "anatomical_groups": len(anatomy_labels),
               "anatomy_kind": config["anatomy_kind"], "occupied_bins": len(bins), "bin_size_mm": BIN_MM,
               "genes": [x["symbol"] for x in gene_entries], "positive_cells": {x["symbol"]: x["n_positive"] for x in gene_entries},
               "selected_section_labels": section_labels,
               "median_section_spacing_estimated_AP_mm": float(np.median(spacing)),
               "min_section_spacing_estimated_AP_mm": float(np.min(spacing)),
               "max_section_spacing_estimated_AP_mm": float(np.max(spacing)),
               "coordinate_note_en": coord_note_en, "coordinate_note_es": coord_note_es,
               "cell_order_identical_to_existing_viewer": True, "raw_counts_exact_cache_match": True,
               "no_interpolation_or_spatial_inference": True, "single_reference_animal": True,
               "count_definition": "Nonnegative integer MERFISH gene-count matrix values; detection means raw count >0; measured zeros retained in every denominator."}
    example = sections.iloc[int(np.argmax(sections.n))]
    summary["example_section"] = {"label": example.label, "n": int(example.n), "CCF_AP_median_mm": float(example.ap_median_mm),
                                  "selection": "Largest number of measured regional cells; independent of gene expression"}
    dump(directory / "summary.json", summary)
    return {"data": result, "cells": joined, "counts": counts, "payloads": gene_payloads,
            "summary": summary, "config": config, "input_paths": [cells_path, cache_path, viewer_path]}


def export_figure(fig, stem, caption_en, caption_es, alt_en, alt_es, manifests):
    records = []
    for ext in ["png", "svg", "pdf"]:
        path = OUT / "figures" / f"{stem}.{ext}"
        fig.savefig(path, dpi=210, facecolor="white", transparent=False)
        records.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha(path)})
    plt.close(fig)
    manifests.append({"stem": stem, "caption_en": caption_en, "caption_es": caption_es,
                      "alt_en": alt_en, "alt_es": alt_es, "exports": records})


def figures(region, bundle, manifests):
    d, cells, config = bundle["data"], bundle["cells"], bundle["config"]
    curated = config["genes"]
    # The example section is selected solely by the number of measured regional cells, not expression.
    section = int(np.argmax([r["n"] for r in d["sections"]]))
    frame = cells.loc[cells.section_code.eq(section)]
    section_name = d["section_labels"][section]
    ap = d["sections"][section]["ap_median_mm"]
    bundle["summary"]["example_section"] = {"label": section_name, "n": len(frame), "CCF_AP_median_mm": ap,
                                            "selection": "Largest number of measured regional cells; independent of gene expression"}
    for spanish in [False, True]:
        T = lambda en, es: es if spanish else en
        suffix = "_es" if spanish else ""
        fig = plt.figure(figsize=(18.0 if region == "HY" else 15.0, 11.5 if region == "HY" else 10.5), layout="constrained")
        grid = fig.add_gridspec(2, 3, height_ratios=[1, .95])
        for k, gene in enumerate(curated):
            ax = fig.add_subplot(grid[0, k])
            payload = bundle["payloads"][gene]
            values = np.asarray(payload["counts"])[frame.index]
            positive = values > 0
            transformed = np.log2(1 + values)
            maximum = max(float(transformed.max()), .1)
            ax.scatter(frame.experimental_x[~positive], frame.experimental_y[~positive], s=3.0, c="#c6cdd0", alpha=.60, rasterized=True)
            scatter = ax.scatter(frame.experimental_x[positive], frame.experimental_y[positive], s=8.0, c=transformed[positive],
                                 cmap="viridis", norm=Normalize(0, maximum), edgecolors="none", rasterized=True)
            ax.set_aspect("equal", adjustable="box")
            ax.invert_yaxis()
            ax.set(xlabel=T("Section x (mm)", "x del corte (mm)"), ylabel=T("Section y (mm)", "y del corte (mm)"),
                   title=f"{chr(65+k)}  {gene} · {positive.sum():,}/{len(frame):,} " + T("detected", "detectadas"))
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.spines[["top", "right"]].set_visible(False)
            fig.colorbar(scatter, ax=ax, shrink=.60, pad=.02, label=T("log₂(1 + count)", "log₂(1 + conteo)"))
        ax = fig.add_subplot(grid[1, :2])
        labels = d["anatomy_labels"]
        anatomy_n = np.bincount(cells.anatomy_code)
        heat = np.array([bundle["payloads"][g]["anatomy"]["detection_fraction"] for g in curated]) * 100
        heat_maximum = float(min(100, max(10, 10 * np.ceil(heat.max() / 10))))
        image = ax.imshow(heat, aspect="auto", cmap="viridis", vmin=0, vmax=heat_maximum, interpolation="nearest")
        ax.set_yticks(range(len(curated)), curated)
        ax.set_xticks(range(len(labels)), [f"{label}\nn={n:,}" for label, n in zip(labels, anatomy_n)], rotation=65 if region == "HY" else 45, ha="right")
        ax.tick_params(axis="x", labelsize=6.7 if region == "HY" else 8)
        ax.set_title("D  " + T("Detection by " + config["anatomy_kind_en"], "Detección por " + config["anatomy_kind_es"]))
        if len(labels) <= 12:
            for i in range(3):
                for j in range(len(labels)):
                    ax.text(j, i, f"{heat[i,j]:.1f}", ha="center", va="center", fontsize=8,
                            color="white" if heat[i,j] < .55 * heat_maximum else "#172f39")
        fig.colorbar(image, ax=ax, shrink=.75, label=T("Measured cells with count >0 (%)", "Células medidas con conteo >0 (%)"))
        ax = fig.add_subplot(grid[1, 2])
        for k, gene in enumerate(curated):
            frac = np.asarray(bundle["payloads"][gene]["section"]["detection_fraction"]) * 100
            # Discrete section points only: no lines connecting unmeasured gaps.
            ax.scatter([x["experimental_z_mm"] for x in d["sections"]], frac, s=29, marker=["o", "s", "^"][k],
                       label=gene, color=["#0072B2", "#D55E00", "#009E73"][k], alpha=.88)
        ax.set(xlabel=T("Estimated slice AP (mm)", "AP estimada del corte (mm)"),
               ylabel=T("Detected measured cells (%)", "Células medidas detectadas (%)"), ylim=(0, 100),
               title="E  " + T("Measured sections; no interpolation", "Cortes medidos; sin interpolación") +
                     f"\nn={min(x['n'] for x in d['sections']):,}–{max(x['n'] for x in d['sections']):,} " + T("cells / section", "células / corte"))
        ax.spines[["top", "right"]].set_visible(False)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.legend(frameon=False, fontsize=9, loc="upper right")
        title = config["name_es"] if spanish else config["name_en"]
        fig.suptitle(title + " · " + T("measured slice and detection evidence", "evidencia medida por corte y detección") +
                     f"\n{section_name} · n={len(frame):,} · " + T("median CCF AP", "AP CCF mediana") + f"={ap:.3f} mm", fontsize=14)
        fig.get_layout_engine().set(rect=(.025, .12, .92, .87), h_pad=.07, w_pad=.05)
        fig.text(.5, .055, T("One unexposed male reference animal. Gray points = measured zeros; positive points = raw count >0. No simulated PM2.5 values.\n"
                             "All cells enter each denominator; anatomy and section summaries use every regional cell, not only the example slice.",
                             "Un único macho de referencia sin exposición. Gris = ceros medidos; puntos positivos = conteo >0. Sin valores PM2.5 simulados.\n"
                             "Todas las células entran en el denominador; los resúmenes anatómicos y por corte usan toda la región, no solo el corte de ejemplo."),
                 ha="center", va="center", fontsize=9, color="#334c55")
        cap_en = (f"Measured {region} reference: {', '.join(curated)} in {section_name}, selected solely as the slice with the most measured regional cells "
                  f"(n={len(frame):,}; median CCF AP {ap:.3f} mm). Within-slice positions use the deposited experimental x/y millimeters, with y displayed downward. "
                  "Gray points retain measured zeros; positive points are colored by log2(1+raw count), with an uncapped per-gene scale. "
                  "Lower panels show raw-count>0 / all measured cells by deposited anatomy and each measured slice, using every regional cell. "
                  "These are descriptive proportions from one unexposed donor, without cell-level tests, smoothing, intervals or exposed measurements. "
                  "ACA labels are layers/subdivisions; HY labels are nuclei/subregions. Unequal group n is shown; section n can be very small and all denominators are downloadable. "
                  f"The anatomy heatmap is scaled 0–{heat_maximum:g}% without clipping; region scales differ. Section points are not connected across gaps.")
        cap_es = (f"Referencia medida {region}: {', '.join(curated)} en {section_name}, elegido solo por contener el mayor número de células regionales medidas "
                  f"(n={len(frame):,}; AP CCF mediana {ap:.3f} mm). Las posiciones dentro del corte usan las x/y experimentales depositadas en mm, con y hacia abajo. "
                  "Gris conserva ceros medidos; los positivos se colorean como log2(1+conteo bruto), con escala propia sin límite percentil por gen. "
                  "Los paneles inferiores muestran conteo>0 / todas las células medidas por anatomía depositada y corte, usando toda la región. "
                  "Son proporciones descriptivas de un donante sin exposición; sin pruebas celulares, suavizado, intervalos ni mediciones expuestas. "
                  "ACA usa capas/subdivisiones; HY núcleos/subregiones. Se muestra el n desigual; puede ser muy pequeño por corte y se descargan todos los denominadores. "
                  f"El mapa anatómico usa una escala 0–{heat_maximum:g}% sin recortar; las escalas regionales difieren. No se conectan cortes a través de huecos.")
        export_figure(fig, f"{region}_sections_evidence{suffix}", cap_en, cap_es,
                      f"Three measured gene scatter maps in one {region} slice, anatomy detection heatmap and discrete slice detection points.",
                      f"Tres mapas medidos por gen en un corte de {region}, mapa de calor de detección anatómica y puntos de detección por corte.", manifests)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-figures", action="store_true")
    args = parser.parse_args()
    for path in [OUT, META, WEB, OUT / "figures"]:
        path.mkdir(parents=True, exist_ok=True)
    metadata_path = ATLAS / "metadata/Zhuang-ABCA-2/20241115/cell_metadata.csv"
    genes_path = ATLAS / "metadata/Zhuang-ABCA-2/20241115/gene.csv"
    priorities_path = ROOT / "results/gene_priorities/prioritized_genes.tsv"
    metadata = pd.read_csv(metadata_path, usecols=["cell_label", "brain_section_label", "x", "y", "z"], dtype={"cell_label": str})
    metadata = metadata.rename(columns={"x": "experimental_x", "y": "experimental_y", "z": "experimental_z"})
    assert metadata.cell_label.is_unique
    all_sections = metadata.groupby("brain_section_label").agg(experimental_z_mm=("experimental_z", "median"), n_metadata_cells=("cell_label", "size")).reset_index().sort_values(["experimental_z_mm", "brain_section_label"]).reset_index(drop=True)
    all_sections["global_section_order"] = np.arange(len(all_sections))
    all_sections["source_section_index"] = all_sections.brain_section_label.str.rsplit(".", n=1).str[-1].astype(int)
    metadata = metadata.merge(all_sections[["brain_section_label", "global_section_order", "source_section_index"]], on="brain_section_label", how="left", sort=False, validate="many_to_one")
    write_table(META / "all_brain_sections.tsv", all_sections)
    genes = pd.read_csv(genes_path)
    priorities = pd.read_csv(priorities_path, sep="\t")
    bundles = {region: build_region(region, config, metadata, genes, priorities) for region, config in REGIONS.items()}
    manifest = []
    if not args.skip_figures:
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
                             "pdf.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white"})
        for region, bundle in bundles.items():
            figures(region, bundle, manifest)
            dump(OUT / region / "summary.json", bundle["summary"])
        dump(OUT / "figure_manifest.json", {item["stem"]: {**item, "files": [x["path"] for x in item["exports"]]} for item in manifest})
    sources = [source(metadata_path, "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zhuang-ABCA-2/20241115/cell_metadata.csv"),
               source(genes_path, "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zhuang-ABCA-2/20241115/gene.csv"), source(priorities_path)]
    # Record the original atlas origins as well as exact local cache/row-order inputs.
    for relative in ["metadata/Zhuang-ABCA-2-CCF/20230830/ccf_coordinates.csv",
                     "metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv",
                     "expression_matrices/Zhuang-ABCA-2/20230830/Zhuang-ABCA-2-raw.h5ad"]:
        sources.append(source(ATLAS / relative, "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/" + relative))
    sources.append(source(ROOT / "metadata/spatial_cortex/source_manifest.json"))
    for bundle in bundles.values():
        sources.extend(source(p) for p in bundle["input_paths"])
    methods = {"created_at": datetime.now(timezone.utc).isoformat(), "scope": "Measured, unexposed baseline only; one male ABCA-2 donor",
               "cell_join": "Exact cell_label; one-to-one; no reordering versus existing regional viewer; no dropped regional cells",
               "gene_selection": "Existing curated measured-top10 groups from gene_priorities; no new significance tests; all raw counts retained",
               "coordinate_definition_en": bundles["ACA"]["data"]["coordinate_note_en"],
               "coordinate_definition_es": bundles["ACA"]["data"]["coordinate_note_es"],
               "detection": "Raw imaging count >0. n_positive/n_all measured regional cells in each explicitly identified group, including zeros.",
               "anatomical_strata": "ACA: all deposited ten layer/subdivision substructures; HY: deposited nuclear/subregion substructures, not cortical layers",
               "binning": {"width_mm": BIN_MM, "edges": "Global multiples of0.1mm in deposited experimental slice x/y; floor(x/0.1),floor(y/0.1), half-open intervals [lower,upper)",
                           "grouping": "Separate bins per exact physical brain_section_label; never pool across sections",
                           "denominator": "All measured regional cells in that occupied bin; empty bins omitted as unobserved, not zero", "not_area_density": True,
                           "no_smoothing_or_neighbor_graph": True},
               "section_example_selection": "Region slice with largest measured cell n, independent of gene expression; ties resolved by AP order",
               "section_spacing": "Estimated experimental slice AP positions are deposited registration/sectioning estimates; nominal~200µm does not mean isotropic 3D sampling",
               "uncertainty": "No binomial cell-level confidence intervals or animal-level inference from one donor; anatomy/section proportions are descriptive",
               "limitations": ["Cell counts are not independent animals", "Detection also reflects assay sensitivity and deposited filtering/registration coverage",
                               "ACA/HY are registered subsets of one hemisphere, not complete bilateral measurements", "These maps do not identify PM2.5 responses by cell type or location",
                               "No RNA imputation, new cells, relocation, jitter, mirrored hemisphere or values between measured sections"],
               "software": {"Python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__}}
    dump(META / "source_manifest.json", {"sources": sources, "atlas_paper_doi": "https://doi.org/10.1038/s41586-023-06808-9",
                                        "official_coordinate_description": "https://github.com/AllenInstitute/abc_atlas_access/blob/main/descriptions/Zhuang-ABCA-2.md"})
    dump(META / "methods.json", methods)
    checks = {"all14582_ACA_and37811_HY_cells_retained": all(len(bundles[r]["cells"]) == n for r, n in [("ACA", 14582), ("HY", 37811)]),
              "cell_order_exact_existing_viewers": True, "section_exact_one_to_one_cell_ID_join": True,
              "all_source_counts_nonnegative_integer": True, "no_coordinates_moved": True,
              "experimental_and_CCF_frames_separate": True, "ten_annotated_measured_genes_per_region": True,
              "all_zeros_in_detection_denominators": True, "group_denominators_partition_cells": True,
              "bins_never_cross_physical_sections": True, "empty_bins_unobserved_not_zero": True,
              "no_expression_imputation_interpolation_or_exposed_cells": True}
    for region, bundle in bundles.items():
        data = bundle["data"]
        for name, payload in bundle["payloads"].items():
            j = next(x["source_gene_index"] for x in data["genes"] if x["symbol"] == name)
            assert np.array_equal(payload["counts"], bundle["counts"][:, j])
            path = WEB / f"genes-{region}" / f"{payload['id']}.js"
            encoded = path.read_text()
            decoded = json.loads(encoded.split(",", 2)[2].strip().removesuffix(");"))
            assert np.array_equal(decoded["counts"], bundle["counts"][:, j])
            for kind in ["section", "anatomy", "section_anatomy", "bins"]:
                a = payload[kind]
                assert sum(a["n"]) == len(bundle["cells"])
                assert sum(a["positive"]) == np.count_nonzero(payload["counts"])
                assert sum(a["sum_count"]) == sum(payload["counts"])
                assert np.all((a["detection_fraction"] >= 0) & (a["detection_fraction"] <= 1))
    checks["all20_gene_payloads_roundtrip_exact"] = True
    checks["every_aggregate_conserves_cells_positive_cells_and_molecules"] = True
    dump(OUT / "validation.json", {"checks": checks, "all_checks_passed": all(checks.values()), "passed": all(checks.values()),
                                   "scope": "Numerical/source/serialization invariants; not evidence for exposure effects or unbiased tissue composition"})
    summaries = {region: x["summary"] for region, x in bundles.items()}
    dump(OUT / "summary.json", summaries)
    narrative_en = "# Measured sections and anatomical detection\n\n"
    narrative_es = "# Cortes medidos y detección anatómica\n\n"
    for region, bundle in bundles.items():
        s = bundle["summary"]
        narrative_en += (f"{REGIONS[region]['name_en']}: all {s['cells']:,} registered measured cells from {s['sections']} original physical sections, "
                         f"with {s['anatomical_groups']} deposited anatomical strata. Ten annotated genes retain exact raw counts and original coordinates. "
                         f"Median estimated slice AP separation is {s['median_section_spacing_estimated_AP_mm']*1000:.1f} µm within this regional subset.\n\n")
        narrative_es += (f"{REGIONS[region]['name_es']}: se conservan las {s['cells']:,} células medidas registradas de {s['sections']} cortes físicos originales, "
                         f"con {s['anatomical_groups']} estratos anatómicos depositados. Diez genes anotados conservan los conteos brutos exactos y las coordenadas originales. "
                         f"La separación AP estimada mediana de cortes es {s['median_section_spacing_estimated_AP_mm']*1000:.1f} µm dentro de este subconjunto regional.\n\n")
        narrative_en += "Baseline regional detection: " + "; ".join(
            f"{gene} {s['positive_cells'][gene]:,}/{s['cells']:,} ({100*s['positive_cells'][gene]/s['cells']:.2f}%)"
            for gene in REGIONS[region]["genes"]) + ". These fractions describe measured baseline detection, not the magnitude or location of a PM2.5 response.\n\n"
        narrative_es += "Detección regional basal: " + "; ".join(
            f"{gene} {s['positive_cells'][gene]:,}/{s['cells']:,} ({100*s['positive_cells'][gene]/s['cells']:.2f}%)"
            for gene in REGIONS[region]["genes"]) + ". Estas fracciones describen detección basal medida, no la magnitud ni la localización de una respuesta a PM2.5.\n\n"
    narrative_en += ("Detection is the number of cells with raw count >0 divided by all measured cells in the displayed anatomical group, slice or occupied 100×100 µm slice bin, including measured zeros. "
                     "Absent bins have no measurements and are not filled with zero. The single donor supplies descriptive proportions; cells are not independent biological replicates. "
                     "Section maps use deposited experimental x/y, distinct from nonlinear CCF registration. There is no interpolation, relocation, smoothing or inferred PM2.5 effect.\n")
    narrative_es += ("La detección es el número de células con conteo bruto >0 dividido por todas las células medidas del grupo anatómico, corte o cuadrícula ocupada de 100×100 µm, incluidos los ceros medidos. "
                     "Las cuadrículas ausentes carecen de mediciones y no se rellenan con cero. El único donante aporta proporciones descriptivas; las células no son réplicas biológicas independientes. "
                     "Los mapas por corte usan x/y experimentales depositadas, distintas del registro CCF no lineal. No hay interpolación, reubicación, suavizado ni efecto PM2.5 inferido.\n")
    (OUT / "narrative_en.md").write_text(narrative_en)
    (OUT / "narrative_es.md").write_text(narrative_es)
    print(json.dumps(clean(summaries), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
