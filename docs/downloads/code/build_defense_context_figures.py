#!/usr/bin/env python3
"""Slide-sized views of saved project results; no model is fitted here."""
from pathlib import Path
import hashlib
import json
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/defense_slides"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
INK, TEAL, RUST, GOLD = "#123d46", "#376e75", "#b65e3c", "#987126"
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":18,
                     "axes.labelsize":20, "axes.titlesize":22,
                     "xtick.labelsize":17, "ytick.labelsize":17,
                     "legend.fontsize":16, "axes.spines.top":False,
                     "axes.spines.right":False, "svg.fonttype":"none",
                     "pdf.fonttype":42, "axes.labelcolor":INK,
                     "text.color":INK, "axes.edgecolor":INK})
MANIFEST = {}

def export(fig, stem, lang, caption, sources):
    name = stem + ("_es" if lang == "es" else "")
    files = []
    for ext in ("png", "svg", "pdf"):
        p = FIG / f"{name}.{ext}"
        fig.savefig(p, dpi=180, facecolor="white")
        files.append({"path":str(p.relative_to(ROOT)), "sha256":hashlib.sha256(p.read_bytes()).hexdigest()})
    plt.close(fig)
    MANIFEST[name] = {"language":lang, "caption":caption,
                      "source_files":[{"path":s,"sha256":hashlib.sha256((ROOT/s).read_bytes()).hexdigest()} for s in sources],
                      "files":files,"no_model_refit":True}

def workflow(stem, lang, caption):
    """Two evidence branches, with an explicit assumption-based bulk transfer."""
    en = lang == "en"
    fig, ax = plt.subplots(figsize=(15, 6.8))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    fig.subplots_adjust(left=.015, right=.985, bottom=.035, top=.98)

    def txt(x, y, value, size=17, weight="normal", color=INK, **kwargs):
        return ax.text(x, y, value, fontsize=size, fontweight=weight,
                       color=color, va="top", linespacing=1.18, **kwargs)

    def arrow(points, color=TEAL, dashed=False):
        # Orthogonal segments keep the two inputs to the scenario explicit.
        style = (0, (4, 3)) if dashed else "solid"
        if len(points) > 2:
            ax.plot(*zip(*points[:-1]), color=color, lw=1.7,
                    linestyle=style, solid_capstyle="butt", zorder=1)
        ax.add_patch(FancyArrowPatch(points[-2], points[-1],
                    arrowstyle="-|>", mutation_scale=15, lw=1.7,
                    color=color, linestyle=style, zorder=1))

    txt(.02, .975, "EXPERIMENTAL BRAIN STUDIES" if en else "ESTUDIOS EXPERIMENTALES DE CEREBRO",
        14, "bold", TEAL)
    ax.plot([.02, .985], [.535, .535], color="#cfdcdd", lw=1.1)
    txt(.02, .50, "INDEPENDENT UNEXPOSED REFERENCE" if en else "REFERENCIA INDEPENDIENTE SIN EXPOSICIÓN",
        14, "bold", TEAL)

    # Experimental source and its two analysis branches.
    ax.plot([.023, .023], [.61, .925], color=TEAL, lw=3)
    txt(.04, .91, "Four source studies" if en else "Cuatro estudios fuente", 20, "bold")
    txt(.04, .835,
        "HY and ACA bulk RNA-seq\nHippocampal miRNA array\nHippocampal O₃ RNA-seq\n24 paired-end O₃ libraries" if en else
        "RNA-seq bulk de HY y ACA\nArray miRNA de hipocampo\nRNA-seq hipocampal de O₃\n24 bibliotecas O₃ paired-end",
        16.5)
    txt(.44, .915, "Original-study reproduction" if en else "Reproducción de los estudios", 20, "bold")
    txt(.44, .845,
        "Author code, thresholds and source animals\nExact recovery or documented discrepancies" if en else
        "Código, umbrales y animales de los autores\nRecuperación exacta o diferencias documentadas",
        17)
    txt(.44, .705, "New study-specific analyses" if en else "Análisis nuevos por estudio", 20, "bold")
    txt(.44, .635,
        "RNA/miRNA models, GO and interactions\nQC and animal-omission sensitivity" if en else
        "Modelos RNA/miRNA, GO e interacciones\nControl de calidad y omisión de animales",
        17)
    arrow([(.305, .775), (.36, .775)])
    ax.plot([.36, .36], [.675, .875], color=TEAL, lw=1.7)
    arrow([(.36, .875), (.425, .875)])
    arrow([(.36, .675), (.425, .675)])
    ax.plot(.36, .775, "o", ms=4, color=TEAL)

    # The atlas is a separate donor and never an exposed-cell measurement.
    ax.plot([.023, .023], [.145, .435], color=TEAL, lw=3)
    txt(.04, .425, "Measured MERFISH atlas" if en else "Atlas MERFISH medido", 20, "bold")
    txt(.04, .35,
        "One unexposed male donor\nOriginal ACA + HY cells\nDeposited coordinates" if en else
        "Un macho sin exposición\nCélulas ACA + HY originales\nCoordenadas depositadas",
        16.5)
    txt(.365, .385,
        "Measured anatomy\nand detection" if en else
        "Anatomía medida\ny detección", 20, "bold")
    txt(.365, .245,
        "Physical sections and layers/nuclei\nPositive / all measured cells" if en else
        "Cortes físicos y capas/núcleos\nPositivas / células medidas",
        16.5)
    arrow([(.305, .29), (.345, .29)])

    txt(.685, .43,
        "Guarded Bayesian\nvalidation" if en else
        "Validación bayesiana\ncon separación", 19, "bold")
    txt(.685, .292,
        "Whole sections + neighbour guards\nIntervals and simple baselines" if en else
        "Cortes completos + vecinos\nIntervalos y referencias simples",
        15.8)
    txt(.685, .17,
        "Illustrative regional scenario" if en else
        "Escenario regional ilustrativo", 18, "bold", RUST)
    txt(.685, .105,
        "Uniform bulk-effect multiplier\nPositions and zeros retained" if en else
        "Multiplicador bulk uniforme\nPosiciones y ceros conservados",
        15.8)
    # The atlas informs both modeled reference estimates and the scenario.
    arrow([(.605, .285), (.64, .285)])
    ax.plot([.64, .64], [.13, .39], color=TEAL, lw=1.7)
    arrow([(.64, .39), (.67, .39)])
    arrow([(.64, .13), (.67, .13)])
    ax.plot(.64, .285, "o", ms=4, color=TEAL)

    # Bulk effects enter only the illustrative scenario. They do not train the
    # unexposed spatial Bayesian models or become measured cell responses.
    arrow([(.84, .59), (.98, .59), (.98, .13), (.95, .13)], RUST, True)
    txt(.977, .485, "Bulk effect only" if en else "Solo efecto bulk", 13, "normal", RUST,
        rotation=90, ha="left", bbox={"facecolor":"white", "edgecolor":"none", "pad":2})
    txt(.02, .035,
        "Solid: analysis flow     Dashed: assumed regional transfer, not an observed cellular response" if en else
        "Continua: flujo analítico     Discontinua: transferencia regional supuesta, no respuesta celular observada",
        12.5, color="#526d73")

    # Catch figure-edge clipping; arrow/text relationships also receive visual QA.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    clipped = [t.get_text() for t in ax.texts
               if not bounds.contains(*t.get_window_extent(renderer).get_points()[0])
               or not bounds.contains(*t.get_window_extent(renderer).get_points()[1])]
    if clipped:
        raise ValueError("Workflow text outside figure: " + repr(clipped))
    export(fig, stem, lang, caption,
           ["results/ozone/methods.json", "results/hypothalamus/methods.json",
            "metadata/spatial_sections/methods.json", "metadata/spatial_bayes/methods.json",
            "metadata/spatial_detection_bayes/methods.json"])

def go_plot(lang):
    p="results/advanced_omics/GO_multi_contrast_NES_source.tsv"
    d=pd.read_csv(ROOT/p,sep="\t")
    # The fixed existing display selection, reduced only for slide legibility.
    selected=list(d["GO_ID"].drop_duplicates())[:8]
    order=["hypothalamus","cortex","ozone_M3_E3","ozone_M3_E4","ozone_M17_E3","ozone_M17_E4"]
    sub=d[d.GO_ID.isin(selected)]
    values=sub.pivot(index="GO_ID",columns="contrast",values="NES").reindex(index=selected,columns=order)
    eligible=sub.pivot(index="GO_ID",columns="contrast",values="eligible_tested").reindex(index=selected,columns=order)
    q=sub.pivot(index="GO_ID",columns="contrast",values="algorithmic_FDR").reindex(index=selected,columns=order)
    names=sub.drop_duplicates("GO_ID").set_index("GO_ID")["term_name"]
    arr=values.to_numpy(float); arr[~eligible.fillna(False).to_numpy(bool)]=np.nan
    cmap=plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#e4e7e7")
    vmax=max(1,float(np.nanmax(abs(arr))))
    fig,ax=plt.subplots(figsize=(14,7.1))
    im=ax.imshow(arr,cmap=cmap,vmin=-vmax,vmax=vmax,aspect="auto")
    ax.set_yticks(range(len(selected)),[textwrap.fill(names[g],42) for g in selected],fontsize=16)
    ax.set_xticks(range(6),["HY","ACA","O₃\n3M/E3","O₃\n3M/E4","O₃\n17M/E3","O₃\n17M/E4"],fontsize=17)
    for j in range(len(selected)):
        for k in range(6):
            if np.isfinite(arr[j,k]) and q.iloc[j,k]<.05:
                ax.text(k,j,"•",ha="center",va="center",color="white" if abs(arr[j,k])>vmax*.55 else "black",fontsize=24)
    fig.colorbar(im,ax=ax,pad=.025,label="NES")
    ax.set_title("Coordinated ranked-gene patterns" if lang=="en" else "Patrones coordinados en el orden génico",pad=18)
    fig.subplots_adjust(left=.37,right=.94,bottom=.14,top=.89)
    cap=("First 8 terms in the project's fixed 22-term display selection. Dot: algorithmic GSEA FDR<0.05. Gray: ineligible/missing. Shared color limits, no new enrichment test; correlated terms do not represent independent mechanisms." if lang=="en" else "Primeros 8 términos de la selección gráfica fija de 22 términos del proyecto. Punto: FDR algorítmica GSEA<0.05. Gris: no elegible/ausente. Escala común, sin repetir enriquecimiento; los términos correlacionados no representan mecanismos independientes.")
    export(fig,"context_go_patterns",lang,cap,[p,"metadata/advanced_omics/methods.json"])

def animal_plot(lang):
    p="results/advanced_omics/ozone_interaction_animal_profiles.tsv"
    d=pd.read_csv(ROOT/p,sep="\t")
    fig,axs=plt.subplots(1,2,figsize=(14,6),sharey=False)
    groups=["M3_E3","M3_E4","M17_E3","M17_E4"]
    for ax,gene in zip(axs,["Ryr3","Rps24"]):
        sub=d[d.gene_symbol==gene]
        for k,exposure in enumerate(["FA","O3"]):
            for j,group in enumerate(groups):
                z=sub[(sub.stratum==group)&(sub.exposure==exposure)].sort_values("sample_id")
                x=j+(-.14 if k==0 else .14)
                offsets=np.linspace(-.035,.035,len(z))
                ax.scatter(x+offsets,z.log2_1_plus_normalized_count,s=95,marker="o" if k==0 else "s",color=TEAL if k==0 else RUST,label=exposure if j==0 else None,zorder=3)
                ax.plot([x-.06,x+.06],[z.log2_1_plus_normalized_count.mean()]*2,color=INK,lw=2)
        ax.set_xticks(range(4),["3M\nE3","3M\nE4","17M\nE3","17M\nE4"])
        ax.set_title(gene,fontstyle="italic",pad=15)
        ax.set_ylabel("log₂(1 + normalized count)" if lang=="en" else "log₂(1 + conteo normalizado)")
        ax.grid(axis="y",alpha=.15); ax.legend(frameon=False)
    fig.subplots_adjust(left=.09,right=.98,bottom=.17,top=.88,wspace=.28)
    cap=("All 24 sequenced animals per gene. Small fixed horizontal offsets separate overlapping symbols; they do not encode another variable. Horizontal segments are descriptive means of plotted transformed values, without pairing or inferential error bars. These genes belong to the eight age-dependent exposure LRT discoveries." if lang=="en" else "Los 24 animales secuenciados por gen. Los pequeños desplazamientos horizontales fijos separan símbolos superpuestos; no codifican otra variable. Los segmentos son medias descriptivas de los valores transformados, sin pareamiento ni barras inferenciales. Estos genes pertenecen a los ocho hallazgos LRT de respuesta a exposición dependiente de edad.")
    export(fig,"context_ozone_animals",lang,cap,[p,"results/advanced_omics/age_response_FDR_genes.tsv"])

def main():
    for lang in ("en","es"):
        if lang=="en":
            cap="Two distinct evidence branches. Four experimental studies feed original-study reproduction and new study-specific analyses. A separate unexposed MERFISH donor feeds measured anatomy/detection and guarded Bayesian validation. A dashed assumption-based regional bulk-effect transfer enters only the illustrative scenario. Solid arrows denote computational flow, not a demonstrated biological causal pathway. O₃ input comprises 24 paired-end libraries."
        else:
            cap="Dos ramas de evidencia diferenciadas. Cuatro estudios experimentales alimentan la reproducción y los análisis nuevos por estudio. Un donante MERFISH independiente sin exposición alimenta la anatomía/detección medidas y la validación bayesiana con separación. Una transferencia supuesta y discontinua del efecto bulk regional entra solo al escenario ilustrativo. Las flechas continuas representan flujo computacional, no una vía causal biológica demostrada. La entrada O₃ incluye 24 bibliotecas paired-end."
        workflow("context_workflow",lang,cap)
        go_plot(lang); animal_plot(lang)
    (OUT/"context_figure_manifest.json").write_text(json.dumps(MANIFEST,indent=2,ensure_ascii=False)+"\n")
    print(f"Exported {len(MANIFEST)} language-specific slide figures")

if __name__=="__main__": main()
