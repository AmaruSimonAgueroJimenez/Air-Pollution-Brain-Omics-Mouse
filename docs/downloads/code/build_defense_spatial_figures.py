#!/usr/bin/env python3
"""Presentation-size spatial figures from frozen, previously computed outputs.

No observations, models, priors, folds or effect estimates are changed/refitted.
Outputs are bilingual PNG/SVG/PDF plus source-file hashes and scoped captions.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.ticker import MultipleLocator, PercentFormatter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/defense_slides/figures'
MANIFEST = ROOT / 'results/defense_slides/spatial_figure_manifest.json'
DARK, TEAL, RUST = '#123d46', '#376E75', '#B65E3C'
GRAY, PALE = '#CDD5D8', '#D4E4E8'
STYLE = {
    'font.family': 'DejaVu Sans', 'font.size': 20,
    'axes.labelsize': 20, 'axes.titlesize': 23, 'axes.titleweight': 'bold',
    'xtick.labelsize': 18, 'ytick.labelsize': 18, 'legend.fontsize': 18,
    'text.color': DARK, 'axes.labelcolor': DARK, 'xtick.color': DARK,
    'ytick.color': DARK, 'axes.edgecolor': '#60767B',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.facecolor': 'white', 'figure.facecolor': 'white',
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'svg.hashsalt': 'brain-defense-spatial-20261028', 'savefig.bbox': None,
}

def read_json(path):
    return json.loads((ROOT / path).read_text())

def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()

def export(fig, stem, es, caption_en, caption_es, sources, alt_en, alt_es, manifest):
    name = stem + ('_es' if es else '')
    files = {}
    # Redraw into the wider evidence canvas, without changing observations.
    fig.set_size_inches(16, 7.2)
    for label in list(fig.texts):
        label.remove()
    fig.subplots_adjust(top=.83, bottom=max(.31, fig.subplotpars.bottom))
    fig.canvas.draw()
    for ax in fig.axes:
        for axis, limits in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            low, high = sorted(limits)
            for tick in axis.get_major_ticks():
                if tick.get_loc() < low-1e-9 or tick.get_loc() > high+1e-9:
                    tick.label1.set_visible(False)
                    tick.label2.set_visible(False)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    from matplotlib.text import Text
    overflows = []
    for artist in fig.findobj(Text):
        if not artist.get_visible() or not artist.get_text():
            continue
        box = artist.get_window_extent(renderer)
        width, height = fig.canvas.get_width_height()
        if box.x0 < -1 or box.y0 < -1 or box.x1 > width+1 or box.y1 > height+1:
            overflows.append(artist.get_text())
    if overflows:
        print(name, 'TEXT OVERFLOW', overflows, flush=True)
    for ext in ('png', 'svg', 'pdf'):
        path = OUT / f'{name}.{ext}'
        metadata = {'CreationDate': None, 'ModDate': None} if ext == 'pdf' else ({'Date': None} if ext == 'svg' else {})
        fig.savefig(path, dpi=140, metadata=metadata)
        files[ext] = path.relative_to(ROOT).as_posix()
    manifest[name] = {
        'stem': name, 'language': 'es' if es else 'en',
        'caption_en': caption_en, 'caption_es': caption_es,
        'alt_en': alt_en, 'alt_es': alt_es,
        'source_paths': sources,
        'source_sha256': {p: sha(p) for p in sources},
        'files': files, 'sha256': {ext: sha(path) for ext, path in files.items()},
        'text_overflow': overflows, 'figure_size_inches': [16, 7.2], 'png_dimensions': [2240, 1008],
        'analysis_scope': 'Focused presentation of frozen outputs; no refit, sampling or new inference',
    }
    plt.close(fig)


def footer(fig, text, y=.025):
    fig.text(.5, y, text, ha='center', va='bottom', fontsize=16, color=DARK)


def measured(region, es, manifest):
    T = lambda en, spa: spa if es else en
    gene = 'Fosb' if region == 'ACA' else 'Tll1'
    summary = read_json(f'results/spatial_sections/{region}/summary.json')
    section = summary['example_section']['label']
    cellfile = f'results/spatial_sections/{region}/cells_with_sections.tsv.gz'
    binfile = f'results/spatial_sections/{region}/detection_by_bins.tsv.gz'
    cells = pd.read_csv(ROOT / cellfile, sep='\t')
    cells = cells.loc[cells.brain_section_label.eq(section)].copy()
    bins = pd.read_csv(ROOT / binfile, sep='\t')
    bins = bins.loc[bins.gene.eq(gene) & bins.section_label.eq(section)].copy()
    assert len(cells) == summary['example_section']['n']
    assert int(bins.n.sum()) == len(cells)
    assert np.allclose(bins.detection_fraction, bins.positive / bins.n)
    anatomy = sorted(cells.substructure.unique())
    cmap = plt.get_cmap('tab10')
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    fig.subplots_adjust(left=.08, right=.875, bottom=.31, top=.78, wspace=.38)
    region_name = T('Anterior cingulate cortex', 'Corteza cingulada anterior') if region == 'ACA' else T('Hypothalamus', 'Hipotálamo')
    fig.suptitle(f'{region_name}  ·  {gene}', y=.96, fontsize=27, weight='bold')
    fig.text(.5, .90, f'{section}  ·  n = {len(cells):,}  ·  ' + T('one unexposed donor', 'un donante sin exposición'), ha='center', fontsize=21)
    for j, label in enumerate(anatomy):
        a = cells.loc[cells.substructure.eq(label)]
        axes[0].scatter(a.experimental_x, a.experimental_y, s=9, c=[cmap(j % 10)], linewidths=0, rasterized=True)
    axes[0].set_title(T('Measured anatomy', 'Anatomía medida'), pad=15)
    rects = [Rectangle((r.x0_mm, r.y0_mm), .096, .096) for r in bins.itertuples()]
    colors = [GRAY if r.positive == 0 else plt.get_cmap('viridis')(r.detection_fraction) for r in bins.itertuples()]
    pc = PatchCollection(rects, facecolors=colors, edgecolors='white', linewidths=.2)
    axes[1].add_collection(pc)
    axes[1].set_title(T('Observed detection\n100 µm occupied bins', 'Detección observada\nCasillas ocupadas de 100 µm'), pad=15, fontsize=21)
    norm = Normalize(0, 1)
    sm = plt.cm.ScalarMappable(norm=norm, cmap='viridis')
    cb = fig.colorbar(sm, ax=axes[1], fraction=.046, pad=.035)
    cb.set_label(T('Positive / measured', 'Positivas / medidas'), labelpad=12)
    cb.ax.yaxis.set_major_formatter(PercentFormatter(1))
    xmin = min(cells.experimental_x.min(), bins.x0_mm.min()) - .08
    xmax = max(cells.experimental_x.max(), (bins.x0_mm + .1).max()) + .08
    ymin = min(cells.experimental_y.min(), bins.y0_mm.min()) - .08
    ymax = max(cells.experimental_y.max(), (bins.y0_mm + .1).max()) + .08
    for ax in axes:
        ax.set(xlim=(xmin, xmax), ylim=(ymax, ymin),
               xlabel=T('Deposited section X (mm)', 'X del corte depositado (mm)'),
               ylabel=T('Deposited section Y (mm)', 'Y del corte depositado (mm)'))
        ax.set_aspect('equal', adjustable='box')
        ax.xaxis.set_major_locator(MultipleLocator(.5))
        ax.yaxis.set_major_locator(MultipleLocator(.5))
        ax.grid(alpha=.12)
    handles = [Line2D([], [], marker='o', color=cmap(j % 10), linestyle='', markersize=8, label=str(label)) for j, label in enumerate(anatomy)]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .09), ncol=5, frameon=False, fontsize=17, columnspacing=1.2, handletextpad=.35)
    total_pos = int(bins.positive.sum())
    footer(fig, T(f'Raw count > 0: {total_pos:,}/{len(cells):,}. Gray occupied bins: zero detection. White gaps: unobserved. No interpolation.',
                  f'Conteo bruto > 0: {total_pos:,}/{len(cells):,}. Casillas grises: detección cero. Huecos blancos: sin observación. Sin interpolación.'))
    cap_en = (f'{region_name if not es else ("Anterior cingulate cortex" if region == "ACA" else "Hypothalamus")}, {gene}, deposited physical section {section}, selected solely for the largest measured regional cell n ({len(cells):,}). '
              'Left: every measured cell at its deposited experimental x/y coordinate, colored by deposited anatomy; these colors do not encode RNA. Right: positive raw count / all measured regional cells in each occupied 100 × 100 µm experimental-XY bin, including measured zeros in every denominator. '
              f'{total_pos:,}/{len(cells):,} cells are positive in this section. Detection color scale is fixed 0–100%; gray occupied bins have zero positives and white gaps contain no observed bin. Neither interpolation, density estimation nor spatial imputation is applied. '
              'The displayed Y-axis is downward. Experimental slice coordinates are not CCF ML/DV. One unexposed donor, no per-cell tests or PM2.5 response inference.')
    cap_es = (f'{"Corteza cingulada anterior" if region == "ACA" else "Hipotálamo"}, {gene}, corte físico depositado {section}, seleccionado solo por mayor n de células regionales medidas ({len(cells):,}). '
              'Izquierda: todas las células medidas en sus x/y experimentales depositadas, coloreadas por anatomía depositada; estos colores no codifican RNA. Derecha: conteo bruto positivo / todas las células regionales medidas por casilla XY experimental ocupada de 100 × 100 µm, incluidos todos los ceros en el denominador. '
              f'{total_pos:,}/{len(cells):,} células positivas en este corte. Escala de detección fija 0–100%; casillas grises ocupadas sin positivos y huecos blancos sin casilla observada. Sin interpolación, densidad ni imputación espacial. '
              'Y se muestra hacia abajo. Las coordenadas experimentales no son CCF ML/DV. Un donante sin exposición; sin pruebas celulares ni inferencia de respuesta PM2.5.')
    export(fig, f'spatial_{region.lower()}_measured', es, cap_en, cap_es, [cellfile, binfile, f'results/spatial_sections/{region}/summary.json'],
           f'Two maps of measured {region} slice: anatomy at actual cell coordinates and {gene} detection fractions in occupied 100 µm bins.',
           f'Dos mapas del corte {region}: anatomía en coordenadas celulares reales y fracciones de detección de {gene} en casillas ocupadas de 100 µm.', manifest)


def guarded_folds(es, manifest):
    T = lambda en, spa: spa if es else en
    fig, axes = plt.subplots(2, 1, figsize=(16, 9))
    fig.subplots_adjust(left=.11, right=.975, bottom=.26, top=.84, hspace=.90)
    sources = []
    for ax, region in zip(axes, ['ACA', 'HY']):
        p = f'results/spatial_bayes/{region}/section_folds.tsv'
        a = f'results/spatial_bayes/{region}/fold_audit.json'
        sources.extend([p, a])
        sections = pd.read_csv(ROOT / p, sep='\t').sort_values('order')
        labels = sections.section.tolist()
        audits = read_json(a)
        for row, audit in enumerate(audits):
            for j, label in enumerate(labels):
                test, guard = label in audit['test_sections'], label in audit['guard_sections']
                color = TEAL if test else RUST if guard else '#EDF2F3'
                ax.add_patch(Rectangle((j-.5, row-.5), 1, 1, facecolor=color, edgecolor='white', linewidth=1.5))
                if test or guard:
                    ax.text(j, row, 'T' if test else 'G', ha='center', va='center', color='white', weight='bold', fontsize=18)
        ax.set(xlim=(-.5, len(labels)-.5), ylim=(4.5, -.5))
        ax.set_yticks(range(5), [T(f'Fold {i}', f'Bloque {i}') for i in range(1, 6)])
        ax.set_xticks(range(len(labels)), [s.rsplit('.', 1)[1] for s in labels])
        ax.set_title(f'{region}  ·  {len(labels)} ' + T('physical sections', 'cortes físicos'), loc='left', pad=12)
        ax.set_xlabel(T('Deposited section index, ordered by acquisition AP', 'Índice de corte depositado, ordenado por AP de adquisición'), labelpad=9)
        for spine in ax.spines.values(): spine.set_visible(False)
        ax.tick_params(length=0)
    fig.suptitle(T('Five whole-section folds with neighboring guards', 'Cinco bloques de cortes completos con vecinos de guarda'), fontsize=27, weight='bold', y=.97)
    handles = [Patch(facecolor='#EDF2F3', edgecolor='#72858A', label=T('Training', 'Entrenamiento')),
               Patch(facecolor=TEAL, label=T('T: held-out test', 'T: prueba reservada')),
               Patch(facecolor=RUST, label=T('G: excluded guard', 'G: guarda excluida'))]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .075), ncol=3, frameon=False, fontsize=20)
    footer(fig, T('Every test section is entirely held out. Training preprocessing uses neither test nor guard cells. One donor.',
                  'Cada corte de prueba se reserva completo. El preprocesamiento no usa células de prueba ni guarda. Un donante.'))
    export(fig, 'spatial_guarded_folds', es,
           'Actual fixed five-fold assignments for all 18 ACA and 20 HY physical sections. Whole sections form contiguous test blocks ordered by estimated acquisition AP. The nearest represented section immediately before and after each test block is excluded as a guard. Section indices need not be consecutive because absent regional sections are not filled. Train-only anatomy vocabularies, coordinate scaling and RBF centers exclude test and guard cells. The same folds support Gaussian and detection models. This assesses within-donor transport, not independent animal replication.',
           'Asignaciones reales de cinco bloques para los 18 cortes ACA y 20 HY. Los cortes completos forman bloques contiguos de prueba según AP estimada de adquisición. El corte regional representado más cercano antes y después se excluye como guarda. Los índices no tienen que ser consecutivos; no se rellenan cortes ausentes. Vocabulario anatómico, escalado y centros RBF se aprenden solo con entrenamiento. Ambos modelos comparten bloques. Evalúa transporte dentro de un donante, no replicación animal independiente.',
           sources, 'Training, test and guard membership by actual physical section and fold in ACA and HY.',
           'Pertenencia de entrenamiento, prueba y guarda por corte físico y bloque en ACA e HY.', manifest)


def gaussian_rmse(es, manifest):
    T = lambda en, spa: spa if es else en
    summary = read_json('results/spatial_bayes/summary.json')
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    fig.subplots_adjust(left=.16, right=.965, bottom=.28, top=.81, wspace=.48)
    for ax, region in zip(axes, ['ACA', 'HY']):
        rows = summary['regions'][region]['heldout_results']
        ratios = np.array([r['rmse']/r['anatomy_baseline_rmse'] for r in rows])
        wins = int((ratios < 1).sum())
        ax.axvline(1, color=DARK, lw=1.6, linestyle='--', zorder=1)
        for i, val in enumerate(ratios):
            color = TEAL if val < 1 else RUST
            ax.plot([1, val], [i, i], color='#A9B9BD', linewidth=2, zorder=2)
            ax.scatter(val, i, c=color, marker='o' if val < 1 else 'X', s=120, linewidths=.8, edgecolors='white', zorder=3)
        ax.set_yticks(range(10), [r['gene'] for r in rows])
        ax.set(ylim=(9.6, -.6), xlim=(.95, 1.15))
        ax.set_xticks(np.arange(.95, 1.151, .05))
        ax.grid(axis='x', alpha=.2)
        ax.set_title(f'{region}: {wins}/10 ' + T('improve', 'mejoran'), pad=16)
        ax.set_xlabel(T('Spatial RMSE / anatomy-only RMSE', 'RMSE espacial / RMSE de solo anatomía'), labelpad=13)
    fig.suptitle(T('The spatial basis seldom improves held-out prediction', 'Las bases espaciales rara vez mejoran la predicción reservada'), fontsize=26, weight='bold', y=.96)
    fig.legend(handles=[Line2D([], [], marker='o', color=TEAL, linestyle='', markersize=10, label=T('Lower error than anatomy-only', 'Menor error que solo anatomía')),
                        Line2D([], [], marker='X', color=RUST, linestyle='', markersize=10, label=T('Equal or higher error', 'Error igual o mayor'))],
               loc='lower center', bbox_to_anchor=(.5, .065), ncol=2, frameon=False, fontsize=20)
    footer(fig, T('Reference = 1. Same guarded folds, measured cells and ln(1 + count) response for both predictors.',
                  'Referencia = 1. Mismos bloques protegidos, células medidas y respuesta ln(1 + conteo) para ambos predictores.'))
    export(fig, 'spatial_gaussian_rmse', es,
           'Per-gene held-out RMSE ratio, primary spatial Gaussian baseline divided by anatomy-only Bayesian baseline, both using identical whole-section guard folds and section nuisance priors. Ratios below 1 improve prediction: 4/10 ACA genes and 1/10 HY genes. Point positions retain exact ratios; the nonzero axis range is explicit and the equality reference is shown. Genes remain in frozen source priority order. Error is measured on natural ln(1+raw imaging count), not raw molecule counts or fold change.',
           'Razón de RMSE por gen: referencia gaussiana espacial principal / referencia bayesiana de solo anatomía, con los mismos bloques de cortes protegidos y priors de corte. Razones menores que 1 mejoran: 4/10 genes ACA y 1/10 HY. Los puntos conservan razones exactas; el rango distinto de cero es explícito y se muestra la igualdad. Orden génico fijado por prioridad de origen. Error en ln(1+conteo bruto), no en moléculas ni fold change.',
           ['results/spatial_bayes/summary.json', 'results/spatial_bayes/ACA/heldout_metrics.tsv', 'results/spatial_bayes/HY/heldout_metrics.tsv'],
           'Two dot plots compare spatial Gaussian held-out RMSE with anatomy-only prediction for ten genes per region.',
           'Dos gráficos de puntos comparan RMSE gaussiano espacial con solo anatomía para diez genes por región.', manifest)


def gaussian_diagnostics(es, manifest):
    T = lambda en, spa: spa if es else en
    summary = read_json('results/spatial_bayes/summary.json')
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    fig.subplots_adjust(left=.11, right=.96, bottom=.30, top=.81, wspace=.36)
    for region, marker, color in [('ACA', 'o', TEAL), ('HY', 's', RUST)]:
        rows = summary['regions'][region]['heldout_results']
        overall = np.array([r['predictive_95_coverage'] for r in rows])*100
        positive = np.array([r['coverage_positive'] for r in rows])*100
        negative = np.array([r['negative_predictive_probability_mean'] for r in rows])*100
        axes[0].scatter(overall, positive, marker=marker, color=color, s=95, edgecolors='white', linewidths=.8, label=f'{region} (10 genes)')
        axes[1].scatter(negative, overall, marker=marker, color=color, s=95, edgecolors='white', linewidths=.8)
        example = 0
        name = rows[example]['gene']
        axes[0].annotate(name, (overall[0], positive[0]), xytext=(-38, 18), textcoords='offset points', fontsize=19, color=color,
                         arrowprops={'arrowstyle':'-', 'color':color})
        offset = (-68, -34) if region == 'ACA' else (-65, 19)
        axes[1].annotate(name, (negative[0], overall[0]), xytext=offset, textcoords='offset points', fontsize=19, color=color,
                         arrowprops={'arrowstyle':'-', 'color':color})
    axes[0].axvline(95, color=DARK, lw=1.4, ls='--')
    axes[1].axhline(95, color=DARK, lw=1.4, ls='--')
    axes[0].set(xlim=(88, 101), ylim=(-3, 103),
                xlabel=T('Overall interval inclusion (%)', 'Inclusión global en intervalo (%)'),
                ylabel=T('Positive-only inclusion (%)\nOutcome-conditioned diagnostic', 'Inclusión solo en positivos (%)\nDiagnóstico condicionado al resultado'))
    axes[1].set(xlim=(0, 55), ylim=(80, 101),
                xlabel=T('Predictive mass below zero (%)', 'Masa predictiva bajo cero (%)'),
                ylabel=T('Overall interval inclusion (%)', 'Inclusión global en intervalo (%)'))
    axes[0].set_title(T('Marginal versus positive observations', 'Observaciones globales y positivas'), pad=16)
    axes[1].set_title(T('Impossible negative support', 'Soporte negativo imposible'), pad=16)
    axes[1].set_xticks(np.arange(0, 51, 10))
    for ax in axes: ax.grid(alpha=.18)
    fig.suptitle(T('Gaussian intervals reveal different limitations', 'Los intervalos gaussianos muestran limitaciones distintas'), fontsize=26, weight='bold', y=.96)
    fig.legend(handles=axes[0].get_legend_handles_labels()[0] + [Line2D([], [], color=DARK, ls='--', label=T('Nominal 95%: overall only', '95% nominal: solo global'))],
               loc='lower center', bbox_to_anchor=(.5, .095), ncol=3, frameon=False, fontsize=20)
    footer(fig, T('Each point is one gene. Positive-only inclusion has no nominal 95% target. Negative transformed counts are physically impossible.',
                  'Cada punto es un gen. La inclusión solo en positivos no tiene meta nominal de 95%. Los conteos transformados negativos son imposibles.'), y=.025)
    export(fig, 'spatial_gaussian_diagnostics', es,
           'Every point is one of ten fixed genes in each region. Left compares empirical overall inclusion in nominal 95% new-cell/new-section predictive intervals with inclusion conditioned on a positive realized observation. The conditional positive-only diagnostic is not required to equal 95% even under marginal calibration; the nominal reference therefore appears only for overall inclusion. Right shows the mean model probability below zero against overall interval inclusion; negative support is incompatible with the nonnegative ln(1+count) response. Tll1 and Fosb are labeled examples. No intervals are clipped or recalibrated. Different plotted quantities retain different axis scales. Within-donor descriptive validation, not independent-cell or animal inference.',
           'Cada punto es uno de diez genes fijos por región. Izquierda compara inclusión global observada en intervalos predictivos nominales de 95% para célula/corte nuevos con inclusión condicionada a una observación positiva. Este diagnóstico condicionado no tiene que igualar 95% bajo calibración marginal; la referencia nominal aparece solo para inclusión global. Derecha muestra probabilidad media modelada bajo cero frente a inclusión global; el soporte negativo es incompatible con ln(1+conteo), no negativo. Tll1 y Fosb son ejemplos rotulados. No se truncan ni recalibran intervalos. Las cantidades usan escalas de eje distintas. Validación descriptiva de un donante, sin inferencia de células o animales independientes.',
           ['results/spatial_bayes/summary.json', 'metadata/spatial_bayes/methods.json'],
           'Two diagnostic scatter plots distinguish overall from positive-only interval inclusion and quantify impossible negative predictive probability.',
           'Dos gráficos distinguen inclusión global y condicionada a positivos, y cuantifican probabilidad predictiva negativa imposible.', manifest)


def beta_performance(es, manifest):
    T = lambda en, spa: spa if es else en
    summary = read_json('results/spatial_detection_bayes/summary.json')
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    fig.subplots_adjust(left=.09, right=.97, bottom=.32, top=.80, wspace=.34)
    x = np.arange(2); width = .30
    for j, (key, col, hatch, name) in enumerate([
        ('brier', TEAL, '', T('Anatomical model', 'Modelo anatómico')),
        ('baseline_brier', '#C9D5D8', '//', T('Regional baseline', 'Referencia regional'))]):
        values = [summary[r][key] for r in ['ACA', 'HY']]
        bars = axes[0].bar(x + (j-.5)*width, values, width, color=col, edgecolor=DARK, linewidth=.6, hatch=hatch, label=name)
        for bar, val in zip(bars, values):
            axes[0].text(bar.get_x()+bar.get_width()/2, val+.004, f'{val:.4f}', ha='center', fontsize=19)
    axes[0].set(xticks=x, xticklabels=['ACA', 'HY'], ylim=(0, .175), ylabel=T('Held-out Brier loss', 'Pérdida de Brier reservada'))
    axes[0].set_title(T('Lower Brier is better', 'Un Brier menor es mejor'), pad=16)
    axes[0].yaxis.set_major_locator(MultipleLocator(.05)); axes[0].grid(axis='y', alpha=.18)
    for region, marker, color, offset in [('ACA', 'o', TEAL, (12, -10)), ('HY', 's', RUST, (12, -20))]:
        d=summary[region]; xx=100*d['count_pi_mean_fraction_width']; yy=100*d['count_pi_coverage']
        axes[1].scatter(xx, yy, marker=marker, color=color, s=160, edgecolors='white', linewidths=1)
        axes[1].annotate(f'{region}\n{yy:.1f}%', (xx, yy), xytext=offset, textcoords='offset points', fontsize=20, color=color)
    axes[1].axhline(95, color=DARK, lw=1.4, ls='--')
    axes[1].text(3, 95.4, T('95% nominal', '95% nominal'), fontsize=18)
    axes[1].set(xlim=(0, 86), ylim=(90, 101),
                xlabel=T('Mean count-interval width / n\nPercentage points', 'Ancho medio del intervalo de conteo / n\nPuntos porcentuales'),
                ylabel=T('Held-out count coverage (%)', 'Cobertura de conteos\nreservados (%)'))
    axes[1].set_title(T('High coverage with wide intervals', 'Cobertura alta con intervalos amplios'), pad=16)
    axes[1].set_xticks(np.arange(0, 81, 20))
    axes[1].grid(alpha=.18)
    fig.suptitle(T('A small ACA gain and worse HY prediction', 'Mejora pequeña en ACA y peor predicción HY'), fontsize=25, weight='bold', y=.96)
    fig.legend(handles=axes[0].get_legend_handles_labels()[0], loc='lower center', bbox_to_anchor=(.5, .09), ncol=2, frameon=False, fontsize=20)
    footer(fig, T('All ten genes pooled per region. HY prior-only: 79/186 anatomical section groups per gene, covering 10,416 actual cells.',
                  'Diez genes agregados por región. HY solo con prior: 79/186 grupos corte–anatomía por gen, con 10,416 células reales.'), y=.025)
    export(fig, 'spatial_beta_performance', es,
           'Primary 64×48 hierarchical Beta–Binomial held-out performance pooled across ten genes per region. Left is cell-weighted Brier loss, with zero baseline on the bar axis, compared with the nonanatomical regional model trained on identical allowed sections and fixed priors. ACA model .1372514 versus regional .1405089; HY .0924574 versus .0362057. Right compares mean predictive count-interval width divided by the observed group n (41.9/61.5 percentage points) with unweighted group coverage (97.89/99.35%). Intervals predict positive counts in a new physical section with its actual n. High coverage with wide intervals does not establish precise prediction. All groups are retained, including 79/186 unique HY section×anatomy groups per gene with prior-only prediction (10,416 of37,811 actual cells). Ten-gene aggregate evaluations are correlated and are not independent samples.',
           'Rendimiento reservado principal 64×48 del modelo Beta–Binomial jerárquico, agregado en diez genes por región. Izquierda: Brier ponderado por células, barras desde cero, comparado con referencia regional no anatómica entrenada con los mismos cortes permitidos y priors. ACA .1372514 frente a .1405089; HY .0924574 frente a .0362057. Derecha: ancho medio del intervalo predictivo de conteos dividido por n observado (41.9/61.5 puntos porcentuales) frente a cobertura no ponderada por grupo (97.89/99.35%). Predicen conteos positivos de un corte nuevo con su n real. Cobertura alta e intervalos amplios no establecen predicción precisa. Se conservan todos los grupos, incluidos 79/186 grupos HY corte×anatomía por gen solo con prior (10,416 de37,811 células reales). Las evaluaciones agregadas en diez genes son correlacionadas, no muestras independientes.',
           ['results/spatial_detection_bayes/summary.json', 'metadata/spatial_detection_bayes/methods.json'],
           'Brier comparison for ACA and HY, and predictive count coverage plotted against interval width.',
           'Comparación de Brier ACA/HY y cobertura predictiva de conteos frente al ancho del intervalo.', manifest)


def beta_posterior(region, es, manifest):
    T = lambda en, spa: spa if es else en
    gene = 'Fosb' if region == 'ACA' else 'Tll1'
    summaryfile = f'results/spatial_sections/{region}/summary.json'
    section = read_json(summaryfile)['example_section']['label']
    knownfile = f'results/spatial_detection_bayes/{region}/known_section_posterior.tsv.gz'
    anatomyfile = f'results/spatial_detection_bayes/{region}/anatomy_posterior.tsv.gz'
    known = pd.read_csv(ROOT / knownfile, sep='\t')
    known = known.loc[known.gene.eq(gene) & known.label.eq(section)].sort_values('anatomy_label').copy()
    anat = pd.read_csv(ROOT / anatomyfile, sep='\t')
    anat = anat.loc[anat.gene.eq(gene)].set_index('label').loc[known.anatomy_label].reset_index()
    assert len(known) == 10 and np.allclose(known.observed_fraction, known.positive/known.n)
    ymax = 9.65; upper = .45 if region == 'ACA' else .22
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    fig.subplots_adjust(left=.23, right=.965, bottom=.30, top=.79, wspace=.36)
    ys = np.arange(10)
    axes[0].errorbar(100*known.p_mean, ys, xerr=100*np.vstack([known.p_mean-known.p_ci_low, known.p_ci_high-known.p_mean]),
                     fmt='o', color=TEAL, ecolor=TEAL, elinewidth=2, capsize=4, markersize=7)
    axes[0].scatter(100*known.observed_fraction, ys, color=RUST, marker='D', s=52, zorder=4)
    axes[0].set_yticks(ys, [f'{r.anatomy_label}  ({r.positive}/{r.n})' for r in known.itertuples()])
    for y, row in zip(ys, anat.itertuples()):
        axes[1].plot(100*np.array([row.new_section_p_ci_low, row.new_section_p_ci_high]), [y, y], color=PALE, linewidth=13, solid_capstyle='butt')
        axes[1].plot(100*np.array([row.mu_ci_low, row.mu_ci_high]), [y, y], color=TEAL, linewidth=3)
        axes[1].scatter(100*row.mu_mean, y, color=DARK, marker='s', s=55, zorder=4)
    axes[1].set_yticks(ys, [f'S = {r.n_sections}' for r in anat.itertuples()])
    for ax in axes:
        ax.set(xlim=(0, upper*100), ylim=(ymax, -.65), xlabel=T('Detection probability / fraction (%)', 'Probabilidad / fracción de detección (%)'))
        ax.grid(axis='x', alpha=.18)
        ax.set_xticks(np.arange(0, upper*100+.01, 10 if region=='ACA' else 5))
    axes[0].set_title(T('p of the measured section', 'p del corte medido'), pad=16, fontsize=21)
    axes[1].set_title(T('Anatomy μ and new-section p', 'μ anatómica y p de un corte nuevo'), pad=16, fontsize=21)
    fig.suptitle(f'{region}  ·  {gene}  ·  ' + T('Distinct posterior targets', 'Objetivos posteriores distintos'), fontsize=26, weight='bold', y=.97)
    fig.text(.5, .89, f'{section}  ·  ' + T('positives/n at left; S = physical sections', 'positivos/n a la izquierda; S = cortes físicos'), ha='center', fontsize=19)
    handles=[Line2D([],[],color=TEAL,marker='o',label=T('Known-section p, 95% credible', 'p del corte, 95% creíble')),
             Line2D([],[],color=RUST,marker='D',linestyle='',label=T('Observed positive fraction', 'Fracción positiva observada')),
             Line2D([],[],color=TEAL,marker='s',label=T('Anatomy μ, 95% credible', 'μ anatómica, 95% creíble')),
             Line2D([],[],color=PALE,linewidth=9,label=T('New-section p, 95% predictive', 'p de corte nuevo, 95% predictivo'))]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .075), ncol=2, frameon=False, fontsize=18, columnspacing=2)
    footer(fig, T('Fitted reference summaries, not held-out predictions. Zero observations remain in n. One unexposed donor.',
                  'Resúmenes ajustados de referencia, no predicciones reservadas. Se conservan ceros en n. Un donante sin exposición.'), y=.025)
    cap_en = (f'{region}/{gene}, the same largest-n physical slice {section} used in the measured figure. Left: observed positive fractions (diamonds), fitted known-section probability means and95% credible intervals under the full-reference hyperposterior. Labels give actual positive/n, including zeros. '
              'Right: the corresponding anatomy mean μ and95% credible intervals over all its measured sections, distinct from95% predictive intervals for a new section propensity p; S is the actual physical-section count. These are different posterior targets. Known-section p intervals use8192 seeded mixture draws, μ uses weighted grid quantiles, and full-data new-section p intervals invert the finite-grid Beta-mixture CDF. '
              'A new-section propensity interval does not include binomial count sampling and is not the held-out count interval. Both panels are full-reference fitted summaries, not held-out predictions or exposure effects. The same anatomical rows are retained in each panel. No model refit or re-selection performed for this presentation.')
    cap_es = (f'{region}/{gene}, mismo corte físico de mayor n {section} usado en la figura medida. Izquierda: fracciones positivas observadas (rombos), medias ajustadas de probabilidad del corte conocido e intervalos creíbles95% bajo la hiperposterior de referencia completa. Etiquetas con positivos/n reales, incluidos ceros. '
              'Derecha: media anatómica μ e intervalos creíbles95% con todos sus cortes medidos, distintos de los intervalos predictivos95% para propensión p de un corte nuevo; S cuenta cortes físicos. Son objetivos posteriores diferentes. Los intervalos p de corte conocido usan8192 sorteos con semilla, μ usa cuantiles ponderados de cuadrícula y p nuevo con datos completos invierte la CDF de mezcla Beta discreta. '
              'Un intervalo de propensión nueva no incluye muestreo binomial de conteos ni es el intervalo de conteos reservados. Ambos paneles son resúmenes ajustados con referencia completa, no predicciones reservadas ni efectos de exposición. Mismas filas anatómicas por panel. Sin reajuste ni selección nueva para esta presentación.')
    export(fig, f'spatial_beta_{region.lower()}_posterior', es, cap_en, cap_es,
           [knownfile, anatomyfile, summaryfile, 'metadata/spatial_detection_bayes/methods.json'],
           f'Two forests distinguish observed and fitted known-section detection from anatomy μ uncertainty and new-section p prediction for {region}/{gene}.',
           f'Dos gráficos distinguen detección observada y p ajustada del corte, incertidumbre de μ anatómica y predicción de p nuevo para {region}/{gene}.', manifest)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {}
    with plt.rc_context(STYLE):
        for es in [False, True]:
            measured('ACA', es, manifest)
            measured('HY', es, manifest)
            guarded_folds(es, manifest)
            gaussian_rmse(es, manifest)
            gaussian_diagnostics(es, manifest)
            beta_performance(es, manifest)
            beta_posterior('ACA', es, manifest)
            beta_posterior('HY', es, manifest)
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n')
    print(f'Created {len(manifest)} focused bilingual figures, each PNG/SVG/PDF. No model refits.')

if __name__ == '__main__':
    main()
