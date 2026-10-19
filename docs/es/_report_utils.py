"""Spanish presentation of the shared report outputs; no analysis is rerun."""
from __future__ import annotations
import html
import importlib.util
from pathlib import Path
import pandas as pd
from IPython.display import Markdown, display

DOCS = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('_omics_report_shared', DOCS.parent / '_report_utils.py')
_shared = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_shared)
ROOT = _shared.ROOT
summary = _shared.summary
read_table = _shared.read_table
metric_cards = _shared.metric_cards

LABELS = {
    'Gene ID': 'ID del gen', 'Gene': 'Gen', 'Mean normalized count': 'Conteo normalizado medio',
    'MLE log₂FC': 'log₂FC de máxima verosimilitud', 'MLE SE': 'Error estándar de MLE',
    'Test statistic': 'Estadístico de prueba', 'Shrunken log₂FC': 'log₂FC regularizado',
    'Shrunken posterior SD': 'DE posterior del efecto regularizado', 'Nominal p': 'p nominal',
    'BH FDR': 'FDR de BH', 'Term BH FDR': 'FDR de BH del término',
    'BH FDR across four contrasts': 'FDR de BH entre cuatro contrastes', 'GSEA estimated FDR': 'FDR estimada de GSEA',
    'Mean difference (signal units)': 'Diferencia media (unidades de señal)', 'Mean published signal': 'Señal publicada media',
    '95% CI lower': 'Límite inferior del IC del 95%', '95% CI upper': 'Límite superior del IC del 95%',
    'Probe ID': 'ID de sonda', 'Sample': 'Muestra', 'Animal label': 'Etiqueta del animal', 'Group': 'Grupo',
    'GEO sample': 'Muestra GEO', 'Source sample title': 'Título de la muestra original', 'PFF injection': 'Inyección de PFF',
    'Sample IDs': 'ID de muestras', 'Starting age': 'Edad inicial', 'Genotype': 'Genotipo', 'Exposure': 'Exposición',
    'Atlas class': 'Clase del atlas', 'Substructure': 'Subestructura', 'Cells': 'Células',
    'Positive-count cells': 'Células con conteo positivo', 'Mean imaging count': 'Conteo de imagen medio',
    'Detection (%)': 'Detección (%)', 'Reference cells (%)': 'Células de referencia (%)',
    'GO ID': 'ID de GO', 'GO term': 'Término de GO', 'Leading-edge genes': 'Genes del núcleo de enriquecimiento',
    'Overlapping genes': 'Genes compartidos', 'Fold enrichment': 'Razón de enriquecimiento',
    'Primary BH FDR': 'FDR de BH principal', 'Same direction (of 7 fits)': 'Misma dirección (de 7 ajustes)',
    'FDR < 0.05 (of 7 fits)': 'FDR < 0.05 (de 7 ajustes)',
    'Lowest omitted-animal log₂FC': 'Menor log₂FC al omitir un animal',
    'Highest omitted-animal log₂FC': 'Mayor log₂FC al omitir un animal',
    'Omitted sample': 'Muestra omitida', 'Omitted animal': 'Animal omitido',
}
_shared.DISPLAY_LABELS = {k: LABELS.get(v, v) for k, v in _shared.DISPLAY_LABELS.items()}
_shared.DISPLAY_LABELS.update({
    'library_size': 'Tamaño de biblioteca', 'detected_genes': 'Genes detectados', 'size_factor': 'Factor de tamaño',
    'mitochondrial_count_fraction': 'Fracción de conteos mitocondriales', 'rRNA_count_fraction': 'Fracción de conteos de ARNr',
    'FA_samples': 'Muestras FA', 'PM_samples': 'Muestras PM', 'genes_testable': 'Genes con prueba válida',
    'genes_fdr_below_0_05': 'Genes con FDR < 0.05', 'primary_hits_still_fdr_below_0_05': 'Hallazgos principales con FDR < 0.05',
    'primary_hits_same_lfc_direction': 'Hallazgos principales con la misma dirección',
    'wald_stat_spearman_with_primary': 'Correlación de Spearman con el Wald principal',
    'contrast': 'Contraste', 'n_FA': 'n FA', 'n_O3': 'n O₃', 'n_control': 'n control', 'n_exposed': 'n expuestos',
    'FDR05_features': 'Características con FDR < 0.05', 'FDR05_higher': 'Aumentadas con FDR < 0.05',
    'FDR05_lower': 'Disminuidas con FDR < 0.05', 'FDR05_low_detection_flag': 'Hallazgos con detección baja',
    'SE_moderated': 'Error estándar moderado', 'effect_units': 'Unidades del efecto', 'dose_mg_per_kg': 'Dosis (mg/kg)',
    'features_tested': 'Características analizadas', 'signal_median': 'Mediana de señal', 'signal_min': 'Señal mínima',
    'features_equal_to_sample_minimum': 'Características en el mínimo de la muestra',
    'present_miRNA_features': 'miRNA detectados', 'absent_miRNA_features': 'miRNA no detectados',
    'low_detection_flag': 'Indicador de detección baja', 'source_features': 'Características de origen',
    'mapped_representatives': 'Representantes con correspondencia', 'DE_FDR_selected_features': 'Características seleccionadas por FDR',
    'ORA_universe': 'Universo ORA', 'ORA_up_genes': 'Genes aumentados de ORA', 'ORA_down_genes': 'Genes disminuidos de ORA',
    'ORA_up_significant_terms': 'Términos ORA significativos: aumento', 'ORA_down_significant_terms': 'Términos ORA significativos: disminución',
    'GSEA_rank_genes': 'Genes ordenados de GSEA', 'GSEA_terms_tested': 'Términos de GSEA analizados',
    'GSEA_terms_algorithmic_FDR_0_05': 'Términos GSEA con FDR algorítmica < 0.05',
    'GSEA_positive_NES_significant': 'Términos significativos con NES positivo',
    'GSEA_negative_NES_significant': 'Términos significativos con NES negativo',
    'GSEA_nonfinite_original_pvalue_genes_excluded': 'Genes excluidos: p original no finito',
    'GSEA_finite_p_missing_padj_genes_retained': 'Genes conservados: p finito y p ajustado ausente',
    'mapped_nonfinite_p_excluded': 'Genes con correspondencia excluidos: p no finito',
    'finite_p_missing_padj_retained': 'Genes conservados: p finito y p ajustado ausente',
    'rank_member_genes': 'Genes del término en la lista', 'foreground_genes': 'Genes seleccionados', 'direction': 'Dirección',
    'dataset': 'Conjunto de datos', 'animals': 'Animales', 'qc_sampled_pairs': 'Pares del control inicial',
    'q30_rate': 'Fracción Q30', 'adapter_trimmed_read_fraction': 'Fracción de lecturas con adaptador recortado',
    'trim_for_quant': 'Recorte para cuantificación', 'raw_pairs': 'Pares originales', 'after_adapter_pairs': 'Pares tras recortar adaptadores',
    'removed_pairs_after_adapter_processing': 'Pares retirados al procesar adaptadores',
    'Raw Q30 (%)': 'Q30 original (%)', 'Raw GC (%)': 'GC original (%)',
    'num_processed_fragments': 'Fragmentos procesados', 'num_mapped_fragments': 'Fragmentos asignados',
    'percent_mapped': 'Asignación (%)', 'library_types': 'Tipos de biblioteca', 'sum_tpm': 'Suma de TPM',
    'quantification_integrity_passed': 'Integridad de cuantificación verificada',
    'total_estimated_fragments': 'Fragmentos estimados totales', 'detected_genes_count_gt_0': 'Genes con conteo > 0',
    'median_normalization_factor': 'Mediana del factor de normalización', 'max_cooks': 'Máxima distancia de Cook',
    'genes_with_nominal_test': 'Genes con prueba nominal', 'FDR05_within_stratum': 'FDR < 0.05 en el estrato',
    'higher_in_O3': 'Aumentados con O₃', 'lower_in_O3': 'Disminuidos con O₃',
    'FDR05_across_four_strata': 'FDR < 0.05 entre cuatro estratos', 'all_zero_in_stratum': 'Todos los conteos cero en el estrato',
    'nonconverged_coefficient_fits': 'Ajustes de coeficientes sin convergencia',
    'nonconverged_shrinkage_fits': 'Ajustes regularizados sin convergencia', 'degrees_of_freedom': 'Grados de libertad',
    'nominal_pvalue': 'p nominal',
})


def fmt_p(value):
    if value is None or pd.isna(value):
        return 'no disponible'
    return f'{value:.2e}' if value < .01 else f'{value:.3f}'


_shared.fmt_p = fmt_p


def show_table(frame, n=20):
    view = frame.copy()
    for column in view.select_dtypes(include=['bool', 'boolean']):
        view[column] = view[column].map(lambda v: 'no disponible' if pd.isna(v) else ('Sí' if v else 'No'))
    if 'direction' in view:
        view['direction'] = view['direction'].replace({'up': 'Aumento', 'down': 'Disminución'})
    _shared.show_table(view, n=n)


def figure(dataset, stem, caption, alt=''):
    root = f'../assets/figures/{dataset}/{stem}'
    display(Markdown('<figure class="research-figure"><a href="' + html.escape(root + '.svg') + '">'
        '<img loading="lazy" src="' + html.escape(root + '.png') + '" alt="' + html.escape(alt or caption) + '"></a>'
        '<figcaption>' + html.escape(caption) + ' <span class="figure-links"><a href="' + html.escape(root + '.svg') +
        '">SVG</a> · <a href="' + html.escape(root + '.pdf') + '">PDF</a></span></figcaption></figure>'))
