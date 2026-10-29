#!/usr/bin/env python3
"""Create presentation-size figures from frozen analytical tables; never refit models.

Outputs are 1600x900 PNG and editable-axis SVG/PDF, in English and Spanish.
All filtering here is display availability (finite coordinates) and is recorded.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import textwrap
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/defense_slides'
FIG=OUT/'figures'
FIG.mkdir(parents=True,exist_ok=True)
COL={'teal':'#176B75','orange':'#BB561D','navy':'#344B75','grey':'#667780','light':'#C6CDD1','red':'#9E3245','ink':'#17333C','grid':'#E3E9EC'}
ORDER=['Control','PM','PFF','PMPFF','PMPFFPM']
GROUPLAB={'Control':'PBS','PM':'PM','PFF':'PFF','PMPFF':'PM-PFF','PMPFFPM':'PM-PFF + PM'}
GROUPCOL=dict(zip(ORDER,[COL['grey'],COL['teal'],COL['orange'],COL['navy'],COL['red']]))
manifest={}
checks={}


def load(path,**kwargs):
 p=ROOT/path
 if p.suffix=='.json': return json.loads(p.read_text())
 return pd.read_csv(p,sep='\t' if p.suffix=='.tsv' else ',',**kwargs)


def clean_text(text):
 if isinstance(text,(list,tuple)):
  return [clean_text(x) for x in text]
 return re.sub(r'\b(All|all|the|The|of|from|among|includes|shows|counts|global|selected|with|los|las|Los|Las|de|entre|muestra|incluye|del|Jaccard)(?=\d)',r'\1 ',text)


def tr(en,es,lang): return clean_text(es if lang=='es' else en)


def base(title,subtitle,lang,ncols=1,width_ratios=None):
 fig,axs=plt.subplots(1,ncols,figsize=(16,6.5),dpi=100,
                       gridspec_kw={'width_ratios':width_ratios} if width_ratios else None)
 fig.subplots_adjust(left=.09,right=.965,bottom=.17,top=.88,wspace=.37)
 # The slide supplies the headline; preserve only data-bearing axes and legends.
 # Full scope and sample context remain in the HTML caption and speaker notes.
 for ax in np.atleast_1d(axs):
  ax.spines[['right','top']].set_visible(False)
  ax.spines[['left','bottom']].set_color(COL['grey'])
  ax.tick_params(labelsize=18)
  ax.grid(axis='y',color=COL['grid'],lw=.7,zorder=0)
  ax.set_axisbelow(True)
 return fig,np.atleast_1d(axs)


def footer(fig,text):
 lines=textwrap.wrap(text,width=105,break_long_words=False,break_on_hyphens=False)
 # Full methodological caption is retained in the manifest, not duplicated inside the image.
 return


def save(fig,stem,lang,cap_en,cap_es,sources,scope,diagnostics):
 name=stem+('_es' if lang=='es' else '')
 cap_en=clean_text(cap_en);cap_es=clean_text(cap_es)
 fig.canvas.draw()
 renderer=fig.canvas.get_renderer()
 clipped=[]
 for artist in fig.findobj(Text):
  if artist.get_visible() and artist.get_text():
   bb=artist.get_window_extent(renderer)
   if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fig.bbox.width+1 or bb.y1 > fig.bbox.height+1:
    clipped.append({'text':artist.get_text(),'bounds':[float(z) for z in bb.bounds]})
 diagnostics['text_outside_canvas']=clipped
 exports={}
 for ext in ('png','svg','pdf'):
  p=FIG/f'{name}.{ext}'
  fig.savefig(p,dpi=100,facecolor='white',transparent=False)
  exports[ext]={'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 plt.close(fig)
 manifest[name]={'language':lang,'caption_en':cap_en,'caption_es':cap_es,
  'caption':cap_es if lang=='es' else cap_en,'files':{e:v['path'] for e,v in exports.items()},
  'exports':exports,'source_paths':sources,'source_hashes':{s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
  'render_size_pixels':[1600,650],'figure_size_inches':[16,6.5],'PNG_dpi':100,'scope':scope,'diagnostics':diagnostics}


def cortex_reproduction(lang):
 summary=load('results/replication_pm/summary.json')['cortex']
 d=load('results/replication_pm/cortex_author_refit_PM.csv.gz')
 valid=np.isfinite(d[['log2FoldChange_author','log2FoldChange_refit']]).all(axis=1)
 q=d.loc[valid]
 assert len(q)==21388
 fig,(a,b)=base(tr('Cortical author results are computationally reproduced','Los resultados corticales se reproducen computacionalmente',lang),
  tr('Same deposited counts and 20 animals across five conditions','Mismos conteos depositados y 20 animales de cinco condiciones',lang),lang,2,[1,1.12])
 groups=['PFF','PM','PMPFF','PMPFFPM'];x=np.arange(4);w=.34
 author=[summary['paper_counts'][g] for g in groups];refit=[summary['refit_counts'][g] for g in groups]
 a.bar(x-w/2,author,w,color='white',edgecolor=COL['grey'],hatch='///',label=tr('Author','Autores',lang),zorder=3)
 a.bar(x+w/2,refit,w,color=COL['teal'],label=tr('Refit','Reajuste',lang),zorder=3)
 for i,n in enumerate(author): a.text(i,n+35,f'{n:,}',ha='center',fontsize=21,fontweight='bold',color=COL['ink'])
 a.set(xticks=x,xticklabels=['PFF','PM','PM-PFF','PM-PFF\n+ PM'],ylim=(0,1500),ylabel=tr('Genes with BH q < 0.05','Genes con BH q < 0.05',lang))
 a.set_title(tr('A  Identical significant-gene sets','A  Conjuntos significativos idénticos',lang),loc='left',fontsize=21,pad=18)
 a.legend(frameon=False,fontsize=17,loc='upper left',ncol=2)
 hb=b.hexbin(q.log2FoldChange_author,q.log2FoldChange_refit,gridsize=52,bins='log',mincnt=1,cmap='cividis',rasterized=True,zorder=3)
 lim=(-5,25);b.plot(lim,lim,color=COL['orange'],ls='--',lw=2,zorder=4)
 b.set(xlim=lim,ylim=lim,xlabel=tr('Author MLE log₂FC','MLE log₂FC de autores',lang),ylabel=tr('Refit MLE log₂FC','MLE log₂FC recalculado',lang))
 b.set_title(tr('B  PM versus PBS effect agreement','B  Efectos PM frente a PBS',lang),loc='left',fontsize=21,pad=18)
 b.text(.04,.96,'n = 21,388\nρ = 1.000000\nmax |Δβ| = 7.93 × 10⁻⁵',transform=b.transAxes,fontsize=18,va='top',bbox={'fc':'white','ec':'none','alpha':.95})
 cax=b.inset_axes([.91,.07,.025,.24]);cb=fig.colorbar(hb,cax=cax)
 cb.set_ticks([1,10**int(np.floor(np.log10(hb.get_array().max())))])
 cb.ax.tick_params(labelsize=15);cb.ax.set_title(tr('Genes/bin','Genes/casilla',lang),fontsize=12,pad=6)
 cb.ax.text(.5,-.20,'log₁₀',transform=cb.ax.transAxes,ha='center',fontsize=13,color=COL['grey'])
 footer(fig,tr('All four FDR sets have Jaccard = 1. This is same-source computational reproduction, not independent biological replication. Dense hexagons include every finite effect; no percentile clipping.','Los cuatro conjuntos FDR tienen Jaccard = 1. Es reproducción computacional de la misma fuente, sin replicación biológica independiente. Los hexágonos incluyen todos los efectos finitos, sin recorte por percentiles.',lang))
 save(fig,'rep_cortex_reproduction',lang,
  'Count-based reproduction of the original five-condition cortex analysis, using 20 animals (four per condition). A: author and refit BH q<0.05 counts are identical for all four contrasts, and each exact Ensembl-ID discovery set has Jaccard1. B: all 21,388 finite PM-versus-PBS maximum-likelihood effects, with identity line and log-count hexagonal binning. No raw-read realignment or independent biological cohort is involved.',
  'Reproducción basada en conteos del análisis cortical original de cinco condiciones, con 20 animales (cuatro por condición). A: los conteos BH q<0.05 de autores y reajuste coinciden en los cuatro contrastes y cada conjunto de IDs Ensembl exactos tiene Jaccard1. B: los 21,388 efectos de máxima verosimilitud PM frente a PBS finitos, con línea de identidad y casillas hexagonales en escala log del número de genes. No incluye realineamiento de lecturas ni una cohorte biológica independiente.',
  ['results/replication_pm/summary.json','results/replication_pm/cortex_author_refit_PM.csv.gz'],
  'Presentation re-expression of frozen computational reproduction outputs',{'finite_effects':len(q),'excluded_nonfinite_effects':int((~valid).sum()),'rho':float(spearmanr(q.log2FoldChange_author,q.log2FoldChange_refit).statistic),'all_four_FDR_sets_identical':summary['all_four_FDR_sets_identical']})


def cortex_markers(lang):
 d=load('results/replication_pm/cortex_marker_expression_by_animal.csv')
 assert len(d)==20 and len(d[d.group=='Control'])==4
 fig,(a,b)=base(tr('Cortical marker counts flag a source-quality question','Los marcadores corticales señalan una duda sobre la fuente',lang),
  tr('All 20 animals retained; raw counts describe a diagnostic, not assigned sex','Se conservan los 20 animales; los conteos originales son un diagnóstico, sin asignar sexo',lang),lang,2,[1.12,1])
 for i,g in enumerate(ORDER):
  z=d[d.group==g];offset=np.linspace(-.19,.19,len(z));y=np.log2(1+z.Y_marker_sum_raw.to_numpy())
  a.scatter(i+offset,y,c=GROUPCOL[g],s=115,edgecolor='white',lw=.7,zorder=4)
  b.scatter(np.log2(1+z.Xist.to_numpy()),np.log2(1+z.Y_marker_sum_raw.to_numpy()),s=115,color=GROUPCOL[g],marker=['o','s','^','D','P'][i],edgecolor='white',lw=.7,zorder=4,label=GROUPLAB[g])
 a.set(xticks=np.arange(5),xticklabels=['PBS','PM','PFF','PM-PFF','PM-PFF\n+ PM'],ylim=(-.7,16),ylabel=tr('log₂(1 + raw Y-marker sum)\nDdx3y + Uty + Kdm5d','log₂(1 + suma Y original)\nDdx3y + Uty + Kdm5d',lang))
 a.set_title(tr('A  Three Y-linked markers','A  Marcadores del cromosoma Y',lang),fontsize=21,loc='left',pad=18)
 b.set(xlabel='log₂(1 + Xist)',ylabel=tr('log₂(1 + raw Y-marker sum)','log₂(1 + suma Y original)',lang),ylim=(-.7,16))
 b.set_title(tr('B  Xist and Y-linked counts','B  Xist y suma de marcadores Y',lang),fontsize=21,loc='left',pad=18)
 b.legend(frameon=False,fontsize=17,loc='upper right')
 footer(fig,tr('All four PBS animals have zero raw counts for all three Y markers. Sex metadata were not recovered. These observations do not assign sample sex, prove confounding or justify sample exclusion.','Los cuatro animales PBS tienen cero conteos originales de los tres marcadores Y. No se recuperaron metadatos de sexo. Estas observaciones no asignan sexo, prueban confusión ni justifican excluir muestras.',lang))
 save(fig,'rep_cortex_y_markers',lang,
  'Descriptive source-count diagnostic for all 20 cortex animals. A: log₂(1+sum of raw Ddx3y,Uty,Kdm5d counts), with fixed within-condition horizontal offsets solely to distinguish animals. B: raw Xist and the same Y-marker sum on log₂(1+count) axes. Shapes and colors distinguish conditions. PBS samples have zero raw counts for each Y marker. No sex classification, exposure-effect test, sample exclusion or covariate was inferred.',
  'Diagnóstico descriptivo de conteos originales de los 20 animales corticales. A: log₂(1+suma de conteos originales Ddx3y,Uty,Kdm5d), con desplazamientos horizontales fijos dentro de cada condición solo para distinguir animales. B: Xist y la misma suma Y en ejes log₂(1+conteo). Formas y colores distinguen condiciones. PBS tiene cero conteos originales en cada marcador Y. No se infiere clasificación de sexo, prueba de exposición, exclusión ni covariable.',
  ['results/replication_pm/cortex_marker_expression_by_animal.csv','results/replication_pm/summary.json'],
  'Descriptive original-count QC; log₂(1+count) transform and deterministic display-only offsets',{'animals':len(d),'PBS_animals':4,'PBS_raw_Y_counts_all_zero':bool((d.loc[d.group=='Control',['Ddx3y','Uty','Kdm5d']]==0).all().all())})


def hypothalamus_sensitivity(lang):
 r=load('results/replication_pm/summary.json')['hypothalamus'];ind=load('results/hypothalamus/independent_deseq2_check/summary.json')
 fig,(a,b)=base(tr('Hypothalamic discoveries depend on the specified analysis','Los hallazgos hipotalámicos dependen del análisis especificado',lang),
  tr('Original RUVg controls and factor count were not recovered','No se recuperaron los controles ni el número de factores RUVg originales',lang),lang,2,[1.35,1])
 names=tr(['Published\nRUVg pipeline\nn=11','Primary\nPyDESeq2\nn=7','Software check\nR DESeq2\nn=7','Declared\nTMM/voom\nn=11'],['Publicado\nflujo RUVg\nn=11','Principal\nPyDESeq2\nn=7','Comprobación\nR DESeq2\nn=7','Declarado\nTMM/voom\nn=11'],lang)
 vals=[r['reported_PM'],r['primary_report_PM'],ind['native_FDR_hits']['R_independent_filter_FALSE'],r['partial_refit_PM']]
 bars=a.bar(np.arange(4),vals,color=[COL['grey'],COL['teal'],COL['navy'],COL['orange']],zorder=3,width=.62)
 for i,n in enumerate(vals): a.text(i,n+.5,str(n),ha='center',fontsize=24,fontweight='bold',color=COL['ink'])
 a.set(xticks=np.arange(4),xticklabels=names,ylim=(0,22),ylabel=tr('Selected genes','Genes seleccionados',lang))
 a.tick_params(axis='x',labelsize=16)
 a.set_title(tr('A  PM versus filtered air','A  PM frente a aire filtrado',lang),fontsize=21,loc='left',pad=18)
 vv=[r['reported_HFD'],r['partial_refit_HFD']];b.bar([0,1],vv,color=[COL['grey'],COL['orange']],width=.55,zorder=3)
 for i,n in enumerate(vv): b.text(i,n+80,f'{n:,}',ha='center',fontsize=24,fontweight='bold',color=COL['ink'])
 b.set(xticks=[0,1],xticklabels=tr(['Published\nRUVg pipeline','Declared\nTMM/voom'],['Publicado\nflujo RUVg','Declarado\nTMM/voom'],lang),ylim=(0,4300),ylabel=tr('Selected genes','Genes seleccionados',lang))
 b.set_title(tr('B  High-fat diet versus filtered air','B  Dieta grasa frente a aire filtrado',lang),fontsize=21,loc='left',pad=18)
 footer(fig,tr('Models, gene universes and abundance/effect gates differ. The declared n=11 sensitivity uses BH<.05, |log₂FC|>.8 and AveExpr>1. Count differences do not establish biological nonreplication.','Difieren modelos, universos génicos y umbrales de abundancia/efecto. La sensibilidad n=11 usa BH<.05, |log₂FC|>.8 y AveExpr>1. Las diferencias no establecen falta de replicación biológica.',lang))
 save(fig,'rep_hypothalamus_sensitivity',lang,
  'Methodological sensitivity of hypothalamus discovery counts. A contrasts the published17 PM-associated genes, primary7-animal PyDESeq2 result11, independent-software R DESeq2 result12, and explicitly unadjusted11-animal TMM/voom result0. B contrasts published761 HFD genes with declared sensitivity3600. Sample sets, models and gene universes differ. The original RUVg control genes and k were not recovered, so these are not matched replications of the original adjustment.',
  'Sensibilidad metodológica de los conteos de hallazgos hipotalámicos. A contrasta17 genes PM publicados,11 del análisis principal PyDESeq2 de 7 animales,12 de DESeq2 en R y0 del TMM/voom explícitamente sin ajuste RUVg de 11 animales. B contrasta761 genes HFD publicados con3600 de la sensibilidad declarada. Difieren muestras, modelos y universos. No se recuperaron controles RUVg ni k, por lo que no son reproducciones equiparables del ajuste original.',
  ['results/replication_pm/summary.json','results/replication_pm/methods.json','results/hypothalamus/independent_deseq2_check/summary.json'],
  'Method comparison across different declared analyses, not effect or causal comparison',{'PM_counts':[int(v) for v in vals],'HFD_counts':[int(v) for v in vv]})


def mirna574(lang):
 s=load('results/mirna/analysis_signal.tsv').set_index('probe_id').loc['mmu-miR-574-5p']
 d=load('results/replication_mirna_ozone/mirna_author_targets.csv')
 q=d[(d.miRNA=='mmu-miR-574-5p')&(d.contrast=='PM5_vs_Control')].iloc[0]
 fig,(a,b)=base(tr('miR-574-5p direction agrees; array-wide evidence remains limited','miR-574-5p concuerda en dirección, con evidencia global limitada',lang),
  tr('Nine hippocampal arrays; deposited normalized signal, without an assumed log₂ scale','Nueve microarreglos hipocampales; señal normalizada depositada, sin suponer escala log₂',lang),lang,2,[1,1.15])
 for i in range(3):
  y=s.iloc[i*3:(i+1)*3].to_numpy(float);offset=np.array([-.13,0,.13])
  a.scatter(i+offset,y,s=130,color=[COL['grey'],COL['teal'],COL['orange']][i],marker=['o','s','^'][i],edgecolor='white',lw=.8,zorder=4)
  mean=y.mean();se=y.std(ddof=1)/np.sqrt(len(y));a.errorbar(i,mean,yerr=se,fmt='_',markersize=30,lw=2.2,color=COL['ink'],capsize=6,zorder=5)
 a.set(xticks=[0,1,2],xticklabels=tr(['Saline','PM 1 mg/kg','PM 5 mg/kg'],['Salina','PM 1 mg/kg','PM 5 mg/kg'],lang),ylabel=tr('Deposited normalized signal','Señal normalizada depositada',lang),xlim=(-.5,2.5))
 a.set_title(tr('A  Every animal, mean ± SE','A  Todos los animales, media ± EE',lang),fontsize=21,loc='left',pad=18)
 for i,(key,c,m) in enumerate([('limma',COL['teal'],'o'),('LSD',COL['orange'],'s')]):
  eff=q.mean_difference;lo=q[f'{key}_CI95_lower'];hi=q[f'{key}_CI95_upper']
  b.errorbar(eff,1-i,xerr=[[eff-lo],[hi-eff]],fmt=m,markersize=11,c=c,lw=3,capsize=6,zorder=4)
  b.text(.99,.84-i*.28,f"p = {q[f'{key}_pvalue']:.4f}\nBH q = {q[f'{key}_BH_FDR']:.3f}",transform=b.transAxes,ha='right',va='top',fontsize=20,color=c)
 b.axvline(0,color=COL['grey'],ls='--',lw=1.8)
 b.set(yticks=[1,0],yticklabels=['limma','LSD'],ylim=(-.6,1.6),xlim=(-.57,.40),xlabel=tr('High dose − saline, signal-unit difference','Dosis alta − salina, diferencia de señal',lang))
 b.set_title(tr('B  Pointwise 95% intervals','B  Intervalos puntuales del 95%',lang),fontsize=21,loc='left',pad=18)
 footer(fig,tr('ANOVA p=.0788: the protected-LSD gate fails. Nominal pairwise evidence does not survive BH across 1,245 features. Means±SE and confidence intervals are not multiplicity-adjusted. qPCR/rescue assays are separate.','ANOVA p=.0788: falla el requisito de LSD protegida. La evidencia nominal por pares no supera BH entre 1,245 características. Media±EE e intervalos no ajustan por multiplicidad. qPCR/rescate son ensayos distintos.',lang))
 save(fig,'rep_mirna574',lang,
  'Prespecified author target mmu-miR-574-5p in all 9 source arrays, three per group. A shows every deposited normalized signal plus group mean±SE. No second log transformation is applied. B shows high-dose minus saline mean difference with limma and unmoderated LSD pointwise 95% intervals, nominal p and separate BH q across 1,245 features. Omnibus ANOVA p=.078773 fails the protected-LSD gate. Direction and nominal pairwise agreement do not establish array-wide FDR confirmation or reproduce independent qPCR/rescue assays.',
  'Diana original preespecificada mmu-miR-574-5p en los 9 microarreglos, tres por grupo. A muestra todas las señales normalizadas depositadas y media±EE. No se aplica un segundo logaritmo. B muestra diferencia dosis alta menos salina con intervalos puntuales95% de limma y LSD no moderada, p nominal y BH q separado entre 1,245 características. ANOVA p=.078773 no cumple el requisito de LSD protegida. La concordancia direccional/nominal no establece confirmación FDR global ni reproduce qPCR/rescate independientes.',
  ['results/mirna/analysis_signal.tsv','results/replication_mirna_ozone/mirna_author_targets.csv','results/replication_mirna_ozone/summary.json'],
  'Display of measured normalized array signals and existing target statistics',{'animals':9,'mean_difference':float(q.mean_difference),'ANOVA_pvalue':float(q.ANOVA_pvalue),'limma_BH_FDR':float(q.limma_BH_FDR),'protected_LSD_nominal':bool(q.protected_LSD_nominal)})


def ozone_concordance(lang):
 x=load('results/replication_mirna_ozone/ozone_concordance_genotype_pooled.csv.gz');y=load('results/replication_mirna_ozone/ozone_concordance_old_ozone_genotype.csv.gz')
 x=x[np.isfinite(x[['log2FoldChange_author','log2FoldChange_ours']]).all(axis=1)]
 y=y[np.isfinite(y[['author_log2FC','log2FoldChange']]).all(axis=1)]
 assert len(x)==15559 and len(y)==858
 fig,(a,b)=base(tr('Ozone-study effect agreement depends on the comparison','La concordancia en el estudio de O₃ depende de la comparación',lang),
  tr('Same-source animals; Salmon reanalysis differs from original HISAT2/featureCounts','Mismos animales de origen; Salmon difiere de HISAT2/featureCounts original',lang),lang,2)
 hb=a.hexbin(x.log2FoldChange_author,x.log2FoldChange_ours,gridsize=48,bins='log',mincnt=1,cmap='cividis',rasterized=True,zorder=3)
 cax=a.inset_axes([.91,.07,.025,.24]);cb=fig.colorbar(hb,cax=cax)
 cb.set_ticks([1,10**int(np.floor(np.log10(hb.get_array().max())))])
 cb.ax.tick_params(labelsize=15);cb.ax.set_title(tr('Genes/bin','Genes/casilla',lang),fontsize=12,pad=6)
 cb.ax.text(.5,-.20,'log₁₀',transform=cb.ax.transAxes,ha='center',fontsize=13,color=COL['grey'])
 b.scatter(y.author_log2FC,y.log2FoldChange,s=15,c=COL['teal'],alpha=.65,edgecolor='none',rasterized=True,zorder=3)
 for ax in (a,b):
  ax.plot([-9.3,9.3],[-9.3,9.3],c=COL['orange'],ls='--',lw=2,zorder=4)
  ax.set(xlim=(-9.3,9.3),ylim=(-9.3,9.3),xlabel=tr('Author MLE log₂FC, E4 / E3','MLE log₂FC original, E4 / E3',lang),ylabel=tr('Candidate MLE log₂FC','MLE log₂FC candidato',lang))
 a.set_title(tr('A  All common genotype effects','A  Efectos de genotipo comunes',lang),fontsize=21,loc='left',pad=18)
 b.set_title(tr('B  Author-selected 17M O₃ list','B  Lista original seleccionada 17M O₃',lang),fontsize=21,loc='left',pad=18)
 a.text(.04,.96,'n = 15,559\nρ = 0.865\n87.9% '+tr('same direction','misma dirección',lang),transform=a.transAxes,fontsize=18,va='top',bbox={'fc':'white','ec':'none','alpha':.95})
 b.text(.04,.96,'858 / 959 '+tr('matched','emparejados',lang)+'\nρ = 0.962\n98.8% '+tr('same direction','misma dirección',lang),transform=b.transAxes,fontsize=18,va='top',bbox={'fc':'white','ec':'none','alpha':.95})
 footer(fig,tr('All finite matched effects are retained and both panels share full-range axes. PanelB is conditional on a published selected list; its correlation is not an unbiased whole-transcriptome replication estimate.','Se conservan todos los efectos finitos emparejados y ambos paneles comparten ejes de rango completo. B condiciona a una lista publicada seleccionada; su correlación no estima reproducción transcriptómica sin sesgo.',lang))
 save(fig,'rep_ozone_concordance',lang,
  'Ozone-study same-source effect concordance. A includes all 15,559 common finite exact-symbol effects from the24-animal genotype-only candidate versus author genotype CSV; hexbin color encodes log gene count. B includes858 of959 author-selected17-month ozone E4-versus-E3 S3 genes with finite exact-symbol effects. All points are retained, with shared full-range axes and identity lines. The98.8% direction agreement andρ=.962 in B are conditional on source selection and do not constitute independent biological replication.',
  'Concordancia de efectos en la misma fuente del estudio O₃. A incluye los 15,559 efectos finitos con símbolos exactos comunes entre candidato de genotipo de 24 animales y CSV original; el color hexagonal codifica log del número de genes. B incluye858 de959 genes S3 seleccionados por autores para E4 frente a E3 bajo O₃ a17 meses. Se retienen todos los puntos, con ejes completos compartidos y líneas de identidad. El98.8% de concordancia direccional yρ=.962 de B condicionan a la selección original y no constituyen replicación biológica independiente.',
  ['results/replication_mirna_ozone/ozone_concordance_genotype_pooled.csv.gz','results/replication_mirna_ozone/ozone_concordance_old_ozone_genotype.csv.gz','results/replication_mirna_ozone/summary.json'],
  'Existing exact-symbol effect comparison; no refit, no selective cropping or omission',{'whole_common_finite_effects':len(x),'S3_common_finite_effects':len(y),'S3_source_selected':959})


def hypothalamus_primary(lang):
 d=load('results/hypothalamus/differential_expression_mle.tsv')
 d=d[d.padj<.05].sort_values('log2FoldChange')
 assert len(d)==11
 fig,(a,)=base(tr('Eleven primary hypothalamic genes pass BH q < 0.05','Once genes hipotalámicos principales superan BH q < 0.05',lang),
  tr('PM versus filtered air in chow-fed mice, n=3 versus 4 animals','PM frente a aire filtrado en dieta estándar, n=3 frente a 4 animales',lang),lang)
 fig.subplots_adjust(left=.20,right=.83)
 yy=np.arange(len(d))
 for i,(_,row) in enumerate(d.iterrows()):
  c=COL['teal'] if row.log2FoldChange<0 else COL['orange']
  a.errorbar(row.log2FoldChange,i,xerr=1.96*row.lfcSE,fmt='o',c=c,ms=10,lw=2.3,capsize=5,zorder=3)
  a.text(1.035,i,f'{row.padj:.4f}',transform=a.get_yaxis_transform(),va='center',fontsize=19,color=COL['ink'])
 a.axvline(0,c=COL['grey'],ls='--',lw=1.6)
 a.set(yticks=yy,yticklabels=d.gene.tolist(),xlabel=tr('MLE log₂(PM / filtered air), pointwise 95% Wald interval','MLE log₂(PM / aire filtrado), intervalo Wald puntual 95%',lang),ylim=(-.75,len(d)-.25))
 a.tick_params(axis='y',labelsize=21)
 a.text(1.035,1.035,'BH q',transform=a.transAxes,fontsize=20,fontweight='bold',color=COL['ink'])
 a.grid(axis='y',visible=False);a.grid(axis='x',color=COL['grid'],lw=.8)
 footer(fig,tr('All 11 primary BH discoveries are displayed. Intervals are MLE ±1.96 SE and do not adjust for gene selection or multiplicity. These are regional bulk effects; no cell-specific or spatial exposure effect is inferred.','Se muestran los 11 hallazgos BH principales. Los intervalos son MLE ±1.96 EE y no ajustan por selección génica ni multiplicidad. Son efectos regionales bulk, sin inferir respuesta celular o espacial de exposición.',lang))
 save(fig,'rep_hypothalamus_primary',lang,
  'All 11 primary hypothalamic BH q<.05 discoveries from the 7-animal PM-versus-filtered-air count analysis, ordered by unshrunken log₂FC. Points are maximum-likelihood coefficients, bars are pointwise±1.96SE Wald intervals, and the separate q column reports BH inference. The displayed intervals do not adjust for selection or multiplicity. No refit, gene omission, or spatial/cell-specific exposure inference is performed.',
  'Los11 hallazgos hipotalámicos principales BH q<.05 del análisis de conteos de 7 animales PM frente a aire filtrado, ordenados por log₂FC sin regularizar. Los puntos son coeficientes de máxima verosimilitud, las barras intervalos Wald puntuales±1.96EE y la columna q informa BH por separado. Los intervalos no ajustan por selección ni multiplicidad. No se reajustan modelos, omiten genes ni infiere respuesta espacial/celular de exposición.',
  ['results/hypothalamus/differential_expression_mle.tsv','results/hypothalamus/summary.json'],
  'Presentation of all preexisting BH discoveries, with pointwise unadjusted Wald intervals',{'discoveries':len(d),'source_effect':'MLE','interval':'pointwise MLE ±1.96SE','full_source_genes':19611})


def cortex_primary(lang):
 d=load('results/cortex/DE_PM_vs_PBS_MLE.tsv');ann=load('results/cortex/gene_annotation.tsv')[['gene_id','gene']]
 d=d.merge(ann,on='gene_id',how='left',validate='1:1')
 valid=np.isfinite(d[['log2FoldChange','pvalue']]).all(axis=1);q=d.loc[valid].copy()
 q['minuslog10p']=-np.log10(np.clip(q.pvalue,1e-300,1))
 sig=q.padj<.05;down=sig&(q.log2FoldChange<0);up=sig&(q.log2FoldChange>=0)
 top=d.sort_values(['padj','pvalue']).head(7).sort_values('log2FoldChange')
 fig,(a,b)=base(tr('Primary cortical contrast yields 1,167 BH discoveries','El contraste cortical principal produce 1,167 hallazgos BH',lang),
  tr('Eight PM/PBS animals only; this differs from the 20-animal author model','Solo ocho animales PM/PBS; difiere del modelo original de 20 animales',lang),lang,2,[1.2,1])
 for mask,c,lab in [(~sig,COL['light'],tr('Other finite tests','Otras pruebas finitas',lang)),(down,COL['teal'],tr('Lower in PM:509','Menor en PM:509',lang)),(up,COL['orange'],tr('Higher in PM:658','Mayor en PM:658',lang))]:
  a.scatter(q.loc[mask,'log2FoldChange'],q.loc[mask,'minuslog10p'],s=9,c=c,alpha=.65,edgecolor='none',rasterized=True,zorder=2 if c==COL['light'] else 3,label=lab)
 a.axvline(0,c=COL['grey'],ls='--',lw=1)
 a.set(xlabel=tr('MLE log₂(PM / PBS)','MLE log₂(PM / PBS)',lang),ylabel=tr('−log₁₀(nominal Wald p)','−log₁₀(p Wald nominal)',lang))
 a.set_title(tr('A  Every finite primary Wald test','A  Todas las pruebas Wald finitas',lang),fontsize=21,loc='left',pad=18)
 a.legend(frameon=False,fontsize=16,loc='upper right')
 yy=np.arange(len(top))
 for i,(_,row) in enumerate(top.iterrows()):
  b.errorbar(row.log2FoldChange,i,xerr=1.96*row.lfcSE,fmt='o',color=COL['teal'] if row.log2FoldChange<0 else COL['orange'],lw=2.1,capsize=5,ms=9,zorder=3)
 b.set(yticks=yy,yticklabels=top.gene.fillna(top.gene_id).tolist(),xlabel=tr('MLE log₂FC, pointwise 95% CI','MLE log₂FC, IC puntual 95%',lang),ylim=(-.6,len(top)-.4))
 b.axvline(0,c=COL['grey'],ls='--',lw=1.4)
 b.set_xticks([-.5,0,.5,1.0])
 b.set_title(tr('B  Seven smallest-BH-q genes','B  Siete genes con menor BH q',lang),fontsize=21,loc='left',pad=18)
 b.grid(axis='y',visible=False)
 footer(fig,tr('All19,225 finite Wald tests are shown without percentile clipping; BH discovery status sets color, nominal p sets height. Top-gene intervals do not adjust for selection. Sex/sample-identity metadata require independent verification.','Se muestran las 19,225 pruebas Wald finitas sin recorte por percentiles; BH define color y p nominal define altura. Los intervalos de genes destacados no ajustan por selección. Sexo/identidad requieren verificación independiente.',lang))
 save(fig,'rep_cortex_primary',lang,
  'Primary8-animal PM-versus-PBS cortical analysis. A displays all 19,225 finite maximum-likelihood effects and nominal Wald p values, with color indicating BH q<.05 status (658 higher,509 lower). Unavailable tests remain in source tables; no display percentile clipping is applied. B labels the 7 smallest-BH-q genes and shows their pointwise unshrunken MLE±1.96SE intervals, without selection/multiplicity adjustment. This primary contrast uses a different fitted sample set than the 20-animal author reproduction.',
  'Análisis cortical principal PM frente a PBS de 8 animales. A muestra los 19,225 efectos de máxima verosimilitud y p Wald nominales finitos, con color para BH q<.05 (658 mayores,509 menores). Las pruebas no disponibles permanecen en las tablas y no se recorta por percentiles. B identifica los 7 genes con menor BH q y sus intervalos puntuales MLE±1.96EE sin ajuste por selección/multiplicidad. El contraste usa un conjunto de muestras distinto de la reproducción original de 20 animales.',
  ['results/cortex/DE_PM_vs_PBS_MLE.tsv','results/cortex/gene_annotation.tsv','results/cortex/summary.json'],
  'Existing primary differential table, full finite p-value domain and transparent selected-gene display',{'finite_Wald_tests':len(q),'unavailable_Wald_tests':int((~valid).sum()),'BH05':int(sig.sum()),'top_genes':top.gene.fillna(top.gene_id).tolist(),'p_floor_for_display':1e-300,'number_pvalues_below_floor':int((q.pvalue<1e-300).sum())})


def ozone_primary(lang):
 d=load('results/ozone/summary.json');keys=['M3_E3_O3_vs_FA','M3_E4_O3_vs_FA','M17_E3_O3_vs_FA','M17_E4_O3_vs_FA']
 values=[d['contrasts'][k]['fdr_005'] for k in keys]
 fig,(a,b)=base(tr('Ozone exposure evidence varies across the planned contrasts','La evidencia de O₃ varía entre los contrastes planificados',lang),
  tr('All 24 libraries, age × APOE × exposure model; three animals per source group','Las 24 bibliotecas, modelo edad × APOE × exposición; tres animales por grupo',lang),lang,2,[1.25,1])
 x=np.arange(4);cols=[COL['grey'],COL['teal'],COL['grey'],COL['navy']]
 a.bar(x,values,color=cols,width=.6,zorder=3)
 for i,n in enumerate(values):a.text(i,n+10,str(n),ha='center',fontsize=24,fontweight='bold',color=COL['ink'])
 a.set(xticks=x,xticklabels=['3M / E3','3M / E4','17M / E3','17M / E4'],ylim=(0,470),ylabel=tr('Genes with BH q <0.05','Genes con BH q <0.05',lang))
 a.tick_params(axis='x',labelsize=17)
 a.set_title(tr('A  Within-stratum O₃ versus air','A  O₃ frente a aire por estrato',lang),fontsize=21,loc='left',pad=18)
 b.bar([0],[d['global_LRT']['fdr_005']],color=COL['orange'],width=.30,zorder=3)
 b.text(0,25.5,str(d['global_LRT']['fdr_005']),ha='center',fontsize=30,fontweight='bold',color=COL['ink'])
 b.set(xticks=[0],xticklabels=tr(['Global exposure LRT\n4 df'],['LRT global de exposición\n4 gl'],lang),xlim=(-.45,1.1),ylim=(0,32),ylabel=tr('Genes with BH q <0.05','Genes con BH q <0.05',lang))
 b.set_title(tr('B  Joint exposure-related terms','B  Términos conjuntos de exposición',lang),fontsize=21,loc='left',pad=18)
 b.text(.98,.70,tr('18,511 valid tests\nSeparate BH family\nNo signed effect','18,511 pruebas válidas\nFamilia BH separada\nSin efecto con signo',lang),transform=b.transAxes,ha='right',va='top',fontsize=18,color=COL['ink'])
 footer(fig,tr('Panel A and B have different y-axis scales and distinct hypotheses. Discovery-count differences do not test response interaction. Ages 3/17M are exposure-onset labels; BH is separate within each planned family.','A y B usan escalas verticales e hipótesis distintas. Diferencias de hallazgos no prueban interacción. Edades 3/17M indican inicio de exposición; BH se aplica por separado dentro de cada familia planificada.',lang))
 save(fig,'rep_ozone_primary',lang,
  'Primary ozone count-model results from all24 original paired-end libraries. A shows within-age/APOE O3-versus-filtered-air BH q<.05 discovery counts0,399,2,3, each contrast adjusted separately. B shows24 global4df exposure-associated likelihood-ratio-test discoveries among18,511 valid tests, a separate BH family. Panels use different y-axis scales. Count differences between strata do not establish statistical interaction or relative susceptibility. No single signed effect is assigned to the global test.',
  'Resultados principales del modelo de conteos O₃ de las 24 bibliotecas pareadas originales. A muestra0,399,2,3 hallazgos BH q<.05 O₃ frente a aire filtrado por edad/APOE, con ajuste separado. B muestra24 hallazgos de la prueba global de exposición de4gl entre18,511 pruebas válidas, familia BH distinta. Los paneles usan escalas verticales diferentes. Diferencias de conteos no establecen interacción ni susceptibilidad relativa. La prueba global no recibe un efecto único con signo.',
  ['results/ozone/summary.json','results/ozone/methods.json'],
  'Existing primary BH discovery summaries, not a new interaction analysis',{'within_stratum_BH05':values,'global_LRT_BH05':d['global_LRT']['fdr_005'],'global_df':4,'all_libraries':24})


def main():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':20,'axes.labelsize':20,'axes.titlesize':22,
  'xtick.labelsize':18,'ytick.labelsize':18,'legend.fontsize':18,'svg.fonttype':'none','pdf.fonttype':42,
  'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':COL['ink'],'text.color':COL['ink'],
  'xtick.color':COL['ink'],'ytick.color':COL['ink'],'savefig.facecolor':'white'})
 builders=[cortex_reproduction,cortex_markers,hypothalamus_sensitivity,mirna574,ozone_concordance,hypothalamus_primary,cortex_primary,ozone_primary]
 for lang in ('en','es'):
  for builder in builders:
   builder(lang)
   print(lang,builder.__name__,flush=True)
 payload={'purpose':'Presentation-size restyling of frozen saved outputs; no statistical refits',
  'script':'scripts/build_defense_replication_figures.py','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'styles':{'pixels':[1600,900],'font':'DejaVu Sans','axis_labels_pt':20,'tick_labels_pt':18,'max_panels':2,'svg_text_editable':True},
  'figures':manifest}
 (OUT/'replication_figure_manifest.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n')
 print('Wrote',len(manifest),'bilingual figures /',len(manifest)*3,'exports')

if __name__=='__main__':main()
