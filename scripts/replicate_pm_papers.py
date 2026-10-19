#!/usr/bin/env python3
"""Source-backed reproduction audit for two mouse-brain PM papers.

Reimplements reviewed statistical steps, never executes downloaded author code.
Raw article text/count downloads remain local under private_inputs.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,re,subprocess,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/replication_pm';META=ROOT/'metadata/replication_pm'
PRIVATE=META/'private_inputs';AUTH=PRIVATE/'cortex_authors/mouse/out'
FIG=OUT/'figures';GROUPS=['Control','PM','PFF','PMPFF','PMPFFPM']
COLORS={'Control':'#0072B2','PM':'#D55E00','PFF':'#009E73','PMPFF':'#CC79A7','PMPFFPM':'#A77C00','FA':'#0072B2','HFD':'#009E73'}
LABELS={'Control':'PBS','PM':'PM2.5','PFF':'PFF','PMPFF':'PM–PFF','PMPFFPM':'PM–PFF + PM2.5','FA':'FA','HFD':'HFD'}
FORMATS=['png','svg','pdf']; records={}
CAPTIONS={
 'cortex_replication_benchmark':(
  'Count-based reproduction of the original five-condition cortex model (20 animals; four per condition). A: TMM logCPM PCA of the 500 most variable genes, scaled to unit variance as in the deposited code. B–D: original and recomputed discovery counts and maximum-likelihood log2 fold changes. Orange lines indicate identity; hexbin color encodes log gene density. E–F: sensitivity to using only the eight PM/PBS animals. The contingency includes only genes with finite adjusted values in both fits; full source tables retain missing values. This is computational reproduction, not independent biological replication.',
  'Reproducción basada en conteos del modelo cortical original de cinco condiciones (20 animales; cuatro por condición). A: PCA de TMM logCPM con los 500 genes más variables, escalados a varianza unitaria como en el código depositado. B–D: conteos de descubrimientos y efectos de máxima verosimilitud originales y recalculados. Las líneas naranjas indican identidad; el color hexagonal representa la densidad génica en escala logarítmica. E–F: sensibilidad al usar solo los ocho animales PM/PBS. La contingencia incluye genes con valores ajustados finitos en ambos modelos; las tablas completas conservan los valores ausentes. Es reproducción computacional, no replicación biológica independiente.'),
 'cortex_multicondition_signature':(
  'A: union of the eight genes with smallest author FDR in each condition. Colors show recomputed apeglm-shrunken log2 effects; black dots mark the original Wald BH FDR <0.05. The top Y-linked genes have large MLEs but near-zero shrunken effects; see the separate marker QC diagnostic before interpreting mechanism. B: exact intersections of reproduced FDR sets, including 582 PM-specific genes and 323 exclusive to the PM/PM–PFF pair. C: the eight smallest-author-FDR genes in the deposited GO:0002443 leading edge (leukocyte-mediated immunity). Pointwise intervals are original-model MLE ±1.96 SE; they do not adjust for gene selection. All contrasts share the same four PBS animals. GO membership was verified from the author table; enrichment was not independently recomputed.',
  'A: unión de los ocho genes con menor FDR de los autores en cada condición. Los colores muestran efectos log2 recalculados y regularizados con apeglm; los puntos negros indican FDR de BH de Wald original <0.05. Los genes ligados al cromosoma Y tienen MLE grandes y efectos regularizados próximos a cero; consulte el diagnóstico de marcadores antes de interpretar mecanismos. B: intersecciones exactas de los conjuntos reproducidos, con 582 genes exclusivos de PM y 323 exclusivos del par PM/PM–PFF. C: ocho genes de menor FDR de los autores en el núcleo de enriquecimiento depositado de GO:0002443 (inmunidad mediada por leucocitos). Los intervalos puntuales corresponden a MLE ±1.96 EE del modelo original y no ajustan por selección de genes. Todos los contrastes comparten cuatro animales PBS. La pertenencia a GO se verificó en la tabla de los autores; el enriquecimiento no se recalculó independientemente.'),
 'hypothalamus_method_sensitivity':(
  'Explicit partial methodological reconstruction; the original RUVg controls and number of unwanted factors were not recovered. All 11 animals are included: FA=4, PM=3, HFD=4. Shared source features were restricted to GENCODE M13 protein-coding symbols; CPM >2 in at least three samples, TMM/voom and limma moderated tests were used. The declared reconstruction gate is BH FDR <0.05, |log2FC| >0.8 and AveExpr >1. A: PCA of 500 variable unadjusted logCPM genes, centered without scaling. B: diet and exposure remain separate contrasts; orange denotes passing either contrast gate. C: union of ten smallest-FDR genes per contrast; row z scores describe expression patterns. D: published 17 PM genes, primary PyDESeq2 11, independent R DESeq2 12, and unadjusted limma 0 use different models and gene universes. These differences cannot establish biological nonreplication or the validity of the unavailable original adjustment.',
  'Reconstrucción metodológica parcial explícita; no se recuperaron los controles RUVg ni el número de factores no deseados originales. Incluye los 11 animales: FA=4, PM=3, HFD=4. Se restringieron los genes comunes a símbolos codificantes de GENCODE M13, con CPM >2 en al menos tres muestras, TMM/voom y pruebas moderadas de limma. El umbral declarado es FDR de BH <0.05, |log2FC| >0.8 y AveExpr >1. A: PCA de 500 genes logCPM variables sin ajuste, centrados sin escalar. B: dieta y exposición son contrastes separados; naranja indica superar los umbrales en cualquiera. C: unión de diez genes de menor FDR por contraste; los valores z por fila describen patrones. D: los 17 genes PM publicados, 11 de PyDESeq2 principal, 12 de DESeq2 en R y 0 de limma sin RUVg usan modelos y universos distintos. Estas diferencias no demuestran falta de replicación biológica ni validan el ajuste original no disponible.')}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(x):
 if isinstance(x,dict):return {k:clean(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v)for v in x]
 if isinstance(x,np.generic):return clean(x.item())
 if isinstance(x,float)and not np.isfinite(x):return None
 return x

def write_json(p,d):Path(p).write_text(json.dumps(clean(d),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def table(d,name):d.to_csv(OUT/name,index=False,na_rep='NA',float_format='%.12g')
def export(fig,name,caption,caption_es,sources):
 caption,caption_es=CAPTIONS.get(name,(caption,caption_es))
 for ext in FORMATS:fig.savefig(FIG/f'{name}.{ext}',dpi=300,facecolor='white')
 records[name]={'caption':caption,'caption_es':caption_es,'source_tables':sources,'size_inches':list(fig.get_size_inches()),'PNG_dpi':300,'sha256':{e:sha(FIG/f'{name}.{e}')for e in FORMATS}}
 plt.close(fig)
def panel(ax,label,title):
 ax.set_title(title,loc='left',fontsize=10.5,pad=10)
 ax.text(-.13,1.06,label,transform=ax.transAxes,fontsize=14,fontweight='bold',va='bottom')
 ax.spines[['top','right']].set_visible(False)
def corr(a,b):
 a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float);keep=np.isfinite(a)&np.isfinite(b)
 return float(spearmanr(a[keep],b[keep]).statistic)
def fetch_sources():
 import requests,gzip
 PRIVATE.mkdir(parents=True,exist_ok=True)
 rec=requests.get('https://zenodo.org/api/records/15374632',timeout=60);rec.raise_for_status()
 write_json(PRIVATE/'zenodo_record.json',rec.json())
 f=next(x for x in rec.json()['files']if x['key']=='Fig5 mouse.zip')
 zpath=PRIVATE/'cortex_authors_mouse.zip'
 if not zpath.exists():
  r=requests.get(f['links']['self'],timeout=120);r.raise_for_status();zpath.write_bytes(r.content)
 with zipfile.ZipFile(zpath)as z:
  for name in z.namelist():
   p=Path(name)
   if p.is_absolute()or '..'in p.parts:raise ValueError('Unsafe archive entry')
   if name.endswith(('.R','.r','.csv','.tsv','.txt')):
    dest=PRIVATE/'cortex_authors'/p;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(name))
 supplement=PRIVATE/'hypothalamus_supplement.pdf'
 if not supplement.exists():
  r=requests.get('https://www.jci.org/articles/view/137315/sd/pdf/render/1',timeout=120);r.raise_for_status()
  if not r.content.startswith(b'%PDF-'):raise ValueError('JCI supplement did not return a PDF')
  supplement.write_bytes(r.content)
 url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M13/gencode.vM13.annotation.gtf.gz'
 path=PRIVATE/'gencode.vM13.annotation.gtf.gz'
 if not path.exists():
  r=requests.get(url,timeout=120);r.raise_for_status();path.write_bytes(r.content)
 symbols=set();n=0
 with gzip.open(path,'rt')as f:
  for line in f:
   if line.startswith('#'):continue
   fields=line.split('\t')
   if fields[2]!='gene':continue
   attrs=dict(re.findall(r'(\w+) "([^"]+)"',fields[8]))
   if attrs.get('gene_type')=='protein_coding':symbols.add(attrs['gene_name']);n+=1
 (PRIVATE/'M13_protein_coding_symbols.txt').write_text('\n'.join(sorted(symbols))+'\n')
 write_json(META/'M13_annotation_source.json',{'URL':url,'bytes':path.stat().st_size,'sha256':sha(path),'protein_coding_gene_entries':n,'unique_protein_coding_symbols':len(symbols),'duplicate_gene_names':n-len(symbols)})

def audit_cortex():
 raw={};metrics=[];sets={};comparison_files=[]
 for g in GROUPS[1:]:
  a=pd.read_csv(AUTH/f'deg_{g}_vs_Control_raw.csv');s=pd.read_csv(AUTH/f'deg_{g}_vs_Control_shrunk.csv');b=pd.read_csv(OUT/f'cortex_{g}_vs_Control.csv')
  assert a.ENSEMBL.is_unique and b.ENSEMBL.is_unique and set(a.ENSEMBL)==set(b.ENSEMBL)
  m=a.merge(b,on='ENSEMBL',suffixes=('_author','_refit'),validate='one_to_one')
  m=m.merge(s[['ENSEMBL','log2FoldChange']].rename(columns={'log2FoldChange':'author_log2FoldChange_shrunk'}),on='ENSEMBL',validate='one_to_one')
  sa=set(a.loc[a.padj.lt(.05),'ENSEMBL']);sb=set(b.loc[b.padj.lt(.05),'ENSEMBL']);sets[g]=sb
  m['author_FDR05']=m.padj_author.lt(.05);m['refit_FDR05']=m.padj_refit.lt(.05)
  m['FDR_call_identical']=m.author_FDR05.eq(m.refit_FDR05)
  finite=np.isfinite(m.pvalue_author)&np.isfinite(m.pvalue_refit)
  assert np.array_equal(m.pvalue_author.isna(),m.pvalue_refit.isna())
  assert np.array_equal(m.padj_author.isna(),m.padj_refit.isna())
  rec={'contrast':g,'genes':len(m),'author_FDR05':len(sa),'refit_FDR05':len(sb),'intersection':len(sa&sb),'symmetric_difference':len(sa^sb),'Jaccard':len(sa&sb)/len(sa|sb),'MLE_Spearman':corr(m.log2FoldChange_author,m.log2FoldChange_refit),'max_abs_MLE_difference':np.nanmax(abs(m.log2FoldChange_author-m.log2FoldChange_refit)),'max_abs_nominal_p_difference':np.max(abs(m.loc[finite,'pvalue_author']-m.loc[finite,'pvalue_refit'])),'max_abs_shrunken_difference':np.nanmax(abs(m.author_log2FoldChange_shrunk-m.log2FoldChange_shrunk)),'shrunken_Spearman':corr(m.author_log2FoldChange_shrunk,m.log2FoldChange_shrunk),'same_pvalue_missingness':True,'same_padj_missingness':True}
  metrics.append(rec);raw[g]=m
  file=f'cortex_author_refit_{g}.csv.gz';table(m,file);comparison_files.append(file)
 metric=pd.DataFrame(metrics);table(metric,'cortex_gene_level_metrics.csv')
 # The current report deliberately fits only8animals; it is a separate estimand/model fit.
 current=pd.read_csv(ROOT/'results/cortex/DE_PM_vs_PBS.tsv',sep='\t')
 current=current.rename(columns={'gene_id':'ENSEMBL'})
 pm=raw['PM'].merge(current[['ENSEMBL','padj','log2FoldChange','log2FoldChange_shrunk']],on='ENSEMBL',how='outer',validate='one_to_one')
 table(pm,'cortex_primary8_vs_author20.csv.gz')
 usable=pm.padj_author.notna()&pm.padj.notna();ct=pd.crosstab(pm.loc[usable,'padj'].lt(.05),pm.loc[usable,'padj_author'].lt(.05)).reindex(index=[False,True],columns=[False,True],fill_value=0)
 table(ct.reset_index().rename(columns={'padj':'primary8_FDR05',False:'author20_not_FDR05',True:'author20_FDR05'}),'cortex_subset_vs_author_contingency.csv')
 all_genes=set.union(*sets.values());pattern=Counter(tuple(g in sets[k]for k in GROUPS[1:])for g in all_genes)
 patterns=pd.DataFrame([{'PM':key[0],'PFF':key[1],'PMPFF':key[2],'PMPFFPM':key[3],'genes':n}for key,n in pattern.items()]).sort_values('genes',ascending=False)
 table(patterns,'cortex_four_condition_intersections.csv')
 return raw,metric,sets,pm,ct,patterns

def benchmark_figure(raw,metric,pm,ct):
 fig,axes=plt.subplots(2,3,figsize=(13.8,8.6),layout='constrained')
 p=pd.read_csv(OUT/'cortex_author_pca.csv');info=json.loads((OUT/'cortex_refit_summary.json').read_text())
 ax=axes[0,0]
 for i,g in enumerate(GROUPS):
  q=p[p.group.eq(g)];ax.scatter(q.PC1,q.PC2,c=COLORS[g],marker=['o','^','s','D','P'][i],s=34,label=LABELS[g])
 ax.set(xlabel=f"PC1 ({info['PCA']['PC1']:.1f}%)",ylabel=f"PC2 ({info['PCA']['PC2']:.1f}%)");ax.legend(fontsize=7.6,ncol=2,frameon=False)
 panel(ax,'A','Author-code PCA · all 20 animals')
 ax=axes[0,1];order=['PFF','PM','PMPFF','PMPFFPM'];m=metric.set_index('contrast').loc[order];x=np.arange(4)
 ax.bar(x-.16,m.author_FDR05,width=.32,color='#667C89',label='Author table');ax.bar(x+.16,m.refit_FDR05,width=.32,color='#0D7475',label='Recomputed')
 for j,n in enumerate(m.refit_FDR05):ax.text(j,n+28,f'{n:,}',ha='center',fontsize=8)
 ax.set(xticks=x,xticklabels=[LABELS[g]for g in order],ylabel='Genes with BH FDR <0.05',ylim=(0,1450));ax.tick_params(axis='x',rotation=18);ax.legend(frameon=False,fontsize=8)
 panel(ax,'B','Exact four-contrast discovery counts')
 for ax,key,label in [(axes[0,2],'PM','C'),(axes[1,0],'PMPFF','D')]:
  d=raw[key];x=d.log2FoldChange_author;y=d.log2FoldChange_refit
  ax.hexbin(x,y,gridsize=65,mincnt=1,bins='log',cmap='viridis',rasterized=True);bounds=[min(x.min(),y.min()),max(x.max(),y.max())];ax.plot(bounds,bounds,'--',color='#D55E00',lw=1)
  ax.set(xlabel='Author MLE log₂FC',ylabel='Recomputed MLE log₂FC')
  maxdiff=float(metric.set_index('contrast').loc[key,'max_abs_MLE_difference'])
  ax.text(.03,.96,f'n={len(d):,}\nρ≈1.000000\nmax |Δβ|={maxdiff:.2e}',transform=ax.transAxes,va='top',fontsize=8,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
  panel(ax,label,f'{LABELS[key]} · gene-level agreement')
 ax=axes[1,1];d=pm.dropna(subset=['log2FoldChange_author','log2FoldChange']);x=d.log2FoldChange_author;y=d.log2FoldChange
 ax.hexbin(x,y,gridsize=65,mincnt=1,bins='log',cmap='cividis',rasterized=True);bounds=[min(x.min(),y.min()),max(x.max(),y.max())];ax.plot(bounds,bounds,'--',color='#D55E00',lw=1)
 ax.set(xlabel='Author 20-animal MLE log₂FC',ylabel='Primary 8-animal MLE log₂FC');ax.text(.03,.96,f'n={len(d):,}; ρ={corr(x,y):.4f}',transform=ax.transAxes,va='top',fontsize=8,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
 panel(ax,'E','Changing the analysis population')
 ax=axes[1,2];im=ax.imshow(np.log10(ct.to_numpy()+1),cmap='Blues',vmin=0,interpolation='nearest')
 for i,j in itertools.product(range(2),repeat=2):ax.text(j,i,f'{ct.iloc[i,j]:,}',ha='center',va='center',fontsize=13,color='white'if im.norm(np.log10(ct.iloc[i,j]+1))>.65 else '#173239')
 ax.set(xticks=[0,1],xticklabels=['No','Yes'],yticks=[0,1],yticklabels=['No','Yes'],xlabel='Author 20: FDR <0.05',ylabel='Primary 8: FDR <0.05')
 panel(ax,'F','Common finite-FDR gene universe')
 export(fig,'cortex_replication_benchmark','Independent count-based reproduction of the published five-condition cortex analysis. A: TMMlogCPM PCA from the500 most variable genes, scaled to unit variance exactly as in deposited author code(n=4 per condition). B–D: original and recomputed counts and MLEs; orange lines are identity, not fitted regressions; hexbin color encodes gene density(log counts). E–F: the previous PM/PBS-only8-animal model differs from the original20-animal model; this comparison is sensitivity to model/population, not independent biological replication. The contingency excludes genes missing an adjusted value in either analysis; full tables retain them.','Reproducción independiente basada en conteos del análisis cortical publicado de cinco condiciones. A: PCA de TMMlogCPM de los500genes más variables, escalados a varianza unitaria como en el código depositado(n=4 por condición). B–D: conteos y MLE originales y recalculados; las líneas naranjas son identidad, no regresiones; el color hexagonal representa densidad génica(log conteos). E–F: el modelo previo de8animales PM/PBS difiere del original de20; es sensibilidad al modelo/población, no replicación biológica. La contingencia excluye genes sin valor ajustado en cualquiera de los análisis; las tablas completas los conservan.',['cortex_author_pca.csv','cortex_gene_level_metrics.csv','cortex_author_refit_PM.csv.gz','cortex_author_refit_PMPFF.csv.gz','cortex_primary8_vs_author20.csv.gz','cortex_subset_vs_author_contingency.csv'])

def signature_figure(raw,sets,patterns):
 # Source-selected immune leading edge, not genes chosen to improve a comparison.
 go=pd.read_csv(AUTH/'gsea_PM_vs_Control.csv');immune=go.loc[go.ID.eq('GO:0002443')].iloc[0]
 candidates=set(str(immune.core_enrichment).split('/'))
 pm=raw['PM'].set_index('ENSEMBL');targets=pm.loc[pm.index.isin(candidates)].sort_values('padj_author').head(8)
 forest=targets[['SYMBOL_author','log2FoldChange_refit','lfcSE_refit','padj_refit']].copy();forest['CI95_lower']=forest.log2FoldChange_refit-1.96*forest.lfcSE_refit;forest['CI95_upper']=forest.log2FoldChange_refit+1.96*forest.lfcSE_refit;table(forest.reset_index(),'cortex_author_immune_leading_edge_targets.csv')
 # Union of8top author-ranked genes per condition, independent of changed-method agreement.
 selected=[]
 for g in GROUPS[1:]:
  selected.extend(raw[g].sort_values('padj_author').ENSEMBL.head(8))
 selected=list(dict.fromkeys(selected))
 effects=pd.DataFrame({g:raw[g].set_index('ENSEMBL').loc[selected,'log2FoldChange_shrunk']for g in GROUPS[1:]})
 qvals=pd.DataFrame({g:raw[g].set_index('ENSEMBL').loc[selected,'padj_refit']for g in GROUPS[1:]})
 annotations=raw['PM'].set_index('ENSEMBL').SYMBOL_author.fillna(pd.Series(raw['PM'].set_index('ENSEMBL').index,index=raw['PM'].set_index('ENSEMBL').index))
 effects.insert(0,'symbol',annotations.loc[selected]);table(effects.reset_index(),'cortex_multicondition_effect_heatmap.csv');table(qvals.reset_index(),'cortex_multicondition_effect_FDR.csv')
 fig=plt.figure(figsize=(13.8,9.0),layout='constrained');gs=fig.add_gridspec(2,3,width_ratios=[1.2,1,1.25],height_ratios=[1.1,1])
 ax=fig.add_subplot(gs[:,0]);arr=effects[GROUPS[1:]].to_numpy();lim=float(np.nanmax(abs(arr)));im=ax.imshow(arr,aspect='auto',norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),cmap='RdBu_r',interpolation='nearest')
 ax.set(xticks=np.arange(4),xticklabels=[LABELS[x]for x in GROUPS[1:]],yticks=np.arange(len(selected)),yticklabels=effects.symbol, ylabel='Union of top 8 author-FDR genes per condition');ax.tick_params(axis='x',rotation=35,labelsize=8);ax.tick_params(axis='y',labelsize=7.5)
 for i,j in zip(*np.where(qvals.to_numpy()<.05)):ax.scatter(j,i,s=5,c='black',marker='o')
 fig.colorbar(im,ax=ax,shrink=.7,label='Shrunken log₂FC vs PBS');panel(ax,'A','Shared and distinct expression signatures')
 sub=gs[0,1:].subgridspec(2,1,height_ratios=[2.2,1],hspace=.05);ax=fig.add_subplot(sub[0]);mat=fig.add_subplot(sub[1],sharex=ax)
 p=patterns.copy();x=np.arange(len(p));ax.bar(x,p.genes,color='#0D7475');
 for i,n in enumerate(p.genes):ax.text(i,n+9,str(n),ha='center',fontsize=7)
 ax.set(ylabel='Genes in exact intersection',ylim=(0,p.genes.max()*1.17));ax.tick_params(axis='x',bottom=False,labelbottom=False);panel(ax,'B','Recomputed four-condition intersections')
 for j,row in enumerate(p.itertuples()):
  on=np.array([row.PM,row.PFF,row.PMPFF,row.PMPFFPM]);mat.scatter([j]*4,range(4),s=24,c=['#173239'if v else '#DCE7E9'for v in on]);indices=np.where(on)[0]
  if len(indices)>1:mat.plot([j,j],[indices.min(),indices.max()],c='#173239',lw=1)
 mat.set(yticks=range(4),yticklabels=[LABELS[g]for g in GROUPS[1:]],xticks=[],ylim=(3.5,-.5));mat.spines[['top','right','bottom','left']].set_visible(False)
 ax=fig.add_subplot(gs[1,1:]);y=np.arange(len(forest));ax.errorbar(forest.log2FoldChange_refit,y,xerr=1.96*forest.lfcSE_refit,fmt='o',color='#D55E00',ecolor='#6D7F84',capsize=3);ax.axvline(0,c='#55686D',ls='--',lw=.8)
 ax.set(yticks=y,yticklabels=forest.SYMBOL_author,xlabel='PM2.5 vs PBS MLE log₂FC ±1.96 SE');ax.invert_yaxis();panel(ax,'C','Author immune leading edge · model uncertainty')
 export(fig,'cortex_multicondition_signature','A: union of the8 lowest-author-FDR genes in each condition; colors are recomputed shrunken log2 effects, black dots indicate original Wald BH FDR<0.05, with one symmetric scale. B: exact gene-set intersections from the reproduced four contrasts; PM-only582 and PM/PM–PFF-only323 match the publication. C:8 smallest-author-FDR genes in the deposited leukocyte-mediated-immunity GO:0002443 leading edge; intervals are pointwise original-model Wald95%intervals, not cell-level intervals, and do not account for displayed-gene selection. All contrasts share a4-animal PBS reference.','A: unión de los8genes de menor FDR de los autores en cada condición; colores: efectos log2regularizados recalculados; puntos negros: FDR de BH de Wald<0.05, con escala simétrica común. B: intersecciones génicas exactas de los cuatro contrastes reproducidos;582exclusivos de PM y323exclusivos de PM/PM–PFF coinciden con la publicación. C:8genes de menor FDR de los autores en el núcleo de enriquecimiento depositado de inmunidad mediada por leucocitosGO:0002443; intervalos puntuales de Wald del95%del modelo original, no celulares ni ajustados por selección. Todos los contrastes comparten4controles PBS.',['cortex_multicondition_effect_heatmap.csv','cortex_multicondition_effect_FDR.csv','cortex_four_condition_intersections.csv','cortex_author_immune_leading_edge_targets.csv'])

def hypothesis_figure():
 pm=pd.read_csv(OUT/'hypothalamus_PM_vs_FA_unadjusted_voom.csv');hfd=pd.read_csv(OUT/'hypothalamus_HFD_vs_FA_unadjusted_voom.csv');d=pm.merge(hfd,on='gene',suffixes=('_PM','_HFD'),validate='one_to_one')
 groups=pd.read_csv(OUT/'hypothalamus_author_like_design.csv');lcpm=pd.read_csv(OUT/'hypothalamus_original_like_logCPM.csv').set_index('gene');sel=list(dict.fromkeys([*pm.sort_values('adj.P.Val').gene.head(10),*hfd.sort_values('adj.P.Val').gene.head(10)]))
 order=groups.sort_values('group',key=lambda s:s.map({'FA':0,'PM':1,'HFD':2}))['sample'].tolist();values=lcpm.loc[sel,order];zs=values.sub(values.mean(axis=1),axis=0).div(values.std(axis=1,ddof=1),axis=0)
 table(zs.reset_index(),'hypothalamus_method_sensitivity_heatmap.csv');table(d,'hypothalamus_PM_HFD_effect_comparison.csv.gz')
 info=json.loads((OUT/'hypothalamus_refit_summary.json').read_text());pca_data=lcpm.loc[lcpm.var(axis=1).nlargest(500).index].T;center=pca_data-pca_data.mean(axis=0);u,s,vt=np.linalg.svd(center.to_numpy(),full_matrices=False);pcs=u[:,:2]*s[:2];variance=s[:2]**2/(s**2).sum()*100;pca=pd.DataFrame(pcs,index=lcpm.columns,columns=['PC1','PC2']).join(groups.set_index('sample'));table(pca.reset_index(),'hypothalamus_unadjusted_pca.csv')
 fig,axes=plt.subplots(2,2,figsize=(12.4,9.6),layout='constrained')
 ax=axes[0,0]
 for k,g in enumerate(['FA','PM','HFD']):
  sub=pca[pca.group.eq(g)];ax.scatter(sub.PC1,sub.PC2,color=COLORS[g],marker=['o','^','s'][k],s=45,label=f'{g} (n={len(sub)})')
  offsets={'PM1':(25,0),'PM2':(25,16),'PM3':(8,23),'FA3':(25,-23),'FA4':(25,-8),'HFD2':(-10,17),'HFD3':(-10,-15)}
  for sample,row in sub.iterrows():
   delta=offsets.get(row.label,(5,5));ax.annotate(row.label,(row.PC1,row.PC2),xytext=delta,textcoords='offset points',fontsize=7.5,arrowprops={'arrowstyle':'-','color':'#9AA8AB','lw':.6}if row.label in offsets else None)
 ax.margins(x=.1,y=.12);ax.set(xlabel=f'PC1 ({variance[0]:.1f}%)',ylabel=f'PC2 ({variance[1]:.1f}%)');ax.legend(frameon=False,fontsize=8);panel(ax,'A','Unadjusted 11-sample reconstruction')
 ax=axes[0,1];ax.scatter(d.logFC_HFD,d.logFC_PM,s=5,c='#899CA1',alpha=.3,rasterized=True);sig=d.original_like_threshold_HFD.fillna(False)|d.original_like_threshold_PM.fillna(False);ax.scatter(d.loc[sig,'logFC_HFD'],d.loc[sig,'logFC_PM'],s=6,c='#D55E00',alpha=.35,rasterized=True)
 ax.axvline(0,c='#7F8C8D',lw=.7);ax.axhline(0,c='#7F8C8D',lw=.7);ax.set(xlabel='HFD vs FA · limma log₂FC',ylabel='PM2.5 vs FA · limma log₂FC');ax.text(.03,.96,f'{len(d):,} coding genes\nOriginal RUVg not recovered',transform=ax.transAxes,va='top',fontsize=8,bbox={'facecolor':'white','alpha':.9,'edgecolor':'none'});panel(ax,'B','Exposure and diet remain distinct contrasts')
 ax=axes[1,0];lim=float(abs(zs).max().max());im=ax.imshow(zs,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),interpolation='nearest');lab=groups.set_index('sample').loc[order,'label'];ax.set(xticks=np.arange(11),xticklabels=lab,yticks=np.arange(len(sel)),yticklabels=sel);ax.tick_params(axis='x',rotation=45,labelsize=8);ax.tick_params(axis='y',labelsize=7.5)
 for i,sample in enumerate(order):ax.plot([i-.49,i+.49],[-.8,-.8],lw=5,color=COLORS[groups.set_index('sample').loc[sample,'group']],clip_on=False)
 fig.colorbar(im,ax=ax,label='Within-gene z score of unadjusted logCPM',shrink=.7);panel(ax,'C','Source-preserving differential patterns')
 ax=axes[1,1];labels=['Published\nRUVg-based','Primary\nPyDESeq2','Independent\nRDESeq2','Unadjusted\nlimma+gates'];numbers=[17,11,12,info['contrasts']['PM_vs_FA']['FDR05_abslogFC_gt08_AveExpr_gt1']];x=np.arange(4);ax.bar(x,numbers,color=['#667C89','#0D7475','#4EA3A4','#D55E00']);ax.set(xticks=x,xticklabels=labels,ylabel='PM-associated genes',ylim=(0,22));ax.tick_params(axis='x',labelsize=8)
 for i,n in enumerate(numbers):ax.text(i,n+.5,str(n),ha='center',fontsize=11)
 ax.text(.02,.95,'Different models/filter families\nNOT a matched replication comparison',transform=ax.transAxes,va='top',fontsize=8,color='#A6431E');panel(ax,'D','Method dependence; exact RUVg unavailable')
 export(fig,'hypothalamus_method_sensitivity','Explicit partial reconstruction of the JCI hypothalamus experiment, not recovery of its unavailable RUVg factors. Common features of all11source files are restricted to GENCODEM13 protein-coding symbols, CPM>2 in>=3samples, TMM/voom, limma~condition and moderated tests; the displayed original-like gate is BH FDR<0.05,|log2FC|>0.8,AveExpr>1. These declared reconstruction choices do not establish the authors exact settings. A:PCA of500 variable unadjustedlogCPMgenes, centered without variance scaling; no sample exclusion. B:two within-model contrasts; orange marks either contrast passing declared gates. C:union of10lowest-FDRgenes percontrast, all11animals, row-zscores are descriptive. D:published17PMgenes vs11primary/12independentDESeq2/0unadjustedlimma; differing gene universes and methods prevent interpreting the differences as biological nonreplication. The original negative-control IDs and RUVgk are unavailable.','Reconstrucción parcial explícita del experimento hipotalámico de JCI, sin recuperar sus factores RUVgno disponibles. Se usan características comunes de los11archivos, símbolos codificantes GENCODEM13, CPM>2en>=3muestras, TMM/voom, limma~condición y pruebas moderadas; umbral reconstruido: FDRBH<0.05,|log2FC|>0.8,AveExpr>1. Estas decisiones declaradas no establecen la configuración exacta de los autores. A:PCAde500geneslogCPMvariables sin ajuste, centrados sin escalar; ningún animal excluido. B:contrastes del mismo modelo; naranjaindica cumplir umbrales en cualquiera. C:unión de10genes de menorFDRporcontraste,11animales;zscoreporfila descriptivo. D:17genesPMpublicados frentea11delanálisisprincipal/12deDESeq2enR/0delimmasinRUVg; universos y métodos diferentes impiden concluir no replicación biológica. No se recuperaron IDsdecontrolesnegativos niRUVgk.',['hypothalamus_original_like_logCPM.csv','hypothalamus_unadjusted_pca.csv','hypothalamus_PM_HFD_effect_comparison.csv.gz','hypothalamus_method_sensitivity_heatmap.csv','hypothalamus_refit_summary.json'])

def marker_diagnostic():
 source=pd.read_csv(PRIVATE/'cortex_authors/mouse/data/mouse_gene_count.txt',sep='\t')
 design=pd.read_csv(OUT/'cortex_author_design.csv');samples=design['sample'].tolist()
 marker_names=['Ddx3y','Uty','Kdm5d','Xist']
 source=source.loc[source.gene_name.isin(marker_names)].set_index('gene_name').loc[marker_names]
 counts=source[samples].T;counts.index.name='sample'
 normalized=pd.read_csv(OUT/'cortex_all20_normalized_counts.csv').set_index('gene_id')
 norm=normalized.loc[source.gene_id,samples].copy();norm.index=marker_names
 observed=counts.join(design.set_index('sample'));observed['Y_marker_sum_raw']=observed[marker_names[:3]].sum(axis=1)
 table(observed.reset_index(),'cortex_marker_expression_by_animal.csv')
 lognorm=np.log2(norm+1);table(lognorm.rename_axis('gene').reset_index(),'cortex_marker_log2_normalized_counts.csv')
 fig,axes=plt.subplots(1,2,figsize=(13.4,5.0),layout='constrained',gridspec_kw={'width_ratios':[1.5,1]})
 ax=axes[0];im=ax.imshow(lognorm,aspect='auto',cmap='viridis',vmin=0,interpolation='nearest')
 labels=[f'{LABELS[g]} {i+1}'for g in ['Control','PFF','PM','PMPFF','PMPFFPM']for i in range(4)]
 ax.set(xticks=np.arange(20),xticklabels=labels,yticks=np.arange(4),yticklabels=marker_names);ax.tick_params(axis='x',rotation=75,labelsize=7)
 for i,row in design.iterrows():ax.plot([i-.49,i+.49],[-.7,-.7],lw=5,color=COLORS[row.group],clip_on=False)
 for x in [3.5,7.5,11.5,15.5]:ax.axvline(x,c='white',lw=1)
 fig.colorbar(im,ax=ax,shrink=.8,label='log₂(DESeq2-normalized count +1)');panel(ax,'A','Observed sex-linked marker heterogeneity')
 ax=axes[1]
 for i,g in enumerate(GROUPS):
  sample=design.loc[design.group.eq(g),'sample'];ax.scatter(np.log2(norm.loc['Xist',sample]+1),np.log2(norm.loc[marker_names[:3],sample].sum(axis=0)+1),c=COLORS[g],marker=['o','^','s','D','P'][i],s=42,label=LABELS[g],alpha=.85)
 ax.set(xlabel='log₂(normalized Xist +1)',ylabel='log₂(normalized sum of three Y markers +1)');ax.legend(frameon=False,fontsize=8,loc='upper right');panel(ax,'B','Descriptive QC; no sex assigned from expression')
 caption='Observed source counts for Ddx3y, Uty, Kdm5d and Xist in every cortex animal. A: normalized expression with condition labels; no animals are omitted. B: Xist versus the sum of three Y-linked markers, without classification or an exposure-effect test. All four PBS animals have zero raw counts for each Y marker, while other conditions contain heterogeneous values. These genes also rank at the top of the original Wald FDR tables but shrink toward zero with apeglm. This is a diagnostic flag for independently verifying sample sex, identity and other technical explanations. Sex cannot be assigned or a confounder established from this plot alone. The original group-only model is preserved.'
 caption_es='Conteos originales observados de Ddx3y, Uty, Kdm5d y Xist en todos los animales de corteza. A: expresión normalizada con etiquetas de condición; no se omiten animales. B: Xist frente a la suma de tres marcadores ligados al cromosoma Y, sin clasificación ni prueba de efecto de exposición. Los cuatro animales PBS tienen cero conteos originales en cada marcador Y; las otras condiciones muestran heterogeneidad. Estos genes están entre los primeros de las tablas FDR de Wald originales, pero apeglm regulariza sus efectos hacia cero. Es una señal diagnóstica para verificar sexo, identidad y explicaciones técnicas con metadatos independientes. La figura por sí sola no permite asignar sexo ni demostrar confusión. Se conserva el modelo original ajustado solo por grupo.'
 export(fig,'cortex_marker_qc',caption,caption_es,['cortex_marker_expression_by_animal.csv','cortex_marker_log2_normalized_counts.csv','cortex_author_design.csv'])
 return {'markers':marker_names,'PBS_raw_Y_counts_all_zero':bool(observed.loc[observed.group.eq('Control'),marker_names[:3]].eq(0).all().all()),'sample_sex_independently_verified':False,'qualification':caption}

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--fetch-sources',action='store_true');parser.add_argument('--refit',action='store_true');args=parser.parse_args()
 OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(exist_ok=True);META.mkdir(parents=True,exist_ok=True)
 if args.fetch_sources:fetch_sources()
 if args.refit:subprocess.run(['Rscript',str(ROOT/'scripts/replicate_pm_papers.R'),'all'],cwd=ROOT,check=True)
 plt.rcParams.update({'font.size':9.5,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
 raw,metrics,sets,pm,ct,patterns=audit_cortex();benchmark_figure(raw,metrics,pm,ct);signature_figure(raw,sets,patterns);hypothesis_figure();marker_qc=marker_diagnostic()
 original=pd.read_csv(AUTH/'DEG_summary_counts.csv');refit=json.loads((OUT/'cortex_refit_summary.json').read_text());hyp=json.loads((OUT/'hypothalamus_refit_summary.json').read_text())
 a=pd.read_csv(PRIVATE/'cortex_authors/mouse/data/mouse_gene_count.txt',sep='\t');b=pd.read_csv(ROOT/'data/GSE297423/processed/GSE297423_mouse_gene_count.txt.gz',sep='\t');assert a.equals(b)
 paper_counts={'PFF':109,'PM':1205,'PMPFF':1222,'PMPFFPM':639};computed={g:refit['contrasts'][g]['FDR05']for g in paper_counts};assert paper_counts==computed
 intersections={'PM_specific':len(sets['PM']-set.union(sets['PFF'],sets['PMPFF'],sets['PMPFFPM'])),'PM_PMPFF_exclusive_pair':len((sets['PM']&sets['PMPFF'])-(sets['PFF']|sets['PMPFFPM']))};assert intersections=={'PM_specific':582,'PM_PMPFF_exclusive_pair':323}
 claims=[]
 for g,n in paper_counts.items():claims.append({'study':'GSE297423','claim':f'{g} vs PBS: DE genes','reported':n,'computed':computed[g],'status':'computationally_reproduced','source':'Science main PDF page 8; Zenodo DEG_summary_counts.csv','qualification':'Same deposited counts and design; not independent biological replication'})
 for k,n in intersections.items():claims.append({'study':'GSE297423','claim':k,'reported':582 if k=='PM_specific' else 323,'computed':n,'status':'computationally_reproduced','source':'Science main PDF page 8; author four-contrast tables','qualification':'Exact exclusive intersection after excluding unavailable FDR'})
 for label,n,v in [('PCA_PC1_percent',29.2,refit['PCA']['PC1']),('PCA_PC2_percent',9.8,refit['PCA']['PC2'])]:claims.append({'study':'GSE297423','claim':label,'reported':n,'computed':v,'status':'reproduced_at_published_rounding','source':'Zenodo mouse/out/pca.pdf','qualification':'Followed code TMM/scaled PCA; prose instead describes filterByExpr'})
 for g,n in [('PFF',13),('PM',365),('PMPFF',266),('PMPFFPM',271)]:
  x=pd.read_csv(AUTH/f'gsea_{g}_vs_Control.csv');assert len(x)==n and x['p.adjust'].lt(.05).all();claims.append({'study':'GSE297423','claim':f'{g}: enriched GO terms','reported':n,'computed':len(x),'status':'author_table_verified_not_independently_recomputed','source':'Science main PDF page 8; Zenodo GSEA tables','qualification':'Original OrgDb/clusterProfiler/fgsea versions and seed unavailable; enrichment not recomputed'})
 for key,n in [('PM_vs_FA',17),('HFD_vs_FA',761)]:claims.append({'study':'GSE145840','claim':f'{key}: DE genes','reported':n,'computed':hyp['contrasts'][key]['FDR05_abslogFC_gt08_AveExpr_gt1'],'status':'partial_method_sensitivity_not_exact_replication','source':'JCI supplement PDF page 10 panels C/D; methods pages 19–20','qualification':'Negative-control IDs and RUVg k not recovered; declared abundance/FC reconstruction; 91 HFD-only features excluded without zero imputation'})
 table(pd.DataFrame(claims),'paper_claim_audit.csv')
 sourcefiles=[PRIVATE/'zenodo_record.json',PRIVATE/'cortex_authors_mouse.zip',PRIVATE/'cortex_authors/mouse/mouse_deg.R',PRIVATE/'cortex_authors/mouse/data/mouse_gene_count.txt',PRIVATE/'hypothalamus_supplement.pdf',PRIVATE/'gencode.vM13.annotation.gtf.gz',ROOT/'data/GSE297423/processed/GSE297423_mouse_gene_count.txt.gz']
 sourcefiles.extend(sorted((ROOT/'data/GSE145840/processed').glob('GSM*_Hypothalamus_*.txt.gz')))
 write_json(META/'source_manifest.json',{'manifest_created_UTC':datetime.now(timezone.utc).isoformat(),'sources':[
  {'source':'https://doi.org/10.1126/science.adu4132','role':'Original cortex paper; main PDF page 8 results and pages 13–14 methods'},
  {'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297423','role':'Public 20-animal cortex count deposit; exact equality with author matrix verified'},
  {'source':'https://doi.org/10.5281/zenodo.15374632','resolved_record':'15374633','download':'https://zenodo.org/api/records/15374633/files/Fig5%20mouse.zip/content','file':'Fig5 mouse.zip','role':'Official author code, count matrix, all four DE tables, PCA and GO exports'},
  {'source':'https://doi.org/10.1172/JCI137315','role':'Original hypothalamus paper'},
  {'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145840','role':'Public hypothalamus source counts and animal metadata'},
  {'source':'https://www.jci.org/articles/view/137315/sd/pdf/render/1','role':'Original supplementary methods pages 19–20 and Figure S9 page 10'},
  {'source':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M13/gencode.vM13.annotation.gtf.gz','role':'Frozen GENCODE M13 protein-coding annotation'}],
  'files':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)}for p in sourcefiles],
  'private_source_inputs':'Publisher text and original source archives stay local under ignored private_inputs. Public artifacts retain curated facts, URLs and hashes.'})
 summary={'created_at_UTC':datetime.now(timezone.utc).isoformat(),'cortex':{'classification':'Successful count-based statistical reproduction, not independent biological replication or raw-read realignment','paper_counts':paper_counts,'refit_counts':computed,'all_four_FDR_sets_identical':bool(metrics.symmetric_difference.eq(0).all()),'gene_level_comparison':metrics.to_dict('records'),'count_matrix_identical_to_GEO':True,'genes_retained':21388,'biological_samples':20,'group_size':4,'PCA':refit['PCA'],'paper_intersections':intersections,'author_GO_tables_verified':{'PFF':13,'PM':365,'PMPFF':266,'PMPFFPM':271},'GSEA_recomputed':False},'hypothalamus':{'classification':'Partial methodological sensitivity only; exact original RUVg fit not recoverable','reported_PM':17,'reported_HFD':761,'partial_refit_PM':hyp['contrasts']['PM_vs_FA']['FDR05_abslogFC_gt08_AveExpr_gt1'],'partial_refit_HFD':hyp['contrasts']['HFD_vs_FA']['FDR05_abslogFC_gt08_AveExpr_gt1'],'primary_report_PM':11,'M13_coding_input':hyp['M13_protein_coding_input'],'tested_genes':hyp['genes_after_abundance_filter'],'shared_source_features':hyp['input_genes'],'excluded_HFD_only_features':91,'not_exact_reason':['Original5000negativecontrolIDsnotprovided; RUVglatent-factorcountknotstated','Supplementdescribes limma-basededgeR and negativebinomialtesting withoutcompleteRcommands','Supplementcontrolselection wording pairs not significantlychanged withFDR<.05, leavingintentambiguous','Originalgene-levelDEtable unavailable; onlypublishedcountbenchmark recovered','MainFig2logCPM>1 differs fromsupplementCPM>2; reconstructionexplicitlydocumentschosenfilter'],'publication_caption_issue':'SupplementS9captionCcalls itmuscle; panelCvisiblylabeledHypothalamus(N11), with12PMspecific+5sharedand756HFDspecific+5shared'},'no_primary_analysis_overwritten':True,'figure_families':list(records)}
 summary['cortex']['marker_QC']=marker_qc
 summary['hypothalamus']['not_exact_reason']=[
  'The original 5000 empirical negative-control gene IDs were not recovered; the RUVg factor count k was not stated.',
  'The supplement describes limma/edgeR and negative-binomial testing without complete reproducible R commands.',
  'Negative-control selection text pairs not significantly changed genes with FDR <0.05, leaving the intended rule ambiguous.',
  'No original gene-level DE table was recovered; only published count benchmarks could be read.',
  'Main Figure 2 logCPM >1 differs from supplement CPM >2; this sensitivity explicitly declares its own abundance rule.']
 summary['hypothalamus']['publication_caption_issue']='Supplement Figure S9 caption C calls the tissue muscle, but panel C is visibly labeled Hypothalamus (N=11): PM 12 specific +5 shared; HFD 756 specific +5 shared.'
 write_json(OUT/'summary.json',summary);write_json(OUT/'figure_manifest.json',records)
 methods={'cortex':'ReviewedreimplementationoftheauthorsZenodoRcode:all20samples~group;>=10countsin>=4samples;DESeq2defaultnormalization/dispersions/Wald/Cooks; results defaultalpha=.1 followedFDR<.05;apeglm;TMMlogCPMtop500genesPCAwithscaleTRUE. Originalscriptwasread,notexecutedwholesale.','cortex_prose_code_discrepancy':'PaperPCAtextfilterByExpr;depositedcodeactuallyuses>=10in4.Paperpresentsselected8sampleexposurecomparisonbutoriginaldispersionfitincludesall20libraries.','hypothalamus':'All11groupsincludedwithconditiondesign;50417sharedfeatures(nozeroimputation),GENCODEM13proteincodingsymbols;CPM>2in>=3samples;TMM/voomlimma~condition;PMvsFAandHFDvsFA;moderatedBHtests,abslogFC>.8,AveExpr>1. NoRUVg/batchcovariatesinvented. This is an explicitly declared sensitivity,not a recovered originalfit.','multiplicity':'CortexBHwithincontrastafterDESeq2independentfilteringatdefaultalpha.1. HypothalamusBHwithincontrastover13398reconstructiongenes. Target-forestintervalspointwiseWaldMLE±1.96SE,selectionnotadjusted.','source_identifiers':'ExactEnsemblIDsforcortexagreement;symbolsmaydifferversuscurrentannotation. Sourcecurrentprimary8sampleanalysiskeptunchanged.','software':refit['software'],'script_sha256':sha(Path(__file__)),'R_script_sha256':sha(ROOT/'scripts/replicate_pm_papers.R')}
 methods.update({
  'cortex':'Reviewed reimplementation of deposited Zenodo R code: all 20 samples, design ~group; retain >=10 counts in >=4 samples; DESeq2 default normalization, dispersion, Wald and Cook filtering; results() default alpha=0.1 for independent filtering, followed by BH FDR <0.05 for discoveries; apeglm shrinkage. PCA uses TMM logCPM, the 500 most variable genes and scale=TRUE. Downloaded source code was reviewed and never executed wholesale.',
  'cortex_prose_code_discrepancy':'PCA prose states filterByExpr; deposited code uses >=10 counts in >=4 samples. The original DESeq2 dispersion fit includes all 20 libraries across five groups.',
  'hypothalamus':'All 11 animals; separate FA, PM and HFD conditions; exact intersection of 50,417 source features with no zero imputation; GENCODE M13 protein-coding symbols; declared CPM >2 in >=3 samples; TMM/voom and limma ~condition moderated tests; PM versus FA and HFD versus FA. Gate: BH FDR <0.05, |log2FC| >0.8 and AveExpr >1. No RUVg or batch factors were invented. This is declared methodological sensitivity, not recovery of the original fit.',
  'multiplicity':'Cortex: BH within each contrast after default alpha=0.1 independent filtering. Hypothalamus: BH within each contrast across 13,398 tested genes. Forest intervals are pointwise Wald MLE +/-1.96 SE and do not adjust for displayed-gene selection.',
  'source_identifiers':'Exact Ensembl IDs establish cortex agreement; displayed symbols come from the deposit. Existing primary eight-sample cortex and seven-sample hypothalamus analyses remain unchanged.',
  'marker_QC':'Descriptive raw and normalized Ddx3y/Uty/Kdm5d/Xist expression in all 20 animals. Independent sex metadata were not recovered. No sex covariate, sample exclusion or mechanistic claim was inferred from the diagnostic.',
  'software':{**refit['software'],**hyp['software'],'Python':__import__('sys').version.split()[0],'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__}})
 write_json(OUT/'methods.json',methods)
 write_json(OUT/'validation.json',{'passed':True,'checks':{'author_and_GEO_counts_identical':True,'all_four_gene_universes_identical':True,'all_four_missing_p_and_padj_masks_identical':True,'all_four_significant_gene_sets_identical':bool(metrics.symmetric_difference.eq(0).all()),'all_four_count_benchmarks_reproduced':True,'paper_exclusive_intersections_reproduced':True,'hypothalamus_source_harmonization_no_zero_imputation':True,'all_figures_exported_in_three_formats':all((FIG/f'{name}.{ext}').exists()for name in records for ext in FORMATS)},'scope':'Independent numerical agreement and export audit. Figure visual review is recorded separately after inspection. No biological replication claim.'})
 print(json.dumps({'cortex_counts':computed,'exact_FDR_sets':summary['cortex']['all_four_FDR_sets_identical'],'hypothalamus_partial_PM':summary['hypothalamus']['partial_refit_PM'],'figures':list(records)},indent=2))

if __name__=='__main__':main()
