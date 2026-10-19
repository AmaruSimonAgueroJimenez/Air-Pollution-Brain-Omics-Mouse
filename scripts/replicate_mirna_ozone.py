#!/usr/bin/env python3
"""Transparent paper-result reproduction candidates for Ku2017 and Nakamya2025.

Uses existing source arrays / completed Salmon count models. It never overwrites
primary analyses. R candidate fits are an explicit separate command/stage.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, re, subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.stats import f_oneway, f as f_dist, t as t_dist, norm, false_discovery_control, spearmanr, pearsonr
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/replication_mirna_ozone';META=ROOT/'metadata/replication_mirna_ozone'
FIG=OUT/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,
 'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
COLORS={'Control':'#287a86','PM1':'#c08022','PM5':'#b24747'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(value):
 if isinstance(value,dict):return{str(k):clean(v)for k,v in value.items()}
 if isinstance(value,(list,tuple)):return[clean(x)for x in value]
 if isinstance(value,np.generic):return clean(value.item())
 if isinstance(value,float)and not np.isfinite(value):return None
 return value

def save_json(path,d):path.write_text(json.dumps(clean(d),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def bh(values):
 p=np.asarray(values,float);q=np.full(p.shape,np.nan);ok=np.isfinite(p)
 if ok.any():q[ok]=false_discovery_control(p[ok],method='bh')
 return q

def save_fig(fig,stem):
 for suffix in ['png','svg','pdf']:fig.savefig(FIG/f'{stem}.{suffix}',dpi=250,facecolor='white')
 plt.close(fig)

def author_sources():
 source=ROOT/'data/E-MTAB-14917/supplements/PMC11942628_supplementary.zip'
 from zipfile import ZipFile
 import io
 with ZipFile(source)as z:
  nested=z.read('ijms-26-02407-s001.zip')
 with ZipFile(io.BytesIO(nested))as z:
  target=META/'private_inputs/ozone_supplements/tables';target.mkdir(parents=True,exist_ok=True)
  for name in z.namelist():
   if name.endswith('.pdf'):
    p=target/Path(name).name;p.write_bytes(z.read(name))
    p.with_suffix('.txt').write_text('\n'.join(x.extract_text()for x in PdfReader(p).pages))
 # S3 has two printed columns of numbered genes. Whitespace within symbols is a PDF extraction artifact.
 p=next(target.glob('Supplementary Table S3*.pdf'))
 text='\n'.join(x.extract_text()for x in PdfReader(p).pages)
 pat=re.compile(r'(?:^|\s)(\d{1,3})\s+([A-Za-z0-9][A-Za-z0-9_. -]*?)\s+(-?\d+\.\d+)\s+(\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)')
 rows=[]
 for page_i,page in enumerate(PdfReader(p).pages,1):
  for line in page.extract_text().splitlines():
   for match in pat.finditer(line):
    row,symbol,lfc,pvalue=match.groups()
    rows.append({'source_row':int(row),'gene_symbol':''.join(symbol.split()),'author_log2FC':float(lfc),
     'author_nominal_p_rounded':float(pvalue),'direction':'up'if float(lfc)>0 else'down','source_pdf_page':page_i,'source_line':line})
 s3=pd.DataFrame(rows);assert len(s3)==959 and s3.gene_symbol.nunique()==959
 assert s3.direction.value_counts().to_dict()=={'up':565,'down':394}
 assert s3.groupby('direction').source_row.apply(lambda x:set(x)==set(range(1,len(x)+1))).all()
 s3.to_csv(OUT/'author_ozone_S3_numeric.csv',index=False)
 lists={};list_rows=[]
 expected={'up':[1490,760,341],'down':[1011,383,207]}
 for code,direction in [('S2a','up'),('S2b','down')]:
  p=next(target.glob(f'Supplementary Table {code}*.pdf'))
  text='\n'.join(x.extract_text()for x in PdfReader(p).pages)
  pats=([r'Genotype E4\s*=1490',r'Age 17M\s*=760',r'Ozone\s*=341',r'Intersection for all']if direction=='up'else
   [r'Genotype E4\s*\(n=1011\)',r'Age 17M\s*\(n=383\)',r'Treatment; Ozone\s*\(n=207\)',r'Intersection for all'])
  matches=[re.search(pattern,text)for pattern in pats];assert all(matches)
  for i,factor in enumerate(['genotype','age','exposure']):
   block=text[matches[i].end():matches[i+1].start()]
   symbols=[''.join(x.split()).rstrip('.')for x in block.split(',')if x.strip()]
   assert len(symbols)==len(set(symbols))==expected[direction][i]
   lists[factor+'_'+direction]=set(symbols)
   list_rows.extend({'factor':factor,'direction':direction,'gene_symbol':s,'source_table':code}for s in symbols)
 pd.DataFrame(list_rows).to_csv(OUT/'author_ozone_S2_factor_gene_lists.csv.gz',index=False)
 author=pd.read_csv(ROOT/'data/E-MTAB-14917/processed/E4_vs_E3_005_total.csv').rename(columns={'Unnamed: 0':'gene_symbol'})
 assert author.gene_symbol.is_unique
 for direction,sign in [('up',1),('down',-1)]:
  selected=author.pvalue.lt(.05)&author.log2FoldChange.abs().gt(.3)&((author.log2FoldChange>0)if sign>0 else(author.log2FoldChange<0))
  assert set(author.loc[selected,'gene_symbol'])==lists['genotype_'+direction]
  assert set(author.loc[author.diffexpressed.eq(direction.upper()),'gene_symbol'])==lists['genotype_'+direction]
 author.to_csv(OUT/'author_ozone_genotype_table.csv.gz',index=False)
 return author,s3,lists

def mirna_statistics():
 x=pd.read_csv(ROOT/'results/mirna/analysis_signal.tsv',sep='\t',index_col=0)
 meta=pd.read_csv(ROOT/'results/mirna/sample_metadata.tsv',sep='\t',index_col=0)
 annotation=pd.read_csv(ROOT/'results/mirna/probe_annotation.tsv',sep='\t')
 groups=['Control','PM1','PM5'];indices=[meta.index[meta.condition.eq(c)]for c in groups]
 assert len(x)==1245 and all(len(i)==3 for i in indices)
 means=np.stack([x.loc[:,i].mean(axis=1)for i in indices],axis=1)
 residual=sum(((x.loc[:,i].to_numpy()-means[:,j,None])**2).sum(axis=1)for j,i in enumerate(indices))
 mse=residual/6;valid=np.isfinite(mse)&(mse>np.finfo(float).eps)
 F=np.full(len(x),np.nan);F[valid]=((3*(means-means.mean(axis=1)[:,None])**2).sum(axis=1)/2/mse)[valid]
 result=pd.DataFrame({'probe_id':x.index,'miRNA':x.index,'mean_Control':means[:,0],'mean_PM1':means[:,1],'mean_PM5':means[:,2],
  'within_group_MSE':mse,'residual_df':6,'valid_residual_variance':valid,'ANOVA_F':F,'ANOVA_pvalue':f_dist.sf(F,2,6)})
 independent=f_oneway(*(x.loc[:,i].to_numpy()for i in indices),axis=1)
 assert np.allclose(result.ANOVA_F,independent.statistic,equal_nan=True,atol=1e-10)
 assert np.allclose(result.ANOVA_pvalue,independent.pvalue,equal_nan=True,atol=1e-10)
 result['ANOVA_BH_FDR']=bh(result.ANOVA_pvalue)
 for j,name in [(1,'PM1'),(2,'PM5')]:
  effect=means[:,j]-means[:,0];se=np.sqrt(mse*2/3);tt=np.full(len(x),np.nan);tt[valid]=effect[valid]/se[valid]
  p=2*t_dist.sf(abs(tt),6)
  result[name+'_difference']=effect;result[name+'_LSD_SE']=se;result[name+'_LSD_t']=tt;result[name+'_LSD_pvalue']=p
  result[name+'_LSD_BH_FDR']=bh(p);result[name+'_LSD_CI95_lower']=effect-t_dist.ppf(.975,6)*se
  result[name+'_LSD_CI95_upper']=effect+t_dist.ppf(.975,6)*se
  result[name+'_protected_LSD_nominal']=np.isfinite(p)&(p<.05)&result.ANOVA_pvalue.lt(.05)
 result.to_csv(OUT/'mirna_author_statistics_candidates.csv.gz',index=False)
 rows=[]
 for name in ['PM1','PM5']:
  limma=pd.read_csv(ROOT/f'results/mirna/DE_{name}_vs_Control.tsv',sep='\t').set_index('probe_id')
  for gene in ['mmu-miR-574-5p','mmu-miR-5119']:
   r=result.set_index('probe_id').loc[gene];l=limma.loc[gene]
   rows.append({'miRNA':gene,'contrast':name+'_vs_Control','mean_difference':r[name+'_difference'],
    'ANOVA_pvalue':r.ANOVA_pvalue,'LSD_pvalue':r[name+'_LSD_pvalue'],'LSD_BH_FDR':r[name+'_LSD_BH_FDR'],
    'protected_LSD_nominal':r[name+'_protected_LSD_nominal'],'LSD_CI95_lower':r[name+'_LSD_CI95_lower'],
    'LSD_CI95_upper':r[name+'_LSD_CI95_upper'],'limma_pvalue':l.pvalue,'limma_BH_FDR':l.padj,
    'limma_CI95_lower':l.CI95_lower,'limma_CI95_upper':l.CI95_upper})
 targets=pd.DataFrame(rows);targets.to_csv(OUT/'mirna_author_targets.csv',index=False)
 counts=[]
 for name in ['PM1','PM5']:
  limma=pd.read_csv(ROOT/f'results/mirna/DE_{name}_vs_Control.tsv',sep='\t')
  counts.append({'contrast':name,'LSD_nominal_p005':int(result[name+'_LSD_pvalue'].lt(.05).sum()),
   'protected_LSD_nominal':int(result[name+'_protected_LSD_nominal'].sum()),
   'LSD_BH_FDR005':int(result[name+'_LSD_BH_FDR'].lt(.05).sum()),'limma_BH_FDR005':int(limma.padj.lt(.05).sum())})
 pd.DataFrame(counts).to_csv(OUT/'mirna_threshold_counts.csv',index=False)
 # New paper-style descriptive heatmap, not a claim to recover the unreported Figure6a row selection.
 selected_samples=list(indices[0])+list(indices[2]);six=x.loc[:,selected_samples]
 order=list(six.var(axis=1,ddof=1).sort_values(ascending=False).head(20).index)
 if 'mmu-miR-574-5p'not in order:order.append('mmu-miR-574-5p')
 z=six.loc[order];z=z.sub(z.mean(axis=1),axis=0).div(z.std(axis=1,ddof=1),axis=0)
 z=z.loc[z.index[leaves_list(linkage(z.to_numpy(),method='average',metric='euclidean'))]]
 z.to_csv(OUT/'mirna_paper_style_heatmap_zscores.tsv',sep='\t')
 return x,meta,result,targets,pd.DataFrame(counts),z

def representative(frame):
 return frame.dropna(subset=['gene_symbol']).sort_values(['baseMean','gene_id'],ascending=[False,True]).drop_duplicates('gene_symbol')

def comparison_stats(merged,author_lfc='author_log2FC',ours_lfc='log2FoldChange',author_selected=None):
 valid=merged[[author_lfc,ours_lfc]].notna().all(axis=1)&np.isfinite(merged[[author_lfc,ours_lfc]]).all(axis=1)
 d=merged.loc[valid];same=np.sign(d[author_lfc])==np.sign(d[ours_lfc])
 result={'common_exact_symbol_rows':len(merged),'common_finite_effects':len(d),
  'effect_Spearman_rho':spearmanr(d[author_lfc],d[ours_lfc]).statistic,
  'effect_Pearson_r':pearsonr(d[author_lfc],d[ours_lfc]).statistic,
  'same_effect_direction':int(same.sum()),'direction_fraction':same.mean(),
  'median_absolute_log2FC_difference':(d[author_lfc]-d[ours_lfc]).abs().median()}
 return result

def ozone_comparisons(author,s3,lists):
 fits=json.loads((OUT/'ozone_candidate_models.json').read_text())
 tables={k:representative(pd.read_csv(OUT/f'ozone_candidate_{k}.csv.gz'))for k in fits['models']}
 summaries={};joins={}
 for key in ['genotype_pooled','genotype_additive']:
  merge=author.merge(tables[key],on='gene_symbol',suffixes=('_author','_ours'),validate='one_to_one')
  merge['author_selected']=merge.pvalue_author.lt(.05)&merge.log2FoldChange_author.abs().gt(.3)
  merge['ours_selected']=merge.author_criteria_selected
  merge.to_csv(OUT/f'ozone_concordance_{key}.csv.gz',index=False)
  v=merge[['log2FoldChange_author','log2FoldChange_ours','pvalue_author','pvalue_ours']].notna().all(axis=1)
  common=merge.loc[v];a=set(common.loc[common.author_selected,'gene_symbol']);b=set(common.loc[common.ours_selected,'gene_symbol'])
  summary=comparison_stats(merge,'log2FoldChange_author','log2FoldChange_ours')
  summary.update(author_full_selected=2501,author_full_genes=len(author),ours_full_selected=fits['models'][key]['nominal_selected'],
   finite_both_tests=len(common),author_selected_in_joint_testable_universe=len(a),ours_selected_in_joint_testable_universe=len(b),
   overlapping_selected_in_joint_testable_universe=len(a&b),selected_Jaccard=len(a&b)/len(a|b),
   author_selection_recall_in_joint_testable_universe=len(a&b)/len(a),ours_selection_precision_in_joint_testable_universe=len(a&b)/len(b))
  summaries[key]=summary;joins[key]=merge
 merge=s3.merge(tables['old_ozone_genotype'],on='gene_symbol',validate='one_to_one')
 merge.to_csv(OUT/'ozone_concordance_old_ozone_genotype.csv.gz',index=False)
 summary=comparison_stats(merge);summary.update(author_full_selected=959,ours_full_selected=fits['models']['old_ozone_genotype']['nominal_selected'],author_up=565,author_down=394,
  exact_symbol_coverage=len(merge)/959,author_genes_retaining_our_nominal_selection=int(merge.author_criteria_selected.sum()),
  author_genes_retaining_our_BH_FDR05=int(merge.BH_FDR05.sum()))
 summaries['old_ozone_genotype']=summary;joins['old_ozone_genotype']=merge
 # The exact-symbol universe is fixed before intersection calculation and shared across both sources.
 universe=set(author.gene_symbol)
 for key in ['genotype_pooled','age_pooled','exposure_pooled']:
  d=tables[key];universe &= set(d.loc[d.pvalue.notna()&np.isfinite(d.log2FoldChange),'gene_symbol'])
 factors=['genotype','age','exposure'];members=pd.DataFrame({'gene_symbol':sorted(universe)})
 for direction in ['up','down']:
  for factor,key in zip(factors,['genotype_pooled','age_pooled','exposure_pooled']):
   d=tables[key];calls=set(d.loc[d.author_criteria_selected&((d.log2FoldChange>0)if direction=='up'else(d.log2FoldChange<0)),'gene_symbol'])
   members[f'author_{factor}_{direction}']=members.gene_symbol.isin(lists[factor+'_'+direction])
   members[f'ours_{factor}_{direction}']=members.gene_symbol.isin(calls)
 members.to_csv(OUT/'ozone_factor_membership_common_universe.csv.gz',index=False)
 patterns=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)];counts=[]
 for direction in ['up','down']:
  for source in ['author','ours']:
   flags=members[[f'{source}_{factor}_{direction}'for factor in factors]].astype(int).to_numpy()
   for pattern in patterns:counts.append({'direction':direction,'source':source,'pattern':''.join(map(str,pattern)),
    'genes':int((flags==pattern).all(axis=1).sum()),'universe':len(universe)})
 counts=pd.DataFrame(counts);counts.to_csv(OUT/'ozone_factor_upset_counts.csv',index=False)
 # Planned biological targets come from the original S1/S3/Figure4 discussion, not from current smallest p-values.
 target_genes=['Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Cacnb2','Tfrc','Bcas1','Sfrp5','Folr1','Bok','Pon3','Barhl1','Casp1','Apaf1']
 targets=merge.loc[merge.gene_symbol.isin(target_genes)].copy();targets.to_csv(OUT/'ozone_author_target_effects.csv',index=False)
 for key,d in tables.items():
  # Independently verify BH family, Wald tails and final candidate selection.
  original=pd.read_csv(OUT/f'ozone_candidate_{key}.csv.gz');ok=np.isfinite(original.pvalue)
  assert np.allclose(original.padj,bh(original.pvalue),equal_nan=True,atol=1e-10)
  assert np.allclose(original.loc[ok,'pvalue'],2*norm.sf(abs(original.loc[ok,'stat'])),atol=1e-10)
  selected=ok&original.pvalue.lt(.05)&original.log2FoldChange.abs().gt(fits['models'][key]['absolute_log2FC_cutoff'])
  assert np.array_equal(selected,original.author_criteria_selected)
  assert int(selected.sum())==fits['models'][key]['nominal_selected']
 return fits,tables,summaries,joins,members,counts,targets

def mirna_plot(x,meta,targets,counts,z,spanish=False):
 tr=lambda a,b:b if spanish else a
 fig=plt.figure(figsize=(12,9.5),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1,1.18],height_ratios=[.85,1.15])
 ax=fig.add_subplot(gs[0,0]);gene='mmu-miR-574-5p'
 for i,group in enumerate(['Control','PM1','PM5']):
  vals=x.loc[gene,meta.index[meta.condition.eq(group)]].to_numpy();shape=['o','s','^'][i]
  ax.scatter(i+np.array([-.08,0,.08]),vals,color=COLORS[group],marker=shape,s=35,zorder=3)
  ax.errorbar(i,vals.mean(),yerr=vals.std(ddof=1)/np.sqrt(3),color='black',capsize=5,fmt='_',markersize=17,zorder=4)
 ax.set(xticks=[0,1,2],xticklabels=[tr('Control','Control'),'1 mg/kg','5 mg/kg'],ylabel=tr('Published signal units','Unidades de señal publicadas'),title='A  miR-574-5p · n = 3 / '+tr('group','grupo'))
 ax=fig.add_subplot(gs[0,1]);r=targets[targets.miRNA.eq(gene)].copy();labels=[]
 for j,(_,row)in enumerate(r.iterrows()):
  for method,offset,color,marker in [('limma',-.10,'#087f8c','o'),('LSD',.10,'#aa5533','D')]:
   lo,hi=row[method+'_CI95_lower'],row[method+'_CI95_upper'];y=j+offset
   ax.errorbar(row.mean_difference,y,xerr=[[row.mean_difference-lo],[hi-row.mean_difference]],fmt=marker,color=color,capsize=3,markersize=5,label=method if j==0 else None)
  labels.append('1 mg/kg'if j==0 else'5 mg/kg')
 ax.axvline(0,color='#777777',ls=':',lw=1);ax.set(yticks=[0,1],yticklabels=labels,xlabel=tr('Difference versus control (signal units)','Diferencia frente a control (unidades de señal)'),title=tr('B  Pointwise 95% intervals; not FDR intervals','B  Intervalos puntuales del 95%; no de FDR'))
 ax.set_ylim(1.85,-.55);ax.legend(loc='upper right',fontsize=8)
 high=r.iloc[1]
 ax.text(.02,.03,tr(f'High dose: LSD p={high.LSD_pvalue:.4f}; limma q={high.limma_BH_FDR:.3f}\nOmnibus ANOVA p={high.ANOVA_pvalue:.4f}: protected LSD does not pass',f'Dosis alta: p LSD={high.LSD_pvalue:.4f}; q limma={high.limma_BH_FDR:.3f}\nANOVA global p={high.ANOVA_pvalue:.4f}: no supera la LSD protegida'),transform=ax.transAxes,fontsize=7.7)
 ax=fig.add_subplot(gs[1,0]);xx=np.arange(2);width=.24
 for offset,(col,label,color,hatch) in enumerate([('LSD_nominal_p005',tr('LSD nominal p < 0.05','LSD: p nominal < 0.05'),'#cc8430',None),('LSD_BH_FDR005',tr('LSD BH FDR < 0.05','LSD: FDR BH < 0.05'),'#725a94','//'),('limma_BH_FDR005',tr('limma BH FDR < 0.05','limma: FDR BH < 0.05'),'#087f8c','..')]):
  vals=counts[col].to_numpy();positions=xx+(offset-1)*width
  ax.bar(positions,vals,width,color=color,label=label,hatch=hatch,edgecolor='white')
  for a,v in zip(positions,vals):ax.annotate(str(v),(a,v),xytext=(0,4),textcoords='offset points',ha='center',fontsize=8,color=color)
 ax.set(xticks=xx,xticklabels=['1 mg/kg','5 mg/kg'],ylabel=tr('Selected miRNA features','Características de miRNA seleccionadas'),ylim=(0,max(counts.LSD_nominal_p005.max()*1.35,3)),title=tr('C  Model and multiplicity alter discovery counts','C  El modelo y la multiplicidad cambian los hallazgos'))
 ax.legend(fontsize=7.4,loc='upper right')
 ax=fig.add_subplot(gs[1,1]);im=ax.imshow(z,vmin=-2.5,vmax=2.5,cmap='RdBu_r',aspect='auto',interpolation='nearest')
 ax.set(yticks=range(len(z)),yticklabels=z.index,xticks=range(6),xticklabels=['C1','C2','C3','H1','H2','H3'],title=tr('D  New paper-style heatmap; original row selection unknown','D  Nuevo mapa de calor; selección original no descrita'))
 ax.tick_params(axis='y',labelsize=7);ax.tick_params(axis='x',labelsize=8)
 for i in range(6):ax.add_patch(plt.Rectangle((i-.5,-.5),1,.35,color=COLORS['Control'if i<3 else'PM5'],clip_on=False))
 fig.colorbar(im,ax=ax,shrink=.65,label=tr('Within-feature signal z score','Puntuación z de señal por característica'))
 save_fig(fig,'mirna_author_reproduction'+('_es'if spanish else''))

def ozone_plot(joins,targets,fits,spanish=False):
 tr=lambda a,b:b if spanish else a
 fig,axes=plt.subplots(2,2,figsize=(12,10),layout='constrained')
 specs=[('genotype_pooled','log2FoldChange_author','log2FoldChange_ours',tr('A  Global genotype · candidate ~genotype','A  Genotipo global · candidato ~genotype')),
 ('old_ozone_genotype','author_log2FC','log2FoldChange',tr('B  17M O3 · E4 versus E3 · author S3','B  17M O3 · E4 frente a E3 · S3 original'))]
 for ax,(key,a,b,title)in zip(axes[0],specs):
  d=joins[key].dropna(subset=[a,b]);v=d[[a,b]].to_numpy();same=np.sign(v[:,0])==np.sign(v[:,1])
  ax.scatter(v[same,0],v[same,1],s=5,alpha=.28,color='#287a86',rasterized=True,label=tr('Same direction','Misma dirección'))
  ax.scatter(v[~same,0],v[~same,1],s=7,alpha=.45,color='#c08022',marker='x',rasterized=True,label=tr('Opposite direction','Dirección opuesta'))
  lim=max(np.abs(v).max()*1.05,1);ax.plot([-lim,lim],[-lim,lim],ls=':',color='#555555',lw=1)
  ax.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel=tr('Author unshrunken log₂FC','log₂FC original de los autores'),ylabel=tr('Candidate unshrunken log₂FC','log₂FC del candidato'),title=title)
  rho=spearmanr(v[:,0],v[:,1]).statistic;ax.text(.04,.94,f'n = {len(d):,}; ρ = {rho:.3f}',transform=ax.transAxes,va='top',fontsize=9)
  ax.legend(fontsize=7.8,loc='lower right')
 ax=axes[1,0]
 genes=['Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Cacnb2','Tfrc','Bcas1'];t=targets.set_index('gene_symbol').loc[genes]
 for i,(gene,row)in enumerate(t.iterrows()):
  ax.scatter(row.author_log2FC,i-.12,marker='D',s=27,color='#826196',label=tr('Published S3 point estimate','Estimación publicada S3')if i==0 else None)
  ax.errorbar(row.log2FoldChange,i+.12,xerr=[[row.log2FoldChange-row.CI95_lower],[row.CI95_upper-row.log2FoldChange]],fmt='o',markersize=4.5,color='#287a86',capsize=2,label=tr('Candidate MLE ± 1.96 SE','MLE candidata ± 1.96 EE')if i==0 else None)
 ax.set(yticks=range(len(genes)),yticklabels=genes,xlabel=tr('E4/E3 log₂FC in 17M O3 mice','log₂FC E4/E3 en ratones 17M O3'),title=tr('C  Prespecified paper targets and uncertainty','C  Dianas del artículo e incertidumbre'))
 ax.invert_yaxis();ax.axvline(0,color='#777777',ls=':',lw=1);ax.legend(fontsize=7.7,loc='lower left')
 ax=axes[1,1];labels=[tr('Global genotype','Genotipo global'),tr('17M O3 genotype','Genotipo 17M O3')];xx=np.arange(2)
 source=[2501,959];ours=[fits['models'][k]['nominal_selected'] for k in ['genotype_pooled','old_ozone_genotype']]
 ax.bar(xx-.18,source,.36,color='#826196',label=tr('Author nominal criteria','Criterios nominales de los autores'))
 ax.bar(xx+.18,ours,.36,color='#287a86',label=tr('Candidate nominal criteria','Criterios nominales del candidato'),hatch='//',edgecolor='white')
 for x0,v in zip(xx-.18,source):ax.text(x0,v+45,f'{v:,}',ha='center',fontsize=8)
 for x0,v in zip(xx+.18,ours):ax.text(x0,v+45,f'{v:,}',ha='center',fontsize=8)
 ax.set(xticks=xx,xticklabels=labels,ylabel=tr('Selected gene features','Características génicas seleccionadas'),ylim=(0,2900),title=tr('D  Full lists; different quantified/filter universes','D  Listas completas; universos de cuantificación/filtro distintos'))
 ax.legend(fontsize=7.5,loc='upper right')
 save_fig(fig,'ozone_author_concordance'+('_es'if spanish else''))

def upset_plot(counts,spanish=False):
 tr=lambda a,b:b if spanish else a
 fig=plt.figure(figsize=(12,8),layout='constrained');outer=fig.add_gridspec(2,2)
 patterns=['100','010','001','110','101','011','111']
 for i,(direction,source)in enumerate([('up','author'),('up','ours'),('down','author'),('down','ours')]):
  sub=outer[i//2,i%2].subgridspec(2,1,height_ratios=[3,1],hspace=.03);ax=fig.add_subplot(sub[0]);mat=fig.add_subplot(sub[1],sharex=ax)
  d=counts[counts.direction.eq(direction)&counts.source.eq(source)].set_index('pattern').loc[patterns];c='#287a86'if direction=='up'else'#c08022'
  ax.bar(range(7),d.genes,color=c,hatch='//'if source=='ours'else None,edgecolor='white',linewidth=.5)
  for j,value in enumerate(d.genes):ax.text(j,value+max(d.genes.max()*.03,1),str(value),ha='center',fontsize=8)
  ax.set(ylim=(0,counts.loc[counts.direction.eq(direction),'genes'].max()*1.24+1),ylabel=tr('Genes in intersection','Genes en la intersección'),title=f'{"ABCD"[i]}  '+(tr('Published','Publicado')if source=='author'else tr('Candidate','Candidato'))+' · '+(tr('higher','aumentados')if direction=='up'else tr('lower','disminuidos')))
  ax.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
  for j,p in enumerate(patterns):
   on=[k for k,a in enumerate(p)if a=='1'];mat.scatter([j]*3,range(3),s=19,color='#dddddd');mat.plot([j]*len(on),on,color='#333333',lw=1);mat.scatter([j]*len(on),on,s=24,color='#333333')
  mat.set(yticks=[0,1,2],yticklabels=[tr('Genotype','Genotipo'),tr('Age','Edad'),tr('Exposure','Exposición')],ylim=(2.5,-.5),xticks=[])
  mat.spines[['left','bottom']].set_visible(False);mat.tick_params(axis='y',length=0,labelsize=8)
 save_fig(fig,'ozone_factor_upset_reproduction'+('_es'if spanish else''))

def outputs(author,s3,lists,mirna,ozone):
 x,meta,stats,targets,counts,z=mirna;fits,tables,summaries,joins,members,upset,otargets=ozone
 high=targets[(targets.miRNA=='mmu-miR-574-5p')&(targets.contrast=='PM5_vs_Control')].iloc[0]
 selected_s3=summaries['old_ozone_genotype'];global_s=summaries['genotype_pooled']
 audits=[
 {'study':'GSE93967','source_claim':'Figure 6: miR-574-5p is lower after high-dose PM2.5','original_assay':'Deposited microarray subset: 3 control / 3 high-dose animals; source qPCR is separate','attempt':'Mean-signal contrast, original-like general ANOVA/LSD candidate and primary limma','observed':f'Mean difference {high.mean_difference:.6f}; LSD p={high.LSD_pvalue:.6f}; limma p={high.limma_pvalue:.6f}, BH q={high.limma_BH_FDR:.6f}','status':'Direction and nominal pairwise evidence reproduced; not FDR-wide confirmation','limitation':'Omnibus ANOVA p=.0788; protected-LSD gate fails. GeneSpring configuration and array feature-selection rule are not reported.'},
 {'study':'GSE93967','source_claim':'miR-574-5p/BACE1/NF-kB mechanism and rescue','original_assay':'qPCR, protein, reporter, ChIP, behavioral intervention','attempt':'None: raw independent assay measurements not present in deposited microarray','observed':'Not tested with expression arrays','status':'Not independently reproducible from this deposit','limitation':'Array association neither establishes nor refutes causal rescue experiments.'},
 {'study':'E-MTAB-14917','source_claim':'S2 global genotype: 1490 up + 1011 down genes','original_assay':'HISAT2/featureCounts RNA-seq; nominal p < .05 and |log2FC| > .3','attempt':'Exact internal reconstruction of deposited author CSV labels; candidate Salmon/tximport genotype-only and additive models','observed':f'AuthorCSV labels exactly match S2 gene sets and threshold; candidate {global_s["ours_full_selected"]} genes; common finite effects {global_s["common_finite_effects"]}, Spearman rho = {global_s["effect_Spearman_rho"]:.4f}','status':'Author-table selection exactly reproduced; biological-count reanalysis shows partial concordance','limitation':'Original formula/filter/annotation release unavailable; current quantifier and DESeq2 version differ. Comparison uses same animals, not independent biological replication.'},
 {'study':'E-MTAB-14917','source_claim':'S3 old 17M O3 E4 vs E3: 565 up + 394 down','original_assay':'HISAT2/featureCounts RNA-seq; nominal p < .05 and |log2FC| > .1','attempt':'Six-animal genotype comparison, refit dispersions, fixed global count filter and inherited length offsets','observed':f'{selected_s3["common_finite_effects"]}/959 exact-symbol finite effects; Spearman rho = {selected_s3["effect_Spearman_rho"]:.4f}; direction agreement = {selected_s3["direction_fraction"]:.2%}','status':'Strong effect-direction concordance within the matched published gene list; selected-list size not identical','limitation':'Published S3 contains only selected genes; correlation is conditional on that selection and same-source animals. Rounded p values cannot be compared as exact raw tests.'},
 {'study':'E-MTAB-14917','source_claim':'Genotype-related hemoglobin differences support a BBB mechanism','original_assay':'RNA-seq, selected qPCR and biological interpretation','attempt':'S3 effect/pointwiseCI comparison of Hba-a1,Hba-a2,Hbb-bs,Hbb-bt','observed':'All four have negative E4/E3 effects and BH FDR < .05 in candidate; remaining paper targets have mixed nominal/FDR support','status':'Hemoglobin expression direction/significance reproduced; BBB function is not inferred','limitation':'Bulk cell composition, residual vascular material and expression regulation cannot be separated by these data; functional interpretation needs independent measurements.'},
 ]
 pd.DataFrame(audits).to_csv(OUT/'paper_replication_audit.csv',index=False)
 sources=[ROOT/'papers/pdf/GSE93967_Ku_2017_Particle_Fibre_Toxicology.pdf',ROOT/'papers/pdf/E-MTAB-14917_Nakamya_2025_IJMS.pdf',
 ROOT/'data/E-MTAB-14917/processed/E4_vs_E3_005_total.csv',ROOT/'data/E-MTAB-14917/supplements/PMC11942628_supplementary.zip',ROOT/'results/mirna/analysis_signal.tsv',ROOT/'results/ozone/dds_full_wald.rds']
 source_manifest={'original_papers':[{'doi':'10.1186/s12989-017-0215-3','methods_local_pages':[14,16],'figure6_local_page':8,'software':'Agilent Feature Extraction 10.7 / GeneSpring 11.0','general_statistics':'One-way ANOVA then LSD; omnibus protection and microarray-specific selection/multiplicity not specified'},
 {'doi':'10.3390/ijms26062407','methods_local_page':15,'statistical_criteria':{'global_genotype_and_age':'nominal p < .05 and |log2FC| > .3','global_exposure_and_individual_comparisons':'nominal p < .05 and |log2FC| > .1'},'software_reported':'HISAT2/featureCounts/GRCm38; DESeq2 1.440 as printed; R 4.4.1','software_used':'Salmon 2.5.1 / tximport 1.38.2 / GENCODE M25; DESeq2 1.50.2'}],
 'private_source_files':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)}for p in sorted((META/'private_inputs').rglob('*'))if p.is_file()],
 'private_source_inputs':'Full-text extracts, publisher artwork and unpacked supplementary PDFs/text are retained locally under metadata/replication_mirna_ozone/private_inputs and excluded from public publishing. Derived numeric comparison tables are public.',
 'source_files':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)}for p in sources],
 'source_urls':['https://doi.org/10.1186/s12989-017-0215-3','https://doi.org/10.3390/ijms26062407','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93967','https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14917'],
 'supplement_parsing':'Source nested s001.zip is already preserved; S3 PDF numeric rows require removing PDF-inserted whitespace within gene names; all 959 row indices (565 up / 394 down) verified. S2componentlists exactly1490/1011 genotype760/383 age341/207 exposure; globalgenotype sourceCSV selected sets match S2 exactly.'}
 save_json(META/'source_manifest.json',source_manifest)
 methods={'purpose':'Attempt method/target/result reproduction, not independent biological replication','miRNA':{'input':'Deposited quantile-normalized continuous signal; same 1245 features as the primary analysis; 9 arrays, 3 per group','candidate':'One-way three-group ANOVA with residual df = 6; pairwise LSD uses pooled within-group MSE; the protected gate is reported separately','multiplicity':'Raw pairwise p and BH within 1245 features per contrast exported separately; separate ANOVA family','scale':'Signal-unit mean differences; no invented log2 scale or fold-change ratios','zero_within_group_variance':'Nonpositive or numerically zero residual MSE produces withheld tests','heatmap':'20 largest-variance features among six control/high-dose arrays plus prespecified miR-574-5p if absent; row z scores (ddof = 1), average linkage and Euclidean distance, fixed columns C1–C3/H1–H3; original Figure 6a selection unreported'},
 'ozone':{'candidate_designs':'All nine candidates are prespecified in replicate_mirna_ozone.R; marginal genotype/age/exposure, additive genotype sensitivity, six-animal old O3 genotype and four 6-animal within-stratum exposures','normalization':'Existing effective-length/size factors retained from full 24-library tximport dataset; refit dispersions per model/subset; inheritance differs from unknown original preprocessing','filter':'Fixed primary 18537-gene filter: ≥10 rounded fragments in ≥3 of 24 libraries; original author filter unknown','outliers':'DESeq2 default minReplicatesForReplace = 7 for author-style candidates; flagged gene counts recorded; no sample removals','inference':'Convergence checked, Cook filter enabled, independent filtering disabled for complete BH families; source-specific nominal .3/.1 effect criteria applied to unshrunken MLEs','gene_matching':'Exact symbols only; duplicated candidate symbols represented by highest baseMean (then gene ID), independently of p/effect; nonmatches retained in coverage accounting','figure1_style':'UpSet counts use one exact-symbol shared universe testable in all 3 candidate factor models and present in the author full genotype gene index; not pooled causal or interaction inference','limitations':['Quantifier/annotation and DESeq2 versions differ','Original exact design/filter/GeneSpring project not deposited','Author S3 is a selected list; agreement conditional on selection','Same-source reanalysis is not independent biological replication','No inferred BBB function, cognition or cell-specific response from bulk expression']},
 'reproduce':['Rscript scripts/replicate_mirna_ozone.R','.venv-analysis/bin/python scripts/replicate_mirna_ozone.py --reuse-models'],
 'figures':'All individual miRNA animals shown; pointwise 95% intervals distinguished from FDR; full scatter ranges and exact-match denominators; all 7 nonempty UpSet factor patterns, including zeros'}
 save_json(OUT/'methods.json',methods)
 summary={'complete':True,'mirna':{'features':len(stats),'valid_residual_variance':int(stats.valid_residual_variance.sum()),'author_target_high_dose':high.to_dict(),'threshold_counts':counts.to_dict('records')},
 'ozone':{'candidate_models':fits['models'],'author_concordance':summaries,'factor_upset_common_universe':len(members),'source_S2_full_counts':{key:len(v)for key,v in lists.items()},'source_S3_rows':len(s3)},
 'replication_levels':['Exact reconstruction of source table selection labels','Targeted direction/nominal evidence assessment','Conditional published-effect concordance','Methodological candidates with explicit quantifier, filter, design and version limits'],
 'all_primary_results_preserved':True,'no_new_spatial_exposure_measurements':True}
 save_json(OUT/'summary.json',summary)
 captions={
 'mirna_author_reproduction':{
 'en':'Four-panel microarray reproduction audit. A: all nine source animals for the prespecified miR-574-5p target; black marks show mean ± SE. B: signal differences with limma and unmoderated LSD pointwise 95% intervals. C: nominal and BH-adjusted discovery counts across 1,245 fixed features, with each contrast adjusted separately. D: a new descriptive control/high-dose heatmap of the 20 highest-variance features plus the author target, row z scored and hierarchically ordered (average linkage, Euclidean distance); C1–C3 = control, H1–H3 = high dose. The original Figure 6a row selection is unreported. qPCR and functional experiments are separate assays.',
 'es':'Auditoría del microarreglo en cuatro paneles. A: los nueve animales de origen para miR-574-5p; las marcas negras indican media ± EE. B: diferencias de señal con intervalos puntuales del 95% de limma y LSD no moderada. C: hallazgos nominales y con ajuste BH entre 1.245 características fijas; ajuste separado de cada contraste. D: nuevo mapa descriptivo de control/dosis alta con las 20 características de mayor varianza más la diana original, puntuaciones z por fila y orden jerárquico (enlace promedio, distancia euclídea); C1–C3 = control, H1–H3 = dosis alta. No se describe la selección de filas de la Figura 6a original. qPCR y experimentos funcionales son ensayos separados.'},
 'ozone_author_concordance':{
 'en':'Reproduction candidates versus source author results. A: all common finite exact-symbol genotype effects. B: 17-month O3 E4/E3 effects conditional on author S3 selection. C: prespecified targets, published point estimates and candidate MLE ± 1.96 SE pointwise intervals. D: full nominal selection counts from different quantified and filtered universes; counts alone do not measure replication. Global genotype uses nominal p < 0.05 and |log2FC| > 0.3; the 17-month O3 subgroup uses p < 0.05 and |log2FC| > 0.1. Current Salmon quantification differs from the source HISAT2/featureCounts workflow. These comparisons use the same animals, rather than an independent cohort.',
 'es':'Candidatos de reproducción frente a resultados originales. A: efectos de genotipo finitos y con símbolos exactos en ambas fuentes. B: efectos E4/E3 a los 17 meses bajo O3, condicionados a la selección publicada S3. C: dianas preespecificadas, estimaciones publicadas e intervalos puntuales candidatos MLE ± 1,96 EE. D: conteos nominales completos de universos de cuantificación y filtro distintos; los conteos no miden por sí solos la reproducción. Genotipo global: p nominal < 0,05 y |log2FC| > 0,3; subgrupo de 17 meses bajo O3: p < 0,05 y |log2FC| > 0,1. Salmon difiere del flujo HISAT2/featureCounts original. Se comparan los mismos animales, no una cohorte independiente.'},
 'ozone_factor_upset_reproduction':{
 'en':f'Factor-intersection comparison inspired by original Figure 1, restricted to the same {len(members):,} exact-symbol genes testable in all three candidate models. Original S2 lists and candidate marginal models are shown separately for higher and lower expression, with matched vertical scales within each direction. All seven nonempty genotype/age/exposure membership patterns are shown, including zero counts. Global genotype and age use nominal p < 0.05 and |log2FC| > 0.3; exposure uses p < 0.05 and |log2FC| > 0.1. These correlated list memberships are not interaction tests or evidence of relative pollutant toxicity.',
 'es':f'Comparación de intersecciones inspirada en la Figura 1 original, restringida a los mismos {len(members):,} genes con símbolos exactos y pruebas válidas en los tres modelos candidatos. Listas S2 originales y modelos marginales candidatos, con aumentos y disminuciones separados y escalas verticales comunes dentro de cada dirección. Se muestran los siete patrones de pertenencia a genotipo/edad/exposición, incluidos los conteos cero. Genotipo y edad globales: p nominal < 0,05 y |log2FC| > 0,3; exposición: p < 0,05 y |log2FC| > 0,1. Estas pertenencias correlacionadas no son pruebas de interacción ni evidencias de toxicidad relativa.'}}
 save_json(OUT/'figure_captions.json',captions)
 en=f"""# Author-result reproduction: miRNA and ozone

The deposited miR-574-5p signal is lower at both PM doses. At 5 mg/kg, the mean difference is {high.mean_difference:.4f} source signal units, with nominal LSD p = {high.LSD_pvalue:.4f} and primary limma p = {high.limma_pvalue:.4f}, but BH q = {high.limma_BH_FDR:.3f}. The direction and nominal pairwise evidence agree with the array claim; omnibus ANOVA p = {high.ANOVA_pvalue:.4f} fails a protected LSD gate. The paper does not specify that gate or its array-specific selection and multiplicity settings. These arrays do not reproduce the separate qPCR, BACE1, ChIP, rescue or behavioral experiments.

The ozone source genotype CSV selection is exactly reconstructed: 1,490 higher and 1,011 lower genes at nominal p < 0.05 and |log2FC| > 0.3; these sets equal supplementary S2. Applying candidate DESeq2 models to the same 24 libraries quantified with Salmon gives partial concordance. The genotype-only candidate has {global_s['common_finite_effects']:,} common finite effects and Spearman ρ = {global_s['effect_Spearman_rho']:.3f}; it selects {global_s['ours_full_selected']:,} genes versus 2,501 published genes. Within the jointly testable exact-symbol universe, {global_s['overlapping_selected_in_joint_testable_universe']:,} selected genes overlap (Jaccard = {global_s['selected_Jaccard']:.3f}).

For E4/E3 among 17-month O3 mice, {selected_s3['common_finite_effects']} of 959 published S3 genes have finite exact-symbol estimates, with ρ = {selected_s3['effect_Spearman_rho']:.3f} and {100*selected_s3['direction_fraction']:.1f}% direction agreement. Of the matched genes, {selected_s3['author_genes_retaining_our_nominal_selection']} retain the candidate nominal criterion and {selected_s3['author_genes_retaining_our_BH_FDR05']} retain BH FDR < 0.05. The candidate's full subgroup list contains {selected_s3['ours_full_selected']:,} genes. This selected-list comparison is conditional on the authors' selection and uses the same animals, rather than independent biological replication.

All four hemoglobin genes retain negative E4/E3 effects and BH significance. Several other highlighted targets retain their direction without equivalent nominal or FDR evidence. Bulk composition, vascular material and expression regulation cannot be distinguished here; BBB function is not inferred. Original HISAT2/featureCounts, the exact design and filter, annotation release and historical DESeq2 version are not fully reproduced. Candidate pooled models use DESeq2's default count-replacement setting, with affected-gene counts recorded and no animals removed. The new UpSet panels compare component-list membership on one shared universe; they are not interaction tests.
"""
 es=f"""# Reproducción de resultados originales: miRNA y ozono

La señal depositada de miR-574-5p disminuye con ambas dosis de PM. A 5 mg/kg, la diferencia media es {high.mean_difference:.4f} unidades de señal originales, con p nominal LSD = {high.LSD_pvalue:.4f} y p de limma = {high.limma_pvalue:.4f}, pero q BH = {high.limma_BH_FDR:.3f}. La dirección y la evidencia nominal por pares concuerdan con la afirmación del microarreglo; ANOVA global p = {high.ANOVA_pvalue:.4f} no supera la condición de LSD protegida. El artículo no especifica esa condición ni las opciones de selección y multiplicidad del microarreglo. Estos datos no reproducen los ensayos separados de qPCR, BACE1, ChIP, rescate o conducta.

Se reconstruye exactamente la selección del CSV original de genotipo de ozono: 1.490 genes aumentados y 1.011 disminuidos con p nominal < 0,05 y |log2FC| > 0,3; los conjuntos coinciden con S2. Los modelos DESeq2 candidatos sobre las mismas 24 bibliotecas cuantificadas con Salmon muestran concordancia parcial. El candidato de genotipo tiene {global_s['common_finite_effects']:,} efectos finitos compartidos y ρ de Spearman = {global_s['effect_Spearman_rho']:.3f}; selecciona {global_s['ours_full_selected']:,} genes frente a 2.501 publicados. En el universo de símbolos exactos con pruebas válidas en ambas fuentes, se comparten {global_s['overlapping_selected_in_joint_testable_universe']:,} genes seleccionados (Jaccard = {global_s['selected_Jaccard']:.3f}).

En E4/E3 a los 17 meses bajo O3, {selected_s3['common_finite_effects']} de 959 genes S3 tienen estimaciones finitas y símbolos exactos, con ρ = {selected_s3['effect_Spearman_rho']:.3f} y {100*selected_s3['direction_fraction']:.1f}% de concordancia direccional. De los genes emparejados, {selected_s3['author_genes_retaining_our_nominal_selection']} mantienen el criterio nominal candidato y {selected_s3['author_genes_retaining_our_BH_FDR05']} mantienen FDR BH < 0,05. La lista candidata completa del subgrupo contiene {selected_s3['ours_full_selected']:,} genes. La comparación está condicionada a la selección original y utiliza los mismos animales; no constituye replicación biológica independiente.

Los cuatro genes de hemoglobina conservan efectos E4/E3 negativos y significación BH. Otras dianas conservan la dirección sin evidencia nominal o de FDR equivalente. Aquí no se pueden separar composición celular, material vascular y regulación de expresión; no se infiere función de la barrera hematoencefálica. No se reproduce totalmente HISAT2/featureCounts, el diseño y filtro exactos, la versión de anotación ni DESeq2 histórico. Los modelos globales candidatos usan la sustitución de conteos predeterminada de DESeq2, registran los genes afectados y no eliminan animales. Los nuevos paneles UpSet comparan pertenencias a listas en un universo compartido; no son pruebas de interacción.
"""
 (OUT/'narrative_en.md').write_text(en);(OUT/'narrative_es.md').write_text(es)
 validation={'passed':True,'source_S3_959_unique_rows_565up394down':True,'source_S2_full_factor_counts_reconciled':True,'source_genotype_CSV_labels_and_S2_gene_sets_exact':True,'nine_candidate_Wald_BH_criteria_recomputed':True,'source_sample_mapping_verified':True,'miRNA_ANOVA_independent_scipy_f_oneway_match':True,'no_primary_outputs_modified':True,'all_nine_candidate_models_converged':all(v['nonconverged_nonzero_genes']==0 for v in fits['models'].values()),'scope':'Arithmetic/source-matching and convergence checks, not independent biological replication','figure_review':'Pending actual image inspection'}
 review_file=OUT/'figure_review.json'
 if review_file.exists():
  review=json.loads(review_file.read_text())
  matches=all((ROOT/f['path']).exists()and sha(ROOT/f['path'])==f['sha256']for f in review.get('figures',[]))
  if matches and len(review.get('figures',[]))==6:
   validation['figure_review']=review['status'];validation['figure_review_record']='results/replication_mirna_ozone/figure_review.json'
 save_json(OUT/'validation.json',validation)


def sample_audit():
 sdrf=pd.read_csv(ROOT/'metadata/E-MTAB-14917/E-MTAB-14917.sdrf.txt',sep='\t')
 samples=pd.read_csv(ROOT/'results/ozone/samples.csv')
 cols=['Source Name','Comment[ENA_RUN]','Characteristics[age]','Characteristics[genotype]','Characteristics[environmental stress]','Characteristics[organism part]']
 src=sdrf[cols].drop_duplicates();assert len(src)==24
 d=samples.merge(src,left_on='source_name',right_on='Source Name',validate='one_to_one')
 assert len(d)==24 and (d.run==d['Comment[ENA_RUN]']).all()
 assert (d.age==('M'+d['Characteristics[age]'].astype(str))).all()
 assert (d.genotype==d['Characteristics[genotype]'].str.replace('Apo','',regex=False)).all()
 assert (d.exposure==d['Characteristics[environmental stress]'].map({'Filtered Air':'FA','Ozone':'O3'})).all()
 assert (d.tissue==d['Characteristics[organism part]']).all()
 d.to_csv(OUT/'ozone_sample_source_mapping.csv',index=False)
 mirna=pd.read_csv(ROOT/'results/mirna/sample_metadata.tsv',sep='\t')
 soft=(ROOT/'metadata/GSE93967/GSE93967_family.soft').read_text()
 for _,r in mirna.iterrows():
  section=soft.split('^SAMPLE = '+r['sample'],1)[1].split('^SAMPLE = ',1)[0]
  assert ('!Sample_title = '+r.title)in section
 assert len(mirna)==9 and mirna.groupby('condition').size().to_dict()=={'Control':3,'PM1':3,'PM5':3}
 save_json(OUT/'sample_mapping_validation.json',{'ozone_24_SDRF_source_run_age_genotype_exposure_tissue_exact':True,'mirna_9_GEO_SOFT_titles_exact':True,'ozone_SDRF_sha256':sha(ROOT/'metadata/E-MTAB-14917/E-MTAB-14917.sdrf.txt'),'miRNA_SOFT_sha256':sha(ROOT/'metadata/GSE93967/GSE93967_family.soft')})


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--reuse-models',action='store_true');args=parser.parse_args()
 for p in [OUT,META,FIG]:p.mkdir(parents=True,exist_ok=True)
 if not args.reuse_models:subprocess.run(['Rscript',str(ROOT/'scripts/replicate_mirna_ozone.R')],cwd=ROOT,check=True)
 sample_audit()
 author,s3,lists=author_sources();mirna=mirna_statistics();ozone=ozone_comparisons(author,s3,lists)
 for spanish in [False,True]:
  mirna_plot(mirna[0],mirna[1],mirna[3],mirna[4],mirna[5],spanish)
  ozone_plot(ozone[3],ozone[6],ozone[0],spanish);upset_plot(ozone[5],spanish)
 outputs(author,s3,lists,mirna,ozone)
 schemas=[]
 for p in sorted(OUT.iterdir()):
  if p.name.endswith(('.csv','.csv.gz','.tsv')):
   d=pd.read_csv(p,sep='\t'if p.suffix=='.tsv'else',')
   schemas.append({'file':str(p.relative_to(ROOT)),'rows':len(d),'columns':list(d.columns)})
 save_json(OUT/'table_schemas.json',schemas)
 print('Paper reproduction candidates, concordance tables, six bilingual figure families and metadata complete.')

if __name__=='__main__':main()
