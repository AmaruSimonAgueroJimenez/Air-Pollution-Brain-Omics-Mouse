#!/usr/bin/env python3
"""Secondary model tests and publication-style figures from existing brain omics outputs.

No primary output or source data is modified. R inference and descriptive figure
steps can be executed separately so figure revisions never require a model refit.
"""
from __future__ import annotations
import argparse, hashlib, importlib.metadata, json, subprocess, textwrap
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import chi2, false_discovery_control, norm

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/advanced_omics';META=ROOT/'metadata/advanced_omics'
FIG=OUT/'figures';SEED=20261019
KEYS=['M3_E3','M3_E4','M17_E3','M17_E4']
COLORS=['#0072B2','#D55E00','#009E73','#CC79A7'];MARKERS=['o','s','^','D']
LABELS=['M3 · APOE3','M3 · APOE4','M17 · APOE3','M17 · APOE4']
EKEYS=['hypothalamus','cortex','ozone_M3_E3','ozone_M3_E4','ozone_M17_E3','ozone_M17_E4']
ELABELS=['PM · HY\n4 FA / 3 PM','PM · ACC\n4 PBS / 4 PM','O₃ · M3 E3\n3 FA / 3 O₃','O₃ · M3 E4\n3 FA / 3 O₃','O₃ · M17 E3\n3 FA / 3 O₃','O₃ · M17 E4\n3 FA / 3 O₃']
LKEYS=['genotype_dependent_exposure','age_dependent_exposure','three_way_interaction']
WKEYS=['genotype_response_difference_at_M3','genotype_response_difference_at_M17','age_response_difference_in_E3','age_response_difference_in_E4']
MANIFEST={}

def t(en,es,lang):return es if lang=='es' else en

def save(tab,name):
    tab.to_csv(OUT/name,index=False,sep='\t',float_format='%.12g',lineterminator='\n')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def read(path):return pd.read_csv(ROOT/path,sep='\t' if '.tsv' in str(path)else ',')

def panel(ax,label):ax.text(-.10,1.05,label,transform=ax.transAxes,fontweight='bold',fontsize=13,va='bottom')

def export(fig,stem,caption,alt,lang):
    name=stem+('_es' if lang=='es' else '')
    files={}
    for ext in ['png','svg','pdf']:
        p=FIG/f'{name}.{ext}';fig.savefig(p,dpi=300,facecolor='white',bbox_inches='tight',pad_inches=.13);files[ext]={'path':str(p.relative_to(ROOT)),'sha256':sha(p)}
    MANIFEST[name]={'caption':caption,'alt':alt,'language':lang,'files':files,'figure_inches':list(fig.get_size_inches()),'png_dpi':300}
    plt.close(fig)

def bh(values):
    a=np.asarray(values,float);good=np.isfinite(a);out=np.full(len(a),np.nan)
    out[good]=false_discovery_control(a[good],method='bh');return out

def validate_interactions():
    summary=json.loads((OUT/'interaction_summary.json').read_text());checks={}
    all_lrt=read('results/advanced_omics/interaction_LRT_all_tests.csv')
    all_wald=read('results/advanced_omics/interaction_Wald_all_tests.csv')
    checks['all_24_unique_animals_and_8_rank_design']=summary['samples']==24 and all(v['full_rank']==8 and v['nested']for v in summary['designs'].values())
    checks['dispersions_and_effective_length_offsets_preserved']=summary['dispersions_preserved']and summary['normalization_factor_offsets_preserved']
    for key in LKEYS:
        d=read(f'results/advanced_omics/interaction_LRT_{key}.csv');p=d.pvalue.to_numpy();valid=np.isfinite(p)
        checks[key+'_chi_square_tail']=np.allclose(p[valid],chi2.sf(d.loc[valid,'stat'],d.loc[valid,'degrees_of_freedom']),atol=1e-12,rtol=1e-7)
        checks[key+'_BH']=np.allclose(d.padj,bh(p),equal_nan=True,atol=1e-11)
        checks[key+'_summary']=int(d.padj.lt(.05).sum())==summary['LRT'][key]['FDR_005']
        bad=~d.full_beta_converged|~d.reduced_beta_converged
        checks[key+'_convergence_mask']=d.loc[bad,'pvalue'].isna().all()
        checks[key+'_fixed_gene_universe']=len(d)==18537 and d.gene_id.is_unique
    for key in summary['Wald']:
        d=read(f'results/advanced_omics/interaction_Wald_{key}.csv');p=d.pvalue.to_numpy();valid=np.isfinite(p)
        checks[key+'_normal_tail']=np.allclose(p[valid],2*norm.sf(np.abs(d.loc[valid,'stat'])),atol=1e-12,rtol=1e-7)
        checks[key+'_BH']=np.allclose(d.padj,bh(p),equal_nan=True,atol=1e-11)
        checks[key+'_pointwise_interval']=np.allclose(d.CI95_low,d.log2FoldChange-norm.ppf(.975)*d.lfcSE,equal_nan=True,atol=1e-10)and np.allclose(d.CI95_high,d.log2FoldChange+norm.ppf(.975)*d.lfcSE,equal_nan=True,atol=1e-10)
        bad=(~d.beta_converged)|d.one_required_stratum_all_zero
        checks[key+'_validity_mask']=d.loc[bad,['pvalue','log2FoldChange','lfcSE']].isna().all().all()
        checks[key+'_summary']=int(d.padj.lt(.05).sum())==summary['Wald'][key]['FDR_005']
    checks['joint_LRT_BH']=np.allclose(all_lrt.padj_across_three_LRTs,bh(all_lrt.pvalue),equal_nan=True,atol=1e-11)
    checks['joint_Wald_BH']=np.allclose(all_wald.padj_across_five_Wald_tests,bh(all_wald.pvalue),equal_nan=True,atol=1e-11)
    assert all(checks.values()),[k for k,v in checks.items()if not v]
    return {k:bool(v)for k,v in checks.items()}

def select_genes():
    global_lrt=read('results/ozone/de_global_ozone_LRT.csv')
    top=global_lrt.loc[global_lrt.gene_type.eq('protein_coding')&global_lrt.baseMean.ge(20)&global_lrt.padj.notna()].sort_values(['padj','pvalue','gene_id']).head(6)
    selected=top[['gene_id','gene_symbol','baseMean','padj']].copy();selected['selection']='Six lowest global ozone LRT q-values among protein-coding genes with baseMean >=20'
    for symbol in ['Fos','Gfap','Aif1','Trem2']:
        candidate=global_lrt.loc[global_lrt.gene_symbol.eq(symbol)].sort_values('baseMean',ascending=False).head(1)
        if len(candidate)and candidate.gene_id.iloc[0]not in set(selected.gene_id):
            row=candidate[['gene_id','gene_symbol','baseMean','padj']].copy();row['selection']='Named interpretability display target; selected for secondary figure, not a primary confirmatory gene set';selected=pd.concat([selected,row],ignore_index=True)
    save(selected,'ozone_display_gene_selection.tsv')
    age=read('results/advanced_omics/interaction_LRT_age_dependent_exposure.csv')
    age_genes=age.loc[age.padj.lt(.05)].sort_values(['padj','pvalue','gene_id'])
    # Keep all eight discoveries, including weakly expressed Gm21586; no abundance-based deletion.
    save(age_genes,'age_response_FDR_genes.tsv')
    return selected,age_genes

def primary_effects(selected,lang):
    tables=[read(f'results/ozone/de_{key}_O3_vs_FA.csv').set_index('gene_id').loc[selected.gene_id].reset_index()for key in KEYS]
    source=pd.concat([d.assign(display_stratum=k)for d,k in zip(tables,KEYS)],ignore_index=True);save(source,'ozone_four_stratum_effect_display.tsv')
    diagnostics=[]
    for key in KEYS:
        entire=read(f'results/ozone/de_{key}_O3_vs_FA.csv');valid=entire.log2FoldChange.notna()&entire.log2FoldChange_shrunk.notna()
        diagnostics.append({'contrast':key,'genes_with_both_estimates':int(valid.sum()),'median_abs_MLE':float(entire.loc[valid,'log2FoldChange'].abs().median()),'median_abs_apeglm_mode':float(entire.loc[valid,'log2FoldChange_shrunk'].abs().median()),'median_posterior_SD':float(entire.loc[valid,'lfcSE_shrunk'].median()),'FDR_genes_with_abs_MLE_gt_1_and_abs_mode_lt_0_001':int((entire.padj.lt(.05)&entire.log2FoldChange.abs().gt(1)&entire.log2FoldChange_shrunk.abs().lt(.001)).sum())})
    save(pd.DataFrame(diagnostics),'ozone_shrinkage_diagnostics.tsv')
    n=len(selected);fig=plt.figure(figsize=(15.3,7.8),layout='constrained');gs=fig.add_gridspec(1,5,width_ratios=[1.15,1,1,1,1.05])
    limits=np.concatenate([d[['wald_CI95_low','wald_CI95_high']].to_numpy().ravel()for d in tables]);finite=limits[np.isfinite(limits)];lim=max(abs(finite))*1.08
    axes=[]
    for j,(d,label,color,marker)in enumerate(zip(tables,LABELS,COLORS,MARKERS)):
        ax=fig.add_subplot(gs[0,j]);axes.append(ax);ax.axvline(0,color='#65767a',ls='--',lw=.8)
        for i,row in d.iterrows():
            if np.isfinite(row.log2FoldChange)and np.isfinite(row.lfcSE):
                ax.plot([row.wald_CI95_low,row.wald_CI95_high],[i,i],color=color,lw=1.1)
                ax.scatter(row.log2FoldChange,i,s=31,marker=marker,facecolor=color if row.padj<.05 else 'white',edgecolor=color,lw=1,zorder=3)
            else:ax.text(0,i,'×',color='#67767a',ha='center',va='center')
        ax.set(xlim=(-lim,lim),ylim=(n-.4,-.6),yticks=range(n),title=label,xlabel=t('MLE log₂FC ± 95% CI','MLE log₂FC ± IC del 95%',lang))
        ax.set_yticklabels(selected.gene_symbol if j==0 else []);ax.grid(axis='x',alpha=.16)
    panel(axes[0],'A')
    ax=fig.add_subplot(gs[0,4]);mat=np.column_stack([d.log2FoldChange_shrunk for d in tables]);v=np.nanmax(abs(mat));cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e3e7e8')
    im=ax.imshow(mat,aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=-v,vcenter=0,vmax=v));ax.set(xticks=range(4),yticks=range(n));ax.set_yticklabels([]);ax.set_xticklabels(KEYS,rotation=45,ha='right');ax.set_title(t('Regularized effects','Efectos regularizados',lang));panel(ax,'B')
    fig.colorbar(im,ax=ax,shrink=.65,label=t('apeglm log₂FC','log₂FC de apeglm',lang))
    fig.suptitle(t('Ozone responses: effect size, uncertainty and regularization','Respuestas al ozono: tamaño del efecto, incertidumbre y regularización',lang),fontsize=15,fontweight='bold')
    handles=[Line2D([0],[0],marker='o',color='#52686e',ls='',markerfacecolor='#52686e',label=t('Within-stratum BH FDR <0.05','FDR de BH del estrato <0.05',lang)),Line2D([0],[0],marker='o',color='#52686e',ls='',markerfacecolor='white',label=t('FDR ≥0.05','FDR ≥0.05',lang))]
    fig.legend(handles=handles,loc='outside lower center',ncol=2,frameon=False)
    caption=t('A: Unshrunken ozone-versus-air log2 fold changes and pointwise95% Wald intervals;3 animals per exposure in each stratum, with dispersion estimation shared across24 animals. Filled symbols mark within-stratum BH FDR<0.05; intervals are not simultaneous. B: apeglm-regularized effects on the same genes. Gray/missing marks mean unavailable estimates, not zero. Six genes were selected by the global ozone LRT among protein-coding genes with mean normalized count>=20; Fos,Gfap,Aif1,Trem2 are named secondary display targets. Selection is post hoc and does not establish gene-specific replication.',
    'A: Cambios log2 de ozono frente a aire sin regularizar e intervalos puntuales de Wald del95%;3 animales por exposición en cada estrato, con dispersión estimada en24 animales. Símbolos rellenos:FDR de BH del estrato<0.05;los intervalos no son simultáneos. B:Efectos regularizados por apeglm de los mismos genes. Gris o ausente indica estimaciones no disponibles, no cero. Seis genes se seleccionaron por la LRT global entre genes codificantes con conteo medio>=20;Fos,Gfap,Aif1,Trem2 son objetivos secundarios de visualización. La selección es posterior al análisis y no establece replicación génica.',lang)
    caption+=t(' Regularization is estimated separately for each contrast and can contract a mode almost to zero despite a significant unshrunken Wald result;for example,Fos in M17/APOE4 has MLE=-2.271,q=0.0212,but mode≈0. The matrix is a regularization diagnostic,not evidence that the original exposure effect is absent.',' La regularización se estima por separado en cada contraste y puede llevar el modo casi a cero pese a una prueba de Wald significativa sin regularizar;por ejemplo,Fos en M17/APOE4 tiene MLE=-2.271,q=0.0212,pero modo≈0. La matriz diagnostica la regularización,no demuestra ausencia del efecto original de exposición.',lang)
    export(fig,'ozone_effect_uncertainty',caption,t('Four forest plots and a regularized-effect heatmap for ten selected ozone-study genes.','Cuatro gráficos de efectos e intervalos y un mapa de efectos regularizados de diez genes seleccionados del estudio de ozono.',lang),lang)

def interaction_evidence(age_genes,lang):
    fig=plt.figure(figsize=(14.5,10),layout='constrained');gs=fig.add_gridspec(2,6,height_ratios=[1,1.2]);tables=[read(f'results/advanced_omics/interaction_LRT_{k}.csv')for k in LKEYS]
    names=[t('Genotype-dependent response ·2df','Respuesta dependiente del genotipo ·2gl',lang),t('Age-dependent response ·2df','Respuesta dependiente de la edad ·2gl',lang),t('Three-way interaction ·1df','Interacción triple ·1gl',lang)]
    ymax=max(-np.log10(d.padj.clip(lower=1e-300)).max()for d in tables)*1.14
    for j,(d,name)in enumerate(zip(tables,names)):
        ax=fig.add_subplot(gs[0,2*j:2*j+2]);good=d.padj.notna();sig=d.padj.lt(.05)
        ax.scatter(np.log10(1+d.loc[good&~sig,'baseMean']),-np.log10(d.loc[good&~sig,'padj']),s=4,c='#9eacb0',alpha=.35,rasterized=True)
        ax.scatter(np.log10(1+d.loc[sig,'baseMean']),-np.log10(d.loc[sig,'padj']),s=27,c='#9a3d2e',marker='D',edgecolors='white',lw=.3,zorder=3)
        ax.axhline(-np.log10(.05),color='#7d4b3a',ls='--',lw=.8)
        ax.set(title=name,xlabel=t('log₁₀(1+mean normalized count)','log₁₀(1+conteo normalizado medio)',lang),ylabel='−log₁₀(BH FDR)',ylim=(-.02,max(1.55,ymax)))
        ax.text(.03,.97,f"{int(sig.sum())} / {int(good.sum()):,}"+t(' genes\nFDR <0.05',' genes\nFDR <0.05',lang),transform=ax.transAxes,va='top',fontsize=9)
        if j==0:panel(ax,'A')
    ax=fig.add_subplot(gs[1,:2]);counts=[int(d.padj.lt(.05).sum())for d in tables];bars=ax.barh(range(3),counts,color=['#536f79','#0072B2','#95819d'],height=.58)
    ax.set(yticks=range(3),yticklabels=[t('Genotype × exposure','Genotipo × exposición',lang),t('Age × exposure','Edad × exposición',lang),t('Age × genotype × exposure','Edad × genotipo × exposición',lang)],xlabel=t('Genes with test-specific BH FDR <0.05','Genes con FDR de BH de cada prueba <0.05',lang),xlim=(0,max(counts+[1])*1.25));ax.invert_yaxis()
    for bar,c in zip(bars,counts):ax.text(c+.12,bar.get_y()+bar.get_height()/2,str(c),va='center',fontweight='bold')
    ax.set_title(t('Secondary omnibus discoveries','Hallazgos ómnibus secundarios',lang));panel(ax,'B');ax.spines[['top','right']].set_visible(False)
    ax=fig.add_subplot(gs[1,2:]);ids=age_genes.gene_id.to_list()
    mats=[];qs=[];source=[]
    for key in WKEYS:
        d=read(f'results/advanced_omics/interaction_Wald_{key}.csv').set_index('gene_id').loc[ids];mats.append(d.log2FoldChange);qs.append(d.padj);source.append(d.reset_index())
    matrix=np.column_stack(mats);q=np.column_stack(qs);lim=np.nanmax(abs(matrix));cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e3e7e8')
    im=ax.imshow(matrix,aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim))
    for i,j in zip(*np.where(q<.05)):ax.plot(j,i,'k.',ms=5)
    ax.set(yticks=range(len(ids)),yticklabels=age_genes.gene_symbol,xticks=range(4));ax.set_xticklabels(['ΔE4−E3\nM3','ΔE4−E3\nM17','ΔM17−M3\nE3','ΔM17−M3\nE4']);ax.set_title(t('Signed differences in ozone response','Diferencias con signo en la respuesta al ozono',lang));panel(ax,'C');fig.colorbar(im,ax=ax,shrink=.8,label='Δ log₂FC')
    fig.suptitle(t('Ozone response heterogeneity: formal interaction tests','Heterogeneidad de respuesta al ozono: pruebas formales de interacción',lang),fontsize=15,fontweight='bold')
    save(pd.concat(source,ignore_index=True),'interaction_evidence_heatmap_source.tsv')
    export(fig,'ozone_interaction_evidence',t('A: Fixed-universe nested-model likelihood-ratio tests;each gene is one point. A separate BH family is used for each test. All24 animals and effective-length offsets are retained;the full/reduced ranks are8/6,8/6,and8/7. B: Eight age-dependent-response genes, zero genotype-omnibus and zero three-way genes meet BH FDR<0.05. A null test is not equivalence. C: Unshrunken differences of ozone log2 fold changes for all eight age-LRT discoveries;black dots indicate the corresponding conditional Wald-family FDR<0.05. These one-df conditional tests and omnibus tests answer different hypotheses;their discovery sets need not agree. All analyses are secondary and exploratory.',
    'A:LRT de modelos anidados con universo génico fijo;cada gen es un punto. Cada prueba usa una familia de BH separada. Se conservan24 animales y los ajustes por longitud efectiva;los rangos completo/reducido son8/6,8/6 y8/7. B:Ocho genes muestran respuesta dependiente de la edad;ninguno supera FDR<0.05 en la prueba ómnibus de genotipo ni en la interacción triple. Un resultado nulo no implica equivalencia. C:Diferencias sin regularizar de los cambios log2 por ozono en los ocho hallazgos de edad;los puntos negros indican FDR<0.05 en la familia condicional de Wald correspondiente. Las pruebas condicionales de1gl y las ómnibus responden hipótesis distintas;no requieren los mismos hallazgos. Son análisis secundarios exploratorios.',lang),t('Three interaction evidence plots, discovery counts and a signed response-difference heatmap.','Tres gráficos de evidencia de interacción, conteos de hallazgos y mapa de diferencias de respuesta con signo.',lang),lang)

def animal_profiles(age_genes,lang):
    selected=age_genes.head(6);counts=read('results/ozone/normalized_counts.csv.gz').set_index('gene_id');meta=read('results/ozone/samples.csv');rng=np.random.default_rng(SEED)
    fig,axes=plt.subplots(2,3,figsize=(14.5,8.5),layout='constrained');sources=[]
    for panel_index,((_,gene),ax)in enumerate(zip(selected.iterrows(),axes.flat)):
        for j,(key,color,marker)in enumerate(zip(KEYS,COLORS,MARKERS)):
            age,gt=key.split('_');sub=meta.loc[meta.age.eq(age)&meta.genotype.eq(gt)]
            means=[]
            for k,exp in enumerate(['FA','O3']):
                rows=sub.loc[sub.exposure.eq(exp)];vals=counts.loc[gene.gene_id,rows.sample_id].to_numpy(float);trans=np.log2(1+vals);x=j*2.8+k
                ax.scatter(x+rng.uniform(-.11,.11,len(vals)),trans,color=color,marker=marker,s=30,edgecolors='white',lw=.4,zorder=3)
                means.append(float(trans.mean()))
                for sid,val,z in zip(rows.sample_id,vals,trans):sources.append({'gene_id':gene.gene_id,'gene_symbol':gene.gene_symbol,'sample_id':sid,'stratum':key,'exposure':exp,'normalized_count':val,'log2_1_plus_normalized_count':z})
            ax.plot([j*2.8,j*2.8+1],means,color=color,lw=1.4,zorder=2);ax.plot([j*2.8-.16,j*2.8+.16],[means[0]]*2,color=color,lw=2);ax.plot([j*2.8+.84,j*2.8+1.16],[means[1]]*2,color=color,lw=2)
        ticks=[j*2.8+k for j in range(4)for k in range(2)]
        ax.set(xticks=ticks,xticklabels=['FA','O₃']*4,ylabel=t('log₂(1+normalized count)','log₂(1+conteo normalizado)',lang),title=f"{gene.gene_symbol} · "+t('age-response q=','q de respuesta por edad=',lang)+f'{gene.padj:.3g}')
        ax.margins(y=.18)
        for j,key in enumerate(KEYS):ax.text(j*2.8+.5,-.15,key.replace('_','\n'),transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8)
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);panel(ax,chr(65+panel_index))
    fig.suptitle(t('Individual animals underlying the strongest age-response findings','Animales individuales tras los principales hallazgos de respuesta por edad',lang),fontsize=15,fontweight='bold')
    save(pd.DataFrame(sources),'ozone_interaction_animal_profiles.tsv')
    export(fig,'ozone_interaction_animals',t('All24 independent animals per gene are shown,3 per exposure/age/genotype cell. Counts use the original DESeq2 gene-length normalization factors, then log2(1+count) for display. Point-position plots use gene-specific nonzero y-ranges,identical across groups within a panel. Connected bars represent arithmetic means of displayed transformed values;they do not connect paired animals and are not inferential intervals. Six lowest-q age-dependent-exposure LRT genes are displayed, preserving secondary selection. Count-model tests were not fitted on these log-transformed values. Gene names and grouping do not establish cell-specific expression.',
    'Se muestran los24 animales independientes por gen,3 por combinación exposición/edad/genotipo. Los conteos usan los factores originales de normalización por longitud génica de DESeq2 y después log2(1+conteo) para visualizar. Los puntos usan rangos verticales específicos por gen que no comienzan en cero,iguales entre grupos dentro de cada panel. Las barras conectadas son medias aritméticas de los valores transformados;no conectan animales pareados ni son intervalos inferenciales. Se muestran los seis menores valores q de la LRT de respuesta dependiente de la edad,como selección secundaria. Las pruebas del modelo de conteos no se ajustaron sobre estos valores logarítmicos. Los nombres génicos y los grupos no establecen expresión específica de células.',lang),t('Six panels showing all24 individual animal values for selected age-response genes.','Seis paneles con los valores de los24 animales individuales en genes seleccionados de respuesta por edad.',lang),lang)

def gmt(path):
    result={}
    for line in path.read_text().splitlines():
        fields=line.split('\t');result[fields[0]]=set(fields[2:])
    return result

def go_selection():
    tables={k:read(f'results/enrichment/{k}/GSEA_all_terms.tsv')for k in EKEYS};sets={}
    for k in EKEYS:
        for go,genes in gmt(ROOT/f'results/enrichment/{k}/GO_BP_tested_gene_sets.gmt').items():sets.setdefault(go,set()).update(genes)
    selected=[];origins=[]
    for key in EKEYS:
        d=tables[key]
        for sign in [1,-1]:
            candidate=d.loc[d.fdr_qvalue.lt(.05)&d.NES.mul(sign).gt(0)].assign(absNES=lambda x:x.NES.abs()).sort_values(['fdr_qvalue','absNES','GO_ID'],ascending=[True,False,True])
            added=0
            for _,row in candidate.iterrows():
                go=row.GO_ID
                if go in selected:continue
                if any(len(sets[go]&sets[other])/len(sets[go]|sets[other])>=.75 for other in selected):continue
                selected.append(go);origins.append({'GO_ID':go,'term_name':row.term_name,'selection_contrast':key,'selection_direction':sign,'selection_q':row.fdr_qvalue,'selection_NES':row.NES});added+=1
                if added==2:break
    save(pd.DataFrame(origins),'GO_heatmap_term_selection.tsv')
    long=[]
    for go in selected:
        for key in EKEYS:
            rows=tables[key].loc[tables[key].GO_ID.eq(go)]
            if len(rows):
                row=rows.iloc[0];long.append({'GO_ID':go,'contrast':key,'term_name':row.term_name,'NES':row.NES,'algorithmic_FDR':row.fdr_qvalue,'eligible_tested':True,'rank_member_genes':row.rank_member_genes})
            else:long.append({'GO_ID':go,'contrast':key,'term_name':next(x['term_name']for x in origins if x['GO_ID']==go),'NES':np.nan,'algorithmic_FDR':np.nan,'eligible_tested':False,'rank_member_genes':np.nan})
    save(pd.DataFrame(long),'GO_multi_contrast_NES_source.tsv')
    return tables,sets,selected,pd.DataFrame(long)

def go_heatmap(selected,source,lang):
    mat=source.pivot(index='GO_ID',columns='contrast',values='NES').reindex(index=selected,columns=EKEYS);q=source.pivot(index='GO_ID',columns='contrast',values='algorithmic_FDR').reindex_like(mat)
    labels=source.drop_duplicates('GO_ID').set_index('GO_ID').term_name;fig,ax=plt.subplots(figsize=(13.8,max(8.8,len(selected)*.34)),layout='constrained');cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dce2e4');lim=np.nanmax(abs(mat.to_numpy()))
    im=ax.imshow(mat.to_numpy(),aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim))
    ax.set(xticks=range(6),yticks=range(len(selected)));ax.set_xticklabels(ELABELS,fontsize=9);ax.set_yticklabels([textwrap.fill(labels[g],51)+'  '+g for g in selected],fontsize=8)
    for i,j in zip(*np.where(q.to_numpy()<.05)):ax.plot(j,i,marker='o',ms=3.4,markerfacecolor='white',markeredgecolor='black',markeredgewidth=.6)
    for i,j in zip(*np.where(~np.isfinite(mat.to_numpy()))):ax.text(j,i,'×',ha='center',va='center',color='#586a70')
    ax.set_xticks(np.arange(-.5,6,1),minor=True);ax.set_yticks(np.arange(-.5,len(selected),1),minor=True);ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',bottom=False,left=False)
    ax.set_title(t('GO Biological Process associations across six distinct contrasts','Asociaciones de procesos biológicos GO en seis contrastes distintos',lang),fontsize=14,fontweight='bold',pad=16);fig.colorbar(im,ax=ax,shrink=.75,label=t('Normalized enrichment score (NES)','Puntuación de enriquecimiento normalizada (NES)',lang))
    export(fig,'GO_multi_contrast_heatmap',t('Signed NES values are displayed for a deterministic, redundancy-reduced selection:up to2 FDR<0.05 terms per direction per contrast, ordered by estimated q then absolute NES;new terms require Jaccard<0.75 with already selected terms using the union of rank-eligible GO members. Outlined circles mark that cell\'s GSEA algorithmic FDR<0.05;these are not BH p-values. Gray× means the term was not eligible/tested, never zero enrichment. Gene-set permutations do not preserve gene-gene correlation or sample uncertainty. Contrasts have distinct tissues/models/assays and are not pooled;differences in NES or q do not constitute a formal between-study test. Source GO names are retained.',
    'Se muestran NES con signo en una selección determinista con reducción de redundancia:hasta2 términos con FDR<0.05 por dirección y contraste, ordenados por q estimado y |NES|;un término nuevo requiere Jaccard<0.75 respecto de los ya elegidos, usando la unión de miembros elegibles. Los círculos delineados indican FDR algorítmica de GSEA<0.05 en esa celda;no son valores p de BH. Gris× significa no elegible/no analizado, nunca enriquecimiento cero. Las permutaciones de conjuntos no conservan la correlación génica ni la incertidumbre muestral. Los contrastes tienen tejidos/modelos/ensayos distintos y no se combinan;diferencias de NES o q no son pruebas formales entre estudios. Se conservan los nombres originales de GO.',lang),t('A signed NES matrix across six contrasts with explicit missing eligibility and algorithmic FDR markers.','Matriz de NES con signo en seis contrastes, con elegibilidad ausente explícita y marcas de FDR algorítmica.',lang),lang)

def go_network(tables,lang):
    key='ozone_M3_E4';d=tables[key];members=gmt(ROOT/f'results/enrichment/{key}/GO_BP_tested_gene_sets.gmt');parts=[]
    for sign in [1,-1]:parts.append(d.loc[d.fdr_qvalue.lt(.05)&d.NES.mul(sign).gt(0)].assign(absNES=lambda x:x.NES.abs()).sort_values(['fdr_qvalue','absNES','GO_ID'],ascending=[True,False,True]).head(10))
    chosen=pd.concat(parts,ignore_index=True);graph=nx.Graph();nodes=[];edges=[]
    for i,row in chosen.iterrows():graph.add_node(row.GO_ID);nodes.append({'node':i+1,'GO_ID':row.GO_ID,'term_name':row.term_name,'NES':row.NES,'algorithmic_FDR':row.fdr_qvalue,'rank_member_genes':len(members[row.GO_ID])})
    for i,a in enumerate(chosen.GO_ID):
        for b in chosen.GO_ID.iloc[i+1:]:
            overlap=len(members[a]&members[b]);union=len(members[a]|members[b]);jaccard=overlap/union
            if jaccard>=.25:graph.add_edge(a,b,weight=jaccard);edges.append({'source':a,'target':b,'overlap_genes':overlap,'union_genes':union,'Jaccard':jaccard})
    save(pd.DataFrame(nodes),'GO_overlap_network_nodes.tsv');save(pd.DataFrame(edges,columns=['source','target','overlap_genes','union_genes','Jaccard']),'GO_overlap_network_edges.tsv')
    # Layout connected components separately: no node/label occlusion, no claim
    # that visual distance is molecular distance or an inferred interaction.
    index={go:i for i,go in enumerate(chosen.GO_ID)}
    components=sorted(nx.connected_components(graph),key=lambda c:(-len(c),min(index[g]for g in c)));pos={}
    for j,component in enumerate(components):
        center=np.array([(j%3)*3.,-(j//3)*3.])
        ordered=sorted(component,key=index.get)
        for k,go in enumerate(ordered):
            angle=2*np.pi*k/len(ordered)
            pos[go]=center if len(ordered)==1 else center+.82*np.array([np.cos(angle),np.sin(angle)])
    fig=plt.figure(figsize=(15.5,9.8),layout='constrained');gs=fig.add_gridspec(1,2,width_ratios=[1,.95]);ax=fig.add_subplot(gs[0]);table=fig.add_subplot(gs[1]);table.axis('off')
    vals=chosen.NES.to_numpy();lim=max(abs(vals));normer=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim)
    nx.draw_networkx_edges(graph,pos,ax=ax,width=[.5+4*graph.edges[e]['weight']for e in graph.edges],edge_color='#9aa8ad',alpha=.65)
    coll=nx.draw_networkx_nodes(graph,pos,ax=ax,nodelist=chosen.GO_ID.to_list(),node_size=4*chosen.rank_member_genes.to_numpy(),node_color=vals,cmap='RdBu_r',vmin=-lim,vmax=lim,edgecolors='#334c54',linewidths=.7)
    nx.draw_networkx_labels(graph,pos,{row['GO_ID']:str(row['node'])for row in nodes},ax=ax,font_size=8,font_weight='bold',font_color='white');ax.axis('off');ax.margins(.16)
    ax.set_title(t('A  Shared GO membership · not gene regulation','A  Miembros GO compartidos · no regulación génica',lang),fontsize=12,fontweight='bold')
    table.set_title(t('B  Selected ozone M3/APOE4 terms','B  Términos seleccionados de ozono M3/APOE4',lang),fontsize=12,fontweight='bold',loc='left')
    y=.98;lineheight=.044
    for row in nodes:
        label=f"{row['node']:>2}. "+textwrap.shorten(row['term_name'],width=58,placeholder='…')
        table.text(0,y,label,transform=table.transAxes,fontsize=8.6,va='top');q='0*'if row['algorithmic_FDR']==0 else f"{row['algorithmic_FDR']:.2g}"
        table.text(.02,y-.018,f"{row['GO_ID']}   NES={row['NES']:.2f}   q={q}   n={row['rank_member_genes']}",transform=table.transAxes,fontsize=7,color='#586b72',va='top');y-=lineheight
    table.text(0,.02,t('Node area ∝ set members in the rank\nEdges: Jaccard ≥0.25; width increases with overlap\n*0 is a finite-permutation estimate, not probability zero.','Área del nodo ∝ miembros del conjunto en la lista\nAristas:Jaccard ≥0.25;grosor creciente con solapamiento\n*0 es una estimación de permutaciones finitas, no probabilidad cero.',lang),transform=table.transAxes,fontsize=8,va='bottom',linespacing=1.5)
    fig.colorbar(plt.cm.ScalarMappable(norm=normer,cmap='RdBu_r'),ax=ax,shrink=.65,label='NES')
    fig.suptitle(t('Enrichment redundancy: twenty data-selected GO terms','Redundancia del enriquecimiento: veinte términos GO seleccionados de los datos',lang),fontsize=15,fontweight='bold')
    export(fig,'GO_overlap_network',t('Ten lowest-q positive and ten lowest-q negative GO-BP terms for ozone M3/APOE4, with |NES| tie-breaks. Nodes encode gene-set size in the eligible rank and NES, not cell abundance. Edges join sets with Jaccard membership overlap>=0.25;they do not represent coexpression, protein interaction, causal regulation or independent pathway discoveries. Connected components use a deterministic grid and circular positions; visual distances are graphical, not anatomical or molecular. Selection is exploratory and restricted to one contrast;the original algorithmic GSEA FDR estimates are retained. Zero q values reflect finite permutations. Full names, memberships, node and edge tables are available.',
    'Diez términos GO-BP de NES positivo y diez de NES negativo con menor q para ozono M3/APOE4, desempatados por |NES|. Los nodos representan tamaño del conjunto en la lista elegible y NES, no abundancia celular. Las aristas unen conjuntos con Jaccard>=0.25;no representan coexpresión, interacción proteica, regulación causal ni vías independientes. Los componentes se distribuyen en una cuadrícula determinista con posiciones circulares; las distancias son gráficas, no anatómicas ni moleculares. La selección es exploratoria y corresponde a un único contraste;se conservan las estimaciones originales de FDR algorítmica. Los q iguales a cero reflejan permutaciones finitas. Se proporcionan nombres completos, miembros y tablas de nodos/aristas.',lang),t('A membership-overlap network of twenty selected GO terms with a numbered term table.','Red de solapamiento de miembros de veinte términos GO seleccionados con tabla numerada.',lang),lang)
    return {'contrast':key,'nodes':len(nodes),'edges':len(edges),'Jaccard_threshold':.25,'layout':'Deterministic connected-component grid with circular within-component positions'}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--refit',action='store_true',help='Run the secondary R interaction models before figures');args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(exist_ok=True);META.mkdir(parents=True,exist_ok=True)
    if args.refit:subprocess.run(['Rscript',str(ROOT/'scripts/analyze_advanced_omics.R')],cwd=ROOT,check=True)
    if not(OUT/'interaction_summary.json').exists():raise RuntimeError('Run with --refit once or execute analyze_advanced_omics.R first.')
    checks=validate_interactions();selected,age_genes=select_genes();tables,sets,terms,source=go_selection()
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False}):
        for lang in ['en','es']:
            primary_effects(selected,lang);interaction_evidence(age_genes,lang);animal_profiles(age_genes,lang);go_heatmap(terms,source,lang);network=go_network(tables,lang)
    summary=json.loads((OUT/'interaction_summary.json').read_text())
    pooled_lrt=read('results/advanced_omics/interaction_LRT_all_tests.csv');pooled_wald=read('results/advanced_omics/interaction_Wald_all_tests.csv')
    for key,item in summary['LRT'].items():item['FDR_005_across_three_LRTs']=int((pooled_lrt.test.eq(key)&pooled_lrt.padj_across_three_LRTs.lt(.05)).sum())
    for key,item in summary['Wald'].items():item['FDR_005_across_five_Wald_tests']=int((pooled_wald.test.eq(key)&pooled_wald.padj_across_five_Wald_tests.lt(.05)).sum())
    summary['figure_families']=5;summary['figure_language_versions']=2;summary['GO_heatmap_terms']=len(terms);summary['GO_overlap_network']=network;summary['display_genes']=selected.gene_symbol.tolist();summary['age_response_gene_symbols']=age_genes.gene_symbol.tolist();summary['complete']=True
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(OUT/'figure_manifest.json').write_text(json.dumps(MANIFEST,indent=2,ensure_ascii=False)+'\n')
    inputs=[ROOT/'results/ozone/dds_full_wald.rds',ROOT/'results/ozone/normalized_counts.csv.gz',ROOT/'results/ozone/de_global_ozone_LRT.csv']+[ROOT/f'results/ozone/de_{k}_O3_vs_FA.csv'for k in KEYS]+[ROOT/f'results/enrichment/{k}/GSEA_all_terms.tsv'for k in EKEYS]
    methods={'created_at_UTC':datetime.now(timezone.utc).isoformat(),'secondary_exploratory':True,'source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'code_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__).resolve(),ROOT/'scripts/analyze_advanced_omics.R']},
    'methods':'Three nested ozone exposure-heterogeneity LRTs and five signed covariance-aware Wald response-difference tests; all24 animals and original18537-gene universe. Effective-length normalization offsets and fitted dispersions preserved; no source re-quantification. Cook filter on, independent filtering off, numerical validity masks before BH. Conditional signed contrasts additionally exclude all-zero required strata.',
    'multiplicity':'BH within each LRT and each Wald contrast; separate supplementary BH across all three LRT families/all five Wald families. Discovery sets across families are not required to be nested.',
    'intervals':'Unshrunken pointwise95% Wald intervals MLE±1.9599639845SE; not simultaneous; not cell-level. Primary-effect forest uses the existing stratum tables. Regularized-effect matrix uses apeglm modes only, not changed hypothesis tests.',
    'GO_display_selection':'Up to2 significant terms per sign per contrast; q then |NES|, skips terms at Jaccard>=0.75 to already selected sets, using union of rank-eligible members. Missing term eligibility remains NA. Original GO term names retained.',
    'GO_network':'20 terms from ozone M3/APOE4(10 per NES sign), ordered by q,then|NES|;edges Jaccard>=0.25 of rank-eligible gene-set membership;node area proportional set size;deterministic component-grid/circular layout. Not regulation/coexpression/anatomy.',
    'animal_profiles':'All24 independent animals for each of6 selected age-LRT genes;log2(1+original DESeq2 normalized count), reproducible tiny horizontal jitter, connected means of displayed values, no pairing or inferential errorbars.',
    'graphics':'PNG300dpi and editable SVG/vector PDF; dense MA points rasterized in SVG/PDF, preserving all points. Both English and Spanish labels; original taxonomy/GO names retained.',
    'documentation':[{'url':'https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions','accessed_date':'2026-10-19','used_for':'Interaction coefficients and nested-model LRT interpretation; installed DESeq2 version is separately pinned.'}],
    'versions':{x:importlib.metadata.version(x)for x in ['numpy','pandas','matplotlib','scipy','networkx']}}
    (META/'methods.json').write_text(json.dumps(methods,indent=2)+'\n')
    validation={'passed':all(checks.values()),'scope':'Independent SciPy p-value/BH/interval arithmetic and documented model invariants; figure visual review separately recorded. Not biological replication.','checks':checks,'figure_exports':len(MANIFEST)*3}
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps({'complete':True,'LRT':summary['LRT'],'GO_terms':len(terms),'network':network,'figure_exports':len(MANIFEST)*3},indent=2))

if __name__=='__main__':main()
