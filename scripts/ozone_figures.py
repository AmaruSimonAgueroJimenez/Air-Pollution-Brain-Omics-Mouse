#!/usr/bin/env python3
"""Export source-backed ozone RNA-seq QC and DE figures, SVG/PDF/PNG."""
from pathlib import Path
import json, platform,subprocess,sys,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
import adjustText
from adjustText import adjust_text

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ozone';FIG=OUT/'figures'
FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
 'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
COLORS={'M3_E3':'#0072B2','M3_E4':'#009E73','M17_E3':'#D55E00','M17_E4':'#CC79A7'}
captions={}

def save(fig,name,caption):
    for ext in ('png','svg','pdf'):fig.savefig(FIG/f'{name}.{ext}',dpi=300)
    plt.close(fig);captions[name]=caption

def quantification_qc(samples):
    records=[];full_records=[]
    txmap=pd.read_csv(ROOT/'data/references/ozone/tx2gene.tsv',sep='\t')
    valid_transcripts=set(txmap.transcript_id)
    for s in samples.to_dict('records'):
        d=json.loads((OUT/'quant'/s['sample_id']/'aux_info/meta_info.json').read_text())
        checkpoint=OUT/'quant'/s['sample_id']/'ozone_quant_complete.json'
        receipt=json.loads(checkpoint.read_text())
        if 'elapsed_pipeline_time' in receipt:
            receipt['completed_at_unix']=receipt.pop('elapsed_pipeline_time')
        receipt['salmon_reported_elapsed_seconds']=d.get('total_time_seconds')
        receipt['metadata_schema']=2
        checkpoint.write_text(json.dumps(receipt,indent=2))
        assert d.get('meta_info_complete') is True and d.get('em_converged') is True
        assert not d.get('quant_errors')
        sf=pd.read_csv(OUT/'quant'/s['sample_id']/'quant.sf',sep='\t')
        assert sf.Name.is_unique and set(sf.Name).issubset(valid_transcripts)
        assert np.isfinite(sf[['Length','EffectiveLength','TPM','NumReads']].to_numpy()).all()
        assert (sf[['Length','EffectiveLength','TPM','NumReads']].to_numpy()>=0).all()
        assert sf.NumReads.sum()>0 and np.isclose(sf.TPM.sum(),1e6,rtol=1e-5)
        assert np.isclose(sf.NumReads.sum(),d['num_mapped'],rtol=1e-6,atol=1)
        record={**s,'num_processed_fragments':d.get('num_processed'),
          'num_mapped_fragments':d.get('num_mapped'),'percent_mapped':d.get('percent_mapped'),
          'library_types':';'.join(d.get('library_types',[])),
          'num_decoy_fragments':d.get('num_decoy_fragments'),
          'num_targets':d.get('num_targets'),'salmon_version':d.get('salmon_version'),
          'em_converged':d.get('em_converged'),'em_iterations':d.get('num_em_iterations'),
          'sum_estimated_transcript_fragments':float(sf.NumReads.sum()),'sum_tpm':float(sf.TPM.sum()),
          'tx2gene_coverage':1.0,'quantification_integrity_passed':True}
        fp=OUT/'trimmed'/s['sample_id']/'fastp.json'
        if fp.exists():
            f=json.loads(fp.read_text());before=f['summary']['before_filtering'];after=f['summary']['after_filtering']
            assert after['total_reads']==2*d['num_processed']
            archive=OUT/'full_qc';archive.mkdir(parents=True,exist_ok=True)
            shutil.copy2(fp,archive/f"{s['sample_id']}.json")
            if fp.with_suffix('.html').exists():shutil.copy2(fp.with_suffix('.html'),archive/f"{s['sample_id']}.html")
            record.update({'raw_pairs_full_qc':before['total_reads']//2,
              'trimmed_pairs_full_qc':after['total_reads']//2,'raw_q30_full_qc':before['q30_rate'],
              'adapter_read_fraction_full_qc':f.get('adapter_cutting',{}).get('adapter_trimmed_reads',0)/before['total_reads']})
            adapters=f.get('adapter_cutting',{})
            full_record={k:s[k]for k in ['sample_id','run','age','genotype','exposure']}
            for prefix,metric in [('raw',before),('after_adapter',after)]:
                full_record.update({prefix+'_pairs':metric['total_reads']//2,
                  prefix+'_bases':metric['total_bases'],prefix+'_q20_rate':metric['q20_rate'],
                  prefix+'_q30_rate':metric['q30_rate'],prefix+'_gc_fraction':metric['gc_content'],
                  prefix+'_r1_mean_length':metric['read1_mean_length'],
                  prefix+'_r2_mean_length':metric.get('read2_mean_length')})
            full_record.update({'adapter_trimmed_reads':adapters.get('adapter_trimmed_reads',0),
              'adapter_trimmed_bases':adapters.get('adapter_trimmed_bases',0),
              'adapter_trimmed_read_fraction':record['adapter_read_fraction_full_qc'],
              'removed_pairs_after_adapter_processing':(before['total_reads']-after['total_reads'])//2,
              'qc_scope':'all original read pairs; no sampling'})
            full_records.append(full_record)
        records.append(record)
    qc=pd.DataFrame(records);qc.to_csv(OUT/'quantification_qc.csv',index=False)
    assert len(full_records)==len(samples)
    pd.DataFrame(full_records).to_csv(OUT/'full_qc_summary.csv',index=False)
    return qc

def main():
    samples=pd.read_csv(OUT/'samples.csv');samples['stratum']=samples.age+'_'+samples.genotype
    raw=pd.read_csv(OUT/'raw_qc_summary.csv');quant=quantification_qc(samples)
    q=samples[['sample_id','stratum','exposure']].merge(raw,on='sample_id',suffixes=('','_qc')).merge(
       quant[['sample_id','percent_mapped','num_processed_fragments']],on='sample_id')
    colors=[COLORS[x] for x in q.stratum]
    fig,axs=plt.subplots(2,2,figsize=(13,7),layout='constrained')
    metrics=[('q30_rate',100,'Q30 in sampled raw-read QC (%)'),
      ('adapter_trimmed_read_fraction',100,'Reads with adapters in sampled QC (%)'),
      ('percent_mapped',1,'Fragments assigned by Salmon (%)'),
      ('num_processed_fragments',1e-6,'Pairs processed by Salmon (millions)')]
    for ax,(field,mult,label) in zip(axs.flat,metrics):
        ax.bar(np.arange(len(q)),q[field]*mult,color=colors,edgecolor='black',linewidth=.3)
        ax.set(xticks=np.arange(len(q)),xticklabels=q.sample_id,ylabel=label)
        ax.tick_params(axis='x',labelrotation=60,labelsize=8);ax.set_ylim(bottom=0)
        if field in ('q30_rate','percent_mapped'):ax.set_ylim(0,100)
    fig.suptitle('O₃ · Hippocampus · 24 animals, 3 per experimental combination',fontsize=14)
    fig.legend(handles=[Line2D([0],[0],marker='s',linestyle='',color=c,
      label=k.replace('M17','17 months').replace('M3','3 months').replace('_',' · '))for k,c in COLORS.items()],
      loc='outside lower center',ncol=4,frameon=False,fontsize=9)
    save(fig,'quality_control','Quality of the first 200,000 pairs per library and quantification of all pairs. Colors denote starting age/genotype; each bar is one animal. raw_qc_summary.csv and quantification_qc.csv contain all values. No animals were excluded.')
    pca=pd.read_csv(OUT/'pca.csv');variance=pd.read_csv(OUT/'pca_variance.csv')
    fig,ax=plt.subplots(figsize=(8,6.5),layout='constrained')
    texts=[]
    for s in pca.to_dict('records'):
        color=COLORS[s['age']+'_'+s['genotype']];marker='o' if s['exposure']=='FA' else '^'
        ax.scatter(s['PC1'],s['PC2'],c=color,marker=marker,s=65,edgecolor='black',linewidth=.6)
        texts.append(ax.text(s['PC1'],s['PC2'],s['sample_id'],fontsize=7))
    ax.set(xlabel=f"PC1 ({100*variance.variance_explained.iloc[0]:.1f}%)",
      ylabel=f"PC2 ({100*variance.variance_explained.iloc[1]:.1f}%)",title='PCA · 500 highest-variance genes after VST')
    handles=[Line2D([0],[0],marker='s',linestyle='',color=c,label=k.replace('_',' · ')) for k,c in COLORS.items()]
    handles += [Line2D([0],[0],marker=m,linestyle='',color='black',label=l)for m,l in [('o','Filtered air'),('^','O₃')]]
    ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.11),ncol=3,frameon=False,fontsize=9)
    np.random.seed(20261018)
    adjust_text(texts,x=pca.PC1,y=pca.PC2,ax=ax,iter_lim=500,
      expand=(1.25,1.35),force_text=(.2,.4),force_static=(.15,.3),
      arrowprops={'arrowstyle':'-','color':'#666666','lw':.5})
    save(fig,'pca','PCA of the 500 highest-variance genes after a variance-stabilizing transformation estimated without using labels. Each point is one animal; color denotes starting age/genotype and shape exposure. PCA is descriptive and does not determine sample exclusions.')
    keys=['M3_E3_O3_vs_FA','M3_E4_O3_vs_FA','M17_E3_O3_vs_FA','M17_E4_O3_vs_FA']
    tables={key:pd.read_csv(OUT/f'de_{key}.csv') for key in keys}
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for ax,(key,d) in zip(axs.flat,tables.items()):
        valid=d.padj.notna()&d.log2FoldChange.notna();p=d.loc[valid]
        sig=p.padj<.05;cl=np.where(sig,np.where(p.log2FoldChange>0,'#D55E00','#0072B2'),'#888888')
        ax.scatter(p.log2FoldChange,-np.log10(np.clip(p.padj,1e-300,1)),c=cl,s=8,alpha=.7,rasterized=True)
        ax.axhline(-np.log10(.05),color='black',linestyle='--',linewidth=.8)
        ax.set(xlabel='Log₂ fold change O₃ / filtered air (MLE)',ylabel='−log₁₀ FDR BH',
         title=key.replace('_O3_vs_FA','').replace('_',' · ')+f' · {int(sig.sum())} genes FDR < 0.05')
    save(fig,'volcano','Wald tests of O₃ versus filtered air, three animals per group within each stratum. The full factorial model provides MLE log fold changes and FDR. Orange/up and blue/down indicate FDR<0.05; other genes are gray. Missing-FDR genes are omitted from plots but retained in tables. These plots do not establish spatial localization of the exposure effect.')
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for ax,(key,d) in zip(axs.flat,tables.items()):
        p=d[d.baseMean.gt(0)&d.log2FoldChange_shrunk.notna()];sig=p.padj.lt(.05)
        ax.scatter(p.baseMean,p.log2FoldChange_shrunk,c=np.where(sig,'#0072B2','#999999'),s=7,alpha=.65,rasterized=True)
        ax.set(xscale='log',xlabel='Mean normalized fragments (24 libraries; log scale)',
         ylabel='Shrunken log₂ fold change (apeglm)',title=key.replace('_O3_vs_FA','').replace('_',' · '))
        ax.axhline(0,color='black',linewidth=.7)
    save(fig,'ma','MA plots use apeglm-shrunken effects for visualization. Blue indicates Wald FDR<0.05. Shrinkage stabilizes weakly informed estimates and does not replace contrast evidence or the original p-values.')
    gl=pd.read_csv(OUT/'de_global_ozone_LRT.csv');sel=gl[gl.padj.lt(.05)].head(24)
    exploratory=sel.empty
    if exploratory:sel=gl[gl.pvalue.notna()].head(24)
    vst=pd.read_csv(OUT/'vst_counts.csv.gz',index_col=0)
    matrix=vst.loc[sel.gene_id].to_numpy();matrix=(matrix-matrix.mean(1,keepdims=True))/np.maximum(matrix.std(1,keepdims=True),1e-8)
    ordered=samples.assign(age_order=samples.age.map({'M3':0,'M17':1})).sort_values(
      ['age_order','genotype','exposure','sample_id'])
    order=ordered.sample_id
    ix=[list(vst.columns).index(s) for s in order]
    fig,ax=plt.subplots(figsize=(12,9),layout='constrained')
    lim=max(2,float(np.nanmax(np.abs(matrix)))) if matrix.size else 2
    sns.heatmap(matrix[:,ix],ax=ax,cmap='RdBu_r',vmin=-lim,vmax=lim,center=0,
      xticklabels=order,yticklabels=sel.gene_symbol,cbar_kws={'label':'VST standardized per gene (z)'} )
    ax.set(title=('Exploratory · smallest global LRT p-values' if exploratory else 'Genes with a global O₃ signal · LRT FDR < 0.05'),
      xlabel='Animal (ordered by starting age, genotype, exposure)',ylabel='Gene')
    group_axis=ax.secondary_xaxis('top')
    group_axis.set_xticks(np.arange(1.5,24,3))
    group_axis.set_xticklabels([f"{s.age[1:]}M {s.genotype}\n{s.exposure}"for s in ordered.iloc[::3].itertuples()],fontsize=8)
    group_axis.tick_params(length=0)
    for boundary in range(3,24,3):ax.axvline(boundary,color='black',linewidth=.5)
    ax.tick_params(axis='x',labelsize=8,rotation=60);ax.tick_params(axis='y',labelsize=8)
    save(fig,'top_genes_heatmap','Heatmap of up to 24 genes with the lowest FDR from the global four-degree-of-freedom LRT, standardized per gene after VST. If no genes meet FDR<0.05, it is explicitly exploratory. Columns are 24 animals, not cells. gene_filter.csv and de_global_ozone_LRT.csv document selection and complete results.')
    summary=json.loads((OUT/'summary.json').read_text())
    summary['raw_qc']={'sampled_pairs_per_library':200000,'q30_range':[float(raw.q30_rate.min()),float(raw.q30_rate.max())],
      'adapter_fraction_range':[float(raw.adapter_trimmed_read_fraction.min()),float(raw.adapter_trimmed_read_fraction.max())]}
    summary['quantification']={'libraries':len(quant),'processed_pairs_total':int(quant.num_processed_fragments.sum()),
      'percent_mapped_range':[float(quant.percent_mapped.min()),float(quant.percent_mapped.max())]}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False))
    manifest={'purpose':'general scientific exploratory figures, no target-journal requirements claimed',
      'software':{'python':platform.python_version(),'matplotlib':matplotlib.__version__,'pandas':pd.__version__,'seaborn':sns.__version__,'adjustText':adjustText.__version__},
      'captions':captions,'formats':['png','svg','pdf'],'png_dpi':300,
      'sources':['results/ozone/raw_qc_summary.csv','results/ozone/quantification_qc.csv','results/ozone/pca.csv',
       'results/ozone/de_all_four_strata.csv','results/ozone/de_global_ozone_LRT.csv','results/ozone/vst_counts.csv.gz']}
    (FIG/'figure_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
    (FIG/'captions.txt').write_text('\n\n'.join(k+'\n'+v for k,v in captions.items()))
    subprocess.run([sys.executable,ROOT/'scripts/ozone_validate.py'],check=True)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
