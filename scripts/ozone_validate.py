#!/usr/bin/env python3
"""Independent, source-output numerical validation of the completed O3 analysis."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control,norm,chi2
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ozone'

def bh(p):
    p=np.asarray(p,dtype=float);valid=np.isfinite(p);q=np.full(len(p),np.nan)
    q[valid]=false_discovery_control(p[valid],method='bh');return q

def compare(actual,expected):
    a=np.asarray(actual,dtype=float);e=np.asarray(expected,dtype=float)
    assert np.array_equal(np.isnan(a),np.isnan(e))
    assert np.allclose(a,e,atol=1e-11,rtol=1e-7,equal_nan=True)

def logical(s):return s.astype(str).str.lower().eq('true').to_numpy()

def main():
    summary=json.loads((OUT/'summary.json').read_text());samples=pd.read_csv(OUT/'samples.csv')
    assert summary['design_rank']==8
    assert len(samples)==24 and samples.sample_id.is_unique
    assert (samples.groupby(['age','genotype','exposure']).size()==3).all()
    counts=pd.read_csv(OUT/'analysis_gene_counts.csv.gz',index_col=0)
    assert counts.index.is_unique and set(counts.columns)==set(samples.sample_id)
    assert len(counts)==summary['n_retained_genes']
    values=counts.to_numpy();assert np.isfinite(values).all() and (values>=0).all()
    assert np.array_equal(values,np.rint(values)) and ((values>=10).sum(1)>=3).all()
    all_tables=pd.read_csv(OUT/'de_all_four_strata.csv')
    assert set(all_tables.contrast)=={f'{a}_{g}_O3_vs_FA'for a in ['M3','M17']for g in ['E3','E4']}
    assert len(all_tables)==4*len(counts)
    assert not all_tables[['contrast','gene_id']].duplicated().any()
    compare(all_tables.padj_across_four_strata,bh(all_tables.pvalue))
    diagnostics={};files=[]
    for age in ['M3','M17']:
        for genotype in ['E3','E4']:
            key=f'{age}_{genotype}_O3_vs_FA';path=OUT/f'de_{key}.csv';d=pd.read_csv(path)
            assert d.gene_id.is_unique and set(d.gene_id)==set(counts.index)
            pooled=all_tables.loc[all_tables.contrast.eq(key)].set_index('gene_id')
            assert len(pooled)==len(d) and set(pooled.index)==set(d.gene_id)
            linked=pooled.loc[d.gene_id]
            for field in ['baseMean','log2FoldChange','stat','pvalue','padj','padj_across_four_strata','log2FoldChange_shrunk']:
                compare(d[field],linked[field])
            compare(d.padj,bh(d.pvalue))
            good=np.isfinite(d.pvalue);compare(d.loc[good,'pvalue'],2*norm.sf(np.abs(d.loc[good,'stat'])))
            selected=samples.loc[samples.age.eq(age)&samples.genotype.eq(genotype),'sample_id']
            zero=counts.loc[d.gene_id,selected].sum(1).eq(0).to_numpy()
            assert np.array_equal(zero,logical(d.stratum_all_zero))
            assert d.loc[zero,'pvalue'].isna().all() and d.loc[zero,'log2FoldChange_shrunk'].isna().all()
            numerical_bad=~logical(d.beta_converged)
            assert d.loc[numerical_bad,'pvalue'].isna().all()
            eligible=logical(d.effect_eligible_for_projection)
            assert not (eligible&(zero|numerical_bad|~logical(d.apeglm_converged))).any()
            assert np.isfinite(d.loc[eligible,['pvalue','padj']].to_numpy()).all()
            assert np.isfinite(d.loc[eligible,['log2FoldChange','lfcSE','log2FoldChange_shrunk','lfcSE_shrunk']].to_numpy()).all()
            n_sig=int(d.padj.lt(.05).sum());assert n_sig==summary['contrasts'][key]['fdr_005']
            diagnostics[key]={'tested':int(good.sum()),'FDR_005':n_sig,'all_zero_stratum':int(zero.sum()),
              'wald_nonconverged':int(numerical_bad.sum()),'apeglm_nonconverged':int((~logical(d.apeglm_converged)).sum())}
            files.append(path)
    path=OUT/'de_global_ozone_LRT.csv';g=pd.read_csv(path)
    assert g.gene_id.is_unique and set(g.gene_id)==set(counts.index)
    assert g.degrees_of_freedom.eq(4).all();compare(g.padj,bh(g.pvalue))
    good=np.isfinite(g.pvalue);compare(g.loc[good,'pvalue'],chi2.sf(g.loc[good,'stat'],4))
    bad=~logical(g.full_beta_converged)|~logical(g.reduced_beta_converged)
    assert g.loc[bad,'pvalue'].isna().all()
    assert int(g.padj.lt(.05).sum())==summary['global_LRT']['fdr_005']
    quant=pd.read_csv(OUT/'quantification_qc.csv')
    assert len(quant)==24 and logical(quant.quantification_integrity_passed).all()
    assert logical(quant.em_converged).all() and np.allclose(quant.sum_tpm,1e6,rtol=1e-5)
    assert np.allclose(quant.sum_estimated_transcript_fragments,quant.num_mapped_fragments,rtol=1e-6,atol=1)
    full_qc=pd.read_csv(OUT/'full_qc_summary.csv')
    assert len(full_qc)==24 and full_qc.sample_id.is_unique
    linked=full_qc.merge(quant[['sample_id','num_processed_fragments']],on='sample_id',validate='one_to_one')
    assert len(linked)==24 and linked.after_adapter_pairs.eq(linked.num_processed_fragments).all()
    assert linked.raw_pairs.ge(linked.after_adapter_pairs).all()
    for field in ['raw_q20_rate','raw_q30_rate','raw_gc_fraction','after_adapter_q20_rate','after_adapter_q30_rate','after_adapter_gc_fraction']:
        assert linked[field].between(0,1).all()
    files += [path,OUT/'de_all_four_strata.csv',OUT/'analysis_gene_counts.csv.gz',OUT/'samples.csv']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files}
    report={'passed':True,'study':'E-MTAB-14917','samples':24,'model_rank':8,'retained_genes':len(counts),
      'checks':['24 unique libraries, 3 per stratum','full model rank is eight','integer count input and fixed filter',
       'pooled and separate contrast tables match; global LRT shares the count universe',
       'per-contrast BH independently reproduced with SciPy','pooled-four BH independently reproduced',
       'Wald normal-tail p-values','global four-df chi-square LRT p-values',
       'all-zero stratum and numerical-convergence exclusions','posterior effects valid for projection',
       'all-read QC covers 24 libraries and retained pair totals match quantification input',
       'all 24 Salmon EM optimizations converged and transcript TPM sums verified'],
      'contrasts':diagnostics,'global_LRT_nonconverged':int(bad.sum()),'source_output_sha256':hashes}
    (OUT/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
