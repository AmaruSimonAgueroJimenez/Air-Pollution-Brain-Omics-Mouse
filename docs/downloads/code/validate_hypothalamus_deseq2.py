#!/usr/bin/env python3
"""Run one fixed-model R DESeq2 fit and compare it with primary PyDESeq2 MLEs."""
from __future__ import annotations
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "results/hypothalamus"
OUT = PRIMARY / "independent_deseq2_check"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [clean(v) for v in value]
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def correlation(a, b, method="pearson"):
    return float(pd.Series(a).corr(pd.Series(b), method=method))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-R-fit", action="store_true", help="Compare the saved fixed-model R results without refitting")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = [PRIMARY/name for name in ("counts.tsv.gz", "sample_metadata.tsv", "differential_expression_mle.tsv", "qc_samples.tsv", "methods.json")]
    before = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    if not args.reuse_R_fit:
        with (OUT / "R_fit.log").open("w") as log:
            subprocess.run(["Rscript", "--vanilla", str(ROOT / "scripts/validate_hypothalamus_deseq2.R")], cwd=ROOT,
                           env=os.environ.copy(), stdout=log, stderr=subprocess.STDOUT, check=True)
    manifest = json.loads((OUT / "R_fit_manifest.json").read_text())
    primary_methods = json.loads((PRIMARY / "methods.json").read_text())
    py = pd.read_csv(PRIMARY/"differential_expression_mle.tsv", sep="\t", index_col="gene")
    r = pd.read_csv(OUT/"R_DESeq2_no_independent_filter.tsv", sep="\t", index_col="gene")
    r_default = pd.read_csv(OUT/"R_DESeq2_default_independent_filter.tsv", sep="\t", index_col="gene")
    r_uncensored = pd.read_csv(OUT/"R_DESeq2_no_Cooks_filter_diagnostic.tsv", sep="\t", index_col="gene")
    assert set(py.index) == set(r.index) and len(py)==19611
    r, r_default, r_uncensored = [x.loc[py.index] for x in (r, r_default, r_uncensored)]
    common = py.index[np.isfinite(py.pvalue) & np.isfinite(r.pvalue)]
    comparison = py.add_suffix("_PyDESeq2").join(r.add_suffix("_R_DESeq2"))
    comparison["padj_R_default_independent_filter"] = r_default.padj
    comparison["in_common_finite_pvalue_universe"] = comparison.index.isin(common)
    comparison["padj_common_PyDESeq2"] = np.nan
    comparison["padj_common_R_DESeq2"] = np.nan
    comparison.loc[common, "padj_common_PyDESeq2"] = multipletests(py.loc[common,"pvalue"], method="fdr_bh")[1]
    comparison.loc[common, "padj_common_R_DESeq2"] = multipletests(r.loc[common,"pvalue"], method="fdr_bh")[1]
    comparison["log2FC_difference_R_minus_Python"] = r.log2FoldChange-py.log2FoldChange
    comparison["SE_ratio_R_over_Python"] = r.lfcSE/py.lfcSE
    comparison["FDR05_Python_primary"] = py.padj.lt(.05)
    comparison["FDR05_R_no_independent_filter"] = r.padj.lt(.05)
    comparison["FDR05_R_default_independent_filter"] = r_default.padj.lt(.05)
    comparison["FDR05_common_Python"] = comparison.padj_common_PyDESeq2.lt(.05)
    comparison["FDR05_common_R"] = comparison.padj_common_R_DESeq2.lt(.05)
    comparison.to_csv(OUT/"all_gene_software_comparison.tsv", sep="\t", na_rep="NA", float_format="%.12g")
    discordant = comparison[comparison.FDR05_Python_primary.ne(comparison.FDR05_R_no_independent_filter)]
    discordant.to_csv(OUT/"discordant_FDR_genes.tsv", sep="\t", na_rep="NA", float_format="%.12g")
    comparison.loc[common].to_csv(OUT/"common_finite_pvalue_family_comparison.tsv", sep="\t", na_rep="NA", float_format="%.12g")
    py_qc = pd.read_csv(PRIMARY/"qc_samples.tsv", sep="\t", index_col="sample")
    r_sf = pd.read_csv(OUT/"R_size_factors.tsv", sep="\t", index_col="sample").loc[py_qc.index]
    sf = pd.DataFrame({"size_factor_PyDESeq2": py_qc.size_factor,"size_factor_R_DESeq2":r_sf.size_factor})
    sf["relative_difference_R_vs_Python"] = sf.size_factor_R_DESeq2/sf.size_factor_PyDESeq2-1
    sf.to_csv(OUT/"size_factor_comparison.tsv", sep="\t", float_format="%.12g")
    py_cooks = py.pvalue.isna() & np.isfinite(py.stat)
    r_cooks = r.pvalue.isna() & np.isfinite(r_uncensored.pvalue)
    cooks = pd.DataFrame({"pvalue_filtered_PyDESeq2":py_cooks,"Cooks_filtered_R_DESeq2":r_cooks})
    cooks.to_csv(OUT/"Cooks_filter_comparison.tsv", sep="\t")
    delta = np.abs(r.loc[common,"log2FoldChange"]-py.loc[common,"log2FoldChange"])
    se_ratio = r.loc[common,"lfcSE"]/py.loc[common,"lfcSE"]
    logp_py = -np.log10(py.loc[common,"pvalue"].clip(lower=1e-300))
    logp_r = -np.log10(r.loc[common,"pvalue"].clip(lower=1e-300))
    metrics = {"size_factor_max_absolute_relative_difference":float(sf.relative_difference_R_vs_Python.abs().max()),
               "common_universe_log2FC_Pearson":correlation(py.loc[common,"log2FoldChange"],r.loc[common,"log2FoldChange"]),
               "common_universe_log2FC_Spearman":correlation(py.loc[common,"log2FoldChange"],r.loc[common,"log2FoldChange"],"spearman"),
               "median_absolute_log2FC_difference":float(delta.median()), "p99_absolute_log2FC_difference":float(delta.quantile(.99)),
               "maximum_absolute_log2FC_difference":float(delta.max()), "median_SE_ratio_R_over_Python":float(se_ratio.median()),
               "p01_SE_ratio_R_over_Python":float(se_ratio.quantile(.01)),"p99_SE_ratio_R_over_Python":float(se_ratio.quantile(.99)),
               "common_universe_Wald_statistic_Pearson":correlation(py.loc[common,"stat"],r.loc[common,"stat"]),
               "common_universe_minus_log10_pvalue_Spearman":correlation(logp_py,logp_r,"spearman")}
    tolerances = {"size_factors_relative_difference_at_most":1e-6,"log2FC_Pearson_at_least":.995,
                  "median_absolute_log2FC_difference_at_most":.01,"p99_absolute_log2FC_difference_at_most":.10,
                  "median_SE_ratio_between":[.95,1.05],"Wald_statistic_Pearson_at_least":.99,
                  "minus_log10_pvalue_Spearman_at_least":.99}
    flags = {"size_factors":metrics["size_factor_max_absolute_relative_difference"]<=1e-6,
             "log2FC_global_correlation":metrics["common_universe_log2FC_Pearson"]>=.995,
             "log2FC_median_difference":metrics["median_absolute_log2FC_difference"]<=.01,
             "log2FC_99th_percentile_difference":metrics["p99_absolute_log2FC_difference"]<=.10,
             "median_SE_ratio":.95<=metrics["median_SE_ratio_R_over_Python"]<=1.05,
             "Wald_statistic_global_correlation":metrics["common_universe_Wald_statistic_Pearson"]>=.99,
             "pvalue_rank_correlation":metrics["common_universe_minus_log10_pvalue_Spearman"]>=.99}
    differences = []
    if not all(flags.values()): differences.append("At least one descriptive numerical-agreement tolerance was not met; inspect the gene-level table before interpreting agreement.")
    if len(discordant): differences.append(f"{len(discordant)} gene(s) change native BH FDR<0.05 status between primary PyDESeq2 and R DESeq2 without independent filtering.")
    if int((py_cooks != r_cooks).sum()): differences.append(f"{int((py_cooks != r_cooks).sum())} genes differ in software Cook's/p-value filtering status.")
    tll = comparison.loc["Tll1"].to_dict() if "Tll1" in comparison.index else {"status":"not in input"}
    if "Tll1" in comparison.index and tll["FDR05_Python_primary"] != tll["FDR05_R_no_independent_filter"]:
        differences.append("Tll1 changes FDR significance across implementations; its effect should be described as software-sensitive.")
    for p in inputs:
        assert sha(p)==before[str(p.relative_to(ROOT))], f"Primary input changed: {p}"
    summary = {"purpose":"Independent software cross-check of the fixed primary hypothesis, not biological replication or a substitute primary analysis",
               "primary_analysis_unchanged":True,"input_SHA256":before,"R_fit":manifest,
               "primary_PyDESeq2_version":primary_methods["versions"]["pydeseq2"],
               "validation_script_SHA256":{"R":sha(ROOT/"scripts/validate_hypothalamus_deseq2.R"),"Python":sha(Path(__file__))},
               "fitted_genes":len(py),"Python_finite_pvalues":int(py.pvalue.notna().sum()),"R_finite_pvalues":int(r.pvalue.notna().sum()),
               "common_finite_pvalue_universe":len(common),"same_gene_universe":set(py.index)==set(r.index),
               "native_FDR_hits":{"PyDESeq2_primary":int(py.padj.lt(.05).sum()),"R_independent_filter_FALSE":int(r.padj.lt(.05).sum()),"R_default_independent_filter_TRUE":int(r_default.padj.lt(.05).sum())},
               "shared_family_BH_FDR_hits":{"PyDESeq2":int(comparison.padj_common_PyDESeq2.lt(.05).sum()),"R_DESeq2":int(comparison.padj_common_R_DESeq2.lt(.05).sum())},
               "Cooks_filter":{"Python_genes":int(py_cooks.sum()),"R_genes":int(r_cooks.sum()),"same_status_genes":int(py_cooks.eq(r_cooks).sum()),"discordant_genes":int(py_cooks.ne(r_cooks).sum())},
               "numerical_agreement_metrics":metrics,"descriptive_tolerances":tolerances,"tolerance_checks":flags,
               "all_descriptive_tolerances_met":all(flags.values()),"FDR_discordant_genes":discordant.index.tolist(),
               "FDR_discordant_gene_details":discordant.reset_index().to_dict("records"),
               "Tll1":tll,"material_differences":differences,
               "interpretation":"Numerical and threshold agreement assess implementation sensitivity only. FDR crossings do not identify the preferred software, and agreement does not compensate for seven biological samples or control-group heterogeneity. No sample, gene threshold or model was adjusted to increase significant hits.",
               "completed_at":datetime.now(timezone.utc).isoformat()}
    (OUT/"summary.json").write_text(json.dumps(clean(summary),indent=2,allow_nan=False)+'\n')
    print(json.dumps(clean({"native_FDR_hits":summary["native_FDR_hits"],"common_genes":len(common),"metrics":metrics,"Tll1":tll,"material_differences":differences}),indent=2),flush=True)


if __name__=="__main__":
    main()
