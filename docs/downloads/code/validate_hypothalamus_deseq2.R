#!/usr/bin/env Rscript
# Independent software cross-check; does not replace the primary PyDESeq2 fit.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
script <- normalizePath(gsub("~+~", " ", script_arg, fixed=TRUE))
root <- dirname(dirname(script))
.libPaths(c(file.path(root, ".r-library"), .libPaths()))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(jsonlite))
register(SerialParam())
options(digits=17)
output <- file.path(root, "results/hypothalamus/independent_deseq2_check")
dir.create(output, recursive=TRUE, showWarnings=FALSE)
counts <- as.matrix(read.table(gzfile(file.path(root, "results/hypothalamus/counts.tsv.gz")),
                              header=TRUE, row.names=1, sep="\t", check.names=FALSE))
meta <- read.table(file.path(root, "results/hypothalamus/sample_metadata.tsv"),
                   header=TRUE, row.names=1, sep="\t", check.names=FALSE, quote="", comment.char="")
stopifnot(identical(colnames(counts), rownames(meta)), nrow(counts)==50417L, ncol(counts)==7L,
          !anyNA(counts), all(counts>=0), all(counts==floor(counts)), max(counts)<.Machine$integer.max)
meta$condition <- factor(meta$condition, levels=c("FA", "PM"))
stopifnot(sum(meta$condition=="FA")==4L, sum(meta$condition=="PM")==3L)
keep <- rowSums(counts>=10)>=3
counts <- counts[keep,,drop=FALSE]
stopifnot(nrow(counts)==19611L)
storage.mode(counts) <- "integer"
dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta["condition"], design=~condition)
message("Fitting the fixed PM-versus-FA model in DESeq2 R with one CPU.")
dds <- DESeq(dds, test="Wald", fitType="parametric", sfType="ratio", betaPrior=FALSE,
             minReplicatesForReplace=7, parallel=FALSE, quiet=TRUE)
contrast <- c("condition", "PM", "FA")
res_common <- results(dds, contrast=contrast, alpha=.05, independentFiltering=FALSE)
res_default <- results(dds, contrast=contrast, alpha=.05, independentFiltering=TRUE)
res_no_cooks <- results(dds, contrast=contrast, alpha=.05, independentFiltering=FALSE, cooksCutoff=FALSE)
save_result <- function(result, filename) {
  frame <- data.frame(gene=rownames(result), as.data.frame(result), row.names=NULL, check.names=FALSE)
  write.table(frame, file.path(output, filename), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}
save_result(res_common, "R_DESeq2_no_independent_filter.tsv")
save_result(res_default, "R_DESeq2_default_independent_filter.tsv")
save_result(res_no_cooks, "R_DESeq2_no_Cooks_filter_diagnostic.tsv")
write.table(data.frame(sample=colnames(dds), condition=dds$condition, size_factor=sizeFactors(dds)),
            file.path(output, "R_size_factors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cooks <- assays(dds)[["cooks"]]
num_coefs <- ncol(model.matrix(design(dds), colData(dds)))
cutoff <- qf(.99, num_coefs, ncol(dds)-num_coefs)
gene_diagnostics <- data.frame(gene=rownames(dds), dispersion=dispersions(dds),
   dispersion_gene_estimate=mcols(dds)$dispGeneEst, dispersion_fitted=mcols(dds)$dispFit,
   beta_converged=mcols(dds)$betaConv, max_Cooks_distance=apply(cooks, 1, max),
   Cooks_filtered_pvalue=is.na(res_common$pvalue) & is.finite(res_no_cooks$pvalue))
write.table(gene_diagnostics, file.path(output, "R_gene_fit_diagnostics.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write.table(data.frame(sample=colnames(dds), raw_Cooks_exceedances=colSums(cooks>cutoff),
                       max_Cooks_distance=apply(cooks, 2, max), raw_Cooks_cutoff=cutoff),
            file.path(output, "R_sample_Cooks_diagnostics.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
manifest <- list(
  purpose="Independent DESeq2 implementation cross-check, not biological replication or a replacement primary analysis",
  software=list(R=R.version.string, DESeq2=as.character(packageVersion("DESeq2"))),
  condition_reference="FA", contrast=contrast, design="~condition", samples=7,
  groups=list(FA=4, PM=3), prefilter="count >=10 in at least 3 of 7 samples before size-factor estimation",
  input_genes=50417, fitted_genes=nrow(dds), size_factor_method="ratio", fit_type="parametric",
  test="two-sided Wald", beta_prior=FALSE, minimum_replicates_for_count_replacement=7,
  replacement_assay_present="replaceCounts" %in% assayNames(dds), parallel=FALSE, CPU_limit=1,
  default_Cooks_filter_enabled=TRUE, Cooks_cutoff=cutoff,
  finite_pvalues_no_independent_filter=sum(is.finite(res_common$pvalue)),
  Cooks_filtered_genes=sum(gene_diagnostics$Cooks_filtered_pvalue),
  FDR_hits_no_independent_filter=sum(res_common$padj<.05, na.rm=TRUE),
  FDR_hits_default_independent_filter=sum(res_default$padj<.05, na.rm=TRUE),
  independent_filter_threshold=as.numeric(metadata(res_default)$filterThreshold),
  fitted_dispersion_trend_type=attr(dispersionFunction(dds), "fitType"),
  completed_at=format(Sys.time(), tz="UTC", format="%Y-%m-%dT%H:%M:%SZ"))
write_json(manifest, file.path(output, "R_fit_manifest.json"), auto_unbox=TRUE, pretty=TRUE, na="null")
capture.output(sessionInfo(), file=file.path(output, "R_session_info.txt"))
message("Independent fixed-model DESeq2 fit completed.")
