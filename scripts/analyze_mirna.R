#!/usr/bin/env Rscript
# Limma on published quantile-normalized signals; no assumed log2 scale.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Supply the project root as the only argument")
root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(file.path(root, ".r-library"), .libPaths()))
suppressPackageStartupMessages(library(limma))
out <- file.path(root, "results", "mirna")
x <- as.matrix(read.delim(file.path(out, "analysis_signal.tsv"), row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(out, "sample_metadata.tsv"), row.names = 1, check.names = FALSE)
if (!identical(colnames(x), rownames(meta))) stop("Expression and metadata sample order differs")
if (any(!is.finite(x))) stop("Nonfinite signals cannot enter limma")
meta$condition <- factor(meta$condition, levels = c("Control", "PM1", "PM5"))
design <- model.matrix(~0 + condition, data = meta)
colnames(design) <- levels(meta$condition)
contrast <- makeContrasts(PM1_vs_Control = PM1 - Control,
                          PM5_vs_Control = PM5 - Control, levels = design)
fit <- lmFit(x, design)
fit_c <- contrasts.fit(fit, contrast)
fit_c <- eBayes(fit_c, trend = TRUE, robust = FALSE)
for (name in colnames(contrast)) {
    tab <- topTable(fit_c, coef = name, number = Inf, adjust.method = "BH", sort.by = "P")
    # limma calls the coefficient logFC generically. Its units here are the
    # published normalized signal units because log2 preprocessing is unverified.
    names(tab)[names(tab) == "logFC"] <- "effect_mean_difference"
    names(tab)[names(tab) == "AveExpr"] <- "published_signal_mean"
    names(tab)[names(tab) == "P.Value"] <- "pvalue"
    names(tab)[names(tab) == "adj.P.Val"] <- "padj"
    tab$probe_id <- rownames(tab)
    position <- match(rownames(tab), rownames(fit_c$coefficients))
    tab$SE_moderated <- fit_c$stdev.unscaled[position, name] * sqrt(fit_c$s2.post[position])
    tab$df_moderated <- fit_c$df.total[position]
    tab$CI95_lower <- tab$effect_mean_difference - qt(.975, tab$df_moderated) * tab$SE_moderated
    tab$CI95_upper <- tab$effect_mean_difference + qt(.975, tab$df_moderated) * tab$SE_moderated
    write.table(tab, file.path(out, paste0("limma_", name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}
omnibus <- topTable(fit_c, coef = c(1, 2), number = Inf, adjust.method = "BH", sort.by = "F")
names(omnibus)[names(omnibus) == "PM1_vs_Control"] <- "effect_PM1_vs_Control"
names(omnibus)[names(omnibus) == "PM5_vs_Control"] <- "effect_PM5_vs_Control"
names(omnibus)[names(omnibus) == "AveExpr"] <- "published_signal_mean"
names(omnibus)[names(omnibus) == "P.Value"] <- "pvalue"
names(omnibus)[names(omnibus) == "adj.P.Val"] <- "padj"
omnibus$probe_id <- rownames(omnibus)
write.table(omnibus, file.path(out, "limma_omnibus.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(design, file.path(out, "limma_design.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(contrast, file.path(out, "limma_contrasts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
writeLines(c(paste0("R: ", R.version.string), paste0("limma: ", packageVersion("limma")),
             "Design: ~0+condition; Control, PM1, PM5 (three samples each)",
             "Empirical Bayes: trend=TRUE, robust=FALSE", "No additional normalization or log transform",
             capture.output(sessionInfo())), file.path(out, "R_session_info.txt"))
message("Limma completed for ", nrow(x), " miRNA probe-level rows and ", ncol(x), " samples")
