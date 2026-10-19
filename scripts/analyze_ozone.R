#!/usr/bin/env Rscript
# Gene-level inference from all 24 original E-MTAB-14917 sequencing libraries.
# The sample, not each read or atlas cell, is the unit of inference.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grep('^--file=', args)])
# Rscript on this macOS host encodes spaces in --file= as the literal ~+~ token.
script <- gsub('~+~', ' ', script, fixed=TRUE)
root <- normalizePath(file.path(dirname(script), '..'),mustWork=TRUE)
setwd(root)
.libPaths(c(file.path(root, '.r-library'), .libPaths()))
suppressPackageStartupMessages({library(tximport);library(DESeq2);library(apeglm);library(jsonlite)})
set.seed(20261018)
out <- file.path(root, 'results/ozone');dir.create(out, recursive=TRUE, showWarnings=FALSE)
meta <- read.csv('metadata/ozone_analysis/samples.csv', stringsAsFactors=FALSE)
rownames(meta) <- meta$sample_id
meta$age <- factor(meta$age, levels=c('M3','M17'))
meta$genotype <- factor(meta$genotype, levels=c('E3','E4'))
meta$exposure <- factor(meta$exposure, levels=c('FA','O3'))
meta$group <- factor(paste(meta$age,meta$genotype,meta$exposure,sep='_'))
stopifnot(nrow(meta)==24, all(table(meta$age,meta$genotype,meta$exposure)==3))
files <- setNames(file.path(out,'quant',meta$sample_id,'quant.sf'),meta$sample_id)
stopifnot(all(file.exists(files)))
for(sid in meta$sample_id){
 info<-fromJSON(file.path(out,'quant',sid,'aux_info/meta_info.json'))
 stopifnot(isTRUE(info$meta_info_complete),isTRUE(info$em_converged),
           length(info$quant_errors)==0,info$num_mapped>0)
}
txmap <- read.delim('data/references/ozone/tx2gene.tsv', stringsAsFactors=FALSE)
genes <- read.delim('data/references/ozone/genes.tsv', stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(txmap$transcript_id),!anyDuplicated(genes$gene_id))
message('Importing Salmon quantified fragments and effective transcript lengths')
txi <- tximport(files, type='salmon', tx2gene=txmap, countsFromAbundance='no')
saveRDS(txi,file.path(out,'tximport.rds'))
write.csv(meta,file.path(out,'samples.csv'),row.names=FALSE)
write_matrix <- function(x,path){
 con<-if(grepl('\\.gz$',path))gzfile(path,'wt')else file(path,'wt')
 on.exit(close(con));write.csv(data.frame(gene_id=rownames(x),x,check.names=FALSE),con,row.names=FALSE)
}
write_matrix(txi$counts,file.path(out,'estimated_gene_fragments.csv.gz'))
write_matrix(txi$abundance,file.path(out,'gene_tpm.csv.gz'))
dds <- DESeqDataSetFromTximport(txi, colData=meta, design=~age*genotype*exposure)
# A conservative, treatment-independent filter fixed before examining DE results.
keep <- rowSums(counts(dds)>=10)>=3
filter_table <- data.frame(gene_id=rownames(dds), retained=keep,
 total_estimated_fragments=rowSums(counts(dds)), n_samples_ge10=rowSums(counts(dds)>=10))
write.csv(filter_table,file.path(out,'gene_filter.csv'),row.names=FALSE)
n_input <- nrow(dds); dds <- dds[keep,]
write_matrix(counts(dds),file.path(out,'analysis_gene_counts.csv.gz'))
message('Fitting full age * genotype * exposure model: ',nrow(dds),' genes')
dds <- DESeq(dds, test='Wald', betaPrior=FALSE, quiet=FALSE)
saveRDS(dds,file.path(out,'dds_full_wald.rds'))
write_matrix(counts(dds,normalized=TRUE),file.path(out,'normalized_counts.csv.gz'))
vsd <- varianceStabilizingTransformation(dds,blind=TRUE)
write_matrix(assay(vsd),file.path(out,'vst_counts.csv.gz'))
rv <- matrixStats::rowVars(assay(vsd)); selected <- head(order(rv,decreasing=TRUE),min(500,length(rv)))
pc <- prcomp(t(assay(vsd)[selected,,drop=FALSE]),center=TRUE,scale.=FALSE)
pca <- data.frame(sample_id=rownames(pc$x),PC1=pc$x[,1],PC2=pc$x[,2],PC3=pc$x[,3])
write.csv(merge(pca,meta,by='sample_id',sort=FALSE),file.path(out,'pca.csv'),row.names=FALSE)
write.csv(data.frame(component=paste0('PC',seq_along(pc$sdev)),
 variance_explained=pc$sdev^2/sum(pc$sdev^2)),file.path(out,'pca_variance.csv'),row.names=FALSE)
write.csv(as.matrix(dist(t(assay(vsd)))),file.path(out,'sample_distances.csv'))
cooks <- assays(dds)[['cooks']]
qc <- data.frame(sample_id=colnames(dds),total_estimated_fragments=colSums(counts(dds)),
 detected_genes=colSums(counts(dds)>0),median_normalization_factor=apply(normalizationFactors(dds),2,median),
 max_cooks=apply(cooks,2,max,na.rm=TRUE))
write.csv(qc,file.path(out,'expression_qc.csv'),row.names=FALSE)
annotate <- function(tab){
 tab$gene_id <- rownames(tab)
 idx <- match(tab$gene_id,genes$gene_id)
 tab$gene_symbol <- genes$gene_symbol[idx]
 tab$gene_type <- genes$gene_type[idx]
 tab[,c('gene_id','gene_symbol','gene_type',setdiff(names(tab),c('gene_id','gene_symbol','gene_type')))]
}
tables <- list();summaries <- list()
for(age in levels(meta$age))for(genotype in levels(meta$genotype)){
 key <- paste(age,genotype,'O3_vs_FA',sep='_')
 fa <- paste(age,genotype,'FA',sep='_');ozone <- paste(age,genotype,'O3',sep='_')
 # This reparameterization has the same model space and fixed dispersions as the
 # full factorial model, and makes each desired LFC a single apeglm coefficient.
 d <- dds;d$group <- relevel(meta$group,ref=fa);design(d)<-~group
 d <- nbinomWaldTest(d,betaPrior=FALSE,maxit=1000,quiet=FALSE)
 coef <- paste0('group_',ozone,'_vs_',fa)
 stopifnot(coef %in% resultsNames(d))
 # The fixed count prefilter defines the universe; exclude numerical failures
 # before BH rather than letting an unreliable p-value affect other genes' FDR.
 res <- results(d,name=coef,alpha=0.05,independentFiltering=FALSE,cooksCutoff=TRUE)
 original <- as.data.frame(res)
 original$beta_converged <- mcols(d)$betaConv
 original$dispersion_outlier <- mcols(d)$dispOutlier
 stratum_samples <- rownames(meta)[meta$age==age&meta$genotype==genotype]
 original$stratum_total_count <- rowSums(counts(d)[,stratum_samples,drop=FALSE])
 original$stratum_all_zero <- original$stratum_total_count==0
 shrink_bundle <- lfcShrink(d,coef=coef,res=res,type='apeglm',returnList=TRUE,quiet=FALSE)
 shrink <- shrink_bundle$res
 stopifnot(identical(rownames(original),rownames(shrink)))
 original$log2FoldChange_shrunk <- shrink$log2FoldChange
 original$lfcSE_shrunk <- shrink$lfcSE
 original$apeglm_convergence_code <- shrink_bundle$fit$diag[rownames(original),'conv']
 original$apeglm_converged <- original$apeglm_convergence_code==0
 invalid_shrink <- is.na(original$apeglm_converged)|!original$apeglm_converged|
                   !is.finite(original$log2FoldChange_shrunk)|!is.finite(original$lfcSE_shrunk)
 original$log2FoldChange_shrunk[invalid_shrink] <- NA_real_
 original$lfcSE_shrunk[invalid_shrink] <- NA_real_
 original$wald_CI95_low <- original$log2FoldChange - qnorm(.975)*original$lfcSE
 original$wald_CI95_high <- original$log2FoldChange + qnorm(.975)*original$lfcSE
 original$contrast <- key
 # Numerical failures do not constitute interpretable evidence of a gene effect.
 invalid <- is.na(original$beta_converged)|!original$beta_converged|original$stratum_all_zero
 original$wald_pvalue_before_convergence_filter <- original$pvalue
 original$wald_log2FoldChange_before_validity_filter <- original$log2FoldChange
 original$wald_lfcSE_before_validity_filter <- original$lfcSE
 original$pvalue[invalid] <- NA_real_
 original$padj <- p.adjust(original$pvalue,method='BH')
 for(field in c('log2FoldChange','lfcSE','log2FoldChange_shrunk','lfcSE_shrunk','wald_CI95_low','wald_CI95_high')){
  original[[field]][invalid] <- NA_real_
 }
 original$effect_eligible_for_projection <- !invalid&!invalid_shrink&
   is.finite(original$pvalue)&is.finite(original$padj)&
   is.finite(original$log2FoldChange)&is.finite(original$lfcSE)&
   is.finite(original$log2FoldChange_shrunk)&is.finite(original$lfcSE_shrunk)
 tables[[key]] <- annotate(original)
}
# Explicit supplementary multiplicity family spanning the four prespecified contrasts.
all_tables <- do.call(rbind,tables)
all_tables$padj_across_four_strata <- p.adjust(all_tables$pvalue,method='BH')
for(key in names(tables)){
 tab <- all_tables[all_tables$contrast==key,]
 tab <- tab[order(tab$padj,tab$pvalue,na.last=TRUE),]
 write.csv(tab,file.path(out,paste0('de_',key,'.csv')),row.names=FALSE,na='')
 summaries[[key]] <- list(genes=nrow(tab),tested=sum(!is.na(tab$pvalue)),
 nonconverged=sum(!tab$beta_converged,na.rm=TRUE),
 shrinkage_nonconverged=sum(!tab$apeglm_converged,na.rm=TRUE),
 all_zero_in_stratum=sum(tab$stratum_all_zero),
 fdr_005=sum(tab$padj<.05,na.rm=TRUE),up_fdr_005=sum(tab$padj<.05&tab$log2FoldChange>0,na.rm=TRUE),
 down_fdr_005=sum(tab$padj<.05&tab$log2FoldChange<0,na.rm=TRUE),
 across_four_fdr_005=sum(tab$padj_across_four_strata<.05,na.rm=TRUE),
 n_FA=3,n_O3=3,top_genes=head(tab$gene_symbol[!is.na(tab$padj)],10))
}
write.csv(all_tables,file.path(out,'de_all_four_strata.csv'),row.names=FALSE,na='')
message('Global ozone likelihood-ratio test, full vs age * genotype (4 df)')
lrt <- nbinomLRT(dds,reduced=~age*genotype,maxit=1000,quiet=FALSE)
global <- as.data.frame(results(lrt,alpha=.05,independentFiltering=FALSE,cooksCutoff=TRUE))
global <- global[,c('baseMean','stat','pvalue','padj')]
global$full_beta_converged <- mcols(lrt)$fullBetaConv
global$reduced_beta_converged <- mcols(lrt)$reducedBetaConv
global$lrt_pvalue_before_convergence_filter <- global$pvalue
invalid_lrt <- is.na(global$full_beta_converged)|!global$full_beta_converged|
               is.na(global$reduced_beta_converged)|!global$reduced_beta_converged
global$pvalue[invalid_lrt] <- NA_real_
global$padj <- p.adjust(global$pvalue,method='BH')
global$degrees_of_freedom <- 4
global <- annotate(global)
write.csv(global[order(global$padj,global$pvalue,na.last=TRUE),],file.path(out,'de_global_ozone_LRT.csv'),row.names=FALSE,na='')
saveRDS(lrt,file.path(out,'dds_global_lrt.rds'))
methods <- list(study='E-MTAB-14917',region='hippocampus',samples=24,
 sex='male, as reported by the primary study',
 biological_replicates='three sequenced animals in each of eight age/genotype/exposure strata',
 age='3 or 17 months at exposure onset, not necessarily at tissue collection',
 read_processing='Initial fastp inspection used the first 200,000 pairs per library and justified consistent paired-end adapter trimming. All reads were then processed with fastp; quality, length, poly-G and duplicate filtering were disabled, while detected adapter dimers were removed by fastp. Full before/after QC is retained per library.',
 quantification='Salmon 2.5.1 selective alignment against GENCODE M25 with whole GRCm38 primary-genome decoys; automatic library type; sequence and GC bias models',
 gene_summarization='tximport: estimated fragments and sample-specific abundance-weighted effective gene lengths; countsFromAbundance=no; DESeqDataSetFromTximport includes length offsets',
 model='negative-binomial GLM, full age * genotype * exposure; group reparameterizations for equivalent stratified coefficients',
 primary_contrasts='O3 versus filtered air within each age x genotype stratum',
 global_test='LRT full age * genotype * exposure versus reduced age * genotype, 4 df',
 filtering='at least ten rounded estimated fragments in at least three of all 24 animals; within-stratum all-zero genes are flagged and their exposure effect/p-values withheld before BH; original coefficient statistics retained in audit columns',
 multiplicity='Fixed count prefilter; adaptive independent filtering disabled. Nonconverged fits are removed before BH within each prespecified contrast and the global LRT separately; supplementary BH across all four valid Wald p-value families',
 shrinkage='apeglm posterior mode and posterior SD for visualization, original Wald p-values preserved',
 sample_handling='All 24 deposited libraries were retained; each sequenced animal is the experimental unit',
 metadata_limits='No library batch or sample RIN covariates deposited; one strain entry differs from all others; human APOE targeted replacement is not allele-specifically quantified by the mouse reference',
 software=list(R=R.version.string,tximport=as.character(packageVersion('tximport')),DESeq2=as.character(packageVersion('DESeq2')),apeglm=as.character(packageVersion('apeglm'))),
 sources=c('https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14917',
 'https://doi.org/10.3390/ijms26062407','https://www.gencodegenes.org/mouse/release_M25.html',
 'https://doi.org/10.1038/nmeth.4197','https://doi.org/10.12688/f1000research.7563.2',
 'https://doi.org/10.1186/s13059-014-0550-8','https://doi.org/10.1093/bioinformatics/bty895'))
write_json(methods,file.path(out,'methods.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
summary <- list(study='E-MTAB-14917',complete=TRUE,n_samples=24,n_input_genes=n_input,
 n_retained_genes=nrow(dds),design_rank=qr(model.matrix(design(dds),as.data.frame(colData(dds))))$rank,
 contrasts=summaries,global_LRT=list(df=4,tested=sum(!is.na(global$pvalue)),fdr_005=sum(global$padj<.05,na.rm=TRUE),
  full_nonconverged=sum(!global$full_beta_converged,na.rm=TRUE),
  reduced_nonconverged=sum(!global$reduced_beta_converged,na.rm=TRUE)),
 pca=list(top_variable_genes=length(selected),PC1_variance=pc$sdev[1]^2/sum(pc$sdev^2),PC2_variance=pc$sdev[2]^2/sum(pc$sdev^2)),
 no_spatial_response_measurement=TRUE)
write_json(summary,file.path(out,'summary.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
message('Ozone differential-expression analysis complete')
