#!/usr/bin/env Rscript
# Secondary ozone response-heterogeneity analysis; preserves every primary output.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grep('^--file=', args)])
script <- gsub('~+~', ' ', script, fixed=TRUE)
root <- normalizePath(file.path(dirname(script),'..'),mustWork=TRUE); setwd(root)
.libPaths(c(file.path(root,'.r-library'),.libPaths()))
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
out <- 'results/advanced_omics';dir.create(out,recursive=TRUE,showWarnings=FALSE)
dds <- readRDS('results/ozone/dds_full_wald.rds')
meta <- as.data.frame(colData(dds))
stopifnot(nrow(dds)==18537,ncol(dds)==24,all(table(meta$age,meta$genotype,meta$exposure)==3))
stopifnot(identical(levels(dds$age),c('M3','M17')),identical(levels(dds$genotype),c('E3','E4')),
          identical(levels(dds$exposure),c('FA','O3')))
full <- ~age*genotype*exposure
stopifnot(identical(deparse(design(dds)),deparse(full)))
X <- model.matrix(full,meta);rank_full<-qr(X)$rank;stopifnot(rank_full==8)
write.csv(data.frame(sample_id=rownames(X),X,check.names=FALSE),file.path(out,'full_design_matrix.csv'),row.names=FALSE)
base_disp <- dispersions(dds);base_nf <- normalizationFactors(dds)
annotation <- read.delim('data/references/ozone/genes.tsv',stringsAsFactors=FALSE)
annotate <- function(tab){
 tab$gene_id<-rownames(tab);idx<-match(tab$gene_id,annotation$gene_id)
 tab$gene_symbol<-annotation$gene_symbol[idx];tab$gene_type<-annotation$gene_type[idx]
 tab[,c('gene_id','gene_symbol','gene_type',setdiff(names(tab),c('gene_id','gene_symbol','gene_type')))]
}
formulas <- list(
 genotype_dependent_exposure=~age*genotype+exposure+age:exposure,
 age_dependent_exposure=~age*genotype+exposure+genotype:exposure,
 three_way_interaction=~(age+genotype+exposure)^2)
expected_df<-c(genotype_dependent_exposure=2,age_dependent_exposure=2,three_way_interaction=1)
summaries<-list();lrts<-list();designs<-list()
for(key in names(formulas)){
 reduced<-formulas[[key]];R<-model.matrix(reduced,meta);rank_reduced<-qr(R)$rank
 stopifnot(rank_full-rank_reduced==expected_df[[key]],qr(cbind(X,R))$rank==rank_full)
 message('Secondary LRT: ',key,'; df=',rank_full-rank_reduced)
 fit<-nbinomLRT(dds,full=full,reduced=reduced,maxit=1000,quiet=FALSE)
 stopifnot(identical(dispersions(fit),base_disp),identical(normalizationFactors(fit),base_nf))
 tab<-as.data.frame(results(fit,alpha=.05,independentFiltering=FALSE,cooksCutoff=TRUE))[,c('baseMean','stat','pvalue','padj')]
 tab$full_beta_converged<-mcols(fit)$fullBetaConv;tab$reduced_beta_converged<-mcols(fit)$reducedBetaConv
 tab$pvalue_before_numerical_filter<-tab$pvalue
 invalid<-is.na(tab$full_beta_converged)|!tab$full_beta_converged|is.na(tab$reduced_beta_converged)|!tab$reduced_beta_converged
 tab$pvalue[invalid]<-NA_real_;tab$padj<-p.adjust(tab$pvalue,method='BH')
 tab$degrees_of_freedom<-rank_full-rank_reduced;tab$test<-key
 tab<-annotate(tab);lrts[[key]]<-tab
 write.csv(tab[order(tab$padj,tab$pvalue,na.last=TRUE),],file.path(out,paste0('interaction_LRT_',key,'.csv')),row.names=FALSE,na='')
 summaries[[key]]<-list(df=rank_full-rank_reduced,genes=nrow(tab),tested=sum(is.finite(tab$pvalue)),
  FDR_005=sum(tab$padj<.05,na.rm=TRUE),full_nonconverged=sum(!tab$full_beta_converged,na.rm=TRUE),
  reduced_nonconverged=sum(!tab$reduced_beta_converged,na.rm=TRUE))
 designs[[key]]<-list(full=deparse(full),reduced=deparse(reduced),full_rank=rank_full,reduced_rank=rank_reduced,
  nested=TRUE,df=rank_full-rank_reduced)
 rm(fit);gc()
}
all_lrt<-do.call(rbind,lrts);all_lrt$padj_across_three_LRTs<-p.adjust(all_lrt$pvalue,method='BH')
write.csv(all_lrt,file.path(out,'interaction_LRT_all_tests.csv'),row.names=FALSE,na='')
# Use the saved full Wald fit; linear combinations carry their fitted covariance.
coefnames<-resultsNames(dds)
stopifnot(all(c('ageM17.exposureO3','genotypeE4.exposureO3','ageM17.genotypeE4.exposureO3')%in%coefnames))
contrasts<-list(
 genotype_response_difference_at_M3=c('genotypeE4.exposureO3'),
 genotype_response_difference_at_M17=c('genotypeE4.exposureO3','ageM17.genotypeE4.exposureO3'),
 age_response_difference_in_E3=c('ageM17.exposureO3'),
 age_response_difference_in_E4=c('ageM17.exposureO3','ageM17.genotypeE4.exposureO3'),
 three_way_response_difference=c('ageM17.genotypeE4.exposureO3'))
strata<-paste(meta$age,meta$genotype,sep='_')
zero<-sapply(unique(strata),function(g)rowSums(counts(dds)[,strata==g,drop=FALSE])==0)
relevant<-list(genotype_response_difference_at_M3=c('M3_E3','M3_E4'),
 genotype_response_difference_at_M17=c('M17_E3','M17_E4'),age_response_difference_in_E3=c('M3_E3','M17_E3'),
 age_response_difference_in_E4=c('M3_E4','M17_E4'),three_way_response_difference=colnames(zero))
wtab<-list();wsum<-list();vectors<-list()
for(key in names(contrasts)){
 coefs<-contrasts[[key]];v<-as.numeric(coefnames%in%coefs);names(v)<-coefnames
 res<-results(dds,contrast=v,alpha=.05,independentFiltering=FALSE,cooksCutoff=TRUE)
 tab<-as.data.frame(res);tab$beta_converged<-mcols(dds)$betaConv
 tab$one_required_stratum_all_zero<-apply(zero[,relevant[[key]],drop=FALSE],1,any)
 tab$pvalue_before_numerical_filter<-tab$pvalue
 tab$MLE_before_validity_filter<-tab$log2FoldChange
 tab$SE_before_validity_filter<-tab$lfcSE
 invalid<-is.na(tab$beta_converged)|!tab$beta_converged|tab$one_required_stratum_all_zero
 tab$pvalue[invalid]<-NA_real_;tab$padj<-p.adjust(tab$pvalue,method='BH')
 tab$log2FoldChange[invalid]<-NA_real_;tab$lfcSE[invalid]<-NA_real_
 tab$CI95_low<-tab$log2FoldChange-qnorm(.975)*tab$lfcSE
 tab$CI95_high<-tab$log2FoldChange+qnorm(.975)*tab$lfcSE
 tab$test<-key;tab<-annotate(tab);wtab[[key]]<-tab
 write.csv(tab[order(tab$padj,tab$pvalue,na.last=TRUE),],file.path(out,paste0('interaction_Wald_',key,'.csv')),row.names=FALSE,na='')
 wsum[[key]]<-list(tested=sum(is.finite(tab$pvalue)),FDR_005=sum(tab$padj<.05,na.rm=TRUE),
  required_stratum_all_zero=sum(tab$one_required_stratum_all_zero),nonconverged=sum(!tab$beta_converged,na.rm=TRUE))
 vectors[[key]]<-v
}
all_wald<-do.call(rbind,wtab);all_wald$padj_across_five_Wald_tests<-p.adjust(all_wald$pvalue,method='BH')
write.csv(all_wald,file.path(out,'interaction_Wald_all_tests.csv'),row.names=FALSE,na='')
write.csv(data.frame(test=names(vectors),do.call(rbind,vectors),check.names=FALSE),file.path(out,'Wald_contrast_vectors.csv'),row.names=FALSE)
write_json(list(secondary_exploratory=TRUE,samples=24,genes=18537,full_model=deparse(full),
 original_Wald_results_names=coefnames,designs=designs,LRT=summaries,Wald=wsum,
 dispersions_preserved=TRUE,normalization_factor_offsets_preserved=TRUE,
 multiplicity='BH separately within each three LRT and each five Wald families; supplementary all-three LRT/all-five Wald BH columns form different families.',
 validity='Cook filtering enabled; adaptive independent filtering disabled. Nonconverged fits have p-values withheld before BH. Signed Wald contrasts additionally withhold estimates when a required six-animal age/genotype stratum is all-zero.',
 intervals='Pointwise unshrunken Wald MLE +/-1.9599639845 SE; covariance-aware linear combinations from the full fitted model. Not simultaneous intervals and not per-cell intervals.',
 scope='Secondary response-heterogeneity tests, not replication of original author methods. All24 animals, fixed18537-gene universe and original effective-length offsets; no dispersion re-estimation.',
 documentation='https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions',
 R=R.version.string,DESeq2=as.character(packageVersion('DESeq2'))),file.path(out,'interaction_summary.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
message('Advanced ozone interactions complete')
