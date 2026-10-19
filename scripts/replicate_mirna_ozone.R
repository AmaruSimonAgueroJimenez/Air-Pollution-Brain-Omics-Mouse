#!/usr/bin/env Rscript
# Prespecified candidate reproductions of published O3 comparisons.
# Quantification differs from the original HISAT2/featureCounts pipeline.
args<-commandArgs(trailingOnly=FALSE)
script<-gsub('~+~',' ',sub('^--file=','',args[grep('^--file=',args)]),fixed=TRUE)
root<-normalizePath(file.path(dirname(script),'..'),mustWork=TRUE);setwd(root)
.libPaths(c(file.path(root,'.r-library'),.libPaths()))
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
out<-'results/replication_mirna_ozone';dir.create(out,recursive=TRUE,showWarnings=FALSE)
base<-readRDS('results/ozone/dds_full_wald.rds')
meta<-as.data.frame(colData(base));meta$age<-factor(meta$age,levels=c('M3','M17'))
meta$genotype<-factor(meta$genotype,levels=c('E3','E4'));meta$exposure<-factor(meta$exposure,levels=c('FA','O3'))
genes<-read.delim('data/references/ozone/genes.tsv',stringsAsFactors=FALSE)
models<-list(
 genotype_pooled=list(formula=~genotype,variable='genotype',numerator='E4',denominator='E3',subset=rep(TRUE,24),lfc=.3),
 genotype_additive=list(formula=~age+exposure+genotype,variable='genotype',numerator='E4',denominator='E3',subset=rep(TRUE,24),lfc=.3),
 age_pooled=list(formula=~age,variable='age',numerator='M17',denominator='M3',subset=rep(TRUE,24),lfc=.3),
 exposure_pooled=list(formula=~exposure,variable='exposure',numerator='O3',denominator='FA',subset=rep(TRUE,24),lfc=.1),
 old_ozone_genotype=list(formula=~genotype,variable='genotype',numerator='E4',denominator='E3',subset=meta$age=='M17'&meta$exposure=='O3',lfc=.1))
for(age in c('M3','M17'))for(gt in c('E3','E4')){
 key<-paste(age,gt,'exposure_subset',sep='_')
 models[[key]]<-list(formula=~exposure,variable='exposure',numerator='O3',denominator='FA',subset=meta$age==age&meta$genotype==gt,lfc=.1)
}
summaries<-list()
for(key in names(models)){
 message('Author-style candidate: ',key)
 spec<-models[[key]];selected<-spec$subset
 d<-DESeqDataSetFromMatrix(countData=counts(base)[,selected,drop=FALSE],colData=droplevels(meta[selected,,drop=FALSE]),design=spec$formula)
 # Reuse the validated length/size normalization factors; refit dispersions per candidate model.
 normalizationFactors(d)<-normalizationFactors(base)[,selected,drop=FALSE]
 d<-DESeq(d,test='Wald',betaPrior=FALSE,minReplicatesForReplace=7,quiet=TRUE)
 rr<-results(d,contrast=c(spec$variable,spec$numerator,spec$denominator),alpha=.05,independentFiltering=FALSE,cooksCutoff=TRUE)
 t<-as.data.frame(rr);t$gene_id<-rownames(t);ix<-match(t$gene_id,genes$gene_id);t$gene_symbol<-genes$gene_symbol[ix]
 t$beta_converged<-mcols(d)$betaConv;t$all_zero<-rowSums(counts(d))==0
 t$pvalue_before_convergence_check<-t$pvalue
 invalid<-is.na(t$beta_converged)|!t$beta_converged|t$all_zero
 t$pvalue[invalid]<-NA_real_;t$padj<-p.adjust(t$pvalue,method='BH')
 for(field in c('log2FoldChange','lfcSE','stat'))t[[field]][invalid]<-NA_real_
 t$CI95_lower<-t$log2FoldChange-qnorm(.975)*t$lfcSE;t$CI95_upper<-t$log2FoldChange+qnorm(.975)*t$lfcSE
 t$author_criteria_selected<-!is.na(t$pvalue)&t$pvalue<.05&abs(t$log2FoldChange)>spec$lfc
 t$BH_FDR05<-!is.na(t$padj)&t$padj<.05
 t$comparison<-key
 t<-t[order(t$pvalue,na.last=TRUE),]
 con<-gzfile(file.path(out,paste0('ozone_candidate_',key,'.csv.gz')),'wt');write.csv(t,con,row.names=FALSE,na='');close(con)
 replaced<-if('replace'%in%names(mcols(d)))sum(mcols(d)$replace,na.rm=TRUE)else 0L
 summaries[[key]]<-list(formula=paste(deparse(spec$formula),collapse=''),contrast=paste(spec$numerator,'vs',spec$denominator),samples=colnames(d),n=ncol(d),genes=nrow(d),finite_pvalue=sum(is.finite(t$pvalue)),count_filter='Inherited >=10 rounded estimated fragments in >=3 of all24 source libraries',normalization='Inherited sample-specific effective-length/size normalization factors',dispersion='Reestimated for each candidate design/subset',DESeq2_minReplicatesForReplace=7,genes_replaced_by_DESeq2=replaced,all_zero_genes=sum(t$all_zero),nonconverged_nonzero_genes=sum(!t$beta_converged&!t$all_zero,na.rm=TRUE),nominal_p_cutoff=.05,absolute_log2FC_cutoff=spec$lfc,nominal_selected=sum(t$author_criteria_selected),nominal_up=sum(t$author_criteria_selected&t$log2FoldChange>0,na.rm=TRUE),nominal_down=sum(t$author_criteria_selected&t$log2FoldChange<0,na.rm=TRUE),BH_FDR05=sum(t$BH_FDR05))
 rm(d);gc()
}
write_json(list(study='E-MTAB-14917',complete=TRUE,models=summaries,original_formula_not_reported=TRUE,original_quantifier_not_replicated=TRUE,DESeq2_version=as.character(packageVersion('DESeq2'))),file.path(out,'ozone_candidate_models.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
writeLines(capture.output(sessionInfo()),file.path(out,'R_session_info.txt'))
