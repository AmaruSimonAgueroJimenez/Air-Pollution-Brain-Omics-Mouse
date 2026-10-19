#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
f <- sub('^--file=', '', grep('^--file=', args, value = TRUE)[1])
f <- gsub('~\\+~', ' ', f)
ROOT <- normalizePath(file.path(dirname(f), '..'))
.libPaths(c(file.path(ROOT, '.r-library'), .libPaths()))
suppressPackageStartupMessages({library(DESeq2);library(apeglm);library(edgeR);library(limma);library(jsonlite)})
OUT <- file.path(ROOT, 'results/replication_pm');dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SOURCE <- file.path(ROOT,'metadata/replication_pm/private_inputs')
write_csv <- function(x,name)write.csv(x,file.path(OUT,name),row.names=FALSE,na='NA')
mode <- commandArgs(trailingOnly=TRUE)
if(length(mode)==0)mode <- 'all'
if(mode %in% c('all','cortex')) {
 message('Cortex: reviewed author-code DESeq2 reimplementation, all20libraries')
 source <- read.delim(gzfile(file.path(ROOT,'data/GSE297423/processed/GSE297423_mouse_gene_count.txt.gz')),check.names=FALSE)
 columns <- grep('^(Control|PFFA|PMA|PMPFFA|PMPFFPMA)',names(source),value=TRUE)
 stopifnot(length(columns)==20)
 counts <- as.matrix(source[,columns]);rownames(counts)<-source$gene_id
 group<-factor(rep(c('Control','PFF','PM','PMPFF','PMPFFPM'),each=4),levels=c('Control','PM','PFF','PMPFF','PMPFFPM'))
 meta<-data.frame(sample=columns,group=group)
 write_csv(meta,'cortex_author_design.csv')
 dds<-DESeqDataSetFromMatrix(countData=counts,colData=data.frame(group=group,row.names=columns),design=~group)
 dds<-dds[rowSums(counts(dds)>=10)>=4,]
 dds<-DESeq(dds,quiet=TRUE,parallel=FALSE)
 # Exact reviewed author settings: default results alpha=0.1 controls independent filtering,
 # followed by discovery threshold padj<0.05; these are distinct quantities.
 summaries<-list()
 for(g in c('PFF','PM','PMPFF','PMPFFPM')) {
  res<-results(dds,contrast=c('group',g,'Control'))
  shr<-lfcShrink(dds,coef=paste0('group_',g,'_vs_Control'),type='apeglm',quiet=TRUE)
  raw<-data.frame(ENSEMBL=rownames(res),as.data.frame(res),check.names=FALSE)
  raw$SYMBOL<-source$gene_name[match(raw$ENSEMBL,source$gene_id)]
  raw$beta_converged<-mcols(dds)$betaConv
  raw$log2FoldChange_shrunk<-shr$log2FoldChange
  raw$lfcSE_shrunk<-shr$lfcSE
  write_csv(raw,paste0('cortex_',g,'_vs_Control.csv'))
  summaries[[g]]<-list(retained_genes=nrow(raw),finite_pvalue=sum(is.finite(raw$pvalue)),finite_padj=sum(is.finite(raw$padj)),FDR05=sum(raw$padj<.05,na.rm=TRUE),up=sum(raw$padj<.05 & raw$log2FoldChange>0,na.rm=TRUE),down=sum(raw$padj<.05 & raw$log2FoldChange<0,na.rm=TRUE),beta_nonconverged=sum(!raw$beta_converged),independent_filter_alpha=.1)
 }
 # Reproduce the actual deposited PCA code, recording its differences from prose methods.
 dge<-DGEList(counts=counts,group=group)
 dge<-dge[rowSums(dge$counts>=10)>=4,,keep.lib.sizes=FALSE]
 dge<-calcNormFactors(dge);lcpm<-cpm(dge,log=TRUE)
 vv<-apply(lcpm,1,var);idx<-order(vv,decreasing=TRUE)[1:500]
 pca<-prcomp(t(lcpm[idx,]),scale.=TRUE)
 write_csv(data.frame(sample=columns,group=group,PC1=pca$x[,1],PC2=pca$x[,2]),'cortex_author_pca.csv')
 write_csv(data.frame(gene=rownames(lcpm)[idx],variance=vv[idx]),'cortex_author_pca_features.csv')
 write_csv(data.frame(gene_id=rownames(dds),counts(dds,normalized=TRUE),check.names=FALSE),'cortex_all20_normalized_counts.csv')
 write_json(list(contrasts=summaries,genes=nrow(dds),PCA=list(PC1=100*pca$sdev[1]^2/sum(pca$sdev^2),PC2=100*pca$sdev[2]^2/sum(pca$sdev^2),top_variable_genes=500,scale=TRUE,normalization='TMM log2CPM'),software=list(R=R.version.string,DESeq2=as.character(packageVersion('DESeq2')),apeglm=as.character(packageVersion('apeglm')),edgeR=as.character(packageVersion('edgeR')))),file.path(OUT,'cortex_refit_summary.json'),auto_unbox=TRUE,pretty=TRUE)
 saveRDS(dds,file.path(OUT,'cortex_author_fit.rds'))
}
if(mode %in% c('all','hypothalamus')) {
 message('Hypothalamus: explicitly unadjusted original-like limma-voom sensitivity, NOT recovered RUVg pipeline')
 meta<-read.csv(file.path(ROOT,'metadata/GSE145840/samples.csv'))
 meta<-meta[match(paste0('GSM',4337093:4337103),meta$gsm),]
 pieces<-lapply(meta$gsm,function(g){p<-Sys.glob(file.path(ROOT,'data/GSE145840/processed',paste0(g,'_Hypothalamus_*.txt.gz')));read.table(gzfile(p),col.names=c('gene','count'))})
 stopifnot(all(vapply(pieces,function(x)!anyDuplicated(x$gene),logical(1))))
 shared_genes<-Reduce(intersect,lapply(pieces,function(x)x$gene))
 union_genes<-Reduce(union,lapply(pieces,function(x)x$gene))
 write_csv(data.frame(gene=setdiff(union_genes,shared_genes),reason='Not measured in every source file; excluded without zero imputation'),'hypothalamus_source_union_only_genes.csv')
 y<-do.call(cbind,lapply(pieces,function(x)x$count[match(shared_genes,x$gene)]));rownames(y)<-shared_genes;colnames(y)<-meta$gsm
 stopifnot(!anyNA(y))
 group<-factor(c(rep('HFD',4),rep('FA',4),rep('PM',3)),levels=c('FA','PM','HFD'))
 # M13 source coding-gene symbols are a separately frozen reference; duplicate gene symbols
 # have a retained mapping decision and are not silently selected by inferential outcome.
 coding_file<-file.path(SOURCE,'M13_protein_coding_symbols.txt')
 coding<-if(file.exists(coding_file))readLines(coding_file) else stop('Missing frozen M13 protein coding annotation')
 eligible<-rownames(y)%in%coding
 write_csv(data.frame(gene=rownames(y),M13_protein_coding=eligible),'hypothalamus_M13_coding_filter.csv')
 y<-y[eligible,]
 # Fixed abundance prefilter uses expression, not group labels or differential evidence.
 dge<-DGEList(counts=y,group=group)
 keep<-rowSums(cpm(dge)>2)>=3
 dge<-dge[keep,,keep.lib.sizes=FALSE];dge<-calcNormFactors(dge)
 design<-model.matrix(~0+group);colnames(design)<-levels(group)
 v<-voom(dge,design,plot=FALSE)
 fit<-lmFit(v,design);ct<-makeContrasts(PM_vs_FA=PM-FA,HFD_vs_FA=HFD-FA,levels=design)
 fit<-eBayes(contrasts.fit(fit,ct))
 summaries<-list()
 for(key in colnames(ct)) {
  tab<-topTable(fit,coef=key,number=Inf,sort.by='none')
  tab<-data.frame(gene=rownames(tab),tab,check.names=FALSE)
  tab$CI95_lower<-tab$logFC-qt(.975,fit$df.total)*fit$stdev.unscaled[,key]*sqrt(fit$s2.post)
  tab$CI95_upper<-tab$logFC+qt(.975,fit$df.total)*fit$stdev.unscaled[,key]*sqrt(fit$s2.post)
  tab$original_like_threshold<-is.finite(tab$adj.P.Val)&tab$adj.P.Val<.05&abs(tab$logFC)>.8&tab$AveExpr>1
  write_csv(tab,paste0('hypothalamus_',key,'_unadjusted_voom.csv'))
  summaries[[key]]<-list(FDR05=sum(tab$adj.P.Val<.05,na.rm=TRUE),FDR05_abslogFC_gt08_AveExpr_gt1=sum(tab$original_like_threshold),tested=nrow(tab))
 }
 write_csv(data.frame(gene=rownames(v$E),v$E,check.names=FALSE),'hypothalamus_original_like_logCPM.csv')
 write_csv(data.frame(sample=colnames(y),label=sub('HYPO-','',meta$title),group=group),'hypothalamus_author_like_design.csv')
 write_json(list(classification='Partial methodological sensitivity only; no original RUVg controls or k recovered',input_genes=length(shared_genes),source_union_genes=length(union_genes),source_union_only_genes=length(setdiff(union_genes,shared_genes)),M13_protein_coding_input=sum(eligible),genes_after_abundance_filter=sum(keep),prefilter='CPM>2 in>=3 of11samples, chosen transparent reconstruction; not specified by original authors',no_RUVg=TRUE,contrasts=summaries,software=list(limma=as.character(packageVersion('limma')),edgeR=as.character(packageVersion('edgeR')))),file.path(OUT,'hypothalamus_refit_summary.json'),auto_unbox=TRUE,pretty=TRUE)
}
capture.output(sessionInfo(),file=file.path(OUT,'R_session_info.txt'))
