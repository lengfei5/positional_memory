##########################################################################
##########################################################################
# Project: positional memory
# Script purpose: scripts for revision
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Mar  6 13:29:22 2023
##########################################################################
##########################################################################
rm(list = ls())

require(ggplot2)
require(DESeq2)
require(dplyr)
require(gridExtra)
require(RColorBrewer)

figureDir = paste0('~/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/',
                   'DEVELOPMENTAL CELL/Review/plots_jingkui/') 
tableDir = paste0('~/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/SupTables/')

##########################################
# Reviewer 1 (#1) gene examples of smartseq2 and microarray data 
##########################################
res = readRDS(file = paste0("../results/microarray/Rdata/", 
                            'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_',
                            'normalized_geneSummary_limma.DE.stats_RARB.rds'))

## manual change the gene annotation MEIS3
rownames(res)[grep('AMEX60DD024424', rownames(res))]
rownames(res)[grep('AMEX60DD024424', rownames(res))] = 'MEIS3_AMEX60DD024424' 

ggs = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})

qv.cutoff = 0.05
logfc.cutoff = 1
select = which(res$fdr.max> -log10(qv.cutoff) & abs(res$logFC.max)> 0)

select = which(res$adj.P.Val_mHand.vs.mLA < qv.cutoff & abs(res$logFC_mHand.vs.mLA) > logfc.cutoff|
                 res$adj.P.Val_mHand.vs.mUA < qv.cutoff & abs(res$logFC_mHand.vs.mUA) > logfc.cutoff |
                 res$adj.P.Val_mLA.vs.mUA < qv.cutoff & abs(res$logFC_mLA.vs.mUA) > logfc.cutoff )
# cat(length(select), ' positional genes found \n')
cat(length(select), ' DE genes selected \n')
ggs.sel = ggs[select]

### plot gene expression of microarray data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

gene_examples = c('HOXA13', 'HOXA9', 'HOXD13', 'HOXD9', 'MEIS1', 'MEIS3', 'SHOX2')

for(n in 1:length(gene_examples)){
 
  # n = 1
  kk = which(ggs == gene_examples[n])
  if(length(kk) == 1){
    test = res[kk, c(1:9)]
    test = data.frame(segs = sapply(names(test), function(x){unlist(strsplit(as.character(x), '_'))[1]}), 
                      reps = sapply(names(test), function(x){unlist(strsplit(as.character(x), '_'))[2]}), 
                      signals = as.numeric(test))
    test$segs = factor(test$segs, levels = c('mUA', 'mLA', 'mHand'))
    
    test = data_summary(test, varname = 'signals', groupnames = 'segs')
    
    # Default bar plot
    p<- ggplot(test, aes(x=segs, y=signals, fill=segs)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=signals-sd, ymax=signals+sd), width=.2,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, size = 14), 
            axis.text.y = element_text(angle = 0, size = 14), 
            axis.title =  element_text(size = 14),
            legend.text = element_text(size=12),
            legend.title = element_text(size = 14)
      )+
      labs(x = "") +
      ggtitle(gene_examples[n])
    
    print(p)
    
    ggsave(paste0(figureDir, 'microarray_boxplot.errorbar_geneExample_', gene_examples[n], '.pdf'), 
           width=8, height = 6)
    
    
  }else{
    cat(length(kk), ' row found for gene: ', gene_examples[n], '\n')
  }
}

##########################################
# Reviewer 1 (#2) check the markers of cell types idenetified in Fig2B and FigS4 
# in smartseq2 and also in microarray data
##########################################
Use.microarray.mUA = FALSE
if(Use.microarray.mUA){
  res = readRDS(file = paste0("../results/microarray/Rdata/",
                              'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_',
                              'normalized_geneSummary_limma.DE.stats_RARB.rds'))
  
  ## manual change the gene annotation MEIS3
  rownames(res)[grep('AMEX60DD024424', rownames(res))]
  rownames(res)[grep('AMEX60DD024424', rownames(res))] = 'MEIS3_AMEX60DD024424' 
  
  ggs = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})
  
  res = data.frame(genes = ggs, res[, c(1:3)], stringsAsFactors = FALSE)
  
}
Use.Smartseq.mUA = FALSE
if(Use.Smartseq.mUA){
  
  
}


##########################################
# Reviewer 1 (#3) 
##########################################
source('functions_chipSeq.R')
source('Functions_atac.R')
require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require(pheatmap)
library(tictoc)
library(ggrepel)
library(dplyr)
library(tibble)
library(reshape2)
library(tidyverse)

res = readRDS(paste0("../results/CT_merged_20220328/Rdata/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest",
                     "_all3histM_CT_merged_20220328.rds"))

source('Functions_histM.R')

## select the significant peaks with all three marks
fdr.cutoff = 0.1; logfc.cutoff = 1

select1 = 
  which(((res$adj.P.Val.mLA.vs.mUA.H3K27me3 < fdr.cutoff & abs(res$logFC.mLA.vs.mUA.H3K27me3) > logfc.cutoff) |
          (res$adj.P.Val.mHand.vs.mUA.H3K27me3 < fdr.cutoff & abs(res$logFC.mHand.vs.mUA.H3K27me3) > logfc.cutoff)|
          (res$adj.P.Val.mHand.vs.mLA.H3K27me3 < fdr.cutoff & abs(res$logFC.mHand.vs.mLA.H3K27me3) > logfc.cutoff)
         ) & res$max_rpkm.H3K27me3> 0.6)

select2 = 
  which(((res$adj.P.Val.mLA.vs.mUA.H3K4me1 < fdr.cutoff & abs(res$logFC.mLA.vs.mUA.H3K4me1) > logfc.cutoff) |
          (res$adj.P.Val.mHand.vs.mLA.H3K4me1 < fdr.cutoff & abs(res$logFC.mHand.vs.mLA.H3K4me1) > logfc.cutoff) |
          (res$adj.P.Val.mHand.vs.mUA.H3K4me1 < fdr.cutoff & abs(res$logFC.mHand.vs.mUA.H3K4me1) > logfc.cutoff)
  ) & res$max_rpkm.H3K4me1 > 0.6)

select3 = 
  which(((res$adj.P.Val.mHand.vs.mLA.H3K4me3 < fdr.cutoff & abs(res$logFC.mHand.vs.mLA.H3K4me3) > logfc.cutoff) |
          (res$adj.P.Val.mLA.vs.mUA.H3K4me3 < fdr.cutoff & abs(res$logFC.mLA.vs.mUA.H3K4me3) > logfc.cutoff)|
          (res$adj.P.Val.mHand.vs.mUA.H3K4me3 < fdr.cutoff & abs(res$logFC.mHand.vs.mUA.H3K4me3) > logfc.cutoff)
  ) & res$max_rpkm.H3K4me3 > 0.6)

select = unique(c(select1, select2, select3))


cat(length(select), ' DE H3K27me3, H3K4me1 or H3K4me3 ',  ' \n')

yy0 = res[select, ]
range <- 2.0

conds_histM = c("H3K4me3", "H3K27me3", "H3K4me1")

yy = matrix(NA, nrow = nrow(yy0), ncol = length(conds_histM)*3)
rownames(yy) = rownames(yy0)
nms = c()

for(n in 1:length(conds_histM)) # transform the data
{
  # n = 1
  jj0 = grep(paste0(conds_histM[n], '_m'), colnames(yy0))
  
  test = yy0[, jj0]
  test = test[, grep('_rRep', colnames(test), invert = TRUE)]
  
  test = cal_sample_means(test, conds = paste0(conds_histM[n], c('_mUA', '_mLA', '_mHand')))
  test = t(apply(test, 1, cal_centering))
  test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
  yy[, c((3*n-2):(3*n))] = test
  nms = c(nms, paste0(conds_histM[n], c('_mUA', '_mLA', '_mHand')))
  
  #yy0[ ,jj0] = t(apply(yy0[,jj0], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
  
}

colnames(yy) = nms

saveTables = FALSE
if(saveTables){
  test = data.frame(geneID = get_geneID(rownames(yy)), 
                    gene = get_geneName(rownames(yy)), 
                    yy, 
                    yy0[, c(9:20)], stringsAsFactors = FALSE)
  
  write.csv2(test, file = paste0(tableDir, 'histM_DE_geneCentric_fdr0.1_log2fc.1_rpkm.max.0.6.csv'), 
             row.names = FALSE)
  
}


df = as.data.frame(sapply(colnames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
colnames(df) = 'segments'
rownames(df) = colnames(yy)

sample_colors = c('springgreen4', 'steelblue2', 'gold2')
annot_colors = list(segments = sample_colors)
gaps_col = c(3, 6)

# reorder each cluster
library(dendextend)
library(ggplot2)
source('Functions_histM.R')
nb_clusters = 5
my_hclust_gene <- hclust(dist(yy), method = "complete")
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)

callback = function(hc, mat){
  sv = abs(svd(t(mat))$v[,4])
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

o1 = c()
diffs = yy[, c(4)] - yy[, 6]
cc1 = which(my_gene_col == 2)
cc1 = cc1[order(-diffs[cc1])]
o1 = c(o1, cc1)

cc2 = which(my_gene_col == 1)
cc2 = cc2[order(diffs[cc2])]
o1 = c(o1, cc2)

cc1 = which(my_gene_col == 3)
cc1 = cc1[order(-diffs[cc1])]
o1 = c(o1, cc1)

yy = yy[o1, ]

pheatmap(yy, cluster_rows=TRUE,
         #cutree_rows = 6,
         show_rownames=FALSE, fontsize_row = 4,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
         show_colnames = FALSE,
         scale = 'none',
         cluster_cols=FALSE, 
         annotation_col=df,
         gaps_col = gaps_col,
         legend = TRUE,
         treeheight_row = 15,
         annotation_legend = FALSE, 
         #annotation_colors = annot_colors,
         clustering_callback = callback,
         #breaks = seq(-2, 2, length.out = 8),
         clustering_method = 'complete', 
         cutree_rows = nb_clusters,
         breaks = seq(-range, range, length.out = 8),
         #gaps_row =  c(22, 79),
         legend_labels = FALSE,
         width = 4, height = 12, 
         filename = paste0(figureDir, 'heatmap_histoneMarker_geneCentric_DE_v4.pdf'))

ggs = get_geneName(rownames(yy))
xx = data.frame(names(ggs), ggs, my_gene_col[o1], stringsAsFactors = FALSE)
colnames(xx) = c('gene', 'geneSymbol', 'clusters')
#saveRDS(xx, file = paste0(RdataDir, '/genes_DE.H3K27me3.rds'))
rownames(yy) = ggs

pheatmap(yy, 
         cluster_rows=TRUE,
         #cutree_rows = 4,
         show_rownames=TRUE, fontsize_row = 6,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
         show_colnames = FALSE,
         scale = 'none',
         cluster_cols=FALSE, annotation_col=df,
         gaps_col = gaps_col,
         legend = TRUE,
         treeheight_row = 15,
         annotation_legend = FALSE, 
         #annotation_colors = annot_colors,
         clustering_callback = callback,
         #breaks = seq(-2, 2, length.out = 8),
         clustering_method = 'complete', 
         cutree_rows = nb_clusters,
         breaks = seq(-range, range, length.out = 8),
         #gaps_row =  c(22, 79),
         legend_labels = FALSE,
         width = 6, height = 16, 
         filename = paste0(figureDir, 'heatmap_histoneMarker_geneCentric_DE_geneSymbols_v4.pdf'))


########################################################
########################################################
# Section : gene-centric histone mark analysis for regeneration samples 
# 
########################################################
########################################################
library(edgeR)
library(qvalue)
require(corrplot)
require(pheatmap)
require(RColorBrewer)

gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
gene = GenomicFeatures::genes(amex)

RdataDir = "../results/CT_merged_20220328/Rdata"

conds_histM = c('H3K27me3','H3K4me1', 'H3K4me3')
ll = width(gene)

for(n_histM in 1:length(conds_histM))
{
  # n_histM = 1
  design.sel = readRDS(file = paste0("../results/CT_merged_20220328/Rdata", 
                                     '/design.sels_bc_TMM_combat_TSSgenebody', 
                                     conds_histM[n_histM], '_CT_merged_20220328.rds'))
  
  cpm = readRDS(file = paste0("../results/CT_merged_20220328/Rdata", 
                              '/fpm_bc_TMM_combat_TSSgenebody', conds_histM[n_histM],
                              '_CT_merged_20220328.rds'))
  
  ### select the samples and extract sample means
  conds = c("mUA", "BL5days", "BL9days", "BL13days.prox", "BL13days.dist")
  sample.sels = c();  
  cc = c()
  sample.means = c()
  
  for(n in 1:length(conds)) 
  {
    # n = 1
    #kk = grep(conds[n], colnames(cpm))
    kk = which(design.sel$sample == conds[n] & (design.sel$batch == 'rRep1'| design.sel$batch == 'rRep2'))
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
    if(length(kk)>1) {
      sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
    }else{
      sample.means = cbind(sample.means, cpm[, kk])
    }
    
  }
  colnames(sample.means) = conds
  
  cpm = cpm[, sample.sels]
  design.sel = design.sel[sample.sels, ]
  
  logCPM = cpm
  f = factor(cc, levels= conds)
  
  mod = model.matrix(~ 0 + f)
  colnames(mod) = conds
  
  #To make all pair-wise comparisons between the three groups one could proceed
  fit <- lmFit(logCPM, mod)
  
  contrast.matrix <- makeContrasts(BL5days - mUA, 
                                   BL9days - mUA, 
                                   BL13days.prox - mUA, 
                                   BL13days.dist - mUA,
                                   levels=mod)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  res = data.frame(fit2$p.value)
  colnames(res) = paste0(c('BL5days.vs.mUA', 'BL9days.vs.mUA', 'BL13days.pro.vs.mUA', 'BL13days.dist.vs.mUA'), '.pval')
  
  xx = topTable(fit2, coef = 1, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.BL5days.vs.mUA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  xx = topTable(fit2, coef = 2, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.BL9days.vs.mUA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  xx = topTable(fit2, coef = 3, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.BL13days.pro.vs.mUA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  xx = topTable(fit2, coef = 4, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.BL13days.dist.vs.mUA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  res$pval.mean = apply(as.matrix(res[, grep('P.Value.', colnames(res))]), 1, function(x) return(mean(-log10(x))))
  res$fdr.mean = apply(as.matrix(res[, grep('adj.P.Val', colnames(res))]), 1, function(x) return(mean(-log10(x))))
  res$logFC.mean =  apply(as.matrix(res[, grep('logFC', colnames(res))]), 1, function(x) return(mean(abs(x))))
  
  res$log2fc = apply(sample.means, 1, function(x) max(x) - min(x))
  res$maxs = apply(sample.means, 1, max)
  res$mins = apply(sample.means, 1, min)
  
  source('Functions_histM.R')
  res$length = width(gene)[match(get_geneID(rownames(res)), gene$gene_id)]
  res$length = res$length + 5000
  res$max_rpkm = res$maxs + log2(10^3/res$length)
  
  #
  xx = data.frame(cpm[match(rownames(res), rownames(cpm)), ],  res, stringsAsFactors = FALSE) 
  res = xx
  
  ## select the significant peaks
  fdr.cutoff = 0.05; logfc.cutoff = 1
  
  select = which((res$adj.P.Val.BL5days.vs.mUA < fdr.cutoff & abs(res$logFC.BL5days.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val.BL9days.vs.mUA < fdr.cutoff & abs(res$logFC.BL9days.vs.mUA) > logfc.cutoff)|
                   (res$adj.P.Val.BL13days.pro.vs.mUA < fdr.cutoff & abs(res$logFC.BL13days.pro.vs.mUA) > logfc.cutoff)|
                   (res$adj.P.Val.BL13days.dist.vs.mUA < fdr.cutoff & abs(res$logFC.BL13days.dist.vs.mUA) > logfc.cutoff)
  )
  
  select = select[which(res$max_rpkm[select]>0.)]
  cat(length(select), ' DE ', conds_histM[n_histM],  ' \n')
  
  #res = res[order(res$adj.P.Val.mHand.vs.mUA), ]
  
  # plot_individual_histMarker_withinATACpeak(res)
  
  saveRDS(res, file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_', 
                             conds_histM[n_histM], '_regeneration.rds'))
  
}


##########################################
# plot segment-specific histone marks overlapped with segement-specific atac
# and segment-specific histone marks overlapped with stable atac
# subclustering was performed 
##########################################
Plot_histMarkers_for_positionalATAC = FALSE
if(Plot_histMarkers_for_positionalATAC){
  
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  require(pheatmap)
  require(RColorBrewer)
  library(khroma)
  source('Functions_histM.R')
  
  ## call function for heamtap 
  source('Functions_plots.R')
  conds_histM = c('H3K27me3', 'H3K4me1', 'H3K4me3')
  # subclustering.postional.histM.postioinalAtacPeaks()
  
  ##########################################
  # merge first the histM tables
  ##########################################
  for(n in 1:length(conds_histM)){
    # n = 1
    if(n == 1){
      res = readRDS(file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_', 
                                  conds_histM[n], '_regeneration.rds'))
      colnames(res) = paste0(colnames(res), '.', conds_histM[n])
      
    }else{
      xx = readRDS(file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_', 
                                 conds_histM[n], '_regeneration.rds'))
      colnames(xx) = paste0(colnames(xx), '.', conds_histM[n])
      xx = xx[match(rownames(res), rownames(xx)),]
      res = data.frame(res, xx, stringsAsFactors = FALSE)
    }
  }
  
  saveRDS(res, file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_all3histM_regeneration.rds'))
  
  ##########################################
  # select the genes based on the H3K27me3 
  ##########################################
  res = readRDS(file =paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_all3histM_regeneration.rds'))
  source('Functions_histM.R')
  
  ## select dynamics peaks only based on H3K27me3
  fdr.cutoff = 0.1; logfc.cutoff = 1
  res = res[, grep('.H3K27me3', colnames(res))]
  #colnames(res) = gsub('.H3K27me3', '', colnames(res))
  
  select = which((res$adj.P.Val.BL5days.vs.mUA.H3K27me3 < fdr.cutoff & 
                    abs(res$logFC.BL5days.vs.mUA.H3K27me3) > logfc.cutoff) |
                   (res$adj.P.Val.BL9days.vs.mUA.H3K27me3 < fdr.cutoff & 
                      abs(res$logFC.BL9days.vs.mUA.H3K27me3) > logfc.cutoff)|
                   (res$adj.P.Val.BL13days.pro.vs.mUA.H3K27me3 < fdr.cutoff & 
                      abs(res$logFC.BL13days.pro.vs.mUA.H3K27me3) > logfc.cutoff)|
                   (res$adj.P.Val.BL13days.dist.vs.mUA.H3K27me3 <fdr.cutoff &
                      abs(res$logFC.BL13days.dist.vs.mUA.H3K27me3) > logfc.cutoff)
  )
  
  select = select[which(res$max_rpkm.H3K27me3[select]>0.6)]
  cat(length(select), ' DE H3K27me3 ',  ' \n')
  
  yy0 = res[select, ]
  
  range <- 2.
  
  #conds_histM = c("H3K4me3", "H3K27me3", "H3K4me1")
  conds_histM = c("H3K27me3")
  conds = c("mUA", "BL5days", "BL9days", "BL13days.prox", "BL13days.dist")
  yy = matrix(NA, nrow = nrow(yy0), ncol = length(conds_histM)*length(conds))
  rownames(yy) = rownames(yy0)
  nms = c()
  
  for(n in 1:length(conds_histM)) # transform the data
  {
    # n = 1
    jj0 = grep(paste0('^', conds_histM[n], '_'), colnames(yy0))
    
    test = yy0[, jj0]
    test = test[, grep('_rRep', colnames(test), invert = FALSE)]
    
    test = cal_sample_means(test, conds = paste0(conds_histM[n], '_', conds))
    test = t(apply(test, 1, cal_centering))
    test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
    yy[, c((length(conds)*(n-1) +1):(length(conds)*n))] = test
    nms = c(nms, paste0(conds_histM[n], "_", conds))
    #yy0[ ,jj0] = t(apply(yy0[,jj0], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
    
  }
  
  colnames(yy) = nms
  
  if(saveTables){
    test = data.frame(geneID = get_geneID(rownames(yy)), 
                      gene = get_geneName(rownames(yy)), 
                      yy, 
                      yy0[, c(11:26)], stringsAsFactors = FALSE)
    
    write.csv2(test, file = paste0(figureDir, 'H3K27me3_geneCentric_regeneration_fdr0.1_log2fc.1_rpkm.max.0.6.csv'), 
               row.names = FALSE)
    
  }
  
  
  df = as.data.frame(sapply(colnames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'condition'
  rownames(df) = colnames(yy)
  
  sample_colors = sample_colors = c('darkblue', 'springgreen', 'springgreen3', 'gold2', 'red')
  annot_colors = list(segments = sample_colors)
  #gaps_col = c(3, 6)
  
  # reorder each cluster
  library(dendextend)
  library(ggplot2)
  source('Functions_histM.R')
  nb_clusters = 5
  my_hclust_gene <- hclust(dist(yy), method = "complete")
  my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
  
  callback = function(hc, mat){
    sv = abs(svd(t(mat))$v[,4])
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  # o1 = c()
  # diffs = yy[, c(4)] - yy[, 6]
  # cc1 = which(my_gene_col == 2)
  # cc1 = cc1[order(-diffs[cc1])]
  # o1 = c(o1, cc1)
  # 
  # cc2 = which(my_gene_col == 1)
  # cc2 = cc2[order(diffs[cc2])]
  # o1 = c(o1, cc2)
  # 
  # cc1 = which(my_gene_col == 3)
  # cc1 = cc1[order(-diffs[cc1])]
  # o1 = c(o1, cc1)
  # 
  # yy = yy[o1, ]
  library(khroma)
  nb_breaks = 7
  sunset <- colour("sunset")
  PRGn <- colour("PRGn")
  cols = rev(PRGn(nb_breaks))
  
  pheatmap(yy, cluster_rows=TRUE,
           #cutree_rows = 4,
           show_rownames=FALSE, fontsize_row = 6,
           #color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
           color = cols,
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           #gaps_col = gaps_col,
           legend = TRUE,
           treeheight_row = 15,
           annotation_legend = FALSE, 
           #annotation_colors = annot_colors,
           #clustering_callback = callback,
           #breaks = seq(-2, 2, length.out = 8),
           clustering_method = 'complete', cutree_rows = nb_clusters,
           breaks = seq(-range, range, length.out = 8),
           gaps_row =  c(22, 79),
           legend_labels = FALSE,
           width = 4, height = 10, 
           filename = paste0(figureDir, 'heatmap_H3K27me3_geneCentric_regeneration_v1.pdf'))
  
  
  ggs = get_geneName(rownames(yy))
  #xx = data.frame(names(ggs), ggs, my_gene_col[o1], stringsAsFactors = FALSE)
  #colnames(xx) = c('gene', 'geneSymbol', 'clusters')
  #saveRDS(xx, file = paste0(RdataDir, '/genes_DE.H3K27me3.rds'))
  rownames(yy) = ggs
  
  pheatmap(yy, cluster_rows=TRUE,
           #cutree_rows = 4,
           show_rownames=TRUE, fontsize_row = 4,
           #show_rownames=FALSE, fontsize_row = 6,
           #color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
           color = cols,
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           #gaps_col = gaps_col,
           legend = TRUE,
           treeheight_row = 15,
           annotation_legend = FALSE, 
           #annotation_colors = annot_colors,
           #clustering_callback = callback,
           #breaks = seq(-2, 2, length.out = 8),
           clustering_method = 'complete', cutree_rows = nb_clusters,
           breaks = seq(-range, range, length.out = 8),
           gaps_row =  c(22, 79),
           legend_labels = FALSE,
           width = 5, height = 20, 
           filename = paste0(figureDir, 'heatmap_H3K27me3_geneCentric_regeneration_geneSymbols_v1.pdf'))
}


