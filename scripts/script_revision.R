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



