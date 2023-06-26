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
require(tibble)

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
  
  #gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
  #amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
  #gene = GenomicFeatures::genes(amex)
  #ll = width(gene)
  res = readRDS(file = paste0("../results/RNAseq_data_used/Rdata", 
                              "/matureSamples_cpm_DEgenes_8selectedSamples.batch4_v47.hox.patch.rds"))
  
  res$gene[which(res$geneID == 'AMEX60DD024424')] = 'MEIS3'
  res$gene[which(res$geneID == 'AMEX60DD002658')] = 'OTOS'
  #xx = res[, 3]
  
  ll = read.table(paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/',
                         'Data/R14026_Quantseq/featurecounts_Q10/204810_S6_R1_001.sort_umiDedup_featureCounts.txt'),
                  header = TRUE)
  res$length = ll$Length[match(res$geneID, ll$Geneid)]
  
  markers = openxlsx::read.xlsx('../data/mUA_celltypes_markerGenes.xlsx')
  markers = markers[, c(6, 2)]
  
  markers = rbind(markers, c('Periskeletal', 'COL8A1'))
  markers = rbind(markers, c('skeletal', 'CNMD'))
  markers = rbind(markers, c('skeletal', 'LECT1'))
  #markers = rbind(markers, c('skeletal', 'SULF2'))
  markers = rbind(markers, c('skeletal', 'OTOS'))
  markers = rbind(markers, c('cycling_cells', 'CDK1'))
  markers = rbind(markers, c('cycling_cells', 'CCNB1'))
  
  
  markers = markers[which(!is.na(match(markers$Gene, res$gene))), ]
  
  markers = markers[which(markers$Gene != 'PINLYP'), ]
  markers = markers[which(markers$Gene != 'VWDE'), ]
  markers = markers[which(markers$Gene != 'THBS4'), ]
  markers = markers[which(markers$Gene != 'HPD'), ]
  
  mm = match(markers$Gene, res$gene)
  xx = res[mm, c(3:4, 10)]
  for(n in 1:nrow(xx))
  {
    xx[n, 1] = xx[n, 1] + log2(10^3/xx[n, 3])
  }
  
  markers$rpkm = apply(xx[, c(1,2)], 1, mean)
  markers$sd = apply(xx[, c(1,2)], 1, sd)
  markers = markers[, c(2:4, 1)]
  
  # select representative markers
  gene_examples = c('COL4A2', 'COL4A1', 
                    "KLF5", 'IGFBP3',
                    'TWIST2', 'F13A1','COL1A1',
                    'TNMD', "SULF2", 'COL2A1', 'COL11A1',
                    'CNMD', 'OTOS', 'CDK1', 'CCNB1')
  
  markers = markers[!is.na(match(markers$Gene, gene_examples)), ]
  
  table(markers$Gene) # some markers for multiple groups: COL2A1, COL11A1, TWIST2
  markers$celltypes[which(markers$Gene == 'COL2A1')] = 'Tendon.Periskeletal'
  markers$celltypes[which(markers$Gene == 'COL11A1')] = 'Tendon.Periskeletal'
  markers$celltypes[which(markers$Gene == 'TWIST2')] = 'fCT_IV.V'
  markers = markers[match(unique((markers$Gene)), markers$Gene), ]
  
  markers$celltypes = factor(markers$celltypes, levels = c('fCT_I', "fCT_II", 'fCT_III', "fCT_IV.V", 
                                                           'fCT_V', 'Tendon', 'Periskeletal', 'Tendon.Periskeletal',
                                                           'skeletal', 'cycling_cells'))
  
  markers$Gene = factor(markers$Gene, levels = unique(markers$Gene[order(markers$celltypes)]))
  #markers$position = c(1:nrow(markers))
  
  
  as_tibble(markers) %>% 
    group_by_(~ celltypes) %>%
    top_n(n = 2) %>%
    ggplot(aes(x=Gene, y=rpkm, fill = celltypes)) + 
    geom_bar(stat="identity", color="black", position ="dodge") +
    #geom_errorbar(aes(ymin=rpkm-sd, ymax=rpkm+sd), width=.2,
    #              position=position_dodge(.9)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 12), 
          axis.text.y = element_text(angle = 0, size = 12), 
          axis.title =  element_text(size = 14),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14)
    ) +
    scale_fill_manual(values=c(colorRampPalette(rev(brewer.pal(n = 8, name ="OrRd")))(5),
                               "springgreen", "blue", "#E69F00", "#56B4E9", "#999999")) + 
    labs(x = "", y = 'log2(RPKM)')
  
  ggsave(paste0(figureDir, 'main_markerGenes_cellcomposition_mature.CT.pdf'), 
         width=8, height = 6)
  
  ##########################################
  # ### compare Tobie's pseudo bulk RNA-seq vs Akane's smartseq2 mUA 
  ##########################################
  load(file = paste0("../results/Rxxxx_R10723_R11637_R12810_atac/Rdata", 
                     '/Gerber_2018_Fluidigm_C1.Rdata'))
  
  annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  
  counts = scRNAseq.counts[, -1]
  rownames(counts) = scRNAseq.counts$ENSEMBL_ID
  rm(scRNAseq.counts)
  
  # convert gene names to gene symbols
  mm = match(rownames(counts), annot$geneID)
  ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
  rownames(counts)[!is.na(mm)] = ggs[!is.na(mm)]
  
  accs = sapply(colnames(counts), function(x) unlist(strsplit(as.character(x), '_'))[1])
  
  mm = match(metadata$metadata.run_accession, accs)
  metadata = metadata[which(!is.na(mm)), ]
  counts = counts[, mm[which(!is.na(mm))]]
  colnames(counts) = as.character(metadata$metadata.run_accession)
  
  metadata$condition = metadata$timepoint
  metadata$sample = as.character(metadata$metadata.run_accession)
  
  #raw = as.matrix(counts[ ,-1])
  #rownames(raw) = counts$ENSEMBL_ID
  
  metadata$batch = sapply(metadata$cell_id, function(x) unlist(strsplit(as.character(x), '[.]'))[1])
  rownames(metadata) = metadata$cell_id
  
  colnames(counts) = rownames(metadata)
  
  require(DESeq2)
  raw = counts
  rm(counts)
  
  conds = c('0dpa', '3dpa', '5dpa', '8dpa', '11dpa', '18dpa', '1apa',  'Stage40', 'Stage44')
  pseudo = c()
  cc = c()
  plates = c()
  for(n in 1:length(conds))
  {
    cat(n, ' : ', conds[n], ' -- ')
    jj = which(metadata$condition == conds[n])
    
    mm = match(metadata$cell_id[jj], colnames(raw))
    if(length(which(is.na(mm)))>0) {
      cat('some cells lost \n')
    }else{
      cat(length(mm), ' cell found \n')
      
      bcs = unique(metadata$batch[jj])
      for(bc in bcs){
        cat('-- ', bc, '--')
        xx = as.matrix(raw[, jj[which(metadata$batch[jj] == bc)]])
        cat(ncol(xx), ' cells \n')
        xx[is.na(xx)] = 0
        pseudo = cbind(pseudo, apply(xx, 1, sum))
        cc = c(cc, conds[n])
        plates = c(plates, bc)
      }
    }
  }
  
  rownames(pseudo) = rownames(raw)
  design = data.frame(condition = cc, plates = plates, stringsAsFactors = FALSE)
  colnames(pseudo) = paste0(design$condition, '_', design$plates)
  
  dds <- DESeqDataSetFromMatrix(pseudo, DataFrame(design), design = ~ condition)
  
  # saveRDS(dds, file = paste0(Rdata.smartseq2, 'Geber_pooledscRNA_muA_regeneration_dev.rds'))
  
  jj = which(design$condition == '0dpa')
  yy = as.matrix(pseudo[, jj])
  yy = apply(yy, 1, sum)
  
  res$pseudo = yy[match(rownames(res), names(yy))]
  res$pseudo = log2(res$pseudo/sum(res$pseudo, na.rm = TRUE)*10^6 + 2^-6)
  
  mm = match(gene_examples, res$gene)
  plot(res$mUA_161512[mm], res$pseudo[mm])
  abline(0, 1, col = 'red')
  
  mm = match(markers$Gene, res$gene)
  markers$mUA = res$mUA_161512[mm]
  markers$pseudo = res$pseudo[mm]
  require("ggrepel")
  markers = as_tibble(markers) 
    #group_by_(~ celltypes) %>%
    #top_n(n = 2) %>%
  ggplot(data = markers, aes(x=mUA, y=pseudo, color=celltypes, label = Gene)) + 
    #geom_bar(stat="identity", color="black", position ="dodge") +
    #geom_errorbar(aes(ymin=rpkm-sd, ymax=rpkm+sd), width=.2,
    #              position=position_dodge(.9)) +
    geom_point(size = 5)+ 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14), 
          axis.text.y = element_text(angle = 0, size = 14), 
          axis.title =  element_text(size = 14),
          legend.text = element_text(size=14),
          legend.title = element_text(size = 14)
    ) +
    geom_text_repel( box.padding = unit(0.5, "lines"), point.padding = unit(0.4, "lines")) +
    scale_color_manual(values=c(colorRampPalette(rev(brewer.pal(n = 8, name ="OrRd")))(5),
                               "springgreen", "blue", "#E69F00", "#56B4E9", "#999999")) + 
    labs(x = "Smart-seq2 (mUA, log2cpm)", y = 'Pseudo-bulk of Gerber et al.') +
    geom_abline(intercept = 0, slope = 1, color = 'black')
    
  ggsave(paste0(figureDir, 'mUA_Smartseq2_vs_Gerber.pseudobulk_markerGenes.pdf'), 
         width=8, height = 5.5)
  
  
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


##########################################
# Reviewer 2 (#4) gene examples of smartseq2 in segments and in regeneration
##########################################
gene_examples = c('HOXA13', 'HOXA9', 'HOXD13', 'HOXD9', 'MEIS1', 'MEIS2', 'MEIS3', 'SHOX2', 'SHOX')

res = readRDS(file = paste0("../results/RNAseq_data_used/Rdata", 
                            "/matureSamples_cpm_DEgenes_8selectedSamples.batch4_v47.hox.patch.rds"))

jj = grep('AMEX60DD024424', rownames(res))
rownames(res)[jj] = 'MEIS3_AMEX60DD024424'
res$gene[jj] = 'MEIS3'

mm = match(gene_examples, res$gene)
aa = res[mm, c(8:9, 1:4)]

res = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/', 
                            'smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked_',
                            'RNAseq_data_used_20220408', '.rds'))

res = res[, c(1:12)]
jj = grep('AMEX60DD024424', rownames(res))
rownames(res)[jj] = 'MEIS3_AMEX60DD024424'

mm = match(rownames(aa), rownames(res))
aa = data.frame(aa, res[mm, ])

res = aa
# standardized the sample names
colnames(res)[7:ncol(res)] = gsub('Mature_UA', 'mUA', colnames(res)[7:ncol(res)])
colnames(res)[7:ncol(res)] = sapply(colnames(res)[7:ncol(res)], function(x){
  test = unlist(strsplit(as.character(x), '_'));
  test = test[1:(length(test)-1)]
  return(paste0(paste0(test[1:(length(test)-1)], collapse = '.'), '_', test[length(test)]))
})

ll = read.table(paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/',
                       'Data/R14026_Quantseq/featurecounts_Q10/204810_S6_R1_001.sort_umiDedup_featureCounts.txt'),
                header = TRUE)
res$length = ll$Length[match(res$geneID, ll$Geneid)]

for(n in 1:nrow(res))
{
  res[n, 3:18] = res[n, 3:18] + log2(10^3/res$length[n])
}

res = res[, -which(colnames(res) == 'length')]

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

library(scales)
show_col(hue_pal()(5))
cols = hue_pal()(6)
cols = cols[c(2:1,3:6)]

for(n in 1:length(gene_examples)){
  
  # n = 1
  cat(n, '--', gene_examples[n], '\n')
  kk = which(res$gene == gene_examples[n])
  
  if(length(kk) == 1){
    test = res[kk, c(3:4, 7:ncol(res))]
    test = data.frame(segs = sapply(names(test), function(x){unlist(strsplit(as.character(x), '_'))[1]}), 
                      reps = sapply(names(test), function(x){unlist(strsplit(as.character(x), '_'))[2]}), 
                      signals = as.numeric(test))
    test$segs = factor(test$segs, levels = c('mHand', 'mUA', 'BL.UA.5days', 'BL.UA.9days', 'BL.UA.13days.distal',
                                             'BL.UA.13days.proximal'))
    
    test = data_summary(test, varname = 'signals', groupnames = 'segs')
    
    # Default bar plot
    p<- ggplot(test, aes(x=segs, y=signals, fill=segs)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=signals-sd, ymax=signals+sd), width=.2,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, size = 14), 
            axis.text.y = element_text(angle = 0, size = 14), 
            axis.title =  element_text(size = 14),
            legend.text = element_text(size=12),
            legend.title = element_text(size = 14)
      ) + 
      scale_fill_manual(values=cols)+
      labs(x = "", y = 'log2(FPKM)') +
      ggtitle(gene_examples[n])
    
    print(p)
    
    ggsave(paste0(figureDir, 'smartseq_mature.segments.regeneration_boxplot.errorbar_', 
                  gene_examples[n], '.pdf'), 
           width=8, height = 6)
  }else{
    cat(length(kk), ' row found for gene: ', gene_examples[n], '\n')
  }
}

########################################################
########################################################
# Review 3 (#2) 
# gene-centric histone mark (H3K27me3, H3K4me1 and H3K4me3) analysis for regeneration samples 
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

######################
## plot histone marks
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
  RdataDir = "../results/CT_merged_20220328/Rdata"
  
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
  #res = res[, grep('.H3K27me3', colnames(res))]
  #colnames(res) = gsub('.H3K27me3', '', colnames(res))
  
  select1 = which((res$adj.P.Val.BL5days.vs.mUA.H3K27me3 < fdr.cutoff & 
                    abs(res$logFC.BL5days.vs.mUA.H3K27me3) > logfc.cutoff) |
                   (res$adj.P.Val.BL9days.vs.mUA.H3K27me3 < fdr.cutoff & 
                      abs(res$logFC.BL9days.vs.mUA.H3K27me3) > logfc.cutoff)|
                   (res$adj.P.Val.BL13days.pro.vs.mUA.H3K27me3 < fdr.cutoff & 
                      abs(res$logFC.BL13days.pro.vs.mUA.H3K27me3) > logfc.cutoff)|
                   (res$adj.P.Val.BL13days.dist.vs.mUA.H3K27me3 <fdr.cutoff &
                      abs(res$logFC.BL13days.dist.vs.mUA.H3K27me3) > logfc.cutoff)
  )
  select1 = select1[which(res$max_rpkm.H3K27me3[select1]>0.6)]
  
  select2 = which((res$adj.P.Val.BL5days.vs.mUA.H3K4me1 < fdr.cutoff & 
                     abs(res$logFC.BL5days.vs.mUA.H3K4me1) > logfc.cutoff) |
                    (res$adj.P.Val.BL9days.vs.mUA.H3K4me1 < fdr.cutoff & 
                       abs(res$logFC.BL9days.vs.mUA.H3K4me1) > logfc.cutoff)|
                    (res$adj.P.Val.BL13days.pro.vs.mUA.H3K4me1 < fdr.cutoff & 
                       abs(res$logFC.BL13days.pro.vs.mUA.H3K4me1) > logfc.cutoff)|
                    (res$adj.P.Val.BL13days.dist.vs.mUA.H3K4me1 <fdr.cutoff &
                       abs(res$logFC.BL13days.dist.vs.mUA.H3K4me1) > logfc.cutoff)
  )
  select2 = select2[which(res$max_rpkm.H3K4me1[select2]>0.6)]
  
  select3 = which((res$adj.P.Val.BL5days.vs.mUA.H3K4me3 < fdr.cutoff & 
                     abs(res$logFC.BL5days.vs.mUA.H3K4me3) > logfc.cutoff) |
                    (res$adj.P.Val.BL9days.vs.mUA.H3K4me3 < fdr.cutoff & 
                       abs(res$logFC.BL9days.vs.mUA.H3K4me3) > logfc.cutoff)|
                    (res$adj.P.Val.BL13days.pro.vs.mUA.H3K4me3 < fdr.cutoff & 
                       abs(res$logFC.BL13days.pro.vs.mUA.H3K4me3) > logfc.cutoff)|
                    (res$adj.P.Val.BL13days.dist.vs.mUA.H3K4me3 <fdr.cutoff &
                       abs(res$logFC.BL13days.dist.vs.mUA.H3K4me3) > logfc.cutoff)
  )
  select3 = select3[which(res$max_rpkm.H3K4me3[select3]>0.6)]
  
  select = unique(c(select1, select2, select3))
  cat(length(select), ' DE histone marks ',  ' \n')
  
  yy0 = res[select, ]
  
  range <- 2.
  
  conds_histM = c("H3K4me3", "H3K27me3", "H3K4me1")
  #conds_histM = c("H3K27me3")
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
    
    write.csv2(test, file = paste0(figureDir, 'All.three.histoneMarks_geneCentric_regeneration_',
                                   'fdr0.1_log2fc.1_rpkm.max.0.6.csv'), 
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
  nb_clusters = 6
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
  #cols = rev(PRGn(nb_breaks))
  cols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(7)
  
  pheatmap(yy, cluster_rows=TRUE,
           #cutree_rows = 4,
           show_rownames=FALSE, fontsize_row = 6,
           #color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
           color = cols,
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           #gaps_col = gaps_col,
           gaps_col = seq(5, ncol(yy), by = 5),
           legend = TRUE,
           treeheight_row = 20,
           annotation_legend = FALSE, 
           #annotation_colors = annot_colors,
           #clustering_callback = callback,
           #breaks = seq(-2, 2, length.out = 8),
           clustering_method = 'complete', cutree_rows = nb_clusters,
           breaks = seq(-range, range, length.out = 8),
           #gaps_row =  c(22, 79),
           legend_labels = FALSE,
           width = 4, height = 10, 
           filename = paste0(figureDir, 'heatmap_all.three.histoneMarks_geneCentric_regeneration_v1.pdf'))
  
  ggs = get_geneName(rownames(yy))
  #xx = data.frame(names(ggs), ggs, my_gene_col[o1], stringsAsFactors = FALSE)
  #colnames(xx) = c('gene', 'geneSymbol', 'clusters')
  #saveRDS(xx, file = paste0(RdataDir, '/genes_DE.H3K27me3.rds'))
  rownames(yy) = ggs
  
  pheatmap(yy, cluster_rows=TRUE,
           #cutree_rows = 4,
           show_rownames=TRUE, fontsize_row = 3,
           #show_rownames=FALSE, fontsize_row = 6,
           #color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
           color = cols,
           show_colnames = FALSE,
           gaps_col = seq(5, ncol(yy), by = 5),
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           #gaps_col = gaps_col,
           legend = TRUE,
           treeheight_row = 20,
           annotation_legend = FALSE, 
           #annotation_colors = annot_colors,
           #clustering_callback = callback,
           #breaks = seq(-2, 2, length.out = 8),
           clustering_method = 'complete', cutree_rows = nb_clusters,
           breaks = seq(-range, range, length.out = 8),
           gaps_row =  c(22, 79),
           legend_labels = FALSE,
           width = 8, height = 25, 
           filename = paste0(figureDir, 'heatmap_all.three.histoneMarks_geneCentric_regeneration_geneSymbols_v1.pdf'))
  
  
  ##########################################
  # test GO enrichment for each clusters 
  ##########################################
  ggs = get_geneName(names(my_gene_col))
  
  jj = grep('AMEX60', ggs, invert = TRUE)
  
  genes_clusters = my_gene_col[jj]
  ggs = ggs[jj]
  
  library(enrichplot)
  library(clusterProfiler)
  library(ggplot2)
  library(stringr)
  #library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  clean_geneNames = function(gg.expressed)
  {
    gg.expressed = unique(unlist(lapply(gg.expressed, 
                                        function(x) { x = unlist(strsplit(as.character(x), '_'));  
                                        return(x[length(x)])})))
    gg.expressed = unique(annot$gene.symbol.toUse[match(gg.expressed, annot$geneID)])
    gg.expressed = gg.expressed[which(gg.expressed != '' & gg.expressed != 'N/A' & !is.na(gg.expressed))]
    
    return(gg.expressed)
  }
  
  # background
  bgs0 = unique(clean_geneNames(rownames(res)))
  
  xx0 = unique(annot$gene.symbol.toUse)
  xx0 = xx0[which(xx0 != '' & xx0 != 'N/A' & !is.na(xx0))]
  bgs = unique(xx0)
  
  for(cluster.id in c(1:3)){
    # cluster.id = 4
    gg.expressed = unique(ggs[which(genes_clusters == cluster.id)])
    cat('# of genes --', length(gg.expressed), '\n')
    
    gg.expressed = firstup(tolower(gg.expressed))
    bgs = firstup(tolower(bgs))
    bgs0 = firstup(tolower(bgs0))
    
    gene.df <- bitr(gg.expressed, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)
    head(gene.df)
    
    #bgs.df <- bitr(bgs, fromType = "SYMBOL",
    #               toType = c("ENSEMBL", "ENTREZID"),
    #               OrgDb = org.Mm.eg.db)
    #head(bgs.df)
    
    bgs0.df <- bitr(bgs0, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)
    
    #pval.cutoff = 0.05
    ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                     universe     = bgs0.df$ENSEMBL,
                     #universe     = bgs.df$ENSEMBL,
                     #OrgDb         = org.Hs.eg.db,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1, 
                     minGSSize = 3)
    
    #head(ego)
    #barplot(ego) + ggtitle("Go term enrichment for positional genes")
    pdfname = paste0(figureDir, 'GOterm_DE_histoneMarkers_cluster_', cluster.id,  '.pdf')
    pdf(pdfname, width = 12, height = 6)
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
    
    dotplot(ego, showCategory=10) + ggtitle(paste0("genes in cluster ", cluster.id))
    
    dev.off()
    
    
  }
  
  
}


##########################################
# change the color of Fig S1B 
##########################################
Process.deeptools.heatmapTable = FALSE
if(Process.deeptools.heatmapTable){
  dpt = read.table(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/', 
                                 'Data/atacseq_using/heatmap_deeptools/heatmaps_all/heatmap_allchrFeatures_matrix4save.txt'),
                   comment.char = "#", sep = '\t', header = TRUE, skip = 2)
  
  peaks = c(rep(1, 76), 
            rep(2, 712), 
            rep(3, 292),
            rep(4, 103), 
            rep(5, 31), 
            rep(6, 32))
  
  samples = colnames(dpt)
  samples = sapply(samples, function(x){x = unlist(strsplit(as.character(x), '_')); paste0(x[1:3], collapse = '_')})          
  
  cluster_order = c(6, 1, 5, 3, 4, 2)
  
  index = c()
  gaps.row = c() 
  for(n in 1:length(cluster_order)) ## compute row gaps as atac-seq peaks
  {
    jj = which(peaks == cluster_order[n])
    index = c(index, jj)
    if(n == 1)  {gaps.row = c(gaps.row, length(jj))
    }else{
      if(n < length(cluster_order)){
        gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(peaks == cluster_order[n])))
      }
    }
  }
  
  dpt = dpt[index, ]
  
  ##########################################
  # select representative samples 
  ##########################################
  library(ComplexHeatmap)
  library("RColorBrewer")
  library("circlize")
  
  ## select the sample to use for plotting
  idx =  c(grep('136164', samples), 
           grep('161521', samples), 
           grep('74940', samples), 
           intersect(c(grep('H3K27me3', samples), grep('H3K4me1', samples), grep('H3K4me3', samples)), grep('mRep2', samples))
  )
  samples = samples[idx]
  dpt = dpt[, idx]
  
  ## sort by H3K27me3 for each cluster
  index = grep('H3K27me3_mUA_mRep2', samples)
  ss = apply(as.matrix(dpt[, index]), 1, mean)
  new_order = c()
  for(c in cluster_order)
  {
    kk = which(peaks==c)
    new_order = c(new_order, kk[order(-ss[kk])])
  }
  
  features = c('atac','H3K4me3', 'H3K27me3', 'H3K4me1')
  
  for(n in 1:length(features))
  {
    # n = 3
    cat(n, ' -- ', features[n], '\n')
    
    index = grep(features[n], samples)
    test = as.matrix(dpt[new_order, index])*20
    
    splits = factor(peaks[new_order], levels = cluster_order)
    sample.sel = samples[index]
    sample.sel = sapply(sample.sel, function(x) unlist(strsplit(as.character(x), '_'))[2])
    ha <- HeatmapAnnotation(samples = sample.sel)
    
    unique(sample.sel)
    
    #col_fun = colorRamp2(c(0, 0.5, 1), c("#377EB8", "white", "#E41A1C"))
    quantile(test, c(0, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 0.8, 0.85,  0.90,  0.95, 0.99, 1))
    
    if(features[n] == 'atac') {range = 2.5; 
    #cols = colorRamp2(c(0, 0.5, 1.5, range), (brewer.pal(n=4, name="Greys")))
    cols = colorRamp2(c(0, range), colors = c('white', 'black'))
    #cols = colorRamp2(c(0, 2, 4), c("blue", "white", "red"));
    }
    if(features[n] == 'H3K4me1'){range = 7; cols = colorRamp2(c(0, range), colors = c('white', 'darkgoldenrod2'))}
    #if(features[n] == 'H3K27me3') {range = 5; cols = colorRamp2(c(0, range), colors = c('white', 'deepskyblue1'))}
    if(features[n] == 'H3K27me3') {range = 5; cols = colorRamp2(c(0, range), colors = c('white', 'dodgerblue2'))}
    
    #if(features[n] == 'H3K4me3'){range = 5; cols = colorRamp2(c(0, range), c("white", "violetred1"))
    if(features[n] == 'H3K4me3'){range = 5; cols = colorRamp2(c(0, range), c("white", "darkorchid1"))
    #cols = colorRamp2(c(0.2, range), c("white", "#3794bf"))
    }
    
    pdf(paste0(figureDir, "/positional_peaks_intensity_heatmap_", features[n], "_test_v4.pdf"),
        width = 4, height = 10) # Open a new pdf file
    Heatmap(test, 
            cluster_rows = FALSE,
            cluster_columns = FALSE, 
            show_row_names = FALSE,
            show_column_names = FALSE,
            row_split = splits, 
            cluster_column_slices = FALSE,
            column_split = factor(sample.sel, levels = c('mUA', 'mLA','mHand')),
            top_annotation = ha,
            show_heatmap_legend = TRUE,
            use_raster = TRUE,
            #raster_resize_mat = FALSE,
            raster_by_magick = FALSE,
            #col = colorRamp2(seq(0, range, length.out = breaks), rev(brewer.pal(n=breaks, name="RdBu")))
            col = cols
            #name = "mtcars", #title of legend
            #column_title = "Variables", row_title = "Samples",
            #row_names_gp = gpar(fontsize = 7) # Text size for row names
    )
    
    dev.off()
    
  }
  
  
}
