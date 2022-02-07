##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Nov 15 13:39:26 2021
##########################################################################
##########################################################################
rm(list = ls())

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0(figureDir, 'tables4plots/')

require(ggplot2)
library(VennDiagram)
library(Vennerable)
require(GenomicFeatures)
require(ChIPpeakAnno)
source('functions_chipSeq.R')
source('Functions_atac.R')

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')


########################################################
########################################################
# Figure 1: ATAC-seq positional peaks 
# 
########################################################
########################################################

##########################################
# Fig 1B peak overlapping between UA, LA and Hand excluding Head control peaks 
##########################################
version.analysis = 'Rxxxx_R10723_R11637_R12810_atac'
resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')

pval.cutoff = 6

#load(file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 'Mature.Rdata'))
load( file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 'Mature.Rdata'))

## exclude head peaks
ua = ua[!overlapsAny(ua, hc)]
la = la[!overlapsAny(la, hc)]
hd = hd[!overlapsAny(hd, hc)]

pmature = union(ua, la)
pmature = union(pmature, hd)

ol.peaks <- makeVennDiagram(list(ua, la, hd), NameOfPeaks=c('mUA', 'mLA', 'mHand'), connectedPeaks="keepAll", by = 'region', 
                            plot = TRUE)
v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

pdfname = paste0(figureDir, 'Fig1B_matureSamples_atacseqPeak_comparison_venndiagram.pdf')
pdf(pdfname, width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

dev.off()

#load(file = paste0(RdataDir, '/peaks_set_union_conditions_pval', pval.cutoff, 'version_', version.analysis, '.Rdata'))
pp = peak.merged

##########################################
# Fig 1C: feature distribution of mature peaks
##########################################
amex = makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))

peakAnnots = annotatePeak(pmature, TxDb=amex, tssRegion = c(-2000, 2000))

pdfname = paste0(figureDir, "Fig1C_matureSample_peak_featureAssignment_distribution.pdf")
pdf(pdfname, width = 8, height = 6)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

plotAnnoPie(peakAnnots)

dev.off()


### plot the partitions of Peaks into different genomics features and distance to TSS
pdfname = paste0(figureDir, "peak_to_TSS_distance_distribution.pdf")
pdf(pdfname, width = 8, height = 4)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

print(plotDistToTSS(peakAnnots))

dev.off()


pdfname = paste0(figureDir, "peak_featureAssignment_distribution_vennpie.pdf")
pdf(pdfname, width = 8, height = 4)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

vennpie(peakAnnots)

dev.off()


# load(file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 'regeneration.Rdata'))
# #ol.peaks <- makeVennDiagram(list(bld5, bld9, bld13.p, bld13.d), 
# #                            NameOfPeaks=c('UA.BL.d5', 'UA.BL.d9', 'UA.BL.d13.p', 'UA.BL.d13.d'), connectedPeaks="keepAll", by = 'region', 
# #                            plot = TRUE, fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
# 
# #v <- venn_cnt2venn(ol.peaks$vennCounts)
# #try(plot(v))
# preg = union(bld5, bld9)
# preg = union(preg, bld13.p)
# preg = union(preg, bld13.d)
# 
# 
# pdfname = paste0(figureDir, 'UA_BL_samples_peak_comparison.pdf')
# pdf(pdfname, width = 10, height = 8)
# par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
# 
# makeVennDiagram(list(bld5, bld9, bld13.p, bld13.d), 
#                 NameOfPeaks=c('UA.BL.d5', 'UA.BL.d9', 'UA.BL.d13.p', 'UA.BL.d13.d'), connectedPeaks="keepAll", by = 'region', 
#                 plot = TRUE, fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
# 
# dev.off()
# 
# 
# load(file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 
#                     'embryoStage.Rdata'))
# 
# pemb = union(es40, es44.d)
# pemb = union(pemb, es44.p)
# 
# ol.peaks <- makeVennDiagram(list(pmature, preg, pemb), NameOfPeaks=c('mature', 'reg', 'embryo'), connectedPeaks="keepAll", by = 'region', 
#                             plot = TRUE)
# 
# v <- venn_cnt2venn(ol.peaks$vennCounts)
# try(plot(v))
# 
# pdfname = paste0(figureDir, 'mature_regeneration_embryo_peak_comparison.pdf')
# pdf(pdfname, width = 10, height = 8)
# par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
# 
# v <- venn_cnt2venn(ol.peaks$vennCounts)
# try(plot(v))
# 
# dev.off()
# 

##########################################
# Fig 1D: heatmap of positional peaks 
##########################################
library(ggplot2)

load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl', version.analysis, '.Rdata'))
design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.2022.rds'))

conds = c("Mature_UA", "Mature_LA", "Mature_Hand", 'HEAD')

sample.sels = c();  cc = c()
for(n in 1:length(conds)) {
  kk = which(design$conds == conds[n]) 
  sample.sels = c(sample.sels, kk)
  cc = c(cc, rep(conds[n], length(kk)))
}

df <- data.frame(cc)
rownames(df) = colnames(keep)
colnames(df) = 'segments'
ii.gaps = c(5, 9, 12)

annot_colors = c('springgreen4', 'steelblue2', 'gold2', 'darkgray')
names(annot_colors) = c('Mature_UA', 'Mature_LA', 'Mature_Hand', 'HEAD')
annot_colors = list(segments = annot_colors)

out = pheatmap(keep, cluster_rows=TRUE, kmeans_k = NA, cutree_rows = 3,
               show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, 
         annotation_colors = annot_colors, 
         filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_rmPeaks.head.pdf'), 
         width = 6, height = 8)

hc <- out$tree_row

lbl <- cutree(hc, 3) # you'll need to change '5' to the number of gene-groups you're interested in


if(saveTable){
  yy = data.frame(keep, xx, stringsAsFactors = FALSE)
  write.csv(yy, file = paste0(tableDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head.csv'), 
            quote = FALSE, row.names = TRUE)
  
}

##########################################
# select top peaks or top promoter peaks
##########################################
load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl', version.analysis, '.Rdata'))
yy = xx

yy = yy[grep('Promoter', yy$annotation), ] # peak close to promoters
yy[grep('HOXA13|MEIS|SHOX|HOXC', yy$transcriptId),]
#yy = xx
#yy = yy[c(1:50), ]
#yy = yy[which(yy$logFC.mean>1.5), ]
# yy = yy[which(yy$max > 4 ), ] # max residual square < 0.1
keep = keep[match(rownames(yy), rownames(keep)), ]

#keep = fpm[!is.na(match(rownames(fpm), rownames(yy))), sample.sels]
gg = yy$geneId
grep('HOXA13', gg)
rownames(keep) = paste0(rownames(keep), '_', gg)
#rownames(keep) = gg
keep = as.matrix(keep)

# jj1 = grep('Embryo_Stage40', colnames(keep))
# jj2 = grep('Embryo_Stage44', colnames(keep))
# jj3 = grep('Mature_UA', colnames(keep))
# mean1 = apply(keep[, jj1], 1, mean)
# mean2 = apply(keep[, jj2], 1, mean)
# mean3 = apply(keep[, jj3], 1, mean)
# 
# kk = which(mean1< 1. & mean2 < 1. & mean3 < 1.)
# 
# # make sure max above the threshold
# nb.above.threshold = apply(keep, 1, function(x) length(which(x>3)))
# keep = keep[which(nb.above.threshold >=3), ] 
gg = rownames(keep)
gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '_'))[2])
gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '[|]'))[1])
#keep = keep[1:50, ]
rownames(keep) = gg
#kk = grep('Mature', cc)
#df <- data.frame(condition = cc[kk])
#keep = keep[,kk]
#rownames(df) = colnames(keep)
#ii.gaps = c(4, 6)
pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, 
         annotation_colors = annot_colors, 
         filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_top.promoters.pdf'), 
         width = 8, height = 6)

if(saveTable){
  write.csv(data.frame(keep, yy, stringsAsFactors = FALSE), 
            file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_top50_promoterPeaks.csv'), 
            quote = FALSE, row.names = TRUE)
  
}


##########################################
# MARA plots 
##########################################
Make.plot.motif.analysis = FALSE
if(Make.plot.motif.analysis){
  r = readRDS(file = paste0(RdataDir, '/MARA_Bayesian_ridge.rds'))
  
  zz = r$Zscore
  zz = apply(as.matrix(zz), 1, max)
  r$max.Zscore = zz
  
  # = sort(r$combined.Zscore, decreasing=TRUE)[1:50]
  sort(r$combined.Zscore, decreasing=TRUE)[1:20]
  
  sort(r$max.Zscore, decreasing=TRUE)[1:20]
  sort(r$max.Zscore, decreasing=TRUE)[1:30]
  
  top20 = sort(r$max.Zscore, decreasing = TRUE)[1:30]
  motif.names = rownames(r$Zscore)
  bb = r$Zscore[match(names(top20), motif.names), ]
  
  df <- data.frame(condition = colnames(bb))
  rownames(df) = colnames(bb)
  
  pheatmap(bb, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, 
           scale = 'none', cluster_cols=FALSE, main = paste0("Inferred z-scores (motif activity) by MARA"), 
           na_col = "white", fontsize_row = 12, annotation_col = df, 
           filename = paste0(figureDir, 'positional_peaks_MARA_ridge.pdf'), 
           width = 12, height = 10)
  
}

##########################################
# make heatmap for regeneration peaks 
##########################################
Make.Heatmap.regeneration.peaks = FALSE
if(Make.Heatmap.regeneration.peaks){
  
  conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal")
  
  sample.sels = c(); cc = c()
  for(n in 1:length(conds)) {
    kk = which(design$conds == conds[n])
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
  }
  
  res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test_v4.rds'))
  
  # select the temporal dynamic peaks
  length(which(res$prob.M0<0.05))
  length(which(res$prob.M0<0.05 & res$log2FC > 1))
  length(which(res$prob.M0<0.01 & res$log2FC > 1))
  length(which(res$prob.M0<0.01 & res$log2FC > 1.5))
  length(which(res$prob.M0<0.01 & res$log2FC > 2))
  
  jj = which(res$prob.M0 < 0.01 & res$log2FC > 1.5 )
  
  xx = res[c(jj), ]
  xx = xx[order(-xx$log2FC), ]
  #xx = xx[which(xx$min < 1), ]
  
  source('Functions_atac.R')
  keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
  keep = as.matrix(keep)
  
  kk = c(grep('Embryo_Stage40', colnames(keep)), 
         grep('Embryo_Stage44', colnames(keep)))
  kk = c(setdiff(c(1:ncol(keep)), kk), kk)
  
  keep = keep[, kk]
  df <- data.frame(cc[kk])
  rownames(df) = colnames(keep)
  colnames(df) = 'condition'
  
  ii.gaps = c(4, 8, 10, 12, 16) + 1
  
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps,
           filename = paste0(figureDir, '/heatmap_DBpeaks_fdr0.01_log2FC.1.5.pdf'), 
           width = 10, height = 16)
  
  
  ##########################################
  # highligh potential regeneration peaks, those not found in mUA, mLA, mHand and embryo stages, only in regeneration process
  # not in head control samples either
  ##########################################
  cpms = fpm[match(rownames(keep), rownames(fpm)), ]
  
  bg.cutoff = 1.0
  kk = which(apply(cpms[, grep('Embryo_Stage40', colnames(cpms))], 1, mean) < bg.cutoff & 
               apply(cpms[, grep('Embryo_Stage44_proximal', colnames(cpms))], 1, mean) < bg.cutoff &
               apply(cpms[, grep('Embryo_Stage44_distal', colnames(cpms))], 1, mean) < bg.cutoff &
               apply(cpms[, grep('Mature_UA', colnames(cpms))], 1, mean) < bg.cutoff &
               apply(cpms[, grep('Mature_LA', colnames(cpms))], 1, mean) < bg.cutoff & 
               apply(cpms[, grep('Mature_Hand', colnames(cpms))], 1, mean) < bg.cutoff 
               & cpms[, grep('HEAD', colnames(cpms))] < bg.cutoff
               )
  
  
  pheatmap(keep[kk, ], cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps,
           filename = paste0(figureDir, '/heatmap_regenerationPeaks_fdr0.01_log2FC.1.5_regeneartion.specific.pdf'), 
           width = 8, height = 10)
  
  if(saveTable){
    yy = keep[kk, ]
    yy = data.frame(yy, xx[match(rownames(yy), rownames(xx)), ], stringsAsFactors = FALSE)
    yy = yy[order(-yy$log2FC), ]                
    
    write.csv(yy, 
              file = paste0(tableDir, 'potential_regenerativepeaks.csv'), 
              quote = FALSE, row.names = TRUE)
    
  }
  
  ##########################################
  # if wannt to show something else, e.g. there is no regeneration enhancers
  ##########################################
  conds = c("Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", '')
  sample.sels = c(); cc = c()
  for(n in 1:length(conds)) {
    kk = which(design$conds == conds[n])
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
  }
  
  df <- data.frame(cc[kk])
  rownames(df) = colnames(keep)
  colnames(df) = 'condition'
  
}


##########################################
# motif ananlysis for temporal peaks 
##########################################
plot.MARA.temporal.peaks = FALSE
if(plot.MARA.temporal.peaks){
  r =  readRDS(file = paste0(RdataDir, '/MARA_Bayesian_ridge_regenerationPeaks.rds'))
  
  # exclude the zscore from mUA
  zz = r$Zscore[, -1]
  zz = apply(as.matrix(zz), 1, max)
  r$max.Zscore = zz
  
  # = sort(r$combined.Zscore, decreasing=TRUE)[1:50]
  sort(r$combined.Zscore, decreasing=TRUE)[1:20]
  
  sort(r$max.Zscore, decreasing=TRUE)[1:20]
  sort(r$max.Zscore, decreasing=TRUE)[1:30]
  sort(r$max.Zscore, decreasing=TRUE)[1:40]
  sort(r$max.Zscore, decreasing=TRUE)[1:50]
  
  topMotifs = sort(r$max.Zscore, decreasing = TRUE)[1:30]
  motif.names = rownames(r$Zscore)
  bb = r$Zscore[match(names(topMotifs), motif.names), ]
  
  df <- data.frame(condition = colnames(bb))
  rownames(df) = colnames(bb)
  
  
  pheatmap(bb, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, 
           scale = 'none', cluster_cols=FALSE, main = paste0("Inferred z-scores (motif activity) by MARA"), 
           na_col = "white", fontsize_row = 12, annotation_col = df, 
           filename = paste0(figureDir, '/MARA_bayesianRidge_temporalPeaks.pdf'), 
           width = 12, height = 10)
  
  
}
