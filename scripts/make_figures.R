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

##########################################
# compare ua, la, hd, hc
##########################################
#ua = ua[!overlapsAny(ua, hc)]
#la = la[!overlapsAny(la, hc)]
#hd = hd[!overlapsAny(hd, hc)]

annot_colors = c('springgreen4', 'steelblue2', 'gold2', 'darkgray')

pmature = union(ua, la)
pmature = union(pmature, hd)
pmature = union(pmature, hc)

library(Vennerable)
#ol.peaks <- makeVennDiagram(list(ua, la, hd, hc), NameOfPeaks=c('mUA', 'mLA', 'mHand', 'Head'), connectedPeaks="merge", by = 'region', 
#                            plot = TRUE, fill=annot_colors)
#v <- venn_cnt2venn(ol.peaks$vennCounts)
#try(plot(v))

pdfname = paste0(figureDir, 'FigS1B_matureSamples_Head_atacseqPeak_comparison_venndiagram.pdf')
pdf(pdfname, width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

#v <- venn_cnt2venn(ol.peaks$vennCounts)
#try(plot(v))
ol.peaks <- makeVennDiagram(list(ua, la, hd, hc), NameOfPeaks=c('mUA', 'mLA', 'mHand', 'Head'), connectedPeaks="keepAll", by = 'region', 
                            plot = TRUE, fill=annot_colors)

dev.off()

#load(file = paste0(RdataDir, '/peaks_set_union_conditions_pval', pval.cutoff, 'version_', version.analysis, '.Rdata'))
# pp = peak.merged


##########################################
# Fig 1C: feature distribution of mature peaks
##########################################
amex = makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))

peakAnnots = annotatePeak(pmature, TxDb=amex, tssRegion = c(-2000, 2000))

library(ggplot2)

pdfname = paste0(figureDir, "FigS1C_matureSample.Head_peak_featureAssignment_distribution.pdf")
pdf(pdfname, width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

plotPeakAnnot_piechart(peakAnnots)

dev.off()


##########################################
# compare ua, la, hd excluding hc
##########################################
## exclude head peaks
ua = ua[!overlapsAny(ua, hc)]
la = la[!overlapsAny(la, hc)]
hd = hd[!overlapsAny(hd, hc)]

pmature = union(ua, la)
pmature = union(pmature, hd)


pdfname = paste0(figureDir, 'FigS1C_matureSamples_atacseqPeak_comparison_venndiagram.pdf')
pdf(pdfname, width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

ol.peaks <- makeVennDiagram(list(ua, la, hd), NameOfPeaks=c('mUA', 'mLA', 'mHand'), connectedPeaks="keepAll", by = 'region', 
                            plot = TRUE, fill=annot_colors[1:3])
#v <- venn_cnt2venn(ol.peaks$vennCounts)
#v <- venn_cnt2venn(ol.peaks$vennCounts)
#try(plot(v))
dev.off()

#load(file = paste0(RdataDir, '/peaks_set_union_conditions_pval', pval.cutoff, 'version_', version.analysis, '.Rdata'))
pp = peak.merged

##########################################
# Fig 1C: feature distribution of mature peaks
##########################################
peakAnnots = annotatePeak(pmature, TxDb=amex, tssRegion = c(-2000, 2000))

pdfname = paste0(figureDir, "Fig1C_matureSample.noHead_peak_featureAssignment_distribution.pdf")
pdf(pdfname, width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

plotPeakAnnot_piechart(peakAnnots = peakAnnots)

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


##########################################
# Fig 1D: heatmap of positional peaks 
##########################################
library(dendextend)
library(ggplot2)

load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl', version.analysis, '.Rdata'))
yy = data.frame(apply(as.matrix(keep[,grep(conds[1], colnames(keep))]), 1, mean), 
                apply(as.matrix(keep[,grep(conds[2], colnames(keep))]), 1, mean),
                apply(as.matrix(keep[,grep(conds[3], colnames(keep))]), 1, mean)
)
colnames(yy) = conds[c(1:3)]  

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

yy <- t(apply(yy, 1, cal_z_score))

my_hclust_gene <- hclust(dist(yy), method = "complete")

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 4)
my_gene_col <- data.frame(cluster = my_gene_col)

xx = data.frame(xx, my_gene_col, stringsAsFactors = FALSE)

save(xx, keep, file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl_clustered.4groups_', 
                             version.analysis, '.Rdata'))


conds = c("Mature_UA", "Mature_LA", "Mature_Hand")

df <- data.frame(conds)
rownames(df) = conds
colnames(df) = 'segments'

annot_colors = c('springgreen4', 'steelblue2', 'gold2')
names(annot_colors) = c('Mature_UA', 'Mature_LA', 'Mature_Hand')
annot_colors = list(segments = annot_colors)

pheatmap(yy, annotation_row = my_gene_col, annotation_col = df, show_rownames = FALSE, scale = 'none', 
         show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = FALSE,  clustering_method = 'complete',
         annotation_colors = annot_colors, cutree_rows = 4, 
         filename = paste0(figureDir, '/Fig1D_heatmap_positionalPeaks_fdr0.01_log2FC.1_rmPeaks.head_1000peaks.pdf'), 
         width = 8, height = 10)


if(saveTable){
  yy = data.frame(keep, xx, stringsAsFactors = FALSE)
  write.csv(yy, file = paste0(tableDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head.csv'), 
            quote = FALSE, row.names = TRUE)
  
  # save bed files of different peak groups
  load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl_clustered.4groups_', 
                     version.analysis, '.Rdata'))
  
  for(c in unique(xx$cluster))
  {
    # c = '1'
    ss = xx[which(xx$cluster == c), ]
    
    ss = data.frame(ss$seqnames, ss$start, ss$end, rownames(ss), rep(0, nrow(ss)), rep('*', nrow(ss)))
    
    write.table(ss, file = 
      paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/makeHeatmaps/peak_groups/', 
             'peak_group_', c, '.bed'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
  }
  
}

##########################################
# feature distribution of groups of positional peaks 
##########################################
amex = makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))

pp = data.frame(t(sapply(rownames(xx), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
pp$strand = '*'

pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")

pp.annots = annotatePeak(pp[which(xx$cluster == 2)], TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

pdfname = paste0(figureDir, "Fig1E_positional_peak_feature_distribution_group2.pdf")
pdf(pdfname, width = 8, height = 6)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

plotAnnoPie(pp.annots)

dev.off()

pp.annots = annotatePeak(pp[which(xx$cluster == 3)], TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

pdfname = paste0(figureDir, "Fig1E_positional_peak_feature_distribution_group3.pdf")
pdf(pdfname, width = 8, height = 6)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

plotAnnoPie(pp.annots)

dev.off()


pp.annots = annotatePeak(pp[which(xx$cluster == 1| xx$cluster == 4)], TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

pdfname = paste0(figureDir, "Fig1E_positional_peak_feature_distribution_group1_4.pdf")
pdf(pdfname, width = 8, height = 6)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

plotAnnoPie(pp.annots)

dev.off()

##########################################
# select top peaks or top promoter peaks
##########################################
load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl_clustered.4groups_', 
                   version.analysis, '.Rdata'))
conds = c("Mature_UA", "Mature_LA", "Mature_Hand")

yy = xx

yy = yy[grep('Promoter', yy$annotation), ] # peak close to promoters
yy[grep('HOXA13|MEIS|SHOX|HOXC', yy$transcriptId),]

keep = keep[match(rownames(yy), rownames(keep)), ]


keep = data.frame(apply(as.matrix(keep[,grep(conds[1], colnames(keep))]), 1, mean), 
                apply(as.matrix(keep[,grep(conds[2], colnames(keep))]), 1, mean),
                apply(as.matrix(keep[,grep(conds[3], colnames(keep))]), 1, mean)
)
colnames(keep) = conds

gg = yy$geneId
grep('HOXA13', gg)
rownames(keep) = paste0(rownames(keep), '_', gg)
#rownames(keep) = gg
keep = as.matrix(keep)
keep = keep[c(1:20),]
xx = xx[c(1:20), ]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

keep <- t(apply(keep, 1, cal_z_score))


#keep = fpm[!is.na(match(rownames(fpm), rownames(yy))), sample.sels]


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
pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'none', show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col = df, 
         annotation_colors = annot_colors,
         filename = paste0(figureDir, '/Figure1F_heatmap_positionalPeaks_fdr0.01_log2FC.1_top.promoters.pdf'), 
         width = 8, height = 6)

if(saveTable){
  write.csv(data.frame(keep,  stringsAsFactors = FALSE), 
            file = paste0(tableDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_top20_promoterPeaks.csv'), 
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
