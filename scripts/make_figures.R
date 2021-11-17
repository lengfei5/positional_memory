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

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

########################################################
########################################################
# Section : ATAC-seq peak summary
# 
########################################################
########################################################
RdataDir = paste0('../results/R10723_Rxxxx_R11637_atacseq_R11876_CutTag/Rdata')
version.analysis = 'R10723_Rxxxx_R11637_atacseq_R11876_CutTag'
pval.cutoff = 4

load(file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 'Mature.Rdata'))
pmature = union(ua, la)
pmature = union(pmature, hd)

ol.peaks <- makeVennDiagram(list(ua, la, hd), NameOfPeaks=c('mUA', 'mLA', 'mHand'), connectedPeaks="keepAll", by = 'region', 
                            plot = TRUE)

v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))


pdfname = paste0(figureDir, 'mature_samples_peak_comparison.pdf')
pdf(pdfname, width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

dev.off()


load(file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 'regeneration.Rdata'))
#ol.peaks <- makeVennDiagram(list(bld5, bld9, bld13.p, bld13.d), 
#                            NameOfPeaks=c('UA.BL.d5', 'UA.BL.d9', 'UA.BL.d13.p', 'UA.BL.d13.d'), connectedPeaks="keepAll", by = 'region', 
#                            plot = TRUE, fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))

#v <- venn_cnt2venn(ol.peaks$vennCounts)
#try(plot(v))
preg = union(bld5, bld9)
preg = union(preg, bld13.p)
preg = union(preg, bld13.d)


pdfname = paste0(figureDir, 'UA_BL_samples_peak_comparison.pdf')
pdf(pdfname, width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

makeVennDiagram(list(bld5, bld9, bld13.p, bld13.d), 
                NameOfPeaks=c('UA.BL.d5', 'UA.BL.d9', 'UA.BL.d13.p', 'UA.BL.d13.d'), connectedPeaks="keepAll", by = 'region', 
                plot = TRUE, fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))

dev.off()


load(file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 
                    'embryoStage.Rdata'))

pemb = union(es40, es44.d)
pemb = union(pemb, es44.p)

ol.peaks <- makeVennDiagram(list(pmature, preg, pemb), NameOfPeaks=c('mature', 'reg', 'embryo'), connectedPeaks="keepAll", by = 'region', 
                            plot = TRUE)

v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

pdfname = paste0(figureDir, 'mature_regeneration_embryo_peak_comparison.pdf')
pdf(pdfname, width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

dev.off()



load(file = paste0(RdataDir, '/peaks_set_union_conditions_pval', pval.cutoff, 'version_', version.analysis, '.Rdata'))
pp = peak.merged

amex = makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))

peakAnnots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000))


### plot the partitions of Peaks into different genomics features and distance to TSS
pdfname = paste0(figureDir, "peak_to_TSS_distance_distribution.pdf")
pdf(pdfname, width = 8, height = 4)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

print(plotDistToTSS(peakAnnots))

dev.off()

pdfname = paste0(figureDir, "peak_featureAssignment_distribution.pdf")
pdf(pdfname, width = 8, height = 6)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

plotAnnoPie(peakAnnots)

dev.off()

pdfname = paste0(figureDir, "peak_featureAssignment_distribution_vennpie.pdf")
pdf(pdfname, width = 8, height = 4)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

vennpie(peakAnnots)

dev.off()

########################################################
########################################################
# Section : positional peaks from mature samples
# 
########################################################
########################################################
RdataDir = "../results/atac_rna_chipseq_analysis_20211007/Rdata"

load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical_v1.Rdata'))
fpm = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))

# prepare the background distribution
fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
rownames(fpm) = gsub('_', '-', rownames(fpm))

hist(fpm.bg, breaks = 100, main = 'background distribution')
abline(v = 1, col = 'red', lwd = 2.0)
quantile(fpm.bg, c(0.95, 0.99))

##########################################
## make Granges and annotate peaks
##########################################
Make.Granges.and.peakAnnotation = TRUE
if(Make.Granges.and.peakAnnotation){
  require(ChIPpeakAnno)
  require(ChIPseeker)
  
  pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  # annotation from ucsc browser ambMex60DD_genes_putative
  amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
  pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  #plotAnnoBar(pp.annots)
  
  pp.annots = as.data.frame(pp.annots)
  rownames(pp.annots) = rownames(fpm)
  
  promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
  
}

##########################################
# make heatmap
##########################################
Make.Heatmap.positional.peaks = FALSE
if(Make.Heatmap.positional.peaks){
  conds = c("Mature_UA", "Mature_LA", "Mature_Hand")
  
  sample.sels = c();  cc = c()
  for(n in 1:length(conds)) {
    #kk = which(design$conds == conds[n] & design$SampleID != '136159')
    kk = which(design$conds == conds[n]) 
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
  }
  
  res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_v6.rds'))
  
  # select the positional peaks with 
  fdr.cutoff = 0.01; logfc.cutoff = 1
  jj = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & res$logFC.mLA.vs.mUA > logfc.cutoff) |
               (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & res$logFC.mHand.vs.mUA > logfc.cutoff)|
               (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & res$logFC.mHand.vs.mLA > logfc.cutoff)
  )
  
  jj1 = which(res$prob.M0<0.01 & res$log2FC>1)
  jj2 = which(res$pval.lrt < 0.001 & res$log2FC > 1)
  
  xx = res[c(jj), ]
  #xx = xx[order(-xx$log2FC.mature), ]
  xx[grep('HOXA13', xx$transcriptId), ]
  
  # filter the peaks from head control sample
  Filtering.peaks.in.Head.samples = TRUE
  if(Filtering.peaks.in.Head.samples){
    p0 = pp[match(rownames(xx), names(pp))]
    
    ctl = fpm[, grep('HEAD', colnames(fpm))]
    ctl = ctl[which(ctl>2.5)]
    p.ctl = pp[match(names(ctl), names(pp))]
    
    non.overlap = !overlapsAny(p0, p.ctl)
    
    xx = xx[non.overlap, ]
    
  }
  
  xx[grep('HOXA13', xx$transcriptId), ]
  
  # sort positional peaks with logFC
  #xx = xx[order(-xx$logFC.mean), ]
  xx = xx[order(-xx$log2FC), ]
  
  keep = fpm[match(rownames(xx), rownames(fpm)), sample.sels]
  keep = as.matrix(keep)
  #keep = keep[match(rownames(xx), rownames(xx)), ]
  
  library(ggplot2)
  df <- data.frame(cc)
  rownames(df) = colnames(keep)
  
  ii.gaps = c(5, 8)
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, 
           filename = paste0(figureDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_rmPeaks.head.pdf'), 
           width = 8, height = 12)
  
  if(saveTable){
    yy = data.frame(keep, xx, stringsAsFactors = FALSE)
    write.csv(yy, file = paste0(tableDir, 'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head.csv'), 
              quote = FALSE, row.names = TRUE)
    
  }
  
  rm.UA.oldBatches = FALSE
  if(rm.UA.oldBatches){
    keep = fpm[match(rownames(xx), rownames(fpm)), sample.sels]
    keep = as.matrix(keep)
    
    jj = grep('102655|74938', colnames(keep))
    keep = keep[, -jj]
    df <- data.frame(cc[-jj])
    rownames(df) = colnames(keep)
    colnames(df) = 'segments'
    
    ii.gaps = c(5, 8) -2
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, 
             filename = paste0(figureDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_rmPeaks.head_rm.oldUA.pdf'), 
             width = 8, height = 12)
    
  }
  
  plot.top.peaks.within.promoters = FALSE
  if(plot.top.peaks.within.promoters){
    
    yy = xx
    
    yy = yy[grep('Promoter', yy$annotation), ] # peak close to promoters
    yy[grep('HOXA13', yy$transcriptId),]
        
    yy = yy[which(yy$max > 3 & yy$min < 2), ]
    
    keep = fpm[match(rownames(yy), rownames(fpm)), sample.sels]
    
    
    gg = res$geneId[match(rownames(keep), rownames(res))]
    grep('HOXA13', gg)
    rownames(keep) = paste0(rownames(keep), '_', gg)
    #rownames(keep) = gg
    keep = as.matrix(keep)
    
    gg = rownames(keep)
    gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '_'))[2])
    gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '[|]'))[1])
    #keep = keep[1:50, ]
    #rownames(keep) = gg
    #kk = grep('Mature', cc)
    #df <- data.frame(condition = cc[kk])
    #keep = keep[,kk]
    #rownames(df) = colnames(keep)
    jj = grep('102655|74938', colnames(keep))
    keep = keep[, -jj]
    df <- data.frame(cc[-jj])
    rownames(df) = colnames(keep)
    colnames(df) = 'segments'
    
    ii.gaps = c(5, 8) -2
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, fontsize_row = 11, gaps_col = ii.gaps,
             filename = paste0(figureDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_top50_promoter.pdf'), 
             width = 16, height = 12)
    
    
    if(saveTable){
      write.csv(data.frame(keep, yy, stringsAsFactors = FALSE), 
                file = paste0(tableDir, 'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_top50_promoterPeaks.csv'), 
                quote = FALSE, row.names = TRUE)
      
    }
       
  }
  
}

##########################################
# MARA plots 
##########################################
Make.plot.motif.analysis = FALSE
if(Make.plot.motif.analysis){
  r = readRDS(file = paste0(resDir, '/MARA_Bayesian_ridge.rds'))
  
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
  
  pheatmap(bb, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, 
           scale = 'none', cluster_cols=FALSE, main = paste0("Inferred z-scores (motif activity) by MARA"), 
           na_col = "white", fontsize_row = 12, 
           filename = paste0(figureDir, 'positional_peaks_MARA_ridge.pdf'), 
           width = 8, height = 10) 
  
}



