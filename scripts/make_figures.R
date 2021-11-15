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
# Section : differentially binding peaks
# 
########################################################
########################################################




