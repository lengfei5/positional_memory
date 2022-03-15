##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: analyze the histone markers to test the chromatin states invovled in postional memory and regeneration
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Dec 14 10:04:12 2021
##########################################################################
##########################################################################
rm(list=ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('functions_chipSeq.R')
source('Functions_atac.R')

version.analysis = 'CT_analysis_20220311'
#peakDir = "Peaks/macs2_broad"
saveTable = TRUE

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

# figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
# tableDir = paste0(figureDir, 'tables4plots/')

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/'
gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')

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

########################################################
########################################################
# Section I : sequencing quality controls
# fragment size distribution 
# sequence saturation analysis
########################################################
########################################################
Process.design.stats = FALSE
if(Process.design.stats){
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R12965_CT/'
  design = read.table(file = paste0(dataDir, 'sampleInfos_v2_parsed.txt'), header = TRUE)
  
  design$sample = sapply(design$fileName, function(x) unlist(strsplit(as.character(x), '_'))[1])
  design$marks = sapply(design$fileName, function(x) unlist(strsplit(as.character(x), '_'))[2])
  
  design = design[, -2]
  
  #design = rbind(design, read.table(paste0(dataDir, 'R12810_cuttag/sampleInfos_conditions_parsed.txt'), header = TRUE))
  stats = read.table(file = paste0(dataDir, 'countStatTable.txt'), header = TRUE)
  stats$ids = sapply(stats$Name, function(x) unlist(strsplit(as.character(x), '_'))[1])
  stats = stats[match(design$sampleID, stats$ids), ]
  
  design = data.frame(design, stats, stringsAsFactors = FALSE)
  
  #stats = rbind(stats, read.table(file = paste0(dataDir, 'R12810_cuttag/nf_out/result/countStatTable.txt'), header = TRUE))
  
  colnames(design)[c(4, 6)] = c('fileName', 'trimmed')
  design = design[, -ncol(design)]
  
  for(n in 5:ncol(design))
  {
    design[,n] = as.numeric(design[,n])/10^6
  }
  
  design$percent.align = design$unique/ design$rmChrM
  design$percent.dup = 1 - design$unique.rmdup/design$unique
  design$percent.usable = design$unique.rmdup/design$Total
  
  design = design[order(design$marks, design$sample), ]
  
  write.csv(design, file = paste0(resDir, '/R12965_CT_QC_stats.csv'), row.names = FALSE, quote = FALSE)
  
  
  
  index = c()
  for(n in 1:nrow(design))
  {
    # n = 1;
    cat(n, '\n')
    #cc.files = cnts[grep(design$sampleID[n], cnts)]
    ii = grep(design$sampleID[n], stats$fileName)
    if(length(ii) == 1){
      index = c(index, ii)
    }else{
      index = c(index, NA)
      cat(length(ii), 'fileName Found for ', design$sampleID[n], '\n')
    }
    
  }
  
  stats = data.frame(design, stats[index, ], stringsAsFactors = FALSE)
  
  colnames(stats) = c('sampleID', 'samples', 'fileName', 'total',  'adapter.trimmed', 'mapped', 'chrM.rm', 'unique', 'unique.rmdup')
  
  stats$trimming.pct = as.numeric(stats$adapter.trimmed)/as.numeric(stats$total)
  stats$mapped.pct =  as.numeric(stats$mapped)/as.numeric(stats$adapter.trimmed)
  stats$mito.pct = 1.0 - as.numeric(stats$chrM.rm)/as.numeric(stats$adapter.trimmed)
  stats$multimapper.pct = 1- as.numeric(stats$unique) / as.numeric(stats$mapped)
  stats$dup.rate = 1.0 - as.numeric(stats$unique.rmdup)/as.numeric(stats$unique)
  stats$pct.usable = stats$unique.rmdup / stats$total
  
  stats$usable = stats$unique.rmdup/10^6
  
  sels = match(c('163648', '163651', '163656', '163659', 
                 '182553', '182554', '182555', '182556'), design$sampleID)
  
  write.csv(stats, file = paste0(resDir, '/', version.analysis,  '_QCs_stats.csv'), row.names = FALSE)
  
}else{
  
  design = read.csv(file = paste0(resDir, '/R11876_R12810_CT_QCs_stats.csv'), header = TRUE)
  design = design[grep('AL1', design$samples, invert = TRUE), ]
  # 163648, 163651, 163656, 163659 technical replicates of 182553-182556
  design = design[grep('163648|163651|163656|163659', design$sampleID, invert = TRUE), ] 
  
  xx = read.csv(file = paste0(resDir, '/R12965_CT_QC_stats.csv'), header = TRUE)
  
  design$sample = sapply(design$samples, function(x) unlist(strsplit(as.character(x), '_'))[1])
  design$marks = sapply(design$samples, function(x) unlist(strsplit(as.character(x), '_'))[2])
  design = design[,c(1, 17, 18, 3:16)]
  
  design = design[, c(1:10, 12, 15:16)]
  colnames(design) = colnames(xx)
  
  design = rbind(design, xx)
  design$sample[which(design$sample == 'Hand')] = 'mHand'
  design$sample[which(design$sample == 'LA')] = 'mLA'
  design$sample[which(design$sample == 'UA')] = 'mUA'
  design$sample[which(design$sample == 'BL13days.d')] = 'BL13days.dist'
  design$sample[which(design$sample == 'BL13days.p')] = 'BL13days.prox'
  
  design$marks[which(design$marks == 'K27ac')] = 'H3K27ac'
  design$marks[which(design$marks == 'K27me3et')] = 'H3K27me3'
  design$marks[which(design$marks == 'K4me1')] = 'H3K4me1'
  design$marks[which(design$marks == 'K4me3')] = 'H3K4me3'
  
  design$condition = paste0(design$marks, '_', design$sample)
  
  saveRDS(design, file = paste0(RdataDir, '/design_matrix_R11876_R12810_R12965_', version.analysis, '.rds'))
  
  write.csv(design, file = paste0(resDir, '/R11876_R12810_R12965_', version.analysis,  '_QCs_stats.csv'), row.names = FALSE)
  
}


##########################################
# saturation analysis 
##########################################
source('Functions_utility.R')
#Sequence.Saturation.Analysis(design)

########################################################
########################################################
# Section II: peak processing and find consensus peaks
# narrow peaks were called with macs2 for H3K4me3 and H3K4me1
# H3K27ac is too noise to call peaks 
# H3K27me3 have broad domains, which are not the main focus for now
# So the peak consensus will be 
# - H3K4me3/1 narrow peaks and ATAC-seq peaks (to calculate the scaling factors for tracks)
# - limb expressed gene promoters and non-expressed gene promoters (negative controls)
########################################################
########################################################
source('functions_chipSeq.R')
design = read.csv(file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/",
                                "histMod_CT_using/R11876_R12810_R12965_CT_analysis_20220217_QCs_AK.csv"))

design = design[!is.na(design$sampleID), ]
design$fileName = paste0(design$condition, '_', design$sampleID)

##########################################
# process called peaks for H3K4me3/1
##########################################
peakDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/calledPeaks/macs2'
peak.files = list.files(path = peakDir,
                        pattern = '*_peaks.xls', full.names = TRUE)

design.sel = design[which(design$marks == 'H3K4me1'| design$marks == 'H3K4me3'), ]

# select peak files
index  = c()
for(n in 1:nrow(design.sel))
{
  test = grep(design.sel$sampleID[n], peak.files)
  if(length(test) != 1) {
    cat(length(test), 'peak files Found \n')
  }else{
    index = c(index, test)
  }
}

peak.files = peak.files[index]

##########################################
# manually check peak overlapping between replicates and also define consensus peaks
##########################################
Manually.identify.peak.consensus = FALSE
pval.cutoff = 6 # the 10^-6 seems to be better when checking peak overlapping in replicates

if(Manually.identify.peak.consensus){
  peaks = c()
  
  for(n in 1:length(peak.files)) 
  {
    cat(n, '\n')
    p = readPeakFile(peak.files[n], as = "GRanges");
    #eval(parse(text = paste0("p = pp.", k)));
    with.p.values = "X.log10.pvalue." %in% colnames(mcols(p))
    if(with.p.values) {
      p <- p[mcols(p)[,"X.log10.pvalue."] > pval.cutoff];
      p = GenomicRanges::reduce(p);
      #peaks10= c(peaks10, p10);
    }else{ 
      cat("no p values conlumn found for -- ", design.matrix$file.name[k], "\n");
      PLOT.p10 = FALSE;
    }
    #p = reduce(p)
    peaks= c(peaks, p)
  }
  
  names(peaks) = design.sel$fileName
  
  saveRDS(peaks, file = paste0(RdataDir, '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_pval.', pval.cutoff, '.rds'))
  save(peaks, design, design.sel, 
       file = paste0(RdataDir, '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_pval.', pval.cutoff, '.Rdata'))
  
  ##########################################
  #  # try to merge BL time series
  ##########################################
  peaks = readRDS(file = paste0(RdataDir, '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_pval.', pval.cutoff, '.rds'))
  design.sel$condition = droplevels(design.sel$condition)
  
  source('functions_chipSeq.R')
  
  conds = as.character(unique(design.sel$condition))
  
  ##########################################
  # H3K4me3
  ##########################################
  ## H3Kme3 UA
  n = 1
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  
  ol.peaks <- makeVennDiagram(peaks[kk[c(2, 3,4)]], NameOfPeaks=names(peaks)[kk[c(2, 3, 4)]], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, '2_p', pval.cutoff, '.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  xx1 = intersect(peaks[[kk[1]]], peaks[[kk[2]]])
  xx2 = intersect(peaks[[kk[1]]], peaks[[kk[3]]])
  xx3 = intersect(peaks[[kk[1]]], peaks[[kk[4]]])
  xx4 = intersect(peaks[[kk[2]]], peaks[[kk[3]]])
  xx5 = intersect(peaks[[kk[2]]], peaks[[kk[4]]])
  xx6 = intersect(peaks[[kk[3]]], peaks[[kk[4]]])
  
  xx = union(xx1, union(xx2, union(xx3, union(xx4, union(xx5, xx6)))))
  
  peaks_H3K4me3_mUA = GenomicRanges::reduce(xx)
  
  peaks.merged = xx
  
  
  ## H3K4me3_mLA
  n = 3
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  
  ol.peaks <- makeVennDiagram(peaks[kk[c(2, 3,4)]], NameOfPeaks=names(peaks)[kk[c(2, 3, 4)]], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, '2_p', pval.cutoff, '.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  peaks_H3K4me3_mLA = peaks[[kk[2]]]
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, peaks_H3K4me3_mLA))
  
  ## H3K4me3_mHand
  n = 4
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk[c(2, 3,4)]], NameOfPeaks=names(peaks)[kk[c(2, 3, 4)]], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, '2_p', pval.cutoff, '.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  xx = peaks[[kk[2]]]
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  
  
  # H3K4me3_BL5days
  n = 8
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  xx = peaks[[kk[1]]]
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  cat(length(xx), ' new peaks \n')
  cat(length(peaks.merged), ' merged peaks \n')
  
  # H3K4me3_BL9days
  n = 10
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  
  xx = intersect(peaks[[kk[1]]], peaks[[kk[2]]])
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  cat(length(xx), ' new peaks \n')
  cat(length(peaks.merged), ' merged peaks \n')
  
  
  # H3K4me3_BL13days.prox
  n = 12
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  
  xx = intersect(peaks[[kk[1]]], peaks[[kk[2]]])
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  cat(length(xx), ' new peaks \n')
  cat(length(peaks.merged), ' merged peaks \n')
  
  # H3K4me3_BL13days.dist
  n = 14
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  xx = intersect(peaks[[kk[1]]], peaks[[kk[2]]])
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  cat(length(xx), ' new peaks \n')
  cat(length(peaks.merged), ' merged peaks \n')
  
  
  peaks.merged.H3K4me3 =  peaks.merged
  saveRDS(peaks.merged.H3K4me3, file = paste0(RdataDir,  
                                              '/histMarkers_macs2peaks_H3K4me3_consensusPeaks_pval.', pval.cutoff, '.rds'))
  
  ##########################################
  # H3K4me1
  ##########################################
  ## peaks_H3K4me1_mLA
  n = 2 
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  
  ol.peaks <- makeVennDiagram(peaks[kk[c(2, 3,4)]], NameOfPeaks=names(peaks)[kk[c(2, 3, 4)]], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, '2_p', pval.cutoff, '.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  peaks_H3K4me1_mLA = peaks[[kk[2]]]
  
  # peaks.merged = peaks_H3K4me1_mLA
  peaks.merged = union(peaks.merged, peaks_H3K4me1_mLA)
  peaks.merged = GenomicRanges::reduce(peaks.merged)
  
  # H3K4me1_mUA
  n = 5
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  
  ol.peaks <- makeVennDiagram(peaks[kk[c(2, 3,4)]], NameOfPeaks=names(peaks)[kk[c(2, 3, 4)]], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, '2_p', pval.cutoff, '.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  xx1 = intersect(peaks[[kk[2]]], peaks[[kk[3]]])
  xx2 = intersect(peaks[[kk[2]]], peaks[[kk[4]]])
  xx3 = intersect(peaks[[kk[3]]], peaks[[kk[4]]])
  
  xx = union(xx1, union(xx2, xx3))
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  
  # H3K4me1_mHand
  n = 6
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  xx = peaks[[kk[2]]]
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  
  # H3K4me1_BL5days
  n = 7
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  xx = peaks[[kk[1]]]
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  cat(length(xx), ' new peaks \n')
  cat(length(peaks.merged), ' merged peaks \n')
  
  # H3K4me1_BL9days
  n = 9
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  
  xx = intersect(peaks[[kk[1]]], peaks[[kk[2]]])
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  cat(length(xx), ' new peaks \n')
  cat(length(peaks.merged), ' merged peaks \n')
  
 
  # H3K4me1_BL13days.prox
  n = 11
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  
  xx = intersect(peaks[[kk[1]]], peaks[[kk[2]]])
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  cat(length(xx), ' new peaks \n')
  cat(length(peaks.merged), ' merged peaks \n')
  
  
  # H3K4me1_BL13days.dist
  n = 13
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  xx = intersect(peaks[[kk[1]]], peaks[[kk[2]]])
  eval(parse(text = paste0('peaks_', cc, ' = xx')))
  peaks.merged = GenomicRanges::reduce(union(peaks.merged, xx))
  cat(length(xx), ' new peaks \n')
  cat(length(peaks.merged), ' merged peaks \n')
  
  peaks.merged.H3K4me1 =  peaks.merged
  saveRDS(peaks.merged.H3K4me1, file = paste0(RdataDir,  
                                              '/histMarkers_macs2peaks_H3K4me1_consensusPeaks_pval.', pval.cutoff, '.rds'))
  
  peaks.all = union(peaks.merged.H3K4me3, peaks.merged.H3K4me1)
  
  save(peaks.all, file = paste0(RdataDir,  
                                   '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_pval.', pval.cutoff, '.Rdata'))
  
  sum(overlapsAny(peaks.merged.H3K4me3, peaks.merged.H3K4me1))
 
  
  ol.peaks <- makeVennDiagram(list(peaks.merged.H3K4me1, peaks.merged.H3K4me3), 
                              NameOfPeaks=c('peaks.merged.H3K4me1', 'peaks.merged.H3K4me3'), connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_H3K4me3_vs_H3K4me1_', 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
}

##########################################
# prepare the consensus peaks for read quantification 
# histoneMarker peaks
# atac-seq peaks
# tss regions
##########################################
## compare histM peaks and atacseq peaks
pval.cutoff = 6
#load(file = paste0(RdataDir, '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_pval.', pval.cutoff, '.Rdata'))
peaks.all = readRDS(file = paste0(RdataDir, 
                                  '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_pval.', pval.cutoff, '.rds'))

# export(object = peaks.all,  con = paste0(resDir, "/histM_peaks_all_440k.bed"), format = 'bed')
pp = readRDS(file = "../results/Rxxxx_R10723_R11637_R12810_atac/Rdata/ATACseq_peak_consensus.rds")

sum(overlapsAny(pp, peaks.all))

ol.peaks <- makeVennDiagram(list(peaks.all, pp), 
                            NameOfPeaks=c('peaks.histM', 'peak.atac'), connectedPeaks="keepAll", main='histM_vs_atac_peaks')
v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

pdf(paste0(resDir, '/manualCheck_peakOverlapping_H3K4me1.3_vs_atacseq_', 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
try(plot(v))
dev.off()


peaks.all = union(peaks.all, pp)
#peaks.all = GenomicRanges::reduce(peaks.all)
saveRDS(peaks.all, file = paste0(RdataDir, 
                                 '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_and_ATACseq_consensusPeaks.rds'))

export(object = peaks.all,  con = paste0(resDir, "/histM_peaks_plus_atacseqPeaks_450k.bed"), format = 'bed')

##########################################
# compare the histM peaks considered with all annotated TSS
##########################################
peaks.all = readRDS(file = paste0(RdataDir, 
                                  '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_and_ATACseq_consensusPeaks.rds'))

annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))

tss_all = data.frame(annot$chr_transcript, annot$start_transcript, (as.numeric(annot$start_transcript) + 1), 
                     annot$strand_transcript, annot$transcriptID, annot$geneID)
colnames(tss_all) = c('X1', 'X2', 'X3', 'strand', 'transcriptID', 'geneID')
tss_all$strand = '*'
pp = makeGRangesFromDataFrame(tss_all, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand", keep.extra.columns=TRUE)
xx = resize(pp, width = 1000, fix = 'center', use.names = TRUE)

pp = xx

saveRDS(pp, file = paste0(RdataDir, '/all_181985_TSS_transcriptID_geneID_with.1kb.width.rds'))

#pp = GenomicRanges::reduce(pp)

sum(overlapsAny(xx, peaks.all,  ignore.strand=TRUE))

tss_missed = pp[!overlapsAny(pp, peaks.all, ignore.strand=TRUE)]

saveRDS(tss_missed, file = paste0(RdataDir, '/missedTSS_140k_inHistM_atacseq_peaks_transcriptID_geneID_with.1kb.width.rds'))

xx = GenomicRanges::reduce(tss_missed) 
saveRDS(xx, file = paste0(RdataDir, '/missedTSS_dupReduced_90k_inHistM_atacseq_peaks_transcriptID_geneID_with.1kb.width.rds'))

peaks_histM_atac_tss = union(peaks.all, tss_missed)

saveRDS(peaks_histM_atac_tss, file = paste0(RdataDir, '/Peaks_histMarkers_ATACseq_missedTSS.dupReduced.rds'))

export(object = peaks_histM_atac_tss,  con = paste0(resDir, "/Peaks_histMarkers_ATACseq_missedTSS.dupReduced.bed"), format = 'bed')

## too slow for some reasons
# ol.peaks <- makeVennDiagram(list(peaks.all, pp), 
#                             NameOfPeaks=c('peaks.histM', 'all.annotated.tss'), connectedPeaks="keepAll", 
#                             main='histM_vs_allTSS')
# 
# v <- venn_cnt2venn(ol.peaks$vennCounts)
# 
# try(plot(v))
# 
# pdf(paste0(resDir, '/manualCheck_peakOverlapping_H3K4me1.3_vs_atacseq_', 'p', pval.cutoff, '.pdf'),  height = 10, width = 10)
# try(plot(v))
# dev.off()
# 

##########################################
#  make saf files
##########################################
peaks = readRDS(file = paste0(RdataDir, '/Peaks_histMarkers_ATACseq_missedTSS.dupReduced.rds'))
tss = readRDS(file = paste0(RdataDir, '/missedTSS_dupReduced_90k_inHistM_atacseq_peaks_transcriptID_geneID_with.1kb.width.rds'))

kk = which(overlapsAny(peaks, tss) == TRUE)
# clean peaks
peaks = data.frame(peaks)
colnames(peaks)[c(1:3)] = c("chr", "start", "end")
dim(peaks)

peaks$peak.name = paste0(peaks$chr,  "_", peaks$start, "_", peaks$end)
peaks$peak.name[kk] = paste0('tss.', peaks$peak.name[kk])

# remove duplications
jj = match(unique(peaks$peak.name), peaks$peak.name)
peaks = peaks[jj, ];
dim(peaks)

# remove contigs
peaks = peaks[grep('chr', peaks$chr), ]
dim(peaks)

# remove chrM
peaks = peaks[which(peaks$chr != 'chrM'), ]
dim(peaks)


peakDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using'
# preapre SAF input for featureCounts
require(Rsubread)

#counts = quantify.signals.within.peaks(peaks = peak.merged, bam.list=bam.list, rpkm.normalization = FALSE, isPairedEnd = FALSE)
SAF = data.frame(GeneID=peaks$peak.name, 
                 Chr=peaks$chr, 
                 Start=peaks$start, 
                 End=peaks$end, 
                 Strand=peaks$strand, stringsAsFactors = FALSE)

write.table(SAF, file = paste0(peakDir, '/Peaks_histMarkers_ATACseq_missedTSS.dupReduced.saf'), sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE) 

##########################################
# read counting within promoters
##########################################
Test.readCounting.for.promoters = FALSE
if(Test.readCounting.for.promoters){
  promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
  
  # clean peaks
  peaks = promoters
  peaks = data.frame(peaks)
  colnames(peaks)[c(1:3)] = c("chr", "start", "end")
  
  dim(peaks)
  
  peaks$peak.name = paste0(peaks$geneID, '_', peaks$geneSymbol, '_', peaks$chr,  "_", peaks$start, "_", peaks$end)
  
  jj = match(unique(peaks$peak.name), peaks$peak.name)
  peaks = peaks[jj, ];
  dim(peaks)
  
  peaks = peaks[which(peaks$chr != 'chrM'), ]
  
  dim(peaks)
  
  
  peakDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R12965_CT'
  # preapre SAF input for featureCounts
  require(Rsubread)
  
  #counts = quantify.signals.within.peaks(peaks = peak.merged, bam.list=bam.list, rpkm.normalization = FALSE, isPairedEnd = FALSE)
  SAF = data.frame(GeneID=peaks$peak.name, 
                   Chr=peaks$chr, 
                   Start=peaks$start, 
                   End=peaks$end, 
                   Strand=peaks$strand, stringsAsFactors = FALSE)
  
  write.table(SAF, file = paste0(peakDir, '/amex6_TSS_all.saf'), sep = '\t', row.names = FALSE, 
              col.names = TRUE, quote = FALSE) 
  
}

########################################################
########################################################
# Section III : processing read counts within consensus peaks 
# 
########################################################
########################################################
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

design = read.csv(file = paste0(dataDir, "R11876_R12810_R12965_CT_analysis_20220217_QCs_AK.csv"))
design = design[!is.na(design$sampleID), ]
design$fileName = paste0(design$condition, '_', design$sampleID)
design = design[, c(1:3, 15, 5, 4, 6:14)]
colnames(design)[5] = 'batch'


xlist<-list.files(path=paste0(dataDir, '/featurecounts_peaks.Q30'),
                  pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

colnames(design)[1] = 'SampleID'
#design = design[grep('90392|90393|108072|108073|108070|108071', design$SampleID, invert = TRUE),]

counts = process.countTable(all=all, design = design[, c(1, 4)])

save(design, counts, 
     file = paste0(RdataDir, '/histoneMarkers_samplesDesign_readCounts_peaks_histMarkers_ATACseq_missedTSS.dupReduced.Rdata'))

##########################################
#  signal normalization and QC with PCA
##########################################
load(file = paste0(RdataDir, '/histoneMarkers_samplesDesign_readCounts_peaks_histMarkers_ATACseq_missedTSS.dupReduced.Rdata'))
design$unique.rmdup = as.numeric(design$unique.rmdup)
jj = which(as.numeric(design$unique.rmdup) > 1000)
design$unique.rmdup[jj] = design$unique.rmdup[jj]/10^6

conds = c('H3K4me3', 'H3K4me1', 'H3K27me3', 'H3K27ac')
save.scalingFactors.for.deeptools = TRUE
Make.pca.plots = TRUE

for(n in 1:length(conds))
{
  # n = 2
  cat(n, ' --', conds[n], '\n')
  sels = which(design$marks == conds[n])
  
  design.sel = design[sels, ]
  counts.sel = counts[, c(1, (sels+1))]
  
  ss = apply(as.matrix(counts.sel[, -1]), 1, mean)
  
  par(mfrow=c(1,2))
  hist(log10(as.matrix(counts.sel[, -1])), breaks = 100, xlab = 'log10(nb of reads within peaks)', main = 'distribution')
  plot(ecdf(log10(as.matrix(counts.sel[, -1]) + 0.1)), xlab = 'log10(nb of reads within peaks)', main = 'cumulative distribution')
  
  ss = apply(as.matrix(counts.sel[, -1]), 2, sum)
  design.sel$usable.reads.withinPeaks = ss
  
  rownames(counts.sel) = counts.sel$gene
  dds <- DESeqDataSetFromMatrix(as.matrix(counts.sel[, -1]), DataFrame(design.sel), design = ~ condition)
  
  ss = rowSums(counts(dds))
  par(mfrow=c(1,1))
  hist(log10(ss), breaks = 100)
  hist(log10(ss[grep('tss', names(ss))]), breaks = 100, add = TRUE, col = 'darkgray')
  length(which(ss > quantile(ss[grep('tss', names(ss))], probs = 0.8)))
  
  dd0 = dds[ss > 100, ]
  dd0 = estimateSizeFactors(dd0)
  sizefactors.UQ = sizeFactors(dd0)
  
  sizeFactors(dds) <- sizefactors.UQ
  fpm = fpm(dds, robust = TRUE)
  #colnames(fpm) = design$condition
  
  if(save.scalingFactors.for.deeptools){
    if(n == 1){
      xx = data.frame(sampleID = design.sel$SampleID,  
                      scalingFactor = design.sel$unique.rmdup/(sizeFactors(dds)*median(as.numeric(design.sel$unique.rmdup))),
                      stringsAsFactors = FALSE)
    }else{
      xx0 = data.frame(sampleID = design.sel$SampleID,  
                         scalingFactor = design.sel$unique.rmdup/(sizeFactors(dds)*median(design.sel$unique.rmdup)),
                         stringsAsFactors = FALSE)
      xx = rbind(xx, xx0)
      
    }
    
  }
  
  if(Make.pca.plots){
    library(ggrepel)
    library(dplyr)
    library(tibble)
    
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
    
    pca=plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = FALSE)
    print(pca)
    
    pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE, ntop = 3000))
    pca2save$name = paste0(design.sel$condition, '_', design.sel$SampleID)
    
    ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
      geom_point(size=3) + 
      #geom_text_repel(size = 2.0, color = 'darkblue')
      geom_text(hjust = 0.2, nudge_y = 0.2, size=2.5)
    
    ggsave(paste0(resDir, "/histM_PCA_", conds[n], ".pdf"), width=16, height = 10)
  }
  
  
}

write.table(xx, file = paste0(resDir, '/histMarkers_DESeq2_scalingFactor_forDeeptools.txt'), sep = '\t',
            col.names = FALSE, row.names = FALSE, quote = FALSE)


##########################################
# test/explore histone markers for postional-related genes
# and regeneration-related gened
##########################################
fpm = readRDS(file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb.rds'))
res = data.frame(log2(fpm + 1), stringsAsFactors = FALSE)

res$gene = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[2]})


res$x1 = res$Hand_K4me3
res$x2 = res$UA_K4me3
rr = res$x1 - res$x2

examples.sel = which(abs(rr)>2 & (res$x1>5|res$x2>5))
examples.sel = unique(c(examples.sel, grep('HOXA11|HOXA9|HOXD13|HOXD11|HOXD9|MEIS', res$gene)))

ggplot(data=res, aes(x=x1, y=x2, label = gene)) +
  geom_point(size = 0.6) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=5, col='darkgray') +
  geom_hline(yintercept=5, col="darkgray") +
  labs(x = "Hand_H3K4me3", y= 'UA_H3K4me3')

ggsave(paste0(figureDir, "histMarker_H3K4me3_scatterplot_Hand.vs.Hand.pdf"), width=12, height = 8)

res$x1 = res$Hand_K27me3
res$x2 = res$UA_K27me3
rr = res$x1 - res$x2

examples.sel = which(abs(rr)>2 & (res$x1>5|res$x2>5))
examples.sel = unique(c(examples.sel, grep('HOXA11|HOXA9|HOXD13|HOXD11|HOXD9|MEIS', res$gene)))

ggplot(data=res, aes(x=x1, y=x2, label = gene)) +
  geom_point(size = 0.6) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=5, col='darkgray') +
  geom_hline(yintercept=5, col="darkgray") +
  labs(x = "Hand_H3K27me3", y= 'UA_H3K27me3')

ggsave(paste0(figureDir, "histMarker_H3K27me3_scatterplot_Hand.vs.Hand.pdf"), width=12, height = 8)

########################################################
########################################################
# Section : Consider the grouping: 
# inactive promoters, active promoter, house-keeping genes, mature-specific genes,
# regeneration genes, positional genes
########################################################
########################################################
fpm = readRDS(file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb.rds'))
res = data.frame(log2(fpm + 2^0), stringsAsFactors = FALSE)
res$gene = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[2]})
res$geneID  = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})

##########################################
# clean TSS, for each gene, only keep the TSS with sigals if there are mulitple ones
##########################################
Clean.TSS = TRUE

if(Clean.TSS){
  ggs = unique(res$geneID)
  sels = c()
  for(n in 1:length(ggs))
  {
    # n = 1
    jj = which(res$geneID == ggs[n])  
    if(length(jj)>0){
      ss = apply(as.matrix(res[jj, c(1:6)]), 1, sum)
      sels = c(sels, jj[which.max(ss)])
    } 
    
  }
}

res = res[sels, ]

saveRDS(res, file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb_TSSfiltered.rds'))

##########################################
# first try to define groups 
##########################################
res = readRDS(file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb_TSSfiltered.rds'))

aa = readRDS(file =  "../results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/TestStat_regeneration_RNAseq.rds") 

aa$gene = sapply(rownames(aa), function(x) unlist(strsplit(as.character(x), '_'))[1])
aa$geneID = sapply(rownames(aa), function(x) {test = unlist(strsplit(as.character(x), '_')); return(test[length(test)])})

aa[grep('DBP|SALL', rownames(aa)), grep('Mature_UA', colnames(aa))]

aa$expr.mUA = apply(aa[, grep('Mature_UA', colnames(aa))], 1, mean)
aa$expr.BLday5 = apply(aa[, grep('BL_UA_5days', colnames(aa))], 1, mean)
aa$expr.BLday9 = apply(aa[, grep('BL_UA_9days', colnames(aa))], 1, mean)
aa$expr.BLday13 = apply(aa[, grep('BL_UA_13days_proximal', colnames(aa))], 1, mean)
aa$expr.BLday13.distal = apply(aa[, grep('BL_UA_13days_distal', colnames(aa))], 1, mean)
aa = aa[, c(5:11, 23:29)]
aa = aa[match(res$geneID, aa$geneID), ]
res = data.frame(res, aa, stringsAsFactors = FALSE)

res$x1 = res$UA_K4me3
res$x2 = res$UA_K27me3
res$ratio = res$x1 - res$x2

examples.sel = c()
examples.sel = unique(c(examples.sel, grep('HOXA13|HOXA11|HOXA9|HOXD13|HOXD11|HOXD9|MEIS|SALL|DBP', res$gene)))


ggplot(data=res, aes(x=x1, y=x2, label = gene)) +
  geom_point(size = 0.25) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=4, col='darkgray') +
  geom_hline(yintercept=4, col="darkgray") +
  labs(x = "UA_H3K4me3", y= 'UA_H3K27me3')


ggplot(data=res, aes(x=ratio, y=expr, label = gene)) +
  geom_point(size = 0.25) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=4, col='darkgray') +
  geom_hline(yintercept=4, col="darkgray") +
  labs(x = "UA_H3K4me3/UA_H3K27me3", y= 'Expr')

load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
hkgs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
nonexp = controls.tissue$geneIDs[which(controls.tissue$tissues != 'housekeeping')]

res$groups = 'limb'
res$groups[!is.na(match(res$geneID, hkgs))] = 'house_keep'
res$groups[!is.na(match(res$geneID, nonexp))] = 'other_tissues'

xx = res[order(-res$expr), ]
head(xx[which(xx$groups != 'house_keep'), ])

ggplot(data=res, aes(x=expr.mUA, y=ratio, label = gene, color = groups)) +
  geom_point(size = 0.25) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=4, col='darkgray') +
  geom_hline(yintercept=4, col="darkgray") +
  labs(y = "UA_H3K4me3/UA_H3K27me3", x= 'Expr')

res[,c(1, 4, 7:ncol(res))] %>% 
  pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
  ggplot(aes(x = factor(groups, levels = c('other_tissues', 'house_keep', 'limb')), y=value, fill=markers)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 1.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14)) +
  labs(x = "", y= 'normalized data (log2)')

##########################################
# redefine different groups with RNA-seq data 
##########################################
Redefine.gene.groups.with.RNAseq = TRUE
if(Redefine.gene.groups.with.RNAseq){
  
  table(res$groups)
  kk = which(res$groups == 'limb')
  
  ss = apply(res[kk, grep('expr', colnames(res))], 1, max)
  
  select = kk[which(ss< -1 | is.na(ss))]
  res$groups[select] = 'lowlyExpr_limb'
  table(res$groups)
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  #jj = which(res$x1>4 & res$x2 >4)
  #xx = res[jj, ]
    
  kk = which(res$groups == 'limb')
  
  diffs = res$expr.BLday5[kk] - res$expr.mUA[kk]
  select = kk[which(res$expr.mUA[kk] < 0 & diffs >2)]
  
  res$groups[select] = 'upregulated.d5'
  table(res$groups)
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) >2)]
  res$groups[select] = 'downregulated.d5'
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 'upregulated.d5', 'downregulated.d5', 'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  
  kk = which(res$groups == 'limb')
  
  select = kk[which(res$expr.mUA[kk] < 0 & (res$expr.BLday5[kk] - res$expr.mUA[kk]) < 2 &  (res$expr.BLday9[kk] - res$expr.mUA[kk]) > 2)]
  res$groups[select] = 'upregulated.d9'
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) < 2 & (res$expr.mUA[kk] - res$expr.BLday9[kk]) >2)]
  res$groups[select] = 'downregulated.d9'
  
  
  kk = which(res$groups == 'limb')
  
  select = kk[which(res$expr.mUA[kk] < 0 & (res$expr.BLday5[kk] - res$expr.mUA[kk]) < 2 &  (res$expr.BLday9[kk] - res$expr.mUA[kk]) < 2 &
                    (res$expr.BLday13[kk] - res$expr.mUA[kk]) >2)]
  res$groups[select] = 'upregulated.d13'
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) < 2 & (res$expr.mUA[kk] - res$expr.BLday9[kk]) < 2 &
                (res$expr.mUA[kk] - res$expr.BLday9[kk]) > 2) ]
  res$groups[select] = 'downregulated.d13'
  
  
  
  
  table(res$groups)
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 
                                             'upregulated.d5', 'upregulated.d9','upregulated.d13', 
                                             'downregulated.d5',  'downregulated.d9',
                                             'downregulated.d13',
                                             'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  
  ggsave(paste0(figureDir, "histMarker_H3K27me3_H3K4me3_up.downregulatedGenes.UA.BL.days.pdf"), width=12, height = 8)
  
  fdr.cutoff = 0.05
  select = which(aa$padj < fdr.cutoff & aa$log2FC >1 & aa$log2FC.mUA.vs.others >0)
  res$groups[!is.na(match(res$geneID, aa$geneID[select]))] = 'mature_highlyExp'
  table(res$groups)
  
  select = which(aa$padj < fdr.cutoff & aa$log2FC >1 & aa$log2FC.mUA.vs.others < 0)
  res$groups[!is.na(match(res$geneID, aa$geneID[select]))] = 'regeneration'
  table(res$groups)
  
  select = which((aa$padj >= fdr.cutoff | aa$log2FC <=1) & log2(aa$baseMean) <4) 
  res$groups[!is.na(match(res$geneID, aa$geneID[select])) & res$groups == 'limb_static'] = 'limb_lowlyExp'
  
    
}


ggplot(data=res, aes(x=x1, y=x2, label = gene, color = groups)) +
  geom_point(size = 0.4) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=4, col='darkgray') +
  geom_hline(yintercept=4, col="darkgray") +
  labs(x = "UA_H3K4me3", y= 'UA_H3K27me3')


#xx = melt(res[], id.vars = c('UA_K4me3', 'UA_K27me3'), variable_name = 'markers')
res[,c(1, 4, 7:13)] %>% 
  pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
ggplot(aes(x = factor(groups, levels = c('other_tissues', 'house_keep', 'limb_lowlyExp',
                                          'mature_highlyExp', 'regeneration', 'limb_static')), y=value, fill=markers)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 1.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14)) +
  labs(x = "", y= 'normalized data (log2)')

  
