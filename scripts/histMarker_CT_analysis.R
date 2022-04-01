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

version.analysis = 'CT_merged_20220328'
#peakDir = "Peaks/macs2_broad"
saveTable = TRUE

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0(figureDir, 'tables4plots/')

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/'
#gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')

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
# update the design with the usable reads after merging the resequenced samples
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
design$usable = NA

xlist = list.files(path = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R13191_CT/QCs/BamStat',
                   pattern = '*stat.txt', full.names = TRUE)

for(n in 1:nrow(design))
{
  kk = grep(design$sampleID[n], xlist)
  cat(n, ' -- sample - ', design$sampleID[n], ' : ', basename(xlist[kk]), '\n')
  if(length(kk) == 1){
    tmp = read.table(xlist[kk], sep = '\t', header = TRUE)
    design$usable[n] = tmp$rmdup_uniq[1]/10^6
  }else{
    stop('Error')
  }
}

saveRDS(design, file = paste0(RdataDir, '/histM_CT_design_info.rds'))

##########################################
# process called peaks for H3K4me3 and H3K4me1
##########################################
design = readRDS(file = paste0(RdataDir, '/histM_CT_design_info.rds'))

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
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
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
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
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
  # H3K4me1 peak consensus
  ##########################################
  ## peaks_H3K4me1_mLA
  n = 2 
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, '2_p', pval.cutoff, '.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  peaks_H3K4me1_mLA = peaks[[kk[2]]]
  
  peaks.merged = peaks_H3K4me1_mLA
  #peaks.merged = union(peaks.merged, peaks_H3K4me1_mLA)
  peaks.merged = GenomicRanges::reduce(peaks.merged)
  
  # H3K4me1_mUA
  n = 5
  cc = conds[n]
  kk = which(design.sel$condition == cc)
  
  cat(n, ' -- ', cc, '--', length(kk), 'replicates\n')
  
  ol.peaks <- makeVennDiagram(peaks[kk[c(1, 2, 3)]], NameOfPeaks=names(peaks)[kk[c(1, 2, 3)]], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, '_1_p', pval.cutoff, '.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  ol.peaks <- makeVennDiagram(peaks[kk[c(2, 3,4)]], NameOfPeaks=names(peaks)[kk[c(2, 3, 4)]], connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_', cc, '_2_p', pval.cutoff, '.pdf'), 
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
  
  ##########################################
  # merge peaks from K4me3 and K4me1
  ##########################################
  peaks.all = union(peaks.merged.H3K4me3, peaks.merged.H3K4me1)
  
  save(peaks.all, file = paste0(RdataDir,  
                                   '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_pval.', pval.cutoff, '.Rdata'))
  
  saveRDS(peaks.all, file = paste0(RdataDir, 
                                   '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_pval.', pval.cutoff, '.rds'))
  
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

export(object = peaks.all,  con = paste0(resDir, "/histM_peaks_all_565k.bed"), format = 'bed')

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

cat(length(peaks.all), ' consensus peaks by adding atac-seq peaks \n')

#peaks.all = GenomicRanges::reduce(peaks.all)
saveRDS(peaks.all, file = paste0(RdataDir, 
                                 '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_and_ATACseq_consensusPeaks.rds'))

export(object = peaks.all,  con = paste0(resDir, "/histM_peaks_plus_atacseqPeaks_570k.bed"), format = 'bed')

##########################################
# compare the histM peaks considered with all annotated TSS
##########################################
peaks.all = readRDS(file = paste0(RdataDir, 
                                  '/histMarkers_macs2peaks_H3K4me3.H3K4me1_32samples_consensusPeaks_and_ATACseq_consensusPeaks.rds'))

cat(length(peaks.all), ' consensus peaks currently \n')

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

cat(length(pp), ' promoters in total \n')

saveRDS(pp, file = paste0(RdataDir, '/all_181985_TSS_transcriptID_geneID_with.1kb.width.rds'))

#pp = GenomicRanges::reduce(pp)

sum(overlapsAny(xx, peaks.all,  ignore.strand=TRUE))

tss_missed = pp[!overlapsAny(pp, peaks.all, ignore.strand=TRUE)]
cat(length(tss_missed), ' tss missed in consensus peaks \n')

# saveRDS(tss_missed, file = paste0(RdataDir, '/missedTSS_140k_inHistM_atacseq_peaks_transcriptID_geneID_with.1kb.width.rds'))

xx = GenomicRanges::reduce(tss_missed) 
saveRDS(xx, file = paste0(RdataDir, '/missedTSS_dupReduced_90k_inHistM_atacseq_peaks_transcriptID_geneID_with.1kb.width.rds'))

peaks_histM_atac_tss = union(peaks.all, tss_missed)

saveRDS(peaks_histM_atac_tss, file = paste0(RdataDir, '/Peaks_histMarkers_ATACseq_missedTSS.dupReduced_662k.rds'))

export(object = peaks_histM_atac_tss,  con = paste0(resDir, "/Peaks_histMarkers_ATACseq_missedTSS.dupReduced_662k.bed"), 
       format = 'bed')

##########################################
#  make saf files
##########################################
peaks = readRDS(file = paste0(RdataDir, '/Peaks_histMarkers_ATACseq_missedTSS.dupReduced_662k.rds'))

tss = readRDS(file = paste0(RdataDir, '/missedTSS_dupReduced_90k_inHistM_atacseq_peaks_transcriptID_geneID_with.1kb.width.rds'))

kk = which(overlapsAny(peaks, tss) == TRUE)

# clean peaks
peaks = data.frame(peaks)
colnames(peaks)[c(1:3)] = c("chr", "start", "end")
dim(peaks)

length(grep('tss', peaks$peak.name))

peaks$peak.name = paste0(peaks$chr,  "_", peaks$start, "_", peaks$end)
peaks$peak.name[kk] = paste0('tss.', peaks$peak.name[kk])

cat(nrow(peaks), ' peaks considered in total \n')

length(grep('tss', peaks$peak.name))

# remove duplications
jj = match(unique(peaks$peak.name), peaks$peak.name)
peaks = peaks[jj, ];
dim(peaks)

# remove contigs and many tss from contigs 
peaks = peaks[grep('chr', peaks$chr), ]
dim(peaks)

length(grep('tss', peaks$peak.name))

# remove chrM
peaks = peaks[which(peaks$chr != 'chrM'), ]
dim(peaks)

length(grep('tss', peaks$peak.name))

peakDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using'
# preapre SAF input for featureCounts
require(Rsubread)

#counts = quantify.signals.within.peaks(peaks = peak.merged, bam.list=bam.list, rpkm.normalization = FALSE, isPairedEnd = FALSE)
SAF = data.frame(GeneID=peaks$peak.name, 
                 Chr=peaks$chr, 
                 Start=peaks$start, 
                 End=peaks$end, 
                 Strand=peaks$strand, stringsAsFactors = FALSE)

length(grep('tss', SAF$GeneID))

write.table(SAF, file = paste0(peakDir, '/Peaks_histMarkers_ATACseq_missedTSS.dupReduced_v2.saf'), sep = '\t', row.names = FALSE, 
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
design = readRDS(file = paste0(RdataDir, '/histM_CT_design_info.rds'))
#design = read.csv(file = paste0(dataDir, "R11876_R12810_R12965_CT_analysis_20220217_QCs_AK.csv"))
#design = design[!is.na(design$sampleID), ]
design$fileName = paste0(design$condition, '_', design$sampleID)
design = design[, c(1:3, 15, 5, 16, 4, 6:14)]
colnames(design)[5] = 'batch'

## read counts within peaks
xlist<-list.files(path=paste0(dataDir, '/featurecounts_peaks.Q30'),
                  pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

colnames(design)[1] = 'SampleID'

counts = process.countTable(all=all, design = design[, c(1, 4)])

save(design, counts, 
     file = paste0(RdataDir, '/histoneMarkers_samplesDesign_readCounts_peaks_histMarkers_ATACseq_missedTSS.dupReduced.Rdata'))

##########################################
## First, for all samples together 
# - filter peaks with low number of reads but keep some background and the one overlapping atac-seq peaks 
# 
# Then for each marker: H3K4me3, H3K4me1, H3K27me3, H3K27ac and IgG
# - normalization,
# - QC with PCA and pairwise plot, filter bad samples
# - batch correction with Combat
# - save the batch-corrected matrix for mature sample comparison and regeneration comparison
# 
##########################################
load(file = paste0(RdataDir, '/histoneMarkers_samplesDesign_readCounts_peaks_histMarkers_ATACseq_missedTSS.dupReduced.Rdata'))
design$unique.rmdup = as.numeric(design$unique.rmdup)
jj = which(as.numeric(design$unique.rmdup) > 1000)
design$unique.rmdup[jj] = design$unique.rmdup[jj]/10^6

##########################################
# filter peaks below certain thrshold of read counts
##########################################
if(Filtering.peaks.with.lowReads){
  
  rownames(counts) = counts$gene
  dds <- DESeqDataSetFromMatrix(as.matrix(counts[, -1]), DataFrame(design), design = ~ condition)
    
  index_bgs = grep('tss', rownames(dds))
  
  ss = rowMaxs(counts(dds))
  
  par(mfrow=c(1,1))
  hist(log10(ss), breaks = 100, main = 'log10(max(reads within peaks across samples))')
  hist(log10(ss[index_bgs]), breaks = 100, add = TRUE, col = 'darkgray') 
  length(which(ss > quantile(ss[index_bgs], probs = 0.95)))
  
  # check the peak length
  peakNames = rownames(dds)
  peakNames = gsub('tss.', '', peakNames)
  pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  
  # rownames(dds) = paste0(pp$X1, ':', pp$X2, '-', pp$X3)
  # rownames(counts) = rownames(dds)
  
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  ll = width(pp)
  
  index_keep = c()
  
  # max read counts normalized by length (nb of reads within 500bp) 
  ss = rowMaxs(counts(dds))/ll*500 
  hist(log10(ss), breaks = 100, col = 'blue', main = 'log2(max of read counts within peaks) ')
  hist(log10(ss[index_bgs]), breaks = 100, add = TRUE, col = 'darkgray')
  abline(v = log10(c(10, 20, 30, 50)), col = 'red', lwd = 2.0)
  
  ## manually check the which cutoff is better
  #ii_test = rownames(dds)[which(ss>9.9 & ss<10.1)]
  #ii_test = rownames(dds)[which(ss>19.9 & ss<20.1)]
  ii_test = rownames(dds)[which(ss>29.9 & ss<30.1)] # cutoff.peak = 30 look good for H3K4me3
  
  cutoff.peak = 30 # cutoff.peak = 30 for H3K4me3  (atac-seq cutoff.peak > 40 for >= 2 samples) 
  #cutoff.bg = 20
  cat(length(which(ss > cutoff.peak)), 'peaks selected with minimum read of the highest peak -- ', cutoff.peak,  '\n')
  #cat(length(which(ss < cutoff.bg)), 'peaks selected with minimum read of the highest peak -- ', cutoff.bg,  '\n')
  
  hist(log10(ss), breaks = 100, col = 'blue', main = 'log2(max of read counts within peaks) ')
  hist(log10(ss[index_bgs]), breaks = 100, add = TRUE, col = 'darkgray')
  abline(v= log10(cutoff.peak), col = 'red', lwd = 2.0)
  
  #nb.above.threshold = apply(counts(dds), 1, function(x) length(which(x>cutoff.peak)))
  index_keep = which(ss > cutoff.peak)
  cat(length(index_keep), ' peaks will be retained \n')
  
  ### keep also peaks if they are overlapped by atac-seq peaks
  atacseq_peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/ATAC_allUsed_20220328/Rdata/',
                                        'ATACseq_peak_consensus_filtered_64k.rds'))
  
  index_keep = unique(c(index_keep, which(overlapsAny(pp, atacseq_peaks) == TRUE)))
  cat(length(index_keep), ' peaks will be retained \n')
  
  if(Select.background.for.peaks){
    #ii.bg = which(ss < cutoff.bg)
    #ii.bg = sample(index_bgs, size = 5000, replace = FALSE)
    #rownames(dds)[ii.bg] = paste0('tss.', rownames(dds)[ii.bg])
    
    index_keep = unique(c(index_keep, index_bgs))
    cat(length(index_keep), ' peaks will be retained \n')
    
    dds = dds[index_keep, ]
    ll.sels = ll[index_keep]
    
  }else{
    dds <- dds[ii, ]
    ll.sels = ll[ss >= cutoff.peak]
  }
  
  save(design, dds, file = paste0(RdataDir, 
                                  '/histoneMarkers_samplesDesign_ddsPeaksMatrix_filtered_incl.atacPeak.missedTSS.bgs_324k.Rdata'))
   
}


##########################################
# for each markers, normalization, sample filtering, and batch correction
##########################################
load(file = paste0(RdataDir, 
                   '/histoneMarkers_samplesDesign_ddsPeaksMatrix_filtered_incl.atacPeak.missedTSS.bgs_324k.Rdata'))

SaveHistM.peaks.afterFiltering = FALSE
if(SaveHistM.peaks.afterFiltering){
  peakNames = rownames(dds)
  peakNames = gsub('tss.', '', peakNames)
  pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  export(object = pp,  con = paste0(resDir, "/histM_peaks_all_afterFiltering_323k.bed"), format = 'bed')
  
  
}

conds = c('H3K4me3', 'H3K4me1', 'H3K27me3', 'H3K27ac', 'IgG')

save.scalingFactors.for.deeptools = TRUE
Make.pca.plots = TRUE
Filtering.peaks.with.lowReads = TRUE
Select.background.for.peaks = TRUE

for(n in 1:length(conds))
{
  # n = 5
  sels = which(dds$marks == conds[n] & dds$SampleID != '163655' & dds$SampleID != '185743')
  
  cat(n, ' --', conds[n], '--', length(sels), 'samples\n')
  
  ddx = dds[, sels]
  design.sel = design[sels, ]
  
  ddx$condition = droplevels(ddx$condition)
  ddx$batch = droplevels(ddx$batch)
  
  ## normalization  
  ss = rowSums(counts(ddx))
  
  jj = grep('tss', rownames(dds))
  hist(log10(ss[-jj]), breaks = 100);
  hist(log10(ss[jj]), breaks = 100, col = 'darkgray', add = TRUE);
  abline(v = c(log10(50), 2, log10(500), 3), col = 'red', lwd = 2.0)
  
  dd0 = ddx[ss > 0, ]
  dd0 = estimateSizeFactors(dd0)
  sizeFactors(ddx) <- sizeFactors(dd0)
  
  fpm = fpm(ddx, robust = TRUE)
  log2(range(fpm[which(fpm>0)]))
  
  fpm = log2(fpm + 2^-3)
  colnames(fpm) = paste0(design.sel$condition, '_', design.sel$SampleID, '_', design.sel$batch)
  
  ## double check the sample qualities with PCA and scatterplots of replicates
  if(Make.pca.plots){
    library(ggrepel)
    library(dplyr)
    library(tibble)
    
    vsd <- varianceStabilizingTransformation(ddx, blind = FALSE)
    
    pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE, ntop = 3000))
    pca2save$name = paste0(design.sel$sample, '_', design.sel$SampleID)
    
    ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
      geom_point(size=3) + 
      geom_text_repel(size = 4.0)
      #geom_text(hjust = 0.4, nudge_y = 1.0, size=3.5)
    
    ggsave(paste0(resDir, "/histM_PCA_peakFiltered_", conds[n], ".pdf"), width=10, height = 6)
    
    pdfname = paste0(resDir, "/histM_peakFiltered_", conds[n], ".pdf")
    pdf(pdfname, width = 10, height = 10)
    
    plot.pair.comparison.plot(fpm[c(1:20000), grep('_mUA', colnames(fpm))], linear.scale = FALSE)
    plot.pair.comparison.plot(fpm[c(1:20000), grep('_mLA|_mHand', colnames(fpm))], linear.scale = FALSE)
    plot.pair.comparison.plot(fpm[c(1:20000), grep('_BL5|_BL9', colnames(fpm))], linear.scale = FALSE)
    plot.pair.comparison.plot(fpm[c(1:20000), grep('_BL13', colnames(fpm))], linear.scale = FALSE)
    
    dev.off()
    
  }
  
  ## batch correction
  Batch.Correction.combat = TRUE
  if(Batch.Correction.combat){
    require(edgeR)
    require(sva)
    source('Functions_atac.R')
    
    d <- DGEList(counts=counts(ddx), group=ddx$condition)
    if(conds[n] == 'H3K27ac' | conds[n] == 'IgG'){
      tmm <- calcNormFactors(d, method='TMMwsp') # TMMwsp used for H3K27ac
    }else{
      tmm <- calcNormFactors(d, method='TMM')
    }
    tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.01)
    
    # ii_bgs = grep('tss', rownames(tmm))
    # tmm = tmm[-ii_bgs, ]
    
    design.sel$condition = droplevels(design.sel$condition)
    design.sel$batch = droplevels(design.sel$batch)
    
    bc = as.factor(design.sel$batch)
    mod = model.matrix(~ as.factor(condition), data = design.sel)
    tmm = as.matrix(tmm)
    # vars = apply(tmm, 1, var)
    # sums = 
    # tmm = tmm[which(vars>0.5), ]
    # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
    # because there is no better batche other others 
    #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
    fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
    
    #design.tokeep<-model.matrix(~ 0 + conds,  data = design.sels)
    #cpm.bc = limma::removeBatchEffect(tmm, batch = bc, design = design.tokeep)
    # plot(fpm.bc[,1], tmm[, 1]);abline(0, 1, lwd = 2.0, col = 'red')
    make.pca.plots(tmm, ntop = 5000, conds.sel = as.character(unique(design.sel$condition)))
    ggsave(paste0(resDir, "/histMarker_", conds[n], "_beforeBatchCorrect_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    make.pca.plots(fpm.bc, ntop = 5000, conds.sel = as.character(unique(design.sel$condition)))
    ggsave(paste0(resDir, "/histMarker_", conds[n], "_afterBatchCorrect_",  version.analysis, ".pdf"), width = 10, height = 6)
    
    fpm = fpm.bc
    rm(fpm.bc)
    
    saveRDS(fpm, file = paste0(RdataDir, '/fpm_bc_TMM_combat_', conds[n], '_', version.analysis, '.rds'))
    saveRDS(design.sel, file = paste0(RdataDir, '/design.sels_bc_TMM_combat_', conds[n], '_', version.analysis, '.rds'))
    
  }
  
  ## save scaling factor for deeptools to make bigwig
  if(save.scalingFactors.for.deeptools){
    if(n == 1){
      sfs = data.frame(sampleID = design.sel$SampleID,  
                      scalingFactor = design.sel$usable/(sizeFactors(ddx)*median(as.numeric(design.sel$usable))),
                      stringsAsFactors = FALSE)
      
    }else{
      sfs0 = data.frame(sampleID = design.sel$SampleID,  
                       scalingFactor = design.sel$usable/(sizeFactors(ddx)*median(as.numeric(design.sel$usable))),
                       stringsAsFactors = FALSE)
      sfs = rbind(sfs, sfs0)
    }
  }
  
}

write.table(sfs, file = paste0(resDir, '/histMarkers_DESeq2_scalingFactor_forDeeptools.txt'), sep = '\t',
            col.names = FALSE, row.names = FALSE, quote = FALSE)

##########################################
# subtrat the IgG signals when there is a peaks there in IgG for each markers and each condition  
##########################################
load(file = paste0(RdataDir, 
                   '/histoneMarkers_samplesDesign_ddsPeaksMatrix_filtered_incl.atacPeak.missedTSS.bgs_324k.Rdata'))
igg = readRDS(file = paste0(RdataDir, '/fpm_bc_TMM_combat_', 'IgG', '_', version.analysis, '.rds'))

conds_histM = unique(design$sample)

ctl.means = c()
for(ii in 1:length(conds_histM)) 
{
  kk = grep(conds_histM[ii], colnames(igg))
  if(length(kk)>1) {
    ctl.means = cbind(ctl.means, apply(igg[, kk], 1, mean))
  }else{
    ctl.means = cbind(ctl.means, igg[, kk])
  }
}

colnames(ctl.means) = conds_histM

### Reason : reduce the ctl.mean with a fixed global factor and then fix the IgG signals
### scale each marker data based on the common IgG background 
ctl.means = ctl.means - 3.5
igg_peak = 4.0 - 3.5 # cutoff for peaks in IgG
#plot.pair.comparison.plot(ctl.means[c(1:20000), ], linear.scale = FALSE)

markers = c('H3K4me3', 'H3K4me1', 'H3K27me3', 'H3K27ac')

for(n in 1:length(markers))
{
  # n = 4
  cpm = readRDS(file = paste0(RdataDir, '/fpm_bc_TMM_combat_', markers[n], '_', version.analysis, '.rds'))
  design.sel = readRDS(file = paste0(RdataDir, '/design.sels_bc_TMM_combat_', markers[n], '_', version.analysis, '.rds'))
  cpm0 = cpm
  
  inputs = ctl.means[match(rownames(cpm), rownames(ctl.means)), ]
  
  ### global cutoff is probably more robust
  marker.means = c()
  for(ii in 1:length(conds_histM)) 
  {
    kk = grep(conds_histM[ii], colnames(cpm))
    if(length(kk)>1) {
      marker.means = cbind(marker.means, apply(cpm[, kk], 1, mean))
    }else{
     marker.means = cbind(marker.means, cpm[, kk])
    }
  }
  colnames(marker.means) = conds_histM
  
  marker.means = apply(marker.means, 1, mean)
  inputs.means = apply(inputs, 1, mean)
  
  #hist(marker.means, breaks = 100, xlim = range(c(marker.means, inputs.means)), ylim = c(0, 2*10^4))
  #hist(inputs.means, breaks = 100, col = 'gray', add = TRUE)
  
  ## find global scaling factors for between this marker and IgG
  ## by assuming that only few 100 peaks solely due to IgG background instead of real peaks
  ratios = inputs.means - marker.means
  
  plot(inputs.means, ratios, cex = 0.2)
  abline(v = c(2, 2.5, 3.5, 4), col = 'red')
  
  igg_cutoff = 2.5 # global factor for IgG
  cat(length(which(inputs.means> igg_cutoff)), ' peaks considered solely due to IgG \n')
  
  global_scaling = median(ratios[which(inputs.means > igg_cutoff)])
  cat('global scaling factor -- ', global_scaling, '\n')
  abline(h = global_scaling, col = 'blue')
  
  #bgs = inputs.means - global_scaling
  #igg_peak_new = igg_peak - median(ratios[which(inputs.means > igg_cutoff)])
  marker.means_new = marker.means + global_scaling
  ratios_new = inputs.means - marker.means_new
  cat('after sacling the markers ratios becomes : ', median(ratios_new[which(inputs.means > igg_cutoff)]))
  
  hist(marker.means_new, breaks = 100)
  hist(inputs.means, breaks = 100, add = TRUE, col = 'red')
  abline(v = igg_peak, col = 'blue')
  
  cpm = cpm + global_scaling
  
  for(m in 1:length(conds_histM))
  {
    # m = 1
    jj = grep(conds_histM[m], colnames(cpm))
    
    bgs = inputs[, which(colnames(inputs) == conds_histM[m])]
    # bgs = bgs - global_scaling
    j_cor = which(bgs > igg_peak)
    
    cat(markers[n], '-',  as.character(conds_histM[m]), '-',
        length(j_cor), ' peaks will subtract IgG signals  \n')
    
    for(j in jj) 
    {
      cpm[j_cor, j] = cpm[j_cor, j] - bgs[j_cor]
    }
    
  }
  
  saveRDS(cpm, file = paste0(RdataDir, '/fpm_bc_TMM_combat_', markers[n], '_IgG.subtrated_', 
                        version.analysis, '.rds'))
  
}

########################################################
########################################################
# Section IV : start to characterize the mUA to characterize the mature samples
# The underlying hypothesise is that the difference between UA, LA and Hand are minor 
# so once we have good idea of the global picture of mUA, it would be much easier to understand the segment-specific histone markers
# also to describe the changes during regneration
########################################################
########################################################
conds = c('H3K4me3', 'H3K4me1', 'H3K27me3',   'H3K27ac')
sample.sel = 'mUA'

keep = c()
for(n in 1:length(conds))
{
  # n = 1
  cpm = readRDS(file = paste0(RdataDir, '/fpm_bc_TMM_combat_', conds[n], '_', version.analysis, '.rds'))
  design.sel = readRDS(file = paste0(RdataDir, '/design.sels_bc_TMM_combat_', conds[n], '_', version.analysis, '.rds'))
  
  ### extract sample means  
  kk = grep(sample.sel, colnames(cpm))
  if(length(kk) == 0) cat('No sample found ! ERROR !!! \n')
  cat(length(kk), ' samples found for ', conds[n], ' - ', sample.sel, '\n')
  
  if(length(kk)>1) {
    keep = cbind(keep, apply(cpm[, kk], 1, mean))
  }else{
    keep = cbind(keep, cpm[, kk])
  }
  
  library(edgeR)
  library(qvalue)
  
}  

colnames(keep) = conds

##########################################
# focus on the atac-seq overlapped peaks
##########################################
atacseq_peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rxxxx_R10723_R11637_R12810_atac/Rdata/',
                                      'ATACseq_peak_consensus_filtered_55k.rds'))
ii_bgs = grep('tss.', rownames(keep))
rownames(keep) = gsub('tss.', '', rownames(keep))

pp = data.frame(t(sapply(rownames(keep), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
pp$strand = '*'
pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")
lls = width(pp)

for(n in 1:ncol(keep))
{
  keep[,n] = keep[,n] + log2(1000/lls)
}

cpm_bgs = keep[ii_bgs, ]
keep = keep[-ii_bgs, ]
pp = pp[-ii_bgs]


ii_overlap = which(overlapsAny(pp, atacseq_peaks) == TRUE)
cpm_nonoverlap = keep[-ii_overlap, ]
keep = keep[ii_overlap, ]


### scale the hisone marker signals to 0 and 1
Transform.histMarkers.to.mitigate.different.background.dynamicRanges = FALSE
if(Transform.histMarkers.to.mitigate.different.background.dynamicRanges){
  cal_transform_histM = function(x, cutoff.min = 3, cutoff.max = 6)
  {
    # x = keep[,1];cutoff.min = 3.5; cutoff.max = 6
    x[which(x<cutoff.min)] = cutoff.min
    x[which(x>cutoff.max)] = cutoff.max
    x = (x - cutoff.min)/(cutoff.max - cutoff.min)
    
  }
  
  xx = keep
  ## c('H3K4me3', 'H3K4me1', 'H3K27me3',   'H3K27ac')
  n = 1
  hist(cpm_nonoverlap[,n], breaks = 100, col = 'darkred'); 
  hist(xx[,n], breaks = 100, col = 'darkblue', add = TRUE);
  hist(cpm_bgs[,n], breaks = 50, add = TRUE, col = 'darkgray')
  abline(v = c(2, 1.5), col = 'blue')
  
  xx[, n] = cal_transform_histM(keep[, n], cutoff.min = 2., cutoff.max = 4) 
  
  n = 2
  hist(cpm_nonoverlap[,n], breaks = 100, col = 'darkred'); 
  hist(xx[,n], breaks = 100, col = 'darkblue', add = TRUE);
  hist(cpm_bgs[,n], breaks = 50, add = TRUE, col = 'darkgray')
  
  xx[, 2] = cal_transform_histM(keep[, 2], cutoff.min = 1.5, cutoff.max = 3.5) 
  
  n = 3
  hist(cpm_nonoverlap[,n], breaks = 100, col = 'darkred'); 
  hist(xx[,n], breaks = 100, col = 'darkblue', add = TRUE);
  hist(cpm_bgs[,n], breaks = 50, add = TRUE, col = 'darkgray')
  
  xx[, 3] = cal_transform_histM(keep[, 3], cutoff.min = 2., cutoff.max = 4) 
  
  n = 4
  hist(cpm_nonoverlap[,n], breaks = 100, col = 'darkred'); 
  hist(xx[,n], breaks = 100, col = 'darkblue', add = TRUE);
  hist(cpm_bgs[,n], breaks = 50, add = TRUE, col = 'darkgray')
  
  xx[, 4] = cal_transform_histM(keep[, 4], cutoff.min = 1.5, cutoff.max = 3.5) 
  
}

##########################################
# visualize mUA chromatin states defined by 4 markers
##########################################
yy = xx[c(1:25000), ]
nb_clusters = 8

library(dendextend)
my_hclust_gene <- hclust(dist(yy), method = "complete")

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)

my_gene_col <- data.frame(cluster =  paste0('cluster_', my_gene_col))
rownames(my_gene_col) = rownames(yy)

df <- data.frame(conds)
rownames(df) = colnames(yy)
colnames(df) = 'histMarker'

col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
          "#33a02c", "#fb9a99", "#e31a1c",
          "#fdbf6f", "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99", "#b15928")

sample_colors = c('springgreen4', 'steelblue2', 'red', 'darkgreen')
names(sample_colors) = conds
cluster_col = col3[1:nb_clusters]
names(cluster_col) = paste0('cluster_', c(1:nb_clusters))
annot_colors = list(
  histMarker = sample_colors,
  cluster = cluster_col)

#gaps.col = c()
#col<- colorRampPalette(c("steelblue", "white", "darkred"))(8)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8)
col = colorRampPalette(c('darkgray',  "orange"))(20)
plt = pheatmap(yy, annotation_row = my_gene_col, 
               annotation_col = df, show_rownames = FALSE, scale = 'none', 
               color = col, 
               show_colnames = FALSE,
               cluster_rows = TRUE, cluster_cols = FALSE,  
               clustering_method = 'complete', cutree_rows = nb_clusters, 
               annotation_colors = annot_colors
               #gaps_col = gaps.col
               ) 


row_order = data.frame(plt$tree_row$order, plt$tree_row$labels, stringsAsFactors = FALSE)
colnames(row_order) = c('index', 'name')
row_order$name = row_order$name[row_order$index]
row_order$cluster = my_gene_col$cluster[match(row_order$name, rownames(my_gene_col))]

callback = function(hc, mat){
  #sv = svd(t(mat))$v[,1]
  #clusters_order = 
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
  
  # o1 = c()
  # for(co in paste0('cluster_', c(1, 6, 5, 4, 3, 2))) {
  #   o1 = c(o1, row_order$index[which(row_order$cluster == co)])
  #   #row_order = row_order[o1, ]
  # }
  # dend = reorder(as.dendrogram(hc), wts = o1)
  # as.hclust(dend)
  
}


pheatmap(yy, annotation_row = my_gene_col, 
         annotation_col = df, show_rownames = FALSE, scale = 'none', 
         color = col, 
         show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = FALSE,  
         clustering_method = 'complete', cutree_rows = nb_clusters, 
         annotation_colors = annot_colors, 
         clustering_callback = callback,
         #gaps_col = gaps.col, 
         filename = paste0(figureDir, '/heatmap_mUA_chromatinStates_histoneMarkers.pdf'), 
         width = 6, height = 12)


#write.csv(xx, file = paste0(saveDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters', 
#                            nb_clusters, '.csv'), quote = FALSE, row.names = TRUE)

#col = colorRampPalette(c("black",  "orange"))(5)
#col = colorRampPalette((brewer.pal(n = 7, name ="YlOrRd")))(10)
pheatmap(xx[c(1:20000), ], show_rownames = FALSE, show_colnames = TRUE, 
         cluster_cols = FALSE, 
         col = col,
         scale = 'none')

########################################################
########################################################
# Section V : segment-specific histone markers
# 
########################################################
########################################################
library(edgeR)
library(qvalue)
require(corrplot)
require(pheatmap)
require(RColorBrewer)

atacseq_peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rxxxx_R10723_R11637_R12810_atac/Rdata/',
                                      'ATACseq_peak_consensus_filtered_55k.rds'))

conds_histM = c('H3K4me3', 'H3K27me3',  'H3K4me1', 'H3K27ac')

for(n_histM in 1:length(conds_histM))
{
  # n_histM = 4
  cpm = readRDS(file = paste0(RdataDir, '/fpm_bc_TMM_combat_', markers[n], '_IgG.subtrated_', version.analysis, '.rds'))
  design.sel = readRDS(file = paste0(RdataDir, '/design.sels_bc_TMM_combat_', conds_histM[n_histM], '_', version.analysis, '.rds'))
  
  ### select the samples and extract sample means
  conds = c("mUA", "mLA", "mHand")
  sample.sels = c();  
  cc = c()
  sample.means = c()
  
  for(n in 1:length(conds)) 
  {
    kk = grep(conds[n], colnames(cpm))
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
  f = factor(cc, levels= c('mUA', 'mLA', 'mHand'))
  
  mod = model.matrix(~ 0 + f)
  colnames(mod) = c('mUA', 'mLA', 'mHand')
  
  #To make all pair-wise comparisons between the three groups one could proceed
  fit <- lmFit(logCPM, mod)
  contrast.matrix <- makeContrasts(mLA - mUA, 
                                   mHand - mUA, 
                                   mHand - mLA, levels=mod)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  res = data.frame(fit2$p.value)
  colnames(res) = paste0(c('mLA.vs.mUA', 'mHand.vs.mUA', 'mHand.vs.mLA'), '.pval')
  
  xx = topTable(fit2, coef = 1, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.mLA.vs.mUA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  xx = topTable(fit2, coef = 2, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.mHand.vs.mUA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  xx = topTable(fit2, coef = 3, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.mHand.vs.mLA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  res$fdr.mean = apply(as.matrix(res[, grep('adj.P.Val', colnames(res))]), 1, function(x) return(mean(-log10(x))))
  res$logFC.mean =  apply(as.matrix(res[, grep('logFC', colnames(res))]), 1, function(x) return(mean(abs(x))))
  
  res$log2fc = apply(sample.means, 1, function(x) max(x) - min(x))
  res$maxs = apply(sample.means, 1, max)
  res$mins = apply(sample.means, 1, min)
  
  xx = data.frame(cpm[match(rownames(res), rownames(cpm)), ],  res, stringsAsFactors = FALSE) 
  res = xx
  
  saveRDS(res, file = paste0(RdataDir, '/fpm_bc_TMM_combat_DBedgeRtest_', conds_histM[n_histM], '_', version.analysis, '.rds'))
  
  #### select the significant peaks
  fdr.cutoff = 0.05; logfc.cutoff = 1
  select = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
               (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
               (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)
  )
  cat(length(select), ' DE ', conds_histM[n_histM],  ' \n')
  res = res[order(-res$log2fc), ]
  sample.means = sample.means[match(rownames(res), rownames(sample.means)), ]
  
  ####### focus on the atac-seq peak sets
  Keep.only.overlapped.by.atacseq.peaks = TRUE
  if(Keep.only.overlapped.by.atacseq.peaks){
    ii_bgs = grep('tss.', rownames(res))
    
    rownames(res) = gsub('tss.', '', rownames(res))
    
    pp = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    lls = width(pp)
    
    #for(n in 1:ncol(keep))
    #{
    #  keep[,n] = keep[,n] + log2(1000/lls)
    #}
    #cpm_bgs = keep[ii_bgs, ]
    #keep = keep[-ii_bgs, ]
    #pp = pp[-ii_bgs]
    
    ii_overlap = which(overlapsAny(pp, atacseq_peaks) == TRUE)
    #cpm_nonoverlap = keep[-ii_overlap, ]
    #keep = keep[ii_overlap, ]
    kk = match(rownames(res), names(pp)[ii_overlap])
    res = res[!is.na(kk), ]
    sample.means = sample.means[!is.na(kk), ]
      
    fdr.cutoff = 0.05; logfc.cutoff = 1
    select = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
                     (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
                     (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)
    )
    cat(length(select), ' DE ', conds_histM[n_histM],  ' \n')
   
  }
  
  cal_transform_histM = function(x, cutoff.min = 1.5, cutoff.max = 4, toScale = FALSE)
  {
    # x = keep[,1];cutoff.min = 3.5; cutoff.max = 6
    x[which(x<cutoff.min)] = cutoff.min
    x[which(x>cutoff.max)] = cutoff.max
    if(toScale) x = (x - cutoff.min)/(cutoff.max - cutoff.min)
    return(x)
    
  }
  
  yy = res[select, grep(conds_histM[n_histM], colnames(res))]
  #yy = apply(yy, 2, cal_transform_histM)
  
  df = as.data.frame(sapply(colnames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'segments'
  rownames(df) = colnames(yy)
  
  sample_colors = c('springgreen4', 'steelblue2', 'gold2')
  annot_colors = list(segments = sample_colors)
  
  pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
           show_colnames = FALSE,
           scale = 'row',
           cluster_cols=FALSE, annotation_col=df,
           #annotation_colors = annot_colors,
           width = 6, height = 12, 
           filename = paste0(figureDir, '/heatmap_histoneMarker_', conds_histM[n_histM], '_overlappedWithAtacseqPeak.pdf'))
  
  
}  


########################################################
########################################################
# Section VI : regeneration-related histone markers
# 
########################################################
########################################################
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


