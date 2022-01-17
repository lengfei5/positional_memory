##########################################################################
##########################################################################
# Project: positional memory
# Script purpose: QCs and processing for ATAC-seq data before the real analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jan 28 10:56:54 2021
##########################################################################
##########################################################################
rm(list=ls())

version.analysis = 'Rxxxx_R10723_R11637_R12810_atac'
#peakDir = "Peaks/macs2_broad"

resDir = paste0("../results/", version.analysis)

RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/'

source('Functions_atac.R')
require(ggplot2)
require(GenomicFeatures)

########################################################
########################################################
# Section : sequencing quality controls
# fragment size distribution 
# sequence saturation analysis
########################################################
########################################################
design_file = paste0(dataDir, 'design_sampleInfo_all.txt')


Collect.design.stat.nf.out = TRUE

if(Collect.design.stat.nf.out){
  
  #design = read.table(paste0(dataDir, 'sampleInfo_parsed.txt'), sep = '\t', header = TRUE)
  design = read.delim(design_file, sep = '\t', header = TRUE)
  
  # Rxxxx, R10723 and R11637 technical replicates merged 
  stats = read.table(paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11637_atac/', 
                            'nf_out/result/countStatTable.txt'), sep = '\t', header = TRUE)
  
  # 2022 new data R12810
  xx = read.table(paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R12810_atac/', 
                         'nf_out/result/countStatTable.txt'), sep = '\t', header = TRUE)
  
  stats = rbind(stats, xx)
  colnames(stats)[c(1, 3)] = c('fileName', 'trimmed')
  
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
  colnames(stats)[c(2,3)] = c('condition', 'samples')
  stats$samples = paste0(stats$condition, '_', stats$sampleID)
  
  stats = stats[!is.na(stats$total), ]
  write.csv(stats, file = paste0(resDir, '/R11637_R12810_atac_QCs_stats.csv'), row.names = FALSE)
  
  plot(stats$usable, stats$mapped); text(stats$usable, stats$mapped, stats$samples)
  
  #save(stats, file = paste0(RdataDir, '/R11876_CutTag_samples_design_stats.Rdata'))
  save(stats, file = paste0(RdataDir, '/R11637_atacseq_samples_design_stats.Rdata'))
  
}else{
  
  extract.stat.for.samples.manual(design_file, dataDir)
  
}

########################################################
########################################################
# Section : Quantify peak signals
# 1) consensus peaks
# 2) quantify the read counts within merged peaks
# 3) overview sample quality with PCA
########################################################
########################################################
source('functions_chipSeq.R')

peakDir = paste0(dataDir,  'calledPeaks/macs2')
peak.files = list.files(path = peakDir,
                        pattern = '*_peaks.xls', full.names = TRUE)

design = readRDS(file = paste0(RdataDir, '/design_merged_technicalReplicates_Rxxxx_R10723_R11637.rds'))
#mm = match(design$sampleID, as.character(c(161523:161526)))
#design = design[is.na(mm), ]
#saveRDS(design, file = paste0(RdataDir, '/design_merged_technicalReplicates_Rxxxx_R10723_R11637.rds'))

load(file = paste0(RdataDir, '/R11637_atacseq_samples_design_stats.Rdata'))
stats = stats[is.na(match(stats$sampleID, design$sampleID)), ]
stats = stats[, c(1, 2, 3, 8, 9)]
colnames(stats) = colnames(design)

stats$unique = stats$unique/10^6
stats$usable = stats$usable/10^6

design = rbind(design, stats)
design = data.frame(design, stringsAsFactors = FALSE)
saveRDS(design, file = paste0(RdataDir, '/design_merged_technicalReplicates_Rxxxx_R10723_R11637_R12810.rds'))


index  = c()
for(n in 1:nrow(design))
{
  test = grep(design$sampleID[n], peak.files)
  if(length(test) != 1) {
    cat(length(test), 'peak files Found \n')
  }else{
    index = c(index, test)
  }
}
peak.files = peak.files[index]


##########################################
# Manually identify consensus peaks across replicates taking into account of different batches and sequencing depth
##########################################
Manually.identify.peak.consensus = FALSE

# union of all peaks from all replicates
peak.merged = merge.peaks.macs2(peak.files, pcutoff = 6)

if(Manually.identify.peak.consensus){
  peaks = c()
  pval.cutoff = 4
  
  for(n in 1:length(peak.files)) 
  {
    cat(n, '\n')
    p = readPeakFile(peak.files[n], as = "GRanges");
    #eval(parse(text = paste0("p = pp.", k)));
    with.p.values = "X.log10.pvalue." %in% colnames(mcols(p))
    if(with.p.values) {
      p <- p[mcols(p)[,"X.log10.pvalue."] > pval.cutoff];
      p = reduce(p);
      #peaks10= c(peaks10, p10);
    }else{ 
      cat("no p values conlumn found for -- ", design.matrix$file.name[k], "\n");
      PLOT.p10 = FALSE;
    }
    #p = reduce(p)
    peaks= c(peaks, p)
  }
  
  names(peaks) = design$fileName
  
  saveRDS(peaks, file = paste0(RdataDir, '/macs2_peaks_mergedTechnialReps_34samples_pval.', pval.cutoff, '.rds'))
  
  peaks = readRDS(file = paste0(RdataDir, '/macs2_peaks_mergedTechnialReps_30samples_pval.', pval.cutoff, '.rds'))
  
  #xx = readRDS(file = paste0(RdataDir, '/macs2_peaks_mergedTechnialReps_30samples.rds'))
  #peaks = readRDS(file = paste0(RdataDir, '/macs2_peaks_mergedTechnialReps_30samples_pval0.001.rds'))
  #peaks = readRDS(file = paste0(RdataDir, '/macs2_peaks_mergedTechnialReps_30samples_pval0.001.rds'))
  
  ##########################################
  #  # try to merge BL time series
  ##########################################
  kk = which(design$condition == 'BL_UA_13days_proximal')
  
  source('functions_chipSeq.R')
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_proximal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_BL.UA.D13.proximal.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  #bld13.p = peaks[[3]][overlapsAny(peaks[[3]], peaks[[4]])]
  bld13.p = GenomicRanges::intersect(peaks[[3]], peaks[[4]], ignore.strand=TRUE)
  
  kk = which(design$condition == 'BL_UA_13days_distal')
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_distal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_BL_UA_13days_distal.pdf'),  height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  #bld13.d = peaks[[1]][overlapsAny(peaks[[1]], peaks[[2]])]
  bld13.d = GenomicRanges::intersect(peaks[[1]], peaks[[2]])
  
  kk = which(design$condition == 'BL_UA_9days')
  
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_distal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if (length(dev.list())!=0) {dev.off()}
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_BL_UA_9days.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  #bld9 = peaks[[10]][overlapsAny(peaks[[10]], peaks[[9]])]
  bld9 = GenomicRanges::intersect(peaks[[10]], peaks[[9]])
  
  kk = which(design$condition == 'BL_UA_5days')
  
  kk = c(7, 8)
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_distal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if (length(dev.list())!=0) {dev.off()}
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_BL_UA_5days_olds.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  #bld5 = peaks[[7]]
  #bld5 = bld5[overlapsAny(bld5, peaks[[8]])]
  bld5 = GenomicRanges::intersect(peaks[[7]], peaks[[8]])
  
  kk = c(5, 6)
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_distal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if (length(dev.list())!=0) {dev.off()}
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_BL_UA_5days_new.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  bld52 = GenomicRanges::intersect(peaks[[5]], peaks[[6]])
  
  
  ol.peaks <- makeVennDiagram(list(bld5, bld52), NameOfPeaks=c('day5.old', 'day5.new'), connectedPeaks="keepAll", main='BL_UA_day13_distal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if (length(dev.list())!=0) {dev.off()}
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_BL_UA_5days_new_vs_old.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  #length(bld5[overlapsAny(bld5, bld52)])
  bld5.old = bld5
  bld5 = bld52
  
  save(bld5, bld5.old, bld9, bld13.p, bld13.d, 
       file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 'regeneration.Rdata'))
  
  ##########################################
  # manual define the overlapping peaks for embryo stages
  ##########################################
  kk = which(design$condition == 'Embryo_Stage40')
  kk = c(11, 12)
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_proximal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Embryo_Stage40_new.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  es40.new = GenomicRanges::intersect(peaks[[11]], peaks[[12]])
 
  kk = c(13, 14)
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_proximal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Embryo_Stage40_old.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  es40.old = GenomicRanges::intersect(peaks[[13]], peaks[[14]])
  
  ol.peaks <- makeVennDiagram(list(es40.old, es40.new), NameOfPeaks=c('es40.old', 'es40.new'), connectedPeaks="keepAll")
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Embryo_Stage40_new_vs_old.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  
  # here for Embryo.stage40, new batche and old merged
  es40 = union(es40.new, es40.old)
  #es40 = es40.new
  
  kk = which(design$condition == 'Embryo_Stage44_proximal')
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_proximal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Embryo_Stage44_proximal.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  es44.p = intersect(peaks[[17]], peaks[[18]])
  
  kk = which(design$condition == 'Embryo_Stage44_distal')
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll")
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Embryo_Stage44_distal.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
    
  es44.d = intersect(peaks[[15]], peaks[[16]])
  
  save(es40, es44.d, es44.p, es40.old, es40.new, 
       file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 
                     'embryoStage.Rdata'))
  
  ##########################################
  # mature samples
  ##########################################
  kk = which(design$condition == 'Mature_UA')
  kk = c(27, 30)
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_distal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Mature_UA.old.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  ua0 = GenomicRanges::intersect(peaks[[30]], peaks[[27]])
  
  kk = c(26, 28, 29)
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll", main='BL_UA_day13_distal')
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Mature_UA.new.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  
  # for the moment, ua peaks overlapped by two replicates will be saved here 
  #ua = peaks[[29]][overlapsAny(peaks[[29]], peaks[[26]])|overlapsAny(peaks[[29]], peaks[[28]])]
  ua1 = intersect(peaks[[26]], peaks[[28]])
  ua2 = intersect(peaks[[26]], peaks[[29]])
  ua3 = intersect(peaks[[29]], peaks[[28]])
  
  ua = union(union(ua1, ua2), ua3)
  ua = reduce(ua)
  
  # LA
  kk = which(design$condition == 'Mature_LA')
  
  kk = c(23, 24, 25)
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll")
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Mature_LA.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  kk = c(25, 32, 34)
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll")
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Mature_LA_R12810_3.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  
  # for the moment, mLA peak consensus are peaks covered by > 2 out of three replicates
  la1 = intersect(peaks[[23]], peaks[[24]])
  la2 = intersect(peaks[[25]], peaks[[24]])
  la3 = intersect(peaks[[23]], peaks[[25]])
  
  la = union(union(la1, la2), la3)
  la = reduce(la)
  
  # for the moment, mHand peak consensus are peaks covered by > 2 out of three replicates
  kk = which(design$condition == 'Mature_Hand')
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll")
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Mature_Hand.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  hd1 = intersect(peaks[[20]], peaks[[21]])
  hd2 = intersect(peaks[[22]], peaks[[21]])
  hd3 = intersect(peaks[[20]], peaks[[22]])
  
  hd = union(union(hd1, hd2), hd3)
  hd = reduce(hd)
  
  # for the moment head control only one sample
  kk = which(design$condition == 'HEAD')
  ol.peaks <- makeVennDiagram(peaks[kk], NameOfPeaks=names(peaks)[kk], connectedPeaks="keepAll")
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  if(length(dev.list())!=0) dev.off()
  pdf(paste0(resDir, '/manualCheck_peakOverlapping_betweenReplicates_Mature_Head.pdf'),  height = 10, width = 10)
  try(plot(v)); dev.off()
  
  hc = peaks[[kk]]
  
  ## here  es40, es44.d, es44.p, bld5, bld9, bld13.p, bld13.d data are complete
  ## in contrast, ua (not sure), la, hd, ua, hc and hand regeneration time points will have new data
  ## when new data comes, the peak consensue will start from here
  save(ua, la, hd, hc, 
       file = paste0(RdataDir, '/consensus_peaks_intersectReplicates_pval', pval.cutoff, 'version_', version.analysis, 
                     'Mature.Rdata'))
  
  
  ##########################################
  # define the peak union across conditions
  ##########################################
  peak.merged = union(bld5, bld9)
  peak.merged = union(peak.merged, bld13.p)
  peak.merged = union(peak.merged, bld13.d)
  peak.merged = union(peak.merged, es40)
  peak.merged = union(peak.merged, es44.d)
  peak.merged = union(peak.merged, es44.p)
  peak.merged = union(peak.merged, ua)
  peak.merged = union(peak.merged, la)
  peak.merged = union(peak.merged, hd)
  
  #peak.merged = union(peak.merged, hc) # probalby the peak from the negative controls heads are not included
  peak.merged = dropSeqlevels(peak.merged, 'chrM', pruning.mode=c("coarse"))
  
  save(peak.merged, file = paste0(RdataDir, '/peaks_set_union_conditions_pval', pval.cutoff, 'version_', version.analysis, '.Rdata'))
  
  
  ##########################################
  # compare the peak overlapping in the UA regeneration
  ##########################################
  Compare.peak.overlapping.UA.regeneration = FALSE
  if(Compare.peak.overlapping.UA.regeneration){
    ol.peaks <- makeVennDiagram(list(ua, bld5), NameOfPeaks=c('mUA', 'BL.day5'), connectedPeaks="keepAll", main='BL_UA_day13_distal')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    ol.peaks <- makeVennDiagram(list(ua, bld9), NameOfPeaks=c('mUA', 'BL.day9'), connectedPeaks="keepAll", main='BL_UA_day13_distal')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    ol.peaks <- makeVennDiagram(list(ua, bld5, bld9), NameOfPeaks=c('mUA', 'BL.day5', 'BL.day9'), connectedPeaks="keepAll", main='BL_UA_day13_distal')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    ol.peaks <- makeVennDiagram(list(bld5, bld9), NameOfPeaks=c('BL.day5', 'BL.day9'), connectedPeaks="keepAll", main='BL_UA_day13_distal')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    ol.peaks <- makeVennDiagram(list(bld5, bld13.p), NameOfPeaks=c('BL.day5', 'BL.d13.p'), connectedPeaks="keepAll", main='BL_UA_day13_distal')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    ol.peaks <- makeVennDiagram(list(bld5, bld13.d), NameOfPeaks=c('BL.day5', 'BL.d13.d'), connectedPeaks="keepAll", main='BL_UA_day13_distal')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    pdfname = paste0(resDir, '/compare_peakOverlapping_betweenReplicates.pdf')
    pdf(pdfname, width = 16, height = 8)
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
    
    #Comparison.overlapping.peaks(design.matrix = design, peaks.list = peak.files, toCompare = 'condition', pval = 3, PLOT.p10 = FALSE)
    
    dev.off()
    
  }
 
}
  
# clean peaks

peaks = peak.merged
peaks = data.frame(peaks)
colnames(peaks)[c(1:3)] = c("chr", "start", "end")

dim(peaks)

peaks$peak.name = paste0(peaks$chr, ":", peaks$start, "_", peaks$end)

jj = match(unique(peaks$peak.name), peaks$peak.name)
peaks = peaks[jj, ];
dim(peaks)

peaks = peaks[which(peaks$chr != 'chrM'), ]

dim(peaks)


# preapre SAF input for featureCounts
require(Rsubread)

#counts = quantify.signals.within.peaks(peaks = peak.merged, bam.list=bam.list, rpkm.normalization = FALSE, isPairedEnd = FALSE)
SAF = data.frame(GeneID=peaks$peak.name, 
                 Chr=peaks$chr, 
                 Start=peaks$start, 
                 End=peaks$end, 
                 Strand=peaks$strand, stringsAsFactors = FALSE)

write.table(SAF, file = paste0(peakDir, 'merge_peak.saf'), sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE) 

Save.Peaklist = FALSE
if(Save.Peaklist){
  
  df = peaks
  
  bb = data.frame(df[ ,c(1, 2, 3, 6, 5)])
  
  write.table(bb, file = paste0(peakDir, 'merge_peak_usedBams_pval.6.bed'), sep = '\t', row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
  saveRDS(bb, file = paste0(RdataDir, '/merged_peaks_usedBams_pval.6.rds'))
  
  
  bb = readRDS(file = paste0(RdataDir, '/merged_peaks_usedBams_pval.6.rds'))
  
  bb = readRDS(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/
                             positional_memory/results/R10723_Rxxxx_atacseq_bowtie2.newParam_mtDNA_picardrmdup_20210208/
                             Rdata/merged_peaks_usedBams_pval.6.rds'))
  
  
  source('functions_chipSeq.R')
  peaks = annotatePeak.curateAxolotl(bb) 
  
  kk = which(!is.na(peaks$promoters))
  peaks.promoter = peaks[kk, c(1:5)]
  peaks.nonpromters = peaks[-kk, c(1:5)]
  
  write.table(peaks.promoter, file = paste0(peakDir, 'merge_peak_usedBams_pval.6_promoters.1000bpUp.200bpDown.bed'), 
              sep = '\t', row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
  write.table(peaks.nonpromters, file = paste0(peakDir, 'merge_peak_usedBams_pval.6_nonPromoters.1000bpUp.200bpDown.bed'), 
              sep = '\t', row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
}


##########################################
# run DESeq2 for QC and detect DE peaks
##########################################
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

#load(file = paste0(RdataDir, '/R11637_atacseq_samples_design_stats.Rdata'))
#design = stats
design = readRDS(file = paste0(RdataDir, '/design_merged_technicalReplicates_Rxxxx_R10723_R11637_R12810.rds'))
mm = match(design$sampleID, as.character(c(161523:161526)))
design = design[is.na(mm), ]

design = rbind(design, c('177595.177597', 'HEAD', 'HEAD_177595.177597', NA, NA))
design = rbind(design, c('177596.177598', 'Mature_LA', 'Mature_LA_177596.177598', NA, NA))

xlist<-list.files(path=paste0(dataDir, 'featurecounts_peaks.Q30'),
                  pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

colnames(design)[1] = 'SampleID'
#design = design[grep('90392|90393|108072|108073|108070|108071', design$SampleID, invert = TRUE),]

counts = process.countTable(all=all, design = design[, c(1,2)])

#design$conds = design$condition

#index = c()
# for(n in 1:nrow(design))
# {
#   index = c(index, grep(design$SampleID[n], stats$samples))
# }
# 
# design = data.frame(design, stats, stringsAsFactors = FALSE)

save(design, counts, file = paste0(RdataDir, 
                                   '/samplesDesign_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical_', 
                                   version.analysis, '.Rdata'))


Compare.with.Old.mature.resequencing = FALSE
if(Compare.with.Old.mature.resequencing){
  
  xlist = list.files(path=paste0(dataDir, 'featurecounts_peaks.Q30'),
                     pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge
  
  all = cat.countTable(xlist, countsfrom = 'featureCounts')
  
  design = colnames(all)[-1]
  design = gsub('_uniq_rmdup_featureCounts.txt', '', design)
  design = data.frame(sample = design, stringsAsFactors = FALSE)
  design$sampleID = sapply(design$sample, function(x) {test = unlist(strsplit(as.character(x), '_')); return(test[length(test)])})
  design$condition = sapply(design$sample, function(x) {
    test = unlist(strsplit(as.character(x), '_')); 
    return(paste0(test[-length(test)], collapse = '_'))})
  
  design = design[, c(2,3, 1)]
  
  colnames(design)[1] = 'SampleID'
  
  counts = process.countTable(all=all, design = design[, c(1,2)])
  
  design$conds = design$condition
  
  design0 = design
  counts0 = counts
  
  save(design0, counts0, file = paste0(RdataDir, '/samplesDesign_readCounts.withinPeaks.pval6_oldMature.Embryo.Regeneration.Rdata'))
  
}

##########################################
#  peak signal normalization
##########################################
load(file = paste0(RdataDir, 
                   '/samplesDesign_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical_', 
                   version.analysis, '.Rdata'))
design$sampleID = design$SampleID

design$batch = 'old'
design$batch[grep('1775', design$sampleID)] = '2022'
sels = grep('HEAD|Mature', design$condition)
design = design[sels, ]
counts = counts[, c(1, sels + 1)]

ss = apply(as.matrix(counts[, -1]), 1, mean)

par(mfrow=c(1,2))
hist(log10(as.matrix(counts[, -1])), breaks = 100, xlab = 'log10(nb of reads within peaks)', main = 'distribution')
plot(ecdf(log10(as.matrix(counts[, -1]) + 0.1)), xlab = 'log10(nb of reads within peaks)', main = 'cumulative distribution')

ss = apply(as.matrix(counts[, -1]), 2, sum)
design$usable.reads.withinPeaks = ss

design$pct.reads.in.peaks = design$usable.reads.withinPeaks/10^6/as.numeric(as.character(design$usable))

par(mfrow=c(1,1))

hist(design$pct.reads.in.peaks, main = 'distribution of pct of usable reads within merged peaks', xlab = 'pct', breaks = 10)

#norms = as.numeric(as.character(unlist(design$unique.rmdup)))
#norms = apply(counts[, -1], 2, sum)
#norms = norms/median(norms)

ss = apply(as.matrix(counts[, -1]), 1, max)

cutoff = 50
hist(log10(ss), breaks = 200)
abline(v = log10(cutoff), col = 'red')
kk = which(ss>cutoff)
length(which(ss>cutoff))

#pdfname = paste0(resDir, "/atacseq_normalization_assessment_allSamples.pdf")
#pdf(pdfname, width = 12, height = 10)

#Check.RNAseq.Quality(read.count=counts[kk, -1], design.matrix = data.frame(design$SampleID, design$conds), norms = norms)
require(ggplot2)
require(DESeq2)

rownames(counts) = counts$gene
dds <- DESeqDataSetFromMatrix(as.matrix(counts[kk, -1]), DataFrame(design), design = ~ condition)

#dds = dds[, grep('1361', design$SampleID)]
#dds$condition = droplevels(dds$condition)

#ss0 = rowMaxs(counts(dds))
#dds = dds[ss0 > 50, ]
#dds = dds[ss > cutoff, ]
ss = rowSums(counts(dds))
length(which(ss > quantile(ss, probs = 0.6)))

dd0 = dds[ss > quantile(ss, probs = 0.6) , ]
dd0 = estimateSizeFactors(dd0)
sizefactors.UQ = sizeFactors(dd0)

sizeFactors(dds) <- sizefactors.UQ
fpm = fpm(dds, robust = TRUE)

dds <- estimateDispersions(dds, fitType = 'parametric')
# plotDispEsts(dds, ymin = 10^-3, main = 'RA')
# 
# dds <- nbinomWaldTest(dds)
# resultsNames(dds)  
# 
# res1 = results(dds, contrast=c("condition", 'BL_UA_5days', 'Mature_UA'), alpha = 0.1)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = colnames(design)[2], ntop = 3000, returnData = FALSE)
print(pca)


pca2save = as.data.frame(plotPCA(vsd, intgroup = colnames(design)[c(2, 7)], returnData = TRUE, ntop = 3000))
pca2save$name = paste0(design$condition, '_', design$SampleID, '_', design$batch)
#pca2save$batch = 'old'
#pca2save$batch[grep('1361|1373', pca2save$name)] = 'new'
pca2save$batch = as.factor(pca2save$batch)

ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape = batch)) + 
  geom_point(size=3) + 
  geom_text(hjust = 0.2, nudge_y = 0.5, size=3)

plot(ggp) 
ggsave(paste0(resDir, "/PCA_allatacseq_new.old.merged_ntop3000_mature.pdf"), width = 16, height = 10)

plot(sizeFactors(dds), design$mapped, log = '')
plot(sizeFactors(dds), design$usable, log = 'xy')
text(sizeFactors(dds), design$usable, labels = design$fileName, cex = 0.7)

#plot(sizeFactors(dds), design$mapped, log = 'xy')
#text(sizeFactors(dds), design$mapped, labels = design$samples, cex = 0.7)

#dev.off()

save.scalingFactors.for.deeptools = FALSE
if(save.scalingFactors.for.deeptools){
  design$usable[which(design$sampleID == '177595.177597')] = 
    sum(as.numeric(as.character(design$usable[which(design$SampleID == '177595'| design$SampleID == '177597')])))
  
  design$usable[which(design$sampleID == '177596.177598')] = 
    sum(as.numeric(as.character(design$usable[which(design$SampleID == '177596'| design$SampleID == '177598')])))
  
  xx = data.frame(sampleID = design$sampleID,  
                  scalingFactor = as.numeric(design$usable)/(sizeFactors(dds)*median(as.numeric(design$usable))),
                  stringsAsFactors = FALSE)
  xx = xx[c(17:18), ]
  write.table(xx, file = paste0(dataDir, '/DESeq2_scalingFactor_forDeeptools.txt'), sep = '\t',
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
}
