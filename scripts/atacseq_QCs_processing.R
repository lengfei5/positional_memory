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

version.analysis = 'R10723_Rxxxx_R11637_atacseq_R11876_CutTag'
#peakDir = "Peaks/macs2_broad"

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

source('Functions_atac.R')

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11637_atac/'
#dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11876_cut.run/'

design_file = paste0(dataDir, 'design_sampleInfo_all.txt')

########################################################
########################################################
# Section : sequencing quality controls
# fragment size distribution 
# sequence saturation analysis
########################################################
########################################################
Collect.design.stat.nf.out = TRUE
if(Collect.design.stat.nf.out){
  design = read.table(paste0(dataDir, 'sampleInfo_parsed.txt'), sep = '\t', header = TRUE)
  
  stats = read.table(paste0(dataDir, 'nf_out/result/countStatTable.txt'), sep = '\t', header = TRUE)
  colnames(stats)[c(1, 3)] = c('fileName', 'trimmed')
  
  #cnts = list.files(path = '../Data//R10723_atac/QCs/cnt_raw', pattern = '*.txt', full.names = TRUE)
  
  #design = design[order(design$fileName), ]
  
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
    
    # total = 0
    # for(m in 1:length(cc.files))
    # {
    #   #cat(m, '\n')
    #   total = total + read.table(cc.files[m], sep = '\t', header = FALSE)
    # }
    # 
    # ss = read.table(files.stat[grep(design$sampleID[n], files.stat)], sep = '\t', header = TRUE)
    # samples = c(samples, as.character(ss[1, 1]))
    # stats = rbind(stats, c(total, ss[1, -1]))
    
  }
  
  
  stats = data.frame(design, stats[index, ], stringsAsFactors = FALSE)
  colnames(stats) = c('sampleID', 'samples', 'fileName', 'total',  'adapter.trimmed', 'mapped', 'chrM.rm', 'unique', 'unique.rmdup')
  
  stats$trimming.pct = as.numeric(stats$adapter.trimmed)/as.numeric(stats$total)
  stats$mapped.pct = as.numeric(stats$mapped)/as.numeric(stats$adapter.trimmed)
  stats$mito.pct = as.numeric(stats$chrM.rm)/as.numeric(stats$adapter.trimmed)
  stats$multimapper.pct = 1- as.numeric(stats$unique) / as.numeric(stats$mapped)
  stats$dup.rate = 1.0 - as.numeric(stats$unique.rmdup)/as.numeric(stats$unique)
  stats$pct.usable = stats$unique.rmdup / stats$total
  
  #stats$sample = gsub('_sorted', '', stats$sample)
  #stats = stats[, c(1, 2, 3, 6, 4, 7, 5, 8)]
  #colnames(stats)[c(5, 7)] = c('uniq.mapped', 'uniq.mapped.rmdup')
  # library(ggplot2)
  # library(dplyr)
  # 
  # xx = stats[, c(8:12)]
  # 
  # # Create data
  # data <- data.frame(
  #   name=c(rep(colnames(xx), each = nrow(xx))),
  #   value=as.numeric(unlist(xx)) %>% round(2)
  # )
  # 
  # 
  # ggplot(data, aes(x=name, y=value, fill=name)) + 
  #   geom_violin()
  # 
  # #boxplot(stats[, c(8:12)])
  # df <- apply(stats,2,as.character)
  #stats = as.data.frame(stats)
  
  write.csv(stats, file = paste0(resDir, '/R11876_CutTag_QCs_stats.csv'), row.names = FALSE)
  
  stats$usable = stats$unique.rmdup/10^6
  colnames(stats)[c(2,3)] = c('condition', 'samples')
  stats$samples = paste0(stats$condition, '_', stats$sampleID)
  
  plot(stats$usable, stats$mapped); text(stats$usable, stats$mapped, stats$samples)
  
  #save(stats, file = paste0(RdataDir, '/R11876_CutTag_samples_design_stats.Rdata'))
  save(stats, file = paste0(RdataDir, '/R11637_atacseq_samples_design_stats.Rdata'))
  
}else{
  
  extract.stat.for.samples.manual(design_file, dataDir)
  
}

##########################################
# fragment size distribution
##########################################
Make.fragsize.distribution = FALSE
if(Make.fragsize.distribution){
  files.insertion = list.files(path = '../Data/R10723_atac/QCs/frag_sizes', pattern = '*.txt', full.names = TRUE)
  
  design = design[order(design$fileName), ]
  
  conds = unique(design$fileName)
  
  
  pdfname = paste0(resDir, '/fragment_size_distribution.pdf')
  pdf(pdfname, width = 12, height = 10)
  
  for(n in 1:length(conds))
  {
    # n = 1
    cat(n, '--', as.character(conds[n]), '\n')
    
    jj = which(design$fileName == conds[n])
    
    par(mfrow=c(2,2))
    #if(length(jj) == 4) 
    #if(length(jj) == 2) par(mfrow = c(1, 2))
    
    for(m in jj){
      #cat(m, '\n')
      ff = files.insertion[grep(design$sampleID[m], files.insertion)]
      frag = read.delim(ff, sep = '\t', header = FALSE, comment.char = "#")[-c(1:3), c(1, 2)]
      
      xx = as.numeric(as.character(frag$V1))
      yy = as.numeric(as.character(frag$V2))
      plot(xx, yy, type = 'l', lwd =2.0,  col = 'red', main = paste0(design$fileName[m], '_',  design$sampleID[m]), 
           xlab = 'fragment size (bp)', ylab = 'counts of usable reads')
      
    }
    
  }
  dev.off()
  
}

########################################################
########################################################
# Section : saturation curve for sequencing depth 
# 
########################################################
########################################################
Sequence.Saturation.Analysis = FALSE
if(Sequence.Saturation.Analysis){
  library("ChIPseeker");
  library("rtracklayer")
  load(file = paste0(RdataDir, '/samples_design_stats.Rdata'))
  
  sampleUsed = 'mergedRep_downsample.Picard'
  Dir.downsampledBam = '../Data/R10723_atac/saturation/bams_downsampled_picard'
  Dir.called.peaks = '../Data/R10723_atac/saturation/calledPeaks_downsampledPicard/macs2'
  
  peak.files = list.files(path = Dir.called.peaks, 
                          pattern = '*macs2_peaks.xls', full.names = TRUE)
  total.files = list.files(path = Dir.downsampledBam , 
                           pattern = '*bam.counts.txt', full.names = TRUE)
  
  macs.peaks.stats = function(f)
  {
    p = readPeakFile(f, as = "GRanges")
    p10 <- p[mcols(p)[,"X.log10.pvalue."] > 10]
    return(c(length(p), sum(width(p)), length(p10), sum(width(p10))))
  }
  
  #test = sapply(peak.files[1:10], macs.peaks.stats)
  test = matrix(NA, ncol = 4, nrow = length(peak.files))
  for(n in 1:length(peak.files)){
    # n =1
    cat(n, '\n')
    #p = readPeakFile(peak.files[n], as = "GRanges")
    #peak.nbs = c(peak.nbs, length(p))
    #widths = c(widths, sum(width(p)))
    test[n, ] = macs.peaks.stats(peak.files[n])
  }
  
  colnames(test) = c('peaks.nb', 'peaks.width', 'peaks.nb.p10', 'peaks.width.p10')
  save(test, file = paste0(RdataDir, '/saturation_data_peak.nbs_widths_', sampleUsed, '.Rdata'))
  
  sat = data.frame(gsub('_macs2_peaks.xls', '', basename(peak.files)), test,  stringsAsFactors = FALSE)
  colnames(sat)[1] = c('samples')
  
  sat$peaks.nb = sat$peaks.nb/10^3
  sat$peaks.nb.p10 = sat$peaks.nb.p10/10^3
  sat$peaks.width = sat$peaks.width/10^6
  sat$peaks.width.p10 = sat$peaks.width.p10/10^6
  
  #pct.downsample = gsub('uniq_rmdup_downsampled.', '', sat$samples)
  
  sat$ss = NA
  sat$pcts = NA
  sat$reads = NA
  
  for(n in 1:nrow(sat))
  {
    # n = 20
    cat(n, '\n')
    sample.pcts = gsub('downsampled.', '', sat$samples[n])
    sample.pcts = unlist(strsplit(as.character(sample.pcts), '_'))
    
    sat$pcts[n] = as.numeric(sample.pcts[length(sample.pcts)])
    sat$ss[n] = paste0(sample.pcts[-length(sample.pcts)], collapse = "_")
    sat$reads[n] = as.numeric(read.table(total.files[grep(paste0(sat$samples[n], '.bam.counts'), total.files)])[1, 1])/10^6
    
  }
  
  save(sat, file = paste0(RdataDir, '/saturation_data_peak.nbs_widths_downsample.ptc_', sampleUsed, '.Rdata'))
  
  load(file = paste0(RdataDir, '/saturation_data_peak.nbs_widths_downsample.ptc_', sampleUsed, '.Rdata'))
  
  sample.uniq = unique(sat$ss)
  
  pdfname = paste0(resDir, '/saturation_curve_sequencing_depth_mergedReplicates_', sampleUsed, '.pdf')
  pdf(pdfname, width = 20, height = 16)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  par(mfrow=c(2,2))
  
  span = 0.5
  
  for(n in 1:length(sample.uniq))
  {
    # n = 1
    cat(n, '\n')
    
    # saturation curve with nb of peaks 
    # usable.reads = stats$unique.rmdup[which(stats$samples == sample.uniq[n])][[1]]
    kk = which(sat$ss == sample.uniq[n])
    xlims = c(0, max(105, max(sat$reads[kk])))
    
    plot(sat$reads[kk], sat$peaks.nb[kk], type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', 
         xlab = 'nb of usable reads (Million)', ylab = 'nb of peaks (K)', xlim = xlims)
    sat1 = data.frame(nb.reads = sat$reads[kk], nb.peaks = sat$peaks.nb[kk])
    loessMod <- loess(nb.peaks ~ nb.reads, data=sat1, span=span)
    smoothed <- predict(loessMod)
    lines(smoothed, x=sat1$nb.reads, col="red", lwd = 2.0)
    #points(sat$pcts[kk] * usable.reads/10^6, sat$nb.peaks[kk]/10^3, type = 'p', col = 'darkblue', pch =1, cex = 1.5)
    abline(v = 100, col = 'blue', lwd = 2.0)
    
    sat2 = data.frame(nb.reads = sat$reads[kk], peak.width = sat$peaks.width[kk])
    loessMod2 <- loess(peak.width ~ nb.reads, data=sat2, span=span)
    smoothed2 <- predict(loessMod2) 
    
    plot(sat2$nb.reads, sat2$peak.width, type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', 
         xlab = 'nb of usable reads (Million)', ylab = 'total peak width (M)', xlim = xlims)
    lines(smoothed2, x=sat2$nb.reads, col="red", lwd = 2.0)
    abline(v = 100, col = 'blue', lwd = 2.0)
    
    
    plot(sat$reads[kk], sat$peaks.nb.p10[kk], type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', 
         xlab = 'nb of usable reads (Million)', ylab = 'nb of peaks.p10 (K)', xlim = xlims)
    sat1 = data.frame(nb.reads = sat$reads[kk], nb.peaks = sat$peaks.nb.p10[kk])
    loessMod <- loess(nb.peaks ~ nb.reads, data=sat1, span=span)
    smoothed <- predict(loessMod)
    lines(smoothed, x=sat1$nb.reads, col="red", lwd = 2.0)
    #points(sat$pcts[kk] * usable.reads/10^6, sat$nb.peaks[kk]/10^3, type = 'p', col = 'darkblue', pch =1, cex = 1.5)
    abline(v = 100, col = 'blue', lwd = 2.0)
    
    sat2 = data.frame(nb.reads = sat$reads[kk], peak.width = sat$peaks.width.p10[kk])
    loessMod2 <- loess(peak.width ~ nb.reads, data=sat2, span=span)
    smoothed2 <- predict(loessMod2) 
    
    plot(sat2$nb.reads, sat2$peak.width, type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', 
         xlab = 'nb of usable reads (Million)', ylab = 'total peak.p10 width (M)', xlim = xlims)
    lines(smoothed2, x=sat2$nb.reads, col="red", lwd = 2.0)
    abline(v = 100, col = 'blue', lwd = 2.0)
    
    
  }
  
  dev.off()
  
  
  pdfname = paste0(resDir, '/saturation_curve_sequencing_depth_mergedReplicates_', sampleUsed, '_all_stringentPeaks.pdf')
  pdf(pdfname, width = 16, height = 8)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  sample.uniq = unique(sat$ss)
  span = 0.75
  # saturation curve with nb of peaks
  xlims = c(0, 300)
  ylims = range(sat$peaks.nb.p10)
  library(RColorBrewer)
  cols = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(length(sample.uniq))
  plot(0, 0, xlim = xlims, ylim = ylims, type ='n', xlab = 'nb of usable reads (Million)', 
       ylab = 'nb of peaks (K)', main = paste0('saturation curve of old and new samples - ', sampleUsed, '- p10 peaks'))
  abline(v = c(100, 150, 200), col = 'blue', lwd = 1.0, lty =2)
  
  #legend('topleft', legend = sample.uniq, col = cols, bty = 'n', lwd = 2.0, cex = 0.7)
  
  for(n in 1:length(sample.uniq))
  {
    # n = 1
    cat(n, '\n')
    kk = which(sat$ss == sample.uniq[n])
    satt = data.frame(nb.reads = sat$reads[kk], nb.peaks = sat$peaks.nb.p10[kk])
    
    points(satt[,1], satt[,2], type= 'p', col = cols[n])
    loessMod <- loess(nb.peaks ~ nb.reads, data=satt, span=span)
    smoothed <- predict(loessMod)
    text(satt[length(kk), 1], smoothed[length(smoothed)], labels = sample.uniq[n], cex = 0.7, pos = 4, offset = 0.2)
    lines(smoothed, x=satt$nb.reads, col=cols[n], lwd = 3.0)
    
  }
  
  dev.off()
  
}

########################################################
########################################################
# Section : Quantify peak signals
# 1) merge peaks
# 2) quantify the read counts within merged peaks
# 3) overview sample quality with PCA
########################################################
########################################################
peakDir = paste0(dataDir,  'calledPeaks/macs2')
peak.files = list.files(path = peakDir,
                        pattern = '*_peaks.xls', full.names = TRUE)
design = readRDS(file = paste0(RdataDir, '/design_merged_technicalReplicates.rds'))

mm = match(design$sampleID, as.character(c(161523:161526)))
design = design[is.na(mm), ]

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

# try to merge BL time series
kk = which(design$condition == 'BL_UA_13days_proximal')
bld13.p = peaks[[3]][overlapsAny(peaks[[3]], peaks[[4]])]

kk = which(design$condition == 'BL_UA_13days_distal')
bld13.d = peaks[[1]][overlapsAny(peaks[[1]], peaks[[2]])]

kk = which(design$condition == 'BL_UA_9days')
bld9 = peaks[[9]]
bld9 = bld9[overlapsAny(bld9, peaks[[10]])]

kk = which(design$condition == 'BL_UA_5days')
bld5 = peaks[[7]]
bld5 = bld5[overlapsAny(bld5, peaks[[8]])]
bld52 = peaks[[5]]
bld52 = bld52[overlapsAny(bld52, peaks[[6]])]

length(bld5[overlapsAny(bld5, bld52)])
bld5 = bld52



pdfname = paste0(resDir, '/compare_peakOverlapping_betweenReplicates.pdf')
pdf(pdfname, width = 16, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

source('functions_chipSeq.R')

#Comparison.overlapping.peaks(design.matrix = design, peaks.list = peak.files, toCompare = 'condition', pval = 3, PLOT.p10 = FALSE)

dev.off()

#peak.files = peak.files[grep('90392|90393|108072|108073|108070|108071', peak.files, invert = TRUE)]

#bam.list = list.files(path = '../Data/R10723_atac/alignments/BAMs_uniq_rmdup', pattern = '*.bam$', full.names = TRUE)
peak.merged = merge.peaks.macs2(peak.files, pcutoff = 6)

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
design = readRDS(file = paste0(RdataDir, '/design_merged_technicalReplicates.rds'))

xlist<-list.files(path=paste0(dataDir, '/featurecounts_peaks.Q30'),
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

save(design, counts, file = paste0(RdataDir, '/samplesDesign_readCounts.withinPeaks.pval6_mergedTechnical.Rdata'))

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
load( file = paste0(RdataDir, '/samplesDesign_readCounts.withinPeaks.pval6_mergedTechnical.Rdata'))
design$batch = 1
design = design[, c(1, 2, 6, 3:5)]


if(Compare.with.Old.mature.resequencing){
  
   load(file = paste0(RdataDir, '/samplesDesign_readCounts.withinPeaks.pval6_oldMature.Embryo.Regeneration.Rdata'))
   design = design[, c(1,2)]
   design$batch = 3
   
   colnames(counts)[-1] = paste0(design[, 2], '_', design[,1], '_', design[, 3])
   
   design0 = design0[, c(1, 2)]
   design0$batch = 2
   colnames(counts0)[-1] = paste0(design0[, 2], '_', design0[,1], '_', design0[, 3])
   
   mm = match(counts$gene, counts0$gene)
   counts = cbind(counts, counts0[mm, -1])
   
   colnames(design0) = colnames(design)
   design = rbind(design, design0)
    
   design$condition[grep('EMbryo_Stage40', design$condition)] = 'Embryo_Stage40'
   
   cc = unique(design$condition[which(design$batch == 3)])
   
   kk = which(!is.na(match(design$condition, cc)))
   
   xx = design[, c(1, 2)]
   xx = xx[match(unique(xx$SampleID), design$SampleID), ]
   colnames(xx) = c('sampleID', 'fileName')
   write.table(xx, file = paste0(dataDir, 'design_sampleInfo_all.txt'), sep = '\t', col.names = TRUE, row.names = FALSE,
               quote = FALSE)
   
      
}

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
hist(log10(ss), breaks = 200)

cutoff = 50
kk = which(ss>cutoff)
length(which(ss>cutoff))

#pdfname = paste0(resDir, "/atacseq_normalization_assessment_allSamples.pdf")
#pdf(pdfname, width = 12, height = 10)

#Check.RNAseq.Quality(read.count=counts[kk, -1], design.matrix = data.frame(design$SampleID, design$conds), norms = norms)
require(ggplot2)
require(DESeq2)

rownames(counts) = counts$gene
dds <- DESeqDataSetFromMatrix(as.matrix(counts[kk, -1]), DataFrame(design), design = ~ condition)

#dds = dds[ss > cutoff, ]
ss = rowSums(counts(dds))
length(which(ss > quantile(ss, probs = 0.65)))

dd0 = dds[ss > quantile(ss, probs = 0.65) , ]
dd0 = estimateSizeFactors(dd0)
sizefactors.UQ = sizeFactors(dd0)

sizeFactors(dds) <- sizefactors.UQ
fpm = fpm(dds, robust = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = colnames(design)[2], ntop = 3000, returnData = FALSE)
print(pca)


pca2save = as.data.frame(plotPCA(vsd, intgroup = colnames(design)[c(2, 3)], returnData = TRUE, ntop = 3000))
pca2save$name = paste0(design$condition, '_', design$SampleID, '_', design$batch)
#pca2save$batch = 'old'
#pca2save$batch[grep('1361|1373', pca2save$name)] = 'new'
pca2save$batch = as.factor(pca2save$batch)

ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape = batch)) + 
  geom_point(size=3) + 
  geom_text(hjust = 0.2, nudge_y = 0.5, size=3)

plot(ggp) + ggsave(paste0(resDir, "/PCA_allatacseq_new.old.merged_ntop3000.pdf"), width = 16, height = 10)

plot(sizeFactors(dds), design$mapped, log = '')
plot(sizeFactors(dds), design$usable, log = 'xy')
text(sizeFactors(dds), design$usable, labels = design$fileName, cex = 0.7)

#plot(sizeFactors(dds), design$mapped, log = 'xy')
#text(sizeFactors(dds), design$mapped, labels = design$samples, cex = 0.7)

#dev.off()

save.scalingFactors.for.deeptools = FALSE
if(save.scalingFactors.for.deeptools){
  xx = data.frame(sampleID = design$SampleID,  
                  scalingFactor = design$usable/(sizeFactors(dds)*median(design$usable)),
                  stringsAsFactors = FALSE)
  
  write.table(xx, file = paste0(dataDir, '/DESeq2_scalingFactor_forDeeptools.txt'), sep = '\t',
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

##########################################
# cost estimate for Batch 1 
# experiment design by Elly (Feb 22th 2021)
# Batch 0 (done): mature.UA, Blastema.UA.day5, Blastema.UA.day9, Blastema.UA.day13.proxomial, Blastema.UA.day13.distal 
# Batch 1:  mature.UA,  mature.LA,  mature.Hand, head (control)
# Batch 2:  mature.Hand, Blastema.Hand.day5, Blastema.Hand.day9, Blastema.Hand.day13 
##########################################
Cost.estimation = FALSE
if(Cost.estimation){
  
  xx = design
  xx = xx[grep('Mature', xx$fileName), ]
  xx = xx[xx$batch == '2020', ]
  
  xx$nb.reads.to.add = (90*10^6 - xx$unique.rmdup)/xx$pct.usable/10^6
  
  sum(xx$nb.reads.to.add)
}

