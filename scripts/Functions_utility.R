##########################################################################
##########################################################################
# Project: Utility functions for this projects
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Oct  5 12:29:45 2021
##########################################################################
##########################################################################
##########################################
# plot fragment size distribution
##########################################
Make.fragsize.distribution = function()
{
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
# make the saturation curve for sequencing depth 
# 
########################################################
########################################################
Sequence.Saturation.Analysis = function(design)
{
  library("ChIPseeker");
  library("rtracklayer")
  
  #load(file = paste0(RdataDir, '/samples_design_stats.Rdata'))
  
  sampleUsed = 'mergedRep_downsample.Picard'
  
  Dir.downsampledBam = paste0(dataDir, 'R12810_cuttag/bam_downsampling/QCs/BamStat')
  Dir.called.peaks = paste0(dataDir, 'R12810_cuttag/bam_downsampling/calledPeaks/macs2_broad')
  
  peak.files = list.files(path = Dir.called.peaks, 
                          pattern = '*peaks.xls', full.names = TRUE)
  total.files = list.files(path = Dir.downsampledBam, 
                           pattern = '*_stat.txt', full.names = TRUE)
  
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
  
  sat = data.frame(gsub('_macs2_broad_fdr_0.05_peaks.xls', '', basename(peak.files)), test,  stringsAsFactors = FALSE)
  colnames(sat)[1] = c('samples')
  
  sat$peaks.nb = sat$peaks.nb/10^3
  sat$peaks.nb.p10 = sat$peaks.nb.p10/10^3
  sat$peaks.width = sat$peaks.width/10^6
  sat$peaks.width.p10 = sat$peaks.width.p10/10^6
  
  #pct.downsample = gsub('uniq_rmdup_downsampled.', '', sat$samples)
  
  sat$ss = NA
  sat$downsample = NA
  sat$total = NA
  sat$usable = NA
  sat$dup = NA
  
  for(n in 1:nrow(sat))
  {
    # n = 1
    cat(n, '\n')
    sample.pcts = gsub('downsampled.', '', sat$samples[n])
    sample.pcts = unlist(strsplit(as.character(sample.pcts), '_'))
    
    sat$downsample[n] = as.numeric(sample.pcts[length(sample.pcts)])
    sat$ss[n] = paste0(sample.pcts[-length(sample.pcts)], collapse = "_")
    ss = read.table(total.files[grep(sat$samples[n], total.files)], sep = '\t', header = TRUE)
    sat$total[n] = as.numeric(ss$total)/10^6
    sat$usable[n] = as.numeric(ss$rmdup_uniq)/10^6
    sat$dup[n] = 1 - as.numeric(ss$rmdup_uniq)/as.numeric(ss$uniq)
    
  }
  
  save(sat, file = paste0(RdataDir, '/saturation_data_peak.nbs_widths_downsample.ptc_', sampleUsed, '.Rdata'))
  
  write.csv(sat, file = paste0(resDir, '/downsampling_states.csv'))
  
  load(file = paste0(RdataDir, '/saturation_data_peak.nbs_widths_downsample.ptc_', sampleUsed, '.Rdata'))
  
  sample.uniq = unique(sat$ss)
  
  pdfname = paste0(resDir, '/saturation_curve_sequencing_depth_mergedReplicates_', sampleUsed, '.pdf')
  pdf(pdfname, width = 12, height = 10)
  #par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  par(mfrow=c(2,2))
  
  span = 0.5
  
  for(n in 1:length(sample.uniq))
  {
    # n = 1
    cat(n, '\n')
    
    # saturation curve with nb of peaks 
    # usable.reads = stats$unique.rmdup[which(stats$samples == sample.uniq[n])][[1]]
    kk = which(sat$ss == sample.uniq[n])
    xlims = c(0, max(105, max(sat$total[kk])))
    
    plot(sat$total[kk], sat$peaks.nb[kk], type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', 
         xlab = 'nb of total reads (Million)', ylab = 'nb of peaks (K)', xlim = xlims, cex = 1.2)
    sat1 = data.frame(nb.reads = sat$total[kk], nb.peaks = sat$peaks.nb[kk])
    loessMod <- loess(nb.peaks ~ nb.reads, data=sat1, span=span)
    smoothed <- predict(loessMod)
    lines(smoothed, x=sat1$nb.reads, col="red", lwd = 2.0)
    #points(sat$pcts[kk] * usable.reads/10^6, sat$nb.peaks[kk]/10^3, type = 'p', col = 'darkblue', pch =1, cex = 1.5)
    abline(v = 100, col = 'blue', lwd = 2.0)
    
    plot(sat$total[kk], sat$dup[kk], type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '',  cex = 1.2, 
         xlab = 'nb of total reads (Million)', ylab = 'dupliocation rates')
    
    plot(sat$usable[kk], sat$peaks.nb[kk], type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', cex = 1.2,
         xlab = 'nb of usable reads (Million)', ylab = 'nb of peaks (K)')
    
    
    plot(sat$total[kk], sat$dup[kk], type= 'n', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '',
         xlab = 'nb of total reads (Million)', ylab = 'dupliocation rates', cex = 1.2)
    
    # sat2 = data.frame(nb.reads = sat$reads[kk], peak.width = sat$peaks.width[kk])
    # loessMod2 <- loess(peak.width ~ nb.reads, data=sat2, span=span)
    # smoothed2 <- predict(loessMod2) 
    # 
    # plot(sat2$nb.reads, sat2$peak.width, type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', 
    #      xlab = 'nb of usable reads (Million)', ylab = 'total peak width (M)', xlim = xlims)
    # lines(smoothed2, x=sat2$nb.reads, col="red", lwd = 2.0)
    # abline(v = 100, col = 'blue', lwd = 2.0)
    
    # 
    # plot(sat$reads[kk], sat$peaks.nb.p10[kk], type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', 
    #      xlab = 'nb of usable reads (Million)', ylab = 'nb of peaks.p10 (K)', xlim = xlims)
    # sat1 = data.frame(nb.reads = sat$reads[kk], nb.peaks = sat$peaks.nb.p10[kk])
    # loessMod <- loess(nb.peaks ~ nb.reads, data=sat1, span=span)
    # smoothed <- predict(loessMod)
    # lines(smoothed, x=sat1$nb.reads, col="red", lwd = 2.0)
    # #points(sat$pcts[kk] * usable.reads/10^6, sat$nb.peaks[kk]/10^3, type = 'p', col = 'darkblue', pch =1, cex = 1.5)
    # abline(v = 100, col = 'blue', lwd = 2.0)
    # 
    # sat2 = data.frame(nb.reads = sat$reads[kk], peak.width = sat$peaks.width.p10[kk])
    # loessMod2 <- loess(peak.width ~ nb.reads, data=sat2, span=span)
    # smoothed2 <- predict(loessMod2) 
    # 
    # plot(sat2$nb.reads, sat2$peak.width, type= 'p', col = 'blue', lwd = 2.0,  main = sample.uniq[n], log = '', 
    #      xlab = 'nb of usable reads (Million)', ylab = 'total peak.p10 width (M)', xlim = xlims)
    # lines(smoothed2, x=sat2$nb.reads, col="red", lwd = 2.0)
    # abline(v = 100, col = 'blue', lwd = 2.0)
    # 
    
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

##########################################
# cost estimate for Batch 1 
# experiment design by Elly (Feb 22th 2021)
# Batch 0 (done): mature.UA, Blastema.UA.day5, Blastema.UA.day9, Blastema.UA.day13.proxomial, Blastema.UA.day13.distal 
# Batch 1:  mature.UA,  mature.LA,  mature.Hand, head (control)
# Batch 2:  mature.Hand, Blastema.Hand.day5, Blastema.Hand.day9, Blastema.Hand.day13 
##########################################
Cost.estimation = function()
{
  
  xx = design
  xx = xx[grep('Mature', xx$fileName), ]
  xx = xx[xx$batch == '2020', ]
  
  xx$nb.reads.to.add = (90*10^6 - xx$unique.rmdup)/xx$pct.usable/10^6
  
  sum(xx$nb.reads.to.add)
}


########################################################
########################################################
# Section : test if CT (cut & tag) controls are similar or not
# Is it necessary to make control samples each time point and each condition
########################################################
########################################################
Compare.CutandRun.controls = function()
{
  rm(list=ls())
  
  version.analysis = 'R10723_Rxxxx_R11637_atacseq_R11876_CutTag'
  #peakDir = "Peaks/macs2_broad"
  
  resDir = paste0("../results/", version.analysis)
  RdataDir = paste0(resDir, '/Rdata')
  if(!dir.exists(resDir)) dir.create(resDir)
  if(!dir.exists(RdataDir)) dir.create(RdataDir)
  
  source('Functions_atac.R')
  
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11876_cut.run/'
  
  design_file = paste0(dataDir, 'sampleInfo_parsed.txt')
  design = read.table(paste0(dataDir, 'sampleInfo_parsed.txt'), sep = '\t', header = TRUE)
  
  ##########################################
  # collect stats from QCs 
  ##########################################
  stats = read.table(paste0(dataDir, 'nf_out/result/countStatTable.txt'), sep = '\t', header = TRUE)
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
  
  design = stats
  
  ##########################################
  # compare peak overlapping
  ##########################################
  source('functions_chipSeq.R')
  peakDir = paste0(dataDir,  'nf_out/peaks_macs2')
  peak.files = list.files(path = peakDir,
                          pattern = '*_peaks.xls', full.names = TRUE)
  
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
  
  peak.merged = merge.peaks.macs2(peak.files, pcutoff = 6)
  # peaks = c()
  # pval.cutoff = 
  # for(n in 1:length(peak.files)) 
  # {
  #   cat(n, '\n')
  #   p = readPeakFile(peak.files[n], as = "GRanges");
  #   #eval(parse(text = paste0("p = pp.", k)));
  #   with.p.values = "X.log10.pvalue." %in% colnames(mcols(p))
  #   if(with.p.values) {
  #     p <- p[mcols(p)[,"X.log10.pvalue."] > pval.cutoff];
  #     p = reduce(p);
  #     #peaks10= c(peaks10, p10);
  #   }else{ 
  #     cat("no p values conlumn found for -- ", design.matrix$file.name[k], "\n");
  #     PLOT.p10 = FALSE;
  #   }
  #   #p = reduce(p)
  #   peaks= c(peaks, p)
  # }
  # 
  # design$fileName = paste0(design$samples, '_', design$sampleID)
  # names(peaks) = design$fileName
  
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
  
  ##########################################
  # compare the normalized peak signals 
  ##########################################
  RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
  RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
  source(RNA.functions)
  source(RNA.QC.functions)
  
  xlist<-list.files(path=paste0(dataDir, 'nf_out/featurecounts_peaks.Q30'),
                    pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge
  
  all = cat.countTable(xlist, countsfrom = 'featureCounts')
  
  colnames(design)[1] = 'SampleID'
    
  counts = process.countTable(all=all, design = design[, c(1,2)])
  
  save(design, counts, file = paste0(RdataDir, '/samplesDesign_readCounts.within_CutAndRUn.Rdata'))
  
  ##
  load(file = paste0(RdataDir, '/samplesDesign_readCounts.within_CutAndRUn.Rdata'))
  
  require(ggplot2)
  require(DESeq2)
  
  rownames(counts) = counts$gene
  colnames(design)[2] = 'condition'
  dds <- DESeqDataSetFromMatrix(as.matrix(counts[, -1]), DataFrame(design), design = ~ condition)
  
  jj = grep('IgG', design$condition)
  dds = dds[, jj]
  
  #design = design[kk, ]
  #counts = counts[]
  
  ss = rowMaxs(counts(dds))
  
  cutoff = 10
  hist(log10(ss), breaks = 200)
  abline(v = log10(cutoff), col = 'red')
  length(which(ss>cutoff))
  
  kk = which(ss>cutoff)
  
  dds = dds[kk,]
  fpm = fpm(dds, robust = FALSE)
  
  RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
  RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
  source(RNA.functions)
  source(RNA.QC.functions)
  
  plot.pair.comparison.plot(fpm, linear.scale = TRUE) 
  
}




