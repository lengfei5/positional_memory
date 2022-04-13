##########################################################################
##########################################################################
# Project: 
# Script purpose: functions of histone marker analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Apr  4 09:51:22 2022
##########################################################################
##########################################################################

##########################################
# utility functions 
##########################################
cal_transform_histM = function(x, cutoff.min = 1., cutoff.max = 5, toScale = FALSE, centering = FALSE)
{
  # x = keep[,1];cutoff.min = 3.5; cutoff.max = 6
  x[which(x<cutoff.min)] = cutoff.min
  x[which(x>cutoff.max)] = cutoff.max
  if(centering) x = x - mean(x)
  if(toScale) x = (x - cutoff.min)/(cutoff.max - cutoff.min)
  return(x)
  
}

cal_centering <- function(x){
  (x - mean(x, na.rm =TRUE))
}

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

## quantile normalization for different histone markers
Quantile.normalize.histMarker = function(shm)
{
  library(preprocessCore)
  xx = normalize.quantiles(shm)
  colnames(xx) = colnames(shm)
  rownames(xx) = rownames(shm)
  shm = xx
  
  return(shm)
  
}

cal_sample_means = function(cpm, conds = c("mUA", "mLA", "mHand") )
{
  sample.means = c()
  for(n in 1:length(conds)) 
  {
    kk = grep(conds[n], colnames(cpm))
    if(length(kk)>1) {
      sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
    }else{
      sample.means = cbind(sample.means, cpm[, kk])
    }
    
  }
  colnames(sample.means) = conds
  
  return(sample.means)
  
}


##########################################
# try to subtract the input IgG 
##########################################
Manually_define_histM_peakConsensus = function()
{
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
  
  
}


Subtract.IgG.inputs_viaGlobalScaling = function()
{
  load(file = paste0(RdataDir, 
                     '/histoneMarkers_samplesDesign_ddsPeaksMatrix_filtered_incl.atacPeak.missedTSS.bgs_324k.Rdata'))
  igg = readRDS(file = paste0(RdataDir, '/fpm_bc_TMM_combat_', 'IgG', '_', version.analysis, '.rds'))
  
  conds_histM = unique(design$sample)
  
  ctl.means = c()
  for(ii in 1:length(conds_histM)) 
  {
    kk = grep(conds_histM[ii], colnames(igg))
    if(length(kk)>1) {
      ctl.means = cbind(ctl.means, apply(igg[, kk], 1, median))
    }else{
      ctl.means = cbind(ctl.means, igg[, kk])
    }
  }
  
  colnames(ctl.means) = conds_histM
  
  ### Reason : reduce the ctl.mean with a fixed global factor and then fix the IgG signals
  ### scale each marker data based on the common IgG background 
  ctl.means = ctl.means - 3.5
  igg_peak = 4.5 - 3.5 # cutoff for peaks in IgG
  #plot.pair.comparison.plot(ctl.means[c(1:20000), ], linear.scale = FALSE)
  
  markers = c('H3K4me3', 'H3K4me1', 'H3K27me3', 'H3K27ac')
  
  for(n in 1:length(markers))
  {
    # n = 1
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
    abline(v = c(2.5,  3., 4), col = 'red')
    
    igg_cutoff = 3 # global factor for IgG
    cat(length(which(inputs.means> igg_cutoff)), ' peaks considered solely due to IgG \n')
    
    global_scaling = median(ratios[which(inputs.means > igg_cutoff)])
    cat('global scaling factor -- ', global_scaling, '\n')
    abline(h = global_scaling, col = 'blue')
    
    #bgs = inputs.means - global_scaling
    #igg_peak_new = igg_peak - median(ratios[which(inputs.means > igg_cutoff)])
    marker.means_new = marker.means + global_scaling
    ratios_new = inputs.means - marker.means_new
    cat('after sacling the markers ratios becomes : ', median(ratios_new[which(inputs.means > igg_cutoff)]), '\n')
    
    hist(marker.means_new, breaks = 100)
    hist(inputs.means, breaks = 100, add = TRUE, col = 'red')
    abline(v = igg_peak, col = 'blue')
    
    cpm = cpm + global_scaling
    
    for(m in 1:length(conds_histM))
    {
      # m = 1
      cat('condition -- ', conds_histM[m], '\n')
      jj = grep(conds_histM[m], colnames(cpm))
      
      bgs = inputs[, which(colnames(inputs) == conds_histM[m])]
      # bgs = bgs - global_scaling
      ss_mean = apply(cpm[,jj], 1, mean)
      j_cor = which(bgs > igg_peak & ss_mean > bgs + 1)
      
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
}


plot_individual_histMarker_withinATACpeak = function(res)
{
  # res = res[order(-res$log2fc), ]
  # sample.means = sample.means[match(rownames(res), rownames(sample.means)), ]
  # xx = res[select, ]

  # ####### focus on the atac-seq peak sets
  # Keep.only.overlapped.by.atacseq.peaks = TRUE
  # if(Keep.only.overlapped.by.atacseq.peaks){
  #   ii_bgs = grep('tss.', rownames(res))
  #
  #   rownames(res) = gsub('tss.', '', rownames(res))
  #
  #   pp = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  #
  #   rownames(res) = paste0(pp[,1], ':', pp[,2], '-', pp[,3])
  #   rownames(pp) = rownames(res)
  #
  #   pp$strand = '*'
  #   pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
  #                                 start.field="X2", end.field="X3", strand.field="strand")
  #   lls = width(pp)
  #
  #   #for(n in 1:ncol(keep))
  #   #{
  #   #  keep[,n] = keep[,n] + log2(1000/lls)
  #   #}
  #   #cpm_bgs = keep[ii_bgs, ]
  #   #keep = keep[-ii_bgs, ]
  #   #pp = pp[-ii_bgs]
  #
  #   ii_overlap = which(overlapsAny(pp, atacseq_peaks) == TRUE)
  #   #cpm_nonoverlap = keep[-ii_overlap, ]
  #   #keep = keep[ii_overlap, ]
  #   kk = match(rownames(res), names(pp)[ii_overlap])
  #
  #   res = res[!is.na(kk), ]
  #
  #   sample.means = sample.means[!is.na(kk), ]
  #
  #   fdr.cutoff = 0.05; logfc.cutoff = 1;
  #   marker.cutoff = 2;
  #   select = which(((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
  #                    (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
  #                    (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)) &
  #                    res$maxs > marker.cutoff
  #   )
  #   cat(length(select), ' DE ', conds_histM[n_histM],  ' \n')
  #
  #   xx = res[select, ]
  #   xx = xx[order(-xx$log2fc), ]
  #
  # }
  #
  # ii_bgs = grep('tss', rownames(cpm))
  # hist(cpm[-ii_bgs, 1], breaks = 100)
  # hist(cpm[ii_bgs, 1], breaks = 50, col = 'red', add = TRUE)
  #
  #
  # yy = res[select, grep(conds_histM[n_histM], colnames(res))]
  # #yy = t(apply(yy, 1, cal_transform_histM, cutoff.min = 0.0, cutoff.max = 6))
  #
  # df = as.data.frame(sapply(colnames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  # colnames(df) = 'segments'
  # rownames(df) = colnames(yy)
  #
  # sample_colors = c('springgreen4', 'steelblue2', 'gold2')
  # annot_colors = list(segments = sample_colors)

  # pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
  #          color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8),
  #          show_colnames = FALSE,
  #          scale = 'row',
  #          cluster_cols=FALSE, annotation_col=df,
  #          #annotation_colors = annot_colors,
  #          width = 6, height = 12,
  #          filename = paste0(figureDir, '/heatmap_histoneMarker_', conds_histM[n_histM], '_overlappedWithAtacseqPeak.pdf'))
  #
  # pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
  #          color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8),
  #          show_colnames = FALSE,
  #          scale = 'none',
  #          cluster_cols=FALSE, annotation_col=df,
  #          #annotation_colors = annot_colors,
  #          width = 6, height = 12,
  #          filename = paste0(figureDir, '/heatmap_histoneMarker_', conds_histM[n_histM], '_overlappedWithAtacseqPeak_nonscaled.pdf'))

}


##########################################
# compare MACS2 and SEACR peak caller 
##########################################
Compare.MACS2.vs.SEACR = function()
{
  library('rtracklayer')
  library(GenomicRanges)
  library('GenomicFeatures')
  
  # MACS2 differnet thresholds (10^-3, 10^-4, 10^-5, 10^-6)
  # SEACR relaxed, stringent, top0.01
  
  ## process Akane's CENs and putative enhancers
  enhancers = read.delim(file = '../data/HoxAD_enhancer_fromAkane.bed.txt', header = FALSE)
  enhancers = data.frame(enhancers, stringsAsFactors = FALSE)
  enhancers$score = 100
  enhancers$strand = '*'
  
  write.table(enhancers, file = paste0('../data/HoxAD_enhancer_fromAkane.bed'), col.names = FALSE, row.names = FALSE, quote = FALSE,
              sep = '\t')
  
  
  ## different stringencies of macs2
  peakDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/calledPeaks/macs2'
  peak.files = list.files(path = peakDir,
                          pattern = '*_peaks.xls', full.names = TRUE)
  peak.files = peak.files[grep('185746', peak.files)]
  
  p = readPeakFile(peak.files, as = "GRanges");
  #eval(parse(text = paste0("p = pp.", k)));
  with.p.values = "X.log10.pvalue." %in% colnames(mcols(p))
  
  for(cutoff in c(3, 4, 5, 6))
  {
    cat('pval cutoff -- ', cutoff, '\n')
    p_sel <- p[mcols(p)[,"X.log10.pvalue."] > cutoff];
    p_sel = GenomicRanges::reduce(p_sel);
    
    export(p_sel, 
           con = paste0('/Volumes//groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/MACS2_vs_SEACR/MACS2/', 
                        gsub('_macs2_peaks.xls', '', basename(peak.files)),  'macs2_pval_', cutoff, '.bed'), format = 'bed')
    
  }
  
  p3 = p
  p6 = p_sel
  p4 <- p[mcols(p)[,"X.log10.pvalue."] > 4];
  p5 = p[mcols(p)[,"X.log10.pvalue."] > 5];
  
  seacrDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/MACS2_vs_SEACR/SEACR/'
  
  stringent = readPeakFile(paste0(seacrDir, 'seacr_H3K4me3_BL5days_rRep1_185746_withIgG_strignet.stringent.bed'))
  relaxed = readPeakFile(paste0(seacrDir, 'seacr_H3K4me3_BL5days_rRep1_185746_withIgG_relaxed.relaxed.bed'))
  tops = readPeakFile(paste0(seacrDir, 'seacr_H3K4me3_BL5days_rRep1_185746_nonIgG_top0.01.stringent.bed'))
  
  counts = c(length(p3), length(p4), length(p5), length(p6), length(stringent), length(relaxed), length(tops))
  names(counts) = c('macs2.p3', 'macs2.p4', 'macs2.p5', 'macs2.p6', 'seacr.stringent', 'seacr.relaxed', 'seacr.top0.01')
  
  counts = data.frame(counts, names = names(counts), stringsAsFactors = FALSE)
  ggplot(data=counts, aes(x=names, y=counts)) +
    geom_bar(stat="identity", fill="steelblue")+
    geom_text(aes(label=counts), vjust=-0.3, size=3.5)+
    theme_minimal()
  ggsave(paste0(resDir, "/peak_numbers_macs2_seacr.pdf"), width=8, height = 6)
  
  
  ol.peaks <- makeVennDiagram(list(p6, stringent, relaxed), 
                              NameOfPeaks=c('macs_p6', 'seacr_stringent_IgG', 'seacr_relaxed_IgG'), connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/peakOverlapping_macs2.p6_seacr.IgG.stringent.relaxed.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  
  ol.peaks <- makeVennDiagram(list(p4, p6, tops), 
                              NameOfPeaks=c('macs_p4', 'macs_p6', 'seacr_top0.01'), connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/peakOverlapping_macs2.p3.p4_seacr.top0.01.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  ol.peaks <- makeVennDiagram(list(p6, tops), 
                              NameOfPeaks=c('macs_p6', 'seacr_top0.01'), connectedPeaks="keepAll", main=cc)
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  pdf(paste0(resDir, '/peakOverlapping_macs2.p6_seacr.top0.01.pdf'), 
      height = 10, width = 10)
  try(plot(v))
  dev.off()
  
  ss = tops[!overlapsAny(tops, p6)]
  ms = p6[!overlapsAny(p6, tops)]
  
  xx = data.frame(width = c(width(ms), width(ss)), peakCaller  = c(rep('macs2.p6', length(ms)), rep('seacr.top0.01', length(ss))))
  xx$width = log10(xx$width)
  
  ggplot(xx, aes(x=width, color=peakCaller)) +
    geom_histogram(fill="white", bins = 100) +
    scale_color_brewer(palette="Dark2") +
    labs(x = ' peak width log10 bp')
  ggsave(paste0(resDir, "/peak_width_distribution_macs2_seacr_specificPeaks.pdf"), width=8, height = 6)
    
  gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
  amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
  
  pdfname = paste0(resDir, "/feature_distribution_MACS2_peaks.pdf")
  pdf(pdfname, width = 6, height = 4)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  pp.annots = annotatePeak(ms, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  
  plotPeakAnnot_piechart(pp.annots)
  
  dev.off()
  
  pdfname = paste0(resDir, "/feature_distribution_SEACR_peaks.pdf")
  pdf(pdfname, width = 6, height = 4)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  #gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
  #amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
  pp.annots = annotatePeak(ss, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  
  plotPeakAnnot_piechart(pp.annots)
  
  dev.off()
  
  
  export(ms, 
  con = paste0('/Volumes//groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/MACS2_vs_SEACR/',
               'macs2_p6_detected.peaks.bed'), format = 'bed')
  
  export(ss, 
         con = paste0('/Volumes//groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/MACS2_vs_SEACR/',
                      'seacr_top0.01_detected.peaks.bed'), format = 'bed')
  
  
}

