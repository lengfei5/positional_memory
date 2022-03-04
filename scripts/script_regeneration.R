##########################################################################
##########################################################################
# Project: Akane's positional memory project
# Script purpose: regeneration time points analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Feb 22 10:46:43 2022
##########################################################################
##########################################################################
########################################################
########################################################
# Section IV : regeneration peaks
# or temporal-peaks test
########################################################
########################################################
##########################################
grouping.temporal.peaks = FALSE
if(grouping.temporal.peaks){
  
  fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration_embryoStages.rds'))
  design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration_embryoStages.rds'))
  
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
    
    save.peak.bed = FALSE
    if(save.peak.bed){
      bed = data.frame(pp[, c(1:3)], rownames(pp), 0, pp[, 4], stringsAsFactors = FALSE)
      write.table(bed, file = paste0(resDir, '/peakset_', version.analysis, '.bed'), col.names = FALSE, row.names = FALSE,
                  sep = '\t', quote = FALSE)
    }
    
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
  
  
  conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal")
  # examples to test
  test.examples = c('HAND2', 'FGF8', 'KLF4', 'Gli3', 'Grem1')
  #test.examples = c('Hoxa13')
  ii.test = which(overlapsAny(pp, promoters[which(!is.na(match(promoters$geneSymbol, test.examples)))]))
  #ii.Hox = which(overlapsAny(pp, Hoxs))
  #ii.test = unique(c(ii.test, ii.Hox))
  
  sample.sels = c(); cc = c()
  for(n in 1:length(conds)) {
    kk = which(design$conds == conds[n])
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
  }
  
  
  Run.temporal.peak.test = FALSE
  if(Run.temporal.peak.test)
  {
    source('Functions_atac.R')
    cpm = fpm[, sample.sels]
    
    # stringent here, with signal > 2.5 in >=3 samples 
    nb.above.threshold = apply(as.matrix(cpm), 1, function(x) length(which(x> 4.5)))
    hist(nb.above.threshold, breaks = c(-1:ncol(cpm)))
    peak.sels = which(nb.above.threshold>=2)
    cat(length(peak.sels), 'peaks after filtering for mature samples\n')
    
    cpm = cpm[peak.sels, ]
    
    tic()
    ## define the dynamic enhancers with mature UA and BL.UA and check them if embryo samples
    sels = grep('Embryo', cc, invert = TRUE) 
    res = t(apply(cpm[, sels], 1, temporal.peaks.test, c = cc[sels]))
    toc()
    
    xx = data.frame(res, pp.annots[match(rownames(cpm), rownames(pp.annots)), ],  stringsAsFactors = FALSE)
    
    res = xx
    
    saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test_v5.rds'))
    
  }
  
  res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test_v5.rds'))
  res = res[order(-res$log2FC), ]
  
  # select the temporal dynamic peaks
  length(which(res$prob.M0<0.05))
  length(which(res$prob.M0<0.05 & res$log2FC > 1))
  length(which(res$prob.M0<0.01 & res$log2FC > 1))
  length(which(res$prob.M0<0.01 & res$log2FC > 1.5))
  length(which(res$prob.M0<0.01 & res$log2FC > 2))
  
  #length(which(res$prob.M0<0.001 & res$log2FC > 2))
  
  #jj = which(res$prob.M0 < 0.05 & res$log2FC >1 )
  
  jj = which(res$prob.M0 < 0.05 & res$log2FC > 1 )
  
  xx = res[c(jj), ]
  xx = xx[order(-xx$log2FC), ]
  
  xx = xx[which(xx$min < 4.5), ]
  
  source('Functions_atac.R')
  
  ##########################################
  # heatmap for all regeneration dynamic peaks
  ##########################################
  keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), ]
  keep = as.matrix(keep)
  
  conds = c("Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", 'BL_UA_13days_distal',
            "Embryo_Stage40", "Embryo_Stage44_proximal", 'Embryo_Stage44_distal')
  
  sample.sels = c(); cc = c()
  for(n in 1:length(conds)) {
    kk = which(design$conds == conds[n])
    if(length(unique(design$batch[kk])) > 1) kk = kk[which(design$batch[kk] == '2021')]
    sample.sels = c(sample.sels, kk)
    cc = c(cc, as.character(design$conds[kk]))
  }
  kk = sample.sels
  
  keep = keep[, kk]
  df <- data.frame(cc)
  rownames(df) = colnames(keep)
  colnames(df) = 'regeneration time'
  
  ii.gaps = seq(2, nrow(df), by = 2)
  
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps,
           filename = paste0(resDir, '/heatmap_regenerationPeaks_M0prob0.05_log2FC.1.pdf'), 
           width = 10, height = 14)
  
  if(saveTable){
    write.csv(data.frame(xx, keep, stringsAsFactors = FALSE), 
              file = paste0(resDir, '/regeneration_peaks_all.csv'), 
              quote = FALSE, row.names = TRUE)
    
  }
  
  ##########################################
  # highligh potential regeneration peaks, not found in mUA and embryo stages only in regeneration process 
  ##########################################
  jj1 = grep('Embryo_Stage40', colnames(keep))
  jj2 = grep('Embryo_Stage44', colnames(keep))
  jj3 = grep('Mature_UA', colnames(keep))
  mean1 = apply(keep[, jj1], 1, mean)
  mean2 = apply(keep[, jj2], 1, mean)
  mean3 = apply(keep[, jj3], 1, mean)
  
  kk = which(mean1< 4.2 & mean2 < 4.2 & mean3 < 4.2)
  
  pheatmap(keep[kk, ], cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps,
           filename = paste0(resDir, '/heatmap_regenerationPeaks_M0prob0.05_log2FC.1_regeneartion.specific.pdf'), 
           width = 8, height = 10)
  
  yy = keep[kk,]
  
  
  if(saveTable){
    write.csv(data.frame(xx[kk, ], keep[kk, ], stringsAsFactors = FALSE), 
              file = paste0(resDir, '/regeneration_peaks_regeneration.specific.csv'), 
              quote = FALSE, row.names = TRUE)
    
  }
  
  ##########################################
  # first motif activity analysis for temporally dynamic peaks 
  ##########################################
  source('MARA_functions.R')
  xx = run.MARA.atac.temporal(keep, cc)
  
}
s