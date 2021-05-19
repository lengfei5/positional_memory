##########################################################################
##########################################################################
# Project: positional memory
# Script purpose: 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Feb 22 11:02:08 2021
##########################################################################
##########################################################################
rm(list=ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

source('functions_chipSeq.R')
source('Functions.R')

version.analysis = 'atac_rna_seq_integration_analysis_20210323'
#peakDir = "Peaks/macs2_broad"

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require("pheatmap")

##########################################
# some annotations for all analysis 
##########################################
Import.HoxCluster.annotation = TRUE
promoters = readRDS(file = paste0(RdataDir, 
                                  '/axolotl_promoters_Granges_AmexT_v47_putative.full.length_N.terminal_upstream2kb.downstream2kb.rds'))
#promoters =  reduce(promoters, drop.empty.ranges = TRUE, ignore.strand = TRUE)
#promoters.dev = readRDS(file= paste0(RdataDir, 
#        '/axolotl_promoters_Granges_AmexT_v47_putative.full.length_N.terminal_upstream2kb.downstream2kb_rmHouseKeepingGenes.rds'))

if(Import.HoxCluster.annotation){
  HoxA = data.frame(chr = 'chr2p', start = 873085043, end = 884416919, strand = '*', stringsAsFactors = FALSE)
  HoxA = makeGRangesFromDataFrame(HoxA, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  HoxD1 = data.frame(chr = 'chr9q', start = 416423355, end = 427456848, strand = '*', stringsAsFactors = FALSE)
  HoxD1 = makeGRangesFromDataFrame(HoxD1, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  HoxD2 = data.frame(chr = 'chr9q', start = 426557960, end = 435711507, strand = '*', stringsAsFactors = FALSE)
  HoxD2 = makeGRangesFromDataFrame(HoxD2, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  Hoxs = c(HoxA, HoxD1, HoxD2)
}

########################################################
########################################################
# Section I : quatitative analysis of peaks: bineary 0 (no peak) and 1 (peak)  
# compare the peak numbers between development, mature and regeneration 
# 
########################################################
########################################################
Compare.peaks.nb.overlapping = FALSE
if(Compare.peaks.nb.overlapping){
  source('functions_chipSeq.R')
  peakDir = '../Data/R10723_atac/calledPeaks_pval.0.001/'
  
  peak.files = list.files(path = paste0(peakDir, 'macs2'), 
                          pattern = '*macs2_peaks.xls', full.names = TRUE)
  # filter low-quality regeneration samples
  peak.files = peak.files[grep('90392|90393|108072|108073|108070|108071', peak.files, invert = TRUE)]
  
  
  ##########################################
  # three big groups: dev, mature, regeneration
  ##########################################
  compare.dev.mature.reg = FALSE
  if(compare.dev.mature.reg){
    peaks.dev = peak.files[grep('_Stage', peak.files)]
    peaks.mature = peak.files[grep('Mature', peak.files)]
    peaks.reg = peak.files[grep('BL_UA', peak.files)]
    #bam.list = list.files(path = '../Data/R10723_atac/alignments/BAMs_uniq_rmdup', pattern = '*.bam$', full.names = TRUE)
    
    peaks.dev = merge.peaks.macs2(peaks.dev, pcutoff = 6)
    peaks.mature = merge.peaks.macs2(peaks.mature, pcutoff = 6)
    peaks.reg = merge.peaks.macs2(peaks.reg, pcutoff = 6)
    
    
    p1 = peaks.dev[overlapsAny(peaks.dev, promoter, ignore.strand = TRUE)]
    p2 = peaks.mature[overlapsAny(peaks.mature, promoter, ignore.strand = TRUE)]
    p3 = peaks.reg[overlapsAny(peaks.reg, promoter, ignore.strand = TRUE)]
    
    Save.peak.coordinates = FALSE
    if(Save.peak.coordinates){
      saveRDS(p1, file = paste0(RdataDir, '/ATACseq_embryo_peaks_within.axolotl.promoters_Granges_Amex.rds'))
      saveRDS(p2, file = paste0(RdataDir, '/ATACseq_mature_peaks_within.axolotl.promoters_Granges_Amex.rds'))
      
      xx = as.data.frame(p1)
      xx$peak.name = paste0(xx$seqnames, ":", xx$start, "_", xx$end)
      xx = xx[, c(1:3, 6, 5)]
      write.table(xx, file = paste0(peakDir, 'merge_peak_Embryo_pval.6_promoters.1000bpUp.200bpDown.bed'), 
                  sep = '\t', row.names = FALSE, 
                  col.names = FALSE, quote = FALSE)
    }
    
    # all peaks overlap between Dev and Mature 
    pp = list(peaks.dev = peaks.dev, peaks.mature= peaks.mature)
    ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='peak number and overlapping Dev vs Mature')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    # promoter peaks overlap between Dev  
    pp = list(peaks.dev = p1, peaks.mature=p2)
    ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='promoter peak overlapping Dev vs Mature')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    
    # all peaks overlap between Dev, Mature and Regeneration 
    pp = list(peaks.dev = peaks.dev, peaks.mature= peaks.mature, peaks.reg = peaks.reg)
    ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='peak number and overlapping Dev, Mature and Reg')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    # promoter peaks overlap between Dev  
    pp = list(peaks.dev = p1, peaks.mature=p2, peaks.reg = p3)
    ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='promoter peak overlapping Dev vs Mature')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
  }
  
  ##########################################
  # compare mature UA, UA.day5, UA.day9   
  ##########################################
  compare.mature.UAregeneration = FALSE
  if(compare.mature.UAregeneration){
    p1 = merge.peaks.macs2(peak.files[grep('Mature_UA', peak.files)], pcutoff = 6)
    p2 = merge.peaks.macs2(peak.files[grep('BL_UA_5days', peak.files)], pcutoff = 6)
    p3 = merge.peaks.macs2(peak.files[grep('BL_UA_9days', peak.files)], pcutoff = 6)
    
    pp = list(Mature.UA = p1, BL.UA.D5= p2, BL.UA.D9 = p3)
    ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='peak number and overlapping Dev, Mature and Reg')
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    ##########################################
    # compare stage 40, mature.UA, UA5, UA9, UA13P, UA13D (only samples from 2021) 
    ##########################################
    p0 = merge.peaks.macs2(peak.files[grep('Stage40_13', peak.files)], pcutoff = 6)
    p1 = merge.peaks.macs2(peak.files[grep('Mature_UA_13', peak.files)], pcutoff = 6)
    p10 = merge.peaks.macs2(peak.files[grep('Mature_UA_74938|Mature_UA_102655', peak.files)], pcutoff = 6)
    p11 = merge.peaks.macs2(peak.files[grep('Mature_LA', peak.files)], pcutoff = 6)
    p12 = merge.peaks.macs2(peak.files[grep('Mature_Hand', peak.files)], pcutoff = 6)
    
    p2 = merge.peaks.macs2(peak.files[grep('BL_UA_5days_13', peak.files)], pcutoff = 6)
    p3 = merge.peaks.macs2(peak.files[grep('BL_UA_9days_13', peak.files)], pcutoff = 6)
    p4 = merge.peaks.macs2(peak.files[grep('BL_UA_13days_proximal_13', peak.files)], pcutoff = 6)
    p5 = merge.peaks.macs2(peak.files[grep('BL_UA_13days_distal_13', peak.files)], pcutoff = 6)
    
    peaks.ref = list(stage40 = p0, mUA = p1, mUL = p11, mHand = p12,  UA5= p2, UA9 = p3, UA13P = p4, UA13D = p5)
    
    Save.peak.coordinates.HoxClusters = FALSE
    if(Save.peak.coordinates.HoxClusters){
      load(file = paste0(RdataDir, '/samplesDesign_readCounts.withinPeaks.pval6.Rdata'))
      
      pp = data.frame(t(sapply(counts$gene, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
      pp$strand = '*'
      pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                    start.field="X2", end.field="X3", strand.field="strand")
      
      Hoxs =GenomicRanges::union(GenomicRanges::union(HoxA, HoxD1), HoxD2)
      
      pxx = pp[overlapsAny(pp, Hoxs)]
      
      require(ChIPpeakAnno)
      require(ChIPseeker)
      
      peakAnnots = as.data.frame(annotatePeak(pxx, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript'))
      
      xx = matrix(0, ncol = length(peaks.ref), nrow = nrow(peakAnnots))
      colnames(xx) = names(peaks.ref)
      for(n in 1:ncol(xx)){
        xx[overlapsAny(pxx, peaks.ref[[n]]), n] = 1
      }
      
      xx = data.frame(peakAnnots, xx, stringsAsFactors = FALSE)
      
      cne.hoxa = data.frame(read.delim(file = '../AkaneToJingkuiShareFiles/HoxA_Candidates_axCNEs/Hoxa_neiboughr_CNEs.bed',
                                       header = FALSE))
      colnames(cne.hoxa) = c('chr', 'start', 'end', 'name')
      cne.hoxa$name = paste0(cne.hoxa$name, '.HoxA')
      
      cne.hoxd = data.frame(read.delim(file = '../AkaneToJingkuiShareFiles/HoxD_Candidates_axCNEs/Axolotl_HoxD_CNEs_NEW/combined.bed',
                                       header = FALSE))
      colnames(cne.hoxd) = c('chr', 'start', 'end', 'name')
      cne.hoxd$name = paste0(cne.hoxd$name, '.HoxD')
      
      cnes = data.frame(rbind(cne.hoxa, cne.hoxd), stringsAsFactors = FALSE)
      cnes$strand = '*'
      
      
      jj = which(cnes$end <= cnes$start)
      test = cnes$start[jj]
      cnes$start[jj] = cnes$end[jj]
      cnes$end[jj] = test
      cnes$coordinates = paste0(cnes$chr, ':', cnes$start, '-', cnes$end)
      cnes = cnes[match(unique(cnes$coordinates), cnes$coordinates), c(1:5)]
      
      cnes = makeGRangesFromDataFrame(cnes, seqnames.field=c("chr"),
                                      start.field="start", end.field="end", strand.field = 'strand',  keep.extra.columns = TRUE)
      
      kk = GenomicRanges::findOverlaps(pxx, cnes, ignore.strand = TRUE, type = 'any', select = 'all')
      kk = data.frame(kk)
      
      cnes = data.frame(cnes)
      
      xx$cnes = NA
      
      for(n in 1:nrow(kk))
      {
        # n = 1
        if(is.na(xx$cnes[kk[n, 1]])) {
          xx$cnes[kk[n, 1]] = cnes$name[kk[n, 2]]      
        }else{
          xx$cnes[kk[n, 1]] = paste0(xx$cnes[kk[n, 1]], ';', cnes$name[kk[n, 2]])      
        }
        
      }
      
      write.table(xx, file = paste0(resDir, '/HoxCluster.peaks_annotation_CNEs.txt'), 
                  sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
      
      write.table(xx[, c(1, 2, 3, 12, 5)], file = paste0(resDir, '/HoxCluster.peaks_annotation.bed'), sep = '\t', 
                  col.names = FALSE, rown)   
    }
    
    library(Vennerable)
    # did not work; there are too many groups
    ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='all peaks', fill=c(1,2,3,4,5))
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    # annotate peaks and select intron and intergeic regions
    require(ChIPpeakAnno); require(ChIPseeker)
    xx = sapply(pp, annotatePeak, TxDb=amex, tssRegion = c(-2000, 2000))
    
    pxx = pp
    for(n in 1:length(xx)){
      # n = 1
      p0 = pp[[n]]
      x0 = xx[[n]]
      x0 = as.data.frame(xx[[n]])
      
      ii = which(overlapsAny(p0, HoxD2, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', x0$annotation))
      
      pxx[[n]] = p0[ii]
    }
    
    library(Vennerable)
    # did not work; there are too many groups
    ol.peaks <- makeVennDiagram(pxx, NameOfPeaks=names(pxx), connectedPeaks="keepAll", main='peak overlapping of HoxD late TAD ',
                                fill=c(1,2,3,4,5))
    
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
  }
  
  ##########################################
  # genomic features distribution of peaks in all conditiosn 
  ##########################################
  Characterize.genomic.feature.distrubtion.all.samples = FALSE
  if(Characterize.genomic.feature.distrubtion.all.samples){
    p0 = merge.peaks.macs2(peak.files[grep('Stage40_13', peak.files)], pcutoff = 6)
    p00 = merge.peaks.macs2(peak.files[grep('Stage40_933', peak.files)], pcutoff = 6)
    p01 = merge.peaks.macs2(peak.files[grep('Stage44_proximal_933', peak.files)], pcutoff = 6)
    p02 = merge.peaks.macs2(peak.files[grep('Stage44_distal_933', peak.files)], pcutoff = 6)
    
    p1 = merge.peaks.macs2(peak.files[grep('Mature_UA_13', peak.files)], pcutoff = 6)
    p10 = merge.peaks.macs2(peak.files[grep('Mature_UA_74938|Mature_UA_102655', peak.files)], pcutoff = 6)
    p11 = merge.peaks.macs2(peak.files[grep('Mature_LA', peak.files)], pcutoff = 6)
    p12 = merge.peaks.macs2(peak.files[grep('Mature_Hand', peak.files)], pcutoff = 6)
    
    p2 = merge.peaks.macs2(peak.files[grep('BL_UA_5days_13', peak.files)], pcutoff = 6)
    p20 = merge.peaks.macs2(peak.files[grep('BL_UA_5days_89', peak.files)], pcutoff = 6)
    
    p3 = merge.peaks.macs2(peak.files[grep('BL_UA_9days_13', peak.files)], pcutoff = 6)
    p4 = merge.peaks.macs2(peak.files[grep('BL_UA_13days_proximal_13', peak.files)], pcutoff = 6)
    p5 = merge.peaks.macs2(peak.files[grep('BL_UA_13days_distal_13', peak.files)], pcutoff = 6)
    
    pp = list(stage40 = p0, stage40.old = p00,  stage44P = p01, stage44D = p02,
              Mature.UA = p1, Mature.UA.old = p10, Mature.LA = p11, Mature.Hand = p12,  
              UA5= p2, UA5.old = p20,  
              UA9 = p3,  UA13P = p4, UA13D = p5
    )
    
    xx = sapply(pp, annotatePeak, TxDb=amex, tssRegion = c(-2000, 2000))
    #peakAnnots = annotatePeak(p0, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    peakAnnots = xx
    
    
    #xx = data.frame(peakAnnots)
    pdf(paste0(resDir, "/Peaks_genomicFeatures_distance_TSS.pdf"), width = 10, height = 6)
    #par(mfrow=c(1,2))
    print(plotAnnoBar(peakAnnots))
    print(plotDistToTSS(peakAnnots))
    dev.off()
    
    
  }
  
  ##########################################
  # qualitatively characterize position-restricted peaks in mature samples 
  ##########################################
  Characterize.position.restricted.peaks = FALSE
  if(Characterize.position.restricted.peaks){
    source('functions_chipSeq.R')
    
    p1 = merge.peaks.macs2(peak.files[grep('Mature_UA_13', peak.files)], pcutoff = 3, select.overlappingPeaks = TRUE)
    #xx = merge.peaks.macs2(peak.files[grep('Mature_UA_13', peak.files)], pcutoff = 6, select.overlappingPeaks = FALSE)
    p10 = merge.peaks.macs2(peak.files[grep('Mature_UA_74938|Mature_UA_102655', peak.files)], pcutoff = 3, 
                            select.overlappingPeaks = TRUE)
    p11 = merge.peaks.macs2(peak.files[grep('Mature_LA', peak.files)], pcutoff = 3, select.overlappingPeaks = TRUE)
    p12 = merge.peaks.macs2(peak.files[grep('Mature_Hand', peak.files)], pcutoff = 3, select.overlappingPeaks = TRUE)
    
    pp = list(mUA = p1, mUA.old = p10, mLA = p11, mHand = p12)
    
    library(Vennerable)
    ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='peak overlapping in Mature',
                                fill=c(1:length(pp)))
    pxx = pp[c(1:2)]
    makeVennDiagram(pxx, NameOfPeaks=names(pxx), connectedPeaks="keepAll", main='peak overlapping in Mature',
                    fill=c(1:length(pxx)))
    
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    
    # pool peaks in all mature samples and get a matrix indicating where they were detected
    pxx = union(union(union(p1, p10,  ignore.strand=TRUE), p11,  ignore.strand=TRUE), p12, ignore.strand=TRUE)
    pxx = reduce(pxx)
    
    peaks.ref = pp
    
    xx = matrix(0, ncol = length(peaks.ref), nrow = length(pxx))
    colnames(xx) = names(peaks.ref)
    rownames(xx) = paste0(seqnames(pxx), ':', start(pxx), '-', end(pxx))
    
    for(n in 1:ncol(xx)){
      xx[overlapsAny(pxx, peaks.ref[[n]]), n] = 1
    }
    ss = apply(xx, 1, sum)
    
    xx = xx[which(ss < 4), ]
    xx = xx[which(xx[,1] == xx[, 2]), ]
    
    ss = apply(xx, 1, sum)
    xx = as.matrix(xx)
    
    pheatmap(xx[c(1:10000), ], cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = TRUE,
             cluster_cols=FALSE)
     
  }
  
  xx = read.delim(file = 
                    '../results/R10723_Rxxxx_atacseq_bowtie2.newParam_mtDNA_picardrmdup_20210208/HoxCluster.peaks_annotation_CNEs.txt',
                  sep = '\t', header = TRUE)
  
  #xx = xx[which(xx$seqnames == 'chr2p'), ]
  #ss = apply(as.matrix(xx[, c(16, 19:22)]), 1, sum)
  #jj = which(xx$UA13D == 1 & ss == 1)
  #xx = xx[jj, ]
  
  #write.table(xx, file = paste0(resDir, '/peakList_coordinates_BL.UA.D13.Distal.txt'), sep = '\t', row.names = FALSE,
  #            col.names = FALSE)
   
  names = paste0(xx$seqnames, ':', xx$start, '-', xx$end)
  
  source('Functions.R')
  
  pdfname = paste0(resDir, "/peak_profiles.pdf")
  pdf(pdfname, width = 10, height = 6)
  par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(6,3,2,0.2), tcl = -0.3)
  
  plot.peak.profiles(peak.name = names)
  
  dev.off()
  
}

########################################################
########################################################
# Section II: quantitative analysis (normalization, clustering, heatmpap) 
# Embryo 40, Embryo 44 (proximal, distal)
# Mature arms: UA, LA and Hand together   
########################################################
########################################################
load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.withinPeaks.pval6.Rdata'))

##########################################
# Quick check the sample infos
##########################################
Manual.change.sample.infos = FALSE
if(Manual.change.sample.infos){
  colnames(counts) = c('gene', paste0(design$fileName, '_', design$SampleID))
  design$total = unlist(design$total)
  design$unique.rmdup = unlist(design$unique.rmdup)
  design$pct.usable = design$unique.rmdup/design$total
  
  design$batch = '2020'
  design$batch[grep('136|137', design$SampleID)] = '2021'
  
  rownames(counts) = counts$gene
  counts = as.matrix(counts[, -1])
  
  save(counts, design, file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.withinPeaks.pval6.Rdata'))
}

hist(design$unique.rmdup/design$total, breaks = 16, main = 'pct of usable reads/total ')


##########################################
# start the DESeq2 and quantile normalization and other normalization will be tested
##########################################
sels = grep('Mature|Embryo|BL_UA', design$conds)

dds <- DESeqDataSetFromMatrix(as.matrix(counts[, sels]), DataFrame(design[sels, ]), design = ~ conds)

#rm(counts)

# check the peak length
peakNames = rownames(dds)
pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))

pp$strand = '*'
pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")
ll = width(pp)


select.peaks.with.readThreshold = TRUE
select.background.for.peaks = TRUE

if(select.peaks.with.readThreshold){
  #ss = rowMax(counts(dds)[, grep('Embryo_', dds$conds)])
  ss = rowMax(counts(dds))/ll*500
  
  hist(log10(ss), breaks = 200, main = 'log2(sum of read within peaks) ')
  cutoff.peak = 40
  cutoff.bg = 10
  cat(length(which(ss >= cutoff.peak)), 'peaks selected with minimum read of the highest peak -- ', cutoff.peak,  '\n')
  cat(length(which(ss < cutoff.bg)), 'peaks selected with minimum read of the highest peak -- ', cutoff.bg,  '\n')
  abline(v= log10(cutoff.peak), col = 'red', lwd = 2.0)
  abline(v= log10(cutoff.bg), col = 'blue', lwd = 2.0)
  
  if(select.background.for.peaks){
    ii = which(ss >= cutoff.peak)
    ii.bg = which(ss < cutoff.bg)
    ii.bg = sample(ii.bg, size = 1000, replace = FALSE)
    rownames(dds)[ii.bg] = paste0('bg_', rownames(dds)[ii.bg])
    dds = dds[c(ii, ii.bg), ]
    ll.sels = ll[c(ii, ii.bg)]
    
  }else{
    dds <- dds[ss >= cutoff.peak, ]
    ll.sels = ll[ss >= cutoff.peak]
  }
  
}

dds <- estimateSizeFactors(dds)
plot(sizeFactors(dds), colSums(counts(dds))/median(colSums(counts(dds))), log = 'xy')

plot(sizeFactors(dds), design$unique.rmdup, log = 'xy')
text(sizeFactors(dds), design$unique.rmdup, labels = design$samples, cex = 0.7)

save.scalingFactors.for.deeptools = FALSE
if(save.scalingFactors.for.deeptools){
  xx = data.frame(sampleID = design$SampleID,  
                  scalingFactor = design$unique.rmdup/(sizeFactors(dds)*median(design$unique.rmdup)),
                  stringsAsFactors = FALSE)
  
  write.table(xx, file = paste0(resDir, '/DESeq2_scalingFactor_forDeeptools.txt'), sep = '\t',
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

##########################################
# QCs again for embryo and mature samples
# Control the peak quality by checking the house-keeping genes and other-tissue specific genes, 
# GO-enrichment analysis
##########################################
QC.PLOT = FALSE
if(QC.PLOT){
  pdfname = paste0(resDir, "/atacseq_Embryo_Mature_QCs_pval6.pdf")
  pdf(pdfname, width = 12, height = 10)
  
  Check.RNAseq.Quality(read.count=counts(dds)[c(1:12000),], design.matrix = data.frame(design$SampleID, 
                                                                                       design$conds, design$batch))
  dev.off()
  
}

Check.promoter.peak.enrichment = FALSE
if(Check.promoter.peak.enrichment){
  source('Functions.R')
  DoubleCheck.promoter.peaks.enrichment(fpm)

}

##########################################
# test normalization and batch correction of ATAC-seq data
# TMM and combat were selected for normalization and batch correction
##########################################
Test.atac.normalization.batch.correction = FALSE
if(Test.atac.normalization.batch.correction){
  source('Functions.R')
  
  #norms = normalize.batch.correct(dds, design, norm.batch.method ='TMM.combat')
  library(edgeR)
  
  d <- DGEList(counts=counts(dds), group=design$conds)
  tmm <- calcNormFactors(d, method='TMM')
  fpm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  rm(tmm)
  rm(d)
  #fpm = log2(tmm + 1)
  
  require("sva")
  bc = as.factor(design$batch)
  mod = model.matrix(~ as.factor(conds), data = design)
  fpm.bc = ComBat(dat=fpm, batch=bc, mod=mod, par.prior=TRUE, ref.batch = '2021')    
  fpm = fpm.bc
  
  make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'all')
  make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'Dev.Mature')
  
  rm(fpm.bc)
  saveRDS(fpm, file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
  
  # normalize the peak width to have fpkm
  FPKM.normalization = FALSE
  if(FPKM.normalization){
    library(preprocessCore)
    
    peakNames = rownames(fpm)
    peakNames = gsub('bg_', '', peakNames)
    peakNames = gsub('_', '-', peakNames)
    
    pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    ll = width(pp)
    
    fpkm = fpm.bc
    for(n in 1:ncol(fpkm))
    {
      fpkm[,n] = fpkm[,n]/ll*10^3
    }
    
    fpkm.qn = normalize.quantiles(fpkm)
    colnames(fpkm.qn) = colnames(fpkm)
    rownames(fpkm.qn) = rownames(fpkm)
    fpkm = fpkm.qn
    #rm(fpm.qn)
    make.pca.plots(fpkm.qn, ntop = 5000, conds.plot = 'all')
    
    rm(fpkm.qn)
    
    save(fpm, fpkm, file = paste0(RdataDir, '/fpm_TMM_combat_fpkm_quantileNorm.Rdata'))
  }
 
  
}

#dds <- estimateDispersions(dds)
#plotDispEsts(dds, ymin = 10^-4)
########################################################
########################################################
# Section : grouping atac-seq peak profiles
# 1) first identify static peaks (probably most of them) with model selection 
# (linear regression with only intercept and nonlinear fitting with spline or GAM)
# 2) grouping dynamic peaks using DP-GP
# 
########################################################
########################################################
Grouping.atac.peaks = FALSE
if(Grouping.atac.peaks){
  
  load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.withinPeaks.pval6.Rdata'))
  fpm = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
  
  # prepare the background distribution
  jj = grep('bg_', rownames(fpm), invert = TRUE)
  fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
  fpm = fpm[jj, ]
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  creat.Granges.and.peakAnnotation = TRUE
  # make Granges and annotate peaks
  if(creat.Granges.and.peakAnnotation){
    pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    require(ChIPpeakAnno)
    require(ChIPseeker)
    
    amex = GenomicFeatures::makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))
    #amex = makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))
    #amex = readRDS(file = paste0(annotDir, 'TxDb_ax6_UCSC_2021_01_26_genes.putative.full.length.rds')) # reimport object does not work
    pp.annots = as.data.frame(annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript'))
    
  }
  
  ##########################################
  # all-peaks test M0 (static peaks), M1 (dynamic peaks above background), M2 (dyanmic peaks with some condtions below background)
  ##########################################
  make.test.All.Peaks = FALSE
  if(make.test.All.Peaks){
    conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", "Embryo_Stage44_distal",
              "Mature_UA",  "Mature_LA", "Mature_Hand", 
              "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", "BL_UA_13days_distal")
    
    # examples to test
    test.examples = c('HAND2', 'FGF8', 'KLF4', 'Gli3', 'Grem1')
    #test.examples = c('Hoxa13')
    ii.test = which(overlapsAny(pp, promoters[which(!is.na(match(promoters$geneSymbol, test.examples)))]))
    #ii.Hox = which(overlapsAny(pp, Hoxs))
    #ii.test = unique(c(ii.test, ii.Hox))
    
    sample.sels = c()
    cc = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n] & design$SampleID != '136159')
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    library(tictoc)
    ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
    
    source('Functions.R')
    tic() 
    res = t(apply(fpm[ii.test, sample.sels], 1, all.peaks.test, c = cc))
    res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
    toc()
    
    saveRDS(res, file = paste0(RdataDir, '/res_allPeaks_test.rds'))
    
    jj = which(res$prob.M0 < 0.5 & res$log2FC >1)
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC.mature), ]
    
    length(which(xx$min.mature <3.0))
    length(which(xx$min.mature <2.5))
    length(which(xx$min.mature <2.))
    
    
    bgs = data.frame(bg = as.numeric(fpm.bg))
    ggplot(bgs, aes(x=bg)) + geom_histogram(binwidth = 0.25, color="darkblue", fill="lightblue")
    
    #xx = xx[which(xx$min.m < 3 & xx$max >3), ]
    
    #write.table(xx, file = paste0(resDir, '/HoxClusters_spatialDynamic_peaks.txt'), 
    #            col.names = TRUE, row.names = TRUE, sep = '\t', quote = FALSE)
    
    mins = res$min
    
    nb.openloci = c()
    thresholds = seq(2, 4, by = 0.1)
    for(cutoff in thresholds)
    {
      nb.openloci = c(nb.openloci, length(which(mins>cutoff)))
    }
    
    opens = data.frame(thresholds = thresholds, nb.open = nb.openloci/nrow(fpm))
    
    ggplot(opens, aes(x = thresholds, y = nb.open)) +
      geom_point(color = 'blue') + geom_line() +
      geom_vline(xintercept = 3.0, color = 'red') +
      ggtitle('% of peaks above the bg thresholds')
    
    
    xx = res
    xx$annotation[grep('Promoter', xx$annotation)] = 'Promoter'
    xx$annotation[grep('Intron', xx$annotation)] = 'Intron'
    xx$annotation[grep('Exon', xx$annotation)] = 'Exon'
    xx$annotation[grep('Downstream', xx$annotation)] = 'Downstream'
    
    aa = c(table(xx$annotation), table(xx$annotation[which(xx$min>3)]))
    aa = data.frame(names(aa), aa, rep(c('Peak.all', 'Peak.above.bg'), each = 7))
    colnames(aa) = c('peakAnnot', 'nb.peaks', 'groups')
    #aa$peakAnnot = as.factor(as.character(aa$peakAnnot))
    #aa = aa[, -3]
    
    ggplot(data=aa, aes(x=peakAnnot, y=nb.peaks, fill=groups)) +
      geom_bar(stat="identity", color="black", position=position_dodge()) +
      theme_minimal()  + scale_fill_manual(values=c('#999999','#E69F00'))
    
  }
  
  ##########################################
  # temporal-peaks test
  ##########################################
  grouping.temporal.peaks = FALSE
  if(grouping.temporal.peaks){
    conds = c("Embryo_Stage40", "Embryo_Stage44_proximal",
              "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal")
    
    # examples to test
    test.examples = c('HAND2', 'FGF8', 'KLF4', 'Gli3', 'Grem1')
    #test.examples = c('Hoxa13')
    ii.test = which(overlapsAny(pp, promoters[which(!is.na(match(promoters$geneSymbol, test.examples)))]))
    #ii.Hox = which(overlapsAny(pp, Hoxs))
    #ii.test = unique(c(ii.test, ii.Hox))
    
    sample.sels = c()
    cc = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n] & design$SampleID != '136159')
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    library(tictoc)
    ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
    
    source('Functions.R')
    tic() 
    res = t(apply(fpm[ii.test, sample.sels], 1, temporal.peaks.test, c = cc))
    res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
    toc()
    
    saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test.rds'))
    
    res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test.rds'))
    # select the temporal dynamic peaks
    length(which(res$prob.M0<0.05))
    length(which(res$prob.M0<0.05 & res$log2FC > 1))
    jj = which(res$prob.M0 < 0.05 & res$log2FC >1 )
    
    jj = which(res$prob.M0 < 0.01 & res$log2FC >2 )
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC), ]
    
    source('Functions.R')
    keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
    keep = as.matrix(keep)
    
    df <- data.frame(cc)
    rownames(df) = colnames(keep)
    
    ii.gaps = c(6, 9)
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    
    ##########################################
    # first motif activity analysis for temporally dynamic peaks 
    ##########################################
    source('Functions.R')
    xx = run.MARA.atac.temporal(keep, cc)
    
      
  }
  
  ##########################################
  # position-dependent test
  ##########################################
  grouping.position.dependent.peaks = FALSE
  if(grouping.position.dependent.peaks){
    # examples to test 
    test.examples = c('HAND2', 'FGF8', 'KLF4', 'Gli3', 'Grem1')
    #test.examples = c('Hoxa13')
    ii.test = which(overlapsAny(pp, promoters[which(!is.na(match(promoters$geneSymbol, test.examples)))]))
    #ii.Hox = which(overlapsAny(pp, Hoxs))
    #ii.test = unique(c(ii.test, ii.Hox))
    
    #conds = as.character(unique(design$conds))
    conds = c("Mature_UA", "Mature_LA", "Mature_Hand", 
              "Embryo_Stage44_proximal", "Embryo_Stage44_distal",
              "BL_UA_13days_proximal", "BL_UA_13days_distal")
    sample.sels = c()
    cc = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n] & design$SampleID != '136159')
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    
    library(tictoc)
    ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
    source('Functions.R')
    tic() 
    res = t(apply(fpm[ii.test, sample.sels], 1, spatial.peaks.test, c = cc))
    res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
    toc()
    
    saveRDS(res, file = paste0(RdataDir, '/res_position_dependant_test.rds'))
    
    res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test.rds'))
    
    # select the spatially dynamic peaks
    jj = which(res$prob.M0.mature < 0.3 & res$log2FC.mature >1)
    #jj1 = which(res$pval.embryo < 0.01 & res$log2FC.embryo >1)
    #jj2 = which(res$pval.BL < 0.01 & res$log2FC.BL >1)
    
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC.mature), ]
    
    length(which(xx$min.mature <3.0))
    length(which(xx$min.mature <2.5))
    length(which(xx$min.mature <2.))
    
    library(ggplot2)
    bgs = data.frame(bg = as.numeric(fpm.bg))
    ggplot(bgs, aes(x=bg)) + geom_histogram(binwidth = 0.25, color="darkblue", fill="lightblue")
    
    #xx = xx[which(xx$min.m < 3 & xx$max >3), ]
    
    #write.table(xx, file = paste0(resDir, '/HoxClusters_spatialDynamic_peaks.txt'), 
    #            col.names = TRUE, row.names = TRUE, sep = '\t', quote = FALSE)
    
    mins = apply(cbind(res$max.mature, res$min.BL, res$min.embryo), 1, min)
    
    nb.openloci = c()
    thresholds = seq(2, 4, by = 0.1)
    for(cutoff in thresholds)
    {
      nb.openloci = c(nb.openloci, length(which(mins>cutoff)))
    }
    
    opens = data.frame(thresholds = thresholds, nb.open = nb.openloci/nrow(fpm))
    
    ggplot(opens, aes(x = thresholds, y = nb.open)) +
      geom_point(color = 'blue') + geom_line() +
      geom_vline(xintercept = 3.0, color = 'red') +
      ggtitle('% of all peaks above the bg thresholds')
    
    
    # visualize the position-dependent peaks
    source('Functions.R')
    
    keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
    keep = as.matrix(keep)
    
    # for(c in c('Mature', 'Embryo', 'BL_UA'))
    # {
    #   jj = grep(c, cc)
    #   keep[,jj] = t(apply(keep[, jj], 1, scale, scale = FALSE))
    # }
    
    df <- data.frame(cc)
    rownames(df) = colnames(keep)
    
    ii.gaps = c(7, 11)
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    
    pdfname = paste0(resDir, "/peak_profiles_test.static.peaks.pdf")
    pdf(pdfname, width = 10, height = 6)
    par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(6,3,2,0.2), tcl = -0.3)
    
    mains = signif(res, d = 2)
    plot.peak.profiles(peak.name = names, fpm = fpm, mains = mains)
    
    dev.off()
    
  }
  
}


########################################################
########################################################
# Section : first analysis for subgrouped peaks, promoter, enhancers
# 
########################################################
########################################################


##########################################
# overview all peaks with heatmap for subgrouped peaks
##########################################
source('Functions.R')

make.heatmap.atacseq()

sel.promoter = which(overlapsAny(pp, promoters, ignore.strand = TRUE) == TRUE)
sel.promoter.dev = which(overlapsAny(pp, promoters.dev, ignore.strand = TRUE) == TRUE)

ii.promoter = grep('Promoter', peakAnnots$annotation)
sel.intron = grep('Intron ', peakAnnots$annotation)
sel.intergen = grep('Intergenic', peakAnnots$annotation)

cat(length(sel.promoter), ' peaks in promoters\n')
cat(length(sel.promoter.dev), ' peaks in promoters.dev \n')

cat(length(sel.intron), ' peaks in introns\n')
cat(length(sel.intergen), ' peaks in intergenic regions\n')

ii.HoxA = which(overlapsAny(pp, HoxA, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))
ii.HoxD1 = which(overlapsAny(pp, HoxD1, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))
ii.HoxD2 = which(overlapsAny(pp, HoxD2, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))

##########################################
#  select conditions of interest
##########################################
#conds.sel = c('Embryo_Stage40_93',  'Embryo_Stage40_13', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',  'Mature_UA_', 'Mature_LA', 
#              'Mature_Hand', 'BL_UA_5days_13', 'BL_UA_5days_89',  'BL_UA_9days_13', 'BL_UA_13days_proximal', 'BL_UA_13days_distal')

conds.sel = c('Embryo_Stage40_93',  'Embryo_Stage40_13', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
              'Mature_UA_13',  'Mature_LA', 'Mature_Hand', 'BL_UA_5days_13', 'BL_UA_5days_89',  'BL_UA_9days_13', 'BL_UA_13days_proximal',
              'BL_UA_13days_distal')
#conds.sel = c('Mature_UA_13', 'Mature_UA_74938|Mature_UA_102655', 'Mature_LA', 'Mature_Hand')

conds.sel = c('Mature_UA_13', 'BL_UA_5days_13', 'BL_UA_9days_13', 'BL_UA_13days_proximal',  'BL_UA_13days_distal')

sample.sel = c()
ii.gaps = c()
for(n in 1:length(conds.sel)) {
  c = conds.sel[n]
  sample.sel = c(sample.sel, grep(c, colnames(fpm)))
  if(n == 1) {
    ii.gaps = c(ii.gaps, length(grep(c, colnames(fpm))))
  }else{
    if(n != length(conds.sel)) ii.gaps = c(ii.gaps, (ii.gaps[n-1] + length(grep(c, colnames(fpm)))))
  }
}

df <- data.frame(colData(dds)[,c("conds", 'batch')])[sample.sel, ]

test.normalization.quantile = FALSE
if(test.normalization.quantile){
  xx = fpm[, sample.sel]
  plot(xx[, c(1, 3)], cex = 0.1); abline(0, 1, lwd = 2.0, col ='red')
  plot(xx[, c(2, 3)], cex = 0.2); abline(0, 1, lwd = 2.0, col ='red')
  
  library(factoextra)
  res.pca <- prcomp(t(xx), scale = TRUE)
  #res.var <- get_pca_var(res.pca)
  
  fviz_pca_ind(res.pca,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  
}

Split.peaks.into.promoters.intron.integeic = FALSE
i(Split.peaks.into.promoters.intron.integeic){
  
  ##########################################
  # promoter peaks
  ##########################################
  keep =fpm[sel.promoter.dev, sample.sel]
  
  filtering.with.signal = TRUE
  if(filtering.with.signal){
    
    ss = apply((keep), 1, max)
    
    max.cutoff = 3
    
    hist(ss, breaks = 100)
    abline(v = max.cutoff, col='red', lwd =2.0)
    cat(length(which(ss>max.cutoff)), 'peak selected\n')
    
    keep = keep[which(ss > max.cutoff), ]
    
  }
  
  keep = as.matrix(keep)
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
  
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
  
  ##########################################
  # promoter introns
  ##########################################
  keep =fpm[sel.intron, sample.sel]
  
  filtering.with.signal = TRUE
  if(filtering.with.signal){
    
    ss = apply((keep), 1, max)
    
    max.cutoff = 3
    hist(ss, breaks = 100)
    abline(v = max.cutoff, col='red', lwd =2.0)
    cat(length(which(ss>max.cutoff)), '\n')
    
    keep = keep[which(ss>max.cutoff), ]
    
  }
  
  #subsample = sample(c(1:nrow(keep)), 15000)
  #keep = keep[subsample, ]
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
  ##########################################
  # enhancer peaks too many and subseting is necessary
  ##########################################
  keep =fpm[sel.intergen, sample.sel]
  
  filtering.with.signal = TRUE
  if(filtering.with.signal){
    ss = apply((keep), 1, max)
    
    max.cutoff = 3.5
    hist(ss, breaks = 100)
    abline(v = max.cutoff, col='red', lwd =2.0)
    cat(length(which(ss>max.cutoff)), '\n')
    
    keep = keep[which(ss>max.cutoff), ]
    
  }
  
  #subsample = sample(c(1:nrow(keep)), 10000)
  #keep = keep[subsample, ]
  
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
  
  ##########################################
  # peaks in introns and intergenic regions in HoxA and HoxD clusters 
  ##########################################
  keep =fpm[ii.HoxA, sample.sel]
  pheatmap(as.matrix(log2(keep + 1)), cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
  keep =fpm[ii.HoxD1, sample.sel]
  pheatmap(as.matrix(log2(keep + 1)), cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
  
}

##########################################
# position-related mature peaks 
##########################################
Test.position.related.peaks.in.mature = FALSE
if(Test.position.related.peaks.in.mature){
  keep = fpm[, sample.sel]
  
  filtering.with.signal = TRUE
  if(filtering.with.signal){
    ss = apply((keep[, c(1:4)]), 1, mean)
    vars = apply(keep[, c(1:4)], 1, var)
    
    plot(ss, vars, cex = 0.15)
    #abline(h=c(1,2), col = 'red', lwd =2.0)
    abline(v = 1, col = 'red', lwd =2.0)
    jj = which(vars < 1.)
    
    keep = keep[jj, ]
    
    ss = apply(keep[, c(5:6)], 1, mean)
    vars = apply(keep[, c(5:6)], 1, var)
    plot(ss, vars, cex = 0.15)
    abline(h=1, col = 'red', lwd =2.0)
    
    jj = which(vars < 1)
    
    keep = keep[jj, ]
    
    ss = apply(keep[, c(7:8)], 1, mean)
    vars = apply(keep[, c(7:8)], 1, var)
    plot(ss, vars, cex = 0.15)
    
    jj = which(vars < 1)
    
    keep = keep[jj, ]
    
    ss = apply(keep, 1, mean)
    vars = apply(keep, 1, var)
    plot(ss, vars, cex = 0.15)
    
    jj = which(ss>0 & vars >1)
    
    keep = keep[jj, ]
    
    
    max.cutoff = 3.5
    hist(ss, breaks = 100)
    abline(v = max.cutoff, col='red', lwd =2.0)
    cat(length(which(ss>max.cutoff)), '\n')
    
    
  }
  
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
  
}


