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
source('Functions_atac.R')

version.analysis = 'atac_rna_chipseq_analysis_20211007'
#peakDir = "Peaks/macs2_broad"
saveTable = TRUE

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')

require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require("pheatmap")

##########################################
# some annotations for all analysis 
##########################################
Import.HoxCluster.annotation = TRUE

if(Import.HoxCluster.annotation){
  HoxA = data.frame(chr = 'chr2p', start = 873085043, end = 884416919, strand = '*', stringsAsFactors = FALSE)
  HoxA = makeGRangesFromDataFrame(HoxA, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  HoxD1 = data.frame(chr = 'chr9q', start = 416423355, end = 427456848, strand = '*', stringsAsFactors = FALSE)
  HoxD1 = makeGRangesFromDataFrame(HoxD1, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  HoxD2 = data.frame(chr = 'chr9q', start = 426557960, end = 435711507, strand = '*', stringsAsFactors = FALSE)
  HoxD2 = makeGRangesFromDataFrame(HoxD2, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  Hoxs = c(HoxA, HoxD1, HoxD2)
  
}

##########################################
# Here import design matrix and read counts of pooled peaks across conditions (pval < 10^-6)
# in the future the IDR will be used to select the peaks across replicates and and then pool peaks
##########################################
load(file = paste0("../results/R10723_Rxxxx_R11637_atacseq_R11876_CutTag/Rdata",
                   '/samplesDesign_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical.Rdata'))

##########################################
# Quick check the sample infos and manually add batch information
##########################################
Manual.change.sample.infos = FALSE
if(Manual.change.sample.infos){
  # colnames(counts) = c('gene', paste0(design$fileName, '_', design$SampleID))
  #design$total = unlist(design$total)
  #design$unique.rmdup = unlist(design$unique.rmdup)
  #design$pct.usable = design$unique.rmdup/design$total
  
  design$batch = NA
  design$batch[grep('895|933|749|102', design$SampleID)] = '2020'
  design$batch[grep('136|137', design$SampleID)] = '2021'
  design$batch[grep('161|166', design$SampleID)] = '2021S' # 2021 summer
  
  rownames(counts) = counts$gene
  counts = as.matrix(counts[, -1])
  
  save(counts, design, 
       file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical.Rdata'))
  
}

#hist(design$unique.rmdup/design$total, breaks = 16, main = 'pct of usable reads/total ')

xx = design
xx$unique.rmdup = as.numeric(xx$usable)

#xx$unique.rmdup = xx$unique.rmdup
ggplot(data = xx, aes(x = fileName, y = unique.rmdup, color= condition)) +   
  geom_bar(aes(fill = batch), position = "dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + 
  geom_hline(yintercept = c(20, 50, 100)) + ylab("unique.rmdup (M)")


########################################################
########################################################
# Section 0 : quatitative analysis of peaks: bineary 0 (no peak) and 1 (peak)
# QCs and first analysis based on binary data (peak, no peak)
########################################################
########################################################
Binary.peaks.QCs.analysis = FALSE
if(Binary.peaks.QCs.analysis){
  source('Functions.R')
  ATACseq.peaks.binary.analysis()

}

########################################################
########################################################
# Section I : normalization and batch correction
########################################################
########################################################
load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical.Rdata'))

Normalization.BatchCorrect = FALSE
if(Normalization.BatchCorrect){
  design$conds = design$condition
  design$unique.rmdup = design$usable
  colnames(design)[which(colnames(design) == 'fileName')] = 'samples'
  #sels = grep('Mature|Embryo|BL_UA', design$conds)
  sels = c(1:nrow(design))
  
  dds <- DESeqDataSetFromMatrix(as.matrix(counts[, sels]), DataFrame(design[sels, ]), design = ~ conds)
  
  #rm(counts)
  
  # check the peak length
  peakNames = rownames(dds)
  pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  ll = width(pp)
  
  save(counts, design, 
       file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical_v1.Rdata'))
  ##########################################
  # filter peaks below certain thrshold of read counts
  # And also consider those filtered peaks as background
  ##########################################
  #load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical_v1.Rdata'))
  select.peaks.with.readThreshold = TRUE
  select.background.for.peaks = TRUE
  
  if(select.peaks.with.readThreshold){
    #ss = rowMax(counts(dds)[, grep('Embryo_', dds$conds)])
    ss = rowMaxs(counts(dds))/ll*500
    hist(log10(ss), breaks = 200, main = 'log2(max of read counts within peaks) ')
    cutoff.peak = 30 # 30 as peak cutoff looks good
    cutoff.bg = 10
    cat(length(which(ss >= cutoff.peak)), 'peaks selected with minimum read of the highest peak -- ', cutoff.peak,  '\n')
    cat(length(which(ss < cutoff.bg)), 'peaks selected with minimum read of the highest peak -- ', cutoff.bg,  '\n')
    abline(v= log10(cutoff.peak), col = 'red', lwd = 2.0)
    abline(v= log10(cutoff.bg), col = 'blue', lwd = 2.0)
    
    nb.above.threshold = apply(counts(dds), 1, function(x) length(which(x>cutoff.peak)))
    ii = which(ss >= cutoff.peak)
    #ii = which(nb.above.threshold>=2)
    
    if(select.background.for.peaks){
      ii.bg = which(ss < cutoff.bg)
      ii.bg = sample(ii.bg, size = 1000, replace = FALSE)
      rownames(dds)[ii.bg] = paste0('bg_', rownames(dds)[ii.bg])
      dds = dds[c(ii, ii.bg), ]
      ll.sels = ll[c(ii, ii.bg)]
      
    }else{
      dds <- dds[ii, ]
      ll.sels = ll[ss >= cutoff.peak]
    }
    
  }
  
  dds <- estimateSizeFactors(dds)
  plot(sizeFactors(dds), colSums(counts(dds))/median(colSums(counts(dds))), log = 'xy')
  
  plot(sizeFactors(dds), design$usable, log = 'xy')
  text(sizeFactors(dds), design$usable, labels = design$samples, cex = 0.7)
  
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
  
  ##########################################
  # test normalization and batch correction of ATAC-seq data
  # TMM and combat were selected for normalization and batch correction
  ##########################################
  edgeR.normalization.batch.correction = FALSE
  if(edgeR.normalization.batch.correction){
    source('Functions_atac.R')
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
    fpm.bc = ComBat(dat=fpm, batch=bc, mod=mod, par.prior=TRUE, ref.batch = '2021S') # 2021S as reference is better for some reasons    
    fpm = fpm.bc
    
    # fpm.bc = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
    make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
    #make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'Dev.Mature')
    
    rm(fpm.bc)
    saveRDS(fpm, file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
    
  }
  
}

########################################################
########################################################
# Section II : grouping atac-seq peak profiles
# 1) first identify static peaks (probably most of them) with model selection 
# (linear regression with only intercept and nonlinear fitting with spline or GAM)
# 2) grouping dynamic peaks using DP-GP
# 
########################################################
########################################################
Grouping.atac.peaks = FALSE
if(Grouping.atac.peaks){
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
  # all-peaks test M0 (static peaks), M1 (dynamic peaks above background), M2 (dyanmic peaks with some condtions below background)
  # we will also loci that are always open across all conditions
  ##########################################
  make.test.All.Peaks = FALSE
  if(make.test.All.Peaks){
    
    Run.all.peaks.test = FALSE
    if(Run.all.peaks.test){
      
      # all conditions included
      conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", "Embryo_Stage44_distal",
                "Mature_UA",  "Mature_LA", "Mature_Hand", 
                "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", "BL_UA_13days_distal")
      
      # test examples
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
      
    }else{
      res = readRDS(file = paste0(RdataDir, '/res_allPeaks_test.rds'))
      
    }
    
    jj = which(res$prob.M0 < 0.05 & res$log2FC >1)
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC), ]
    
    length(which(xx$min <3.0))
    length(which(xx$min <2.5))
    length(which(xx$min <2.))
    
    bgs = data.frame(bg = as.numeric(fpm.bg))
    
    ggplot(bgs, aes(x=bg)) + geom_histogram(binwidth = 0.25, color="darkblue", fill="lightblue") +
      theme(axis.text.x = element_text(angle = 0, size = 16)) + 
      xlab(" background signals (log2)") 
    
    Test.promoter.openness = FALSE
    if(Test.promoter.openness){
      
      source('Functions.R')
      Test.promoter.openness.enrichment(res, bg) # two different methods to test nonHS promoters have tendancy to be open across conditions
      
      #DoubleCheck.promoter.peaks.enrichment(fpm) # double check those open promoters are not due to sample contamination
      
    }
    
  }
  
  ##########################################
  # Position-dependent test
  # mainly use the mature samples, mUA, mLA and mHand
  # 
  ##########################################
  grouping.position.dependent.peaks = FALSE
  if(grouping.position.dependent.peaks){
    
    #conds = as.character(unique(design$conds))
    # conds = c("Mature_UA", "Mature_LA", "Mature_Hand", 
    #           "Embryo_Stage44_proximal", "Embryo_Stage44_distal",
    #           "BL_UA_13days_proximal", "BL_UA_13days_distal")
    
    conds = c("Mature_UA", "Mature_LA", "Mature_Hand")
              
    sample.sels = c();  cc = c()
    for(n in 1:length(conds)) {
      #kk = which(design$conds == conds[n] & design$SampleID != '136159')
      kk = which(design$conds == conds[n]) 
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    Run.test.spatial.peaks = FALSE
    if(Run.test.spatial.peaks){
      # examples to test
      #test.examples = c('HAND2', 'FGF8', 'KLF4', 'Gli3', 'Grem1')
      #test.examples = c('Hoxa13')
      #ii.test = which(overlapsAny(pp, promoters[which(!is.na(match(promoters$geneSymbol, test.examples)))]))
      #ii.Hox = which(overlapsAny(pp, Hoxs))
      #ii.test = unique(c(ii.test, ii.Hox))
      
      library(tictoc)
      ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
      source('Functions_atac.R')
      
      # select the mature samples
      cpm = fpm[, sample.sels]
      
      # select the peaks that are above background with >3 samples
      quantile(fpm.bg, 0.999)
      
      hist(fpm.bg)
      
      # stringent here, with signal > 2.5 in >=3 samples 
      nb.above.threshold = apply(as.matrix(cpm), 1, function(x) length(which(x> 2.5)))
      hist(nb.above.threshold, breaks = c(-1:ncol(cpm)))
      peak.sels = which(nb.above.threshold>=3)
      cpm = cpm[peak.sels, ]
      
      tic() 
      res = spatial.peaks.test(cpm = cpm, c = cc, test.Dev.Reg = FALSE)
      res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
      toc()
      
      xx = data.frame(res, pp.annots[match(rownames(res), rownames(pp.annots)), ], stringsAsFactors = FALSE)
      
      res = xx
      saveRDS(res, file = paste0(RdataDir, '/res_position_dependant_test_v6.rds'))
      
    }
    
    ##########################################
    # select all positional-dependent loci with below threshold
    ##########################################
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
    
    keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
    keep = as.matrix(keep)
    
    library(ggplot2)
    df <- data.frame(cc)
    rownames(df) = colnames(keep)
    
    ii.gaps = c(5, 8)
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, 
             filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_rmPeaks.head.pdf'), 
             width = 8, height = 12)
    
    if(saveTable){
      yy = data.frame(keep, xx, stringsAsFactors = FALSE)
      write.csv(yy, file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head.csv'), 
                quote = FALSE, row.names = TRUE)
      
    }
    
    ##########################################
    ## select top peaks genome-wide and also promoter peaks
    ##########################################
    #res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_v2.rds'))
    #jj = which(res$prob.M0.mature < 0.01 & res$log2FC.mature > 3 & res$min.mature <1)
    #xx = res[jj, ]
    #xx = xx[order(-xx$log2FC.mature), ]
    #xx = data.frame(xx, pp.annots[match(rownames(xx), rownames(pp.annots)), ], stringsAsFactors = FALSE)
    xx = xx[order(-xx$logFC.mean), ]
    
    yy = xx
    yy = yy[which(yy$max > 3 & yy$min < 1 & yy$res2.max< 0.1), ]
    
    yy = yy[order(yy$logFC.mean), ]
    
    keep = fpm[!is.na(match(rownames(fpm), rownames(yy))), sample.sels]
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
    ii.gaps = c(5, 8)
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, fontsize_row = 8, gaps_col = ii.gaps,
             filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_top300_genomewide.pdf'), 
             width = 16, height = 40)
    
    if(saveTable){
      write.csv(data.frame(keep, yy, stringsAsFactors = FALSE), 
                file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_top300_genomewide.csv'), 
                quote = FALSE, row.names = TRUE)
      
    }
    
    ##########################################
    # select top promoter peaks
    ##########################################
    yy = xx
    
    yy = yy[grep('Promoter', yy$annotation), ] # peak close to promoters
    yy[grep('HOXA13', yy$transcriptId),]
    #yy = xx
    #yy = yy[c(1:50), ]
    #yy = yy[which(yy$logFC.mean>1.5), ]
    yy = yy[which(yy$max > 3 & yy$min < 2 & yy$res2.max< 0.1), ]
    
    keep = fpm[!is.na(match(rownames(fpm), rownames(yy))), sample.sels]
    gg = res$geneId[match(rownames(keep), rownames(res))]
    grep('HOXA13', gg)
    rownames(keep) = paste0(rownames(keep), '_', gg)
    #rownames(keep) = gg
    keep = as.matrix(keep)
    
    # jj1 = grep('Embryo_Stage40', colnames(keep))
    # jj2 = grep('Embryo_Stage44', colnames(keep))
    # jj3 = grep('Mature_UA', colnames(keep))
    # mean1 = apply(keep[, jj1], 1, mean)
    # mean2 = apply(keep[, jj2], 1, mean)
    # mean3 = apply(keep[, jj3], 1, mean)
    # 
    # kk = which(mean1< 1. & mean2 < 1. & mean3 < 1.)
    # 
    # # make sure max above the threshold
    # nb.above.threshold = apply(keep, 1, function(x) length(which(x>3)))
    # keep = keep[which(nb.above.threshold >=3), ] 
    # 
    # nb.below.bg = apply(keep, 1, function(x) length(which(x<2)))
    # keep = keep[which(nb.below.bg >=3), ] 
    
    gg = rownames(keep)
    gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '_'))[2])
    gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '[|]'))[1])
    #keep = keep[1:50, ]
    #rownames(keep) = gg
    #kk = grep('Mature', cc)
    #df <- data.frame(condition = cc[kk])
    #keep = keep[,kk]
    #rownames(df) = colnames(keep)
    ii.gaps = c(5, 8)
    pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, fontsize_row = 11, gaps_col = ii.gaps,
             filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_top50_promoter.pdf'), 
             width = 16, height = 12)
    
    
    if(saveTable){
      write.csv(data.frame(keep, yy, stringsAsFactors = FALSE), 
                file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_top50_promoterPeaks.csv'), 
                quote = FALSE, row.names = TRUE)
      
    }
    ##########################################
    # first motif activity analysis for positional-dependent peaks 
    ##########################################
    source('MARA_functions.R')
    
    saveRDS(keep, file = paste0(RdataDir, '/matrix.saved.for.spatial.MARA.rds'))
    
    xx = run.MARA.atac.spatial(keep, cc)
   
    
  }
  
  ##########################################
  # temporal-peaks test
  ##########################################
  grouping.temporal.peaks = FALSE
  if(grouping.temporal.peaks){
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
    
    library(tictoc)
    # ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
    
    Run.temporal.peak.test = FALSE
    if(Run.temporal.peak.test)
    {
      source('Functions_atac.R')
      cpm = fpm[, sample.sels]
      
      quantile(fpm.bg, 0.99)
      nb.above.threshold = apply(as.matrix(cpm), 1, function(x) length(which(x> 2)))
      
      hist(nb.above.threshold, breaks = c(-1:ncol(cpm)))
      length(which(nb.above.threshold>6))
      
      peak.sels = which(nb.above.threshold>=3)
      cpm = cpm[peak.sels, ]
      
      tic()
      ## define the dynamic enhancers with mature UA and BL.UA and check them if embryo samples
      sels = grep('Embryo', cc, invert = TRUE) 
      res = t(apply(cpm[, sels], 1, temporal.peaks.test, c = cc[sels]))
      toc()
      
      xx = data.frame(res, pp.annots[match(rownames(cpm), rownames(pp.annots)), ],  stringsAsFactors = FALSE)
     
      res = xx
      
      saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test_v3.rds'))
      
    }
    
    
    res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test_v3.rds'))
    
    # select the temporal dynamic peaks
    length(which(res$prob.M0<0.05))
    length(which(res$prob.M0<0.05 & res$log2FC > 1))
    length(which(res$prob.M0<0.01 & res$log2FC > 1))
    length(which(res$prob.M0<0.01 & res$log2FC > 1.5))
    length(which(res$prob.M0<0.01 & res$log2FC > 2))
    
    #length(which(res$prob.M0<0.001 & res$log2FC > 2))
    
    #jj = which(res$prob.M0 < 0.05 & res$log2FC >1 )
    
    jj = which(res$prob.M0 < 0.01 & res$log2FC > 2 )
    
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC), ]
    xx = xx[which(xx$min < 1), ]
    
    source('Functions_atac.R')
    keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
    keep = as.matrix(keep)
    
    kk = c(grep('Embryo_Stage40', colnames(keep)), 
           grep('Embryo_Stage44', colnames(keep)))
    kk = c(setdiff(c(1:ncol(keep)), kk), kk)
    
    keep = keep[, kk]
    df <- data.frame(cc[kk])
    rownames(df) = colnames(keep)
    
    ii.gaps = c(4, 8, 10, 12, 16)+1
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps,
             filename = paste0(resDir, '/heatmap_regenerationPeaks_fdr0.01_log2FC.1.pdf'), 
             width = 12, height = 12)
    
    ##########################################
    # highligh potential regeneration peaks, not found in mUA and embryo stages only in regeneration process 
    ##########################################
    jj1 = grep('Embryo_Stage40', colnames(keep))
    jj2 = grep('Embryo_Stage44', colnames(keep))
    jj3 = grep('Mature_UA', colnames(keep))
    mean1 = apply(keep[, jj1], 1, mean)
    mean2 = apply(keep[, jj2], 1, mean)
    mean3 = apply(keep[, jj3], 1, mean)
    
    kk = which(mean1< 1. & mean2 < 1. & mean3 < 1.)
    
    pheatmap(keep[kk, ], cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps,
             filename = paste0(resDir, '/heatmap_regenerationPeaks_fdr0.01_log2FC.1_regeneartion.specific.pdf'), 
             width = 8, height = 10)
    
    yy = keep[kk,]
    
    ##########################################
    # first motif activity analysis for temporally dynamic peaks 
    ##########################################
    source('MARA_functions.R')
    xx = run.MARA.atac.temporal(keep, cc)
    
  }
  
}


