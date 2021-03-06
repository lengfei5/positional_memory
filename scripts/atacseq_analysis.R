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

gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')

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

xx = design
xx$unique.rmdup = xx$unique.rmdup/10^6
ggplot(data = xx, aes(x = samples, y = unique.rmdup, color= conds)) +   
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
Normalization.BatchCorrect = FALSE
if(Normalization.BatchCorrect){
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
    
    dds <- estimateSizeFactors(dds)
      
  }
  
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
  load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.withinPeaks.pval6.Rdata'))
  fpm = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
  
  # prepare the background distribution
  fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
  fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
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
    pp.annots = as.data.frame(annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript'))
    
  }
  
  promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
  
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
  # position-dependent test
  ##########################################
  grouping.position.dependent.peaks = FALSE
  if(grouping.position.dependent.peaks){
    
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
    
    Run.test.spatial.peaks = FALSE
    if(Run.test.spatial.peaks){
      # examples to test
      test.examples = c('HAND2', 'FGF8', 'KLF4', 'Gli3', 'Grem1')
      #test.examples = c('Hoxa13')
      ii.test = which(overlapsAny(pp, promoters[which(!is.na(match(promoters$geneSymbol, test.examples)))]))
      #ii.Hox = which(overlapsAny(pp, Hoxs))
      #ii.test = unique(c(ii.test, ii.Hox))
      
      library(tictoc)
      ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
      source('Functions.R')
      tic() 
      res = t(apply(fpm[ii.test, sample.sels], 1, spatial.peaks.test, c = cc, test.Dev.Reg = FALSE))
      res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
      toc()
      
      saveRDS(res, file = paste0(RdataDir, '/res_position_dependant_test_v2.rds'))
      
    }
    ##########################################
    # select all positional-dependent loci with below threshold
    ##########################################
    res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_v2.rds'))
    
    # select the spatially dynamic peaks
    jj = which(res$prob.M0.mature < 0.01 & res$log2FC.mature > 1 )
    
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC.mature), ]
    
    length(which(xx$min.mature <1.))
    length(which(xx$min.mature <2.))
    length(which(xx$min.mature <2.5))
    length(which(xx$min.mature < 3.0))
    
    keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
    keep = as.matrix(keep)
    
    library(ggplot2)
    df <- data.frame(cc)
    rownames(df) = colnames(keep)
    
    ii.gaps = c(7, 11)
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    
    ##########################################
    # # select loci that closes at some conditions
    ##########################################
    res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_v2.rds'))
    jj = which(res$prob.M0.mature < 0.01 & res$log2FC.mature > 3 & res$min.mature <1)
    xx = res[jj, ]
    #xx = xx[order(-xx$log2FC.mature), ]
    keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
    gg = res$geneId[match(rownames(keep), rownames(res))]
    grep('HOXA13', gg)
    rownames(keep) = paste0(rownames(keep), '_', gg)
    keep = as.matrix(keep)
    
    kk = grep('Mature', cc)
    df <- data.frame(condition = cc[kk])
    keep = keep[,kk]
    rownames(df) = colnames(keep)
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, fontsize_row = 9)
    
    
    ##########################################
    # first motif activity analysis for positional-dependent peaks 
    ##########################################
    source('MARA_functions.R')
    res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_v2.rds'))
    
    # select the spatially dynamic peaks
    jj = which(res$prob.M0.mature < 0.01 & res$log2FC.mature > 2)
    cat(length(jj), ' peaks selected \n')
    
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC.mature), ]
    
    keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), ]
    keep = as.matrix(keep)
    
    conds = c("Mature_UA", "Mature_LA", "Mature_Hand")
    
    sample.sels = c()
    cc = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n] & design$SampleID != '136159')
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    keep = keep[ , sample.sels]
    
    xx = run.MARA.atac.spatial(keep, cc)
   
    
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
    
    Run.temporal.peak.test = FALSE
    if(Run.temporal.peak.test)
    {
      source('Functions.R')
      tic() 
      res = t(apply(fpm[ii.test, sample.sels], 1, temporal.peaks.test, c = cc))
      res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
      toc()
      
      saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test.rds'))
      
    }
    
    res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test.rds'))
    # select the temporal dynamic peaks
    length(which(res$prob.M0<0.05))
    length(which(res$prob.M0<0.05 & res$log2FC > 1))
    jj = which(res$prob.M0 < 0.05 & res$log2FC >1 )
    
    jj = which(res$prob.M0 < 0.01 & res$log2FC > 2 )
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
    source('MARA_functions.R')
    xx = run.MARA.atac.temporal(keep, cc)
    
  }
  
}


