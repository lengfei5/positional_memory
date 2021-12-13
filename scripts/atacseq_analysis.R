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

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0(figureDir, 'tables4plots/')


annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')

require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require(pheatmap)
library(tictoc)

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
  
  # check the peak length
  peakNames = rownames(dds)
  pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  ll = width(pp)
  
  ##########################################
  # filter peaks below certain thrshold of read counts
  # And also consider those filtered peaks as background
  ##########################################
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
    
    sfs = data.frame(sample = colnames(dds), sf = sizeFactors(dds)*median(colSums(counts(dds))), stringsAsFactors = FALSE)
    
    saveRDS(sfs, file = paste0(RdataDir, '/DESeq2_peaks.based_scalingFactors_forGenomicRanger.rds'))
    
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
    require("sva")
    require(limma)
    
    ##########################################
    # batch correct samples separately for mature samples and regeneration samples  
    ##########################################
    table(design$condition, design$batch)
    
    Split.Mature.Regeneration.samples = TRUE
    if(Split.Mature.Regeneration.samples){
      
      # start with mature samples
      Batch.Correct.matureSamples = FALSE
      if(Batch.Correct.matureSamples){
        sels = grep('Mature|HEAD', design$conds)
        sels = setdiff(sels, which(design$SampleID == '74938'| design$SampleID == '74939'|design$SampleID == '74940'))
        design.sels = design[sels, ]
        
        design.sels$conds = droplevels(design.sels$conds)
        
        #design.sels$batch[which(design.sels$batch == '2021')] = '2021S'
        
        #design.sels$batch[grep('749', design.sels$SampleID)] = '2019'
        
        #design.sels$batch = droplevels(design.sels$batch)
        table(design.sels$conds, design.sels$batch)
        
        ddx = dds[, sels]
        ddx$conds = droplevels(ddx$conds)
        ss = rowSums(counts(ddx))
        # remove low count genes, otherwise combat returns error 
        # 'Error in while (change > conv) { : missing value where TRUE/FALSE needed'
        ddx = ddx[which(ss>5), ] 
        ddx = estimateSizeFactors(ddx)
        vsd <- varianceStabilizingTransformation(ddx, blind = TRUE)
        tmm = assay(vsd)
        
        #d <- DGEList(counts=counts(ddx), group=design.sels$conds)
        #tmm <- calcNormFactors(d, method='TMM')
        #tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
        
        tmm.vars = apply(as.matrix(tmm), 1, var) # row with var = 0 pose problem for ComBat
        tmm = tmm[which(tmm.vars>0 & !is.na(tmm.vars)), ]
        
        bc = as.factor(design.sels$batch)
        mod = model.matrix(~ as.factor(conds), data = design.sels)
        
        # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
        # because there is no better batche other others 
        #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
        fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
        
        #design.tokeep<-model.matrix(~ 0 + conds,  data = design.sels)
        #cpm.bc = limma::removeBatchEffect(tmm, batch = bc, design = design.tokeep)
        # plot(fpm.bc[,1], tmm[, 1]);abline(0, 1, lwd = 2.0, col = 'red')
        make.pca.plots(fpm.bc, ntop = 1000, conds.plot = 'Mature')
        
        make.pca.plots(tmm, ntop = 1000, conds.plot = 'Mature')
        
        
        fpm = fpm.bc
        
        rm(fpm.bc)
        
        saveRDS(fpm, file = paste0(RdataDir, '/fpm.bc_TMM_combat_MatureSamples_batch2020.2021.2021S_rmOldBatch.rds'))
        saveRDS(design.sels, file = paste0(RdataDir, '/design_sels_bc_TMM_combat_MatureSamples_batch2020.2021.2021S_rmOldBatch.rds'))
        
      }
      # regeneration time points and embryo stages
      Batch.Correct.regeneration.embryoStage = FALSE
      if(Batch.Correct.regeneration.embryoStage){
        
        sels = unique(c(setdiff(which(design$batch == '2020'), grep('Mature', design$condition)), 
                        which(design$batch == '2021'), which(design$condition == 'BL_UA_9days')))
        #sels = sels[which(design$batch[sels] != '2020')]
        design.sels = design[sels, ]
        design.sels$conds = droplevels(design.sels$conds)
        
        table(design.sels$conds, design.sels$batch)
        
        design.sels$batch[which(design.sels$batch == '2021S')] = '2021'
        #design.sels$batch = droplevels(design.sels$batch)
        table(design.sels$conds, design.sels$batch)
        
        ddx = dds[, sels]
        ddx$conds = droplevels(ddx$conds)
        ss = rowSums(counts(ddx))
        # remove low count genes, otherwise combat returns error 
        # 'Error in while (change > conv) { : missing value where TRUE/FALSE needed'
        ddx = ddx[which(ss>5), ] 
        ddx = estimateSizeFactors(ddx)
        vsd <- varianceStabilizingTransformation(ddx, blind = TRUE)
        tmm = assay(vsd)
        
        #d <- DGEList(counts=counts(ddx), group=design.sels$conds)
        #tmm <- calcNormFactors(d, method='TMM')
        #tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
        
        tmm.vars = apply(as.matrix(tmm), 1, var) # row with var = 0 pose problem for ComBat
        tmm = tmm[which(tmm.vars>0 & !is.na(tmm.vars)), ]
        
        bc = as.factor(design.sels$batch)
        mod = model.matrix(~ as.factor(conds), data = design.sels)
        
        # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
        # because there is no better batche other others 
        #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
        fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
        
        #design.tokeep<-model.matrix(~ 0 + conds,  data = design.sels)
        #cpm.bc = limma::removeBatchEffect(tmm, batch = bc, design = design.tokeep)
        # plot(fpm.bc[,1], tmm[, 1]);abline(0, 1, lwd = 2.0, col = 'red')
        make.pca.plots(fpm.bc, ntop = 1000, conds.plot = 'all')
        
        make.pca.plots(tmm, ntop = 1000, conds.plot = 'all')
        
        
        fpm = fpm.bc
        
        rm(fpm.bc)
        
        saveRDS(fpm, file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration_embryoStages.rds'))
        saveRDS(design.sels, file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration_embryoStages.rds'))
          
      }
      
      
    }else{
      d <- DGEList(counts=counts(dds), group=design$conds)
      tmm <- calcNormFactors(d, method='TMM')
      tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
      #rm(tmm)
      rm(d)
      #fpm = log2(tmm + 1)
      make.pca.plots(tmm[,which(design$batch == '2020')], ntop = 1000, conds.plot = 'all')
      
      design$batch[grep('749', design$SampleID)] = '2019'
      
      table(design$condition, design$batch)
      
      bc = as.factor(design$batch)
      mod = model.matrix(~ as.factor(conds), data = design)
      
      fpm.bc = ComBat(dat=tmm, batch=bc, mod=mod, par.prior=FALSE, ref.batch = NULL) # No reference is better for some reasons    
      
      
      xx = tmm[, grep('Mature', colnames(fpm.bc))]
      yy = fpm.bc[, grep('Mature', colnames(fpm.bc))]
      
      # fpm.bc = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
      make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
      make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'Dev.Mature')
      make.pca.plots(fpm.bc, ntop = 1000, conds.plot = 'Mature')
      
      fpm = fpm.bc
      
      rm(fpm.bc)
      saveRDS(fpm, file = paste0(RdataDir, '/fpm_TMM_combat.rds')) 
    }
    
  }
  
}


########################################################
########################################################
# Section II: promoter analysis
# grouping atac-seq peak profiles
# 1) first identify static peaks (probably most of them) with model selection 
# (linear regression with only intercept and nonlinear fitting with spline or GAM)
# 2) grouping dynamic peaks using DP-GP
# ##########################################
# all-peaks test M0 (static peaks), M1 (dynamic peaks above background), M2 (dyanmic peaks with some condtions below background)
# we will also loci that are always open across all conditions
########################################################
########################################################
make.test.All.Peaks = FALSE
if(make.test.All.Peaks){
  
  # to change for all peaks
  fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.rds'))
  design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.rds'))
  
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

########################################################
########################################################
# Section III : positional peaks
# ##########################################
# Position-dependent test
# mainly use the mature samples, mUA, mLA and mHand
# 
##########################################
########################################################
########################################################
grouping.position.dependent.peaks = FALSE
if(grouping.position.dependent.peaks){
  
  ##########################################
  # import batch corrected gene expression and design 
  ##########################################
  fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_MatureSamples_batch2020.2021.2021S_rmOldBatch.rds'))
  design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_MatureSamples_batch2020.2021.2021S_rmOldBatch.rds'))
  
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
  
  ##########################################
  # run spatial test for mature samples
  ##########################################
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
    quantile(fpm.bg, 0.99)
    
    hist(fpm.bg)
    
    # stringent here, with signal > 2.5 in >=3 samples 
    nb.above.threshold = apply(as.matrix(cpm), 1, function(x) length(which(x> 4.5)))
    hist(nb.above.threshold, breaks = c(-1:ncol(cpm)))
    peak.sels = which(nb.above.threshold>=2)
    cat(length(peak.sels), 'peaks after filtering for mature samples\n')
    
    cpm = cpm[peak.sels, ]
    
    tic() 
    res = spatial.peaks.test(cpm = cpm, c = cc, test.Dev.Reg = FALSE)
    #res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
    toc()
    
    xx = data.frame(res, pp.annots[match(rownames(res), rownames(pp.annots)), ], stringsAsFactors = FALSE)
    
    res = xx
    saveRDS(res, file = paste0(RdataDir, '/res_position_dependant_test_v9.rds'))
    
  }
  
  ##########################################
  # select all positional-dependent loci with below threshold
  ##########################################
  res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_v9.rds'))
  
  # select the positional peaks with pairwise comparisions 
  # limma logFC is in log2 scale
  fdr.cutoff = 0.01; logfc.cutoff = 0.5
  jj = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
               (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
               (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)
  )
  cat(length(jj), '\n')
  
  jj1 = which(res$prob.M0<0.01 & res$log2FC>logfc.cutoff)
  cat(length(jj1), '\n')
  jj2 = which(res$pval.lrt < 0.001 & res$log2FC > logfc.cutoff)
  cat(length(jj2), '\n')
  
  xx = res[c(jj), ]
  #xx = xx[order(-xx$log2FC.mature), ]
  xx[grep('HOXA13|SHOX', xx$transcriptId), ]
  fpm[which(rownames(fpm) == 'chr2p:873464923-873465440'), ]
  
  # filter the peaks from head control sample
  Filtering.peaks.in.Head.samples = TRUE
  if(Filtering.peaks.in.Head.samples){
    p0 = pp[match(rownames(xx), names(pp))]
    
    ctl = fpm[, grep('HEAD', colnames(fpm))]
    ctl = ctl[which(ctl>4.5)]
    p.ctl = pp[match(names(ctl), names(pp))]
    
    non.overlap = !overlapsAny(p0, p.ctl)
    
    xx = xx[non.overlap, ]
  }
  
  dim(xx)
  
  xx[grep('HOXA13|SHOX|MEIS', xx$transcriptId), ]
  
  cat(nrow(xx), ' peaks left\n')
  # sort positional peaks with logFC
  #xx = xx[order(-xx$logFC.mean), ]
  xx = xx[order(-xx$log2FC), ]
  
  keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
  keep = as.matrix(keep)
  
  library(ggplot2)
  df <- data.frame(cc)
  rownames(df) = colnames(keep)
  colnames(df) = 'segments'
  ii.gaps = c(4, 6)
  
  annot_colors = c('springgreen4', 'steelblue2', 'gold2')
  names(annot_colors) = c('Mature_UA', 'Mature_LA', 'Mature_Hand')
  annot_colors = list(segments = annot_colors)
  
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, 
           annotation_colors = annot_colors, 
           filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC0.5_rmPeaks.head.pdf'), 
           width = 8, height = 10)
  
  if(saveTable){
    yy = data.frame(keep, xx, stringsAsFactors = FALSE)
    write.csv(yy, file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head.csv'), 
              quote = FALSE, row.names = TRUE)
    
  }
  
  ##########################################
  ## select top peaks genome-wide
  ##########################################
  Select.top.peaks = FALSE
  if(Select.top.peaks){
    #res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_v2.rds'))
    #jj = which(res$prob.M0.mature < 0.01 & res$log2FC.mature > 3 & res$min.mature <1)
    #xx = res[jj, ]
    #xx = xx[order(-xx$log2FC.mature), ]
    #xx = data.frame(xx, pp.annots[match(rownames(xx), rownames(pp.annots)), ], stringsAsFactors = FALSE)
    xx = xx[order(-xx$logFC.mean), ]
    
    yy = xx
    yy[grep('HOXA13', yy$geneId), ]
    
    # res2 here is residues in the fitting for each condition, UA, LA and Hand
    yy = yy[which(yy$max > 5 & yy$min < 4 & yy$res2.max< 0.1), ]
    #yy = yy[which(yy$max > 3 & yy$min < 1), ]
    
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
    ii.gaps = c(4, 6)
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, fontsize_row = 8, gaps_col = ii.gaps,
             filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_top300_genomewide.pdf'), 
             width = 16, height = 40)
    
    if(saveTable){
      write.csv(data.frame(keep, yy, stringsAsFactors = FALSE), 
                file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_top300_genomewide.csv'), 
                quote = FALSE, row.names = TRUE)
      
    }
  }
  
  ##########################################
  # select top promoter peaks
  ##########################################
  yy = xx
  
  yy = yy[grep('Promoter', yy$annotation), ] # peak close to promoters
  yy[grep('HOXA13|MEIS|SHOX|HOXC', yy$transcriptId),]
  #yy = xx
  #yy = yy[c(1:50), ]
  #yy = yy[which(yy$logFC.mean>1.5), ]
  yy = yy[which(yy$max > 5), ] # max residual square < 0.1
  
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
  gg = rownames(keep)
  gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '_'))[2])
  gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '[|]'))[1])
  #keep = keep[1:50, ]
  #rownames(keep) = gg
  #kk = grep('Mature', cc)
  #df <- data.frame(condition = cc[kk])
  #keep = keep[,kk]
  #rownames(df) = colnames(keep)
  ii.gaps = c(4, 6)
  pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = df, fontsize_row = 14, gaps_col = ii.gaps,
           annotation_colors = annot_colors,
           filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC0.5_top.promoters.pdf'), 
           width = 18, height = 10)
  
  
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
  
  # prepare step is in the MARA_functions.R
  xx = run.MARA.atac.spatial(keep, cc)
  
}

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
  
  
  # ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
  
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
  
  #xx = xx[which(xx$min < 4), ]
  
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

