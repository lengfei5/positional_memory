##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: process TSS and merge with regeneration TSS
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon May 16 10:38:58 2022
##########################################################################
##########################################################################

########################################################
########################################################
# Section : first process TSS from mature samples 
# batch correction
########################################################
########################################################
process.normalize.atac.histM.allTSS.matureSamples = function()
{
  RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
  RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
  source(RNA.functions)
  source(RNA.QC.functions)
  
  ## import histM sample info of mature samples
  design = readRDS(file = paste0('../results/CT_merged_20220328/Rdata/histM_CT_design_info.rds'))
  #design = read.csv(file = paste0(dataDir, "R11876_R12810_R12965_CT_analysis_20220217_QCs_AK.csv"))
  #design = design[!is.na(design$sampleID), ]
  design$fileName = paste0(design$condition, '_', design$sampleID)
  design = design[, c(1:3, 15, 5, 16, 4, 6:14)]
  colnames(design)[5] = 'batch'
  
  design.histM = design[, c(1:5)]
  rm(design)
  design.histM = design.histM[grep('mRep', design.histM$batch), ] # select only mature samples in mRep1/2 batches
  
  
  ## import atac-seq mature samples
  design = readRDS(file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                                '/design_sels_bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.2022.rds'))
  
  design = design[, c(1, 2, 6)]
  design$sample = design$condition
  design$marks = 'atac'
  design = design[, c(1, 4,5, 2:3)]
  
  # design matrix of atac and histM for mature samples
  colnames(design.histM) = colnames(design)
  design = rbind(design, design.histM)
  
  # design$sample = gsub('BL_UA_13days_distal', 'BL13days.dist', design$sample)
  # design$sample = gsub('BL_UA_13days_proximal', 'BL13days.prox', design$sample)
  # design$sample = gsub('BL_UA_9days', 'BL9days', design$sample)
  # design$sample = gsub('BL_UA_5days', 'BL5days', design$sample)
  design = design[grep('IgG', design$marks, invert = TRUE), ]
  
  design$sample = gsub('Mature_UA', 'mUA', design$sample)
  design$sample = gsub('Mature_LA', 'mLA', design$sample)
  design$sample = gsub('Mature_Hand', 'mHand', design$sample)
  
  #design = design[grep('Embryo', design$sample, invert = TRUE), ]
  design$condition = paste0(design$marks, '_', design$sample)
  
  saveRDS(design, file = paste0(RdataDir, '/design_sampleInfos_TSS_matureSamples.rds'))
  
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_histM_both/'
  xlist = list.files(path=paste0(dataDir, 'featurecounts_peaks.Q30'), pattern = 'featureCounts.txt.summary', 
                     full.names = TRUE)
  
  design$total.reads = NA
  for(n in 1:nrow(design))
  {
    kk = grep(design$SampleID[n], xlist)
    cat(length(kk), '\n')
    test = read.table(file = xlist[kk], sep = '\t', header = TRUE)
    design$total.reads[n] = sum(as.numeric(as.character(test[,2])))
    
  }
  
  saveRDS(design, file = paste0(RdataDir, '/atac_histM_sample_design_matureSamples.rds'))
  
  ## process read counts for selected samples
  design = readRDS(file = paste0(RdataDir, '/atac_histM_sample_design_matureSamples.rds'))
  
  xlist<-list.files(path=paste0(dataDir, 'featurecounts_peaks.Q30'),
                    pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge
  
  sels = c()
  for(n in 1:nrow(design)) sels = c(sels, grep(design$SampleID[n], xlist))
  xlist = xlist[sels]
  
  all = cat.countTable(xlist, countsfrom = 'featureCounts')
  
  colnames(design)[1] = 'SampleID'
  
  counts = process.countTable(all=all, design = design[, c(1,4)])
  
  save(design, counts, file = paste0(RdataDir, '/atac_histM_sample_design_regeneration_counts_TSS_matureSamples.Rdata'))
  
  ## select the tss if multiple ones for one gene using the regeneration TSS
  load(file = paste0(RdataDir, '/atac_histM_sample_design_regeneration_counts_TSS_matureSamples.Rdata'))
  rownames(counts) = counts$gene
  
  design.matures = design
  load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration_counts_withEmbryo.Rdata')) # tss to use
  
  mm = match(tss$transcript, rownames(counts))
  xx = tss # save the tss from regeneration
  tss = counts[mm, ] # same number of TSS as regeneration TSS
  tss = data.frame(tss[, -1], geneID = xx$geneID, transcript = tss$gene, coords = xx$coords, stringsAsFactors = FALSE)
  rm(design)
  rm(xx)
  design = design.matures
  rm(design.matures)
  
  save(tss, design, file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_matureSamples.Rdata'))
  
  ##########################################
  # select genes to consider by removing non-expressed genes in limb fibroblast 
  ##########################################
  load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_matureSamples.Rdata')) 
  design.matures = design
  tss.matures = tss
  
  # import design and tss of regeneration data
  load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration.embryo_counts_expressedGenes_controls.Rdata'))
  mm = match(tss$coords, tss.matures$coords)
  tss.matures = tss.matures[mm, ]
  
  tss = tss.matures
  design = design.matures
  
  save(tss, design, 
       file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_matureSamples_counts_expressedGenes_controls.Rdata'))
  
  
  ##########################################
  # batch correct TSS for atac data and histM have already majority of TSS considered here 
  ##########################################
  Check.overlaps.between.tss.with.atacseq.And.BatchCorrection = FALSE
  if(Check.overlaps.between.tss.with.atacseq.And.BatchCorrection){
    load( file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_matureSamples_counts_expressedGenes_controls.Rdata'))
    
    tp = data.frame(t(sapply(tss$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    tp$strand = '*'
    tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    load(file = paste0(RdataDir, '/regeneration_samples_beforeBatchCorrection.Rdata'))
    
    pp = data.frame(t(sapply(rownames(ddx), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    
    mapping = findOverlaps(tp, pp, ignore.strand=TRUE,  minoverlap=100L)
    jj = (unique(mapping@from))
    missed = setdiff(c(1:nrow(tss)), jj)
    
    cat(length(missed), ' tss missed \n')
    
    ## add missing tss to regeneration atac data
    load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration.embryo_counts_expressedGenes_controls.Rdata'))
    load(file = paste0(RdataDir, '/regeneration_samples_beforeBatchCorrection.Rdata'))
    
    counts = counts(ddx)
    #kk = grep('Embryo_', design.sels$condition, invert = TRUE)
    #design.sels = design.sels[kk, ]
    #counts = counts[, kk]
    
    kk2 = grep('atac', design$marks)
    raw = tss[missed, kk2]
    design.atac = design[kk2, ]
    
    rownames(raw) = tss$coords[missed]
    
    mm = match(design.sels$sampleID, design.atac$SampleID)
    design.atac = design.atac[mm, ]
    raw = raw[ ,mm]
    
    colnames(raw) = colnames(counts)
    counts = rbind(counts, raw)
    
    ## batch correction for those 69k loci 
    table(design.sels$conds, design.sels$batch)
    
    d <- DGEList(counts=counts, group=design.sels$conds)
    tmm <- calcNormFactors(d, method='TMM')
    tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
    
    bc = as.factor(design.sels$batch)
    mod = model.matrix(~ as.factor(conds), data = design.sels)
    
    # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
    # because there is no better batche other others 
    #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
    fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
    
    #design.tokeep<-model.matrix(~ 0 + conds,  data = design.sels)
    #cpm.bc = limma::removeBatchEffect(tmm, batch = bc, design = design.tokeep)
    # plot(fpm.bc[,1], tmm[, 1]);abline(0, 1, lwd = 2.0, col = 'red')
    
    make.pca.plots(tmm, ntop = 3000, conds.plot = 'all')
    ggsave(paste0(resDir, "/regeneration_embryo_TSS_Samples_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
    ggsave(paste0(resDir, "/regeneration_embryo_TSS_Samples_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    fpm = fpm.bc
    
    rm(fpm.bc)
    
    saveRDS(fpm, file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration.dev.TSS_2Batches.R10723_R7977_', version.analysis, '.rds'))
    saveRDS(design.sels, 
            file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration.dev.TSS_2Batches.R10723_R7977',
                          version.analysis, '.rds'))
    
    ## dynamic peak analysis
    fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration.dev.TSS_2Batches.R10723_R7977_', version.analysis, '.rds'))
    design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration.dev.TSS_2Batches.R10723_R7977',
                                   version.analysis, '.rds'))
    
    # prepare the background distribution
    fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
    fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
    # rownames(fpm) = gsub('_', '-', rownames(fpm))
    
    
    conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", 'Embryo_Stage44_distal', 
              "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", 'BL_UA_13days_distal')
    
    sample.sels = c(); cc = c()
    sample.means = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n])
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
      sample.means = cbind(sample.means, apply(fpm[, kk], 1, mean))
    }
    colnames(sample.means) = conds
    
    source('Functions_atac.R')
    cpm = fpm[, sample.sels]
    
    source('Functions_atac.R')
    tic()
    ## define the dynamic enhancers with mature UA and BL.UA and check them if embryo samples
    #sels = grep('Embryo', cc, invert = TRUE) 
    res = temporal.peaks.test(cpm, c = cc)
    
    toc()
    
    res = data.frame(fpm, res, stringsAsFactors = FALSE)
    
    saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks_regeneration.dev.TSS_2Batches.R10723_R7977_peakAnnot_v1.rds'))
    
    
  }
  
  
  
}
