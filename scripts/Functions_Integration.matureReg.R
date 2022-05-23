##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: mainly for integration of mature and regeneration data and process TSS and merge with regeneration TSS
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
  # batch correct TSS in atac data,  histM have already majority of TSS considered here 
  ##########################################
  Check.overlaps.between.tss.with.atacseq.And.BatchCorrection = FALSE
  if(Check.overlaps.between.tss.with.atacseq.And.BatchCorrection){
    load( file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_matureSamples_counts_expressedGenes_controls.Rdata'))
    
    tp = data.frame(t(sapply(tss$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    tp$strand = '*'
    tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    load(file = paste0(RdataDir, '/atac_matureSamples_beforeBatchCorrection.Rdata'))
    
    pp = data.frame(t(sapply(rownames(ddx), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    
    mapping = findOverlaps(tp, pp, ignore.strand=TRUE,  minoverlap=100L)
    jj = (unique(mapping@from))
    missed = setdiff(c(1:nrow(tss)), jj)
    
    cat(length(missed), ' tss missed \n')
    
    ## add missing tss to regeneration atac data
    load( file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_matureSamples_counts_expressedGenes_controls.Rdata'))
    load(file = paste0(RdataDir, '/atac_matureSamples_beforeBatchCorrection.Rdata'))
    
    counts = counts(ddx)
    #kk = grep('Embryo_', design.sels$condition, invert = TRUE)
    #design.sels = design.sels[kk, ]
    #counts = counts[, kk]
    
    kk2 = grep('atac', design$marks)
    raw = tss[missed, kk2]
    design.atac = design[kk2, ]
    
    # rownames(raw) = tss$coords[missed]
    
    mm = match(design.sels$sampleID, design.atac$SampleID)
    design.atac = design.atac[mm, ]
    raw = raw[ ,mm]
    
    colnames(raw) = colnames(counts)
    counts = rbind(counts, raw)
    
    ## batch correction for those 69k loci 
    require(edgeR)
    require(sva)
    source('Functions_atac.R')
    
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
    ggsave(paste0(resDir, "/matureSamples_TSS_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
    ggsave(paste0(resDir, "/matureSamples_TSS_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    fpm = fpm.bc
    
    rm(fpm.bc)
    
    saveRDS(fpm, file = paste0(RdataDir, '/fpm.bc_TMM_combat_matureSamples.TSS_2Batches.R10723_R7977_', version.analysis, '.rds'))
    saveRDS(design.sels, 
            file = paste0(RdataDir, '/design_sels_bc_TMM_combat_matureSamples.TSS_2Batches.R10723_R7977',
                          version.analysis, '.rds'))
    
    ##########################################
    # ## segment-specific peak analysis
    ##########################################
    fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_matureSamples.TSS_2Batches.R10723_R7977_', version.analysis, '.rds'))
    design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_matureSamples.TSS_2Batches.R10723_R7977',
                                   version.analysis, '.rds'))
    
    # prepare the background distribution
    fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
    fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
    rownames(fpm) = gsub('_', '-', rownames(fpm))
    
    conds = c("Mature_UA", "Mature_LA", "Mature_Hand")
    
    sample.sels = c();  cc = c()
    for(n in 1:length(conds)) {
      #kk = which(design$conds == conds[n] & design$SampleID != '136159')
      kk = which(design$conds == conds[n]) 
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    library(tictoc)
    ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
    source('Functions_atac.R')
    
    cpm = fpm[, sample.sels]
    
    source('Functions_atac.R')
    
    tic() 
    res = spatial.peaks.test(cpm = cpm, c = cc, test.Dev.Reg = FALSE)
    #res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
    toc()
    
    res = data.frame(fpm, res, stringsAsFactors = FALSE)
    
    saveRDS(res, file = paste0(RdataDir, '/res_segment_specific_atac_matureSamples.TSS_2Batches.R10723_R797.rds'))
    
  }
}


Add.TSS.chromatinFeatures.in.matureSamples = function()
{
  source('Functions_histM.R')
  
  ## load tss with already regeneration data
  tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM_geneCorrection_v3.rds'))
  tss$gene[which(rownames(tss) == 'AMEX60DD028208')] = 'PROD1'
  tss$gene[which(rownames(tss) == 'AMEX60DD024424')] = NA
  
  ## add analysis results
  #rna = readRDS(file = paste0("../results/RNAseq_data_used/Rdata/", 
  #                            'regeneration_dynamicGeneClusters_allGenes.rds'))
  #rna$geneID = get_geneID(rownames(rna))
  #tss[, c(8:12)] = rna[match(tss$geneID, rna$geneID), c(1:5)]
  #colnames(tss)[8:12] = paste0('smartseq2_', colnames(tss)[8:12])
  #annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
  #                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  # tss$gene = annot$gene.symbol.toUse[match(tss$geneID, annot$geneID)]
  tp = data.frame(t(sapply(tss$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  tp$strand = '*'
  tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  ########
  ## start to add atac peaks overlapping tss
  ## and the histM overlapping the tss
  ########
  aa = readRDS(file = paste0(RdataDir, '/res_segment_specific_atac_matureSamples.TSS_2Batches.R10723_R797.rds'))
  names = rownames(aa)
  names = gsub('bg_', '', names)
  jj = grep('^chr', names, invert = TRUE)
  names[jj] = tss$coords[match(names[jj], tss$transcript)]
  names.uniq = unique(names)
  names.uniq = names.uniq[!is.na(names.uniq)]
  
  aa = aa[match(names.uniq, names), ]
  aa = data.frame(aa, stringsAsFactors = FALSE)
  rownames(aa) = as.character(names.uniq)
  
  mat = aa[, c(1:14)]
  aa = aa[, -c(1:14)]
  
  # mat = mat[, grep('Embryo_', colnames(mat), invert = TRUE)]
  mat = cal_sample_means(mat, conds = c('Mature_UA', 'Mature_LA', 'Mature_Hand', 'HEAD'))
  newcc = c('mUA', 'mLA', 'mHand', 'Head')
  colnames(mat) = newcc
  
  pp = data.frame(t(sapply(rownames(aa), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(tp, pp, ignore.strand=TRUE,  minoverlap=100L)
  jj = (unique(mapping@from))
  missed = setdiff(c(1:nrow(tss)), jj)
  
  mapping = data.frame(mapping) # mapping from gene to 
  
  # jj_sel = mapping$subjectHits[match(c(1:nrow(tss)), mapping$queryHits)] # randomly select one atac-seq peaks
  jj_sels = c() # if multiple atac peaks overlapping with the tss considered, the one with max signals will be picked up
  for(n in 1:nrow(tss))
  {
    kk = mapping$subjectHits[which(mapping$queryHits == n)]
    if(length(kk) == 0){
      cat(n, '-- missing \n')
      jj_sels = c(jj_sels, NA)
    }else{
      if(length(kk) == 1){
        jj_sels = c(jj_sels, kk)
      }else{
        ss = apply(mat[kk, ], 1, mean)
        jj_sels = c(jj_sels, kk[which.max(ss)])
      }
    }
  }
  
  res = data.frame(mat[jj_sels, ], aa[jj_sels, ], stringsAsFactors = FALSE)
  
  fdr.cutoff = 0.05; logfc.cutoff = 1
  select = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
                        (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
                        (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)
  )
  
  cat(length(select), 'DE peaks found !\n')
  res$dynamic = NA
  res$dynamic[select] = 1
  
  colnames(res) = paste0('atac.M_', colnames(res))
  
  tss =  data.frame(tss, res, stringsAsFactors = FALSE)
  
  ## add histM analysis results
  keep = readRDS(file = paste0('../results/CT_merged_20220328/Rdata/matureSamples_combined_4histMarkers_DE_345k.Rdata'))
  
  peakNames = rownames(keep)
  peakNames = gsub('tss.', '', peakNames)
  pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(tp, pp, ignore.strand=TRUE,  minoverlap=100L)
  jj = (unique(mapping@from))
  missed = setdiff(c(1:nrow(tss)), jj)
  mapping = data.frame(mapping)
  
  #jj_sel = mapping$subjectHits[match(c(1:nrow(tss)), mapping$queryHits)]
  jj_sels = c() # if multiple atac peaks overlapping with the tss considered, the one with max H3K4me3 signals will be picked up
  for(n in 1:nrow(tss))
  {
    kk = mapping$subjectHits[which(mapping$queryHits == n)]
    if(length(kk) == 0){
      cat(n, '-- missing \n')
      jj_sels = c(jj_sels, NA)
    }else{
      if(length(kk) == 1){
        jj_sels = c(jj_sels, kk)
      }else{
        ss = apply(keep[kk, grep('H3K4me3_mUA|H3K4me3_mHand|H3K4me3_LA', colnames(keep))], 1, mean)
        jj_sels = c(jj_sels, kk[which.max(ss)])
      }
    }
  }
  
  res = keep[jj_sels,]
  
  tss =  data.frame(tss, res, stringsAsFactors = FALSE)
  
  saveRDS(tss, file = paste0(RdataDir, '/regeneration_matureSamples_tss_perGene_smartseq2_atac_histM_v4.rds'))
  
}


##########################################
# collect putative enhancer regions  
##########################################
find_enhancers_integeticRegions_Introns_H3K4me1 = function()
{
  source('Functions_histM.R')
  
  samples = c('mHand', 'mLA', 'mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
  
  ### enhancer regions must be open at least one condition, meaning overlapping with 55k atac-seq peaks
  ### so the total nubmer of candidates are 55 regions
  atac = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rxxxx_R10723_R11637_R12810_atac/Rdata/',
                                        'ATACseq_peak_consensus_filtered_55k.rds'))
  enhancers = data.frame(coords = names(atac))
  rownames(enhancers) = enhancers$coords
  
  ## add positional atac-seq data
  fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.2022.rds'))
  #fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
  fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  mm = match(rownames(enhancers), rownames(fpm))
  
  fpm = fpm[mm, ]
  fpm.mean = cal_sample_means(fpm, conds = c('Mature_Hand', 'Mature_LA', 'Mature_UA'))
  colnames(fpm.mean) = paste0('atac_', samples[1:3], '_mRep')
  
  enhancers = data.frame(enhancers, fpm.mean, stringsAsFactors = FALSE)  
  
  ## add regeneration atac data
  fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration_dev_2Batches.R10723_R7977_', version.analysis, '.rds'))
  fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  mm = match(rownames(enhancers), rownames(fpm))
  
  fpm = fpm[mm, ]
  fpm.mean = cal_sample_means(fpm, conds = c('Mature_UA', 'BL_UA_5days', 'BL_UA_9days', 'BL_UA_13days_proximal', 'BL_UA_13days_distal'))
  colnames(fpm.mean) = paste0('atac_', samples[3:7], '_rRep')
  
  enhancers = data.frame(enhancers, fpm.mean, stringsAsFactors = FALSE)  
  
  
  ## add histM feature by features
  #load(file = paste0('../results/CT_merged_20220328/Rdata', 
  #                   '/combined_4histMarkers_overlapped55kATACseq_DE_fdr0.05.Rdata'))
  for(n_histM in 2:length(features))
  {
    # n_histM = 2
    cat('--', features[n_histM], '--\n')
    
    cpm = readRDS(file = paste0('../results/CT_merged_20220328/Rdata', 
                                '/fpm_bc_TMM_combat_', features[n_histM], '_', 'CT_merged_20220328', '.rds'))
    
    conds = c("mUA", "mLA", "mHand")
    sm1 = cal_sample_means(cpm[, grep('_mRep', colnames(cpm))], conds = conds)
    colnames(sm1) = paste0(features[n_histM], '_', samples[1:3], '_mRep')
    
    conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')
    sm2 = cal_sample_means(cpm[, grep('_rRep', colnames(cpm))], conds = conds)
    colnames(sm2) = paste0(features[n_histM], '_', samples[3:7], '_rRep')
    sms = cbind(sm1, sm2)
    
    if(n_histM == 2){
      keep = data.frame(sms, stringsAsFactors = FALSE)
    }else{
      keep = data.frame(keep, sms[match(rownames(keep), rownames(sms)), ], stringsAsFactors = FALSE)
    }
  }
  
  saveRDS(keep, file = paste0(RdataDir, '/histM_samleMean_345k_loci.rds'))
  
  ## relate atac peaks to histM 
  tp = data.frame(t(sapply(enhancers$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  tp$strand = '*'
  tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  pp = gsub('tss.', '', rownames(keep))
  pp = data.frame(t(sapply(rownames(keep), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  rownames(pp) = rownames(keep)
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = GenomicRanges::findOverlaps(tp, pp, ignore.strand=TRUE)
  mapping = data.frame(mapping)
  
  index_enhancer = rep(NA, nrow(enhancers))
  for(n in 1:length(index_enhancer))
  {
    # n = 1
    cat(n, '\n')
    kk = mapping$subjectHits[which(mapping$queryHits == n)]
    if(length(kk) == 1){
      index_enhancer[n] = kk
    }
    if(length(kk) >1){
      stop()
    }
  }
  
  enhancers = data.frame(enhancers, keep[index_enhancer, ], histM_coords = rownames(keep)[index_enhancer],  stringsAsFactors = FALSE)
  
  enhancers = data.frame(enhancers, histM_coords = rownames(keep)[index_enhancer],  stringsAsFactors = FALSE)
  saveRDS(enhancers, file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM.rds'))
  
  enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM.rds'))
  
  # save also length information
  tp = data.frame(t(sapply(enhancers$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  tp$strand = '*'
  tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  enhancers$atac_peak_width = width(tp)
  
  tp = as.character((enhancers$histM_coords))
  jj = which(!is.na(tp))
  tp = tp[jj]
  tp = t(sapply(tp, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':'))))
  tp = data.frame((tp))
  tp$strand = '*'
  tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  enhancers$histM_width = NA
  enhancers$histM_width[jj] = width(tp)
  
  test = as.matrix(enhancers[, grep('H3K4me1', colnames(enhancers))])
  ss = apply(test, 1, max)
  ss = ss + log2(10^3/enhancers$histM_width)
  
  enhancers$rpkm_H3K4me1 = ss
  enhancers$enhancer = NA
  enhancers$enhancer[which(ss>1)] = 1
  
  saveRDS(enhancers, file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1.rds'))
  
  enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1.rds'))
  
  ### peak annotation and then we know if they are integeic/exon
  tp = data.frame(t(sapply(enhancers$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  tp$strand = '*'
  tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  # limb fibroblast expressing genes
  gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
  amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
  
  pp.annots = annotatePeak(tp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  pp.annots = data.frame(pp.annots, stringsAsFactors = FALSE)
  
  # shorten the annotation
  pp.annots = pp.annots[, c(6:14)]
  xx = pp.annots
  xx$annotation[grep('Intron', xx$annotation)] = 'Intron'
  xx$annotation[grep('Promoter', xx$annotation)] = 'Promoter'
  xx$annotation[grep('Exon', xx$annotation)] = 'Exon'
  xx$annotation[grep('Downstream', xx$annotation)] = 'Downstream'
  colnames(xx) = paste0(colnames(xx), '_chipseeker')
  
  enhancers = data.frame(enhancers, xx, stringsAsFactors = FALSE)
  
  rm(pp.annots)
  
  saveRDS(enhancers, file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot.rds'))
  
  cat(length(which(enhancers$enhancer == 1 & 
                 (enhancers$annotation_chipseeker == 'Distal Intergenic'| enhancers$annotation_chipseeker == 'Intron'))),
       ' enhancer candidates \n')
  
  
}

########################################################
########################################################
# Section : integration analysis for gene examples
# 
########################################################
########################################################
plot_rna_chromainFeatures_geneExamples = function(tss, 
                                                  geneList = c('FGF10', 'SALL4'), 
                                                  outDir = 'Gene_Examepls',
                                                  incl_Mature = FALSE,
                                                  log2fc = FALSE)
{
  source('Functions_histM.R')
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library("cowplot")
  require(gridExtra)
  library(tidyr)
  require(patchwork)
  
  # outDir = "/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/positional_genes"   
  # geneList = dev.genes
  
  features = c('atac', 'H3K4me3', 'H3K27me3', 'H3K4me1')
  
  for(n in 1:length(geneList))
  {
    # n = 1
    kk = which(tss$gene == geneList[n])
    if(length(kk) == 1){
      
      cat(n, ' -- ', geneList[n], '-- with tss :', tss$coords[kk], '\n')
      
      if(incl_Mature){
        samples = c('mHand', 'mLA', 'mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
        test = matrix(NA, nrow = length(samples), ncol = length(features))
        rownames(test) = samples
        colnames(test) = features
        
        for(m in 1:ncol(test))
        {
          if(m == 1) {
            mm =c(grep(paste0(colnames(test)[m], '.M_mHand$'), colnames(tss)),  
                  grep(paste0(colnames(test)[m], '.M_mLA$'), colnames(tss)),
                  grep(paste0(colnames(test)[m], '.M_mUA$'), colnames(tss)),
                  grep(paste0(colnames(test)[m], '_mUA|', colnames(test)[m], '_X'), colnames(tss)))
            if(log2fc){
              test[,m] = c(as.numeric(tss[kk, mm[1:2]]) - as.numeric(tss[kk, mm[3]]), 0, 
                           as.numeric(tss[kk, mm[5:8]])- as.numeric(tss[kk, mm[4]]))
            }else{
              test[,m] = c(as.numeric(tss[kk, mm[1:2]]), mean(as.numeric(tss[kk, mm[3:4]])), as.numeric(tss[kk, mm[5:8]]))
            }
            
          }
          if(m > 1){
            mm0 = grep(paste0(colnames(test)[m], '_mUA'), colnames(tss))
            mm0 = mm0[grep('Rep1$|Rep2$', colnames(tss)[mm0])]
            mm = c(grep(paste0(colnames(test)[m], '_mHand'), colnames(tss)),
                   grep(paste0(colnames(test)[m], '_mLA'), colnames(tss)),
                   mm0[grep('mRep', colnames(tss)[mm0])], 
                   mm0[grep('rRep', colnames(tss)[mm0])], 
                   grep(paste0(colnames(test)[m], '_BL'), colnames(tss))
            )
            tmp = apply(matrix(as.numeric(tss[kk, mm]), nrow = 2), 2, mean)
            if(log2fc){
              test[ ,m] = c(tmp[c(1:2)] - tmp[3], 0, tmp[5:8] - tmp[4])
            }else{
              test[ ,m] = c(tmp[c(1:2)], mean(tmp[3:4]), tmp[5:8])
            }
          }
          
        }
        
        test = data.frame(test, sample = rownames(test))
        
        #test[which(test[,1]<(-2)),1] = -2
        #test[,1] = (test[, 1] - min(test[, 1]))/(max(test[,1]) - min(test[, 1])) * 3
        ylims = range(test[, c(1:length(features))])
        ylims[1] = min(c(0, ylims[1]))
        as_tibble(test) %>%  gather(features, signals,  1:4) %>% 
          mutate(features = factor(features, levels=c('atac', 'H3K4me3', 'H3K27me3', 'H3K4me1'))) %>%
          mutate(sample = factor(sample, levels = c('mHand','mLA', 'mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
          ggplot(aes(y=signals, x=sample, color = features, group = features)) + 
          geom_line(aes(linetype=features, color=features), size = 1.0) +
          geom_point(size = 4.5) +
          theme_classic() +
          geom_hline(yintercept=c(0), col="gray", size = 0.5) +
          #geom_vline(xintercept=c(3), col="blue", size = 1.2) +
          theme(axis.text.x = element_text(angle = 0, size = 14), 
                axis.text.y = element_text(angle = 0, size = 14), 
                axis.title =  element_text(size = 14),
                legend.text = element_text(size=10),
                legend.title = element_text(size = 10)
          ) +
          scale_color_manual(values=c('springgreen', 'blue', 'red','gold2')) +
          scale_linetype_manual(values=c("solid", "longdash", "solid", "longdash")) + 
          labs( x = '', y = 'feature signals') +
          ylim(ylims[1], ylims[2]) +
          ggtitle(geneList[n])
        
        ggsave(paste0(outDir, "/Regeneration_allFeatures_", geneList[n],  ".pdf"),  width = 8, height = 4) 
        
      }else{
        samples = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
        test = matrix(NA, nrow = length(samples), ncol = length(features))
        rownames(test) = samples
        colnames(test) = features
        
        for(m in 1:ncol(test))
        {
          if(m == 1) {
            mm =c(grep(paste0(colnames(test)[m], '_mUA|', colnames(test)[m], '_X'), colnames(tss)))
            if(log2fc){
              test[,m] = as.numeric(tss[kk, mm]) - as.numeric(tss[kk, mm[1]])
            }else{
              test[,m] = as.numeric(tss[kk, mm])
            }
            
          }
          if(m > 1){
            mm0 = grep(paste0(colnames(test)[m], '_mUA'), colnames(tss))
            mm0 = mm0[grep('Rep1$|Rep2$', colnames(tss)[mm0])]
            mm = c(mm0[grep('rRep', colnames(tss)[mm0])], 
                   grep(paste0(colnames(test)[m], '_BL'), colnames(tss))
            )
            tmp = apply(matrix(as.numeric(tss[kk, mm]), nrow = 2), 2, mean)
            test[ ,m] = tmp 
            if(log2fc){
              test[,m] = test[,m] - test[1,m]
            }
            
          }
          
        }
        
        test = data.frame(test, sample = rownames(test))
        
        #test[which(test[,1]<(-2)),1] = -2
        #test[,1] = (test[, 1] - min(test[, 1]))/(max(test[,1]) - min(test[, 1])) * 3
        ylims = range(test[, c(1:length(features))])
        ylims[1] = min(c(0, ylims[1]))
        as_tibble(test) %>%  gather(features, signals,  1:4) %>% 
          mutate(features = factor(features, levels=c('atac', 'H3K4me3', 'H3K27me3', 'H3K4me1'))) %>%
          mutate(sample = factor(sample, levels = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
          ggplot(aes(y=signals, x=sample, color = features, group = features)) + 
          geom_line(aes(linetype=features, color=features), size = 1.2) +
          geom_point(size = 4.0) +
          theme_classic() +
          geom_hline(yintercept=c(0), col="darkgray", size = 0.7) +
          #geom_vline(xintercept=c(3), col="blue", size = 1.2) +
          theme(axis.text.x = element_text(angle = 0, size = 14), 
                axis.text.y = element_text(angle = 0, size = 14), 
                axis.title =  element_text(size = 14),
                legend.text = element_text(size=12),
                legend.title = element_text(size = 14)
          ) +
          scale_color_manual(values=c('springgreen', 'blue', 'red','gold2')) +
          scale_linetype_manual(values=c("solid", "longdash", "solid", "longdash")) + 
          labs( x = '', y = 'feature signals') +
          ylim(ylims[1], ylims[2]) +
          ggtitle(geneList[n])
        
        ggsave(paste0(outDir, "/Regeneration_allFeatures_", geneList[n],  ".pdf"),  width = 7, height = 4) 
      }
      
    }else{
      cat(n, ' -- ', geneList[n], 'FOUND in tss',  length(kk), '\n')  
    }
  }
  
}

Analysis_TSS_positionalGenes_in_mature_regeneration = function(tss, ids)
{
  features = c('atac', 'H3K4me3', 'H3K27me3', 'H3K4me1', 'H3K27ac')
  
  kks = match(ids, tss$geneID)
  kks = kks[!is.na(kks)]
  
  samples = c('mHand', 'mLA', 'mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
  test = matrix(NA, nrow = length(kks), ncol = length(samples)*length(features))
  rownames(test) = tss$geneID[kks]
  test = data.frame(test, gene = tss$gene[kks], stringsAsFactors = FALSE)
  test$coords = tss$coords[kks]
  #colnames(test) = features
  
  for(n in 1:nrow(test))
  { 
    # n = 1
    cat(n, ' -- ', test$gene[n], '-- with tss :', tss$coords[n], '\n')
    kk = kks[n]
    
    for(m in 1:length(features))
    {
      if(m == 1) {
        mm =c(grep(paste0(features[m], '.M_mHand$'), colnames(tss)),  
              grep(paste0(features[m], '.M_mLA$'), colnames(tss)),
              grep(paste0(features[m], '.M_mUA$'), colnames(tss)),
              grep(paste0(features[m], '_mUA|', features[m], '_X'), colnames(tss)))
        test[n,c(((m-1)*length(samples) +1):(m*length(samples)))] = c(as.numeric(tss[kk, mm[1:2]]) - as.numeric(tss[kk, mm[3]]), 0, 
                     as.numeric(tss[kk, mm[5:8]])- as.numeric(tss[kk, mm[4]]))
      }
      if(m > 1 & m <5){
        mm0 = grep(paste0(features[m], '_mUA'), colnames(tss))
        mm0 = mm0[grep('Rep1$|Rep2$', colnames(tss)[mm0])]
        mm = c(grep(paste0(features[m], '_mHand'), colnames(tss)),
               grep(paste0(features[m], '_mLA'), colnames(tss)),
               mm0[grep('mRep', colnames(tss)[mm0])], 
               mm0[grep('rRep', colnames(tss)[mm0])], 
               grep(paste0(features[m], '_BL'), colnames(tss))
        )
        tmp = apply(matrix(as.numeric(tss[kk, mm]), nrow = 2), 2, mean)
        test[n,c(((m-1)*length(samples) +1):(m*length(samples)))] = c(tmp[c(1:2)] - tmp[3], 0, tmp[5:8] - tmp[4])
      }
      
      if(m >=5){
        mm0 = grep(paste0(features[m], '_mUA'), colnames(tss))
        mm0 = mm0[grep('Rep1$|Rep2$', colnames(tss)[mm0])]
        mm = c(grep(paste0(features[m], '_mHand'), colnames(tss)),
               grep(paste0(features[m], '_mLA'), colnames(tss)),
               mm0[grep('mRep', colnames(tss)[mm0])], 
               mm0[grep('rRep', colnames(tss)[mm0])], 
               grep(paste0(features[m], '_BL'), colnames(tss))
        )
        tmpx = as.numeric(tss[kk, mm])
        tmp = apply(matrix(tmpx[-3], nrow = 2), 2, mean)
        tmp = c(tmp[1], tmpx[3], tmp[2:length(tmp)])
        rm(tmpx)
        test[n,c(((m-1)*length(samples) +1):(m*length(samples)))] = c(tmp[c(1:2)] - tmp[3], 0, tmp[5:8] - tmp[4])
        
      }
      
    }
    
  }
  
  colnames(test)[c(1:(length(features)*length(samples)))] = paste0(rep(features, each = length(samples)), '_', samples)
  saveRDS(test, file = paste0(RdataDir, '/positional_gene_TSS_chromatinFeatures.rds'))
  
  ### add RNA-seq data
  test = readRDS(file = paste0(RdataDir, '/positional_gene_TSS_chromatinFeatures.rds'))
  samples = c('mHand', 'mLA', 'mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
  
  ## here import the microarray data of all genes, becasue Prod1, Tig1 are missing in the list of postional genes identified
  ## but still HOXA9, HOXA11 and HOXD9 are missing in the microarray probes
  ## test the microarray data, however the dynamic ranges are very different from the smartseq2 data in regeneration
  ## so use smart-seq2 mature samples as well
  rna = readRDS( file = paste0("../results/RNAseq_data_used/Rdata/matureSamples_cpm_DEgenes_8selectedSamples.batch4_v47.hox.patch.rds"))
  ids = rna$geneID
  ggs = rna$gene
  
  rna = data.frame(rna$log2FoldChange_mHand.vs.mUA, rna$log2FoldChange_mLA.vs.mUA, rep(0, nrow(rna)))
  mm = match(rownames(test), ids)
  missed = which(is.na(mm))
  #mm_missed = match(rownames(test)[missed], ggs)
  #mm[missed] = mm_missed 
  
  rna = rna[mm, ]
  rownames(rna) = rownames(test)
  colnames(rna) = samples[1:3]
  colnames(rna) = paste0('rna_', colnames(rna))
  #rna = rna[, c(3:1)]
  
  ## add regeneration data
  res = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/', 
          'smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked_RNAseq_data_used_20220408.rds'))
  cpm = res[, grep('log2FoldChange_d', colnames(res))]
  #cpm = cal_sample_means(cpm, conds = c("Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal",  "BL_UA_13days_distal"))
  colnames(cpm) = samples[4:7]
  #cpm = cpm - cpm[,1]
  ids = get_geneID(rownames(cpm))
  mm = match(rownames(rna), ids)
  
  cpm = cpm[mm, ]
  colnames(cpm) = paste0('rna_', colnames(cpm))
  
  rna = data.frame(rna, cpm, stringsAsFactors = FALSE)
  
  test = data.frame(test, rna, stringsAsFactors = FALSE)
  
  saveRDS(test, file = paste0(RdataDir, '/positional_gene_TSS_chromatinFeatures_smartseq2_mature.reg.rds'))
  
  ### SVD analysis
  Test.SVD = FALSE
  if(Test.SVD){
    for(n in 1:nrow(test))
    {
      # n = 7
      cat(n, '--', test$gene[n], '\n')
      xx = t(matrix(as.numeric(test[n, c(1:35)]), nrow = length(samples)))
      rownames(xx) = features
      colnames(xx) = samples
      #xx = t(apply(xx, 1, cal_z_score))
      ss = svd(xx)
      #L = length(ss$d)
      pl = ss$d^2/sum(ss$d^2)
      
      #### first 2 components, first change the first two conlums to positive
      u1 = ss$u[,1]
      #u2 = ss$u[,2]
      v1 = ss$v[,1]
      #v2 = ss$v[,2]
      #nb.positive = length(which(u1>=0))
      #nb.negative = length(which(u1<0))
      #if(nb.positive<nb.negative){u1 = -u1; v1 = -v1;}
      
      xx = rbind(xx, v1)
      
      data.frame(xx, fts = rownames(xx)) %>%
        as_tibble() %>% gather(ss, signals,  1:7) %>%
        mutate(ss = factor(ss, levels = c('mHand', 'mLA', 'mUA', 'X5dpa', 'X9dpa', 'X13dpa.p', 'X13dpa.d'))) %>%
        ggplot(aes(y=signals, x=ss, color = fts, group = fts)) + 
        geom_line(aes(linetype=fts, color=fts), size = 0.7) +
        theme_classic() +
        scale_color_manual(values=c('springgreen', 'gold2', 'black', 'blue', 'orange', 'red')) +
        scale_linetype_manual(values=c("longdash", "longdash", "longdash", "longdash", "longdash", 'solid'))
    }
  }
  
  return(test)
  
}
