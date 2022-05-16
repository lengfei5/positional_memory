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