##########################################################################
##########################################################################
# Project: Akane's positional memory project
# Script purpose: regeneration time points analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Feb 22 10:46:43 2022
##########################################################################
##########################################################################
rm(list=ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('functions_chipSeq.R')
source('Functions_atac.R')

version.analysis = 'Rxxxx_R10723_R11637_R12810_atac'
#peakDir = "Peaks/macs2_broad"

resDir = paste0("../results/", version.analysis)

RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/'
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0(figureDir, 'tables4plots/')

saveTables = FALSE

require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require(pheatmap)
library(tictoc)


########################################################
########################################################
# Section I : normalization and batch correction
# TMM and combat were selected for normalization and batch correction
########################################################
########################################################
source('Functions_atac.R')
library(edgeR)
require("sva")
require(limma)

load(file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata', 
                   '/ATACseq_selected.63k.peaks_cutoff.40.at.least.2sample.Rdata'))

Save.peak.consensus = FALSE
if(Save.peak.consensus){
  pp = data.frame(t(sapply(rownames(dds), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp$strand = '*'

  bed = data.frame(pp[, c(1:3)], rownames(pp), 0, pp[, 4], stringsAsFactors = FALSE)
  write.table(bed, file = paste0(resDir, '/ATAC_peakset_64K_', version.analysis, '.bed'), col.names = FALSE, row.names = FALSE,
              sep = '\t', quote = FALSE)
  
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")
  saveRDS(pp, file = paste0(RdataDir, '/ATACseq_peak_consensus_filtered_64k.rds'))
  
}


table(design$condition, design$batch)

# regeneration time points and embryo stages
Batch.Correct.regeneration.embryoStage = FALSE
if(Batch.Correct.regeneration.embryoStage){
  
  design$batch[which(design$condition == 'BL_UA_9days')] = '2021'
  
  table(design$condition, design$batch)
  
  sels = which(design$batch == '2021'| (design$batch == '2020' & design$condition != 'Mature_LA' & design$condition != 'Mature_Hand' &
                                          design$condition != 'Mature_UA'))
  
  # sels = sels[which(design$SampleID[sels] != '89542' & design$SampleID[sels] != '89543')]
  
  design.sels = design[sels, ]
  design.sels$conds = droplevels(design.sels$conds)
  
  table(design.sels$conds, design.sels$batch)
  
  #design.sels$batch[grep('749', design.sels$SampleID)] = '2019.1'
  #design.sels$batch[grep('1026', design.sels$SampleID)] = '2019.2'
  
  #design.sels$batch[which(design.sels$batch == '2021S')] = '2021'
  #design.sels$batch = droplevels(design.sels$batch)
  table(design.sels$conds, design.sels$batch)
  
  ddx = dds[, sels]
  ddx$conds = droplevels(ddx$conds)
  ss = rowSums(counts(ddx))
  # remove low count genes, otherwise combat returns error 
  # 'Error in while (change > conv) { : missing value where TRUE/FALSE needed'
  #ddx = ddx[which(ss>5), ] 
  #ddx = estimateSizeFactors(ddx)
  #vsd <- varianceStabilizingTransformation(ddx, blind = TRUE)
  #tmm = assay(vsd)
  
  d <- DGEList(counts=counts(ddx), group=design.sels$conds)
  tmm <- calcNormFactors(d, method='TMM')
  tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  
  #tmm.vars = apply(as.matrix(tmm), 1, var) # row with var = 0 pose problem for ComBat
  #tmm = tmm[which(tmm.vars>0 & !is.na(tmm.vars)), ]
  
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
  ggsave(paste0(resDir, "/regeneration_embryo_Samples_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)
  
  make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
  ggsave(paste0(resDir, "/matureSamples_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)
  
  fpm = fpm.bc
  
  rm(fpm.bc)
  
  saveRDS(fpm, file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration_dev_2Batches.R10723_R7977_', version.analysis, '.rds'))
  saveRDS(design.sels, file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration_dev_2Batches.R10723_R7977',
                                     version.analysis, '.rds'))
  
  
}

########################################################
########################################################
# Section II : regeneration peaks
# or temporal-peaks test
########################################################
########################################################
##########################################
grouping.temporal.peaks = FALSE
if(grouping.temporal.peaks){
  require(ChIPpeakAnno)
  require(ChIPseeker)
  
  fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration_dev_2Batches.R10723_R7977_', version.analysis, '.rds'))
  design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration_dev_2Batches.R10723_R7977',
                                 version.analysis, '.rds'))
  
  # prepare the background distribution
  fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
  fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
  # rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  hist(fpm.bg, breaks = 100, main = 'background distribution')
  abline(v = c(1, 2.5, 3), col = 'red', lwd = 2.0)
  quantile(fpm.bg, c(0.95, 0.99))
  
  ##########################################
  ## make Granges and annotate peaks
  ##########################################
  Make.Granges.and.peakAnnotation = TRUE
  if(Make.Granges.and.peakAnnotation){
    
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
    
    ## annotation from ucsc browser ambMex60DD_genes_putative
    ## to update in the future 
    amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
    pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    #plotAnnoBar(pp.annots)
    
    pp.annots = as.data.frame(pp.annots)
    rownames(pp.annots) = rownames(fpm)
    
    promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
    
  }
  
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
  
  Run.temporal.peak.test = FALSE
  if(Run.temporal.peak.test)
  {
    source('Functions_atac.R')
    cpm = fpm[, sample.sels]
    
    # stringent here, with signal > 2.5 in >=3 samples 
    #nb.above.threshold = apply(as.matrix(cpm), 1, function(x) length(which(x> 2.5)))
    #hist(nb.above.threshold, breaks = c(-1:ncol(cpm)))
    maxs = apply(sample.means, 1, max)
    length(which(maxs>2))
    length(which(maxs>2.5))
    length(which(maxs>3))
    
    # no further filtering anymore
    # peak.sels = which(maxs > 2)
    # cat(length(peak.sels), 'peaks after filtering for mature samples\n')
    # cpm = cpm[peak.sels, ]
    # pp.sel = pp[peak.sels]
    
    # examples to test
    test.examples = c('HAND2', 'FGF8', 'KLF4', 'Gli3', 'Grem1')
    ii.test = which(overlapsAny(pp.sel, promoters[which(!is.na(match(promoters$geneSymbol, test.examples)))]))
    
    source('Functions_atac.R')
    tic()
    ## define the dynamic enhancers with mature UA and BL.UA and check them if embryo samples
    #sels = grep('Embryo', cc, invert = TRUE) 
    res = temporal.peaks.test(cpm, c = cc)
    
    toc()
    
    
    saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks_mUA_regeneration_dev_2Batches.R10723_R7977_v8_noPeakAnnot.rds'))
    
    xx = data.frame(cpm, res, pp.annots[match(rownames(cpm), rownames(pp.annots)), ],  stringsAsFactors = FALSE)
    
    res = xx
    
    #res = data.frame(cpm, res, stringsAsFactors = FALSE)
    saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks__mUA_regeneration_dev_2Batches.R10723_R7977_peakAnnot_v8.rds'))
    
    
  }
  
  ##########################################
  # reload the test result and select dynamic peaks
  ##########################################
  library(qvalue)
  res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks__mUA_regeneration_dev_2Batches.R10723_R7977_peakAnnot_v8.rds'))
  
  # select only the dynamic peak results without annotation part
  res = res[, c(1:50)]
  
  res = res[order(-res$log2FC), ]
  qv = qvalue(res$pval.lrt)
  res$fdr.lrt = qv$qvalues
  
  # select the temporal dynamic peaks
  fdr.cutoff = 0.01; logfc.cutoff = 1
  
  # length(which(res$fdr.lrt < fdr.cutoff))
  # length(which(res$padj_LRT<fdr.cutoff & res$log2fc>1))
  # length(which(res$padj_LRT<fdr.cutoff & res$log2fc>2))
  # length(which(res$padj_LRT<fdr.cutoff & res$log2fc>1.5))
  
  select = which(  (res$adj.P.Val_5dpa.vs.mUA < fdr.cutoff & abs(res$logFC_5dpa.vs.mUA) > logfc.cutoff)| 
                   (res$adj.P.Val_9dpa.vs.mUA < fdr.cutoff & abs(res$logFC_9dpa.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_13dpap.vs.mUA < fdr.cutoff & abs(res$logFC_13dpap.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_13dpad.vs.mUA < fdr.cutoff & abs(res$logFC_13dpad.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_s40.vs.mUA < fdr.cutoff & abs(res$logFC_s40.vs.mUA) > logfc.cutoff ) | 
                   (res$adj.P.Val_s44p.vs.mUA < fdr.cutoff & abs(res$logFC_s44p.vs.mUA) > logfc.cutoff ) |
                   (res$adj.P.Val_s44d.vs.mUA < fdr.cutoff & abs(res$logFC_s44d.vs.mUA) > logfc.cutoff ))
  
  cat(length(select), 'DE peaks found !\n')
  
  res = res[select, ]
  res = res[order(-res$log2FC), ]
  
  ##########################################
  # heatmap for all regeneration dynamic peaks
  ##########################################
  res$mins =  apply(sample.means[match(rownames(res), rownames(sample.means)), ], 1, min)
  # xx = xx[which(xx$min < 4.5), ]
  
  res = res[which(res$mins < 3), ]
  cat(nrow(res), 'DE peaks found with at least condition below background !\n')
  
  source('Functions_atac.R')
  conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", 'Embryo_Stage44_distal', 
            "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", 'BL_UA_13days_distal'
            )
  
  keep = as.matrix(res[, c(1:20)])
  sample.sels = c(); cc = c()
  for(n in 1:length(conds)) {
    kk = grep(conds[n], colnames(keep))
    if(conds[n] == 'BL_UA_5days') kk = kk[grep('136', colnames(keep)[kk])]
    if(conds[n] == 'Embryo_Stage40') kk = kk[grep('137', colnames(keep)[kk])]
    #if(length(unique(design$batch[kk])) > 1) kk = kk[which(design$batch[kk] == '2021')]
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
  }
  
  keep = keep[, sample.sels]
  
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  library(dendextend)
  
  yy <- t(apply(keep, 1, cal_z_score))
  
  df <- data.frame(cc)
  rownames(df) = colnames(keep)
  colnames(df) = 'sample'
  
  sample_colors = c('magenta', 'darkblue', 'springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2',
                    'red')[c(1:length(conds))]
  
  names(sample_colors) = conds
  
  nb_clusters = 8
  my_hclust_gene <- hclust(dist(yy), method = "complete")
  
  my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
  my_gene_col <- data.frame(cluster =  paste0('cluster_', my_gene_col))
  rownames(my_gene_col) = rownames(yy)
  
  res$clusters = my_gene_col
  
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
            "#33a02c", "#fb9a99", "#e31a1c",
            "#fdbf6f", "#ff7f00", "#cab2d6",
            "#6a3d9a", "#ffff99", "#b15928")
  cluster_col = col3[1:nb_clusters]
  names(cluster_col) = paste0('cluster_', c(1:nb_clusters))
  annot_colors = list(
    condition = sample_colors,
    cluster = cluster_col)
    
  
  ii.gaps = c(6, 8)
  col = colorRampPalette(c("navy", "white", "red3"))(16)
  
  pheatmap(yy, 
           annotation_row = my_gene_col, 
           annotation_col = df, show_rownames = FALSE, scale = 'none', 
           color = col, 
           show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = FALSE,  
           clustering_method = 'complete', cutree_rows = nb_clusters, 
           annotation_colors = annot_colors, 
           #clustering_callback = callback,
           gaps_col = ii.gaps, 
           filename = paste0(resDir, '/heatmap_regenerationPeaks_M0prob0.05_log2FC.1_scaled.pdf'), 
           width = 8, height = 16)
  
  pheatmap(keep, 
           # annotation_row = my_gene_col, 
           annotation_col = df, show_rownames = FALSE, scale = 'none', 
           color = col, 
           show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = FALSE,  
           #clustering_method = 'complete', cutree_rows = nb_clusters, 
           #annotation_colors = annot_colors, 
           #clustering_callback = callback,
           gaps_col = ii.gaps, 
           filename = paste0(resDir, '/heatmap_regenerationPeaks_M0prob0.05_log2FC.1.pdf'), 
           width = 8, height = 12)
  
  ##########################################
  # highligh potential regeneration peaks, not found in mUA and embryo stages only in regeneration process 
  ##########################################
  means.sel = sample.means[match(rownames(keep), rownames(sample.means)), ]
  
  dev.mature.maxs = apply(means.sel[, grep('Embryo|Mature_UA', colnames(means.sel))] , 1, max)
  bl.maxs = apply(means.sel[ , grep('BL_UA', colnames(means.sel))], 1, max)
  
  kk = which(dev.mature.maxs < 2.5 & bl.maxs > 3)
  
  yy0 = yy[kk, ]
  
  pheatmap(yy0, 
           #annotation_row = my_gene_col, 
           annotation_col = df, show_rownames = FALSE, scale = 'none', 
           color = col, 
           show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = FALSE,  
           #clustering_method = 'complete', cutree_rows = nb_clusters, 
           #annotation_colors = sample_colors,
           #clustering_callback = callback,
           gaps_col = ii.gaps, 
           filename = paste0(resDir, '/heatmap_regenerationPeaks_edgeRtest_fdr0.05_log2FC.1_regeneartion.specific_v2.pdf'), 
           width = 8, height = 6)
    
  if(saveTable){
    write.csv(data.frame(yy0, stringsAsFactors = FALSE), 
              file = paste0(resDir, '/regeneration_peaks_regeneration.specific.csv'), 
              quote = FALSE, row.names = TRUE)
    
  }
  
}


########################################################
########################################################
# Section III : motif analysis
# 
########################################################
########################################################
##########################################
# first motif activity analysis for temporally dynamic peaks 
##########################################
source('MARA_functions.R')
library(qvalue)
res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks__mUA_regeneration_dev_2Batches.R10723_R7977_peakAnnot_v8.rds'))

Select.dynamic.peaks.for.MARA = TRUE
if(Select.dynamic.peaks.for.MARA){
  res = res[order(-res$log2FC), ]
  qv = qvalue(res$pval.lrt)
  res$fdr.lrt = qv$qvalues
  
  # select the temporal dynamic peaks
  fdr.cutoff = 0.01; logfc.cutoff = 1.
  
  # length(which(res$fdr.lrt < fdr.cutoff))
  # length(which(res$padj_LRT<fdr.cutoff & res$log2fc>1))
  # length(which(res$padj_LRT<fdr.cutoff & res$log2fc>2))
  # length(which(res$padj_LRT<fdr.cutoff & res$log2fc>1.5))
  select = which(  (res$adj.P.Val_5dpa.vs.mUA < fdr.cutoff & abs(res$logFC_5dpa.vs.mUA) > logfc.cutoff)| 
                     (res$adj.P.Val_9dpa.vs.mUA < fdr.cutoff & abs(res$logFC_9dpa.vs.mUA) > logfc.cutoff) |
                     (res$adj.P.Val_13dpap.vs.mUA < fdr.cutoff & abs(res$logFC_13dpap.vs.mUA) > logfc.cutoff) |
                     (res$adj.P.Val_13dpad.vs.mUA < fdr.cutoff & abs(res$logFC_13dpad.vs.mUA) > logfc.cutoff) 
                   #(res$adj.P.Val_s40.vs.mUA < fdr.cutoff & abs(res$logFC_s40.vs.mUA) < logfc.cutoff ) | 
                   #(res$adj.P.Val_s44p.vs.mUA < fdr.cutoff & abs(res$logFC_s44p.vs.mUA) < logfc.cutoff ) |
                   #(res$adj.P.Val_s44d.vs.mUA < fdr.cutoff & abs(res$logFC_s44d.vs.mUA) < logfc.cutoff )
  )
  cat(length(select), 'DE peaks found !\n')
  
  res = res[select, ]
  res = res[order(-res$log2FC), ]
  keep = as.matrix(res[, c(1:20)])
  
  conds = c("Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", 'BL_UA_13days_distal')
  
  sample.sels = c(); cc = c()
  for(n in 1:length(conds)) {
    kk = grep(conds[n], colnames(keep))
    if(conds[n] == 'BL_UA_5days') kk = kk[grep('136', colnames(keep)[kk])]
    if(conds[n] == 'Embryo_Stage40') kk = kk[grep('137', colnames(keep)[kk])]
    #if(length(unique(design$batch[kk])) > 1) kk = kk[which(design$batch[kk] == '2021')]
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
  }
  
  keep = keep[, sample.sels]
  
  
}

xx = run.MARA.atac.temporal(keep, cc)



