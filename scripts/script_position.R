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
source('Functions_histM.R')
source('Functions_Integration.matureReg.R')

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
tableDir = paste0('/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/SupTables/')

saveTables = FALSE

require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require(pheatmap)
library(tictoc)
library(tidyr)
library(dplyr)
require(ggplot2)
library("gridExtra")
library("cowplot")
require(ggpubr)

########################################################
########################################################
# Section I: batch correction 
# TMM and combat were selected for normalization and batch correction
########################################################
########################################################
source('Functions_atac.R')
library(edgeR)
require("sva")
require(limma)

load(file = paste0(RdataDir, '/ATACseq_selected.55k.peaks_cutoff.50.at.least.1sample.Rdata'))

Split.Mature.Regeneration.samples = TRUE

if(Split.Mature.Regeneration.samples){
 
  ##########################################
  # batch correct samples separately for mature samples and regeneration samples  
  ##########################################
  table(design$condition, design$batch)
  
  # start with mature samples
  Batch.Correct.matureSamples = FALSE 
  if(Batch.Correct.matureSamples){
    
    sels = grep('Mature|HEAD', design$conds)
    
    #sels = setdiff(sels, which(design$SampleID == '74938'| design$SampleID == '74939'|design$SampleID == '74940'))
    design.sels = design[sels, ]
    
    design.sels$conds = droplevels(design.sels$conds)
    
    design.sels$batch[grep('749', design.sels$SampleID)] = '2019'
    
    #design.sels$batch = droplevels(design.sels$batch)
    table(design.sels$conds, design.sels$batch)
    
    ddx = dds[, sels]
    ddx$conds = droplevels(ddx$conds)
    ss = rowSums(counts(ddx))
    
    save(ddx, design.sels, file = paste0(RdataDir, '/atac_matureSamples_beforeBatchCorrection.Rdata'))
    # remove low count genes, otherwise combat returns error 
    # 'Error in while (change > conv) { : missing value where TRUE/FALSE needed'
    #ddx = ddx[which(ss>5), ] 
    #ddx = estimateSizeFactors(ddx)
    
    #vsd <- varianceStabilizingTransformation(ddx, blind = TRUE)
    
    #tmm = assay(vsd)
    
    #tmm = log2(fpm(ddx) + )
    #tmm.vars = apply(as.matrix(tmm), 1, var) # row with var = 0 pose problem for ComBat
    #tmm = tmm[which(tmm.vars>0 & !is.na(tmm.vars)), ]
    
    d <- DGEList(counts=counts(ddx), group=design.sels$conds)
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
    make.pca.plots(tmm, ntop = 1000, conds.plot = 'Mature')
    ggsave(paste0(resDir, "/matureSamples_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    
    make.pca.plots(fpm.bc, ntop = 1000, conds.plot = 'Mature')
    ggsave(paste0(resDir, "/matureSamples_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    
    fpm = fpm.bc
    
    rm(fpm.bc)
    
    saveRDS(fpm, file = paste0(RdataDir, '/fpm.bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.2022.rds'))
    saveRDS(design.sels, file = paste0(RdataDir, '/design_sels_bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.2022.rds'))
      
  }
      
}
  

########################################################
########################################################
# Section II : positional peaks
# ##########################################
# Position-dependent test
# mainly use the mature samples, mUA, mLA and mHand
# 
##########################################
########################################################
########################################################
require(ChIPpeakAnno)
require(ChIPseeker)

##########################################
# import batch corrected gene expression and design 
##########################################
fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.2022.rds'))
design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.2022.rds'))

# prepare the background distribution
fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
rownames(fpm) = gsub('_', '-', rownames(fpm))

hist(fpm.bg, breaks = 100, main = 'background distribution')
abline(v = 3, col = 'red', lwd = 2.0)
quantile(fpm.bg, c(0.95, 0.99))

##########################################
## make Granges and annotate peaks
##########################################
Make.Granges.and.peakAnnotation = FALSE
if(Make.Granges.and.peakAnnotation){
  pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  # saveRDS(pp, file = paste0(RdataDir, '/ATACseq_peak_consensus_filtered_55k.rds'))
  # export(object = pp,  con = paste0(resDir, "/atacseq_peaks_filtered_55k.bed"), format = 'bed')
  
  # annotation from ucsc browser ambMex60DD_genes_putative
  amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
  pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    
  pp.annots = as.data.frame(pp.annots)
  rownames(pp.annots) = rownames(yy)
  
  
  promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
  
}

##########################################
# run spatial test for mature samples
##########################################
Run.test.spatial.peaks = FALSE
if(Run.test.spatial.peaks){
  conds = c("Mature_UA", "Mature_LA", "Mature_Hand")
  
  sample.sels = c();  cc = c()
  for(n in 1:length(conds)) {
    #kk = which(design$conds == conds[n] & design$SampleID != '136159')
    kk = which(design$conds == conds[n]) 
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
    
  }
  
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
  
  further.peak.signal.filtering = FALSE
  if(further.peak.signal.filtering){
    signal.threshold = 3.0
    # stringent here, with signal > 2.5 in >=3 samples 
    nb.above.threshold = apply(as.matrix(cpm), 1, function(x) length(which(x>signal.threshold )))
    hist(nb.above.threshold, breaks = c(-1:ncol(cpm)))
    peak.sels = which(nb.above.threshold>=2)
    cat(length(peak.sels), 'peaks after filtering for mature samples\n')
    
    cpm = cpm[peak.sels, ]
    
  }
  
  tic() 
  res = spatial.peaks.test(cpm = cpm, c = cc, test.Dev.Reg = FALSE)
  #res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
  toc()
  
  xx = data.frame(res, pp.annots[match(rownames(res), rownames(pp.annots)), ], stringsAsFactors = FALSE)
  
  res = xx
  saveRDS(res, file = paste0(RdataDir, '/res_position_dependant_test_', version.analysis, '_v12.rds'))
  
  rm(xx)
  
}

##########################################
# select all positional-dependent loci with below threshold
##########################################
res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_', version.analysis, '_v12.rds'))
fdr.cutoff = 0.05; logfc.cutoff = 1

Filter.select.positional.atacPeaks = FALSE
if(Filter.select.positional.atacPeaks){
  # select the positional peaks with pairwise comparisions 
  # limma logFC is in log2 scale
  
  jj = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
               (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
               (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)
  )
  cat(length(jj), '\n')
  
  jj1 = which(res$prob.M0< fdr.cutoff & res$log2FC>logfc.cutoff)
  cat(length(jj1), '\n')
  jj2 = which(res$pval.lrt < fdr.cutoff & res$log2FC > logfc.cutoff)
  cat(length(jj2), '\n')
  
  
  xx = res[c(jj), ]
  #xx = xx[order(-xx$log2FC.mature), ]
  xx[grep('HOXA13|SHOX', xx$transcriptId), ]
  # fpm[which(rownames(fpm) == 'chr2p:873464923-873465440'), ]
  
  xx[grep('HOXA13|SHOX|MEIS', xx$transcriptId), ]
  
  cat(nrow(xx), ' peaks left\n')
  # sort positional peaks with logFC
  #xx = xx[order(-xx$logFC.mean), ]
  xx = xx[order(-xx$log2FC), ]
  
  ########
  ## asscociate the signifiant postional peaks with expression matrix
  ########
  conds = c("Mature_UA", "Mature_LA", "Mature_Hand", 'HEAD')
  
  sample.sels = c();  cc = c()
  for(n in 1:length(conds)) {
    kk = which(design$conds == conds[n]) 
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
  }
  
  keep = fpm[(match(rownames(xx), rownames(fpm))), sample.sels]
  keep = as.matrix(keep)
  
  ### filter the peaks from head control sample
  ## either filter soly based on the head signals or based on the logFC between max(mature samples/control) or both
  Filtering.peaks.in.Head.samples = TRUE
  if(Filtering.peaks.in.Head.samples){
    maxs = apply(keep[, grep('Mature_', colnames(keep))], 1, function(x) return(max(c(mean(x[1:5]), mean(x[6:9]), mean(x[10:12])))))
    mins = apply(keep[, grep('Mature_', colnames(keep))], 1, function(x) return(min(c(mean(x[1:5]), mean(x[6:9]), mean(x[10:12])))))
    
    ctl.mean = apply(keep[, grep('HEAD', colnames(keep))], 1, mean)
    
    rr = maxs - ctl.mean
    #p.ctl = pp[match(names(ctl.sels), names(pp))]
    #non.overlap = !overlapsAny(p0, p.ctl)
    plot(maxs, ctl.mean, cex = 0.2);
    abline(0, 1, lwd = 2.0, col = 'red')
    abline(v = 3, lwd = 2.0, col = 'red')
    abline(h = 3, lwd = 2.0, col = 'red')
    
    #xx = xx[non.overlap, ]
    plot(rr, ctl.mean, cex = 0.2);
    abline(h = 3, col = 'red', lwd = 2.0)
    abline(v = c(0.5, 1), col = 'red', lwd = 2.0)
    
    jj = which(maxs > 3 & mins <3)
    
    sels = which(rr>1 & ctl.mean<3 & maxs > 3 & mins <3)
    cat(length(sels), 'peaks selected \n')
    
    xx = xx[sels, ]
    keep = keep[sels, ]
    
  }
  
  dim(xx)
  
  save(xx, keep, file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl', version.analysis, '.Rdata'))
  
}

##########################################
# Positional atac-seq peaks heatmap 
# postional peak clustering using three replicates of each segments
#  one of main overview of  
##########################################
library(dendextend)
library(ggplot2)
source('Functions_histM.R')

### load the previous result: xx (test result and peak annotation) and keep (log2 data of mature samples)
load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl', version.analysis, '.Rdata'))
conds = c("Mature_UA", "Mature_LA", "Mature_Hand")

# saveRDS(keep, file = paste0(RdataDir, '/positional_peakSignals.rds'))
pp = data.frame(t(sapply(rownames(xx), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
pp$strand = '*'
pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")

# make.pca.plots(keep, ntop = 1246, conds.plot = 'Mature')
rep.sels = grep('HEAD|102657|102655|74938', colnames(keep), invert = TRUE)
yy = keep[, rep.sels]
yy <- t(apply(yy, 1, cal_z_score))
#saveRDS(yy, file = paste0(RdataDir, '/positional_atacPeaks_data_3reps_forHeatmap.rds'))

nb_clusters = 6
my_hclust_gene <- hclust(dist(yy), method = "complete")

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
xx$clusters = my_gene_col

my_gene_col <- data.frame(cluster =  paste0('cluster_', my_gene_col))
rownames(my_gene_col) = rownames(yy)

df <- data.frame(rep(conds[1:3], each = 3))
rownames(df) = colnames(yy)
colnames(df) = 'segments'

col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
          "#33a02c", "#fb9a99", "#e31a1c",
          "#fdbf6f", "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99", "#b15928")

sample_colors = c('springgreen4', 'steelblue2', 'gold2')
names(sample_colors) = c('Mature_UA', 'Mature_LA', 'Mature_Hand')
cluster_col = col3[1:nb_clusters]
names(cluster_col) = paste0('cluster_', c(1:nb_clusters))
annot_colors = list(
  segments = sample_colors,
  cluster = cluster_col)

gaps.col = c(3, 6)

col<- colorRampPalette(c("steelblue", "white", "darkred"))(8)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8)

plt = pheatmap(yy, annotation_row = my_gene_col, 
               annotation_col = df, show_rownames = FALSE, scale = 'none', 
               color = col, 
               show_colnames = FALSE,
               cluster_rows = TRUE, cluster_cols = FALSE,  
               clustering_method = 'complete', cutree_rows = nb_clusters, 
               annotation_colors = annot_colors, 
               gaps_col = gaps.col) 

row_order = data.frame(plt$tree_row$order, plt$tree_row$labels, stringsAsFactors = FALSE)
colnames(row_order) = c('index', 'name')
row_order$name = row_order$name[row_order$index]
row_order$cluster = my_gene_col$cluster[match(row_order$name, rownames(my_gene_col))]

callback = function(hc, mat = row_order){
  #sv = svd(t(mat))$v[,1]
  #clusters_order = 
  o1 = c()
  for(co in paste0('cluster_', c(1, 6, 5, 4, 3, 2))) {
    o1 = c(o1, row_order$index[which(row_order$cluster == co)])
    #row_order = row_order[o1, ]
  }
  dend = reorder(as.dendrogram(hc), wts = o1)
  as.hclust(dend)
  
}

col<- colorRampPalette(c("blue4", "white", "darkred"))(8)
col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdGy")))(8)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="PuOr")))(16)
col = palette(hcl.colors(8, "Viridis"))

col = colorRampPalette(c("navy", "white", "red3"))(8)

plt = pheatmap(yy, annotation_row = my_gene_col, 
               annotation_col = df, show_rownames = FALSE, scale = 'none', 
               color = col, 
               show_colnames = FALSE,
               cluster_rows = TRUE, cluster_cols = FALSE,  
               clustering_method = 'complete', cutree_rows = nb_clusters, 
               annotation_colors = annot_colors, 
               clustering_callback = callback,
               gaps_col = gaps.col, 
               treeheight_row = 20,
               annotation_legend = FALSE,
               filename = paste0(figureDir, '/positional_atacPeaks_fdr0.05_log2FC.1_rmHeadPeaks_v2.pdf'), 
               width = 4, height = 12)

saveRDS(plt, file = paste0(RdataDir, '/postional_atacPeaks_heatmap_orderSaved.rds'))

saveRDS(xx, file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))

write.csv(xx, file = paste0(figureDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters', 
                            nb_clusters, '.csv'), quote = FALSE, row.names = TRUE)

## save group bed of each cluster for deeptools heatmap
Save.peak.groups.for.Deeptools = FALSE
if(Save.peak.groups.for.Deeptools){
  res = readRDS(file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))
  mcs = unique(res$clusters)
  
  peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                        'positional_memory/Data//atacseq_using/makeHeatmaps/peak_groups')
  
  for(n in 1:length(mcs)){
    jj = which(res$clusters == mcs[n])
    pp_atac = data.frame(t(sapply(rownames(yy)[jj], function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    pp_atac$name = rownames(yy)[jj]
    pp_atac$score = 0
    pp_atac$strand = '*'
    
    write.table(pp_atac, file = paste0(peakGroupDir, '/peak_group_', mcs[n], '.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
}

##########################################
# Add histone marker with together as atac-seq peaks 
# this part is firstly done in histMarker_CT_analysis.R
# overview of segment-specific chromatin landscape
# histone marker with the same cluster order as atac-seq peaks 
# 
##########################################
Find.histonePeaks.overlapping.nonOverlapping.with.dynamicATACpeaks = FALSE
if(Find.histonePeaks.overlapping.nonOverlapping.with.dynamicATACpeaks)
{
  # load keep object in which all histM peaks overlapping with atac-seq peaks and segement-specific test 
  load(file = paste0('../results/CT_merged_20220328/Rdata', '/combined_4histMarkers_overlapped55kATACseq_DE_fdr0.05.Rdata'))
  design = readRDS(file = paste0('../results/CT_merged_20220328/Rdata', '/histM_CT_design_info.rds'))
  
  yy = keep[, c(53:60, 27:34, 1:8, 79:85)]
  sampleID = sapply(colnames(yy), function(x) unlist(strsplit(as.character(x), '_'))[3])
  mm = match(sampleID, design$sampleID)
  colnames(yy) = paste0(design$condition[mm], '_', design$Batch[mm], '_', design$sampleID[mm])
  #yy = yy[, grep('mRep1|mRep2', colnames(yy))]
  
  #peaks = read.csv(paste0(figureDir, 'positional_peaks_clusters_6/', 
  #                        'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters6.csv'))
  peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rdata/', 
                                'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))
  
  pp_atac = data.frame(t(sapply(rownames(peaks), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp_atac$strand = '*'
  pp_atac = makeGRangesFromDataFrame(pp_atac, seqnames.field=c("X1"),
                                     start.field="X2", end.field="X3", strand.field="strand")
  
  pp_histM = data.frame(t(sapply(rownames(yy), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp_histM$strand = '*'
  pp_histM = makeGRangesFromDataFrame(pp_histM, seqnames.field=c("X1"),
                                      start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(pp_atac, pp_histM)
  
  ## peaks mapped to postional atac-peaks 
  yy1 = yy[mapping@to, ]
  
  ## peaks not mapped to positional atac-peaks but overlapped wiht atac peaks
  yy0 = yy[-unique(mapping@to), ]
  
  ss = apply(DE.locus[, c(1:3)], 1, sum)
  mm = match(names(ss[which(ss>0)]), rownames(yy0))
  length(which(!is.na(mm)))
  
  yy0 = yy0[mm[!is.na(mm)], ] ## histM overlapped with stable atacPeak but significantly changed
  
  #plot.pair.comparison.plot(yy1[, grep('H3K4me1_mUA_', colnames(yy1))], linear.scale = FALSE)
  #plot.pair.comparison.plot(yy1[, grep('H3K4me3_mUA_', colnames(yy1))], linear.scale = FALSE)
  #plot.pair.comparison.plot(yy1[, grep('H3K27me3_mUA_', colnames(yy1))], linear.scale = FALSE)
  saveRDS(yy1, file = paste0(RdataDir, '/peak_signals_atac_4histM_positionalPeaks.rds'))
  saveRDS(yy0, file = paste0(RdataDir, '/peak_signals_atac_4histM_notOverlapped.positionalPeaks.rds'))
  
}

##########################################
# plot segment-specific histone marks overlapped with segement-specific atac
# and segment-specific histone marks overlapped with stable atac
# subclustering was performed 
##########################################
Plot_histMarkers_for_positionalATAC = FALSE
if(Plot_histMarkers_for_positionalATAC){
  
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  require(pheatmap)
  require(RColorBrewer)
  library(khroma)
  source('Functions_histM.R')
  
  conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
  load(file = paste0('../results/CT_merged_20220328/Rdata', '/combined_4histMarkers_overlapped55kATACseq_DE_fdr0.05.Rdata'))
  design = readRDS(file = paste0('../results/CT_merged_20220328/Rdata', '/histM_CT_design_info.rds'))
  
  yy0 = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_notOverlapped.positionalPeaks.rds'))
  yy1 = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_positionalPeaks.rds'))
  
  yy1 = yy1[, grep('mRep', colnames(yy1))]
  yy0 = yy0[, grep('mRep', colnames(yy0))]
  
  mm = match(rownames(yy1), rownames(DE.locus))
  DE.peaks = DE.locus[mm, ]
  ss = apply(DE.peaks, 1, sum)
  length(which(ss>0))
  
  # not consider H3K27ac in the main figure
  yy1 = yy1[, grep('H3K27ac', colnames(yy1), invert = TRUE)]
  
  
  ## call function for heamtap 
  source('Functions_plots.R')
  subclustering.postional.histM.postioinalAtacPeaks()
  
  ##########################################
  # dynamic histM peaks overlapped with stable atac peaks;
  ##########################################
  yy0 = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_notOverlapped.positionalPeaks.rds'))
  yy0 = yy0[, grep('mRep', colnames(yy0))]
  
  # conds_histM = c('H3K4me3', 'H3K27me3', 'H3K4me1', 'H3K27ac')
  conds_histM = c('H3K4me3', 'H3K27me3', 'H3K4me1')
  range <- 2.0
  
  yy = matrix(NA, nrow = nrow(yy0), ncol = length(conds_histM)*3)
  rownames(yy) = rownames(yy0)
  nms = c()
  for(n in 1:length(conds_histM)) # transform the data
  {
    jj0 = grep(conds_histM[n], colnames(yy0))
    
    test = yy0[, jj0]
    
    test = cal_sample_means(test, conds = paste0(conds_histM[n], c('_mUA', '_mLA', '_mHand')))
    test = t(apply(test, 1, cal_centering))
    test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
    yy[, c((3*n-2):(3*n))] = test
    nms = c(nms, paste0(conds_histM[n], c('_mUA', '_mLA', '_mHand')))
    #yy0[ ,jj0] = t(apply(yy0[,jj0], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
    
  }
  colnames(yy) = nms
  
  df = as.data.frame(sapply(colnames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'segments'
  rownames(df) = colnames(yy)
  
  sample_colors = c('springgreen4', 'steelblue2', 'gold2')
  annot_colors = list(segments = sample_colors)
  
  gaps_col = c(3, 6)
  
  callback = function(hc, mat){
    sv = svd(t(mat))$v[,8]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  pheatmap(yy, cluster_rows=TRUE,
           #cutree_rows = 4,
           show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           gaps_col = gaps_col,
           legend = TRUE,
           treeheight_row = 20,
           annotation_legend = TRUE, 
           #annotation_colors = annot_colors,
           clustering_callback = callback,
           width = 4, height = 12, 
           filename = paste0(outDir, '/Fig_S1F_heatmap_histoneMarker_DE_notoverlapped.with.atac.postionalPeaks.pdf'))
  
  
  ##########################################
  # highlight some positional genes and promoters
  ##########################################
  pp = data.frame(t(sapply(rownames(yy), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
  pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-5000, 5000), level = 'transcript')  
  
  ### highlight promoter peaks
  pp.annots = as.data.frame(pp.annots)
  promoter.sels = grep('Promoter', pp.annots$annotation)
  
  yy.sels = yy[promoter.sels, ]
  peaks.sels = pp.annots[promoter.sels, ]
  rownames(peaks.sels) = rownames(yy.sels)
  mm = match(rownames(peaks.sels), rownames(keep))
  peaks.sels$log2fc = apply(keep[mm, setdiff(grep('logFC.mHand.vs.mUA_', colnames(keep)), grep('H3K27ac', colnames(keep)))], 
                            1, function(x){x[which.max(abs(x))]})
  
  o1 = order(-abs(peaks.sels$log2fc))
  peaks.sels = peaks.sels[o1, ]
  yy.sels = yy.sels[o1, ]
  
  #tss$gene[which(rownames(tss) == 'AMEX60DD028208')] = 'PROD1'
  #tss$gene[which(rownames(tss) == 'AMEX60DD024424')] = 'MEIS3'
  
  peaks.sels[grep('HOXA13|MEIS|SHOX', peaks.sels$transcriptId), ]
  
  ## further selection of positional genes with microarray data
  genelists = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/', 
                                    'microarray_positionalGenes_data.rds'))
  ids = get_geneID(rownames(genelists)) 
  #ids = c(tss$geneID[match(positional.genes, tss$gene)], ids)
  ids = unique(ids)
  
  # convert transcriptID to gene ID
  annot = read.delim(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'annot_AmexT_v47_transcriptID_transcriptCotig_geneID_geneSymbol.hs_nr_CDS.txt'))
  
  peaks.sels$geneId = sapply(peaks.sels$geneId, function(x) {x = unlist(strsplit(as.character(x), '[|]')); 
  x[grep('AMEX60', x)]} )
  peaks.sels$geneId = as.character(annot$geneID[match(peaks.sels$geneId, annot$transcriptID)])
  
  kk = match(peaks.sels$geneId, ids)
  kk = which(!is.na(kk))
  kk = unique(c(grep('MEIS', peaks.sels$transcriptId), kk))
  
  #test = peaks.sels[kk, ]
  #test[grep('HOXA13|MEIS|SHOX', test$transcriptId), ]
  peaks.sels = peaks.sels[kk, ]
  yy.sels = yy.sels[kk, ]
  
  segments = c('mUA', 'mLA', 'mHand')
  conds = c('atac', 'H3K4me3','H3K27me3', 'H3K4me1')
  cc = paste(rep(conds, each = length(segments)), segments, sep = "_")
  
  ####### top 30 promoter peaks
  # ntop = 30
  # ntop = nrow(peaks.sels)
  # 
  # o1 = order(-peaks.sels$fdr.mean)
  # peaks.sels = peaks.sels[o1[1:ntop], ]
  # yy.sels = yy.sels[o1[1:ntop], ]
  
  geneSymbols = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                               'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  # jj = which(geneSymbols$geneID == 'AMEX60DD024424')
  # geneSymbols$gene.symbol.toUse[jj] = 'MEIS3'
  # geneSymbols$manual[jj] = 'MEIS3'
  # 
  # saveRDS(geneSymbols, file = paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
  #                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  ggs = peaks.sels$geneId
  mm = match(ggs, geneSymbols$geneID)
  ggs[!is.na(mm)] = geneSymbols$gene.symbol.toUse[mm[!is.na(mm)]]
  
  grep('HOX13|MEIS|SHOX', ggs)
  
  yy.sels = as.matrix(yy.sels)
  
  rownames(yy.sels) = ggs
  test = yy.sels
  
  df = data.frame(segments = sapply(colnames(test), function(x) unlist(strsplit(as.character(x), '_'))[2]))
  colnames(df) = c('seg')
  rownames(df) = colnames(test)
  
  sample_colors = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors) = c('mUA', 'mLA', 'mHand')
  annot_colors = list(segments = sample_colors)
  
  gaps.col = c(3, 6)
  pheatmap(test, 
           annotation_col = df, show_rownames = TRUE, scale = 'none', 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
           show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = FALSE, 
           annotation_colors = annot_colors, 
           gaps_col = gaps.col, fontsize_row = 8,
           #breaks = seq(-2, 2, length.out = 8),
           clustering_method = 'complete', cutree_rows = 4,
           breaks = seq(-range, range, length.out = 8),
           #gaps_row =  gaps.row, 
           treeheight_row = 20,
           legend_labels = FALSE,
           annotation_legend = FALSE,
           filename = paste0(outDir, 'Fig_S1G_heatmap_positional_histMarks_promoters.pdf'), 
           width = 5, height = 15)
  
}

##########################################
# characterize the correlations between atac-seq changes and histone markers
##########################################
library(tidyr)
library(dplyr)
require(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(ggpubr)
library("cowplot")
source('Functions_histM.R')

signals = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_positionalPeaks.rds'))
peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rdata/', 
                              'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))
ps = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rxxxx_R10723_R11637_R12810_atac/Rdata', 
                           '/positional_peakSignals.rds'))
segments = c('mUA', 'mLA', 'mHand')

xx1 = cal_sample_means(ps, conds =  c("Mature_UA", "Mature_LA", "Mature_Hand"))
colnames(xx1) = paste0('atac_', segments)
xx1 = data.frame(clusters = peaks$clusters, xx1, stringsAsFactors = FALSE)

# conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
# levels = c('atac', 'H3K4me3', 'H3K4me1', 'H3K27me3', 'H3K27ac')
conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1')
levels = c('atac', 'H3K4me3', 'H3K4me1', 'H3K27me3')

cc = paste(rep(conds_histM, each = length(segments)), segments, sep = "_")
xx2 = cal_sample_means(signals, conds = cc)
rownames(xx2) = rownames(xx1)

xx1 = data.frame(xx1, xx2, stringsAsFactors = FALSE)

cols = c('Blue', 'darkgreen',  'purple1', 'orange2')

# cluster order : 6, 1, 5, 3, 4, 2
yy = xx1[which(xx1$clusters == 6),]
p1 = as_tibble(yy) %>% 
  gather(marker, signal, 2:ncol(yy)) %>%
  separate(marker, c('markers', 'segment')) %>%
  ggplot(aes(x = segment, y = signal, fill = factor(markers, levels = levels))) +
  #geom_bar(stat = "identity") +
  geom_boxplot() + 
  #scale_fill_brewer(palette="Dark2") +
  scale_fill_manual(values = cols) +
  guides(fill=guide_legend(title="")) +
  theme_classic() +
  ggtitle('cluster 1')

yy = xx1[which(xx1$clusters == 1),]
p2 = as_tibble(yy) %>% 
  gather(marker, signal, 2:ncol(yy)) %>%
  separate(marker, c('markers', 'segment')) %>%
  ggplot(aes(x = segment, y = signal, fill = factor(markers, levels = levels))) +
  #geom_bar(stat = "identity") +
  geom_boxplot() + 
  #scale_fill_brewer(palette="Dark2") +
  scale_fill_manual(values = cols) +
  guides(fill=guide_legend(title="")) +
  theme_classic() +
  ggtitle('cluster 2')


yy = xx1[which(xx1$clusters == 5| xx1$clusters == 3),]
p3 = as_tibble(yy) %>% 
  gather(marker, signal, 2:ncol(yy)) %>%
  separate(marker, c('markers', 'segment')) %>%
  ggplot(aes(x = segment, y = signal, fill = factor(markers, levels = levels))) +
  #geom_bar(stat = "identity") +
  geom_boxplot() + 
  #scale_fill_brewer(palette="Dark2") +
  scale_fill_manual(values = cols) +
  guides(fill=guide_legend(title="")) +
  theme_classic() +
  ggtitle('cluster 3,4')

yy = xx1[which(xx1$clusters == 2| xx1$clusters == 4),]
p4 = as_tibble(yy) %>% 
  gather(marker, signal, 2:ncol(yy)) %>%
  separate(marker, c('markers', 'segment')) %>%
  ggplot(aes(x = segment, y = signal, fill = factor(markers, levels = levels))) +
  #geom_bar(stat = "identity") +
  geom_boxplot() + 
  #scale_fill_brewer(palette="Dark2") +
  scale_fill_manual(values = cols) +
  guides(fill=guide_legend(title="")) +
  theme_classic() +
  ggtitle('cluster 5,6') 
#scale_fill_brewer(palette="BuPu")
#geom_violin(width = 0.8) +
#scale_color_discrete(values=c('darkblue', 'forestgreen',  "red",   'orange', 'black')) + 
#geom_jitter(width = 0.2, size = 0.5) + 

#theme(legend.position = "none")  + 
# theme(axis.text.x = element_text(angle = 0))
plot_list=list()
plot_list[['p1']]=p1[[4]]
plot_list[['p2']]=p2[[4]]
plot_list[['p3']]=p3[[4]]
plot_list[['p4']]=p4[[4]]

pdf(paste0(figureDir, "/Fig_S1H_atac_histMarkers_boxplot_byClusters.pdf"),
    width = 10, height = 6) # Open a new pdf file

plot_grid(p1, p2, p3, p4,
          ncol = 2, nrow = 2)

dev.off()

##########################################
# group based on the fact if histone markers follow atac-seq changes, only with UA and Hand 
##########################################
peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rdata/', 
                              'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))
signals = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_positionalPeaks.rds'))
load(file = paste0(RdataDir, '/combined_4histMarkers_overlapped55kATACseq_DE_fdr0.05.Rdata'))

mm = match(rownames(signals), rownames(keep))
missed = which(is.na(mm))
cat(missed, ' -- ', rownames(signals)[missed], '\n')
grep('chr9p:41110738-41112452', rownames(signals))

rownames(keep)[mm[174]]
mm[missed] = mm[174]

keep = keep[mm, ]

conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1')
res = matrix(0, nrow = nrow(keep), ncol = (1+length(conds_histM)))
colnames(res) = c('atac', conds_histM)
rownames(res) = rownames(peaks)
res = data.frame(res, stringsAsFactors = FALSE)

res$atac = peaks$logFC.mHand.vs.mUA

fdr.cutoff = 0.05;
logfc.cutoff = 0;
marker.cutoff = 1;

for(n in 1:length(conds_histM))
{
  marker = conds_histM[n]
  eval(parse(text = paste0('sels = which(abs(keep$logFC.mHand.vs.mUA_', marker, 
                           ') > logfc.cutoff & keep$adj.P.Val.mHand.vs.mUA_', marker, ' < fdr.cutoff)')))
  
  eval(parse(text = paste0('res$', marker, '[sels] = keep$logFC.mHand.vs.mUA_', marker, '[sels]')))
  #eval(parse(text = paste0('res$', marker, ' = keep$logFC.mHand.vs.mUA_', marker))) 
  
}

yy = res

for(n in 1:ncol(yy)){
  yy[which(yy[,n] > 0),n] = 1
  yy[which(yy[,n] < 0),n] = -1
}

yy = t(yy)

# yy = yy[c(nrow(yy):1), ]

outDir = "/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/Figure_S1/"
pheatmap(yy, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = FALSE,
         color = c('darkred', 'white', 'darkgreen'), treeheight_row = 0, treeheight_col = 0, 
         filename = paste0(outDir, '/Fig_S1H_check_histMarkers_follow_atacPeak.pdf'), 
         width = 6, height = 2)


##########################################
# feature distribution of positional atac-seq peaks
##########################################
library('rtracklayer')
library(GenomicRanges)
library('GenomicFeatures')

## save gtf of axolotl limb fibroblast expressing genes
# annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
# annot = import(paste0(annotDir, 'AmexT_v47_Hox.patch.gtf'))
# 
# gene.expr = readRDS(file = '../data/expressedGenes_list_limb_fibroblast_using_smartseq2.mature.regeneration_pooledscRNAseq.dev.rds')
# id.expr = gene.expr$geneID[which(gene.expr$expressed>0)]
# 
# aa = annot[which(!is.na(match(annot$gene_id, id.expr)))]
# export(aa, con = paste0('../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'), 
#        format = 'gtf')
gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rdata/', 
                              'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))
ps = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rxxxx_R10723_R11637_R12810_atac/Rdata', 
                           '/positional_peakSignals.rds'))

mm = match(rownames(peaks), rownames(ps))
peaks = data.frame(ps[mm, ], peaks, stringsAsFactors = FALSE)

pp = data.frame(t(sapply(rownames(peaks), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
pp$strand = '*'
pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")
# annotation from ucsc browser ambMex60DD_genes_putative
amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

pdfname = paste0(figureDir, "Fig1_feature_distribution_positionalPeaks.pdf")
pdf(pdfname, width = 6, height = 4)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

plotPeakAnnot_piechart(pp.annots)

dev.off()

pp.annots = as.data.frame(pp.annots)
rownames(pp.annots) = rownames(peaks)

## update the peak annoation using limb-fibroblast expressing genes
peaks[, c(36:49)] = pp.annots

saveRDS(peaks, file = paste0(RdataDir, 
                             '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters6_DEtest_peakSignals_peakAnnot.updated.rds'))


##########################################
# Figure 1 promoter peaks and top enhancer peaks with atac-seq and histone markers
##########################################
signals = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_positionalPeaks.rds'))
load(file = paste0(RdataDir, '/combined_4histMarkers_overlapped55kATACseq_DE_fdr0.05.Rdata'))

peaks = readRDS(paste0(RdataDir, 
    '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters6_DEtest_peakSignals_peakAnnot.updated.rds'))
enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual_targets.rds'))

mm = match(rownames(signals), rownames(keep))
missed = which(is.na(mm))
cat(missed, ' -- ', rownames(signals)[missed], '\n')
grep('chr9p:41110738-41112452', rownames(signals))

rownames(keep)[mm[174]]
mm[missed] = mm[174]
keep = keep[mm, ]

grep('chr9p:41110738-41112452', rownames(keep))

xx = peaks[, c(1:14)]
colnames(xx) = paste0('atac_', colnames(xx))
colnames(xx) = gsub('Mature_UA', 'mUA', colnames(xx))
colnames(xx) = gsub('Mature_LA', 'mLA', colnames(xx))
colnames(xx) = gsub('Mature_Hand', 'mHand', colnames(xx))
colnames(xx) = gsub('HEAD', 'mHEAD', colnames(xx))

xx = data.frame(xx, signals, stringsAsFactors = FALSE)

source('Functions_histM.R')
segments = c('mUA', 'mLA', 'mHand')
conds = c('atac', 'H3K4me3','H3K27me3', 'H3K4me1')
cc = paste(rep(conds, each = length(segments)), segments, sep = "_")

yy = cal_sample_means(xx, cc)

### highlight promoter peaks
promoter.sels = grep('Promoter', peaks$annotation)

yy.sels = yy[promoter.sels, ]
peaks.sels = peaks[promoter.sels, ]

peaks.sels[grep('HOXA13|MEIS|SHOX|HOX', peaks.sels$transcriptId), ]

####### top 30 promoter peaks
ntop = 30
#ntop = nrow(peaks.sels)

o1 = order(-peaks.sels$fdr.mean)
peaks.sels = peaks.sels[o1[1:ntop], ]
yy.sels = yy.sels[o1[1:ntop], ]

geneSymbols = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                             'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

ggs = peaks.sels$geneId
mm = match(ggs, geneSymbols$geneID)
ggs[!is.na(mm)] = geneSymbols$gene.symbol.toUse[mm[!is.na(mm)]]

grep('HOX', ggs)

yy.sels = as.matrix(yy.sels)

if(saveTables){
  test = data.frame(peak = rownames(peaks.sels), gene = ggs, geneID = peaks.sels$geneId,  
                    yy.sels, stringsAsFactors = FALSE)
  test = test[, c(1:6)]
  test$specificSegment = 'mHand'
  test$specificSegment[which(test$atac_mUA > test$atac_mHand)] = 'mUA'
  test = test[, -c(4:6)]
  
  write.csv(test, file = paste0(tableDir, 'Figure1_allPromoter_segment_specific_atacPeak.csv'), 
              row.names = FALSE)
  
}

rownames(yy.sels) = ggs
test = yy.sels
for(n in 1:length(conds))
{
  ii = grep(conds[n], colnames(test))
  test[, ii] = t(apply(test[,ii], 1, cal_centering))
}

range <- 2.0
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))

df = data.frame(segments = sapply(colnames(test), function(x) unlist(strsplit(as.character(x), '_'))[2]))
colnames(df) = c('seg')
rownames(df) = colnames(test)

sample_colors = c('springgreen4', 'steelblue2', 'gold2')
names(sample_colors) = c('mUA', 'mLA', 'mHand')
annot_colors = list(segments = sample_colors)

gaps.col = c(3, 6, 9)
pheatmap(test, 
         annotation_col = df, show_rownames = TRUE, scale = 'none', 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
         show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = FALSE, 
         annotation_colors = annot_colors, 
         gaps_col = gaps.col, fontsize_row = 10,
         #breaks = seq(-2, 2, length.out = 8),
         breaks = seq(-range, range, length.out = 8),
         #gaps_row =  gaps.row, 
         treeheight_row = 20,
         legend_labels = FALSE,
         annotation_legend = FALSE,
         filename = paste0(figureDir, '/heatmap_positionalPeaks_top30.promoters.pdf'), 
         width = 6, height = 5)

saveRDS(ggs, file = paste0(RdataDir, '/positional_atac_promoters_top30.rds'))

##########################################
# #### highlight the enhancers
##########################################
peaks$targets = enhancers$targets[match(rownames(peaks), rownames(enhancers))]
sels = intersect(grep('Intergenic|Intron', peaks$annotation), which(!is.na(peaks$targets)))
ntop = 50

yy.sels = yy[sels, ]
peaks.sels = peaks[sels, ]

peaks.sels[grep('HOXA13|MEIS|SHOX|HOX', peaks.sels$transcriptId), ]

## sort the promoters with fdr and takes the top30
o1 = order(-peaks.sels$fdr.mean)
peaks.sels = peaks.sels[o1[1:ntop], ]
yy.sels = yy.sels[o1[1:ntop], ]

## use the updated peak-to-gene assignment
geneSymbols = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                             'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
ggs = peaks.sels$targets
mm = match(ggs, geneSymbols$geneID)
ggs[!is.na(mm)] = geneSymbols$gene.symbol.toUse[mm[!is.na(mm)]]

grep('HOX', ggs)

rownames(yy.sels) = ggs
yy.sels = as.matrix(yy.sels)

test = yy.sels
for(n in 1:length(conds))
{
  ii = grep(conds[n], colnames(test))
  test[, ii] = t(apply(test[,ii], 1, cal_centering))
}

range <- 2.0
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))

df = data.frame(segments = sapply(colnames(yy.sels), function(x) unlist(strsplit(as.character(x), '_'))[2]))
colnames(df) = c('seg')
rownames(df) = colnames(yy.sels)

sample_colors = c('springgreen4', 'steelblue2', 'gold2')
names(sample_colors) = c('mUA', 'mLA', 'mHand')
annot_colors = list(segments = sample_colors)

saveRDS(rownames(test), file = paste0(RdataDir, '/positional_atac_enhancers_top50.rds'))

gaps.col = c(3, 6, 9)
pheatmap(test, 
         annotation_col = df, show_rownames = TRUE, scale = 'none', 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8),  
         show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = FALSE, 
         annotation_colors = annot_colors, 
         gaps_col = gaps.col, fontsize_row = 10,
         breaks = seq(-range, range, length.out = 8),
         treeheight_row = 20,
         legend_labels = FALSE,
         annotation_legend = FALSE,
         #gaps_row =  gaps.row, 
         filename = paste0(figureDir, '/heatmap_positionalPeaks_top50.enhancers.pdf'), 
         width = 6, height = 10)

##########################################
# distance distribution of between positional peaks and positional genes
##########################################
# peaks = readRDS(paste0(RdataDir, 
#     '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters6_DEtest_peakSignals_peakAnnot.updated.rds'))
# 
# pp = data.frame(t(sapply(rownames(peaks), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
# pp$strand = '*'
# pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
#                               start.field="X2", end.field="X3", strand.field="strand")
# 
# gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
# aa = import(gtf.file)
# 
# # load the positional gene list
# genes = readRDS(file = '../data/positional_genes.716_microarray.rds')
# ids = sapply(rownames(genes), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
# kk = unique(c(grep('^HOX', aa$gene_id), which(!is.na(match(aa$gene_id, ids)))))
# 
# aa = aa[kk, ]
# 
# amex = GenomicFeatures::makeTxDbFromGRanges(aa)
# pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
# 
# pp.annots = as.data.frame(pp.annots)
# rownames(pp.annots) = rownames(peaks)
# 
# dists = data.frame(peak.to.expr = peaks$distanceToTSS, peak.to.positional = pp.annots$distanceToTSS, stringsAsFactors = FALSE)
# 
# dists = apply(dists, 2, function(x) {sign(x) * log10(abs(x)+1)})
# 
# library(tidyr)
# library(dplyr)
# require(ggplot2)
# library(gridExtra)
# library(grid)
# library(lattice)
# library(ggpubr)
# library("cowplot")
# as_tibble(dists) %>% 
#   gather(geneGroup, dist, 1:2) %>%
#   ggplot(aes(x=dist, color = geneGroup)) + 
#   geom_density(size = 1.) +
#   theme_classic() +
#   scale_color_brewer(palette="Dark2") +
#   labs(x = 'distance to closest TSS (log10 bp)') +
#   geom_vline(xintercept=c(-7, -5.5, 5.5, 7), col='gray', size = 1.) +
#   theme(legend.text = element_text(size=10),
#         legend.title = element_text(size = 10),
#         legend.position=c(0.5, 0.8),
#         plot.margin = margin()
#         #legend.key.size = unit(1, 'cm')
#         #legend.key.width= unit(1, 'cm')
#   )
# 
# ggsave(paste0(figureDir, "distance_to_closest_genes_positionalPeaks.pdf"), width=4, height = 3)

########################################################
########################################################
# Section III : positiona features in the regenerations 
# TSS of postional genes
# positional enhancers 
########################################################
########################################################
source('Functions_histM.R')
source('Functions_Integration.matureReg.R')
library(ggrepel)
library(dplyr)
library(tibble)
library("cowplot")
require(gridExtra)
library(tidyr)
require(patchwork)

## bivalent analysis here
tss = readRDS(file = paste0(RdataDir, '/regeneration_matureSamples_tss_perGene_smartseq2_atac_histM_v4.rds'))

tss$gene[which(rownames(tss) == 'AMEX60DD028208')] = 'PROD1'
tss$gene[which(rownames(tss) == 'AMEX60DD024424')] = 'MEIS3'

# saveRDS(tss, file =paste0(RdataDir, '/regeneration_matureSamples_tss_perGene_smartseq2_atac_histM_v5.rds'))

kk = which(!is.na(tss$gene))
rownames(tss)[kk] = paste0(tss$gene[kk], '_', rownames(tss)[kk])
tss$gene[-kk] = rownames(tss)[-kk]
positional.genes = c('HOXA13','PROD1', 'RARRES1', 'MEIS1', 'MEIS2', 'SHOX', 'SHOX2',  'HOXA11', 'HOXA9', 'HOXD13',
                     'HOXD11', 'HOXD9')

##########################################
# plot individual gene examples of different features, RNAseq, atac, histone marks around TSS 
##########################################
Plot.gene.examples = FALSE
if(Plot.gene.examples){
  outDir = "/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/Gene_Examples"   
  source('Functions_Integration.matureReg.R')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  plot_rna_chromainFeatures_geneExamples(tss, geneList = positional.genes, outDir = outDir, incl_Mature = TRUE, log2fc = TRUE)
  
}

##########################################
# TSS of positional genes (known and from microarray data) in mature and regeneration samples 
##########################################
source('Functions_Integration.matureReg.R')
source('Functions_histM.R')

genelists = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/', 
                                             'microarray_positionalGenes_data.rds'))

ids = get_geneID(rownames(genelists)) 
ids = c(tss$geneID[match(positional.genes, tss$gene)], ids)
ids = unique(ids)

#test = Analysis_TSS_positionalGenes_in_mature_regeneration(tss, ids)
#test = readRDS(file = paste0(RdataDir, '/positional_gene_TSS_chromatinFeatures.rds'))
test = readRDS(file = paste0(RdataDir, '/positional_gene_TSS_chromatinFeatures_smartseq2_mature.reg.rds'))
xx = as.matrix(test[,c(1:35)])
yy = as.matrix(test[, c(38:ncol(test))])
rownames(xx) = test$gene
rownames(yy) = rownames(xx)

## remove the H3K27ac in the main figure
xx = xx[, grep('H3K27ac', colnames(xx), invert = TRUE)]
range <- 3.0
xx = t(apply(xx, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
range = 6.
yy = t(apply(yy, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
# pheatmap(xx, 
#          #annotation_col = df, 
#          show_rownames = FALSE, scale = 'none', 
#          color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(7),
#          show_colnames = FALSE,
#          cluster_rows = TRUE, 
#          cluster_cols = FALSE, 
#          #annotation_colors = annot_colors, 
#          gaps_col = seq(7, 35, by = 7), 
#          fontsize_row = 8,
#          treeheight_row = 50,
#          cutree_rows = 8,
#          #gaps_row =  gaps.row, 
#          filename = paste0(figureDir, '/heatmap_positionalGens_TSS_mature_regeneration_all.pdf'), 
#          width = 6, height = 15)

maxs = apply(xx, 1, function(x) max(abs(x)))
tops = which(maxs>2)
xx = xx[tops, ]
yy = yy[tops, ]

df = as.data.frame(sapply(colnames(xx), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
colnames(df) = 'samples'
rownames(df) = colnames(xx)

sample_colors = c('gold2',  'steelblue2', 'springgreen4',   'springgreen3', 'magenta', 'darkblue', 'red')
names(sample_colors) = unique(df$samples)
annot_colors = list(samples = sample_colors)

plt = pheatmap(xx, 
         annotation_col = df, 
         show_rownames = FALSE, scale = 'none', 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(7),
         show_colnames = FALSE,
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         annotation_colors = annot_colors, 
         gaps_col = seq(7, 28, by = 7), 
         fontsize_row = 7,
         treeheight_row = 20,
         cutree_rows = 8,
         #gaps_row =  gaps.row, 
         filename = paste0(figureDir, '/heatmap_positionalGens_TSS_mature_regeneration_log2FC.2.pdf'), 
         width = 6, height = 14)

df = as.data.frame(sapply(colnames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
colnames(df) = 'samples'
rownames(df) = colnames(yy)

source('Functions_histM.R')
gaps.row = cal_clusterGaps(plt, nb_clusters = 8)

pheatmap(yy[plt$tree_row$order, ], 
         annotation_col = df, 
         show_rownames = TRUE, 
         scale = 'none', 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12),
         show_colnames = FALSE,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         annotation_colors = annot_colors, 
         #gaps_col = seq(7, 35, by = 7), 
         fontsize_row = 7,
         #treeheight_row = 40,
         #cutree_rows = 8,
         gaps_row =  gaps.row, 
         filename = paste0(figureDir, '/heatmap_positionalGens_smartseq2_mature_regeneration_log2FC.2.pdf'), 
         width = 3.2, height = 14)


## predicted postional genes from above clusters 
outDir = paste0("/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/", 
                "figure/plots_4figures/Gene_Examples")

source('Functions_Integration.matureReg.R')
if(!dir.exists(outDir)) dir.create(outDir)

tss$gene[which(tss$geneID == 'AMEX60DD023963')] = 'AMEX60DD023963'
tss$gene[which(tss$geneID == 'AMEX60DD004776')] = 'AMEX60DD004776'
plot_rna_chromainFeatures_geneExamples(tss, geneList = c('LHX2', 'COL9A2', 'GDF5', 'SALL1', 'SCUBE2', 
                                                         'ECM2', 'CYP2F1', 'CYP2A13'), 
                                       outDir = outDir, incl_Mature = TRUE, log2fc = TRUE)

##########################################
# Search for potential positional genes that are not expressed in mature samples,
# but with positional chromatin features; e.g. SHH in A-P position
# a few criterion should be met:
# 1) Regeneration-response genes
# 2) showing difference in chromatin features (ATAC, histone marks) in mature samples or dpa13 proximal and distal
##########################################
source('Functions_Integration.matureReg.R')
source('Functions_histM.R')

# start with regeneration RNA-seq analysis
geneClusters = readRDS(file = paste0("../results/RNAseq_data_used/Rdata/", 'regeneration_geneClusters.rds'))
geneClusters$gene = rownames(geneClusters)
geneClusters$geneID = sapply(rownames(geneClusters), 
                             function(x) {test = unlist(strsplit(as.character(x), '_')); return(test[length(test)])})

dpcompare = readRDS(file = paste0("../results/RNAseq_data_used/Rdata/", 
                            'smartseq2_regeneration_dpa13dist.vs.dpa13prox.rds'))

ggs = rownames(dpcompare)[which(abs(dpcompare$log2FoldChange)> 1 & dpcompare$padj <0.05)]
ggs = get_geneID(ggs)

genelists = intersect(geneClusters$geneID[which(geneClusters$mUA< (-1))], ggs)

genelist2 = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/', 
                                    'microarray_positionalGenes_data.rds'))
ids = get_geneID(rownames(genelist2))
genelists = unique(setdiff(genelists, ids))
genelists = genelists[grep('HOXA|HOXB|HOXC|HOXD', genelists, invert = TRUE)]

rm(genelist2)
cat(length(genelists), ' gene candidates to search \n')

mm = match(genelists, geneClusters$geneID)
as.character(get_geneName(rownames(geneClusters)[mm]))

fdr.cutoff= 0.01
select = which(tss$atac.M_adj.P.Val.mHand.vs.mLA < fdr.cutoff |
                 tss$atac.M_adj.P.Val.mHand.vs.mUA < fdr.cutoff |
                 tss$H3K4me3.M_adj.P.Val.mLA.vs.mUA < fdr.cutoff |
                 tss$H3K4me3.M_adj.P.Val.mHand.vs.mUA < fdr.cutoff |
                 tss$H3K27me3.M_adj.P.Val.mLA.vs.mUA < fdr.cutoff |
                 tss$H3K27me3.M_adj.P.Val.mHand.vs.mUA < fdr.cutoff |
                 tss$H3K4me1.M_adj.P.Val.mHand.vs.mUA < fdr.cutoff |
                 tss$H3K4me1.M_adj.P.Val.mLA.vs.mUA < fdr.cutoff
                 )
ids = tss$geneID[select]
genelists = intersect(genelists, ids)

mm = match(genelists, geneClusters$geneID)
as.character(get_geneName(rownames(geneClusters)[mm]))

source('Functions_Integration.matureReg.R')

test = Analysis_TSS_positionalGenes_in_mature_regeneration(tss = tss, ids = genelists)
  
#test = Analysis_TSS_positionalGenes_in_mature_regeneration(tss, ids)
#test = readRDS(file = paste0(RdataDir, '/positional_gene_TSS_chromatinFeatures.rds'))
#test = readRDS(file = paste0(RdataDir, '/positional_gene_TSS_chromatinFeatures_smartseq2_mature.reg.rds'))
xx = as.matrix(test[,c(1:28)])

yy = as.matrix(test[, c(38:ncol(test))])
rownames(xx) = test$gene
rownames(yy) = rownames(xx)

quantile(xx, c(0.05, 0.1,  0.90,  0.95, 0.99), na.rm = TRUE)
range <- 2
xx = t(apply(xx, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))

quantile(yy, c(0.05, 0.1,  0.90,  0.95, 0.99), na.rm = TRUE)
range = 5
yy = t(apply(yy, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))

kk = which(rownames(xx) == 'AMEX60DD050556')
xx = xx[-kk, ]
yy = yy[-kk, ]

yy[which(is.na(yy[, 1])), 1] = 0
yy[which(is.na(yy[, 2])), 2] = 0

kk = apply(xx, 1, function(x){return(all(!is.na(x)))})
xx = xx[kk, ]
yy = yy[kk, ]

df = as.data.frame(sapply(colnames(xx), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
colnames(df) = 'samples'
rownames(df) = colnames(xx)

sample_colors = c('gold2',  'steelblue2', 'springgreen4',   'springgreen3', 'magenta', 'darkblue', 'red')
names(sample_colors) = unique(df$samples)
annot_colors = list(samples = sample_colors)

plt = pheatmap(xx, 
               annotation_col = df, 
               show_rownames = FALSE, scale = 'none', 
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(7),
               show_colnames = FALSE,
               cluster_rows = TRUE, 
               cluster_cols = FALSE, 
               annotation_colors = annot_colors, 
               gaps_col = seq(7, ncol(xx), by = 7), 
               fontsize_row = 7,
               treeheight_row = 20,
               cutree_rows = 4,
               #gaps_row =  gaps.row, 
               legend = FALSE,
               annotation_legend = FALSE,
               filename = paste0(figureDir, '/heatmap_positionalGens_notExprmUA_TSS_mature_regeneration_log2FC.2.pdf'), 
               width = 6, height = 5)

pheatmap(xx, 
         annotation_col = df, 
         show_rownames = FALSE, scale = 'none', 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(7),
         show_colnames = FALSE,
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         annotation_colors = annot_colors, 
         gaps_col = seq(7, ncol(xx), by = 7), 
         fontsize_row = 7,
         treeheight_row = 20,
         cutree_rows = 4,
         #gaps_row =  gaps.row, 
         legend = TRUE,
         annotation_legend = FALSE,
         filename = paste0(figureDir, '/heatmap_positionalGens_notExprmUA_TSS_mature_regeneration_log2FC.2_lenged.pdf'), 
         width = 6, height = 5)

df = as.data.frame(sapply(colnames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
colnames(df) = 'samples'
rownames(df) = colnames(yy)

source('Functions_histM.R')
gaps.row = cal_clusterGaps(plt, nb_clusters = 4)

pheatmap(yy[plt$tree_row$order, ], 
         annotation_col = df, 
         show_rownames = TRUE, 
         scale = 'none', 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12),
         show_colnames = FALSE,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         annotation_colors = annot_colors, 
         #gaps_col = seq(7, 35, by = 7), 
         fontsize_row = 7,
         #treeheight_row = 40,
         #cutree_rows = 8,
         gaps_row =  gaps.row, 
         filename = paste0(figureDir, '/heatmap_positionalGens_notExprmUA_smartseq2_mature_regeneration_log2FC.2.pdf'), 
         width = 3.5, height = 5)

##########################################
# positional enhancers in mature and regeneration 
# not sure if it is necessary
##########################################

########################################################
########################################################
# Section : backup analysis
# 
########################################################
########################################################
##########################################
# first motif activity analysis for positional-dependent peaks 
##########################################
source('MARA_functions.R')

saveRDS(keep, file = paste0(RdataDir, '/matrix.saved.for.spatial.MARA.rds'))

# prepare step is in the MARA_functions.R
xx = run.MARA.atac.spatial(keep, cc)

##########################################
# ## test the distribution of features of different groups
##########################################
source('Functions_atac.R')

pdfname = paste0(figureDir, "/Fig1E_positional_peak_feature_distribution_cluster_all.pdf")
pdf(pdfname, width = 8, height = 6)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

stats = plotPeakAnnot_piechart(pp.annots)

dev.off()

colnames(stats)[ncol(stats)] = paste0('all_', pp.annots@peakNum)

scores = c()
for(m in 1:nb_clusters)
{
  # m = 1
  cat('cluster - ',  m, '\n')
  pp.annots = annotatePeak(pp[which(xx$cluster == m)], TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  
  pdfname = paste0(figureDir, "/Fig1E_positional_peak_feature_distribution_cluster_", m, ".pdf")
  pdf(pdfname, width = 8, height = 6)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  stats = data.frame(stats,  plotPeakAnnot_piechart(pp.annots)[, 2])
  colnames(stats)[ncol(stats)] = paste0('cluster', m, '_', pp.annots@peakNum)
  
  test = c()
  test2 = c()
  for(i in 1:nrow(stats))
  {
    total = length(pp); 
    mm = round(total * stats[i, 2]/100)
    nn = total - mm
    qq = round(pp.annots@peakNum * stats[i, ncol(stats)]/100)
    test = c(test, phyper(qq-1, mm, nn, pp.annots@peakNum, lower.tail = FALSE, log.p = FALSE))
    test2 = c(test2, phyper(qq+1, mm, nn, pp.annots@peakNum, lower.tail = TRUE, log.p = FALSE))
  }
  scores = cbind(scores, test, test2)
  colnames(scores)[c(ncol(scores)-1, ncol(scores))] = c(paste0('enrich.pval.cluster', m), paste0('depelet.pval.cluster', m))
  
  dev.off()
}
rownames(scores) = rownames(scores)

stats = data.frame(stats, scores, stringsAsFactors = FALSE)

write.csv(stats, file = paste0(figureDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters', 
                               nb_clusters, '_feature.Enrichment.depeletion.csv'), quote = FALSE, row.names = TRUE)


library(tidyr)
library(dplyr)
require(ggplot2)

as_tibble(stats) %>%  gather(group, freq,  2:ncol(stats)) %>% 
  ggplot(aes(fill=group, y=freq, x=Feature)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, size = 10)) 

ggsave(paste0(figureDir, "/Fig1E_positional_peak_feature_distribution_cluster_comparison.pdf"),  width = 12, height = 8)

##########################################
# compare the ATAC-seq peaks and microarray data
# integrative analysis of positional peaks and positional mRNAs
# either from positional peak, check the potentially regulated target is position dependent
# or from postional mRNA, are there assigned peaks are positional dependent
##########################################
## first check the targeted gene of positional peaks have positio-dependent gene expression
# positonal peaks
#load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl', version.analysis, '.Rdata')) 
#saveRDS(peaks, file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl_signals_peakAnnot_', 
#                             version.analysis, '.rds')) 
peaks = readRDS(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl_signals_peakAnnot_', 
                              version.analysis, '.rds'))

annot = readRDS(paste0(annotDir, 
                       'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))

# positional gene 
rna = readRDS(file = paste0("../results/microarray/Rdata/", 
                            'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_limma.DE.stats.rds'))

matures = data.frame(ua = apply(rna[, c(1:3)], 1, mean), 
                     la = apply(rna[, c(4:6)], 1, mean), 
                     hd = apply(rna[, c(7:9)], 1, mean))

rna$maxs = apply(matures, 1, max)

qv.cutoff.rna = 0.05
logfc.cutoff.rna = 1

peaks$transcript = sapply(peaks$transcriptId, function(x) {x = unlist(strsplit(as.character(x), '[|]')); return(x[length(x)])})

peaks$Gene = annot$geneID[match(peaks$transcript, annot$transcriptID)]


ggs = rownames(rna)
ggs = sapply(ggs, function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})

select = which((rna$adj.P.Val_mHand.vs.mLA < qv.cutoff.rna & abs(rna$logFC_mHand.vs.mLA) > logfc.cutoff.rna|
                  rna$adj.P.Val_mHand.vs.mUA < qv.cutoff.rna & abs(rna$logFC_mHand.vs.mUA) > logfc.cutoff.rna |
                  rna$adj.P.Val_mHand.vs.mLA < qv.cutoff.rna & abs(rna$logFC_mHand.vs.mLA) > logfc.cutoff.rna) )

pos.genes = rna[select, ]
saveRDS(pos.genes, file = paste0('../data/positional_genes.716_microarray.rds'))

## VennDiagram showing the overlapp between positional peak targets and positional genes
library(VennDiagram)
set1 <- unique(peaks$Gene)
set2 = unique(ggs[select])

# Chart
venn.diagram(
  x = list(set1, set2),
  category.names = c("targets of positional peak" , "positional mRNA"),
  filename = paste0(figureDir, '/comparison_postionalPeak.targets_positionalGenes.png'),
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c('darkgreen', 'darkblue'),
  output=TRUE,
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.pos = c(-10, 10)
)


# heatmap showing the comparions of overlapping peaks
mm = match(peaks$Gene, ggs)
### note: some peak targets are not found in microarray data due to some reason, e.g. probe design or inaccurate peak gene assignment
#mm = unique(mm[!is.na(mm)])
yy = peaks[!is.na(mm), ]
mapped = rna[mm[!is.na(mm)], ]
#targets = rownames(rna)[mm]

select = which(mapped$adj.P.Val_mHand.vs.mLA < qv.cutoff.rna & abs(mapped$logFC_mHand.vs.mLA) > logfc.cutoff.rna|
                 mapped$adj.P.Val_mHand.vs.mUA < qv.cutoff.rna & abs(mapped$logFC_mHand.vs.mUA) > logfc.cutoff.rna |
                 mapped$adj.P.Val_mHand.vs.mLA < qv.cutoff.rna & abs(mapped$logFC_mHand.vs.mLA) > logfc.cutoff.rna )

yy = yy[select,  ]
mapped = mapped[select, ]

yy = yy[, c(1:14)]
mapped = mapped[, c(1:9)]

rep.sels = grep('HEAD|102657|102655|74938', colnames(keep), invert = TRUE)
yy = yy[, rep.sels]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

yy <- t(apply(yy, 1, cal_z_score))
mapped = t(apply(mapped, 1, cal_z_score))

yy = data.frame(yy, mapped)  

df <- data.frame(rep(rep(conds[1:3], each = 3), 2))
rownames(df) = colnames(yy)
colnames(df) = 'segments'
df$modality = rep(c('atac', 'mRNA'), each = 9)

col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
          "#33a02c", "#fb9a99", "#e31a1c",
          "#fdbf6f", "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99", "#b15928")

sample_colors = c('springgreen4', 'steelblue2', 'gold2')
names(sample_colors) = c('Mature_UA', 'Mature_LA', 'Mature_Hand')
data_colors = c('darkgreen', 'darkblue')
names(data_colors) = c('atac', 'mRNA')
annot_colors = list(
  segments = sample_colors,
  modality = data_colors)

gaps.col = c(9)

names = sapply(rownames(mapped), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[1])})
rownames(yy) = paste0(names, '_', rownames(yy))

pheatmap(yy, 
         annotation_col = df, show_rownames = TRUE, scale = 'none', 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlGn")))(8), 
         show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = FALSE,  
         clustering_method = 'complete', cutree_rows = 4, 
         annotation_colors = annot_colors, 
         gaps_col = gaps.col, 
         fontsize_row = 8,
         #gaps_row =  gaps.row, 
         filename = paste0(figureDir, '/heatmap_positionalPeaks_positionalGenes_overlapped.pdf'), 
         width = 12, height = 16)



##########################################
# histone markers and atac-seq of positional genes
# 
##########################################
Explore_histM_positional.genes = FALSE
if(Explore_histM_positional.genes){
  fpm = readRDS(file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb.rds'))
  res = data.frame(log2(fpm + 1), stringsAsFactors = FALSE)
  
  res$gene = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[2]})
  
  res$x1 = res$Hand_K4me3
  res$x2 = res$UA_K4me3
  rr = res$x1 - res$x2
  
  examples.sel = which(abs(rr)>2 & (res$x1>5|res$x2>5))
  examples.sel = unique(c(examples.sel, grep('HOXA11|HOXA9|HOXD13|HOXD11|HOXD9|MEIS', res$gene)))
  
  ggplot(data=res, aes(x=x1, y=x2, label = gene)) +
    geom_point(size = 0.6) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
    #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=5, col='darkgray') +
    geom_hline(yintercept=5, col="darkgray") +
    labs(x = "Hand_H3K4me3", y= 'UA_H3K4me3')
  
  ggsave(paste0(figureDir, "histMarker_H3K4me3_scatterplot_Hand.vs.Hand.pdf"), width=12, height = 8)
  
  res$x1 = res$Hand_K27me3
  res$x2 = res$UA_K27me3
  rr = res$x1 - res$x2
  
  examples.sel = which(abs(rr)>2 & (res$x1>5|res$x2>5))
  examples.sel = unique(c(examples.sel, grep('HOXA11|HOXA9|HOXD13|HOXD11|HOXD9|MEIS', res$gene)))
  
  ggplot(data=res, aes(x=x1, y=x2, label = gene)) +
    geom_point(size = 0.6) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
    #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=5, col='darkgray') +
    geom_hline(yintercept=5, col="darkgray") +
    labs(x = "Hand_H3K27me3", y= 'UA_H3K27me3')
  
  ggsave(paste0(figureDir, "histMarker_H3K27me3_scatterplot_Hand.vs.Hand.pdf"), width=12, height = 8)
  
}
