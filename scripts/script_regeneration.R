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

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0(figureDir, 'tables4plots/')

saveTables = FALSE

require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require(pheatmap)
library(tictoc)

# limb fibroblast expressing genes
gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)

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

load(file = paste0(RdataDir, '/ATACseq_selected.55k.peaks_cutoff.50.at.least.1sample.Rdata'))

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
# Section II : identify dynamic peaks in regeneration
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
  
  gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
  
  ##########################################
  ## make Granges and annotate peaks
  ##########################################
  Make.Granges.and.peakAnnotation = TRUE
  if(Make.Granges.and.peakAnnotation){
    
    pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
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
    
    # promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
    
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
    
    # saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks_mUA_regeneration_dev_2Batches.R10723_R7977_v8_noPeakAnnot.rds'))
    
    xx = data.frame(cpm, res, pp.annots[match(rownames(cpm), rownames(pp.annots)), ],  stringsAsFactors = FALSE)
    
    res = xx
    
    #res = data.frame(cpm, res, stringsAsFactors = FALSE)
    saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks__mUA_regeneration_dev_2Batches.R10723_R7977_peakAnnot_v8.rds'))
    
  }
}


########################################################
########################################################
# Section : characterize regeneration chromatin landscape
# - all dynamic peaks (atac, histone markers) 
# - regeneration-specific peaks (atac, histone markers)
########################################################
########################################################
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
source('Functions_histM.R')
source('Functions_atac.R')
conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", 'Embryo_Stage44_distal', 
          "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", 'BL_UA_13days_distal'
)

sample.means = cal_sample_means(res[, c(1:20)], conds = conds)

res$mins =  apply(sample.means[match(rownames(res), rownames(sample.means)), ], 1, min)
res$max = apply(sample.means[match(rownames(res), rownames(sample.means)), ], 1, max)
# xx = xx[which(xx$min < 4.5), ]

#length(which(res$min<3))
length(which(res$max > 3))
res = res[which(res$max > 3), ]
cat(nrow(res), 'DE peaks found with at least condition above background !\n')

# saveRDS(res, file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration.rds'))

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
sample.means = sample.means[match(rownames(keep), rownames(sample.means)), ]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

yy <- t(apply(sample.means, 1, cal_z_score))

df <- data.frame(conds)
rownames(df) = colnames(yy)
colnames(df) = 'sample'

sample_colors = c('magenta', 'darkblue', 'springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2',
                  'red')[c(1:length(conds))]

names(sample_colors) = conds
col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
          "#33a02c", "#fb9a99", "#e31a1c",
          "#fdbf6f", "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99", "#b15928")

Regeneration.peaks.clustering.DPGP = FALSE
if(!Regeneration.peaks.clustering.DPGP)
{
  # specify the nb of clusters
  library(dendextend)
  nb_clusters = 8
  my_hclust_gene <- hclust(dist(yy), method = "complete")
  
  my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
  my_gene_col <- data.frame(cluster =  paste0('cluster_', my_gene_col))
  rownames(my_gene_col) = rownames(yy)
  
  res$clusters = my_gene_col
  
  cluster_col = col3[1:nb_clusters]
  names(cluster_col) = paste0('cluster_', c(1:nb_clusters))
  
  annot_colors = list(
    condition = sample_colors,
    cluster = cluster_col)
  
  ii.gaps = c(3, 4)
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
           filename = paste0(figureDir, 'heatmap_regenerationPeaks_scaled.pdf'), 
           width = 6, height = 12)
  
  plt = pheatmap(yy, 
                 annotation_row = my_gene_col, 
                 annotation_col = df, show_rownames = FALSE, scale = 'none', 
                 color = col, 
                 show_colnames = FALSE,
                 cluster_rows = TRUE, cluster_cols = FALSE,  
                 clustering_method = 'complete', cutree_rows = nb_clusters, 
                 annotation_colors = annot_colors, 
                 #clustering_callback = callback,
                 gaps_col = ii.gaps)
  
  save(yy, plt, file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap.Rdata'))
  
  
}else{
  
  #res2 = process.dynamic.peaks.clustering.GPDP(yy, res)
  # reload the processed GPDP clusters
  res = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds'))
  res = res[which(!is.na(res$clusters)), ]
  yy = yy[match(rownames(res), rownames(yy)), ]
  
  ######
  ### given the order of groups and order the peaks withinin groups using hclust
  ## find the row gap and order subclusters   
  
  #res$clusters = as.character(res$clusters)
  cluster_order = paste0('mc', c(1:8))
  
  gaps.row = c()
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  
  for(n in 1:length(cluster_order))
  {
    kk = which(res$cluster == cluster_order[n])
    cat('clsuter ', cluster_order[n], ' -- ', length(kk), ' peaks \n')
    
    hm_hclust <- hclust(dist(as.matrix(yy[kk,])), method = "complete")
    #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
    peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
    
    if(n == 1)  {gaps.row = c(gaps.row, length(which(res$clusters == cluster_order[n])))
    }else{
      if(n < length(cluster_order)) {
        gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(res$clusters == cluster_order[n])))
      }
    }
  }
  
  new_order = match(peakNm, rownames(yy))
  
  test = yy[new_order, ]
  
  my_gene_col <- data.frame(cluster = res$clusters[match(rownames(test), rownames(res))])
  rownames(my_gene_col) = rownames(test)
  
  cluster_col = col3[1:length(cluster_order)]
  names(cluster_col) = cluster_order
  
  annot_colors = list(
    condition = sample_colors,
    cluster = cluster_col)
  
  ii.gaps = c(3, 4)
  col = colorRampPalette(c("navy", "white", "red3"))(16)
  
  
  pheatmap(test, 
           annotation_row = my_gene_col, 
           annotation_col = df, show_rownames = FALSE, scale = 'none', 
           color = col, 
           show_colnames = FALSE,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           #clustering_method = 'complete', cutree_rows = nb_clusters, 
           annotation_colors = annot_colors, 
           #clustering_callback = callback,
           gaps_col = ii.gaps, 
           gaps_row =  gaps.row, 
           filename = paste0(figureDir, 'heatmap_regenerationPeaks_scaled.pdf'), 
           width = 6, height = 12)
  
  yy = test
  plt = pheatmap(yy, 
                 annotation_row = my_gene_col, 
                 annotation_col = df, show_rownames = FALSE, scale = 'none', 
                 color = col, 
                 show_colnames = FALSE,
                 cluster_rows = FALSE, cluster_cols = FALSE,  
                 #clustering_method = 'complete', cutree_rows = nb_clusters, 
                 annotation_colors = annot_colors, 
                 #clustering_callback = callback,
                 gaps_col = ii.gaps, 
                 gaps_row =  gaps.row
                 #filename = paste0(figureDir, 'heatmap_regenerationPeaks_scaled.pdf'), 
                 #width = 6, height = 12
                 )
   
  save(yy, plt, res, file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters.Rdata'))
  
  
}


## save data for DPGP clustering
## unfornately the result is that there are too many clusters
Save.Matrix.for.DPGP = FALSE
if(Save.Matrix.for.DPGP){
  
  rownames(yy) = gsub(':', '_', rownames(yy))
  #yy = rbind(c(1:ncol(yy)), yy)
  colnames(yy) = c(1:ncol(yy))
  write.table(yy[c(1:1000), ], 
              file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/DPGP_clustering/', 
                            'atac_peakSignals_1000.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  write.table(yy[c(1:2000), ], 
              file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/DPGP_clustering/', 
                            'atac_peakSignals_2000.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  write.table(yy[c(1:5000), ], 
              file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/DPGP_clustering/', 
                            'atac_peakSignals_5000.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  write.table(yy[c(1:10000), ], 
              file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/DPGP_clustering/', 
                            'atac_peakSignals_10000.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(yy, 
              file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/DPGP_clustering/', 
                            'atac_peakSignals_all.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
  
}

##########################################
# save peak list for deeptool heatmap  
##########################################
Save_peak.list_for_deeptools = FALSE
if(Save_peak.list_for_deeptools){
  # import atac-seq peak
  # res = readRDS(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration.rds')) # stat of dynamic peaks
  # z-score of data and heatmap (variable: plt, yy and res)
  load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters.Rdata')) 
  #design_atac = design
  res_atac = res
  yy_atac = yy
  
  Test_cluster_order = FALSE
  if(Test_cluster_order){
    test = yy
    test = test[plt$tree_row$order, ]
    pheatmap(test, 
             nnotation_row = my_gene_col, 
             annotation_col = df, show_rownames = FALSE, scale = 'none', 
             color = col, 
             show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE,  
             #clustering_method = 'complete', cutree_rows = nb_clusters, 
             annotation_colors = annot_colors, 
             #clustering_callback = callback,
             gaps_col = ii.gaps)
    
    ## save group bed of each cluster for deeptools heatmap
    load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters.Rdata')) 
    
    res = res[match(rownames(yy), rownames(res)), ]
    mcs = unique(res$clusters)
    for(n in 1:length(mcs)){
      jj = which(res$clusters == mcs[n])
      pp_atac = data.frame(t(sapply(rownames(yy)[jj], function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
      pp_atac$name = rownames(yy)[jj]
      pp_atac$score = 0
      pp_atac$strand = '*'
      
      peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                            'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_groups')
      write.table(pp_atac, file = paste0(peakGroupDir, '/peak_group_', n, '.bed'), 
                  sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    }
    
    ## save bed files of promoters, enhancers of up- and down-regulated atac-seq peaks
    res = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds'))
    res = res[which(!is.na(res$clusters)), ]
    
    pp = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    pp.annots = as.data.frame(pp.annots)
    
    res$atacseq = 'up'
    res$atacseq[which(res$clusters == 'mc1'| res$clusters == 'mc6')] = 'down'
    
    res$annots = NA
    res$annots[grep('Promoter', pp.annots$annotation)] = 'promoter'
    res$annots[grep('Intergenic|Intron', pp.annots$annotation)] = 'enhancer'
    
    pp_atac = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    pp_atac$name = rownames(res)
    pp_atac$score = 0
    pp_atac$strand = '*'
    
    peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                          'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_promoter_enhancer')
    
    jj = which(res$atacseq == 'up' & res$annots == 'enhancer')
    write.table(pp_atac[jj, ], file = paste0(peakGroupDir, '/peak_enhancer_up.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    jj = which(res$atacseq == 'down' & res$annots == 'enhancer')
    write.table(pp_atac[jj, ], file = paste0(peakGroupDir, '/peak_enhancer_down.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    
    jj = which(res$atacseq == 'up' & res$annots == 'promoter')
    xx = pp_atac[jj, ]
    xx0 =  pp.annots[jj, ]
    xx$strand = xx0$geneStrand # keep tss strand
    kk = which(xx$strand == '1')
    xx$X2 = as.numeric(as.character(xx$X2))
    xx$X3 = as.numeric(as.character(xx$X3))
    xx$X2[kk] = xx0$geneStart[kk]; xx$X2[-kk] = xx0$geneStart[-kk] - 1 
    xx$X3[kk] = xx$X2[kk] + 1; xx$X3[-kk] = xx$X2[-kk]
    xx$strand[kk] = '+'
    xx$strand[-kk] = '-'
    
    write.table(xx, file = paste0(peakGroupDir, '/peak_tss_up.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    jj = which(res$atacseq == 'down' & res$annots == 'promoter')
    xx = pp_atac[jj, ]
    xx0 =  pp.annots[jj, ]
    xx$strand = xx0$geneStrand # keep tss strand
    kk = which(xx$strand == '1')
    xx$X2 = as.numeric(as.character(xx$X2))
    xx$X3 = as.numeric(as.character(xx$X3))
    xx$X2[kk] = xx0$geneStart[kk]; xx$X2[-kk] = xx0$geneStart[-kk] - 1 
    xx$X3[kk] = xx$X2[kk] + 1; xx$X3[-kk] = xx$X2[-kk]
    xx$strand[kk] = '+'
    xx$strand[-kk] = '-'
    
    
    write.table(xx, file = paste0(peakGroupDir, '/peak_tss_down.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    ## save promoters to test histone markers
    peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                          'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_tss')
    
    tss = promoters(amex, upstream=2000, downstream=2000, use.names=TRUE)
    ggs = genes(amex)
    ggs = as.data.frame(ggs)
    
    tss = ggs
    jj = which(tss$strand == '-')
    tss$start[jj] = tss$end[jj]
    tss$end = tss$start + 1
    #tss$strand = '*'
    tss$score = 0
    tss = tss[, c(1, 2, 3, 6, 7, 5)] 
    
    write.table(tss, file = paste0(peakGroupDir, '/tss_expressed.genes.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
}

##########################################
# combine profile plots of histone markers of each cluster  
##########################################
Assembly_histMarkers_profilePlots = FALSE
if(Assembly_histMarkers_profilePlots){
  library(tidyr)
  library(dplyr)
  require(ggplot2)
  library("gridExtra")
  library("cowplot")
  require(ggpubr)
  
  profileDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                      'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/heatmaps_allHist_8clusters/')
  
  # +/-4kb; binSize = 50bp;
  nb_bins = 8000/50
  pfs = read.delim(file = paste0(profileDir, 'allHist_8clusters_data4profile.tab'), header = TRUE)
  pfs = pfs[-1,]
  pfs = pfs[, !is.na(pfs[1,])] 
  
  # select +/-3kb window sizes
  xx = pfs[, c(3:ncol(pfs))]
  colnames(xx) = c(1:ncol(xx))
  xx = xx[, -c(1:20, ((ncol(xx)-19):ncol(xx)))]
  
  pfs = data.frame(pfs[, c(1:2)], xx, stringsAsFactors = FALSE)
  colnames(pfs) = c('marker', 'group',  c(1:ncol(xx)))
  pfs$group = gsub('peak_group_', '', pfs$group) 
  pfs$group =  paste0('mc', gsub('.bed', '', pfs$group))  
  pfs$sample = sapply(pfs$marker, function(x) unlist(strsplit(as.character(x), '_'))[2])
  pfs$marker = sapply(pfs$marker, function(x) unlist(strsplit(as.character(x), '_'))[1])
  pfs = pfs[, c(2, 1, ncol(pfs), 3:122)]
  
  grps = unique(pfs$group)
  markers = unique(pfs$marker)
  markers = markers[grep('H3K', markers)]
  
  plot_histM_8groups = function(xx, ylims)
  {
    plot = as_tibble(xx) %>%  gather(bins, signals,  4:123) %>% 
      mutate(bins = factor(bins, levels=c(1:120))) %>%
      mutate(sample = factor(sample, levels = c('mUA', 'BL5days', 'BL9days', 'BL13days.prox', 'BL13days.dist', 'IgG'))) %>%
      ggplot(aes(y=signals, x=bins, color = sample, group = sample)) + 
      geom_line(size = 1.2) +
      theme_classic() +
      theme(axis.text.y = element_text(angle = 90, size = 10), legend.position = c(0.2, 0.8),
            axis.text.x = element_blank(), axis.ticks = element_blank()) +
      scale_color_manual(values=c('springgreen', 'springgreen2', 'springgreen3', 'gold2','red')) +
      labs( x = '', y = 'signals (rpkm)') +
      ylim(0, ylims) 
    return(plot)
  }
  
  for(m in 1:length(markers))
  {
    # m = 3
    ylims = 5.0
    cat(markers[m], '\n')
    for(n in 1:length(grps))
    {
      #n = 1
      cat(n, ' -- ', grps[n], '\n')
      
      xx = pfs[which((pfs$marker == markers[m]) & pfs$group == grps[n]), ]
      xx[, c(4:123)] = xx[, c(4:123)] *1000/50
      
      eval(parse(text = paste0('p', n, ' = plot_histM_8groups(xx, ylims = ylims)')))
    }
    
    pdf(paste0(figureDir, "dynamics_peaks_averageProfiles_", markers[m], ".pdf"),
        width = 4, height = 20) # Open a new pdf file
    
    grid.arrange(p1, p2 + rremove('legend'), p3 + rremove('legend'), p4+rremove('legend'),
                p5 + rremove('legend'), p6 + rremove('legend'), p7 +rremove('legend'), p8+rremove('legend'), 
                 ncol = 1, nrow = 8)
    #ggsave(paste0(figureDir, "dynamics_peaks_averageProfiles_", markers[m], ".pdf"), width=4, height = 20)
    dev.off()
    
   
  }
}


##########################################
# feature distribution of all peaks and regenration dynamic peaks 
##########################################
Peak.annotation_feature.distribution = FALSE
if(Peak.annotation_feature.distribution){
  library('rtracklayer')
  library(GenomicRanges)
  library('GenomicFeatures')
  source('Functions_atac.R')
  
  # limb fibroblast expressing genes
  gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
  amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
  
  # dynamic peaks
  res = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds'))
  res = res[which(!is.na(res$clusters)), ]
  
  pp = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  
  pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  #pp.annots = as.data.frame(pp.annots)
  
  # pdfname = paste0(figureDir, "feature_distribution_positionalPeaks.pdf")
  # pdf(pdfname, width = 6, height = 4)
  # par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  plotPeakAnnot_piechart(pp.annots)
  
  # dev.off()
  
  # all peaks 
  peak.all = readRDS(file = paste0(RdataDir, '/ATACseq_peak_consensus_filtered_55k.rds'))
                        
  all.annots = annotatePeak(peak.all, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  plotPeakAnnot_piechart(all.annots)
  
  
  ## compare the features distributions and make test: promoter is more stable in popultion
  # pdfname = paste0(saveDir, "/feature_distribution_regenerationPeaks.pdf")
  # pdf(pdfname, width = 8, height = 6)
  # par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  # pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  #dev.off()
  
  stats = plotPeakAnnot_piechart(all.annots)
  
  colnames(stats)[ncol(stats)] = 'all'
  stats = data.frame(stats,  plotPeakAnnot_piechart(pp.annots)[, 2])
  colnames(stats)[ncol(stats)] = 'dynamic'
  
  scores = c()
  test = c()
  test2 = c()
  for(i in 1:nrow(stats))
  {
    total = length(peak.all); 
    mm = round(total * stats[i, 2]/100)
    nn = total - mm
    qq = round(pp.annots@peakNum * stats[i, ncol(stats)]/100)
    test = c(test, phyper(qq-1, mm, nn, pp.annots@peakNum, lower.tail = FALSE, log.p = FALSE))
    test2 = c(test2, phyper(qq+1, mm, nn, pp.annots@peakNum, lower.tail = TRUE, log.p = FALSE))
  }
  
  scores = cbind(scores, test, test2)
  colnames(scores) = c('enrich.pval', 'depelet.pval')
  rownames(scores) = rownames(stats)
  
  stats = data.frame(stats, scores, stringsAsFactors = FALSE)
  
  library(tidyr)
  library(dplyr)
  require(ggplot2)
  
  features = sapply(stats$Feature, function(x){x = unlist(strsplit(as.character(x), ' ')); x = x[-length(x)]; paste0(x, collapse = ' ')})
  stats$Feature = features
  
  as_tibble(stats) %>%  gather(group, freq,  2:3) %>% 
    mutate(Feature = factor(Feature, levels=stats$Feature)) %>%
    ggplot(aes(fill=group, y=freq, x=Feature)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 10)) +
    scale_fill_manual(values=c('#999999','#E69F00')) +
    labs( x = '', y = 'percentage (%)') 
  
  ggsave(paste0(figureDir, "feature_distribution_regenerationPeaks.pdf"),  width = 8, height = 6)
  
  # for(m in 1:nb_clusters)
  # {
  #   # m = 1
  #   cat('cluster - ',  m, '\n')
  #   # pp.annots = annotatePeak(pp[which(xx$cluster == m)], TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  #   
  #   # pdfname = paste0(saveDir, "/Fig1E_positional_peak_feature_distribution_cluster_", m, ".pdf")
  #   # pdf(pdfname, width = 8, height = 6)
  #   # par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  #   
  #   stats = data.frame(stats,  plotPeakAnnot_piechart(pp.annots)[, 2])
  #   colnames(stats)[ncol(stats)] = paste0('cluster', m, '_', pp.annots@peakNum)
  #   
  #   test = c()
  #   test2 = c()
  #   for(i in 1:nrow(stats))
  #   {
  #     total = length(pp); 
  #     mm = round(total * stats[i, 2]/100)
  #     nn = total - mm
  #     qq = round(pp.annots@peakNum * stats[i, ncol(stats)]/100)
  #     test = c(test, phyper(qq-1, mm, nn, pp.annots@peakNum, lower.tail = FALSE, log.p = FALSE))
  #     test2 = c(test2, phyper(qq+1, mm, nn, pp.annots@peakNum, lower.tail = TRUE, log.p = FALSE))
  #   }
  #   scores = cbind(scores, test, test2)
  #   colnames(scores)[c(ncol(scores)-1, ncol(scores))] = c(paste0('enrich.pval.cluster', m), paste0('depelet.pval.cluster', m))
  #   
  #   dev.off()
  # }
  # rownames(scores) = rownames(scores)
  # 
  # stats = data.frame(stats, scores, stringsAsFactors = FALSE)
  # 
  # write.csv(stats, file = paste0(saveDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters', 
  #                                nb_clusters, '_feature.Enrichment.depeletion.csv'), quote = FALSE, row.names = TRUE)
  
  
}

########################################################
########################################################
# Section : Regeneration genes and associated CREs
# specific questions: chromatin states of regeneration genes in mature sample
# is it bivalent ?? or not 
########################################################
########################################################

# prepare all tss to quantify the read counts within those tss 
gtf.all = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/AmexT_v47_Hox.patch.gtf'
amex.all = GenomicFeatures::makeTxDbFromGFF(file = gtf.all)
tss = GenomicFeatures::promoters(amex.all, upstream = 2000, downstream = 2000, use.names = TRUE)

tss = as.data.frame(tss)
tss$tx_name = sapply(tss$tx_name, function(x){x = unlist(strsplit(as.character(x), '[|]')); x[length(x)]})
tss$start[which(tss$start<=1)] = 1

SAF = data.frame(GeneID=tss$tx_name, 
                 Chr=tss$seqnames, 
                 Start=tss$start, 
                 End=tss$end, 
                 Strand=tss$strand, 
                 stringsAsFactors = FALSE)

write.table(SAF, file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_histM_both/',
                               'amex6_TSS_all_56670genes.saf'), 
            sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE) 


## load processed tss atac and histMarker signals
#genes = GenomicFeatures::genes(amex.all, columns="gene_id")
res = readRDS(file = paste0(RdataDir,  '/atac_histoneMarkers_normSignals_axolotlAllTSS.2kb_TSS_genelevel.rds'))
res = data.frame(res)
res$geneID = rownames(res)

res$gene = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[2]})
res$geneID  = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})


##########################################
# first try to define groups 
##########################################
res = readRDS(file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb_TSSfiltered.rds'))

aa = readRDS(file =  "../results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/TestStat_regeneration_RNAseq.rds") 

aa$gene = sapply(rownames(aa), function(x) unlist(strsplit(as.character(x), '_'))[1])
aa$geneID = sapply(rownames(aa), function(x) {test = unlist(strsplit(as.character(x), '_')); return(test[length(test)])})

aa[grep('DBP|SALL', rownames(aa)), grep('Mature_UA', colnames(aa))]

aa$expr.mUA = apply(aa[, grep('Mature_UA', colnames(aa))], 1, mean)
aa$expr.BLday5 = apply(aa[, grep('BL_UA_5days', colnames(aa))], 1, mean)
aa$expr.BLday9 = apply(aa[, grep('BL_UA_9days', colnames(aa))], 1, mean)
aa$expr.BLday13 = apply(aa[, grep('BL_UA_13days_proximal', colnames(aa))], 1, mean)
aa$expr.BLday13.distal = apply(aa[, grep('BL_UA_13days_distal', colnames(aa))], 1, mean)
aa = aa[, c(5:11, 23:29)]
aa = aa[match(res$geneID, aa$geneID), ]
res = data.frame(res, aa, stringsAsFactors = FALSE)

res$x1 = res$UA_K4me3
res$x2 = res$UA_K27me3
res$ratio = res$x1 - res$x2

examples.sel = c()
examples.sel = unique(c(examples.sel, grep('HOXA13|HOXA11|HOXA9|HOXD13|HOXD11|HOXD9|MEIS|SALL|DBP', res$gene)))


ggplot(data=res, aes(x=x1, y=x2, label = gene)) +
  geom_point(size = 0.25) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=4, col='darkgray') +
  geom_hline(yintercept=4, col="darkgray") +
  labs(x = "UA_H3K4me3", y= 'UA_H3K27me3')


ggplot(data=res, aes(x=ratio, y=expr, label = gene)) +
  geom_point(size = 0.25) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=4, col='darkgray') +
  geom_hline(yintercept=4, col="darkgray") +
  labs(x = "UA_H3K4me3/UA_H3K27me3", y= 'Expr')

load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
hkgs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
nonexp = controls.tissue$geneIDs[which(controls.tissue$tissues != 'housekeeping')]

res$groups = 'limb'
res$groups[!is.na(match(res$geneID, hkgs))] = 'house_keep'
res$groups[!is.na(match(res$geneID, nonexp))] = 'other_tissues'

xx = res[order(-res$expr), ]
head(xx[which(xx$groups != 'house_keep'), ])

ggplot(data=res, aes(x=expr.mUA, y=ratio, label = gene, color = groups)) +
  geom_point(size = 0.25) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=4, col='darkgray') +
  geom_hline(yintercept=4, col="darkgray") +
  labs(y = "UA_H3K4me3/UA_H3K27me3", x= 'Expr')

res[,c(1, 4, 7:ncol(res))] %>% 
  pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
  ggplot(aes(x = factor(groups, levels = c('other_tissues', 'house_keep', 'limb')), y=value, fill=markers)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 1.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14)) +
  labs(x = "", y= 'normalized data (log2)')

##########################################
# redefine different groups with RNA-seq data 
##########################################
Redefine.gene.groups.with.RNAseq = TRUE
if(Redefine.gene.groups.with.RNAseq){
  
  table(res$groups)
  kk = which(res$groups == 'limb')
  
  ss = apply(res[kk, grep('expr', colnames(res))], 1, max)
  
  select = kk[which(ss< -1 | is.na(ss))]
  res$groups[select] = 'lowlyExpr_limb'
  table(res$groups)
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  #jj = which(res$x1>4 & res$x2 >4)
  #xx = res[jj, ]
  
  kk = which(res$groups == 'limb')
  
  diffs = res$expr.BLday5[kk] - res$expr.mUA[kk]
  select = kk[which(res$expr.mUA[kk] < 0 & diffs >2)]
  
  res$groups[select] = 'upregulated.d5'
  table(res$groups)
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) >2)]
  res$groups[select] = 'downregulated.d5'
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 'upregulated.d5', 'downregulated.d5', 'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  
  kk = which(res$groups == 'limb')
  
  select = kk[which(res$expr.mUA[kk] < 0 & (res$expr.BLday5[kk] - res$expr.mUA[kk]) < 2 &  (res$expr.BLday9[kk] - res$expr.mUA[kk]) > 2)]
  res$groups[select] = 'upregulated.d9'
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) < 2 & (res$expr.mUA[kk] - res$expr.BLday9[kk]) >2)]
  res$groups[select] = 'downregulated.d9'
  
  
  kk = which(res$groups == 'limb')
  
  select = kk[which(res$expr.mUA[kk] < 0 & (res$expr.BLday5[kk] - res$expr.mUA[kk]) < 2 &  (res$expr.BLday9[kk] - res$expr.mUA[kk]) < 2 &
                      (res$expr.BLday13[kk] - res$expr.mUA[kk]) >2)]
  res$groups[select] = 'upregulated.d13'
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) < 2 & (res$expr.mUA[kk] - res$expr.BLday9[kk]) < 2 &
                      (res$expr.mUA[kk] - res$expr.BLday9[kk]) > 2) ]
  res$groups[select] = 'downregulated.d13'
  
  
  
  
  table(res$groups)
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 
                                             'upregulated.d5', 'upregulated.d9','upregulated.d13', 
                                             'downregulated.d5',  'downregulated.d9',
                                             'downregulated.d13',
                                             'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  
  ggsave(paste0(figureDir, "histMarker_H3K27me3_H3K4me3_up.downregulatedGenes.UA.BL.days.pdf"), width=12, height = 8)
  
  fdr.cutoff = 0.05
  select = which(aa$padj < fdr.cutoff & aa$log2FC >1 & aa$log2FC.mUA.vs.others >0)
  res$groups[!is.na(match(res$geneID, aa$geneID[select]))] = 'mature_highlyExp'
  table(res$groups)
  
  select = which(aa$padj < fdr.cutoff & aa$log2FC >1 & aa$log2FC.mUA.vs.others < 0)
  res$groups[!is.na(match(res$geneID, aa$geneID[select]))] = 'regeneration'
  table(res$groups)
  
  select = which((aa$padj >= fdr.cutoff | aa$log2FC <=1) & log2(aa$baseMean) <4) 
  res$groups[!is.na(match(res$geneID, aa$geneID[select])) & res$groups == 'limb_static'] = 'limb_lowlyExp'
  
  
}

ggplot(data=res, aes(x=x1, y=x2, label = gene, color = groups)) +
  geom_point(size = 0.4) + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
  #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
  #scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=4, col='darkgray') +
  geom_hline(yintercept=4, col="darkgray") +
  labs(x = "UA_H3K4me3", y= 'UA_H3K27me3')


#xx = melt(res[], id.vars = c('UA_K4me3', 'UA_K27me3'), variable_name = 'markers')
res[,c(1, 4, 7:13)] %>% 
  pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
  ggplot(aes(x = factor(groups, levels = c('other_tissues', 'house_keep', 'limb_lowlyExp',
                                           'mature_highlyExp', 'regeneration', 'limb_static')), y=value, fill=markers)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 1.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14)) +
  labs(x = "", y= 'normalized data (log2)')


########################################################
########################################################
# Section : regeneration peaks, not found in mUA and embryo stages only in regeneration process 
# 
########################################################
########################################################
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



########################################################
########################################################
# Section : chromatin dynamics of position-related CREs and genes
# focus on the dynamics of Hand-specific CREs in regeneration process
# 
########################################################
########################################################

##########################################
# process the positional result and regeneration analysis and combine them 
##########################################

## import positional peaks in Figure 1
peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rdata/', 
                              'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))
table(peaks$clusters)
cluster_order = c(6, 1, 5, 3, 4, 2)
peaks$new_clusters = NA

# change the peak names 
for(n in 1:length(cluster_order))
{
  peaks$new_clusters[which(peaks$clusters == cluster_order[n])] = n 
  #peaks$new_clusters[which(peaks$clusters == '1')] = 2 
}

## import regeneration peaks 
library(qvalue)
res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks__mUA_regeneration_dev_2Batches.R10723_R7977_peakAnnot_v8.rds'))

# ignore the current annotation
res = res[, c(1:50)]

res = res[order(-res$log2FC), ]
qv = qvalue(res$pval.lrt)
res$fdr.lrt = qv$qvalues

# select the temporal dynamic peaks
fdr.cutoff = 0.01; logfc.cutoff = 1

select = which(  (res$adj.P.Val_5dpa.vs.mUA < fdr.cutoff & abs(res$logFC_5dpa.vs.mUA) > logfc.cutoff)| 
                   (res$adj.P.Val_9dpa.vs.mUA < fdr.cutoff & abs(res$logFC_9dpa.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_13dpap.vs.mUA < fdr.cutoff & abs(res$logFC_13dpap.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_13dpad.vs.mUA < fdr.cutoff & abs(res$logFC_13dpad.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_s40.vs.mUA < fdr.cutoff & abs(res$logFC_s40.vs.mUA) > logfc.cutoff ) | 
                   (res$adj.P.Val_s44p.vs.mUA < fdr.cutoff & abs(res$logFC_s44p.vs.mUA) > logfc.cutoff ) |
                   (res$adj.P.Val_s44d.vs.mUA < fdr.cutoff & abs(res$logFC_s44d.vs.mUA) > logfc.cutoff ))
cat(length(select), 'DE peaks found !\n')

res$dynamic.peaks = 0
res$dynamic.peaks[select] = 1

res = res[order(-res$log2FC), ]

## add dpgp cluster result
xx = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds'))
res$dpgp.clusters = NA
res$dpgp.clusters[match(rownames(xx), rownames(res))] = xx$clusters

## transfer the positional peak info into the regeneration analysis
rownames(res) = gsub('_', '-', rownames(res))
mm = match(rownames(peaks), rownames(res))

res$positional.peaks = 0
res$positional.clusters = NA
res$positional.peaks[mm] = 1
res$positional.clusters[mm] = peaks$new_clusters

## add histone marker analysis
RdataHistM = '/Users/jiwang/workspace/imp/positional_memory/results/CT_merged_20220328/Rdata'
load(file = paste0(RdataHistM, '/combined_4histMarkers_overlapped55kATACseq_DE_regeneration_v2.Rdata')) # variable (keep and DE.locus)

mm = match(rownames(keep), rownames(DE.locus))
DE.locus = DE.locus[mm, ]

colnames(DE.locus) = paste0('DE_', colnames(DE.locus))

mm = match(rownames(res), keep$atac_peak_coord_H3K4me3)
res = data.frame(res, DE.locus[mm, ], stringsAsFactors = FALSE)

saveRDS(res, file = paste0(RdataDir, '/allPeaks_regeneraton_positional_histM_analysisRes.rds'))
saveRDS(keep, file = paste0(RdataDir, '/allPeaks_regeneraton_positional_histM_analysisRes_histMdata.rds'))

##########################################
# explore the atac-seq peak dyanmic  
##########################################
res = readRDS(file = paste0(RdataDir, '/allPeaks_regeneraton_positional_histM_analysisRes.rds'))

res = res[which(res$positional.peaks == 1 & res$positional.clusters != 1 & res$positional.clusters !=2), ]

ss = apply(res[, grep('DE_H3K', colnames(res))], 1, sum)
length(which(ss>0))

res$DE_histM = 0
res$DE_histM[which(ss>0)] = 1

table(res$dynamic.peaks)
table(res$DE_histM)

table(res$dynamic.peaks, res$DE_histM)
table()



