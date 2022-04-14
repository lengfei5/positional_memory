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
res$mins =  apply(sample.means[match(rownames(res), rownames(sample.means)), ], 1, min)
res$max = apply(sample.means[match(rownames(res), rownames(sample.means)), ], 1, max)
# xx = xx[which(xx$min < 4.5), ]

#length(which(res$min<3))
length(which(res$max > 3))
res = res[which(res$max > 3), ]
cat(nrow(res), 'DE peaks found with at least condition above background !\n')

# saveRDS(res, file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration.rds'))

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
sample.means = sample.means[match(rownames(keep), rownames(sample.means)), ]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

library(dendextend)
yy <- t(apply(sample.means, 1, cal_z_score))

df <- data.frame(conds)
rownames(df) = colnames(yy)
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

##########################################
# combine atac-seq peaks and histone marker to make heatmap 
##########################################
Assembly_histMarkers_togetherWith_ATACseq = FALSE
if(Assembly_histMarkers_togetherWith_ATACseq){
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  require(pheatmap)
  require(RColorBrewer)
  
  # import atac-seq peak
  res = readRDS(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration.rds')) # stat of dynamic peaks
  load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap.Rdata')) # z-score of data and heatmap (variable plt, yy)
  
  design_atac = design
  res_atac = res
  yy_atac = yy
  
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
  
  # import histone marker 
  RdataHistM = '/Users/jiwang/workspace/imp/positional_memory/results/CT_merged_20220328/Rdata'
  load(file = paste0(RdataHistM, '/combined_4histMarkers_overlapped55kATACseq_DE.Rdata')) # variables (keep and DE.locus)
  design = readRDS(file = paste0(RdataHistM, '/histM_CT_design_info.rds'))
  
  yy = keep[, c(53:60, 27:34, 1:8, 79:85)]
  sampleID = sapply(colnames(yy), function(x) unlist(strsplit(as.character(x), '_'))[3])
  mm = match(sampleID, design$sampleID)
  colnames(yy) = paste0(design$condition[mm], '_', design$Batch[mm], '_', design$sampleID[mm])
    
  ## select the regions overlapping with atac-seq
  pp_atac = data.frame(t(sapply(rownames(yy_atac), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
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
  
  ## peaks not mapped to positional atac-peaks
  yy0 = yy[-mapping@to, ]
  
  ss = apply(DE.locus[, c(1:3)], 1, sum)
  mm = match(names(ss[which(ss>0)]), rownames(yy0))
  length(which(!is.na(mm)))
  
  yy0 = yy0[mm[!is.na(mm)], ]
  
  #plot.pair.comparison.plot(yy1[, grep('H3K4me1_mUA_', colnames(yy1))], linear.scale = FALSE)
  #plot.pair.comparison.plot(yy1[, grep('H3K4me3_mUA_', colnames(yy1))], linear.scale = FALSE)
  #plot.pair.comparison.plot(yy1[, grep('H3K27me3_mUA_', colnames(yy1))], linear.scale = FALSE)
  saveRDS(yy1, file = paste0(RdataDir, '/peak_signals_atac_4histM_regenerationPeaks.rds'))
  saveRDS(yy0, file = paste0(RdataDir, '/peak_signals_atac_4histM_notOverlapped.regenerationPeaks.rds'))
  
  yy1 = yy1[, grep('rRep', colnames(yy1))]
  yy0 = yy0[, grep('rRep', colnames(yy0))]
  
  # not consider H3K27ac in the main figure
  yy1 = yy1[, grep('H3K27ac', colnames(yy1), invert = TRUE)]
  
  mm = match(rownames(yy1), rownames(DE.locus))
  DE.peaks = DE.locus[mm, ]
  ss = apply(DE.peaks, 1, sum)
  length(which(ss>0))
  
  source('Functions_histM.R')
  conds_histM = c('H3K4me1', 'H3K27me3', 'H3K4me3', 'H3K27ac')
  for(n in 1:length(conds_histM))
  {
    jj = grep(conds_histM[n], colnames(yy1))
    #yy1[,jj] = t(apply(yy1[,jj], 1, cal_centering))
    yy1[ ,jj] = t(apply(yy1[,jj], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
    
    jj0 = grep(conds_histM[n], colnames(yy0))
    yy0[ ,jj0] = t(apply(yy0[,jj0], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
    
  }
  
  
  ### peaks overlapped with atac-seq 
  df = as.data.frame(sapply(colnames(yy1), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'segments'
  rownames(df) = colnames(yy1)
  
  sample_colors = c('springgreen4', 'steelblue2', 'gold2')
  annot_colors = list(segments = sample_colors)
  
  gaps_col = c(6, 12)
  pheatmap(yy1, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           gaps_col = gaps_col,
           legend = TRUE,
           #annotation_colors = annot_colors,
           width = 4, height = 8, 
           filename = paste0(figureDir, '/heatmap_histoneMarker_1246.postionalPeaks.pdf'))
  
  ## peaks not overlapped with atac
  df = as.data.frame(sapply(colnames(yy0), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'segments'
  rownames(df) = colnames(yy0)
  
  sample_colors = c('springgreen4', 'steelblue2', 'gold2')
  annot_colors = list(segments = sample_colors)
  
  gaps_col = c(6, 12, 18)
  pheatmap(yy0, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           gaps_col = gaps_col,
           legend = TRUE,
           #annotation_colors = annot_colors,
           width = 4, height = 8, 
           filename = paste0(figureDir, '/heatmap_histoneMarker_DE_notoverlapped.with.atac.postionalPeaks.pdf'))
  
  
  
  ##########################################
  #  # reorder the histM according to the atac-seq peak clusters
  # cluster order : 6, 1, 5, 3, 4, 2
  # gaps.row = c(32, 32+76, 31 + 32 + 76, 292 + 31 + 32 + 76, 292 + 31 + 32 + 76 + 103)
  ##########################################
  peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rdata/', 
                                'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))
  table(peaks$clusters)
  cluster_order = c(6, 1, 5, 3, 4, 2)
  
  ## specify row gaps as atac-seq peaks
  gaps.row = c()
  for(n in 1:(length(cluster_order)-1))
  {
    if(n == 1)  {gaps.row = c(gaps.row, length(which(peaks$clusters == cluster_order[n])))
    }else{
      gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(peaks$clusters == cluster_order[n])))
    }
  }
  
  ## refine histone subclusters for each ATACseq peak cluster
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  
  for(ac in cluster_order)
  {
    # ac = 6
    kk = which(peaks$cluster == ac)
    cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
    
    hm_hclust <- hclust(dist(as.matrix(yy1[kk,])), method = "complete")
    #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
    peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
  }
  
  new_order = match(peakNm, rownames(yy1))
  
  
  ## plot the histone marker side by side with replicates
  df_histM = data.frame(segments = sapply(colnames(yy1), function(x) unlist(strsplit(as.character(x), '_'))[2])
                        #markers = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[2]), stringsAsFactors = FALSE
  )
  colnames(df_histM) = c('seg')
  rownames(df_histM) = colnames(yy1)
  sample_colors_histM = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors_histM) = c('mUA', 'mLA', 'mHand')
  #marker_colors_histM = c('blue', 'red', 'deepskyblue2', 'darkgreen')
  #names(marker_colors_histM) = histMs
  annot_colors_histM = list(seg = sample_colors_histM)
  gaps.col_histM = c(2, 4)
  
  kk = c(1:6)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = colorRampPalette((brewer.pal(n = 7, name ="BrBG")))(10)
  p1 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                #color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlGnBu")))(8), 
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col= df_histM_new,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  kk = c(7:12)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  
  cols = rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da"))
  p2 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  
  kk = c(13:18)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = rev(terrain.colors(10))
  p3 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                annotation_legend = FALSE,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  plot_list[['p3']]=p3[[4]]
  
  pdf(paste0(figureDir, "/Segemet_specific_chromatin_landscape_atac_histM.pdf"),
      width = 8, height = 10) # Open a new pdf file
  
  layout = matrix(c(1, 2, 3), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  dev.off()
  
  
  ### collect the four markers with means 
  histMs =  c('H3K4me3','H3K27me3', 'H3K4me1')
  shm = c()
  
  for(n in 1:length(histMs))
  {
    # n = 1
    cat(n, ' -- ', histMs[n], '\n')
    cpm = readRDS(file = paste0(RdataHistM, '/fpm_bc_TMM_combat_', histMs[n], '_', version.analysis.HistM, '.rds'))
    cpm = cpm[ , grep(histMs[n], colnames(cpm))]
    #design.sel = readRDS(file = paste0(RdataHistM, '/design.sels_bc_TMM_combat_DBedgeRtest_', histMs[n], '_', 
    #                                   version.analysis.HistM, '.rds'))
    
    ### select the samples and extract sample means
    conds_histM = c("mUA", "mLA", "mHand")
    sample.sels = c();  
    cc = c()
    
    sample.means = c()
    for(ii in 1:length(conds_histM)) 
    {
      kk = grep(conds_histM[ii], colnames(cpm))
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds_histM[ii], length(kk)))
      if(length(kk)>1) {
        sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
      }else{
        sample.means = cbind(sample.means, cpm[, kk])
      }
    }
    colnames(sample.means) = paste0(conds_histM, '_', histMs[n])
    
    if(n == 1){
      shm = sample.means
    }else{
      shm = cbind(shm, sample.means[match(rownames(shm), rownames(sample.means)), ])
    }
  }    
  
  ## quantile normalization for different histone markers
  library(preprocessCore)
  xx = normalize.quantiles(shm)
  colnames(xx) = colnames(shm)
  rownames(xx) = rownames(shm)
  shm = xx
  
  pp_histM = data.frame(t(sapply(rownames(shm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp_histM$strand = '*'
  pp_histM = makeGRangesFromDataFrame(pp_histM, seqnames.field=c("X1"),
                                      start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(pp, pp_histM, type = 'within')
  
  yy0 = shm[mapping@to, ]
  #rownames(yy0) = rownames(yy)
  
  # reorder the histM according to the atac-seq peak clusters
  yy0 = yy0[plt$tree_row$order, ]
  yy0 = t(apply(yy0, 1, cal_z_score))
  #test = yy[plt$tree_row$order, ]
  #rownames(yy0) = rownames(yy)
  #yy0 = cbind(yy, yy0)
  
  df_histM = data.frame(segments = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[1])
                        #markers = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[2]), stringsAsFactors = FALSE
  )
  colnames(df_histM) = c('seg')
  rownames(df_histM) = colnames(yy0)
  sample_colors_histM = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors_histM) = c('mUA', 'mLA', 'mHand')
  marker_colors_histM = c('blue', 'red', 'deepskyblue2', 'darkgreen')
  names(marker_colors_histM) = histMs
  annot_colors_histM = list(seg = sample_colors_histM)
  gaps.col_histM = seq(1, ncol(yy0), by = 1)
  #clusters <- rownames(subsetDat[heat$tree_row[["order"]],])
  clusters <- sort(cutree(plt$tree_row, k = nb_clusters))
  
  # cluster order : 6, 1, 5, 3, 4, 2
  gaps.row = c(32, 32+76, 31 + 32 + 76, 292 + 31 + 32 + 76, 292 + 31 + 32 + 76 + 103)
  
  #yy0 = yy0[plt$tree_row$order, ]
  p1 = pheatmap(test, 
                nnotation_row = my_gene_col, 
                annotation_col = df, show_rownames = FALSE, scale = 'none', 
                color = col, 
                show_colnames = FALSE,
                cluster_rows = FALSE, cluster_cols = FALSE,  
                #clustering_method = 'complete', cutree_rows = nb_clusters, 
                annotation_colors = annot_colors, 
                #clustering_callback = callback,
                legend = TRUE,
                gaps_col = gaps.col, 
                gaps_row = gaps.row)
  
  p2 = pheatmap(yy0, cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  #grid.arrange(p1, p2,  nrow = 1)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  #plot_list[['p3']]=p2[[4]]
  grid.arrange(grobs=plot_list, ncol=2)
  
  indexs = data.frame(peak = names(clusters), cluster= clusters, stringsAsFactors = FALSE)
  indexs = indexs[match(rownames(test), rownames(indexs)), ]
  c_order = unique(indexs$cluster)
  
  yy1 = yy0
  
  ### transform the histM to account for different backgrounds and scaling
  # for(m in 1:ncol(yy1))
  # {
  #   jj1 = which(yy1[ ,m] < 1.5)
  #   yy1[jj1, m] = 1.5
  #   jj2 = which(yy1[ ,m] > 4)
  #   yy1[jj2, m] = 4
  #    
  # }
  # 
  
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  
  for(ac in c_order)
  {
    # ac = 6
    kk = which(indexs$cluster == ac)
    cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
    
    hm_hclust <- hclust(dist(as.matrix(yy1[kk,])), method = "complete")
    #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
    peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
    
  }
  
  new_order = match(peakNm, rownames(yy1))
  
  p1 = pheatmap(test[new_order, ], 
                nnotation_row = my_gene_col, 
                annotation_col = df, show_rownames = FALSE, scale = 'none', 
                color = col, 
                show_colnames = FALSE,
                cluster_rows = FALSE, cluster_cols = FALSE,  
                #clustering_method = 'complete', cutree_rows = nb_clusters, 
                annotation_colors = annot_colors, 
                #clustering_callback = callback,
                legend = TRUE,
                annotation_legend = FALSE,
                gaps_col = gaps.col, 
                gaps_row = gaps.row)
  
  gaps.col_histM = c(1:2)
  df_histM_new = as.data.frame(df_histM[c(1:3),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(1:3)]
  
  p2 = pheatmap(yy1[new_order, c(1:3)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col= df_histM_new,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  df_histM_new = as.data.frame(df_histM[c(4:6),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(4:6)]
  p3 = pheatmap(yy1[new_order, c(4:6)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  df_histM_new = as.data.frame(df_histM[c(7:9),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(7:9)]
  p4 = pheatmap(yy1[new_order, c(7:9)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                annotation_legend = FALSE,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  plot_list[['p3']]=p3[[4]]
  plot_list[['p4']]=p4[[4]]
  
  pdf(paste0(figureDir, "/Segemet_specific_chromatin_landscape_atac_histM_firstTry.pdf"),
      width = 10, height = 12) # Open a new pdf file
  
  layout = matrix(c(1, 1, 2, 3, 4), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  dev.off()
  
  
}



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

########################################################
########################################################
# Section : Consider the grouping peaks or chromatin states: 
# inactive promoters, active promoter, house-keeping genes, mature-specific genes,
# regeneration genes, positional genes
########################################################
########################################################
fpm = readRDS(file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb.rds'))
res = data.frame(log2(fpm + 2^0), stringsAsFactors = FALSE)
res$gene = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[2]})
res$geneID  = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})

##########################################
# clean TSS, for each gene, only keep the TSS with sigals if there are mulitple ones
##########################################
Clean.TSS = TRUE

if(Clean.TSS){
  ggs = unique(res$geneID)
  sels = c()
  for(n in 1:length(ggs))
  {
    # n = 1
    jj = which(res$geneID == ggs[n])  
    if(length(jj)>0){
      ss = apply(as.matrix(res[jj, c(1:6)]), 1, sum)
      sels = c(sels, jj[which.max(ss)])
    } 
    
  }
}

res = res[sels, ]

saveRDS(res, file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb_TSSfiltered.rds'))

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



