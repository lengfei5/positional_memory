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
source('Functions_histM.R')

version.analysis = 'Rxxxx_R10723_R11637_R12810_atac'
#peakDir = "Peaks/macs2_broad"

resDir = paste0("../results/", version.analysis)

RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/'
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

figureDir = '~/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0('~/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/SupTables/')

saveTables = FALSE

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

# limb fibroblast expressing genes
gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
#amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)

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
  # saveRDS(pp, file = paste0(RdataDir, '/ATACseq_peak_consensus_filtered_64k.rds'))
  
}


##########################################
# test PCA for three samples 
##########################################
PCA_for_Wouter.Monika = FALSE
if(PCA_for_Wouter.Monika){
  sels = which(design$SampleID == '137330'| design$SampleID == '136164'| design$SampleID == '136166')
  
  #design = design[sels, ]
  ddx = dds[, sels]
  ddx$conds = droplevels(ddx$conds)
  ss = rowSums(counts(ddx))
  #dds$batch = droplevels(dds$batch)
  condition <- factor(c("A","A","A"))
  ddx = DESeqDataSetFromMatrix(counts(ddx), DataFrame(condition), ~ 1 )
  
  ddx = estimateSizeFactors(ddx)
  
  vsd <- varianceStabilizingTransformation(ddx, blind = FALSE)
  
  #pca=plotPCA(vsd, intgroup = c('conds'), returnData = FALSE)
  #print(pca)
  
  # modify the PCA plot (https://github.com/mikelove/DESeq2/blob/devel/R/plots.R)
  # calculate the variance for each gene
  rv <- rowVars(assay(vsd))
  
  # select the ntop genes by variance
  ntop = 5000
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(vsd)[select,]), rank. = 10)
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  # assembly the data for the plot
  pcsToUse = 1:3
  pcs <- paste0("PC", pcsToUse)
  d <- data.frame(V1=pca$x[,pcsToUse[1]],
                  V2=pca$x[,pcsToUse[2]],
                  V3=pca$x[,pcsToUse[2]],
                  sample=colnames(vsd))
  colnames(d)[1:3] <- pcs
  
  pca2save = d
  #pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition'), returnData = TRUE, ntop = 5000))
  
  ggplot(data=pca2save, aes(PC1, PC2, label = sample, color= sample))  + 
    geom_point(size=4) + 
    geom_text(hjust = 1, nudge_y = 1, size=4) +
    theme_classic() + 
    theme(legend.text = element_text(size=12),
          legend.title = element_text(size = 14),
          #legend.position=c(0.8, 0.2),
          plot.margin = margin(),
          axis.text.x = element_text(angle = 0, size = 14), 
          axis.text.y = element_text(angle = 0, size = 14), 
          axis.title =  element_text(size = 14)
          #legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    ) 
  
  
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
  
  # save(ddx, design.sels, file = paste0(RdataDir, '/regeneration_samples_beforeBatchCorrection.Rdata'))
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
  
  write.csv2(design, file = paste0(tableDir, 'RegDev_atac_sampelInfos.csv'), row.names = FALSE)
  
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
    # one enhancer defined by Akane: chr2p:880963799-880966329
    gr <- GRanges(
      seqnames = Rle(c("chr2p", "chr7p"), c(1, 1)),
      ranges = IRanges(start = c(880963799, 216589802), end = c(880966329, 216590446), names = c('axe.R', 'Runx1.promoter')),
      strand = Rle(strand(c("*")), c(2)),
      score = rep(1, 2)
    )
    gr 
    
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
    
    overlapsAny(gr, pp)
    pp[overlapsAny(pp, gr)]
    
    loci =   pp[overlapsAny(pp, gr)]
    ## peak name overlapping with Akane's axe.R is chr2p:880964725_880965275
    
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

#mm = match(colnames(keep), design$samples)
#design = design[mm, ]


sample.means = sample.means[match(rownames(keep), rownames(sample.means)), ]

saveRDS(sample.means, paste0(RdataDir, '/sampleMean_regeneration_peaks_beforeScaling_GPDPclustering.rds'))

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

yy <- t(apply(sample.means, 1, cal_z_score))

if(saveTables)
{
  sample.means = readRDS(paste0(RdataDir, '/sampleMean_regeneration_peaks_beforeScaling_GPDPclustering.rds'))

  res = readRDS(file = paste0(RdataDir, 
                              '/renegeration_dynamicPeaks_GPDPclustering.merged.extended_manualBCcorrected.rds'))
  res = res[which(!is.na(res$clusters)), ]
  #yy = yy[match(rownames(res), rownames(yy)), ]
  
  mm = match(rownames(res), rownames(sample.means))
  xx = sample.means[mm, ]
  
  res = res[, -c(21:27,51:52)]
  
  xx = data.frame(xx, stringsAsFactors = FALSE)
  xx$clusters = res$clusters[match(rownames(xx), rownames(res))]
  
  write.csv2(xx, 
             file = paste0(tableDir, 
                           'Dynamic_atacseqPeaks_clustering.csv'), row.names = TRUE)
  
  write.table(res, 
              file = paste0(tableDir, 
                                 'Dynamic_atacseqPeaks_clustering.txt'), 
              row.names = TRUE, col.names = TRUE)
  
  
}

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
  #res = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds'))
  res = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended_manualBCcorrected.rds'))
  res = res[which(!is.na(res$clusters)), ]
  yy = yy[match(rownames(res), rownames(yy)), ]
  
  Annotate.clustered.peaks = FALSE
  if(Annotate.clustered.peaks){
    pp = data.frame(t(sapply(rownames(yy), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
    pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    
    pp.annots = as.data.frame(pp.annots)
    rownames(pp.annots) = rownames(yy)
    res = data.frame(res, pp.annots, stringsAsFactors = FALSE)
    
    kk1 = grep('HAND2|FGF8|SHH|FGF10|GREM1', res$transcriptId)
    kk2 = grep('Distal', res$annotation, invert = TRUE)
    kk1 = intersect(kk1, kk2)
    
    res[kk1, c(53, 59, 66:67)]
    
    
  }
  #res[match(names(loci), rownames(res)),]
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
           filename = paste0(figureDir, 'heatmap_regenerationPeaks_scaled_bgcorrected.pdf'), 
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
   
  save(yy, plt, res, file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters_bgcorrected.Rdata'))
  
  
}

## save data for DPGP clustering
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
# characterize the histone markers dynamics around the dynamic atac-seq peaks 
##########################################
PLOT.global.dynamic.parameters.for.histM = FALSE
if(PLOT.global.dynamic.parameters.for.histM){
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  require(pheatmap)
  require(RColorBrewer)
  library(khroma)
  # import atac-seq peak
  # z-score of data and heatmap (variable: plt, yy and res)
  # # saved variables: yy (data of atac-seq peaks), res (DE test result) and plt (saved the heatmap of dynamic atac-seq peaks)
  #load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters.Rdata')) 
  load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters_bgcorrected.Rdata'))
  #design_atac = design
  res_atac = res
  yy_atac = yy
  res_atac = res_atac[match(rownames(yy_atac), rownames(res_atac)), ]
  rm(yy)
  rm(res)
  
  # import histone marker 
  RdataHistM = '~/workspace/imp/positional_memory/results/CT_merged_20220328/Rdata'
  load(file = paste0(RdataHistM, '/combined_4histMarkers_overlapped55kATACseq_DE_regeneration_v2.Rdata')) # variable (keep and DE.locus)
  design = readRDS(file = paste0(RdataHistM, '/histM_CT_design_info.rds'))
  
  yy = keep[, grep('^H3K4me3|^H3K27me3|^H3K4me1|^H3K27ac', colnames(keep))]
  
  sampleID = sapply(colnames(yy), function(x) unlist(strsplit(as.character(x), '_'))[3])
  mm = match(sampleID, design$sampleID)
  colnames(yy) = paste0(design$condition[mm], '_', design$Batch[mm], '_', design$sampleID[mm])
  
  source('Functions_histM.R')
  conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
  conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')
  cc = paste(rep(conds_histM, each = length(conds)), conds, sep = "_")
  yy = cal_sample_means(yy, conds = cc)
  
  ## select the regions overlapping with atac-seq
  pp_atac = data.frame(t(sapply(rownames(yy_atac), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  pp_atac$strand = '*'
  pp_atac = makeGRangesFromDataFrame(pp_atac, seqnames.field=c("X1"),
                                     start.field="X2", end.field="X3", strand.field="strand")
  
  pp_histM = data.frame(t(sapply(rownames(keep), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp_histM$strand = '*'
  pp_histM = makeGRangesFromDataFrame(pp_histM, seqnames.field=c("X1"),
                                      start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(pp_atac, pp_histM, ignore.strand=TRUE, minoverlap = 100L)
  mapping = data.frame(mapping)
  mapping = mapping[match(unique(mapping$queryHits), mapping$queryHits), ]
  
  yy1 = yy[mapping$subjectHits, ]
  keep1 = keep[mapping$subjectHits, ]
  hde = DE.locus[mapping$subjectHits, ]
  
  yy_atac = yy_atac[mapping$queryHits,]
  res_atac = res_atac[mapping$queryHits, ]
  
  if(saveTables){
    xx = yy_atac
    colnames(xx) = paste0('atac_', colnames(xx))
    xx = data.frame(xx, yy1, stringsAsFactors = FALSE)
    xx$peak_coordinates = rownames(xx)
    xx$clusters = res_atac$clusters
    xx$clusters = gsub('mc', 'c', xx$clusters)
    xx = xx[, c(30, 29, 1:28)]
    write.csv2(xx, 
               file = paste0(tableDir, 
                             'Dynamic_atacseqPeaks_clustering_and_histM.csv'), row.names = TRUE)
     
  }
  
  ## specify row gaps as atac-seq peaks
  cluster_order = paste0('mc', c(1:8))
  gaps.row = c()
  gaps.row = c()
  for(n in 1:(length(cluster_order)-1))
  {
    if(n == 1)  {gaps.row = c(gaps.row, length(which(res_atac$clusters == cluster_order[n])))
    }else{
      gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(res_atac$clusters == cluster_order[n])))
    }
  }
  
  #conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
  conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1')
  conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')
  
  ## heatmap of test significance for 3 histone marks without H3K27ac
  pheatmap(hde[, c(1:3)], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
           color = c('darkgray', 'red'), 
           gaps_row = gaps.row,
           legend = FALSE,
           gaps_col = c(1,2),
           filename = paste0(figureDir, '/regeneration_histM_dynamics_for_dynamicATACpeaks_bgcor.pdf'), 
           width = 3, height = 12)
  
  ## heatmaps of log2FC, time point vs mUA for each histone marks 
  for(n in 1:length(conds_histM))
  {
    # n = 1
    ii.test = intersect(grep('logFC_', colnames(keep1)), grep(conds_histM[n], colnames(keep1)))
    test = keep1[, ii.test]
    fc_sels = c('5dpa.vs.mUA', '9dpa.vs.mUA', '13dpap.vs.mUA', '13dpad')
    
    jj.test = c()
    for(fc in fc_sels) jj.test = c(jj.test, grep(fc, colnames(test)))
    test = test[, jj.test]
    colnames(test) = fc_sels
    
    range <- 3
    test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
    
    nb_breaks = 7
    sunset <- colour("sunset")
    PRGn <- colour("PRGn")
    #highcontrast <- colour("high contrast")
    if(n == 1) cols = (sunset(nb_breaks))
    if(n == 2) cols = rev(PRGn(nb_breaks-1))
    if(n == 3)   cols = colorRampPalette(rev((brewer.pal(n = 8, name ="BrBG"))))(6)
    #cols =  colorRampPalette(colour("high contrast"))(nb_breaks)
    #cols = rev(terrain.colors(10))
    
    pheatmap(test, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
             #color = c('darkgray', 'blue'), 
             #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
             color = cols, 
             breaks = seq(-range, range, length.out = nb_breaks), 
             gaps_row = gaps.row,
             filename = paste0(figureDir, '/regeneration_histM_dynamics_log2fc.vs.mUA_dynamicATACpeaks_', 
                               conds_histM[n], '_bgcor.pdf'), 
             width = 3, height = 12)
  }
  
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
# combine profile plots of histone markers of each cluster from Deeptools output  
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
      mutate(sample = factor(sample, levels = c('mUA', 'BL5days', 'BL9days', 'BL13days.prox', 'BL13days.dist'))) %>%
      ggplot(aes(y=signals, x=bins, color = sample, group = sample)) + 
      geom_line(aes(linetype=sample, color=sample), size = 0.7) +
      theme_classic() +
      theme(axis.text.y = element_text(angle = 90, size = 10), legend.position = c(0.2, 0.8),
            axis.text.x = element_blank(), axis.ticks = element_blank()) +
      scale_color_manual(values=c('springgreen', 'gold2', 'blue', 'red','red')) +
      scale_linetype_manual(values=c("solid", "longdash", "longdash", "solid", "longdash")) + 
      labs( x = '', y = 'signals (rpkm)') +
      ylim(0, ylims) 
    
    # scale_color_manual(values=c('springgreen', 'gold2', 'blue', 'red','red')) +
    #  scale_linetype_manual(values=c('solid', 'longdashed', 'dotted', "solid", "dotted")) + 
    return(plot)
  }
  
  
  for(m in 1:length(markers))
  {
    # m = 2
    ylims = max(pfs[which(pfs$marker == markers[m]), c(4:123)])*1000/50
    cat(markers[m], '\n')
    for(n in 1:length(grps))
    {
      # n = 1
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
  
}

########################################################
########################################################
# Section :
# Regeneration genes and associated CREs
# specific questions: chromatin states of regeneration genes in mature sample
# Bivalency ?? or not 
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

##########################################
# start with the dynamic genes from RNA-seq data, together with expressed genes and some controls
# the collect all chromatin features, atac, histone markers
# bivalent chromatin/promoters
# 
## load processed tss atac and histMarker signals
# Gene groups defined by regeneration RNAseq data, house keeping, 
# limb fibrobalst expressing gene, non-expressing control genes
##########################################
source('Functions_histM.R')
#load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration.embryo_counts_expressedGenes_controls.Rdata'))
tss = readRDS(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration.embryo_counts_expressedGenes_controls.rds'))
tss = tss[,61:72]
tss = tss[, c(1, 6, 2:4, 12, 5, 7:11)]
# limb fibroblast expressing genes
#expressed = readRDS(file = '../data/expressedGenes_list_limb_fibroblast_using_smartseq2.mature.regeneration_pooledscRNAseq.dev.rds')
#expressed = as.character(expressed$geneID[which(expressed$expressed == 1)])
#cat(length(expressed), ' gene expressed to consider \n')

## add rna analysis results
rna = readRDS(file = paste0("../results/RNAseq_data_used/Rdata/", 
                                     'regeneration_dynamicGeneClusters_allGenes.rds'))
rna$geneID = get_geneID(rownames(rna))

tss[, c(8:12)] = rna[match(tss$geneID, rna$geneID), c(1:5)]
colnames(tss)[8:12] = paste0('smartseq2_', colnames(tss)[8:12])
annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

tss$gene = annot$gene.symbol.toUse[match(tss$geneID, annot$geneID)]

tp = data.frame(t(sapply(tss$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
tp$strand = '*'
tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")

########
## start to add atac peaks overlapping tss
## and the histM overlapping the tss
########
aa = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks_regeneration.dev.TSS_2Batches.R10723_R7977_peakAnnot_v1.rds'))
names = rownames(aa)
names = gsub('bg_', '', names)
jj = grep('^chr', names, invert = TRUE)
names[jj] = tss$coords[match(names[jj], tss$geneID)]
names.uniq = unique(names)
aa = aa[match(names.uniq, names), ]
rownames(aa) = names.uniq

mat = aa[, c(1:20)]
aa = aa[, -c(1:20)]
mat = mat[, grep('Embryo_', colnames(mat), invert = TRUE)]
mat = cal_sample_means(mat, conds = c('Mature_UA', 'BL_UA_5days', 'BL_UA_9days', 'BL_UA_13days_proximal', 'BL_UA_13days_distal'))
newcc = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
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

fdr.cutoff = 0.01; logfc.cutoff = 1
select = which(  (res$adj.P.Val_5dpa.vs.mUA < fdr.cutoff & abs(res$logFC_5dpa.vs.mUA) > logfc.cutoff)| 
                   (res$adj.P.Val_9dpa.vs.mUA < fdr.cutoff & abs(res$logFC_9dpa.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_13dpap.vs.mUA < fdr.cutoff & abs(res$logFC_13dpap.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_13dpad.vs.mUA < fdr.cutoff & abs(res$logFC_13dpad.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val_s40.vs.mUA < fdr.cutoff & abs(res$logFC_s40.vs.mUA) > logfc.cutoff ) | 
                   (res$adj.P.Val_s44p.vs.mUA < fdr.cutoff & abs(res$logFC_s44p.vs.mUA) > logfc.cutoff ) |
                   (res$adj.P.Val_s44d.vs.mUA < fdr.cutoff & abs(res$logFC_s44d.vs.mUA) > logfc.cutoff ))
cat(length(select), 'DE peaks found !\n')
res$dynamic = NA
res$dynamic[select] = 1

colnames(res) = paste0('atac_', colnames(res))

tss =  data.frame(tss, res, stringsAsFactors = FALSE)


## add histM analysis results
keep = readRDS(file = paste0('../results/CT_merged_20220328/Rdata/regeneration_combined_4histMarkers_DE_345k.rds'))

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
      ss = apply(keep[kk, grep('H3K4me3_mUA|H3K4me3_BL', colnames(keep))], 1, mean)
      jj_sels = c(jj_sels, kk[which.max(ss)])
    }
  }
}

res = keep[jj_sels,]

tss =  data.frame(tss, res, stringsAsFactors = FALSE)

saveRDS(tss, file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM.rds'))

##########################################
# analysis around tss of regeneration-response genes
# start with which chromatin features or combinations are related to regeneration-response genes
# reminder: use the updataed tss with gene annotations
##########################################
source('Functions_histM.R')
#tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM.rds'))
tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM_geneCorrection_v3.rds'))

kk = which(!is.na(tss$gene))
rownames(tss)[kk] = paste0(tss$gene[kk], '_', rownames(tss)[kk])

dev.example = c('HOXA13', 'HOXA11', 'HOXA9', 'HOXD13','HOXD11', 'HOXD9',
                'SHH', 'FGF8', 'FGF10', 'HAND2', 'BMP4', 'ALX1',
                'ALX4', 'PRRX1', 'GREM1', 'LHX2', 'LHX9', 
                'TBX2', 'TBX4', 'LMX1', 'MEIS1', 'MEIS2', 'SALL4')

# gene list from Akane
# Creb5, Bmp2, Efna1, Stmn2, Gja3, Cpa2, Cbfa2t3, Cdh3, Lmo2, Tfap2b, Jag1, Hey1, shox2, shox1,PRRX1
mature.example = c('COL1A1', 'COL4A1', 'COL4A2', 'COL6A', 'LAMA4', 'TNXB', 
                   'MATN2', 'FBN1', 'FBLN2', 'FBLN5', 'PRELP', 'ELN', 'RSPO1', 
                   'DPT', 'IGFBP3',  'TNMD', 'TWIST2', 'COL3A1', 'COL8A2', 'RARRES1', 'KLF5')

examples.sel = unique(grep(paste0(dev.example, collapse = '|'), rownames(tss)))
tss[examples.sel, c(2,7:12)]

# original code from https://github.com/const-ae/ggupset
library(ggplot2)
library(tidyverse, warn.conflicts = FALSE)
library(ggupset)
require(UpSetR)

## upregulated genes 
DEs = tss[which(tss$groups == 'DE_up'), c(2, grep('cluster|dynamic', colnames(tss)))]
DEs = DEs[which(!is.na(DEs$cluster)), ]

kk = which(!is.na(DEs$gene))
rownames(DEs)[kk] = paste0(DEs$gene[kk], '_', rownames(DEs)[kk])

DEs = DEs[,-c(1:2)]
DEs[is.na(DEs)] = 0
colnames(DEs) = gsub('_dynamic|dynamic_', '', colnames(DEs))

pdfname = paste0(figureDir, 'chromatinFeatures_relatedTo_dynamicGenes_upregulated.pdf')
pdf(pdfname, width = 8, height = 5)
par(cex = 1.0, las = 1, mgp = c(2,2,0), mar = c(3,2,2,0.2), tcl = -0.3)

upset(DEs[, c(1:4)], sets = c("H3K4me3", "atac","H3K4me1", "H3K27me3"), 
      mb.ratio = c(0.65, 0.35), order.by = "freq", keep.order =TRUE, decreasing = TRUE, show.numbers = FALSE,
      point.size = 3., line.size = 1.5, mainbar.y.label = "gene counts", sets.x.label = "gene counts", 
      text.scale = 1.6, 
      main.bar.color = 'blue')

dev.off()

examples.sel = unique(grep(paste0(dev.example, collapse = '|'), rownames(DEs)))
DEs[examples.sel, ]

## heatmap of tss with H3K27me3 changed
ss = apply(DEs, 1, sum)
jj = which(ss==1 & DEs$H3K27me3 == 1)
ids = get_geneID(rownames(DEs)[jj])

test = tss[match(ids, tss$geneID), grep('H3K27me3_logFC_', colnames(tss))]
range(test)

range <- 3.5
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
rownames(test) = rownames(DEs)[jj]

library(khroma)
nb_breaks = 7
sunset <- colour("sunset")
PRGn <- colour("PRGn")
cols = rev(PRGn(nb_breaks))

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(test, 
         #scale = 'row',
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         treeheight_row = 15,
         #color = c('darkgray', 'blue'), 
         #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
         color = cols,
         fontsize_row = 8,
         #breaks = seq(-range, range, length.out = nb_breaks), 
         #gaps_row = gaps.row,
         clustering_callback = callback,
         filename = paste0(figureDir, '/regeneration_H3K27me3_dynamics_log2fc.vs.mUA_upregulatedGenes.pdf'), 
         width = 3, height = 12)

plt = pheatmap(test, 
               #scale = 'row',
               cluster_rows = TRUE, 
               cluster_cols = FALSE, 
               show_rownames = TRUE, 
               show_colnames = FALSE,
               treeheight_row = 20,
               color = cols,
               clustering_callback = callback)

## check  the H3K4me3 changes/ levels
kk = intersect(grep('H3K4me3_', colnames(tss)), grep('rRep', colnames(tss)))
test = tss[match(ids, tss$geneID), kk]
test = cal_sample_means(test, conds = paste0('H3K4me3_', c('mUA', 'BL5days', 'BL9days', 'BL13days.prox', 'BL13days.dist')))

quantile(test, c(0, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 0.8, 0.85,  0.90,  0.95, 0.99, 1))
range <- 6
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
rownames(test) = rownames(DEs)[jj]

nb_breaks = 7
cols = (sunset(nb_breaks))

pheatmap(test[plt$tree_row$order, ], 
         #scale = 'row',
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         #treeheight_row = ,
         #color = c('darkgray', 'blue'), 
         #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
         color = cols,
         #breaks = seq(-range, range, length.out = nb_breaks), 
         #gaps_row = gaps.row,
         clustering_callback = callback,
         filename = paste0(figureDir, '/regeneration_H3K4me3_dynamics_log2fc.vs.mUA_upregulatedGenes.pdf'), 
         width = 3, height = 12)

# check their enhancer H3K4me1
enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual.rds'))
mm = match(enhancers$targets, ids)
mm = which(!is.na(mm))
mm = mm[which(enhancers$annotation_chipseeker[mm] != 'Promoter')]
targets = enhancers$targets[mm]

kk = intersect(grep('H3K4me1_', colnames(enhancers)), grep('rRep', colnames(enhancers)))
test = enhancers[mm, kk]
test = cal_sample_means(test, conds = paste0('H3K4me1_', c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')))
test = test[, c(2:5)] - test[, 1]

xx = test
test = matrix(NA, nrow = length(ids), ncol = ncol(test))
colnames(test) = colnames(xx)
for(n in 1:length(ids))
{
  # n = 1
  kk = which(targets == ids[n])
  if(length(kk) == 1) {
    test[n, ] = xx[kk,]
  }
  if(length(kk) >1) {
    test[n, ] = apply(xx[kk, ], 2, max)
  }
}

quantile(test, c(0.05, 0.1,  0.90,  0.95, 0.99), na.rm = TRUE)
range <- 3
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))


nb_breaks = 7
cols = (sunset(nb_breaks))

pheatmap(test[plt$tree_row$order, ], 
         #scale = 'row',
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         #treeheight_row = ,
         #color = c('darkgray', 'blue'), 
         #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
         color = cols,
         #breaks = seq(-range, range, length.out = nb_breaks), 
         #gaps_row = gaps.row,
         clustering_callback = callback,
         filename = paste0(figureDir, '/regeneration_H3K4me1_enhancers_log2fc.vs.mUA_upregulatedGenes.pdf'), 
         width = 3, height = 12)


## downregulated genes
DEs = tss[which(tss$groups == 'DE_down'), c(2, grep('cluster|dynamic', colnames(tss)))]
DEs = DEs[which(!is.na(DEs$cluster)), ]
kk = which(!is.na(DEs$gene))
rownames(DEs)[kk] = paste0(DEs$gene[kk], '_', rownames(DEs)[kk])

DEs = DEs[,-c(1:2)]
DEs[is.na(DEs)] = 0
colnames(DEs) = gsub('_dynamic|dynamic_', '', colnames(DEs))

pdfname = paste0(figureDir, 'chromatinFeatures_relatedTo_dynamicGenes_downregulated.pdf')
pdf(pdfname, width = 8, height = 5)
par(cex = 1.0, las = 1, mgp = c(2,2,0), mar = c(3,2,2,0.2), tcl = -0.3)

upset(DEs[, c(1:4)], sets = c("H3K4me3", "atac","H3K4me1", "H3K27me3"), 
          mb.ratio = c(0.65, 0.35), order.by = "freq", keep.order =TRUE, decreasing = TRUE, show.numbers = FALSE,
      point.size = 3., line.size = 1.5, mainbar.y.label = "gene counts", sets.x.label = "gene counts", 
      text.scale = 1.6,
      main.bar.color = 'red')

dev.off()

ss = apply(DEs, 1, sum)
jj = which(ss==1 & DEs$H3K27me3 == 1)
ids = get_geneID(rownames(DEs)[jj])

test = tss[match(ids, rownames(tss)), grep('H3K27me3_logFC_', colnames(tss))]
range <- 3
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))

library(khroma)
nb_breaks = 7
sunset <- colour("sunset")
PRGn <- colour("PRGn")
cols = rev(PRGn(nb_breaks-1))

pheatmap(test, 
         #scale = 'row',
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, show_colnames = FALSE,
         treeheight_row = 20,
         #color = c('darkgray', 'blue'), 
         #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
         color = cols,
         #breaks = seq(-range, range, length.out = nb_breaks), 
         #gaps_row = gaps.row,
         filename = paste0(figureDir, '/regeneration_H3K27me3_dynamics_log2fc.vs.mUA_downregulatedGenes.pdf'), 
         width = 3, height = 12)


########################################################
########################################################
# Section :  co-localization of H3K27me3 and H3K4me3 at the promoters
# probably together with other chromatin states
# 
########################################################
########################################################
source('Functions_histM.R')
library(ggrepel)
library(dplyr)
library(tibble)
library("cowplot")
require(gridExtra)
library(tidyr)
require(patchwork)

## bivalent analysis here
tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM_geneCorrection_v3.rds'))
tss$gene[which(rownames(tss) == 'AMEX60DD024424')] = 'MEIS3'
tss$gene[which(rownames(tss) == 'AMEX60DD028208')] = 'PROD1'
kk = which(!is.na(tss$gene))

rownames(tss)[kk] = paste0(tss$gene[kk], '_', rownames(tss)[kk])
tss$gene[-kk] = rownames(tss)[-kk]


res = tss[which(!is.na(tss$groups)), ]
res$x = apply(res[, grep('H3K4me3_mUA', colnames(res))], 1, mean)
res$y = apply(res[, grep('H3K27me3_mUA', colnames(res))], 1, mean)
#res$groups[which(res$groups != 'house_keep'& res$groups!= 'non_expr')] = NA
res = res[which(res$x<7 & res$y <7),]

dev.example = c('HOXA13', 'HOXA11', 'HOXA9', 'HOXD13','HOXD11', 'HOXD9',
                'SHH', 'FGF8', 'FGF10', 'HAND2', 'BMP4', 'ALX1',
                'ALX4', 'PRRX1', 'GREM1', 'LHX2', 'LHX9', 
                'TBX2_', 'TBX4', 'MEIS1', 'MEIS2', 'SALL4', 'MEIS3')

mature.example = c('COL1A1', 'COL4A1', 'COL4A2', 'COL6A1', 'LAMA4', 'TNXB', 
                   'MATN2', 'FBN1', 'FBLN2', 'FBLN5', 'PRELP', 'ELN', 'RSPO1', 
                   'DPT')

examples.sel = unique(grep(paste0(dev.example, collapse = '|'), rownames(res)))
examples.sel = examples.sel[which(rownames(res)[examples.sel] != 'HAND2_AMEX60DD030069')]
examples.sel = examples.sel[which(rownames(res)[examples.sel] != 'MEIS1_AMEX60DD024424')]

matures.sel = unique(grep(paste0(mature.example, collapse = '|'), rownames(res)))

##########################################
# plot the bivalent TSS in mUA
##########################################
jj = which(res$x>1 & res$y > 0)
length(jj)
#tfs.sel = unique(grep(paste0(tfs, collapse = '|'), res$gene))
#tfs.sel = tfs.sel[!is.na(match(tfs.sel, jj))]

ggplot(data=res, aes(x=x, y=y, label = gene)) +
  geom_point(size = 0.1, color = 'darkgray') + 
  #geom_point(size = 0.1 ) +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  #scale_color_manual(values=c("darkred", 'darkorange')) + 
  #geom_point(data=res[res$groups == 'reg_up', ], aes(x=H3K4me3_mUA, y=H3K27me3_mUA),  size=0.7) +
  #geom_point(data=res[res$groups == 'reg_down', ], aes(x=H3K4me3_mUA, y=H3K27me3_mUA),  size=0.7) +
  geom_point(data=res[res$groups == 'house_keep', ], aes(x=x, y=y),  size=0.3, color = 'red') +
  geom_point(data=res[res$groups == 'non_expr', ], aes(x=x, y=y),  size=0.3, color = 'darkorange') +
  geom_point(data=res[examples.sel, ], aes(x=x, y=y),  size=1.5, color = 'black') +
  #geom_text_repel(data= res[examples.sel, ], size = 4.0, color = 'blue') +
  #geom_point(data=res[matures.sel, ], aes(x=x, y=y),  size=1.5, color = 'black') +
  geom_text_repel(data= res[c(examples.sel), ], 
                  aes(x, y),
                  size = 5,
                  color = "blue",
                  #family = 'Times',
                  fontface = 'bold',
                  # Add extra padding around each text label.
                  box.padding = unit(0.3, 'lines'),
                  # Add extra padding around each data point.
                  point.padding = unit(1.6, 'lines')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
        #legend.position=c(0.2, 0.8),
        #plot.margin = margin()
        #legend.key.size = unit(1, 'cm')
        #legend.key.width= unit(1, 'cm')
  ) + 
  geom_vline(xintercept=c(1), col='black') +
  geom_hline(yintercept=c(0), col="black") +
  labs(x = "UA_H3K4me3 (log2 cpm)", y= 'UA_H3K27me3 (log2 cpm)') +
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave(paste0(figureDir, "Bivalent_TSS_mUA_scatterplot.pdf"),  width = 8, height = 6)

## gene expression mUA vs 9dpa
rna = readRDS(file = paste0("../results/RNAseq_data_used/Rdata/",
              "smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked_RNAseq_data_used_20220408.rds"))

res$x = res$smartseq2_mUA
res$y = res$smartseq2_X9dpa
ids = get_geneID(rownames(rna))
res$log2fc.9dpa.vs.mUA = rna$log2FoldChange_dpa9.vs.mUA[match(res$geneID, ids)]
res$fdr.9dpa.vs.mUA = rna$padj_dpa9.vs.mUA[match(res$geneID, ids)]
examples.up = which(res$log2fc.9dpa.vs.mUA > 1 & res$fdr.9dpa.vs.mUA <0.05)
examples.down =  which(res$log2fc.9dpa.vs.mUA < (-1) & res$fdr.9dpa.vs.mUA <0.05)

ggplot(data=res, aes(x=x, y=y, label = gene)) +
  geom_point(size = 0.1, color = 'darkgray') + 
  #geom_point(size = 0.1 ) +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  #geom_point(data=res[res$groups == 'house_keep', ], aes(x=x, y=y),  size=0.3, color = 'red') +
  #geom_point(data=res[res$groups == 'non_expr', ], aes(x=x, y=y),  size=0.3, color = 'darkorange') +
  geom_point(data=res[examples.up, ], aes(x=x, y=y),  size=0.1, color = 'red') +
  geom_point(data=res[examples.down, ], aes(x=x, y=y),  size=0.1, color = 'blue') +
  #geom_text_repel(data= res[examples.sel, ], size = 4.0, color = 'blue') +
  #geom_point(data=res[matures.sel, ], aes(x=x, y=y),  size=1.5, color = 'black') +
  geom_point(data=res[examples.sel, ], aes(x=x, y=y),  size=0.5, color = 'black') +
  geom_text_repel(data= res[c(examples.sel), ], 
                  aes(x, y),
                  size = 5,
                  color = "black",
                  #family = 'Times',
                  fontface = 'bold',
                  # Add extra padding around each text label.
                  box.padding = unit(0.2, 'lines'),
                  # Add extra padding around each data point.
                  point.padding = unit(1., 'lines')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 15), 
        axis.text.y = element_text(angle = 0, size = 15), 
        axis.title =  element_text(size = 18),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)) + 
  geom_abline(slope = 1, intercept = 0, col = 'black') +
  #geom_vline(xintercept=c(1), col='black') +
  #geom_hline(yintercept=c(0), col="black") +
  labs(x = "gene expression in mUA (log2 cpm)", y= 'gene expression in 9dpa (log2 cpm)') +
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave(paste0(figureDir, "geneExpression_comparison_mUA_BLday9.pdf"),  width = 8, height = 6)


ego = read.csv(file = paste0(tableDir, "FigS1E_GO_term_enrichmenet_for_positional_genes_Microarray_geneSymbols.csv"),
               header = TRUE)
gg1 = unique(c(unlist(strsplit(as.character(ego$geneSymbols[2]), ';')),
               unlist(strsplit(as.character(ego$geneSymbols[5]), ';'))))
gg2 = unique(c(unlist(strsplit(as.character(ego$geneSymbols[14]), ';')),
               unlist(strsplit(as.character(ego$geneSymbols[15]), ';'))))

examples.m = which(!is.na(match(res$gene, toupper(gg1))))
examples.r =  which(!is.na(match(res$gene, toupper(gg2))))

ggplot(data=res, aes(x=x, y=y, label = gene)) +
  geom_point(size = 0.1, color = 'darkgray') + 
  #geom_point(size = 0.1 ) +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  #geom_point(data=res[res$groups == 'house_keep', ], aes(x=x, y=y),  size=0.3, color = 'red') +
  #geom_point(data=res[res$groups == 'non_expr', ], aes(x=x, y=y),  size=0.3, color = 'darkorange') +
  geom_point(data=res[examples.m, ], aes(x=x, y=y),  size=0.1, color = 'red') +
  geom_point(data=res[examples.r, ], aes(x=x, y=y),  size=0.1, color = 'darkblue') +
  #geom_text_repel(data= res[examples.sel, ], size = 4.0, color = 'blue') +
  #geom_point(data=res[matures.sel, ], aes(x=x, y=y),  size=1.5, color = 'black') +
  #geom_point(data=res[examples.sel, ], aes(x=x, y=y),  size=0.5, color = 'purple1') +
  #geom_text_repel(data= res[c(examples.sel), ], 
                  # aes(x, y),
                  # size = 5,
                  # color = "purple1",
                  # #family = 'Times',
                  # fontface = 'bold',
                  # # Add extra padding around each text label.
                  # box.padding = unit(0.2, 'lines'),
                  # # Add extra padding around each data point.
                  # point.padding = unit(1., 'lines')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 15), 
        axis.text.y = element_text(angle = 0, size = 15), 
        axis.title =  element_text(size = 18),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)) + 
  geom_abline(slope = 1, intercept = 0, col = 'black') +
  #geom_vline(xintercept=c(1), col='black') +
  #geom_hline(yintercept=c(0), col="black") +
  labs(x = "gene expression in mUA (log2 cpm)", y= 'gene expression in 9dpa (log2 cpm)') +
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave(paste0(figureDir, "geneExpression_comparison_mUA_BLday9.pdf"),  width = 8, height = 6)


##########################################
# plot RNA and chromatin features for gene examples 
##########################################
PLOT.rna.chromatinFeatures_geneExample = FALSE
if(PLOT.rna.chromatinFeatures_geneExample)
{
  tss = readRDS(file = paste0(RdataDir, '/regeneration_matureSamples_tss_perGene_smartseq2_atac_histM_v4.rds'))
  
  ## examples
  tss$gene[which(rownames(tss) == 'AMEX60DD028208')] = 'PROD1'
  tss$gene[which(rownames(tss) == 'AMEX60DD024424')] = NA
  
  kk = which(!is.na(tss$gene))
  rownames(tss)[kk] = paste0(tss$gene[kk], '_', rownames(tss)[kk])
  tss$gene[-kk] = rownames(tss)[-kk]
  
  dev.genes = c('SHH', 'FGF8', 'FGF10', 'HAND2', 'BMP4', 'ALX1',
                'ALX4', 'PRRX1', 'GREM1', 'LHX2', 'LHX9', 'TBX2_', 'TBX4', 'SALL4')
  
  outDir = "/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/Gene_Examples"   
  
  source('Functions_atac.R')
  
  plot_rna_chromainFeatures_geneExamples(tss, geneList = dev.genes, outDir = outDir, incl_Mature = FALSE)
  
}

##########################################
# check the gene expression of those bivalent promoters 
##########################################
# rna = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/',
#                             'smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked_RNAseq_data_used_20220408.rds'))
# rna = rna[, grep('Mature_UA', colnames(rna))]
# rna = data.frame(rna, mUA.mean = apply(rna, 1, median))
# rna$gene = rownames(rna)
# rna$geneID = get_geneID(rna$gene)
cutoff_active = 1
cutoff_repress = 0
jj = which(res$x > cutoff_active & res$y > cutoff_repress)
length(jj)
jj1 = which(res$x > cutoff_active & res$y <= cutoff_repress)
jj2 = which(res$x < cutoff_active & res$y > cutoff_repress)

rna = res[, c(1:2, 7:8)] 
rna$groups = 'absent'
rna$groups[jj] = 'both'
rna$groups[jj1] = 'active'
rna$groups[jj2] = 'repressive'
colnames(rna)[ncol(rna)] = 'mUA'
rna$mUA[which(is.na(rna$mUA))] = -6

library(khroma)
rna = data.frame(rna)

## gene expression distribution in  boxplot
rna%>%
  mutate(groups = factor(groups, levels = c('active', 'both', 'absent', 'repressive'))) %>%
  ggplot(aes(x = groups, y=mUA, fill= groups)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 0.8) +
  #scale_fill_discrete(values=c('black', "orange", 'darkgray',  "red",   'green')) + 
  scale_fill_manual(values=c("#117733",  "blue", 'cyan',  'magenta')) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position=c(0.9, 0.8)) +
  labs(x = "", y= 'normalized gene expression (log2 cpm)')

## gene expression distribution in histogram
rna%>%
  mutate(groups = factor(groups, levels = c('active', 'both', 'absent', 'repressive'))) %>%
  ggplot(aes(x = mUA)) +
  geom_histogram(aes(color = groups, fill = groups), 
                 position = "identity", bins = 20, alpha = 0.6) +
  scale_color_manual(values = c("#117733",  "blue", 'cyan',  'magenta')) +
  scale_fill_manual(values=c("#117733",  "blue", 'cyan',  'magenta')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position=c(0.9, 0.7)) +
  labs(y = "counts", x= 'gene expression in mUA (log2 cpm)')

ggsave(paste0(figureDir, "Bivalent_TSS_mUA_geneExpression.pdf"),  width = 6, height = 4)

##########################################
# go term enrichment for bivalent promoters 
##########################################
library(enrichplot)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(stringr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(tidyr)
library(dplyr)
require(ggplot2)
library("gridExtra")
library("cowplot")
require(ggpubr)


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

grps = c('active', 'both', 'absent', 'repressive')

for(n in 1:length(grps))
{
  # n = 1
  jj = which(rna$groups == grps[n])
  gg.expressed = unique(rna$gene[jj]) # gene set 
  xx0 = rna$gene # background
  
  gg.expressed = gg.expressed[which(gg.expressed != '' & gg.expressed != 'N/A' & !is.na(gg.expressed))]
  bgs = unique(xx0[which(xx0 != '' & xx0 != 'N/A' & !is.na(xx0))])
  
  gg.expressed = firstup(tolower(gg.expressed))
  bgs = firstup(tolower(bgs))
  
  gene.df <- bitr(gg.expressed, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  bgs.df <- bitr(bgs, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Mm.eg.db)
  
  ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                   universe     = bgs.df$ENSEMBL,
                   #OrgDb         = org.Hs.eg.db,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
  
  barplot(ego, showCategory = 20)
  
  kk = grep('deoxyribonucleic', ego@result$Description)
  ego@result$Description[kk] = 'endonuclease activity'
  kk = grep('transferase activity',  ego@result$Description)
  ego@result$Description[kk] = 'transferase activity'
  
  kk = grep('binding transcription factor activity', ego@result$Description)
  ego@result$Description[kk] = 'DNA-binding transcription factor activity'
  
  eval(parse(text = paste0('p', n, '= barplot(ego, showCategory=20) + ggtitle("GO enrichment of ', grps[n], ' genes")')))
  
}

pdfname = paste0(figureDir, 'GoEnrichment_bivalent_TSS_mUA.pdf')
pdf(pdfname, width = 12, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

dev.off()


##########################################
# Check if bivalent genes keep their bivalent state 
# or transit into active or repressive states during regeneration 
##########################################
conds_histM = c('H3K4me3','H3K27me3')
conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')

cutoff_active = 1
cutoff_repress = 0

res_sel = res[which(res$groups == 'DE_down'|res$groups == 'DE_up'|res$groups == 'house_keep'|res$groups == 'non_expr'), ]

namesx = paste0('H3K4me3_', conds)
x = cal_sample_means(as.matrix(res_sel[, grep(paste0(namesx, collapse = '|'), colnames(res_sel))]), conds = namesx)
namesy = paste0('H3K27me3_', conds)
y = cal_sample_means(as.matrix(res_sel[, grep(paste0(namesy, collapse = '|'), colnames(res_sel))]), conds = namesy)

states = data.frame(matrix(NA, nrow = nrow(res_sel), ncol = length(conds)))
colnames(states) = conds
rownames(states) = rownames(res_sel)

for(n in 1:ncol(states))
{
  # n = 1
  states[ ,n] = 0
  states[which(x[,n] < cutoff_active & y[,n] >= cutoff_repress) ,n] = 1
  states[which(x[,n] >=cutoff_active & y[,n] >= cutoff_repress) ,n] = 2
  states[which(x[,n] >=cutoff_active & y[,n] < cutoff_repress) ,n] = 3
  # jj = which(res_sel$x > cutoff_active & res_sel$y > cutoff_repress)
  # length(jj)
  # jj1 = which(res_sel$x > cutoff_active & res_sel$y <= cutoff_repress)
  # jj2 = which(res_sel$x < cutoff_active & res_sel$y > cutoff_repress)
}

states = data.frame(states)
grps = res_sel$groups
states$group = grps
#grps[which(grps == 'house_keep')] = 'highlyExpr_stable'

cluster_order = c('house_keep', 'DE_up', 'DE_down', 'non_expr')
gaps.row = c()
peakNm = c()
library(dendextend)
library(ggplot2)

for(n in 1:length(cluster_order))
{
  kk = which(grps == cluster_order[n])
  cat('clsuter ', cluster_order[n], ' -- ', length(kk), ' peaks \n')
  
  hm_hclust <- hclust(dist(as.matrix(states[kk, c(1:5)])), method = "complete")
  #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
  peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
  
  if(n == 1)  {gaps.row = c(gaps.row, length(kk))
  }else{
    if(n < length(cluster_order)) {
      gaps.row = c(gaps.row,  gaps.row[n-1] + length(kk))
    }
  }
}

new_order = match(peakNm, rownames(states))
states = states[new_order, ]


library(khroma)
muted <- colour("muted")
sunset <- colour("sunset")
pale <- colour("pale")
# prepare the sample info and also the gene groups
df <- data.frame(conds)
rownames(df) = colnames(states)[1:5]
colnames(df) = 'sample'

sample_colors = c('magenta', 'darkblue', 'springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2',
                  'red')[c(1:length(conds))]
names(sample_colors) = conds

my_gene_col <- data.frame(groups =  grps)
rownames(my_gene_col) = rownames(states)
cluster_col = c('red', '#009988', 'black',  "#F99858")
names(cluster_col) = unique(grps)

annot_colors = list(
  condition = sample_colors,
  groups = cluster_col)

#sunset <- colour("sunset")
#highcontrast <- colour("high contrast")
cols = rev(c("#117733",  "blue", 'cyan',  'magenta'))
#cols =  colorRampPalette(colour("high contrast"))(nb_breaks)

pheatmap(as.matrix(states[, c(1:5)]), 
         annotation_row = my_gene_col, 
         annotation_col = df,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         scale = 'none',
         #treeheight_row = 10,
         #color = c('darkgray', 'blue'), 
         #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
         color = cols,
         annotation_colors = annot_colors,
         fontsize_row = 20,
         breaks = seq(0, 3, length.out = 5),
         gaps_row = gaps.row,
         filename = paste0(figureDir, 'chromatinStates_bivalency_regenerations.pdf'), 
         width = 4, height = 10)
######
## check  the H3K4me3 levels and H3K27me3
geneGroup = 'DE_down'

kk = intersect(grep('H3K4me3_', colnames(res_sel)), grep('rRep', colnames(res_sel)))
test = res_sel[match(rownames(states)[which(states$group == geneGroup)], rownames(res_sel)), kk]
test = cal_sample_means(test, conds = paste0('H3K4me3_', c('mUA', 'BL5days', 'BL9days', 'BL13days.prox', 'BL13days.dist')))

quantile(test, c(0, 0.05, 0.1,  0.90,  0.95, 0.99, 1), na.rm = TRUE)
range <- 6
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= cutoff_active)] = cutoff_active; x}))

nb_breaks = 9
cols = sunset(nb_breaks)
#cols = colorRamp2(c(0, range), colors = c('white', cols[nb_breaks]))
#cols = colorRampPalette(c("white", cols[nb_breaks]))(nb_breaks)
cols = cols[3:9]
pheatmap(test, 
         #scale = 'row',
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         #treeheight_row = ,
         #color = c('darkgray', 'blue'), 
         #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
         color = cols,
         #breaks = seq(-range, range, length.out = nb_breaks), 
         #gaps_row = gaps.row,
         #clustering_callback = callback,
         filename = paste0(figureDir, '/state_transition_H3K4me3_', geneGroup, '.pdf'), 
         width = 3, height = 6)

kk = intersect(grep('H3K27me3_', colnames(res_sel)), grep('rRep', colnames(res_sel)))
test = res_sel[match(rownames(states)[which(states$group == geneGroup)], rownames(res_sel)), kk]
test = cal_sample_means(test, conds = paste0('H3K27me3_', c('mUA', 'BL5days', 'BL9days', 'BL13days.prox', 'BL13days.dist')))

quantile(test, c(0, 0.05, 0.1,  0.90,  0.95, 0.99, 1), na.rm = TRUE)
range <- 5
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= cutoff_repress)] = cutoff_repress; x}))

nb_breaks = 9
PRGn <- colour("PRGn")
cols = rev(PRGn(nb_breaks))[4:nb_breaks]

pheatmap(test, 
         #scale = 'row',
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, show_colnames = FALSE,
         #treeheight_row = 20,
         #color = c('darkgray', 'blue'), 
         #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
         color = cols,
         #breaks = seq(-range, range, length.out = nb_breaks), 
         #gaps_row = gaps.row,
         filename = paste0(figureDir, '/state_transition_H3K27me3_', geneGroup, '.pdf'), 
         width = 3, height = 6)


##########################################
# Check the chromatin features that are relevant 
# to predict gene expression patterns (gene groups defined by smartseq2)
##########################################
## collect features
yy = res[, grep('_mUA_', colnames(res))]
yy = cal_sample_means(yy, conds = c('H3K4me3_mUA', 'H3K27me3_mUA', 'H3K4me1_mUA', 'H3K27ac_mUA'))
yy = data.frame(groups = res$groups, 
                gene = res$gene,
                promoters = res$coords, 
                rna_mUA = res$smartseq2_mUA, 
                atac_mUA = res$atac_mUA, 
                yy, 
                stringsAsFactors = FALSE)

# yy = add_CpG_features(yy)
# yy = add_TFmotifs_feature(yy)

## load features importance with RF
imps = readRDS(file = paste0(RdataDir, '/RF_featuresImportance.rds'))
imps = imps[c(1:10), ]

ggplot(data = imps, aes(x = scores, y = rank, label = names)) +   
  geom_point(size = 3.0, color = 'blue') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 16)) +
  #geom_text_repel(data=subset(yy, pvalue_pos > 2), size = 4)
  geom_text_repel(size = 4, box.padding = 0.4) +
  labs(x = "Importance scores (Random Forest)", y= 'feature rank')

ggsave(paste0(figureDir, "promoter_chromatinFeatures_Importance_RF_forUpDownRegulated.pdf"),  width = 6, height = 4)

## make plot for importance features
yy = readRDS(file = paste0(RdataDir, '/chromatin_promoter_features_Tfmotifs_geneGroups.rds')) 

fts = c('rna_mUA')
yy$rna_mUA[is.na(yy$rna_mUA)] = -6

yy %>% 
  pivot_longer(cols = fts, names_to = 'features') %>%
  mutate(features = factor(features, levels = fts)) %>%
  mutate(groups = factor(groups, 
                         levels = c('DE_up', 'DE_down', 'house_keep', 'highlyExpr_stable', 'lowlyExpr_stable',  'non_expr'))) %>% 
  ggplot(aes(x = features, y=value, fill= groups)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 16), 
        axis.text.y = element_text(angle = 0, size = 16), 
        axis.title =  element_text(size = 16),
        legend.text = element_text(size=10),
        legend.title = element_text(size = 10),
        legend.position = 'none'
        ) +
  labs(x = "", y= ' log2 cpm') + 
  scale_fill_brewer(palette = "Dark2")
ggsave(paste0(figureDir, "GeneGroups_features_geneExpression.pdf"),  width = 4, height = 4)

fts = c('cpg.oe')
yy %>% 
  pivot_longer(cols = fts, names_to = 'features') %>%
  mutate(features = factor(features, levels = fts)) %>%
  mutate(groups = factor(groups, levels = c('DE_up', 'DE_down', 
                                            'house_keep', 'highlyExpr_stable', 'lowlyExpr_stable',  'non_expr'))) %>%
  ggplot(aes(x = features, y=value, fill= groups)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 16), 
        axis.text.y = element_text(angle = 0, size = 16), 
        axis.title =  element_text(size = 16),
        legend.text = element_text(size=10),
        legend.title = element_text(size = 10)
        #legend.position=c(0.3, 0.9)
  ) +
  labs(x = "", y= 'CpG ob/exp ratio') +
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_discrete(labels=c("CpG at promoter"))

ggsave(paste0(figureDir, "GeneGroups_features_CpG.pdf"),  width = 6, height = 4)


fts = c('H3K27me3_mUA', 'H3K4me3_mUA', 'atac_mUA')
yy$rna_mUA[is.na(yy$rna_mUA)] = -6
yy %>% 
  pivot_longer(cols = fts, names_to = 'features') %>%
  mutate(features = factor(features, levels = fts)) %>%
  mutate(groups = factor(groups, levels = c('DE_up', 'DE_down', 'house_keep', 'highlyExpr_stable', 
                                            'lowlyExpr_stable',  'non_expr'))) %>% 
  ggplot(aes(x = features, y=value, fill= groups)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 16), 
        axis.text.y = element_text(angle = 0, size = 16), 
        axis.title =  element_text(size = 16),
        legend.text = element_text(size=10),
        legend.title = element_text(size = 10),
        legend.position='none'
  ) +
  labs(x = "", y= 'log2 cpm at promoter') + 
  scale_fill_brewer(palette = "Dark2")
ggsave(paste0(figureDir, "GeneGroups_features_promoterChromatinFeatures.pdf"),  width = 8, height = 4)


## plot non-high-rank features

fts = c('H3K4me1_mUA')
yy %>% 
  pivot_longer(cols = fts, names_to = 'features') %>%
  mutate(features = factor(features, levels = fts)) %>%
  mutate(groups = factor(groups, levels = c('DE_up', 'DE_down', 'house_keep', 
                                            'highlyExpr_stable', 'lowlyExpr_stable',  'non_expr'))) %>% 
  ggplot(aes(x = features, y=value, fill= groups)) + 
  geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 10), 
        axis.text.y = element_text(angle = 0, size = 10), 
        axis.title =  element_text(size = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size = 10)
        #legend.position=c(0.3, 0.9)
  ) +
  labs(x = "", y= ' log2 cpm at promoter') + 
  scale_fill_brewer(palette = "Dark2")

outDir =  "/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/Figure_S3/"
ggsave(paste0(outDir, "Fig_S3E_GeneGroups_features_promoterChromatinFeatures.pdf_lessImportant.pdf"),  width = 6, height = 4)


fts = c("ZFP42_MA1651.1")
yy %>% 
  pivot_longer(cols = fts, names_to = 'features') %>%
  mutate(features = factor(features, levels = fts)) %>%
  mutate(groups = factor(groups, levels = c('DE_up', 'DE_down', 'house_keep', 'highlyExpr_stable', 
                                            'lowlyExpr_stable',  'non_expr'))) %>%
  ggplot(aes(x = features, y=value, fill= groups)) + 
  #geom_boxplot(outlier.alpha = 0.1) + 
  #geom_jitter(width = 0.1)+
  #geom_violin(width = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 10), 
        axis.text.y = element_text(angle = 0, size = 10), 
        axis.title =  element_text(size = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size = 10)
        #legend.position=c(0.3, 0.9)
  ) +
  labs(x = "", y= '') +
  scale_fill_brewer(palette = "Dark2")

ggsave(paste0(figureDir, "CpGscores_geneGroups.pdf"),  width = 6, height = 4)


########################################################
# not used for now
# regeneration-specific enhancers highlight
########################################################
# means.sel = sample.means[match(rownames(keep), rownames(sample.means)), ]
# 
# dev.mature.maxs = apply(means.sel[, grep('Embryo|Mature_UA', colnames(means.sel))] , 1, max)
# bl.maxs = apply(means.sel[ , grep('BL_UA', colnames(means.sel))], 1, max)
# 
# kk = which(dev.mature.maxs < 2.5 & bl.maxs > 3)
# 
# yy0 = yy[kk, ]
# 
# pheatmap(yy0, 
#          #annotation_row = my_gene_col, 
#          annotation_col = df, show_rownames = FALSE, scale = 'none', 
#          color = col, 
#          show_colnames = FALSE,
#          cluster_rows = TRUE, cluster_cols = FALSE,  
#          #clustering_method = 'complete', cutree_rows = nb_clusters, 
#          #annotation_colors = sample_colors,
#          #clustering_callback = callback,
#          gaps_col = ii.gaps, 
#          filename = paste0(resDir, '/heatmap_regenerationPeaks_edgeRtest_fdr0.05_log2FC.1_regeneartion.specific_v2.pdf'), 
#          width = 8, height = 6)
# 
# if(saveTable){
#   write.csv(data.frame(yy0, stringsAsFactors = FALSE), 
#             file = paste0(resDir, '/regeneration_peaks_regeneration.specific.csv'), 
#             quote = FALSE, row.names = TRUE)
#   
# }

########################################################
########################################################
# Section :  Regeneraiton Motif Analysis and footprint analysis
# 
########################################################
########################################################

##########################################
# peak-to-gene assignment method and characterization
# Compare the corrections between features and rna and it turns out that 
# atac-seq data have the best correlation with RNA-seq data
##########################################
load(file = paste0(RdataDir, '/peak_to_gene_assignment_atacseqPeaks_all.Rdata'))

as_tibble(corrMax) %>% 
  gather(features, corr, 1:4) %>%
  ggplot(aes(x = corr, color = features)) +
  geom_density(size = 1.) +
  scale_color_manual(values = c("#117733", 'red', 'cyan', 'black')) +
  #scale_color_brewer(palette="Dark2") +
  #scale_fill_manual(values=c("#117733",  "blue", 'cyan',  'magenta')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 16), 
        axis.text.y = element_text(angle = 0, size = 16), 
        axis.title =  element_text(size = 16),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position=c(0.4, 0.8)) +
  labs(x = "Correlation between features and expression", y= 'density')

ggsave(paste0(figureDir, "peak_to_gene_assignmenet_correlation_chromatinFeatures_rna_minus_H3K27me3.pdf"), 
       width = 6, height = 4)

check.plot.distance.to.TSS = FALSE
if(check.plot.distance.to.TSS){
  xx = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual_targets.rds'))
  #xx = data.frame(xx[, c(54:64)], stringsAsFactors = FALSE)
  #xx = xx[, c(1:2, 9:11)]
  xx = data.frame(chipseeker = xx$distanceToTSS_chipseeker,  corr.withinTAD = xx$distanceToTSS)
  xx = xx[which(!is.na(xx$corr.withinTAD)), ]
  
  for(n in 1:ncol(xx))
  {
    xx[,n] = sign(xx[,n]) * log10(abs(xx[,n] + 1))
  }
  
  as_tibble(xx) %>%
    gather(cond, dist, 1:2) %>%
    ggplot(aes(x=dist, color = cond)) + 
    geom_density(size = 1.) +
    #geom_histogram(size = 1.) +
    theme_classic() +
    scale_color_brewer(palette="Dark2") +
    labs(x = 'distance between peak to gene (bp in log10)') +
    #geom_vline(xintercept=c(-6, -3,  3, 6), col='gray', size = 1.) +
    theme(legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.position=c(0.2, 0.8),
          plot.margin = margin(),
          axis.text.x = element_text(angle = 0, size = 14), 
          axis.text.y = element_text(angle = 0, size = 14),
          axis.title =  element_text(size = 16),
          #legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    ) + 
    scale_x_continuous(breaks=seq(-7,7, by = 1))
  
  ggsave(paste0(figureDir, "distance_peak_to_target.pdf"), width=6, height = 4)
  
}

##########################################
# regeneration dynamic peaks and gene expressoin of targets 
##########################################
source('Functions_histM.R')

## peak-to-gene assignment of 55k atac-seq peaks
peaks = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual_targets.rds'))

res = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds')) # dynamic atac-seq peak
res = res[which(!is.na(res$clusters)), ]
rownames(res) = gsub('_', '-', rownames(res))

mm = match(rownames(res), rownames(peaks))
length(which(is.na(mm)))

peaks = peaks[mm, ]

rna = readRDS(file = paste0(RdataDir, '/TSS_chromatinFeatures_smartseq2_mature.reg.rds'))
rna = rna[, c(1:4,7, grep('rna_', colnames(rna)))]

mm = match(peaks$targets, rna$geneID)
jj = which(!is.na(mm))
mm = mm[jj]

peaks = peaks[jj,]
rna = rna[mm, ]

ss = apply(as.matrix(rna[, grep('rna_', colnames(rna))]), 1, mean)
sels = which(!is.na(ss) & (rna$groups == 'DE_up'| rna$groups == 'DE_down'))
peaks = peaks[sels, ]
rna = rna[sels, ]

## import dpgp clustered atac peaks
# dpgp = readRDS(file = paste0(RdataDir, '/DPGP_clusters_extended_allPeaks.rds'))
# dpgp$gene = gsub('_', '-', dpgp$gene)
# mm = match(rownames(peaks), dpgp$gene)
# length(which(is.na(mm)))
# peaks$clusters = dpgp$cluster[mm]

yy = peaks[, intersect(grep('atac', colnames(peaks)), grep('rRep', colnames(peaks)))] 
yy <- t(apply(yy, 1, cal_z_score))
colnames(yy) = gsub('_rRep', '', gsub('atac_', '', colnames(yy)))

df <- data.frame(colnames(yy))
rownames(df) = colnames(yy)
colnames(df) = 'sample'

sample_colors = c('darkblue', 'springgreen', 'springgreen3', 'gold2', 'red')[1:nrow(df)]
names(sample_colors) = df$sample

annot_colors = list(
  sample = sample_colors)

#ii.gaps = c(3, 4)
col = colorRampPalette(c("navy", "white", "red3"))(10)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

plt = pheatmap(yy, 
               #annotation_row = my_gene_col, 
               annotation_col = df, show_rownames = FALSE, scale = 'none', 
               color = col, 
               show_colnames = FALSE,
               cluster_rows = TRUE, cluster_cols = FALSE,  
               clustering_method = 'complete', cutree_rows = 6, 
               annotation_colors = annot_colors, 
               clustering_callback = callback,
               treeheight_row = 20,
               legend = TRUE, annotation_legend =  FALSE,
               #gaps_col = ii.gaps, 
               filename = paste0(figureDir, 'dynamic_regeneration_peaks.pdf'), 
               width = 3.2, height = 10)

xx = rna[, grep('rna_', colnames(rna))]
xx = xx[, -c(1:2)]
xx <- t(apply(xx, 1, cal_z_score))
colnames(xx) = gsub('rna_', '', colnames(yy))

gaps.row = cal_clusterGaps(plt, nb_clusters = 6)
pheatmap(xx[plt$tree_row$order, ], 
         #annotation_row = my_gene_col, 
         annotation_col = df, show_rownames = FALSE, scale = 'none', 
         color =  colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10), 
         show_colnames = FALSE,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,  
         #clustering_method = 'complete', cutree_rows = nb_clusters, 
         annotation_colors = annot_colors, 
         width = 3., height = 10, 
         legend = TRUE, annotation_legend =  FALSE,
         #clustering_callback = callback,
         #treeheight_row = 30,
         gaps_row = gaps.row,
         filename = paste0(figureDir, '/putative_targets_dynamic_regeneration_peaks.pdf'))

##########################################
# MARA motif activity analysis for temporally dynamic peaks 
##########################################
source('Functions_MARA.R')
source('Functions_histM.R')
res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks__mUA_regeneration_dev_2Batches.R10723_R7977_peakAnnot_v8.rds'))

Select.dynamic.peaks.for.MARA = TRUE
if(Select.dynamic.peaks.for.MARA){
  res = res[order(-res$log2FC), ]
  
  library(qvalue)
  qv = qvalue(res$pval.lrt)
  res$fdr.lrt = qv$qvalues
  
  # select the temporal dynamic peaks
  fdr.cutoff = 0.05; logfc.cutoff = 1.
  
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
  
  # select samples to use
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
  keep = cal_sample_means(keep, conds = conds)
  rownames(keep) = gsub('_', '-', rownames(keep))
  
  
}

xx = run.MARA.atac.temporal(keep)


##########################################
# plot footprint analysis results from TOBIAS 
##########################################
library(tidyr)
library(dplyr)
require(ggplot2)
library("gridExtra")
library("cowplot")
require(ggpubr)

file.list = list.files(path = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/',
                       'Data/atacseq_using/footprinting/outs_footprints_oneRep/plots'), pattern = '*_all.txt', 
                       full.names = TRUE)

motif = 'RUNX2'

f = file.list[grep(motif, file.list)]
ff = read.table(f, as.is = c(1:3))
xx = matrix(NA, nrow = nrow(ff), ncol = 500)
for(n in 1:nrow(xx))
{
  xx[n, ] = sapply(ff[n, 3],  function(x) as.numeric(unlist(strsplit(as.character(x), ','))))
}

## scale the data according to 50bp in right/left sides
ss = apply(xx[, c(1:50, 450:500)], 1, median)
for(n in 1:nrow(xx)) xx[n,] = xx[n,] - ss[n]
apply(xx[, c(1:50, 450:500)], 1, median)

xx = data.frame(sample = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'), xx,  stringsAsFactors = FALSE)
colnames(xx)[2:501] = as.factor(seq(1, 500, by = 1))

as_tibble(xx) %>%
  gather(bins, signals, 2:501) %>%
  mutate(bins = factor(bins, levels=c(1:500))) %>% 
  mutate(sample = factor(sample, levels = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
  ggplot(aes(y=signals, x=bins, color = sample, group = sample)) + 
  #geom_line(aes(linetype=sample, color = sample), size = 0) +
  geom_smooth(span = 0.1, method = 'loess', se = FALSE, size =1.5) + 
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, size = 14), 
        axis.text.x = element_text(angle = 0, size = 14),
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14)
        #axis.text.x = element_blank(), 
        #axis.ticks = element_blank()
        ) +
  #scale_color_manual(values=c('gray60', 'darkgreen', 'blue', 'gray70','gray80')) +
  #scale_linetype_manual(values=c("dotted", "solid", "solid", "dotted", "dotted")) + 
  labs( x = 'bp from center', y = 'Enrichment score') +
  scale_x_discrete(
    breaks=c(50, 150,  250, 350, 450),
    labels=c("-200", '-100', "0", '100', "200")
  )
  
ggsave(paste0(figureDir, "Footprint_RUNX.pdf"),  width = 6, height = 4)

### motif NRF1
motif = 'NRF1'
f = file.list[grep(motif, file.list)]
ff = read.table(f, as.is = c(1:3))
xx = matrix(NA, nrow = nrow(ff), ncol = 500)
for(n in 1:nrow(xx))
{
  xx[n, ] = sapply(ff[n, 3],  function(x) as.numeric(unlist(strsplit(as.character(x), ','))))
}

## scale the data according to 50bp in right/left sides
ss = apply(xx[, c(1:50, 450:500)], 1, median)
for(n in 1:nrow(xx)) xx[n,] = xx[n,] - ss[n]
apply(xx[, c(1:50, 450:500)], 1, median)

xx = data.frame(sample = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'), xx,  stringsAsFactors = FALSE)
colnames(xx)[2:501] = as.factor(seq(1, 500, by = 1))

as_tibble(xx) %>%
  gather(bins, signals, 2:501) %>%
  mutate(bins = factor(bins, levels=c(1:500))) %>% 
  mutate(sample = factor(sample, levels = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
  ggplot(aes(y=signals, x=bins, color = sample, group = sample)) + 
  #geom_line(aes(linetype=sample, color = sample), size = 0) +
  geom_smooth(span = 0.05, method = 'loess', se = FALSE, size =1.5) + 
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, size = 16), 
        axis.text.x = element_text(angle = 0, size = 16),
        axis.title = element_text(angle = 0, size = 16),
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14)
        #axis.text.x = element_blank(), 
        #axis.ticks = element_blank()
  ) +
  #scale_color_manual(values=c('gray60', 'darkgreen', 'blue', 'gray70','gray80')) +
  #scale_linetype_manual(values=c("dotted", "solid", "solid", "dotted", "dotted")) + 
  labs( x = 'bp from center', y = 'Enrichment score') +
  scale_x_discrete(
    breaks=c(50, 150,  250, 350, 450),
    labels=c("-200", '-100', "0", '100', "200")
  )
ggsave(paste0(figureDir, "Footprint_Motif_", motif, ".pdf"),  width = 6, height = 4)


motif = 'HOXA13'
f = file.list[grep(motif, file.list)]

ff = read.table(f, as.is = c(1:3))
xx = matrix(NA, nrow = nrow(ff), ncol = 500)
for(n in 1:nrow(xx))
{
  xx[n, ] = sapply(ff[n, 3],  function(x) as.numeric(unlist(strsplit(as.character(x), ','))))
}

## scale the data according to 50bp in right/left sides
ss = apply(xx[, c(1:50, 450:500)], 1, median)
for(n in 1:nrow(xx)) xx[n,] = xx[n,] - ss[n]
apply(xx[, c(1:50, 450:500)], 1, median)

xx = data.frame(sample = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'), xx,  stringsAsFactors = FALSE)
colnames(xx)[2:501] = as.factor(seq(1, 500, by = 1))

as_tibble(xx) %>%
  gather(bins, signals, 2:501) %>%
  mutate(bins = factor(bins, levels=c(1:500))) %>% 
  mutate(sample = factor(sample, levels = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
  ggplot(aes(y=signals, x=bins, color = sample, group = sample)) + 
  #geom_line(aes(linetype=sample, color = sample), size = 0) +
  geom_smooth(span = 0.1, method = 'loess', se = FALSE, size =1.5) + 
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, size = 16), 
        axis.text.x = element_text(angle = 0, size = 16),
        axis.title = element_text(angle = 0, size = 16),
        legend.position = 'none',
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14)
        #axis.text.x = element_blank(), 
        #axis.ticks = element_blank()
  ) +
  #scale_color_manual(values=c('gray60', 'darkgreen', 'blue', 'gray70','gray80')) +
  #scale_linetype_manual(values=c("dotted", "solid", "solid", "dotted", "dotted")) + 
  labs( x = 'bp from center', y = 'Enrichment score') +
  scale_x_discrete(
    breaks=c(50, 150,  250, 350, 450),
    labels=c("-200", '-100', "0", '100', "200")
  )

ggsave(paste0(figureDir, "Footprint_Motif_", motif, ".pdf"),  width = 6, height = 4)

motif = 'BACH1'
f = file.list[grep(motif, file.list)]
f = f[grep('BACH1_MA1633.2_BACH1_MA1633.2', f)]

ff = read.table(f, as.is = c(1:3))
xx = matrix(NA, nrow = nrow(ff), ncol = 500)
for(n in 1:nrow(xx))
{
  xx[n, ] = sapply(ff[n, 3],  function(x) as.numeric(unlist(strsplit(as.character(x), ','))))
}

## scale the data according to 50bp in right/left sides
ss = apply(xx[, c(1:50, 450:500)], 1, median)
for(n in 1:nrow(xx)) xx[n,] = xx[n,] - ss[n]
apply(xx[, c(1:50, 450:500)], 1, median)

xx = data.frame(sample = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'), xx,  stringsAsFactors = FALSE)
colnames(xx)[2:501] = as.factor(seq(1, 500, by = 1))

as_tibble(xx) %>%
  gather(bins, signals, 2:501) %>%
  mutate(bins = factor(bins, levels=c(1:500))) %>% 
  mutate(sample = factor(sample, levels = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
  ggplot(aes(y=signals, x=bins, color = sample, group = sample)) + 
  #geom_line(aes(linetype=sample, color = sample), size = 0) +
  geom_smooth(span = 0.03, method = 'loess', se = FALSE, size =1.5) + 
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, size = 16), 
        axis.text.x = element_text(angle = 0, size = 16),
        axis.title = element_text(angle = 0, size = 16),
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14)
        #axis.text.x = element_blank(), 
        #axis.ticks = element_blank()
  ) +
  #scale_color_manual(values=c('gray60', 'darkgreen', 'blue', 'gray70','gray80')) +
  #scale_linetype_manual(values=c("dotted", "solid", "solid", "dotted", "dotted")) + 
  labs( x = 'bp from center', y = 'Enrichment score') +
  scale_x_discrete(
    breaks=c(50, 150,  250, 350, 450),
    labels=c("-200", '-100', "0", '100', "200")
  )

ggsave(paste0(figureDir, "Footprint_Motif_", motif, ".pdf"),  width = 6, height = 4)


motif = 'MAFK'
f = file.list[grep(motif, file.list)]
f = f[grep('MAFK_MA0496.3_MAFK_MA0496.3', f)]

ff = read.table(f, as.is = c(1:3))
xx = matrix(NA, nrow = nrow(ff), ncol = 500)
for(n in 1:nrow(xx))
{
  xx[n, ] = sapply(ff[n, 3],  function(x) as.numeric(unlist(strsplit(as.character(x), ','))))
}

## scale the data according to 50bp in right/left sides
ss = apply(xx[, c(1:50, 450:500)], 1, median)
for(n in 1:nrow(xx)) xx[n,] = xx[n,] - ss[n]
apply(xx[, c(1:50, 450:500)], 1, median)

xx = data.frame(sample = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'), xx,  stringsAsFactors = FALSE)
colnames(xx)[2:501] = as.factor(seq(1, 500, by = 1))

as_tibble(xx) %>%
  gather(bins, signals, 2:501) %>%
  mutate(bins = factor(bins, levels=c(1:500))) %>% 
  mutate(sample = factor(sample, levels = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
  ggplot(aes(y=signals, x=bins, color = sample, group = sample)) + 
  #geom_line(aes(linetype=sample, color = sample), size = 0) +
  geom_smooth(span = 0.03, method = 'loess', se = FALSE, size =1.5) + 
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, size = 16), 
        axis.text.x = element_text(angle = 0, size = 16),
        axis.title = element_text(angle = 0, size = 16),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14),
        legend.position = 'none'
        #axis.text.x = element_blank(), 
        #axis.ticks = element_blank()
  ) +
  #scale_color_manual(values=c('gray60', 'darkgreen', 'blue', 'gray70','gray80')) +
  #scale_linetype_manual(values=c("dotted", "solid", "solid", "dotted", "dotted")) + 
  labs( x = 'bp from center', y = 'Enrichment score') +
  scale_x_discrete(
    breaks=c(50, 150,  250, 350, 450),
    labels=c("-200", '-100', "0", '100', "200")
  )

ggsave(paste0(figureDir, "Footprint_Motif_", motif, ".pdf"),  width = 6, height = 4)


motif = 'RELA'
f = file.list[grep(motif, file.list)]

ff = read.table(f, as.is = c(1:3))
xx = matrix(NA, nrow = nrow(ff), ncol = 500)
for(n in 1:nrow(xx))
{
  xx[n, ] = sapply(ff[n, 3],  function(x) as.numeric(unlist(strsplit(as.character(x), ','))))
}

## scale the data according to 50bp in right/left sides
ss = apply(xx[, c(1:50, 450:500)], 1, median)
for(n in 1:nrow(xx)) xx[n,] = xx[n,] - ss[n]
apply(xx[, c(1:50, 450:500)], 1, median)

xx = data.frame(sample = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'), xx,  stringsAsFactors = FALSE)
colnames(xx)[2:501] = as.factor(seq(1, 500, by = 1))

as_tibble(xx) %>%
  gather(bins, signals, 2:501) %>%
  mutate(bins = factor(bins, levels=c(1:500))) %>% 
  mutate(sample = factor(sample, levels = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
  ggplot(aes(y=signals, x=bins, color = sample, group = sample)) + 
  #geom_line(aes(linetype=sample, color = sample), size = 0) +
  geom_smooth(span = 0.05, method = 'loess', se = FALSE, size =1.5) + 
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, size = 16), 
        axis.text.x = element_text(angle = 0, size = 16),
        axis.title = element_text(angle = 0, size = 16),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14),
        legend.position = 'none'
        #axis.text.x = element_blank(), 
        #axis.ticks = element_blank()
  ) +
  #scale_color_manual(values=c('gray60', 'darkgreen', 'blue', 'gray70','gray80')) +
  #scale_linetype_manual(values=c("dotted", "solid", "solid", "dotted", "dotted")) + 
  labs( x = 'bp from center', y = 'Enrichment score') +
  scale_x_discrete(
    breaks=c(50, 150,  250, 350, 450),
    labels=c("-200", '-100', "0", '100', "200")
  )

ggsave(paste0(figureDir, "Footprint_Motif_", motif, ".pdf"),  width = 6, height = 4)

########################################################
########################################################
# Section : cross-species analysis
# 1) shared motif analysis
# 2) shared RUNX1 and AP1 targets
# 3) Runx1 CNEs analysis 
########################################################
########################################################

##########################################
# regulators cross-species:
# some code origined from https://davemcg.github.io/post/lets-plot-scrna-dotplots/
# or from https://github.com/lengfei5/pallium_evo/blob/main/analysis/ComparativeSpeciesAnalysis.Rmd
##########################################
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork)

## collect MARA results
ntop = 135 
bb = readRDS(file =  paste0(RdataDir, '/MARA_output_Motifs_maxZscore_filteredTFsExpr_top', ntop, '.rds'))

#sels1 = which(bb[, 2]> bb[, 1]|bb[, 3] > bb[, 1]|bb[, 4]> bb[,1] |bb[, 5] > bb[, 1])
sels = which(apply(bb[, c(2:5)], 1, function(x) any(abs(x)>1.8)) == TRUE)

bb = bb[sels, ]
keep = as.matrix(bb[, c(2:5)])
rownames(keep) = bb$gene
keep = data.frame(keep, stringsAsFactors = FALSE)

keep$zscore = apply(keep, 1, function(x) max(abs(x)))
keep = keep[order(-keep$zscore), ]
keep$rank = c(1:nrow(keep)) 

keep = keep[, -c(1:4)]
keep$species = 'axolotl'
keep$gene = rownames(keep)
rownames(keep) = NULL

### zebrafish fin
xx = readRDS(file =  paste0('../results/cross_species_20220621/zebrafish_fin/Rdata', 
                            '/MARA_output_Motifs_filteredTFsExpr.rds'))
test = data.frame(xx[, c(2, 4)])
rownames(test) = xx$gene

test$zscore = apply(test, 1,function(x) max(abs(x)))
test = test[order(-test$zscore), ]

test$rank = c(1:nrow(test)) 
test = test[, -c(1, 2)]
test$species = 'zebrafishFin'
test$gene = rownames(test)
rownames(test) = NULL

keep = data.frame(rbind(keep, test), stringsAsFactors = FALSE)

### zebrafish heart
xx = readRDS(file =  paste0('../results/cross_species_20220621/zebrafish_heart/Rdata/', 
                            '/MARA_output_Motifs_filteredTFsExpr.rds'))
xx = xx[which(xx$motifs != 'FOS_MA1951.1'), ]
test = data.frame(xx[, c(2,3)])
rownames(test) = xx$gene

test$zscore = apply(test, 1,function(x) max(abs(x)))
test = test[order(-test$zscore), ]
#test$zscore = apply(test, 1, max)

test$rank = c(1:nrow(test)) 
test = test[, -c(1, 2)]
test$species = 'zebrafishHeart'
test$gene = rownames(test)
keep = data.frame(rbind(keep, test), stringsAsFactors = FALSE)

### acoel
xx = readRDS(file =  paste0('../results/cross_species_20220621/acoel/Rdata/', 
                            '/MARA_output_sorted.rds'))

#kk1 = apply(as.matrix(xx[, c(1:6)]), 1, function(x) which.max(x) > 1 & max(x) > 2)
#kk2 = apply(as.matrix(xx[, c(7:10)]), 1,  function(x) which.max(x) > 1 & max(x) > 2)
#kk = kk1 + kk2 > 0
#xx = xx[kk, ]
sels = apply(xx[, c(2:6, 8:10)], 1, function(x) any(max(x)>1.8))
xx = xx[sels, ]

xx = xx[which(xx$motifs != 'MEIS1_MA1639.1'), ]
test = data.frame(xx[, c(2:6)])
rownames(test) = xx$gene

test$zscore = apply(test, 1,function(x) max(abs(x)))
test = test[order(-test$zscore), ]

test$rank = c(1:nrow(test)) 
test = test[, c((ncol(test)-1):ncol(test))]
test$species = 'acoel'
test$gene = rownames(test)

keep = data.frame(rbind(keep, test), stringsAsFactors = FALSE)

saveRDS(keep, file = paste0(RdataDir, '/motif_collections_crossSpecies_v2.rds'))

###########
### manually merge similar motifs or from the same family and make the plot 
###########
keep = readRDS(file = paste0(RdataDir, '/motif_collections_crossSpecies_v2.rds'))
keep$rank = as.integer(keep$rank)

## manually merge the similar motifs and discard motifs only once present in acoel, zebrafish fin or heart
counts = table(keep$species, keep$gene)

# merge HOXA HOXC
kk = grep('HOXA13|HOXC13', colnames(counts))
counts[, kk[1]] = apply(counts[, kk], 1, sum)
colnames(counts)[kk[1]] = 'HOXA13/HOXC13'
keep$gene[grep('HOXA13|HOXC13', keep$gene)] = 'HOXA13/HOXC13'
keep[grep('HOXA13|HOXC13', keep$gene), ]


# merge CEBPA CEBPG
kk = grep('CEBPA|CEBPG', colnames(counts))
counts[, kk[1]] = apply(counts[, kk], 1, sum)
colnames(counts)[kk[1]] = 'CEBPA/G'
keep$gene[grep('CEBPA|CEBPG', keep$gene)] = 'CEBPA/G'
keep[grep('CEBPA/G', keep$gene), ]

# merge IRFs
kk = grep('IRF1|IRF2|IRF7', colnames(counts))
counts[, kk[1]] = apply(counts[, kk], 1, sum)
colnames(counts)[kk[1]] = 'IRF1/2/7'
keep$gene[grep('IRF1|IRF2|IRF7', keep$gene)] = 'IRF1/2/7'
keep[grep('IRF1/2/7', keep$gene), ]

# merge RELA RELB
kk = grep('REL|RELA|RELB', colnames(counts))
counts[, kk[1]] = apply(counts[, kk], 1, sum)
colnames(counts)[kk[1]] = 'REL/RELA/RELB'
keep$gene[grep('REL|RELA|RELB', keep$gene)] = 'REL/RELA/RELB'
keep[grep('REL/RELA/RELB', keep$gene), ]

# merge RUNX1 and RUNX2
kk = grep('RUNX1|RUNX2', colnames(counts))
counts[, kk[1]] = apply(counts[, kk], 1, sum)
colnames(counts)[kk[1]] = 'RUNX1/2'
keep$gene[grep('RUNX1|RUNX2', keep$gene)] ='RUNX1/2'
keep[grep('RUNX1/2', keep$gene), ]

# merge SMAD4/2/5
kk = grep('SMAD', colnames(counts))
counts[, kk[2]] = apply(counts[, kk], 1, sum)
colnames(counts)[kk[2]] = 'SMAD'
keep$gene[grep('SMAD', keep$gene)] = 'SMAD'
keep[grep('SMAD', keep$gene), ]

## remove the reduncy of keep after manual changes
gene_species = paste0(keep$gene, '_', keep$species)
idx = c()
for(g in unique(gene_species))
{
  jj = which(gene_species == g)
  if(length(jj) == 1){
    idx = c(idx, jj)
  }else{
    idx = c(idx, jj[which.min(keep$rank)])  
  }
}
keep = keep[idx, ]

ss = apply(counts>0, 2, sum)
genes.NC = colnames(counts)[which(ss > 2)]

mat = counts[, !is.na(match(colnames(counts), genes.NC))] >0
clust = hclust(dist(t(mat))) # hclust with distance matrix
ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

orders = clust$labels[clust$order]
ss = apply(mat, 2, sum)
ss = ss[match(orders, names(ss))]
orders = orders[order(-ss)]
ss = ss[match(orders, names(ss))]

# newlevels = c("RUNX1.2","RELA", 
#              "SMAD",  "MAFK", "NFE2L2",  "HIC1", "FOSL1_JUN", "BACH1", "FOS_JUND",   "PRDM1",  "TCF7L2", 
#              "ZNF341",     "ZBTB32",     "ZBTB26",     "ZBTB18", 
# "TFAP4",      "PRDM5",      "KLF17",    "FOSL2_JUNB", "JDP2",  "TEAD3",      "MAF_NFE2",   "CREB1",     
# "FOS" ,       "THAP11",     "TBP",        "OSR2",       "NR2F1",      "LEF1",       "BCL11B",     "HOXC13",    
# "EGR3",       "NFATC1",     "STAT6",      "STAT2",      "PRDM15",     "POU6F1",     "EGR2",       "NFYA",      
# "ZNF384",     "TCF7",       "ZFX")     
# newlevels = newlevels[c(length(newlevels):1)]

newlevels = clust$labels[clust$order]
newlevels = newlevels[c(1:21, 24, 23, 22)] # manually specify orders

keep$rank = as.numeric(keep$rank)

as_tibble(keep) %>% 
  filter(gene %in% genes.NC) %>% 
  # mutate(gene = factor(gene, levels = clust$labels[clust$order])) %>% 
  mutate(gene = factor(gene, levels = newlevels)) %>% 
  ggplot(aes(y=gene, x = factor(species, levels = c('axolotl', 'zebrafishFin', 'zebrafishHeart', 'acoel')), 
             color = zscore, size = reorder(rank, -rank))) + 
  # ggplot(aes(y=gene, x = factor(species, levels = c('axolotl', 'zebrafishFin', 'zebrafishHeart', 'acoel')), 
  #            color = zscore)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,8), oob = scales::squish, name = 'motif activity') +
  theme(axis.text.x = element_text(angle = 60,  size = 18, hjust = 0.4), 
        axis.text.y = element_text(angle = 0, size = 12), 
        axis.title =  element_text(size = 12),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position='right',
        #plot.margin = margin()
        #legend.key.size = unit(1, 'cm')
        #legend.key.width= unit(1, 'cm')
  )

ggsave(paste0(figureDir, "CrossSpecies_shared_Regulators_v4.pdf"),  width = 8, height = 10)

##########################################
##########################################
# gene targets of RUNX1/2, RELA, AP-1, Bach1, MAFK in axolotl
##########################################
##########################################
library(tidyr)
library(dplyr)
require(ggplot2)
library("gridExtra")
library("cowplot")
require(ggpubr)
library(enrichplot)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(stringr)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(tidyr)
library(dplyr)
require(ggplot2)
library("gridExtra")
library("cowplot")
require(ggpubr)
#require(ChIPseeker)

amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM_geneCorrection_v3.rds'))

DEgenes = tss[which(tss$groups == 'DE_up'|tss$groups == 'DE_down'), ]

dir.list = list.dirs(path = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/',
                                     'Data/atacseq_using/footprinting'), recursive = FALSE,  full.names = TRUE)
dir.list = dir.list[grep('ATAC_footprint_BL_UA_5days|ATAC_footprint_BL_UA_9days', dir.list)]
dir.list = dir.list[grep('_13616', dir.list)]

##### RUNX targets
motif = 'RUNX' 
for(m in 1:length(dir.list))
{
  # m = 1
  subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
  subdir = subdir[grep(motif, subdir)]
  subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
  subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
  
  for(n in 1:length(subdir))
  {
    bed.list = list.files(path = subdir[n], pattern = '*.txt', full.names = TRUE)
    bed.file = bed.list[grep('_overview', bed.list)]
    
    bounds = read.table(bed.file, header = TRUE)
    kk = grep('_bound', colnames(bounds))
    bounds = bounds[which(bounds[,kk] == 1), ]
    bounds = makeGRangesFromDataFrame(bounds, seqnames.field=c("TFBS_chr"),
                                           start.field="TFBS_start", end.field="TFBS_end", strand.field="TFBS_strand ")
    
    cat(basename(dir.list[m]), ' -- ', basename(bed.file), ' -- bound site ', length(bounds), '\n')
    if(m == 1 & n == 1){
      footprint = bounds
    }else{
      footprint = union(footprint, bounds)
    }
  }
}
cat('total footprint found -- ', length(footprint), '\n')

pp.annots = annotatePeak(footprint, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

source('Functions_crossSpecies.R')
p = run_enrichGo_axolotl(pp.annots, distanceToTSS = 2000, regulation = 'both', title = 'RUNX targets')

pdfname = paste0(figureDir, 'GoEnrichment_footprinting_targetGenes_', motif, '_axolotl.pdf')
pdf(pdfname, width = 8, height = 5)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
grid.arrange(p, ncol = 1, nrow = 1)
dev.off()

## save the target genes
pp.annots = as.data.frame(pp.annots)
pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 2000), ]
pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
ggs = unique(as.character(ggs))
cat(length(ggs), ' targets \n')

saveRDS(ggs, file = paste0(RdataDir, '/targetGenes_footprint_', motif, '_axolotl.rds'))

### FOS-JUND motif
motif = 'FOS_JUND'
dir.list = list.dirs(path = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/',
                                   'Data/atacseq_using/footprinting'), recursive = FALSE,  full.names = TRUE)
#dir.list = dir.list[grep('ATAC_footprint_BL_UA_5days|ATAC_footprint_BL_UA_9days|Mature_UA', dir.list)]
dir.list = dir.list[grep('ATAC_footprint_BL_UA_5days|ATAC_footprint_BL_UA_9days', dir.list)]
dir.list = dir.list[grep('_13616', dir.list)]

for(m in 1:length(dir.list))
{
  # m = 1
  subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
  subdir = subdir[grep(motif, subdir)]
  subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
  subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
  
  for(n in 1:length(subdir))
  {
    # n = 1
    
    bed.list = list.files(path = subdir[n], pattern = '*.txt', full.names = TRUE)
    bed.file = bed.list[grep('_overview', bed.list)]
    
    bounds = read.table(bed.file, header = TRUE)
    kk = grep('_bound', colnames(bounds))
    bounds = bounds[which(bounds[,kk] == 1), ]
    bounds = makeGRangesFromDataFrame(bounds, seqnames.field=c("TFBS_chr"),
                                      start.field="TFBS_start", end.field="TFBS_end", strand.field="TFBS_strand ")
    
    cat(basename(dir.list[m]), ' -- ', basename(bed.file), ' -- bound site ', length(bounds), '\n')
    if(m == 1 & n == 1){
      footprint = bounds
    }else{
      if(length(grep('Mature_UA', dir.list[m]))==0) {
        footprint = union(footprint, bounds) 
      }else{
        footprint = GenomicRanges::setdiff(footprint, bounds)
      }
     
    }
  }
}

cat('total footprint found -- ', length(footprint), '\n')

pp.annots = annotatePeak(footprint, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

#plotAnnoBar(pp.annots)
# pp.annots = as.data.frame(pp.annots)
# pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 2000), ]
source('Functions_crossSpecies.R')
p = run_enrichGo_axolotl(pp.annots, distanceToTSS = 5000, regulation = 'DE_down', title.plot = 'AP-1 targets')


pdfname = paste0(figureDir, 'GoEnrichment_footprinting_targetGenes_', motif, '_axolotl.pdf')
pdf(pdfname, width = 8, height = 5)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
grid.arrange(p, ncol = 1, nrow = 1)
dev.off()

## save the target genes
pp.annots = as.data.frame(pp.annots)
pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 2000), ]
pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
ggs = unique(as.character(ggs))
cat(length(ggs), ' targets \n')
saveRDS(ggs, file = paste0(RdataDir, '/targetGenes_footprint_', motif, '_axolotl.rds'))


## motif RELA 
motif = 'REL'
dir.list = list.dirs(path = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/',
                                   'Data/atacseq_using/footprinting'), recursive = FALSE,  full.names = TRUE)
dir.list = dir.list[grep('ATAC_footprint_BL_UA_5days|ATAC_footprint_BL_UA_9days', dir.list)]
dir.list = dir.list[grep('_13616', dir.list)]

for(m in 1:length(dir.list))
{
  # m = 1
  subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
  subdir = subdir[grep(motif, subdir)]
  subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
  subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
  
  for(n in 1:length(subdir))
  {
    # n = 1
    
    bed.list = list.files(path = subdir[n], pattern = '*.txt', full.names = TRUE)
    bed.file = bed.list[grep('_overview', bed.list)]
    
    bounds = read.table(bed.file, header = TRUE)
    kk = grep('_bound', colnames(bounds))
    bounds = bounds[which(bounds[,kk] == 1), ]
    bounds = makeGRangesFromDataFrame(bounds, seqnames.field=c("TFBS_chr"),
                                      start.field="TFBS_start", end.field="TFBS_end", strand.field="TFBS_strand ")
    
    cat(basename(dir.list[m]), ' -- ', basename(bed.file), ' -- bound site ', length(bounds), '\n')
    if(m == 1 & n == 1){
      footprint = bounds
    }else{
      footprint = union(footprint, bounds)
    }
  }
}
cat('total footprint found -- ', length(footprint), '\n')

pp.annots = annotatePeak(footprint, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

#plotAnnoBar(pp.annots)
source('Functions_crossSpecies.R')
p = run_enrichGo_axolotl(pp.annots, distanceToTSS = 2000, regulation = 'DE_down', title.plot = 'RELA targets')

pdfname = paste0(figureDir, 'GoEnrichment_footprinting_targetGenes_', motif, '_axolotl.pdf')
pdf(pdfname, width = 8, height = 5)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
grid.arrange(p, ncol = 1, nrow = 1)
dev.off()

## save the target genes
pp.annots = as.data.frame(pp.annots)
pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 2000), ]
pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
ggs = unique(as.character(ggs))
cat(length(ggs), ' targets \n')
saveRDS(ggs, file = paste0(RdataDir, '/targetGenes_footprint_', motif, '_axolotl.rds'))


motif = 'BACH1'
for(m in 1:length(dir.list))
{
  # m = 1
  subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
  subdir = subdir[grep(motif, subdir)]
  subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
  subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
  subdir = subdir[grep('MA1633.1|MAFK', subdir, invert = TRUE)]
  
  for(n in 1:length(subdir))
  {
    # n = 1
    
    bed.list = list.files(path = subdir[n], pattern = '*.txt', full.names = TRUE)
    bed.file = bed.list[grep('_overview', bed.list)]
    
    bounds = read.table(bed.file, header = TRUE)
    kk = grep('_bound', colnames(bounds))
    bounds = bounds[which(bounds[,kk] == 1), ]
    bounds = makeGRangesFromDataFrame(bounds, seqnames.field=c("TFBS_chr"),
                                      start.field="TFBS_start", end.field="TFBS_end", strand.field="TFBS_strand ")
    
    cat(basename(dir.list[m]), ' -- ', basename(bed.file), ' -- bound site ', length(bounds), '\n')
    if(m == 1 & n == 1){
      footprint = bounds
    }else{
      footprint = union(footprint, bounds)
    }
  }
}
cat('total footprint found -- ', length(footprint), '\n')

pp.annots = annotatePeak(footprint, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

source('Functions_crossSpecies.R')
p = run_enrichGo_axolotl(pp.annots, distanceToTSS = 2000, regulation = 'DE_down', title.plot = 'BACH1 targets')

pdfname = paste0(figureDir, 'GoEnrichment_footprinting_targetGenes_', motif, '_axolotl.pdf')
pdf(pdfname, width = 8, height = 5)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
grid.arrange(p, ncol = 1, nrow = 1)
dev.off()

## save the target genes
pp.annots = as.data.frame(pp.annots)
pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 2000), ]
pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
ggs = unique(as.character(ggs))
cat(length(ggs), ' targets \n')
saveRDS(ggs, file = paste0(RdataDir, '/targetGenes_footprint_', motif, '_axolotl.rds'))

motif = 'MAF_NFE2'
for(m in 1:length(dir.list))
{
  # m = 1
  subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
  subdir = subdir[grep(motif, subdir)]
  subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
  subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
  subdir = subdir[grep('MA1633.1', subdir, invert = TRUE)]
  subdir = subdir[grep('BACH1', subdir, invert = TRUE)]
  
  for(n in 1:length(subdir))
  {
    # n = 1
    
    bed.list = list.files(path = subdir[n], pattern = '*.txt', full.names = TRUE)
    bed.file = bed.list[grep('_overview', bed.list)]
    
    bounds = read.table(bed.file, header = TRUE)
    kk = grep('_bound', colnames(bounds))
    bounds = bounds[which(bounds[,kk] == 1), ]
    bounds = makeGRangesFromDataFrame(bounds, seqnames.field=c("TFBS_chr"),
                                      start.field="TFBS_start", end.field="TFBS_end", strand.field="TFBS_strand ")
    
    cat(basename(dir.list[m]), ' -- ', basename(bed.file), ' -- bound site ', length(bounds), '\n')
    if(m == 1 & n == 1){
      footprint = bounds
    }else{
      footprint = union(footprint, bounds)
    }
  }
}
cat('total footprint found -- ', length(footprint), '\n')

pp.annots = annotatePeak(footprint, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')

source('Functions_crossSpecies.R')
p = run_enrichGo_axolotl(pp.annots, distanceToTSS = 2000, regulation = 'DE_up', title.plot = 'MAFK targets')

pdfname = paste0(figureDir, 'GoEnrichment_footprinting_targetGenes_', motif, '_axolotl.pdf')
pdf(pdfname, width = 8, height = 5)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
grid.arrange(p, ncol = 1, nrow = 1)
dev.off()

## save the target genes
pp.annots = as.data.frame(pp.annots)
pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 2000), ]
pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
ggs = unique(as.character(ggs))
cat(length(ggs), ' targets \n')
saveRDS(ggs, file = paste0(RdataDir, '/targetGenes_footprint_', motif, '_axolotl.rds'))

##########################################
# Collect RUNX targets cross-species and find the shared  
# the network plot was original from 
# https://github.com/lengfei5/pallium_evo/blob/main/analysis/GRN_analysis/moo_branch_grns.R
##########################################
source('Functions_crossSpecies.R')
source('Functions_histM.R')

ax = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_axolotl_geneID_coordinate_4bed.rds'))

targets = data.frame(ax$gene, ax$geneId, ax$transcriptId, stringsAsFactors = FALSE)
targets$ortho_zebrafishFin = NA
targets$ortho_zebrafishHeart = NA
targets$ortho_acoel = NA

## start the orthologues in fish fin
annotFish = 
  read.delim('/Volumes/groups/tanaka/People/current/jiwang/Genomes/zebrafish/GRCz11/annotation_ens_biomart_proteinID.txt')

pp1 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_zebrafishFin_geneID_coordinate_4bed.rds'))
pp2 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_zebrafishHeart_geneID_coordinate_4bed.rds'))
#ort1 = read.delim(paste0(outDir_orthofinder, '/8296_Ambystoma_mexicanum__v__7955_Danio_rerio.tsv'), 
#                  header = TRUE)
pp3 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_acoel_geneID_coordinate_4bed.rds'))


# gg1 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_zebrafish_fin.rds'))
# gg2 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_zebrafishHeart.rds'))
# ggs = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_axolotl.rds'))
# gg3 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_acoel.rds'))
gg1 = get_geneName(pp1$gene)
gg2 = get_geneName(pp2$gene)

targets = intersect(gg1, gg2)
targets = targets[order(targets)]

ggs = get_geneName(ax$gene)
ggs = ggs[grep('AMEX|LOC', ggs, invert = TRUE)]
ggs = as.character(ggs[order(ggs)])

intersect(targets, ggs)
#d = rbind(d, data.frame(gene = intersect(targets, ggs), species = rep('')) )

gg3 = unique(pp3$geneSymbols)

#gg3 = gg3[grep('^g', gg3, invert = TRUE)]
#gg3 = gg3[order(gg3)]
intersect(targets, gg3)
intersect(ggs, gg3)


## manually checking the protein in the same family
# targets[grep('HMG', targets)]; ggs[grep('HMG', ggs)]; gg3[grep('HMG', gg3)]
# targets[which(targets == 'HMGB1')] = 'HMGB1/3'
# ggs[which(ggs == 'HMGB3')] = 'HMGB1/3'
# 
# targets[grep('KDM', targets)]; ggs[grep('KDM', ggs)]; gg3[grep('KDM', gg3)]
# 
# ggs[grep('MRP', ggs)]; gg3[grep('MRP', gg3)]
# 
# ggs[grep('PHF', ggs)]; gg3[grep('PHF', gg3)]
# ggs[grep('PSM', ggs)]; gg3[grep('PSM', gg3)];
# 
# targets[grep('RNF', targets)]; ggs[grep('RNF', ggs)]; gg3[grep('RNF', gg3)];
# targets[grep('SALL', targets)];ggs[grep('SALL', ggs)]; gg3[grep('SALL', gg3)];
# 
# ggs[grep('SLC44', ggs)]; gg3[grep('SLC44', gg3)];
# ggs[grep('TMEM', ggs)]; gg3[grep('TMEM', gg3)]; 

##########
## Visualize the shared targets
#########
library(igraph)
library(ggraph)
library(tidygraph)
library(graphlayouts) 
library(RColorBrewer) # This is the color library


#xx = unique(targets, intersect(targets),)
d = data.frame(regulator = rep('RUNX', length(targets)),
               gene =  targets, species = rep('fish', length(targets)), stringsAsFactors = FALSE)

jj1 = which(!is.na(match(d$gene, ggs)))
jj2 = which(!is.na(match(d$gene, gg3)))
jj = intersect(jj1, jj2)

d$species[jj1] = 'fish.axolotl'
d$species[jj2] = 'fish.acoel'
d$species[jj] = 'all'

d  =d[which(d$species != 'fish'), ]

g <- graph_from_data_frame(d, directed = TRUE)
V(g)$species = d$species[match(V(g)$name, d$gene)]
V(g)$species[1] = 'none'

require(khroma)
highcontrast <- colour("high contrast")
got_palette = highcontrast(3)

ggraph(g, layout='tree', circular=TRUE) +
  #geom_edge_link(aes(edge_width = 0.3), edge_colour = s) + 
  #geom_edge_link(edge_colour = "gray80" ) +
  geom_edge_diagonal(width=0.3) +
  # geom_node_label(aes(label=name), size=12/ggplot2::.pt, 
  #                 label.padding=unit(0.1, 'cm'), label.size=0.3, fill = 'black',  color='white') +
  geom_node_label(aes(label=name, fill = species), size=12/ggplot2::.pt, 
                  fontface = "bold",label.padding=unit(0.1, 'cm'), label.size=0.3) +
  #geom_node_label(aes(label=name, fill = species), size=10/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.2) +
  #geom_node_point(aes(fill = species,size=10/ggplot2::.pt), shape = 21) + 
  #geom_node_text(aes(label=name), size=6/ggplot2::.pt, repel=T, family = "serif")+
  #scale_color_manual(values = got_palette)  +
  guides(fill = NULL) + 
  scale_color_brewer(palette="Set2") +
  scale_x_continuous(expand=c(0.1, 0)) +
  scale_y_continuous(expand=c(0.1, 0)) +
  theme_void()

ggsave(paste0(figureDir, "crossSpecies_RUNX_targets.pdf"), width=7, height = 4)


########################################################
########################################################
# Section : addtional utility analysis 
# 
########################################################
########################################################
meta = read.csv2(file = paste0(tableDir, 'metadata/atac_mature_sampleInfo.csv'))
xx = read.csv(file = paste0(resDir, '/R11637_R12810_atac_QCs_stats.csv'))

mm = match(xx$sampleID, meta$SampleID)

yy = readRDS(file = paste0(RdataDir, '/design_merged_technicalReplicates_Rxxxx_R10723_R11637_R12810.rds'))

load((paste0(RdataDir, '/R11637_atacseq_samples_design_stats.Rdata')))

xx = read.csv2(file = paste0(tableDir, 'metadata/CT_histM_sampleInfos.csv'))
xx = xx[,c(1, 2, 3, 5, 7, 8, 10, 12, 13)]

for(n in c(6:9))
{
  test = as.numeric(as.character(xx[,n]))
  test[which(test<10^6)] = test[which(test<10^6)]*10^6
  xx[,n] = test
}

write.csv2(xx, file = paste0(tableDir, 'metadata/CT_histM_sampleInfos_correctedNumbers.csv'))


xx = read.csv2(file = paste0(tableDir, 'metadata/smartseq2_mature_sampleInfo.csv'))
xx2 = read.csv2(file = paste0(tableDir, 'metadata/smartseq2_reg_sampleInfo.csv'))

xx = xx[, c(1:3)]
xx2 = xx2[, c(1,2, 6)]

xx = rbind(xx2, xx)

yy = read.csv(file = paste0('/Users/jingkui.wang/workspace/imp/positional_memory/results/rnaseq_RNAseqSamples_all/',
  'designSampleInfos_QCstat.csv'))

mm = match(xx$SampleID, yy$SampleID)

xx$total.reads = NA
xx$uniqe.aligned = NA
xx$assigned = NA
xx$total.reads = yy$total.reads[mm]
xx$uniqe.aligned = yy$unique.aligned[mm]
xx$assigned = yy$assigned.reads[mm]

yy = read.csv(file = paste0('~/workspace/imp/positional_memory/results/rnaseq_RNAseqSamples_all/',
                            'QCs_stats_R161513.csv'), header = TRUE, sep = ',', dec = '.', quote = "")

colnames(yy)[1] = 'sampleID'
yy$sampleID = gsub('"', '', yy$sampleID)

jj = which(is.na(xx$total.reads))
mm = match(xx$SampleID[jj], yy$sampleID)
yy = yy[mm, ]
yy = yy[, c(1, 2, 7, 10, 12)]

xx$total.reads[jj] = yy$X..total.reads..
xx$uniqe.aligned[jj] = yy$X..unique.aligned..
xx$assigned[jj] = gsub('"', '', yy$X..assigned.reads...)

write.csv2(xx, file = paste0(tableDir, 'metadata/smartseq2_sampleInfos_correctedNumbers.csv'))

