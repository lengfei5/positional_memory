##########################################################################
##########################################################################
# Project: postional memory 
# Script purpose: gene-centric histone mark analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep  5 07:42:37 2022
##########################################################################
##########################################################################
rm(list=ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('functions_chipSeq.R')
source('Functions_atac.R')

version.analysis = 'CT_merged_20220328'
#peakDir = "Peaks/macs2_broad"
saveTable = TRUE

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0(figureDir, 'tables4plots/')

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/histMod_CT_using/'
#gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')

require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require(pheatmap)
library(tictoc)
library(ggrepel)
library(dplyr)
library(tibble)
library(reshape2)
library(tidyverse)

########################################################
########################################################
# Section : prepare gene coordinates from full gtf 
# 
########################################################
########################################################
gtf.all = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/AmexT_v47_Hox.patch.gtf'

amex.all = GenomicFeatures::makeTxDbFromGFF(file = gtf.all)
# tss = GenomicFeatures::promoters(amex.all, upstream = 5000, downstream = 0, use.names = TRUE)
gene = GenomicFeatures::genes(amex.all)
#gene = resize(gene, , fix="end", use.names = TRUE)

gene = as.data.frame(gene)
gene$start[which(gene$strand == '+')] = gene$start[which(gene$strand == '+')] - 5000
gene$end[which(gene$strand == '-')] = gene$end[which(gene$strand == '-')] + 5000

#gene$tx_name = sapply(gene$tx_name, function(x){x = unlist(strsplit(as.character(x), '[|]')); x[length(x)]})
gene$start[which(gene$start<=1)] = 1

SAF = data.frame(GeneID=gene$gene_id, 
                 Chr=gene$seqnames, 
                 Start=gene$start, 
                 End=gene$end, 
                 Strand=gene$strand, 
                 stringsAsFactors = FALSE)

write.table(SAF, file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_histM_both/',
                               'amex6_TSSgene_all_56670genes.saf'), 
            sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE) 


########################################################
########################################################
# Section : process the count data
# 
########################################################
########################################################

##########################################
# process the sample info and count table
##########################################
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

## import sample infos
design = readRDS(file = paste0(RdataDir, '/histM_CT_design_info.rds'))
#design = read.csv(file = paste0(dataDir, "R11876_R12810_R12965_CT_analysis_20220217_QCs_AK.csv"))
#design = design[!is.na(design$sampleID), ]
design$fileName = paste0(design$condition, '_', design$sampleID)
design = design[, c(1:3, 15, 5, 16, 4, 6:14)]
colnames(design)[5] = 'batch'

# design.histM = design[, c(1:5)]
# design.histM = design.histM[grep('rRep', design.histM$batch), ]
# 
# design = readRDS(file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata', 
#                                '/design_sels_bc_TMM_combat_mUA_regeneration_dev_2Batches.R10723_R7977',
#                                'Rxxxx_R10723_R11637_R12810_atac', '.rds'))
# design = design[, c(1, 2, 6)]
# design$sample = design$condition
# design$marks = 'atac'
# design = design[, c(1, 4,5, 2:3)]
# 
# colnames(design.histM) = colnames(design)
# 
# design = rbind(design, design.histM)
# design$sample = gsub('BL_UA_13days_distal', 'BL13days.dist', design$sample)
# design$sample = gsub('BL_UA_13days_proximal', 'BL13days.prox', design$sample)
# design$sample = gsub('BL_UA_9days', 'BL9days', design$sample)
# design$sample = gsub('BL_UA_5days', 'BL5days', design$sample)
# design$sample = gsub('Mature_UA', 'mUA', design$sample)

#design = design[grep('Embryo', design$sample, invert = TRUE), ]
# design = design[grep('IgG', design$marks, invert = TRUE), ]
# design$condition = paste0(design$marks, '_', design$sample)

## read counts within peaks and process the counts
dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_histM_both/'
xlist<-list.files(path=paste0(dataDir, '/featurecounts_peaks.Q30'),
                  pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

xlist<-list.files(path=paste0(dataDir, 'featurecounts_peaks.Q30'),
                  pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

sels = c()
for(n in 1:nrow(design)) sels = c(sels, grep(design$SampleID[n], xlist))
xlist = xlist[sels]

all = cat.countTable(xlist, countsfrom = 'featureCounts')

colnames(design)[1] = 'SampleID'

counts = process.countTable(all=all, design = design[, c(1,4)])

save(design, counts, file = paste0(RdataDir, '/histM_design_mature_regeneration_counts_TSSgenebody.Rdata'))

##########################################
# reload the save design and counts for TSS and gene body
# and select only the limb fibrobalst exprssing gene
##########################################
## select only expressed gene
gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
gene = GenomicFeatures::genes(amex)

load(file = paste0(RdataDir, '/histM_design_mature_regeneration_counts_TSSgenebody.Rdata'))
rownames(counts) = counts$gene
counts = as.matrix(counts[, -1])

mm = match(gene$gene_id, rownames(counts))
counts = counts[mm, ]

## use geneID 
annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
names = rownames(counts)
mm = match(rownames(counts), annot$geneID)
ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
names[!is.na(mm)] = ggs[!is.na(mm)]

rownames(counts) = names

dds <- DESeqDataSetFromMatrix(as.matrix(counts), DataFrame(design[, c(1:5)]), design = ~ condition)

save(dds,  design, file = paste0(RdataDir, '/histM_design_mature_regeneration_counts_TSSgenebody_geneNames.Rdata'))
##########################################
# for each markers, normalization, sample filtering, and batch correction
##########################################
load(file = paste0(RdataDir, '/histM_design_mature_regeneration_counts_TSSgenebody_geneNames.Rdata'))

conds = c('H3K27me3', 'H3K4me3', 'H3K4me1')

Make.pca.plots = FALSE

for(n in 1:length(conds))
{
  # n = 3
  # sels = which(dds$marks == conds[n] & dds$SampleID != '163655' & dds$SampleID != '185743' )
  sels = which(dds$marks == conds[n] & dds$SampleID != '163655')
  
  cat(n, ' --', conds[n], '--', length(sels), 'samples\n')
  
  ddx = dds[, sels]
  design.sel = design[sels, ]
  
  ddx$condition = droplevels(ddx$condition)
  ddx$batch = droplevels(ddx$batch)
  
  ## normalization  
  ss = rowSums(counts(ddx))
  
  hist(log10(ss), breaks = 100,  add = FALSE);
  abline(v = c(log10(50), 2, log10(500), 3), col = 'red', lwd = 2.0)
  
  if(conds[n] == 'H3K27me3')  dd0 = ddx[ss > 100, ]
  if(conds[n] == 'H3K4me3'|conds[n] == 'H3K4me1') dds = ddx[ss>10^3, ]
  dd0 = estimateSizeFactors(dd0)
  sizeFactors(ddx) <- sizeFactors(dd0)
  
  fpm = fpm(ddx, robust = TRUE)
  log2(range(fpm[which(fpm>0)]))
  
  fpm = log2(fpm + 2^-4)
  colnames(fpm) = paste0(design.sel$condition, '_', design.sel$SampleID, '_', design.sel$batch)
  colnames(ddx) = colnames(fpm)
  
  ## double check the sample qualities with PCA and scatterplots of replicates
  if(Make.pca.plots){
    library(ggrepel)
    library(dplyr)
    library(tibble)
    
    vsd <- varianceStabilizingTransformation(ddx, blind = FALSE)
    
    pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE, ntop = 3000))
    pca2save$name = paste0(design.sel$sample, '_', design.sel$SampleID)
    
    ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
      geom_point(size=3) + 
      geom_text_repel(size = 4.0)
    #geom_text(hjust = 0.4, nudge_y = 1.0, size=3.5)
    ggsave(paste0(resDir, "/histM_PCA_TSSgenebody_", conds[n], ".pdf"), width=10, height = 6)
    
    pdfname = paste0(resDir, "/histM_TSSgenebody_", conds[n], ".pdf")
    pdf(pdfname, width = 10, height = 10)
    
    plot.pair.comparison.plot(fpm[c(1:20000), grep('_mUA', colnames(fpm))], linear.scale = FALSE)
    plot.pair.comparison.plot(fpm[c(1:20000), grep('_mLA|_mHand', colnames(fpm))], linear.scale = FALSE)
    plot.pair.comparison.plot(fpm[c(1:20000), grep('_BL5|_BL9', colnames(fpm))], linear.scale = FALSE)
    plot.pair.comparison.plot(fpm[c(1:20000), grep('_BL13', colnames(fpm))], linear.scale = FALSE)
    
    dev.off()
    
  }
  
  ## batch correction
  Batch.Correction.combat = TRUE
  if(Batch.Correction.combat){
    require(edgeR)
    require(sva)
    source('Functions_atac.R')
    
    d <- DGEList(counts=counts(ddx), group=ddx$condition)
    
    if(conds[n] == 'H3K27ac' | conds[n] == 'IgG'){
      tmm <- calcNormFactors(d, method='TMMwsp') # TMMwsp used for H3K27ac
    }else{
      tmm <- calcNormFactors(d, method='TMM')
    }
    
    tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
    
    design.sel$condition = droplevels(design.sel$condition)
    design.sel$batch = droplevels(design.sel$batch)
    
    bc = as.factor(design.sel$batch)
    mod = model.matrix(~ as.factor(condition), data = design.sel)
    
    make.pca.plots(tmm, ntop = 5000, conds.sel = as.character(unique(design.sel$condition)))
    ggsave(paste0(resDir, "/histMarker_TSSgenebody", conds[n], "_beforeBatchCorrect_",  version.analysis, ".pdf"), 
           width = 16, height = 14)
    
    #tmm = as.matrix(tmm)
    #vars = apply(tmm, 1, var)
    #tmm = tmm[which(vars>0.5), ]
    # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
    # because there is no better batche other others 
    #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
    fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
    
    make.pca.plots(fpm.bc, ntop = 5000, conds.sel = as.character(unique(design.sel$condition)))
    ggsave(paste0(resDir, "/histMarker_", conds[n], "_afterBatchCorrect_",  version.analysis, ".pdf"), width = 10, height = 6)
    
    fpm = fpm.bc
    rm(fpm.bc)
    
    saveRDS(fpm, file = paste0(RdataDir, '/fpm_bc_TMM_combat_TSSgenebody', conds[n], '_', version.analysis, '.rds'))
    saveRDS(design.sel, file = paste0(RdataDir, '/design.sels_bc_TMM_combat_TSSgenebody', 
                                      conds[n], '_', version.analysis, '.rds'))
    
  }
}

########################################################
########################################################
# Section III : DE test for segment-specific histone markers
# 
########################################################
########################################################
library(edgeR)
library(qvalue)
require(corrplot)
require(pheatmap)
require(RColorBrewer)

gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
gene = GenomicFeatures::genes(amex)

conds_histM = c('H3K27me3','H3K4me1', 'H3K4me3')
ll = width(gene)

for(n_histM in 1:length(conds_histM))
{
  # n_histM = 1
  design.sel = readRDS(file = paste0(RdataDir, '/design.sels_bc_TMM_combat_TSSgenebody', 
                              conds_histM[n_histM], '_', version.analysis, '.rds'))
  cpm = readRDS(file = paste0(RdataDir, '/fpm_bc_TMM_combat_TSSgenebody', conds_histM[n_histM], '_', version.analysis, '.rds'))
  
  ### select the samples and extract sample means
  conds = c("mUA", "mLA", "mHand")
  sample.sels = c();  
  cc = c()
  sample.means = c()
  
  for(n in 1:length(conds)) 
  {
    kk = grep(conds[n], colnames(cpm))
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
    if(length(kk)>1) {
      sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
    }else{
      sample.means = cbind(sample.means, cpm[, kk])
    }
    
  }
  colnames(sample.means) = conds
  
  cpm = cpm[, sample.sels]
  design.sel = design.sel[sample.sels, ]
  
  logCPM = cpm
  f = factor(cc, levels= c('mUA', 'mLA', 'mHand'))
  
  mod = model.matrix(~ 0 + f)
  colnames(mod) = c('mUA', 'mLA', 'mHand')
  
  #To make all pair-wise comparisons between the three groups one could proceed
  fit <- lmFit(logCPM, mod)
  contrast.matrix <- makeContrasts(mLA - mUA, 
                                   mHand - mUA, 
                                   mHand - mLA, levels=mod)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  res = data.frame(fit2$p.value)
  colnames(res) = paste0(c('mLA.vs.mUA', 'mHand.vs.mUA', 'mHand.vs.mLA'), '.pval')
  
  xx = topTable(fit2, coef = 1, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.mLA.vs.mUA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  xx = topTable(fit2, coef = 2, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.mHand.vs.mUA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  xx = topTable(fit2, coef = 3, number = nrow(res))
  xx = xx[, c(1, 4, 5)]
  colnames(xx) = paste0(colnames(xx), '.mHand.vs.mLA')
  res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
  
  res$pval.mean = apply(as.matrix(res[, grep('P.Value.', colnames(res))]), 1, function(x) return(mean(-log10(x))))
  res$fdr.mean = apply(as.matrix(res[, grep('adj.P.Val', colnames(res))]), 1, function(x) return(mean(-log10(x))))
  res$logFC.mean =  apply(as.matrix(res[, grep('logFC', colnames(res))]), 1, function(x) return(mean(abs(x))))

  res$log2fc = apply(sample.means, 1, function(x) max(x) - min(x))
  res$maxs = apply(sample.means, 1, max)
  res$mins = apply(sample.means, 1, min)
  
  source('Functions_histM.R')
  res$length = width(gene)[match(get_geneID(rownames(res)), gene$gene_id)]
  res$length = res$length + 5000
  res$max_rpkm = res$maxs + log2(10^3/res$length)
  
  #
  xx = data.frame(cpm[match(rownames(res), rownames(cpm)), ],  res, stringsAsFactors = FALSE) 
  res = xx
  
  ## select the significant peaks
  fdr.cutoff = 0.05; logfc.cutoff = 1
  select = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
                   (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)
  )
  select = select[which(res$max_rpkm[select]>0)]
  cat(length(select), ' DE ', conds_histM[n_histM],  ' \n')
  
  res = res[order(res$adj.P.Val.mHand.vs.mUA), ]
  
  # plot_individual_histMarker_withinATACpeak(res)
  
  saveRDS(res, file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_', 
                             conds_histM[n_histM], '_', version.analysis, '.rds'))
    
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
  
  ## call function for heamtap 
  source('Functions_plots.R')
  # subclustering.postional.histM.postioinalAtacPeaks()
  
  ##########################################
  # merge first the histM tables
  ##########################################
  conds_histM = c('H3K27me3', 'H3K4me1', 'H3K4me3')
  for(n in 1:length(conds_histM)){
    # n = 2
    if(n == 1){
      res = readRDS(file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_', 
                                  conds_histM[n], '_', version.analysis, '.rds'))
      colnames(res) = paste0(colnames(res), '.', conds_histM[n])
      
    }else{
      xx = readRDS(file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_', 
                                 conds_histM[n], '_', version.analysis, '.rds'))
      colnames(xx) = paste0(colnames(xx), '.', conds_histM[n])
      xx = xx[match(rownames(res), rownames(xx)),]
      res = data.frame(res, xx, stringsAsFactors = FALSE)
    }
  }
 
  saveRDS(res, file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_all3histM_', version.analysis, '.rds'))
  
  ##########################################
  # select the genes based on the H3K27me3 
  ##########################################
  res = readRDS(file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_all3histM_', version.analysis, '.rds'))
  source('Functions_histM.R')
  
  ## select the significant peaks
  fdr.cutoff = 0.1; logfc.cutoff = 1
  select = which((res$adj.P.Val.mLA.vs.mUA.H3K27me3 < fdr.cutoff & abs(res$logFC.mLA.vs.mUA.H3K27me3) > logfc.cutoff) |
                   (res$adj.P.Val.mHand.vs.mUA.H3K27me3 < fdr.cutoff & abs(res$logFC.mHand.vs.mUA.H3K27me3) > logfc.cutoff)|
                   (res$adj.P.Val.mHand.vs.mLA.H3K27me3 < fdr.cutoff & abs(res$logFC.mHand.vs.mLA.H3K27me3) > logfc.cutoff)
  )
  select = select[which(res$max_rpkm.H3K27me3[select]>0.6)]
  cat(length(select), ' DE H3K27me3 ',  ' \n')
  
  
  yy0 = res[select, ]
  range <- 2.0
  
  yy = matrix(NA, nrow = nrow(yy0), ncol = length(conds_histM)*3)
  rownames(yy) = rownames(yy0)
  nms = c()
  
  for(n in 1:length(conds_histM)) # transform the data
  {
    # n = 1
    jj0 = grep(paste0(conds_histM[n], '_m'), colnames(yy0))
    
    test = yy0[, jj0]
    test = test[, grep('_rRep', colnames(test), invert = TRUE)]
    
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
           show_rownames=TRUE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           gaps_col = gaps_col,
           legend = TRUE,
           treeheight_row = 15,
           annotation_legend = FALSE, 
           #annotation_colors = annot_colors,
           clustering_callback = callback,
           #breaks = seq(-2, 2, length.out = 8),
           clustering_method = 'complete', cutree_rows = 3,
           breaks = seq(-range, range, length.out = 8),
           #gaps_row =  gaps.row, 
           legend_labels = FALSE,
           width = 5, height = 10, 
           filename = paste0(figureDir, 'heatmap_histoneMarker_geneCentric_DE.pdf'))
  
  
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


