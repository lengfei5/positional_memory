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
  #
  xx = data.frame(cpm[match(rownames(res), rownames(cpm)), ],  res, stringsAsFactors = FALSE) 
  res = xx
  
  saveRDS(res, file = paste0(RdataDir, '/TSSgenebody_fpm_bc_TMM_combat_DBedgeRtest_', 
                             conds_histM[n_histM], '_', version.analysis, '.rds'))
  
  
  ## select the significant peaks
  fdr.cutoff = 0.05; logfc.cutoff = 1
  select = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
                   (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
                   (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)
  )
  select = select[which(res$max_rpkm[select]>0)]
  cat(length(select), ' DE ', conds_histM[n_histM],  ' \n')
  
  res = res[order(res$adj.P.Val.mHand.vs.mUA), ]
  
  source('Functions_histM.R')
  res$length = width(gene)[match(get_geneID(rownames(res)), gene$gene_id)]
  res$length = res$length + 5000
  res$max_rpkm = res$maxs + log2(10^3/res$length)
  # plot_individual_histMarker_withinATACpeak(res)
  length(which(res$max_rpkm[select]> 0))
  
  
    
}
