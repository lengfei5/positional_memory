##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: analyze the histone markers to test the chromatin states invovled in postional memory and regeneration
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Dec 14 10:04:12 2021
##########################################################################
##########################################################################

rm(list=ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('functions_chipSeq.R')
source('Functions_atac.R')

version.analysis = 'atac_rna_chipseq_analysis_20211007'
#peakDir = "Peaks/macs2_broad"
saveTable = TRUE

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0(figureDir, 'tables4plots/')

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')

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
# Section : quantify histone markers around TSS/enhancers
# prepare the SAF file for promoters
########################################################
########################################################
promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')


# clean peaks
peaks = promoters
peaks = data.frame(peaks)
colnames(peaks)[c(1:3)] = c("chr", "start", "end")

dim(peaks)

peaks$peak.name = paste0(peaks$geneID, '_', peaks$geneSymbol, '_', peaks$chr,  "_", peaks$start, "_", peaks$end)

jj = match(unique(peaks$peak.name), peaks$peak.name)
peaks = peaks[jj, ];
dim(peaks)

peaks = peaks[which(peaks$chr != 'chrM'), ]

dim(peaks)


peakDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11876_cut.run/nf_out/'
# preapre SAF input for featureCounts
require(Rsubread)

#counts = quantify.signals.within.peaks(peaks = peak.merged, bam.list=bam.list, rpkm.normalization = FALSE, isPairedEnd = FALSE)
SAF = data.frame(GeneID=peaks$peak.name, 
                 Chr=peaks$chr, 
                 Start=peaks$start, 
                 End=peaks$end, 
                 Strand=peaks$strand, stringsAsFactors = FALSE)

write.table(SAF, file = paste0(peakDir, 'amex6_TSS_all.saf'), sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE) 

##########################################
# import counts by feature counts around TSS
##########################################
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

design = read.delim(file = 
        '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11876_cut.run/sampleInfos_R11876_parsed.txt',
        header = TRUE, sep = '\t')

#mm = match(design$sampleID, as.character(c(161523:161526)))
#design = design[is.na(mm), ]
xlist<-list.files(path=paste0(peakDir, '/featurecounts_peaks.Q30'),
                  pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

colnames(design)[1] = 'SampleID'
#design = design[grep('90392|90393|108072|108073|108070|108071', design$SampleID, invert = TRUE),]

counts = process.countTable(all=all, design = design[, c(1,2)])

save(design, counts, file = paste0(RdataDir, '/histoneMarkers_samplesDesign_readCounts_axolotlAllTSS.2kb.Rdata'))

##########################################
#  signal normalization
##########################################
load(file = paste0(RdataDir, '/histoneMarkers_samplesDesign_readCounts_axolotlAllTSS.2kb.Rdata'))
design$condition = design$fileName
design$condition = gsub('K27me3et', 'K27me3', design$condition)

cc.sels = apply(expand.grid(c('UA', 'LA', 'Hand'),  c('K4me3', 'K27me3')), 1, paste, collapse="_")
sels = match(cc.sels, design$condition)

design = design[sels, ]
counts = counts[, c(1, (sels+1))]


ss = apply(as.matrix(counts[, -1]), 1, mean)

par(mfrow=c(1,2))
hist(log10(as.matrix(counts[, -1])), breaks = 100, xlab = 'log10(nb of reads within peaks)', main = 'distribution')
plot(ecdf(log10(as.matrix(counts[, -1]) + 0.1)), xlab = 'log10(nb of reads within peaks)', main = 'cumulative distribution')

ss = apply(as.matrix(counts[, -1]), 2, sum)
design$usable.reads.withinPeaks = ss

rownames(counts) = counts$gene
dds <- DESeqDataSetFromMatrix(as.matrix(counts[, -1]), DataFrame(design), design = ~ condition)

ss = rowSums(counts(dds))
#length(which(ss > quantile(ss, probs = 0.4)))

dd0 = dds[ss > 20, ]
dd0 = estimateSizeFactors(dd0)
sizefactors.UQ = sizeFactors(dd0)

sizeFactors(dds) <- sizefactors.UQ
fpm = fpm(dds, robust = TRUE)
colnames(fpm) = design$condition

saveRDS(fpm, file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb.rds'))
##########################################
# test/explore histone markers for postional-related genes
# and regeneration-related gened
##########################################
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

########################################################
########################################################
# Section : Consider the grouping: 
# inactive promoters, active promoter, house-keeping genes, mature-specific genes,
# regeneration genes, positional genes
########################################################
########################################################
fpm = readRDS(file = paste0(RdataDir,  '/histoneMarkers_normSignals_axolotlAllTSS.2kb.rds'))
res = data.frame(log2(fpm + 1), stringsAsFactors = FALSE)
res$gene = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[2]})
res$geneID  = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})

##########################################
# clean TSS, for each gene, only keep the TSS with sigals if there are mulitple ones
##########################################
res$used = NA

Clean.TSS = FALSE
if(Clean.TSS){
  ggs = unique(res$geneID)
  
  for(n in 1:length(ggs))
  {
    # n = 1
    jj = which(res$geneID == ggs[n])  
    if(length(jj)>0){
      hist.max = apply()
    } 
    
  }
}

##########################################
# first try to define groups 
##########################################
res$x1 = res$UA_K4me3
res$x2 = res$UA_K27me3
#rr = res$x1 - res$x2

#examples.sel = which(abs(rr)>2 & (res$x1>5|res$x2>5))
examples.sel = c()
examples.sel = unique(c(examples.sel, grep('HOXA13|HOXA11|HOXA9|HOXD13|HOXD11|HOXD9|MEIS', res$gene)))

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


load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
hkgs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
nonexp = controls.tissue$geneIDs[which(controls.tissue$tissues != 'housekeeping')]

res$groups = 'limb'
res$groups[!is.na(match(res$geneID, hkgs))] = 'house.keep'
res$groups[!is.na(match(res$geneID, nonexp))] = 'other.tissues'

##########################################
# redefine different groups with RNA-seq data 
##########################################
Redefine.gene.groups.with.RNAseq = TRUE
if(Redefine.gene.groups.with.RNAseq){
  load(file = paste0('~/workspace/imp/positional_memory/results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/RNAseq_design_dds.object.Rdata'))
  annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  dds = dds[,which(dds$condition != 'Mature_HEAD')]
  dds$condition = droplevels(dds$condition)
  
  
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
res[,c(1, 4, 7:12)] %>% 
  pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
ggplot(aes(x = groups, y=value, fill=markers)) + 
  geom_boxplot(outlier.shape = NA) + 
  #geom_jitter(width = 0.1)+
  theme(axis.text.x = element_text(angle = 0, size = 12))







