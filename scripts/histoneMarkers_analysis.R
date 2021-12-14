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
fpm = log2(fpm + 2^-2)





