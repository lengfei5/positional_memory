##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: cross-species analysis (zebrafish fin and heart, worm)
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jun 21 11:39:12 2022
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
#source('functions_chipSeq.R')

version.analysis = 'cross_species_20220621'
resDir = paste0("../results/", version.analysis)

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0('/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/SupTables/')

saveTables = FALSE

require(DESeq2)
require(GenomicRanges)
require(pheatmap)
library(tictoc)
require(rtracklayer)

library(tidyr)
library(dplyr)
require(ggplot2)
library("gridExtra")
library("cowplot")
require(ggpubr)

species = 'zebrafish_fin'

resDir = paste0(resDir, '/', species)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/', 
                 'other_species_atac_rnaseq/zebrafish_fin_Lee2020')
atacDir = paste0(dataDir, '/atac_seq')
rnaDir = paste0(dataDir, '/rna_seq')

atacMetadata = paste0(atacDir, '/PRJNA523011_atacseq_metadata.txt')
gtf.file = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/zebrafish/GRCz11/Danio_rerio.GRCz11.106.gtf'

########################################################
########################################################
# Section : ATAC-seq data processing
# 
########################################################
########################################################
metadata = read.table(atacMetadata, sep = '\t', header = TRUE)

### manual selection of column of metadata
metadata = data.frame(SampleID = metadata$run_accession,  
                      sample = metadata$sample_title, stringsAsFactors = FALSE)
metadata$condition = gsub('_rep1|_rep2', '', metadata$sample)

design = metadata
rm(metadata)

##########################################
# find the peak consensus
##########################################
peak.files = list.files(path = atacDir, pattern = '*_peaks.xls', full.names = TRUE)

index  = c()
for(n in 1:nrow(design))
{
  test = grep(design$SampleID[n], peak.files)
  if(length(test) != 1) {
    cat(length(test), 'peak files Found \n')
  }else{
    index = c(index, test)
  }
}

peak.files = peak.files[index]

source('functions_chipSeq.R')
peak.consensus = merge.peaks.macs2(peak.files,  select.overlappingPeaks.acrossRep = TRUE,
                                   cc = design$condition, pcutoff = 6)

## save bed for motif scanning and saf file for peak signal counting
export(peak.consensus, format = 'bed', con = paste0(atacDir, '/peakConsensus_for_fimo.bed'))

peak = as.data.frame(peak.consensus)
peak$name = paste0(peak$seqnames, ':', peak$start, '-', peak$end)

SAF = data.frame(GeneID=peak$name, 
                 Chr=peak$seqnames, 
                 Start=peak$start, 
                 End=peak$end, 
                 Strand=peak$strand, 
                 stringsAsFactors = FALSE)

write.table(SAF, file = paste0(atacDir, '/peakConsensus.saf'),   sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE) 

saveRDS(design, file = paste0(RdataDir, '/design_sampleInfo.rds'))

########################################################
########################################################
# Section : Differential binding analysis for atac-seq peaks
# 
########################################################
########################################################
design = readRDS(file = paste0(RdataDir, '/design_sampleInfo.rds'))

xlist<-list.files(path=paste0(atacDir, '/featurecounts_peaks.Q30'),
                  pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

colnames(design)[1] = 'SampleID'

counts = process.countTable(all=all, design = design[, c(1, 3)])

save(design, counts, file = paste0(RdataDir, '/samplesDesign_readCounts.within_peakConsensus.Rdata'))

ss = apply(as.matrix(counts[, -1]), 1, mean)

par(mfrow=c(1,2))
hist(log10(as.matrix(counts[, -1])), breaks = 50, xlab = 'log10(nb of reads within peaks)', main = 'distribution')
plot(ecdf(log10(as.matrix(counts[, -1]) + 0.1)), xlab = 'log10(nb of reads within peaks)', main = 'cumulative distribution')

ss = apply(as.matrix(counts[, -1]), 1, max)
cutoff = 50
hist(log10(ss), breaks = 200)
abline(v = log10(cutoff), col = 'red')
kk = which(ss>cutoff)
length(which(ss>cutoff))

require(ggplot2)
require(DESeq2)

rownames(counts) = counts$gene
dds <- DESeqDataSetFromMatrix(as.matrix(counts[kk, -1]), DataFrame(design), design = ~ condition)

ss = rowSums(counts(dds))
length(which(ss > quantile(ss, probs = 0.6)))

dd0 = dds[ss > quantile(ss, probs = 0.6) , ]
dd0 = estimateSizeFactors(dd0)
sizefactors.UQ = sizeFactors(dd0)

sizeFactors(dds) <- sizefactors.UQ
fpm = fpm(dds, robust = TRUE)

dds <- estimateDispersions(dds, fitType = 'parametric')

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = colnames(design)[3], ntop = 3000, returnData = FALSE)
print(pca)

# ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape = batch)) + 
#   geom_point(size=3) + 
#   geom_text(hjust = 0.2, nudge_y = 0.5, size=3)
# 
# plot(ggp)
ggsave(paste0(resDir, "/PCA_allatacseq_new.old.merged_ntop3000_allSamples.pdf"), width = 10, height = 6)


## quick peak annotation 
require(ChIPpeakAnno)
require(ChIPseeker)

pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
pp$strand = '*'
pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")

gtf = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
pp.annots = annotatePeak(pp, TxDb=gtf, tssRegion = c(-2000, 2000), level = 'transcript')
#plotAnnoBar(pp.annots)

pp.annots = as.data.frame(pp.annots) # sometimes there are less peaks annotated 
rownames(pp.annots) = paste0(pp.annots$seqnames, ':', pp.annots$start, '-', pp.annots$end)

##########################################
# DE test
##########################################
dds$condition = droplevels(dds$condition)
dds$condition <- relevel(dds$condition, ref = "ATAC_sp7ne0dpa")

dds <- DESeq(dds, test="Wald", fitType = c("parametric"))

plotDispEsts(dds, ymin = 10^-3)
resultsNames(dds)

#res <- results(dds)
#res = data.frame(res[, c(2, 5, 6)])
#colnames(res) = paste0(colnames(res), '_LRT')

res = results(dds, name="condition_ATAC_sp7ne4dpa_vs_ATAC_sp7ne0dpa", test = 'Wald')
res <- lfcShrink(dds, coef="condition_ATAC_sp7ne4dpa_vs_ATAC_sp7ne0dpa")
colnames(res) = paste0(colnames(res), "_sp7ne_4dpa.vs.0dpa")
res = data.frame(res[, c(2, 5, 6)])

res.ii = results(dds, contrast = c('condition', 'ATAC_sp7po4dpa', 'ATAC_sp7po0dpa'), test = 'Wald')
res.ii <- lfcShrink(dds, contrast = c('condition', 'ATAC_sp7po4dpa', 'ATAC_sp7po0dpa'))
colnames(res.ii) = paste0(colnames(res.ii), "_sp7po_4dpa.vs.0dpa")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res = data.frame(fpm, res, stringsAsFactors = FALSE)
saveRDS(res, file = paste0(RdataDir, '/fpm_DE_binding_lfcShrink_res.rds'))

##########################################
# select DE peaks
##########################################
res = readRDS(file = paste0(RdataDir, '/fpm_DE_binding_lfcShrink_res.rds'))

fdr.cutoff = 0.01; logfc.cutoff = 1

jj = which((res$padj_sp7ne_4dpa.vs.0dpa < fdr.cutoff & abs(res$log2FoldChange_sp7ne_4dpa.vs.0dpa) > logfc.cutoff) |
             (res$padj_sp7po_4dpa.vs.0dpa < fdr.cutoff & abs(res$log2FoldChange_sp7po_4dpa.vs.0dpa) > logfc.cutoff)
)
cat(length(jj), '\n')

res = res[jj, ]
mm = match(rownames(res), rownames(pp.annots))

res = data.frame(res, pp.annots[mm, ], stringsAsFactors = FALSE)

saveRDS(res, file = paste0(RdataDir, '/DE_atacPeaks_fpm_annotatation.rds'))

########################################################
########################################################
# Section : Motif analysis with MARA
# 
########################################################
########################################################
res = readRDS(file = paste0(RdataDir, '/DE_atacPeaks_fpm_annotatation.rds'))

cat(nrow(res), 'DE peaks found !\n')

res = res[order(-res$log2FoldChange_sp7po_4dpa.vs.0dpa), ]

keep = as.matrix(res[, c(1:8)])

conds = unique(design$condition)

# select samples to use
#keep = keep[, sample.sels]
keep = cal_sample_means(keep, conds = conds)
#rownames(keep) = gsub('_', '-', rownames(keep))

### process fimo output 
fimoDir = paste0(dataDir, '/FIMO_atacPeak/fimo_out')



### run MARA analysis


##########################################
# footprint analysis  
##########################################
Run_footprint_analysis = FALSE
if(Run_footprint_analysis){
  
}

########################################################
########################################################
# Section : RNA-seq analysis
# 
########################################################
########################################################


