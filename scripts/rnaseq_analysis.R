##########################################################################
##########################################################################
# Project: position memory project
# Script purpose: analyze RNA-seq samples
# Usage example:  
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Mar  1 14:07:23 2021
##########################################################################
##########################################################################
rm(list = ls())
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

version.Data = 'rnaseq_RNAseqSamples_all';
version.analysis = paste0("_", version.Data, "_20210301")

## Directories to save results 
design.file = "../exp_design/RNAseq_sampleInfos.xlsx"
#dataDir = "../data/quantseq/featurecounts_R7183"

resDir = paste0("../results/", version.Data)
tabDir =  paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}


##################################################
##################################################
## Section: import design matrix and prepare count table
##################################################
##################################################
Processing.design.matrix = TRUE
if(Processing.design.matrix){
  library(openxlsx)
  design = read.xlsx(design.file, sheet = 1)
}

colnames(design)[1] = 'SampleID'
design$conds = paste0(design$conditions, '_', design$SampleID,  '.batch', design$batch)

# prepare the data table for different batches
batch1 = read.delim('../Data/R10724_rnaseq/nf_out_RNAseq/featureCounts/merged_gene_counts.txt', sep = '\t',  header = TRUE)
batch2 = read.delim('../Data/Rxxxx_rnaseq_old/nf_out_RNAseq/featureCounts/merged_gene_counts.txt', sep = '\t',  header = TRUE)

batch3 = batch2[, grep('HLVGMDRXX', colnames(batch2), invert = TRUE)]
batch2 = batch2[, c(1, grep('HLVGMDRXX', colnames(batch2)))]

source(RNA.functions)
xx1 = process.countTable(all=batch1, design = design[which(design$batch == '3'), c(1:2)], ensToGeneSymbol = FALSE)
colnames(xx1)[-1] = paste0(colnames(xx1)[-1], '.batch3')

xx2 = process.countTable(all=batch2, design = design[which(design$batch == '2'), c(1:2)], merge.technicalRep.sameID = TRUE,
                           ensToGeneSymbol = FALSE)
colnames(xx2)[-1] = paste0(colnames(xx2)[-1], '.batch2')

xx3 = process.countTable(all=batch3, design = design[which(design$batch == '1'), c(1:2)], merge.technicalRep.sameID = TRUE,
                         ensToGeneSymbol = FALSE)

colnames(xx3)[-1] = paste0(colnames(xx3)[-1], '.batch1')

all <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(xx1, xx2, xx3))

# make sure that the sample order is the same in design and all matrix
mm = match(design$conds, colnames(all))
all = all[, c(1, mm)]

write.csv(design, file = paste0(resDir, '/designSampleInfos_QCstat.csv'), row.names = FALSE)

save(design, all, file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))

##########################################
# add stats into design 
##########################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))
stat0 = read.csv(file = '../results/R10724_rnaseq_202102019/QCs_stats.csv')
stat1 = read.csv(file = '../results/Rxxxx.rnaseq.old_202102020/QCs_stats.csv')

xx = data.frame(design, matrix(NA, nrow = nrow(design), ncol = 9), stringsAsFactors = FALSE)
colnames(xx)[7:15] = colnames(stat0)[-c(1:3)]

mm = match(stat0$sampleID, design$SampleID)

xx[mm, c(7:15)] = stat0[, -c(1:3)]

for(n in 1:nrow(xx))
{
  if(is.na(xx$pct.duplication[n])){
    kk = which(stat1$sampleID == xx$SampleID[n])
    if(xx$batch[n] == 1) {
      kk = kk[grep('HLVGMDRXX_', stat1$sample[kk], invert = TRUE)]
      if(length(kk) != 1) cat(n, ' -- Error \n')
      xx[n, c(7:15)] = stat1[kk, c(4:12)]
    }
    if(xx$batch[n] == 2) {
      kk = kk[grep('HLVGMDRXX_', stat1$sample[kk])]
      if(length(kk) != 2) cat(n, ' -- Error \n')
      xx[n, c(7:9, 12)] = apply(stat1[kk, c(4:6, 9)], 2, mean)
      xx[n, c(10:11, 13:15)] = apply(stat1[kk, c(7:8, 10:12)], 2, sum)
    }
    
  }
}

xx = data.frame(xx, stringsAsFactors = FALSE)

xx$pct.assigned.to.features = xx$assigned.reads/xx$total.reads *100
xx$alignment.uniq.rate = xx$unique.aligned/xx$trimmed.reads *100
xx = xx[, c(1:11, 13:15, 12, 17, 16)]

design = xx


save(design, all, file=paste0(RdataDir, 'Design_stats_readCounts_', version.analysis, '.Rdata'))

####################
## first QC for cpm 
####################
load(file=paste0(RdataDir, 'Design_stats_readCounts_', version.analysis, '.Rdata'))

gene.mapping = read.delim('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/AmexT_v47.release_rm.contigs_geneID_geneSymbol.mapping.txt', 
                          sep = '\t', header = TRUE)


design$batch = as.factor(design$batch)

Select.genes.having.symbols = TRUE
if(Select.genes.having.symbols){
  gene.mapping = gene.mapping[which(!is.na(gene.mapping$gene.symbol.nr) | !is.na(gene.mapping$gene.symbol.hs)), ]
}

all = all[!is.na(match(all$gene, gene.mapping$gene.id)), ]


QC.for.cpm = TRUE
if(QC.for.cpm){
  
  #index.qc = c(4)
  source(RNA.functions)
  source(RNA.QC.functions)
  
  kk = grep('Mature_UA|Mature_LA|Mature_Hand', design$conditions)
  kk = kk[which(design$batch[kk] == 1)]
  raw = all[, -1]
  raw = raw[, -kk]
  
  ss = apply(as.matrix(raw), 1, sum)
  
  
  raw = raw[which(ss >10), ]
  
  pdfname = paste0(resDir, "/Data_qulity_assessment_remove.MatureBatch1_filterGenes.withoutSymbols_10", version.analysis, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  
  Check.RNAseq.Quality(read.count=raw, design.matrix = design[-kk, c(1, 2, 5)], 
                       lowlyExpressed.readCount.threshold=500)
  
  dev.off()
  
}

##########################################
# because the gene symbol from nr and hs are not consistent sometimes, so we keep gene.id from AMEXDD60
# Dimensionality reduction to visulize the difference between time points
# Here we select only the batch 3 and batch 2
##########################################
require(ggplot2)
require(DESeq2)
library("dplyr")
library("ggplot2")

raw = as.matrix(all[, -1])
rownames(raw) = all$gene

sels = which(design$batch != 1)

design.matrix = design[sels, ]
raw = raw[, sels]

dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ conditions)

ss = rowSums(counts(dds))

hist(log2(ss), breaks = 200, main = 'log2(sum of reads for each gene)')

cutoff.peak = 2^5
cat(length(which(ss > cutoff.peak)), 'peaks selected \n')

dds <- dds[ss > cutoff.peak, ]

# normalization and dimensionality reduction
dds <- estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = TRUE))
ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_batch2.batch3.pdf"), width=12, height = 8)

library("glmpca")
gpca <- glmpca(counts(dds), L=2)

gpca.dat <- gpca$factors
gpca.dat$conditions <- dds$conditions
gpca.dat$batch <- dds$batch
gpca.dat$name = rownames(gpca.dat)

ggplot(gpca.dat, aes(x = dim1, y = dim2, label = name, color = conditions, shape = batch)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA") +
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5) + 
  ggsave(paste0(resDir, "/GPCAplot_batch2.batch3.pdf"), width=12, height = 8)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$conditions, vsd$batch, sep = " - " )
colnames(sampleDistMatrix) <- NULL

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
mds$name = rownames(mds)

ggplot(mds, aes(x = `1`, y = `2`, label = name, color = conditions, shape = batch)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data") +
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5) 

##########################################
# dynamic genes  
##########################################
dds <- estimateDispersions(dds)
plotDispEsts(dds, ymin = 10^-4)

dds <- nbinomLRT(dds, reduced = ~1 )
res <- results(dds)




