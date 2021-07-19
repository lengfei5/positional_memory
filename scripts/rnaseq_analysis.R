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
require(openxlsx)
require(ggplot2)
require(DESeq2)
library("dplyr")
require(gridExtra)

version.Data = 'rnaseq_Rxxxx.old_R10724_R161513';
version.analysis = paste0("_", version.Data, "_20210628")

## Directories to save results 
design.file = "../exp_design/RNAseq_sampleInfos.xlsx"
#dataDir = "../data/quantseq/featurecounts_R7183"

resDir = paste0("../results/", version.Data)
tabDir =  paste0(resDir, "/tables/")
tfDir = '~/workspace/imp/positional_memory/results/motif_analysis'
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}


########################################################
########################################################
# Section I : data processing and sequencing quality controls
# 
########################################################
########################################################

##########################################
# prepare design matrix, count table and QC table
# here the statistics were collected batch by batch (request by request)
# merge in the later steps
##########################################
#dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/Rxxxx_rnaseq_old/'
#design = read.table(paste0(dataDir, 'sampleInfos_parsed.txt'), sep = '\t', header = TRUE)
dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R161513_rnaseq/'
design = read.csv(file = paste0(dataDir, 'sampleInfos.csv'))

colnames(design) = c('sampleID', 'fileName')
stats = read.delim(paste0(dataDir, 'nf_out_RNAseq/MultiQC/multiqc_data/multiqc_general_stats.txt'), sep = '\t', 
                   header = TRUE)
alignment = read.delim(paste0(dataDir, 'nf_out_RNAseq/MultiQC/multiqc_data/multiqc_hisat2.txt'), sep = '\t')

stats = stats[, c(1, 2, 3, 4, 6, 8, 9, 10)]
colnames(stats) = c('sample', 'pct.duplication', 'pct.GC', 'avg.seq.length', 'total.reads', 
                    'pct.assign', 'assigned.reads', 'alignment.rate')

stats = data.frame(stats, alignment[match(stats$sample, alignment$Sample), c(2, 4, 5)], stringsAsFactors = FALSE)
colnames(stats)[c(9:11)] = c('trimmed.reads', 'unique.aligned', 'multimapper')
stats = stats[, c(1:5, 9, 8, 10, 11, 7)]

ii = c()
jj = c()
for(n in 1:nrow(design))
{
  # n = 1;
  cat(n, '\n')
  kk = grep(design$sampleID[n], stats$sample)
  ii = c(ii, rep(n, length(kk)))
  jj = c(jj, kk)
  #kk = c(kk, grep(design$sampleID[n], stats$sample))
  #kk = c(kk, which(design$sampleID == ))
}

xx = data.frame(design[ii, ], stats[jj, ], stringsAsFactors = FALSE)
xx = xx[order(xx$fileName), ]

# write.csv(xx, file = paste0(resDir, '/QCs_stats.csv'), row.names = FALSE)

#xx2 = xx
#xx1 = xx;
#design = xx 


########################################################
# saturation curve from rseqc
# the r code from rseqc output
########################################################
Saturation.test = FALSE
if(Saturation.test){
  
  rseqc.file = list.files('../Data/R10724_rnaseq/saturation_rseqc', pattern = 'junctionSaturation_plot.r', 
                          full.names = TRUE)
  library(stringr)
  
  yy = c()
  for(n in 1:length(rseqc.file))
  {
    cat(n, '\n')
    xx = read.delim(rseqc.file[n])
    xx = xx[grep('y=', xx[, 1 ]), ]
    #xx = gsub('y=c', '', xx)
    # Get the parenthesis and what is inside
    k <- str_extract_all(xx, "\\([^()]+\\)")[[1]]
    # Remove parenthesis
    k <- substring(k, 2, nchar(k)-1)
    #k = gsub('["]', '', k)
    k = as.numeric(unlist(strsplit(as.character(k), ',')))
    yy = rbind(yy, k)
  }
  
  rownames(yy) = gsub('_junction.junctionSaturation_plot.r', '', basename(rseqc.file))
  
  
  pdfname = paste0(resDir, '/saturation_curve_rseqc_knownJunctions.pdf')
  pdf(pdfname, width = 16, height = 8)
  par(cex = 1.0, las = 1,  mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0), tcl = -0.3)
  
  yy = yy[which(rownames(yy) != '136150_TTAACCTTCGAGGCCAGACA_HNF3KDSXY_3_20201223B_20201223'), ]
  span = 0.75
  # saturation curve with nb of peaks
  xlims = c(0, 120)
  ylims = range(yy/10^3)
  frac = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)/100
  
  library(RColorBrewer)
  cols = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(nrow(yy))
  plot(0, 0, xlim = xlims, ylim = ylims, type ='n', xlab = 'nb of TOTAL reads (Million)', 
       ylab = 'nb of known junctions (K)', main = paste0('saturation curve from rseqc'))
  abline(v = c(20, 30, 40, 50), col = 'blue', lwd = 1.0, lty =2)
  
  #legend('topleft', legend = sample.uniq, col = cols, bty = 'n', lwd = 2.0, cex = 0.7)
  
  for(n in 1:nrow(yy))
  {
    # n = 1
    cat(n, '\n')
    
    kk = which(design$sample == rownames(yy)[n])
    
    satt = data.frame(nb.reads = design$total.reads[kk]*frac/10^6, nb.junctions = yy[n, ]/10^3)
    
    points(satt[,1], satt[,2], type= 'p', col = cols[n])
    loessMod <- loess(nb.junctions ~ nb.reads, data=satt, span=span)
    smoothed <- predict(loessMod)
    lines(smoothed, x=satt$nb.reads, col=cols[n], lwd = 3.0)
    
    text(satt[nrow(satt), 1], smoothed[length(smoothed)], labels = paste0(design$fileName[kk], '_', design$sampleID[kk]), 
         cex = 0.7, pos = 4, offset = 0.2)
    
  }
  
  dev.off()
  
}

##################################################
## Import design matrix and prepare count table
##################################################
design = read.xlsx(design.file, sheet = 1) # all samples included in this file

colnames(design)[1] = 'SampleID'
design$conds = paste0(design$conditions, '_', design$SampleID,  '.batch', design$batch)

# prepare the data table for different batches
dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/'
batch1 = read.delim(paste0(dataDir,  'R10724_rnaseq/nf_out_RNAseq/featureCounts/merged_gene_counts.txt'), sep = '\t',  header = TRUE)
batch2 = read.delim(paste0(dataDir, 'Rxxxx_rnaseq_old/nf_out_RNAseq/featureCounts/merged_gene_counts.txt'), sep = '\t',  header = TRUE)
batch4 = read.delim(paste0(dataDir, 'R161513_rnaseq/nf_out_RNAseq/featureCounts/merged_gene_counts.txt'), sep = '\t', header = TRUE)

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

xx4 = process.countTable(all=batch4, design = design[which(design$batch == '4'), c(1:2)], merge.technicalRep.sameID = FALSE,
                         ensToGeneSymbol = FALSE)

colnames(xx4)[-1] = paste0(colnames(xx4)[-1], '.batch4')

all <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(xx1, xx2, xx3, xx4))

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
stat2 = read.csv(file = '../results/rnaseq_Rxxxx.old_R10724_R161513/QCs_stats.csv')

xx = data.frame(design, matrix(NA, nrow = nrow(design), ncol = 9), stringsAsFactors = FALSE)
colnames(xx)[7:15] = colnames(stat0)[-c(1:3)]

mm = match(stat0$sampleID, design$SampleID)

xx[mm, c(7:15)] = stat0[, -c(1:3)]

for(n in 1:nrow(xx))
{
  if(is.na(xx$pct.duplication[n])){
    # batch 1
    if(xx$batch[n] == 1) {
      kk = which(stat1$sampleID == xx$SampleID[n])
      kk = kk[grep('HLVGMDRXX_', stat1$sample[kk], invert = TRUE)]
      if(length(kk) != 1) cat(n, ' -- Error \n')
      xx[n, c(7:15)] = stat1[kk, c(4:12)]
    }
    # batch 2
    if(xx$batch[n] == 2) {
      kk = which(stat1$sampleID == xx$SampleID[n])
      kk = kk[grep('HLVGMDRXX_', stat1$sample[kk])]
      if(length(kk) != 2) cat(n, ' -- Error \n')
      xx[n, c(7:9, 12)] = apply(stat1[kk, c(4:6, 9)], 2, mean)
      xx[n, c(10:11, 13:15)] = apply(stat1[kk, c(7:8, 10:12)], 2, sum)
    }
    # batch 4
    if(xx$batch[n] == 4){
     kk = which(stat2$sampleID == xx$SampleID[n])
     if(length(kk) != 1) cat(n, ' -- Error \n')
     xx[n, c(7:15)] = stat2[kk, c(4:12)]
    }
  }
}

xx = data.frame(xx, stringsAsFactors = FALSE)

xx$pct.assigned.to.features = xx$assigned.reads/xx$total.reads *100
xx$alignment.uniq.rate = xx$unique.aligned/xx$trimmed.reads *100
xx = xx[, c(1:11, 13:15, 12, 17, 16)]

design = xx

save(design, all, file=paste0(RdataDir, 'Design_stats_readCounts_', version.analysis, '.Rdata'))

########################################################
########################################################
# Section II : Analyze the RNA-seq data 
# 
########################################################
########################################################
load(file=paste0(RdataDir, 'Design_stats_readCounts_', version.analysis, '.Rdata'))
design$batch = as.factor(design$batch)
annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

tfs = readRDS(file = paste0('../results/motif_analysis/TFs_annot/curated_human_TFs_Lambert.rds'))

Refine.ax.gene.annot = FALSE
if(Refine.ax.gene.annot){
  annot$gene.symbol.toUse = NA
  
  length(which(is.na(annot$gene.symbol.toUse)))
  kk = which(!is.na(annot$gene.evidence.synteny)) # synteny evidence
  annot$gene.symbol.toUse[kk] = annot$gene.evidence.synteny[kk]
  length(which(is.na(annot$gene.symbol.toUse)))
  
  kk = which(is.na(annot$gene.symbol.toUse) & !is.na(annot$gene.evidence.same.gene.symbol.hs.nr)) # shared names by hs and nr
  annot$gene.symbol.toUse[kk] = annot$gene.evidence.same.gene.symbol.hs.nr[kk]
  length(which(is.na(annot$gene.symbol.toUse)))
  
  kk = which(is.na(annot$gene.symbol.toUse) & !is.na(annot$gene.symbol.hs)) # with hs name
  annot$gene.symbol.toUse[kk] = annot$gene.symbol.hs[kk]
  length(which(is.na(annot$gene.symbol.toUse)))
  
  kk = which(is.na(annot$gene.symbol.toUse) & !is.na(annot$gene.symbol.nr)) # with hs name
  annot$gene.symbol.toUse[kk] = annot$gene.symbol.nr[kk]
  length(which(is.na(annot$gene.symbol.toUse)))
  
  saveRDS(annot, file = paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                               'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
}

##########################################
# convert gene names to gene symbols
##########################################
mm = match(all$gene, annot$geneID)
ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
all$gene[!is.na(mm)] = ggs[!is.na(mm)]

Select.genes.having.symbols = FALSE
if(Select.genes.having.symbols){
  gene.mapping = gene.mapping[which(!is.na(gene.mapping$gene.symbol.nr) | !is.na(gene.mapping$gene.symbol.hs)), ]
  all = all[!is.na(match(all$gene, gene.mapping$gene.id)), ]
}

##########################################
# general QC for RNA-seq
##########################################
QC.for.cpm = FALSE
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
raw = as.matrix(all[, -1])
rownames(raw) = all$gene

# select samples 
sels = which(design$batch != 1 & !(design$SampleID == '136150' & design$batch == 3))
design.matrix = design[sels, ]

raw = raw[, sels]

rm(design)

dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ conditions)

ss = rowSums(counts(dds))

length(which(ss > quantile(ss, probs = 0.85)))

dd0 = dds[ss > quantile(ss, probs = 0.85) , ]
dd0 = estimateSizeFactors(dd0)
sizefactors.UQ = sizeFactors(dd0)

jj = c(1:length(sizefactors.UQ))
jj = which(design.matrix$batch == 4)

plot(sizefactors.UQ[jj], colSums(counts(dds))[jj], log = 'xy')
text(sizefactors.UQ[jj], colSums(counts(dds))[jj], colnames(dd0), cex =0.6)

design.matrix$sizefactor = sizefactors.UQ

hist(log10(ss), breaks = 200, main = 'log2(sum of reads for each gene)')

cutoff.gene = 100
cat(length(which(ss > cutoff.gene)), 'genes selected \n')

dds <- dds[ss > cutoff.gene, ]
#design.matrix = design.matrix[with(design.matrix, order(conditions, SampleID)), ]

save.scalingFactors.for.deeptools = FALSE
if(save.scalingFactors.for.deeptools){
  ss = colSums(counts(dds))
  plot(ss[jj], (design.matrix$alignment.rate*design.matrix$trimmed.reads)[jj])
  
  reads.mapped = design.matrix$trimmed.reads*design.matrix$alignment.rate/100
  xx = data.frame(sampleID = design.matrix$SampleID,  
                  scalingFactor = reads.mapped/(design.matrix$sizefactor*median(reads.mapped)),
                  stringsAsFactors = FALSE)
  xx = xx[jj,]
  
  write.table(xx, file = paste0(resDir, '/DESeq2_scalingFactor_forDeeptools.txt'), sep = '\t',
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}


# normalization and dimensionality reduction
sizeFactors(dds) = sizefactors.UQ

fpm = fpm(dds, robust = TRUE)

save(dds, design.matrix, file = paste0(RdataDir, 'RNAseq_design_dds.object.Rdata'))
save(fpm, design.matrix, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))

kk = intersect(which(design.matrix$batch == 4), grep('Mature_LA|Mature_Hand', design.matrix$conditions))
plot.pair.comparison.plot(fpm[, kk[order(design.matrix$conditions[kk])]])


ii = grep('HOXA13|MEIS2|SOX9', rownames(fpm))
par(mfrow = c(2, 2))
plot(fpm[, c(21, 28)], log='xy', cex = 0.5, main = 'LA before correction');
points(fpm[ii, 21], fpm[ii, 28], cex = 2, col = 'red', pch = 16)
text(fpm[ii, 21], fpm[ii, 28], c('Meis2', 'HOXA13', 'SOX9'), col = 'red', pos = 4, offset = 1, cex = 1.5)
abline(0, 1, col='blue', lwd = 1.5)

plot(fpm[, c(29, 32)], log='xy', cex = 0.5, main = 'Hand before correction');
points(fpm[ii, 29], fpm[ii, 32], cex = 2, col = 'red', pch = 16)
text(fpm[ii, 29], fpm[ii, 32], c('Meis2', 'HOXA13', 'SOX9'), col = 'red', pos = 4, offset = 1, cex = 1.5)
abline(0, 1, col='blue', lwd = 1.5)

plot(fpm[, c(28, 32)], log='xy', cex = 0.5, main = 'corrected LA');
points(fpm[ii, 28], fpm[ii, 32], cex = 2, col = 'red', pch = 16)
text(fpm[ii, 28], fpm[ii, 32], c('Meis2', 'HOXA13', 'SOX9'), col = 'red', pos = 4, offset = 1, cex = 1.5)
abline(0, 1, col='blue', lwd = 1.5)

plot(fpm[, c(21, 29)], log='xy', cex = 0.5, main = 'corrected Hand');
points(fpm[ii, 21], fpm[ii, 29], cex = 2, col = 'red', pch = 16)
text(fpm[ii, 21], fpm[ii, 29], c('Meis2', 'HOXA13', 'SOX9'), col = 'red', pos = 4, offset = 1, cex = 1.5)
abline(0, 1, col='blue', lwd = 1.5)


vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = TRUE))
ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_batch2.3.4.pdf"), width=12, height = 8)

ggp = ggplot(data=pca2save[which(pca2save$batch==4), ], aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_batch4.pdf"), width=12, height = 8)


##########################################
# try to correct batches 
##########################################
require("sva")
sels = c(1:nrow(design.matrix))
cpm = log2(fpm[, sels] + 2^-6)

bc = droplevels(design.matrix$batch[sels])
#bc = levelsdroplevels(bc)
mod = model.matrix(~ as.factor(conditions), data = design.matrix[sels, ])

cpm.bc = ComBat(dat=cpm, batch=bc, mod=mod, par.prior=TRUE, ref.batch = 4, mean.only = FALSE)    

ntop = 500
library(factoextra)
xx = as.matrix(cpm.bc)
vars = apply(xx, 1, var)
xx = xx[order(-vars), ]
xx = xx[1:ntop, ]

#res.pca <- prcomp(t(xx[ ,grep('Mature', colnames(xx))]), scale = TRUE)
res.pca <- prcomp(t(xx), scale = TRUE)
#res.var <- get_pca_var(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
)


Test.glmpca.mds = FALSE
if(Test.glmpca.mds){
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
  
}

########################################################
########################################################
# Section III: check the differential expressed genes, 
# mainly TFs that are changed in mature UA, LA and Hand (head as a control)
# TFs that are changes in regeneration mature UA, BL.day5, BL.day9, BL.day13.proximal (Development as control)
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'RNAseq_design_dds.object.Rdata'))
dds0 = dds
fpm0 = fpm(dds)

# select mature samples
sels = intersect(which(design.matrix$batch == 4), grep('Mature', design.matrix$conditions))
dds = dds[, sels]
cpm = log2(fpm0[, sels] + 2^-6)

dds$conditions = droplevels(dds$conditions)
dds <- estimateDispersions(dds)

plotDispEsts(dds, ymin = 10^-3)

dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

res.ii = results(dds, contrast=c("conditions", 'Mature_Hand', 'Mature_UA'))
colnames(res.ii) = paste0(colnames(res.ii), "_Hand.vs.UA")
res = data.frame(res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("conditions", 'Mature_Hand', 'Mature_LA'))
colnames(res.ii) = paste0(colnames(res.ii), "_Hand.vs.LA")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("conditions", 'Mature_LA', 'Mature_UA'))
colnames(res.ii) = paste0(colnames(res.ii), "_LA.vs.UA")
res = data.frame(res, res.ii[, c(2, 5, 6)])


#dds <- nbinomLRT(dds, reduced = ~1 )
#res0 <- results(dds)

library("pheatmap")
pval.cutoff = 0.01
select = which(res$pvalue_Hand.vs.LA < pval.cutoff | res$pvalue_Hand.vs.UA < pval.cutoff | res$pvalue_LA.vs.UA < pval.cutoff)

df <- as.data.frame(colData(dds)[,c("conditions", 'batch')])
o1 = c(grep('UA', df$conditions), grep('LA', df$conditions), grep('Hand', df$conditions), grep('Head', df$conditions))

yy = cpm[select, o1]
ss = apply(as.matrix(yy), 1, mean)
yy = yy[which(ss>-2), ]

pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, show_colnames = FALSE,
         scale = 'row',
         cluster_cols=FALSE, annotation_col=df[o1, ])


#xx = res[select, ]
#xx = xx[order(xx$pvalue), ]
#xx = cpm
ggs = rownames(yy)
ggs = sapply(ggs, function(x) unlist(strsplit(as.character(x), '_'))[1])

mm = match(tfs[, 3], ggs)
xx = yy[mm[!is.na(mm)], ]

tf.sels = match(rownames(xx), rownames(res))
df <- as.data.frame(colData(dds)[,c("conditions", 'batch')])
o1 = c(grep('UA', df$conditions), grep('LA', df$conditions), grep('Hand', df$conditions), grep('Head', df$conditions))

yy = cpm[tf.sels, o1]
#ss = apply(as.matrix(yy), 1, mean)
#yy = yy[which(ss>-2), ]

pheatmap(yy, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE,
         scale = 'row',
         cluster_cols=FALSE, annotation_col=df[o1, ], fontsize_row = 8)


##########################################
# dynamic TFs in regeneration (development stages as contorls) 
##########################################
load(file = paste0(RdataDir, 'RNAseq_design_dds.object.Rdata'))
dds0 = dds
fpm0 = fpm(dds)

# select mature samples
sels = unique(c(which(design.matrix$batch == 3 & design.matrix$conditions != 'BL_UA_13days_distal'), 
             which(design.matrix$conditions == 'Embryo_Stage46_proximal')))
dds = dds[, sels]
cpm = log2(fpm0[, sels] + 2^-6)

dds$conditions = droplevels(dds$conditions)
dds <- estimateDispersions(dds)

plotDispEsts(dds, ymin = 10^-3)

dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

res.ii = results(dds, contrast=c("conditions", 'BL_UA_13days_proximal', 'Mature_UA'))
colnames(res.ii) = paste0(colnames(res.ii), "_BL.D13.vs.UA")
res = data.frame(res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("conditions", 'BL_UA_5days', 'Mature_UA'))
colnames(res.ii) = paste0(colnames(res.ii), "_BL.D5.vs.LA")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("conditions", 'BL_UA_9days', 'Mature_UA'))
colnames(res.ii) = paste0(colnames(res.ii), "_BL.D9.vs.UA")
res = data.frame(res, res.ii[, c(2, 5, 6)])


#dds <- nbinomLRT(dds, reduced = ~1 )
#res0 <- results(dds)

library("pheatmap")
pval.cutoff = 0.01
select = which(res$pvalue_BL.D13.vs.UA < pval.cutoff | res$pvalue_BL.D5.vs.LA < pval.cutoff | res$pvalue_BL.D9.vs.UA < pval.cutoff)

df <- as.data.frame(colData(dds)[,c("conditions", 'batch')])
o1 = c(grep('Mature_UA', df$conditions), grep('BL_UA_5days', df$conditions), grep('BL_UA_9days', df$conditions), 
       grep('BL_UA_13days_proximal', df$conditions), grep('Embryo_Stage40', df$conditions), 
       grep('Embryo_Stage46_proximal', df$conditions))

yy = cpm[select, o1]
ss = apply(as.matrix(yy), 1, mean)
yy = yy[which(ss>-2), ]

pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, show_colnames = FALSE,
         scale = 'row',
         cluster_cols=FALSE, annotation_col=df[o1, ])


#xx = res[select, ]
#xx = xx[order(xx$pvalue), ]
#xx = cpm
ggs = rownames(yy)
ggs = sapply(ggs, function(x) unlist(strsplit(as.character(x), '_'))[1])

mm = match(tfs[, 3], ggs)
xx = yy[mm[!is.na(mm)], ]

tf.sels = match(rownames(xx), rownames(res))
df <- as.data.frame(colData(dds)[,c("conditions", 'batch')])

yy = cpm[tf.sels, o1]
vars = apply(yy[, c(1:7)],1, var)
ss = apply(yy[, c(1:7)], 1, mean)
yy = yy[which(vars>=1.5), ]
#ss = apply(as.matrix(yy), 1, mean)
#yy = yy[which(ss>-2), ]

pheatmap(yy, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE,
         scale = 'row',
         cluster_cols=FALSE, annotation_col=df[o1, ], fontsize_row = 8, 
         filename = paste0(resDir, '/heatmap_DE_TFs_mUA_regeneration.pdf'),
         width = 12, height = 20)

write.table(yy, file = paste0(resDir, '/DEtfs_mUA_regeneration_dev.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)



