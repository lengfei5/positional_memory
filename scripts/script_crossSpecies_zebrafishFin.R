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
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 

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
annot.file = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/zebrafish/GRCz11/annotation_ens_biomart.txt'
annot = read.delim(annot.file)

tfs = readRDS(file = paste0('../results/motif_analysis/TFs_annot/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)

########################################################
########################################################
# Section : RNA-seq analysis
# 
########################################################
########################################################
count.file = paste0(rnaDir, '/GSE126701_RNA_finRegen_featureCounts.txt')
counts = read.table(count.file, sep = '\t', header = TRUE)

mm = match(counts$Geneid, annot$Gene.stable.ID)
jj = which(annot$Human.gene.name[mm] != '' & !is.na(mm))
mm = mm[jj]
ggs = paste0(annot$Human.gene.name[mm], '_',  annot$Gene.stable.ID[mm])
rownames(counts) = counts$Geneid
rownames(counts)[jj] = ggs
counts = counts[, -1]

conds = colnames(counts)
conds = gsub('_rep1|_rep2', '', conds)
conds = data.frame(condition = conds)
dds <- DESeqDataSetFromMatrix(counts, DataFrame(conds), design = ~ condition)

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100);abline(v = log10(20), lwd = 2.0, col = 'red')
cat(length(which(ss>10)), ' gene selected \n')

dds = dds[which(ss>10), ]

dds = estimateSizeFactors(dds)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('condition'), returnData = FALSE)
print(pca)

##########################################
# DE test
##########################################
dds$condition = droplevels(dds$condition)
dds$condition <- relevel(dds$condition, ref = "sp7po0dpa")

dds <- DESeq(dds, test="Wald", fitType = c("parametric"))

plotDispEsts(dds, ymin = 10^-3)
resultsNames(dds)

fpm = log2(fpm(dds) + 2^-6)

res = results(dds, name="condition_sp7po4dpa_vs_sp7po0dpa", test = 'Wald')
res <- lfcShrink(dds, coef="condition_sp7po4dpa_vs_sp7po0dpa")
colnames(res) = paste0(colnames(res), "_sp7po_4dpa.vs.0dpa")
res = data.frame(res[, c(2, 5, 6)])

res.ii = results(dds, contrast = c('condition', 'sp7ne4dpa', 'sp7ne0dpa'), test = 'Wald')
res.ii <- lfcShrink(dds, contrast = c('condition', 'sp7ne4dpa', 'sp7ne0dpa'))
colnames(res.ii) = paste0(colnames(res.ii), "_sp7ne_4dpa.vs.0dpa")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res = data.frame(fpm, res, stringsAsFactors = FALSE)

saveRDS(res, file = paste0(RdataDir, '/RNAseq_fpm_DEgenes_lfcShrink_res.rds'))

##########################################
# select significant TFs  
##########################################
res = readRDS(file = paste0(RdataDir, '/RNAseq_fpm_DEgenes_lfcShrink_res.rds'))

fdr.cutoff = 0.05; logfc.cutoff = 1 # select only activated genes in regeneration
jj = which((res$padj_sp7ne_4dpa.vs.0dpa < fdr.cutoff & res$log2FoldChange_sp7ne_4dpa.vs.0dpa > logfc.cutoff) |
             (res$padj_sp7po_4dpa.vs.0dpa < fdr.cutoff & (res$log2FoldChange_sp7po_4dpa.vs.0dpa) > logfc.cutoff)
)
cat(length(jj), '\n')

res = res[jj, ]

saveRDS(res, file = paste0(RdataDir, '/RNAseq_fpm_DEgenes_lfcShrink_res_RRgenes.rds'))

ggs = get_geneName(rownames(res))
print(intersect(ggs, tfs))

mm = match(ggs, tfs)
yy = res[unique(c(which(!is.na(mm)))), ]
cat(nrow(yy), ' DE TFs found \n')
yy = yy[order(yy$padj_sp7ne_4dpa.vs.0dpa), ]

grep('RUNX', rownames(yy))
yy[grep('RUNX', rownames(yy)), ]
ggs = as.character(get_geneName(rownames(yy)))

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

# save sample info for footprinting 
saveSampleInfo4Footprint = FALSE
if(saveSampleInfo4Footprint){
  xx = design[, c(1, 3)]
  colnames(xx) = c('sampleID', 'condition')
  
  write.table(xx, file = paste0(dataDir, '/footprinting/sample_infos.txt'), sep = '\t', quote = FALSE, col.names = TRUE,
              row.names = FALSE)
  
}



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
# DE peak test
##########################################
dds$condition = droplevels(dds$condition)
dds$condition <- relevel(dds$condition, ref = "ATAC_sp7ne0dpa")

dds <- DESeq(dds, test="Wald", fitType = c("parametric"))

plotDispEsts(dds, ymin = 10^-3)
resultsNames(dds)

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

########################################################
########################################################
# Section : Motif analysis with MARA
# 
########################################################
########################################################
## process fimo output 
fimo.out = paste0(dataDir, '/FIMO_atacPeak/fimo_out/fimo.tsv')
source('Functions_MARA.R')

prefix = paste0(RdataDir, '/motif_oc_fimo_atacPeaks_jaspar2022.core')
make.motif.oc.matrix.from.fimo.output(fimo.out = fimo.out, 
                                      pval = c(10^-4, 10^-5, 10^-6), 
                                      prefix = prefix)

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

##########################################
# ## prepare the motif oc matrix and response matrix
##########################################
motif.oc = readRDS(file = paste0(RdataDir, '/motif_oc_fimo_atacPeaks_jaspar2022.core_pval_0.0001.rds'))
grep('RUNX', colnames(motif.oc))

res = readRDS(file = paste0(RdataDir, '/DE_atacPeaks_fpm_annotatation.rds'))

cat(nrow(res), 'DE peaks found !\n')

res = res[order(-res$log2FoldChange_sp7po_4dpa.vs.0dpa), ]

## prepare the response matrix
keep = as.matrix(res[, c(1:8)])
conds = unique(sapply(colnames(keep), function(x) unlist(strsplit(as.character(x), '_'))[2]))

## make sure response is logscale 
keep = log2(keep + 2^-4)

# select samples to use
#keep = keep[, sample.sels]
Y = cal_sample_means(keep, conds = conds)
rownames(Y) = sapply(rownames(Y), function(x) {x = unlist(strsplit(gsub(':', '-', as.character(x)), '-'));
paste0(x[1], ':', (as.numeric(x[2]) -1), '-', x[3])}) 
mm = match(rownames(Y), rownames(motif.oc))

cat(length(which(is.na(mm))), ' missed peaks in ater fimo scanning \n')
Y = Y[which(!is.na(mm)), ]
motif.oc = motif.oc[mm[which(!is.na(mm))], ]
grep('RUNX1', colnames(motif.oc))

### run MARA analysis
source('Functions_MARA.R')
aa1 = run.MARA.atac(motif.oc, Y[ ,c(3,4)],  method = 'Bayesian.ridge')
aa2 = run.MARA.atac(motif.oc, Y[ ,c(1,2)],  method = 'Bayesian.ridge')

aa = data.frame(aa1[, c(1:2)], aa2[, c(1:2)], aa1[, c(3:7)], stringsAsFactors = FALSE)
aa$combine.Zscore = apply(as.matrix(aa[, c(1:4)]), 1, function(x) sqrt(mean(x^2)))
aa$maxZscore = apply(as.matrix(aa[, c(1:4)]), 1, function(x){x.abs = abs(x); return(max(x.abs))})
aa$rank = order(aa$combine.Zscore)
aa = aa[order(-aa$combine.Zscore), ]
grep('RUNX', rownames(aa))

saveRDS(aa, file = paste0(RdataDir, '/MARA_output_sorted.rds'))

##########################################
# select motif of interest and add associated TF expressed 
##########################################
bb = readRDS(file = paste0(RdataDir, '/MARA_output_sorted.rds'))
kk = which(bb[, 2]> bb[, 1] | bb[, 4]>bb[, 3] ) ## motif activated in regeneration
bb = bb[kk, ]

kk = which(abs(bb[, 2]) >2| abs(bb[, 4])>2) # 
bb = bb[kk, ]

## import processed RNAseq data
cpm = readRDS(file = paste0(RdataDir, '/RNAseq_fpm_DEgenes_lfcShrink_res.rds'))
cpm = cpm[, c(1:8)]
cpm = cal_sample_means(cpm, conds = c("sp7po0dpa", "sp7po4dpa", "sp7ne0dpa", "sp7ne4dpa"))
ss = apply(cpm, 1, max)

ss[grep('RUNX', names(ss))]

cpm = cpm[which(ss > 0), ]
ggs = get_geneName(rownames(cpm))

## find expression of associated TFs
bb$index = NA
bb$corr.tfs = NA
for(n in 1:nrow(bb)){
  tfs = unique(unlist(strsplit(as.character(bb$gene[n]), '_')))
  if(length(tfs) == 1){
    if(tfs == 'SNAI3') tfs = 'SNAI2'
    if(tfs == 'ZNF701') tfs = 'ZNF419' 
  }
  
  jj = c()
  for(tf in tfs) jj = c(jj, which(ggs == tf))
  
  if(length(jj) == 1) {
    bb$index[n] = jj
    bb$corr.tfs[n] = cor(as.numeric(bb[n, c(1:4)]), as.numeric(cpm[jj, ]))
  }
  if(length(jj) >=2){
    correlation = cor(as.numeric(bb[n, c(1:4)]), t(cpm[jj,]), method = 'pearson')
    bb$index[n] = jj[which.max(abs(correlation))]
    bb$corr.tfs[n] = correlation[which.max(abs(correlation))]
  }
  if(length(jj) == 0){
    cat(n, '--', tfs, ': ')
    cat(length(jj), 'tfs found \n')
  }
}

bb$tfs = ggs[bb$index]
bb$tf_ids = rownames(cpm)[bb$index]

xx = bb[which(!is.na(bb$tfs)), ]

xx$tfs.speciees = annot$Gene.name[match(get_geneID(xx$tf_ids), annot$Gene.stable.ID)]
rnas = cpm[xx$index, ]
colnames(rnas) = paste0(colnames(rnas), '_rna')

xx = data.frame(xx, rnas, stringsAsFactors = FALSE)

saveRDS(xx, file =  paste0(RdataDir, '/MARA_output_Motifs_filteredTFsExpr.rds'))

##########################################
# plot the motif activity and associated TF expression
##########################################
xx = readRDS(file =  paste0(RdataDir, '/MARA_output_Motifs_filteredTFsExpr.rds'))
test = as.matrix(xx[, c(1:4)])
rownames(test) = xx$gene

cat('range of motif activity -- ', range(test), '\n')

range <- 6 ; breaks = 10
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))

df <- data.frame(colnames(test))
rownames(df) = colnames(test)
colnames(df) = 'samples'

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
cols <- colorRampPalette(cols)(breaks)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

plt = pheatmap(test, show_rownames=TRUE, show_colnames = FALSE, 
               color = cols,
               scale = 'none', cluster_cols=FALSE, main = '', 
               na_col = "white", #fontsize_row = 10, 
               breaks = seq(-range, range, length.out = breaks), 
               annotation_col = df, 
               cluster_rows=TRUE,
               treeheight_row = 20,
               #clustering_method = 'complete', cutree_rows = 6, 
               annotation_legend = TRUE,
               clustering_callback = callback,
               filename = paste0(figureDir, 'MARA_bayesianRidge_Jaspar2022_', species,'.pdf'), 
               width = 6, height = 10)

### plot associated TF expression
test = xx[plt$tree_row$order, c(14:18)]
test$tfs.speciees = as.character(test$tfs.speciees)

rownames(test) = test$tfs.speciees

grep('fosl2', test$tfs.speciees)
test$tfs.speciees[grep('fosl2', test$tfs.speciees)]

#test$tfs.speciees[41] = 'ctcf.2'
test$tfs.speciees[58] = 'fosl2.2'
rownames(test) = test$tfs.speciees

test = test[, -1]
test = t(apply(test, 1, cal_centering))

cat('range of centered TF expression --',  range(test), '\n')
range <- 2.5
test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
cols = rev(colour('PRGn')(5))
df <- data.frame(colnames(test))
rownames(df) = colnames(test)
colnames(df) = 'samples'

pheatmap(test, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, 
         color = cols,
         scale = 'none', cluster_cols=FALSE, main = '', 
         na_col = "white",
         #fontsize_row = 10, 
         #breaks = seq(-8, 8, length.out = 8), 
         annotation_col = df, 
         #treeheight_row = 20,
         #clustering_method = 'complete', cutree_rows = 6, 
         annotation_legend = TRUE,
         #clustering_callback = callback,
         filename = paste0(figureDir, 'TFs_associatedMotifs_regeneration_', species, '.pdf'), 
         width = 5, height = 10)


########################################################
########################################################
# Section : search for TF target with footprint analysis  
# 
########################################################
########################################################


motif = 'REL' 
for(m in 1:length(dir.list))
{
  # m = 1
  subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
  subdir = subdir[grep(motif, subdir)]
  subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
  subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
  subdir = subdir[grep('RELB', subdir, invert = TRUE)]
  
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

pp.annots = annotatePeak(footprint, TxDb=gtf, tssRegion = c(-2000, 2000), level = 'transcript')

## save the target genes
pp.annots = as.data.frame(pp.annots)
pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 2000), ]

pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
ggs = unique(as.character(ggs))
cat(length(ggs), ' targets \n')

saveRDS(ggs, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                           '/targetGenes_footprint_', motif, '_', species,  '.rds'))
