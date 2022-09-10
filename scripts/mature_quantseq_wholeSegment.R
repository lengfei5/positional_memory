##########################################################################
##########################################################################
# Project: positional memory 
# Script purpose : analyze whole-segment quant-seq data 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Sep 10 11:22:32 2022
##########################################################################
##########################################################################
rm(list = ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('Functions_rnaseq.R')
require(openxlsx)
require(ggplot2)
require(DESeq2)
require(dplyr)
require(gridExtra)

version.Data = 'QuantSeq_noSortedCells';

# in the version 20220408, counts were done with update annot with Hox patch
version.analysis = paste0("_", version.Data, "_20220910") 
## Directories to save results 
dataDir = "/Volumes//groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R14026_Quantseq"

design.file = paste0(dataDir, '/paramFile.csv')

resDir = paste0("../results/", version.Data)
#tabDir =  paste0(resDir, "/tables/")
tfDir = '~/workspace/imp/positional_memory/results/motif_analysis'
RdataDir = paste0(resDir, "/Rdata/")
#shareDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/AkaneToJingkuiShareFiles/results_rnaseq/positional_genes'
figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
#tableDir = paste0(figureDir, 'tables4plots/')
tableDir = paste0('/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/SupTables/')

annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))


tfs = readRDS(file = paste0('../results/motif_analysis/TFs_annot/curated_human_TFs_Lambert.rds'))
sps = readRDS(file = '~/workspace/imp/organoid_patterning/results/Rdata/curated_signaling.pathways_gene.list_v2.rds')
eps = readRDS(file = paste0('../data/human_chromatin_remodelers_Epifactors.database.rds'))
rbp = readRDS(file = paste0('../data/human_RBPs_rbpdb.rds'))
tfs = unique(tfs$`HGNC symbol`)
sps = toupper(unique(sps$gene))
sps = setdiff(sps, tfs)

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tableDir)){dir.create(tableDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}


########################################################
########################################################
# Section : process count tables 
# 
########################################################
########################################################
design = read.csv(design.file)

xlist = list.files(path=paste0(dataDir, '/featurecounts_Q10'),
                   pattern = "*featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

design = design[, c(3,2,1)]

colnames(design)[1] = 'SampleID'
counts = process.countTable(all=all, design = design[, c(1,2)], merge.technicalRep.sameID = FALSE)


all = counts

mm = match(all$gene, annot$geneID)
ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
all$gene[!is.na(mm)] = ggs[!is.na(mm)]

raw = as.matrix(all[, -1])
rownames(raw) = all$gene
#design$batch = paste0(design$request, '_', design$protocol)
colnames(design)[2:3] = c('condition', 'batch')
dds <- DESeqDataSetFromMatrix(raw, DataFrame(design), design = ~ condition)

save(design, dds, file = paste0(RdataDir, 'design_dds_wholeSegements_Quantseq.Rdata'))


##########################################
# filering and normalization 
##########################################
load(file = paste0(RdataDir, 'design_dds_wholeSegements_Quantseq.Rdata'))

table(design$condition, design$batch)
design$condition = paste0('Mature_', design$condition)
counts = counts(dds)
colnames(counts) = paste0('Mature_', colnames(counts))

dds <- DESeqDataSetFromMatrix(counts, DataFrame(design), design = ~ condition)

#design$protocol = gsub(' ', '', design$protocol)
#design$batch = paste0(design$request, '_', design$protocol)

#dds = DESeqDataSetFromMatrix(raw, DataFrame(design), design = ~ condition)
dds$condition = droplevels(dds$condition)

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100);abline(v = log10(20), lwd = 2.0, col = 'red')
cat(length(which(ss>10)), ' gene selected \n')

dds = dds[which(ss>10), ]

dds$condition = droplevels(dds$condition)

x <- DESeqDataSetFromMatrix(counts(dds), DataFrame(design), design = ~ batch + condition )
#dds$batch = droplevels(dds$batch)

dds = estimateSizeFactors(x)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE, ntop = 3000))
pca2save$condition = as.factor(pca2save$condition)
pca2save$batch = as.factor(pca2save$batch)
ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
  geom_point(size=4) + 
  geom_text(hjust = 1, nudge_y = 1, size=4)

ggsave(paste0(resDir, '/PCA_QuangSeq_wholeSegments.pdf'),  width=12, height = 8)

cpm = fpm(dds)
cpm = log2(cpm + 2^-4)

#plot.pair.comparison.plot(cpm[, c(1:7)], linear.scale = FALSE)

## batch correction with combat
require(sva)
tmm = cpm
bc = as.factor(design$batch)
mod = model.matrix(~ as.factor(condition), data = design)

# if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
# because there is no better batche other others 
#ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE)
#fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 

source('Functions_atac.R')
make.pca.plots(tmm, ntop = 3000, conds.plot = 'Mature')
ggsave(paste0(resDir, "/smartseq2_R10724_R11635_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)


make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'Mature')
ggsave(paste0(resDir, "/smartseq2_R10724_R11635_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)

plot.pair.comparison.plot(fpm.bc[, c(1:6)], linear.scale = FALSE)

#dds$condition = droplevels(dds$condition)
dds$condition <- relevel(dds$condition, ref = "Mature_UA")

#dds$condition = droplevels(dds$condition)
dds <- estimateDispersions(dds, fitType = 'parametric')
plotDispEsts(dds, ymin = 10^-3); abline(h = 0.1, col = 'blue', lwd = 2.0)

dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

res.ii = results(dds, contrast=c("condition", 'Mature_Hand', 'Mature_UA'), alpha = 0.1)
colnames(res.ii) = paste0(colnames(res.ii), "_mHand.vs.mUA")
res = data.frame(res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("condition", 'Mature_Hand', 'Mature_LA'), alpha = 0.1)
colnames(res.ii) = paste0(colnames(res.ii), "_mHand.vs.mLA")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("condition", 'Mature_LA', 'Mature_UA'), alpha = 0.1)
colnames(res.ii) = paste0(colnames(res.ii), "_mLA.vs.mUA")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res$gene = sapply(rownames(res), function(x) unlist(strsplit(as.character(x), '_'))[1])
res$geneID = sapply(rownames(res), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})

res = res[order(res$padj_mHand.vs.mUA), ]
fpm.bc = fpm.bc[match(rownames(res), rownames(fpm.bc)), ]

saveRDS(data.frame(fpm.bc, res, stringsAsFactors = FALSE), 
        file = paste0(RdataDir, "matureSamples_Quantseq_cpm_DEgenes.rds"))

write.csv(data.frame(fpm.bc, res, stringsAsFactors = FALSE), 
          file = paste0(tableDir, 'wholeSegments_matureSamples_Quantseq_cpm_DEgenes.csv'), 
          row.names = TRUE)

