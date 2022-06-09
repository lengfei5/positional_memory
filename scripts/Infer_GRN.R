##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: infer axolotl limb regeneration GRN with bulk-ATAC, bulk-RNA and single cell RNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jun  7 10:46:45 2022
##########################################################################
##########################################################################
rm(list=ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('functions_chipSeq.R')
source('Functions_atac.R')

version.analysis = 'Rxxxx_R10723_R11637_R12810_atac'
#peakDir = "Peaks/macs2_broad"

resDir = paste0("../results/", version.analysis)

RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/'
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0('/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/SupTables/')

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

source('Functions_histM.R')
library(dplyr)
library(Seurat)
library(patchwork)

########################################################
########################################################
# Section : define gene set of interest
# mainly TFs here from bulk RNA-seq and limb-development genes 
########################################################
########################################################
##########################################
# collect TFs of interest
##########################################
# TF annotation from human
tfs = readRDS(file = paste0('../results/motif_analysis/TFs_annot/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)

# limb dev TFs
limb.genes = read.delim(file = '../data/limb_development_GO_term_summary_20220218_162116.txt', skip = 1, header = FALSE)
limb.genes = unique(limb.genes$V2)
limb.genes = toupper(as.character(limb.genes))
limb.genes = limb.genes[which(!is.na(match(limb.genes, tfs)))]
tfs.limb = limb.genes;
rm(limb.genes)

# MARA top TFs
bb = readRDS(file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata', '/MARA_output_top135_maxZscore.rds'))
tfs.mara = unique(bb$gene)
tfs.mara = unique(unlist(lapply(tfs.mara, function(x) unlist(strsplit(as.character(x), '_')))))

## postional TFs from microarray
tfs.positional = readRDS(file = paste0(RdataDir, '/postional_genes_tfs.rds'))

# DE tfs in regeneration samples of smart-seq2 data
rna = readRDS(file = paste0("../results/RNAseq_data_used/Rdata/", 
                            'regeneration_dynamicGeneClusters_allGenes.rds'))
rna$gene = get_geneName(rownames(rna))

tfs.limb[which(is.na(match(tfs.limb, rna$gene)))]
tfs.mara[which(is.na(match(tfs.mara, rna$gene)))]
tfs.positional[which(is.na(match(tfs.positional, rna$gene)))]

rna = rna[which(!is.na(rna$cluster)), ]
rna = rna[which(!is.na(match(rna$gene, tfs))), ]

tfs.rna = unique(rna$gene)

tfs.limb[which(is.na(match(tfs.limb, tfs.rna)))]
tfs.mara[which(is.na(match(tfs.mara, tfs.rna)))]
tfs.positional[which(is.na(match(tfs.positional, tfs.rna)))]

##########################################
# import scRNA-seq data for the dataset to used in GENIE3
##########################################
aa = readRDS(file = paste0(RdataDir, '/Gerber2018_Fluidigm.C1_batches_seuratObj.rds'))
aa = subset(aa, cells = colnames(aa)[which(aa$timepoint != 'Stage28' & aa$timepoint != 'Stage40' & 
                                             aa$timepoint != 'Stage44')])
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(aa)
aa <- ScaleData(aa, features = all.genes)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
ElbowPlot(aa)

aa <- FindNeighbors(aa, dims = 1:15)
aa <- FindClusters(aa, resolution = 0.5)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, reduction = "umap", group.by = 'timepoint')

# add gene names 
annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
mm = match(rownames(aa), annot$geneID)
ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
ggs[is.na(mm)] = rownames(aa)[is.na(mm)]

E = aa@assays$RNA@data
rownames(E) = ggs
genes = get_geneName(rownames(E))

tfs.limb[which(is.na(match(tfs.limb, genes)))]
tfs.mara[which(is.na(match(tfs.mara, genes)))]
tfs.positional[which(is.na(match(tfs.positional, genes)))]
tfs.rna[which(is.na(match(tfs.rna, genes)))]

### gene set of interest
geneSet = unique(c(tfs.limb, tfs.mara, tfs.positional,  tfs.rna))
cat('total TFs of interest : ', length(geneSet), '\n')

mm = match(genes, geneSet)
E = E[!is.na(mm), ]  # many gene IDs share the same gene symbols in the annotation
genes = get_geneName(rownames(E)) 
gg.counts = table(genes)
cat(length(gg.counts), ' TFs was found in scRNA-seq data\n')

cat(names(gg.counts)[which(gg.counts>=4)])

# remove the TFs with more than 5 IDs in the annotations, mainly zfinger proteins
mm = match(genes, names(gg.counts)[which(gg.counts<4)])
E = E[!is.na(mm), ]  
genes = get_geneName(rownames(E)) 
gg.counts = table(genes)

## 346 TFs from scRNA-seq data with 289 unique gene symbols
saveRDS(E, file = paste0(RdataDir, '/GRNinference_scExprMatrix_346geneID_289TFsymbols.rds')) 

########################################################
########################################################
# Section : Prior network from bulk atac-seq
# 
########################################################
########################################################

##########################################
# collect CRE regions for selected TFs
##########################################
Prepare.enhancer.tss.4fimo.scanning = FALSE
if(Prepare.enhancer.tss.4fimo.scanning)
{
  E = readRDS(file = paste0(RdataDir, '/GRNinference_scExprMatrix_346geneID_289TFsymbols.rds'))
  geneID = get_geneID(rownames(E))
  
  # tss used in the analysis
  tss = readRDS(file =paste0(RdataDir, '/regeneration_matureSamples_tss_perGene_smartseq2_atac_histM_v5.rds'))
  tss = data.frame(geneID = tss$geneID, coords = tss$coords)
  missed = which(is.na(match(geneID, tss$geneID)))
  cat(length(missed), ' TFs miss TSS \n')
  
  tss2 = readRDS(file = paste0(RdataDir, '/missed_tss_GRN_inference.rds'))
  tss = rbind(tss, tss2)
  
  missed = which(is.na(match(geneID, tss$geneID)))
  cat(length(missed), ' TFs miss TSS \n')
  
  Add.missed.TSS = FALSE
  if(Add.missed.TSS){
    gtf.all = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/AmexT_v47_Hox.patch.gtf'
    amex.all = GenomicFeatures::makeTxDbFromGFF(file = gtf.all)
    tss = GenomicFeatures::promoters(amex.all, upstream = 2000, downstream = 2000, use.names = TRUE)
    
    tss = as.data.frame(tss)
    tss$tx_name = sapply(tss$tx_name, function(x){x = unlist(strsplit(as.character(x), '[|]')); x[length(x)]})
    tss$start[which(tss$start<=1)] = 1
    
    annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
    'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
    tss$geneID = annot$geneID[match(tss$tx_name, annot$transcriptID)]
    
    kk = match(geneID[missed], tss$geneID)
    
    xx = data.frame(geneID = tss$geneID[kk], coords = paste0(tss$seqnames[kk], ':', tss$start[kk], '-', tss$end[kk]))
    saveRDS(xx, file = paste0(RdataDir, '/missed_tss_GRN_inference.rds'))
    
  }
  
  ## resize the promoters to -2kb--+300bp
  tp = data.frame(t(sapply(tss$coords, 
                           function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  
  tp$X3 = as.numeric(as.character(tp$X3)) - 2000 +300
  #tp$X2 = as.numeric(as.character(tp$X2))
  rownames(tp) = rownames(tss)
  tp$name = rownames(tss)
  tp$strand = '*'
  tp$X3 = as.factor(tp$X3)
  
  # atac-seq peaks
  enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual.rds'))
  
  missed = which(is.na(match(geneID, enhancers$targets)))
  cat(length(missed), ' TFs miss atacPeaks \n')
  
  pp = rownames(enhancers)
  pp = data.frame(t(sapply(pp, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp$name = rownames(pp)
  pp$strand = '*'
  
  pp = rbind(pp, tp)
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  saveDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/motif_analysis/peaks/'
  rtracklayer::export(pp, format = 'bed', con = paste0(saveDir, 'atacPeaks_2kbtss_for_fimo.bed'))
  
}

##########################################
# process matrix for gene-regulator matrix 
##########################################
## CRE-motif-occurence matrix
moc1 = readRDS(file = '../results/motif_analysis/motif_oc_fimo_jaspar2022_pval.0.0001_v1.rds')
moc2 = readRDS(file = '../results/motif_analysis/motif_oc_fimo_jaspar2022_pval.0.0001_regenerationTSS.rds')
cres = c(rownames(moc1), rownames(moc2))

mm = match(colnames(moc1), colnames(moc2))
moc2 = moc2[,mm]

mocs = rbind(moc1, moc2)
rm(moc1); rm(moc2)

# collect CREs for targets defined previously = rownames(E)
# target-to-CRE matrix
targets.ids = get_geneID(rownames(E))
targets = data.frame(geneID = targets.ids, stringsAsFactors = FALSE)
targets$CREs = NA

bed = read.table(file = 
    paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/motif_analysis/FIMO_promoters/', 
           'promoter_tfs/sorted.bed'), header = FALSE)
kk = match(targets$geneID, bed$V4)
mm = which(!is.na(kk))
targets$CREs[mm] = paste0(bed$V1[kk[mm]], ':', bed$V2[kk[mm]], '-', bed$V3[kk[mm]])

enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual.rds'))

mm = match(targets.ids, enhancers$targets)
ii = which(!is.na(mm))
mm = mm[ii]
targets2 = data.frame(geneID = enhancers$targets[mm], CREs = rownames(enhancers)[mm], stringsAsFactors = FALSE)
targets = rbind(targets, targets2)

xx = table(targets$geneID, targets$CREs)

ss = apply(xx, 1, sum)

xx = xx[which(ss>0), ]
mm = match(rownames(xx), targets.ids)
rownames(xx) = rownames(E)[mm]

jj = match(colnames(xx), rownames(mocs))
ii = which(!is.na(jj))
jj = jj[ii]

xx = xx[,ii]
mocs = mocs[jj, ]

target.moc = xx %*% mocs
ss = apply(target.moc, 2, sum)


## motif to TF association matrix
mapping =readRDS(file = 
      paste0('../data/JASPAR2022_CORE_vertebrates_nonRedundant_metadata.rds')) # association between motifs and TF symbols
mapping = mapping[, c(3, 2)]

atm = matrix(0, nrow = ncol(target.moc), ncol = nrow(target.moc))
colnames(atm) = rownames(target.moc)
rownames(atm) = colnames(target.moc)
genes = get_geneName(colnames(atm))

for(n in 1:nrow(atm))
{
  # n = 164
  tfs = unlist(strsplit(as.character(mapping$tfs[which(mapping$name == rownames(atm)[n])]), '_'))
  cat(n, '--', rownames(atm)[n], ' -- tfs --', tfs, '\n')
  
  for(tf in tfs)
  {
    jj = which(genes == tf)
    atm[n, jj] = 1
  }
}

target.tfs = target.moc %*% atm

x = target.tfs[grep('HOXA13', rownames(target.tfs)), ];
x[which(x>0)]

ss = apply(target.tfs, 2, sum) # due to missing motif for associated TFs
target.tfs = target.tfs[, which(ss>0)]

ss = apply(target.tfs>0, 1, sum)
target.tfs = target.tfs[which(ss>0), ]

saveRDS(target.tfs, file = paste0(RdataDir, '/GRN_priorNetwork_target_TFs.rds'))

########################################################
########################################################
# Section : GENIE3 
# 
########################################################
########################################################
source('myGENIE3.R')

E = readRDS(file = paste0(RdataDir, '/GRNinference_scExprMatrix.rds'))
target.tfs = readRDS(file = paste0(RdataDir, '/GRN_priorNetwork_target_TFs.rds'))

mm = match(rownames(target.tfs), rownames(E))
E = E[mm, ]
E = as.matrix(E)

sds = apply(E, 1, sd)
sels = which(sds>10^-4)
E = E[sels, ]
target.tfs = target.tfs[sels, ]

#tic()
#weight.matrix1 = GENIE3(expr.matrix = E, priorRegulators =  target.tfs)
#toc()

tic()
wtm = GENIE3(expr.matrix = E, priorRegulators =  target.tfs, ncore = 8)
saveRDS(wtm, file = paste0(RdataDir, '/first_test_Genie3.rds'))

toc()


wtm = readRDS(file =  paste0(RdataDir, '/first_test_Genie3.rds'))
ggs = colnames(wtm)
genes = get_geneName(ggs)
gene.counts = table(genes)
gene.unique = names(gene.counts)[which(gene.counts==1)]
mm = match(genes, gene.unique)
kk = which(!is.na(mm))
ggs[kk] = genes[kk]
colnames(wtm) = ggs
rownames(wtm) = ggs

cutoff = 0.01;
ss1 = apply(wtm, 1, function(x) length(which(x>cutoff)))
ss2 = apply(wtm, 2, function(x) length(which(x>cutoff)))
sels = which(ss1>=5 & ss2>=5)
wtm = wtm[sels, sels]

link.list = get.link.list(wtm, threshold = cutoff)
dim(link.list)
head(link.list)


