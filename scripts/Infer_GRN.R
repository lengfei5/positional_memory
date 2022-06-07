##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: infer axolotl limb regeneration GRN with bulk-ATAC, bulk-RNA and single cell RNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jun  7 10:46:45 2022
##########################################################################
##########################################################################
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

# smart-seq2 data
rna = readRDS(file = paste0("../results/RNAseq_data_used/Rdata/", 
                            'regeneration_dynamicGeneClusters_allGenes.rds'))
rna$gene = get_geneName(rownames(rna))

tfs.limb[which(is.na(match(tfs.limb, rna$gene)))]
tfs.mara[which(is.na(match(tfs.mara, rna$gene)))]

rna = rna[which(!is.na(rna$cluster)), ]
rna = rna[which(!is.na(match(rna$gene, tfs))), ]

tfs.rna = unique(rna$gene)

tfs.limb[which(is.na(match(tfs.limb, tfs.rna)))]
tfs.mara[which(is.na(match(tfs.mara, tfs.rna)))]

##########################################
# scRNA-seq data 
##########################################
aa = readRDS(file = paste0(RdataDir, '/Gerber2018_Fluidigm.C1_batches_seuratObj.rds'))
p1 = DimPlot(aa, reduction = "umap", group.by = 'timepoint')
#p2 = DimPlot(aa, reduction = "umap", group.by = 'batch')
p1

annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

mm = match(rownames(aa), annot$geneID)
ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
ggs[is.na(mm)] = rownames(aa)[is.na(mm)]
rownames(aa) = ggs

genes = get_geneName(ggs)

tfs.limb[which(is.na(match(tfs.limb, genes)))]
tfs.mara[which(is.na(match(tfs.mara, genes)))]
tfs.rna[which(is.na(match(tfs.rna, genes)))]


