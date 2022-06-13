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
  rownames(tss) = tss$geneID
  missed = which(is.na(match(geneID, tss$geneID)))
  cat(length(missed), ' TFs miss TSS \n')
  
  tss2 = readRDS(file = paste0(RdataDir, '/missed_tss_GRN_inference.rds'))
  rownames(tss2) = tss2$geneID
  
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
  enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual_targets.rds'))
  
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
  rtracklayer::export(pp, format = 'bed', con = paste0(saveDir, 'atacPeaks_2kbtss_for_fimo_v2.bed'))
  
}

##########################################
# process matrix for gene-regulator matrix 
##########################################
## manually correct motif-to-tf association 
Manual_correction_motif_tf_association = FALSE
if(Manual_correction_motif_tf_association){
  E = readRDS(file = paste0(RdataDir, '/GRNinference_scExprMatrix_346geneID_289TFsymbols.rds'))
  genes = unique(get_geneName(rownames(E)))
  
  # association between motifs and TF symbols from JASPAR core and unvalided
  
  mapping =readRDS(file = 
                     paste0('../data/JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata.rds')) 
  
  mapping = rbind(mapping, c('UN0262.1', 'SALL1', 'SALL1_UN0262.1'))
  mapping = rbind(mapping, c('UN0262.1', 'SALL3', 'SALL3_UN0262.1'))
  mapping = rbind(mapping, c('MA0835.2', 'BATF2', 'BATF2_MA0835.2'))
  mapping = rbind(mapping, c('MA1504.1', 'HOXC5', 'HOXC5_MA1504.1'))
  mapping = rbind(mapping, c('MA0103.3', 'ZEB2', NA))
  mapping = rbind(mapping, c('MA1122.1', 'TFDP2', NA))
  mapping = rbind(mapping, c('MA0597.2', 'THAP4', NA))
  mapping = rbind(mapping, c('MA1928.1', 'BNC1', NA))
  mapping = rbind(mapping, c('UN0130.1', 'IRX5', NA))
  mapping = rbind(mapping, c('MA0626.1', 'NPAS3', NA))
  mapping = rbind(mapping, c('UN0116.1', 'DPF3', NA))
  mapping = rbind(mapping, c('MA0833.2', 'ATF5', NA))
  mapping = rbind(mapping, c('MA1989.1', 'BCL11A', NA))
  mapping = rbind(mapping, c('MA0108.2', 'TBPL1', NA))
  mapping = rbind(mapping, c('MA0690.2', 'TBX22', NA))
  mapping = rbind(mapping, c('MA0597.2', 'THAP6', NA))
  mapping = rbind(mapping, c('MA1638.1', 'SCX', NA))
  mapping = rbind(mapping, c('MA1707.1', 'DMRT2', NA))
  
  tfs = c()
  for(n in 1:nrow(mapping))
  {
    tfs = unique(c(tfs, unlist(strsplit(as.character(mapping$tfs[n]), '_'))))
  }
  missed = which(is.na(match(genes, tfs)))
  genes[missed]
  
  saveRDS(mapping, file = '../data/JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata_manual.rds')
  
  ##########################################
  # prioritize the core motifs over unvalided ones
  ##########################################
  mapping = readRDS(file = '../data/JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata_manual.rds')
  jj = which(is.na(mapping$name))
  mapping$name[jj] = paste0(mapping$tfs[jj], '_', mapping$motifs[jj])
  
  kk = grep('^UN', mapping$motifs)
  
  idx = c()
  ggs = c()
  for(n in 1:nrow(mapping))
  {
    test = unlist(strsplit(as.character(mapping$tfs[n]), '_'))
    ggs = c(ggs, test)
    idx = c(idx, rep(n, length(test)))
  }
  
  mapping = data.frame(mapping[idx, ], gene = ggs, stringsAsFactors = FALSE)
  
  idx.rm = c()
  ggs = unique(mapping$gene)
  for(n in 1:length(ggs))
  {
    # n = 1
    jj = which(mapping$gene == ggs[n])
    motifs = mapping$motifs[jj]
    jj1 = jj[grep('MA', motifs)]
    jj2 = jj[grep('UN', motifs)]
    if(length(jj1)>0 & length(jj2)>0) {
      cat(n, '--', ggs[n], '\n')
      idx.rm = c(idx.rm, jj2)  
    }
  }
  
  mapping = mapping[-idx.rm, ]
  saveRDS(mapping, file = '../data/JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata_manual_rmRedundantUNVALIDED.rds')
    
}

### CRE-motif-occurence matrix
mocs = readRDS(file = '../results/motif_analysis/motif_oc_fimo_atacPeaks.2kbTSS_jaspar2022.core.unvalided_pval.0.00001_v1.rds')
cres = unique(rownames(mocs))

### target-to-CRE matrix
E = readRDS(file = paste0(RdataDir, '/GRNinference_scExprMatrix_346geneID_289TFsymbols.rds'))
targets.ids = get_geneID(rownames(E))

#targets = data.frame(geneID = targets.ids, stringsAsFactors = FALSE)
#targets$CREs = NA
bed = read.table(file = 
    paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/motif_analysis/', 
           'FIMO_atacPeak_tss_mediumQ_core.unvalided_64Gmem/peaks_for_fimo_sorted.bed'), header = FALSE)
kk = match(bed$V4, targets.ids)
mm = which(!is.na(kk))
targets = data.frame(geneID = bed$V4[mm], CREs = paste0(bed$V1[mm], ':', bed$V2[mm], '-', bed$V3[mm]), stringsAsFactors = FALSE)

enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual.rds'))

mm = match(enhancers$targets, targets.ids)
mm = which(!is.na(mm))

targets2 = data.frame(geneID = enhancers$targets[mm], CREs = rownames(enhancers)[mm], stringsAsFactors = FALSE)
targets2$CREs = sapply(targets2$CREs, function(x) {x = unlist(strsplit(gsub(':', '-', as.character(x)), '-'));
 paste0(x[1], ':', (as.numeric(x[2]) -1), '-', x[3])}) 
# shift one bp due to the bed format and fasta formation

targets = rbind(targets, targets2)
rm(targets2)

targets = targets[match(unique(targets$CREs), targets$CREs), ]
targets$geneID = as.character(droplevels(targets$geneID))
#targets$CREs =droplevels(targets$CREs)

xx = table(targets$geneID, targets$CREs) ## target-CRE matrix: each row is tf and each column is associated tss and enhancer peaks
ss = apply(xx, 1, sum)

xx = xx[which(ss>0), ]
mm = match(rownames(xx), targets.ids)
rownames(xx) = rownames(E)[mm]

## match the tf-CRE
jj = match(colnames(xx), rownames(mocs))
missed = which(is.na(jj))
cat(length(missed), ' CRE missed \n')

ii = which(!is.na(jj))
jj = jj[ii]

xx = xx[,ii]
mocs = mocs[jj, ]

target.moc = xx %*% mocs
ss = apply(target.moc, 2, sum)
length(which(ss==0))

## motif to TF association matrix
mapping =readRDS(file = '../data/JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata_manual_rmRedundantUNVALIDED.rds')
# mapping = mapping[, c(3, 2)]

mm = match(colnames(target.moc), mapping$name)
target.moc = target.moc[,which(!is.na(mm))] # filter the unvalided motifs if core motifs are present

#atm = table(mapping$name, mapping$gene)
atm = matrix(0, nrow = ncol(target.moc), ncol = nrow(target.moc))
colnames(atm) = rownames(target.moc)
rownames(atm) = colnames(target.moc)
genes = get_geneName(colnames(atm))

for(n in 1:nrow(atm))
{
  # n = 164
  tfs = unlist(strsplit(as.character(mapping$tfs[which(mapping$name == rownames(atm)[n])]), '_'))
  for(tf in tfs)
  {
    
    jj = which(genes == tf)
    if(length(jj)>0){
      cat(n, '--', rownames(atm)[n], ' -- tfs --', tf, '\n')
      atm[n, jj] = 1
    }
    
  }
}

target.tfs = target.moc %*% atm

x = target.tfs[grep('HOXA13', rownames(target.tfs)), ];
x[which(x>0)]

ss1 = apply(target.tfs, 2, sum) # due to missing motif for associated TFs
#target.tfs = target.tfs[, which(ss>0)]
ss2 = apply(target.tfs>0, 1, sum)

#target.tfs = target.tfs[which(ss>0), ]

saveRDS(target.tfs, file = paste0(RdataDir, '/GRN_priorNetwork_target_TFs.rds'))

########################################################
########################################################
# Section : GENIE3 
# 
########################################################
########################################################
source('myGENIE3.R')

E = readRDS(file = paste0(RdataDir, '/GRNinference_scExprMatrix_346geneID_289TFsymbols.rds'))
target.tfs = readRDS(file = paste0(RdataDir, '/GRN_priorNetwork_target_TFs.rds'))

mm = match(rownames(target.tfs), rownames(E))
E = E[mm, ]
E = as.matrix(E)

sds = apply(E, 1, sd)
sels = which(sds>10^-4)
E = E[sels, ]
target.tfs = target.tfs[sels, sels]

#tic()
#weight.matrix1 = GENIE3(expr.matrix = E, priorRegulators =  target.tfs)
#toc()

tic()
wtm = GENIE3(expr.matrix = E, priorRegulators =  target.tfs, ncore = 8)
saveRDS(wtm, file = paste0(RdataDir, '/first_test_Genie3_v2.rds'))

toc()


wtm = readRDS(file =  paste0(RdataDir, '/first_test_Genie3_v2.rds'))
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



