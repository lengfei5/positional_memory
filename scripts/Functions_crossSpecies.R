##########################################################################
##########################################################################
# Project: positional memory 
# Script purpose: set of functions for cross-species analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jul  7 11:00:22 2022
##########################################################################
##########################################################################
##########################################
# GO term enrichment help function
##########################################
library(tidyr)
library(dplyr)
require(ggplot2)
library("gridExtra")
library("cowplot")
require(ggpubr)
library(enrichplot)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(stringr)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(tidyr)
library(dplyr)
require(ggplot2)
library("gridExtra")
library("cowplot")
require(ggpubr)
#require(ChIPseeker)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

run_enrichGo_axolotl = function(pp.annots, distanceToTSS = 2000,  regulation = 'DE_down', 
                                title.plot = 'GoEnrichment')
{
  pp.annots = as.data.frame(pp.annots)
  pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < distanceToTSS), ]
  
  tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM_geneCorrection_v3.rds'))
  bgs = tss$gene
  
  # tss = tss[which(tss$groups == 'DE_up'|tss$groups == 'DE_down'), ]
  
  if(regulation != 'both')  {
    tss = tss[which(tss$groups == regulation), ]
  }else {
    tss = tss[which(tss$groups == 'DE_up'|tss$groups == 'DE_down'), ]
  }
  
  pp = pp.annots[!is.na(match(pp.annots$geneId, tss$geneID)), ]
  ggs = tss$gene[match(pp$geneId, tss$geneID)]
  ggs = unique(as.character(ggs))
  cat(length(ggs), ' targets \n')
  
  gg.expressed = ggs # gene set 
  xx0 = bgs # background
  
  gg.expressed = gg.expressed[which(gg.expressed != '' & gg.expressed != 'N/A' & !is.na(gg.expressed))]
  bgs = unique(xx0[which(xx0 != '' & xx0 != 'N/A' & !is.na(xx0))])
  
  gg.expressed = firstup(tolower(gg.expressed))
  bgs = firstup(tolower(bgs))
  
  gene.df <- bitr(gg.expressed, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  bgs.df <- bitr(bgs, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Mm.eg.db)
  
  ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                   universe     = bgs.df$ENSEMBL,
                   #OrgDb         = org.Hs.eg.db,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.1)
  
  # kk = grep('deoxyribonucleic', ego@result$Description)
  # ego@result$Description[kk] = 'endonuclease activity'
  
  eval(parse(text = paste0('p = barplot(ego, showCategory=20) + ggtitle("', title.plot, '")')))
  plot(p)
  
  return(p)
  
}

########################################################
########################################################
# Section : find conserved targets for conserved regulators
# 
########################################################
########################################################
collect_target_footprint = function(species = 'axolotl', motif = 'RUNX', distanceToTSS = 2000, 
                                    filter.with.DEgene = FALSE) 
{
  if(species == 'zebrafishFin'){
    species = 'zebrafishFin'
    cat(' star the footprinting analysis for zebrafish fin\n')
    
    # gene annotation
    gtf.file = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/zebrafish/GRCz11/Danio_rerio.GRCz11.106.gtf'
    gtf = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
    
    dir.list = list.dirs(path = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/",
                                       "other_species_atac_rnaseq/zebrafish_fin_Lee2020/footprinting"), 
                         recursive = FALSE,  full.names = TRUE)
    
    dir.list = dir.list[grep('sp7ne4dpa_SRR8587716.SRR8587717.merged|sp7po4dpa_SRR8587720.SRR8587721', dir.list)]
    
    ##### RUNX targets
    motif = 'RUNX' 
    for(m in 1:length(dir.list))
    {
      # m = 1
      subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
      subdir = subdir[grep(motif, subdir)]
      subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
      subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
      
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
    saveFootprint = FALSE
    if(saveFootprint){
      
      saveDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/other_species_atac_rnaseq/CNEs_test'
      
      xx = resize(footprint, width = 300, fix = 'center');
      xx = xx[which(start(xx)>0)]
      export(xx,
             con = paste0(saveDir, '/', species, '_', motif, '.bed'), format = 'bed')
      
    }
    
    pp.annots = annotatePeak(footprint, TxDb=gtf, tssRegion = c(-2000, 2000), level = 'transcript')
    
    ## select targets with maximum distance of 5kb
    pp.annots = as.data.frame(pp.annots)
    pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 5000), ]
    cat(length(unique(pp.annots$geneId)), 'targets by promoter footprinting \n')
    
    ### DE genes from RNA-seq data
    res = readRDS(file = paste0('../results/cross_species_20220621/zebrafish_fin/Rdata', 
                                '/RNAseq_fpm_DEgenes_lfcShrink_res.rds'))
    
    fdr.cutoff = 0.1; logfc.cutoff = 0 # select only activated genes in regeneration
    jj = which((res$padj_sp7ne_4dpa.vs.0dpa < fdr.cutoff & res$log2FoldChange_sp7ne_4dpa.vs.0dpa > logfc.cutoff) |
                 (res$padj_sp7po_4dpa.vs.0dpa < fdr.cutoff & (res$log2FoldChange_sp7po_4dpa.vs.0dpa) > logfc.cutoff)
    )
    cat(length(jj), '\n')
    DEgenes = res[jj, ]
    
    DEgenes = data.frame(DEgenes)
    DEgenes$gene = get_geneName(rownames(DEgenes))
    DEgenes$geneID = get_geneID(rownames(DEgenes))
    
    pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
    
    ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
    ggs = unique(as.character(ggs))
    cat(length(ggs), ' targets \n')
    
    saveRDS(ggs, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                               '/targetGenes_footprint_', motif, '_', species,  '.rds'))
    
    ## save both gene ID and gene symbol 
    annot.file = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/zebrafish/GRCz11/annotation_ens_biomart.txt'
    annot = read.delim(annot.file)
    
    pp$gene = pp$geneId
    mm = match(pp$gene, annot$Gene.stable.ID)
    jj = which(annot$Human.gene.name[mm] != '' & !is.na(mm))
    mm = mm[jj]
    pp$gene[jj] = paste0(annot$Human.gene.name[mm], '_',  annot$Gene.stable.ID[mm])
    
    saveRDS(pp, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                              '/targetGenes_footprint_', motif, '_', species,  '_geneID_coordinate_4bed.rds'))
    
  }
  
  if(species == 'zebrafishHeart'){
    
    species = 'zebrafishHeart'
    cat(' star the footprinting analysis \n')
    
    # gene annotation
    gtf.file = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/zebrafish/GRCz11/Danio_rerio.GRCz11.106.gtf'
    gtf = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
    
    
    dir.list = list.dirs(path = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/",
                                       "other_species_atac_rnaseq/zebrafish_heart_Cao2022/footprinting"), 
                         recursive = FALSE,  full.names = TRUE)
    dir.list = dir.list[grep('zebrah_3dpa_|zebrah_7dpa', dir.list)]
    
    ##### RUNX targets
    motif = 'RUNX' 
    for(m in 1:length(dir.list))
    {
      # m = 1
      subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
      subdir = subdir[grep(motif, subdir)]
      subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
      subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
      
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
    
    saveFootprint = FALSE
    if(saveFootprint){
      saveDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/other_species_atac_rnaseq/CNEs_test'
      xx = resize(footprint, width = 300, fix = 'center');
      xx = xx[which(start(xx)>0)]
      export(xx,
             con = paste0(saveDir, '/', species, '_', motif, '.bed'), format = 'bed')
    }
    
    ## save the target genes
    pp.annots = as.data.frame(pp.annots)
    pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 5000), ]
    cat(length(unique(pp.annots$geneId)), 'targets by promoter footprinting \n')
    
    # DE genes from RNA-seq data
    res = readRDS(file = paste0('../results/cross_species_20220621/zebrafish_heart', 
                                '/Rdata/RNAseq_fpm_DEgenes_lfcShrink_res.rds'))
    
    fdr.cutoff = 0.1; logfc.cutoff =  # select only activated genes in regeneration
    jj = which((res$padj_3dpa.vs.0dpa < fdr.cutoff & res$log2FoldChange_3dpa.vs.0dpa > logfc.cutoff) |
                 (res$padj_7dpa.vs.0dpa < fdr.cutoff & (res$log2FoldChange_7dpa.vs.0dpa) > logfc.cutoff)
    )
    cat(length(jj), '\n')
    
    #DEgenes = res[jj, ]
    
    DEgenes = res
        
    DEgenes = data.frame(DEgenes)
    DEgenes$gene = get_geneName(rownames(DEgenes))
    DEgenes$geneID = get_geneID(rownames(DEgenes))
    
    pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
    ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
    ggs = unique(as.character(ggs))
    cat(length(ggs), ' targets \n')
    
    saveRDS(ggs, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                               '/targetGenes_footprint_', motif, '_', species,  '.rds'))
    
    ## save gene ID and gene symbols
    annot.file = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/zebrafish/GRCz11/annotation_ens_biomart.txt'
    annot = read.delim(annot.file)
    
    pp$gene = pp$geneId
    
    mm = match(pp$gene, annot$Gene.stable.ID)
    jj = which(annot$Human.gene.name[mm] != '' & !is.na(mm))
    mm = mm[jj]
    pp$gene[jj] = paste0(annot$Human.gene.name[mm], '_',  annot$Gene.stable.ID[mm])
    
    saveRDS(pp, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                              '/targetGenes_footprint_', motif, '_', species,  '_geneID_coordinate_4bed.rds'))
    
    gg2 = readRDS(paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                         '/targetGenes_footprint_RUNX_zebrafish_fin.rds'))
      
    pp0 = readRDS(paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata/',
                  'targetGenes_footprint_RUNX_zebrafishFin_geneID_coordinate_4bed.rds'))
    
    mm = match(pp0$geneId, pp$geneId)
    pp0 = pp0[which(!is.na(mm)), ]
    
  }
  
  if(species == 'axolotl'){
    
    species = 'axolotl'
    RdataDir = "../results/Rxxxx_R10723_R11637_R12810_atac/Rdata"
    gtf.file =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
    amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
    
    dir.list = list.dirs(path = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/',
                                       'Data/atacseq_using/footprinting'), recursive = FALSE,  full.names = TRUE)
    dir.list = dir.list[grep('ATAC_footprint_BL_UA_5days|ATAC_footprint_BL_UA_9days', dir.list)]
    dir.list = dir.list[grep('_13616', dir.list)]
    
    for(m in 1:length(dir.list))
    {
      # m = 1
      subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
      subdir = subdir[grep(motif, subdir)]
      subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
      subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
      
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
    
    saveFootprint = FALSE
    if(saveFootprint){
      
      saveDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/other_species_atac_rnaseq'
      export(resize(footprint, width = 300, fix = 'center'),
             con = paste0(saveDir, '/', species, '_', motif, '.bed'), format = 'bed')
      
    }
    
    pp.annots = annotatePeak(footprint, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    ## save the target genes
    pp.annots = as.data.frame(pp.annots)
    pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 5000), ]
    cat(length(unique(pp.annots$geneId)), 'targets by promoter footprinting \n')
    
    ## DEgenes by RNA-seq data
    res = readRDS(file = 
                    paste0('../results/RNAseq_data_used/Rdata/', 
            'smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked_RNAseq_data_used_20220408.rds'))
    
    fdr.cutoff = 0.05
    logfc.cutoff = 0
    
    select = which(
      (res$padj_dpa5.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_dpa5.vs.mUA) > logfc.cutoff) | 
        (res$padj_dpa9.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_dpa9.vs.mUA) > logfc.cutoff) |
        (res$padj_dpa13prox.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_dpa13prox.vs.mUA) > logfc.cutoff) |
        (res$padj_dpa13dist.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_dpa13dist.vs.mUA) > logfc.cutoff))
    
    cat(length(select), ' DE genes \n')
    DEgenes = res[select, ]
    DEgenes$gene = get_geneName(rownames(DEgenes))
    DEgenes$geneID = get_geneID(rownames(DEgenes))
    
    pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
    ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
    ggs = unique(as.character(ggs))
    cat(length(ggs), ' targets found for ', motif, ' \n')
    
    saveRDS(ggs, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata', 
                               '/targetGenes_footprint_', motif, '_axolotl.rds'))

    gg1 = readRDS(paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                     '/targetGenes_footprint_RUNX_zebrafish_fin.rds'))
    gg2 = readRDS(paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                         '/targetGenes_footprint_RUNX_zebrafishHeart.rds'))
    
    xx = intersect(gg1, gg2)
    intersect(xx, ggs)
    
    pp$gene = rownames(DEgenes)[match(pp$geneId, DEgenes$geneID)]
    
    saveRDS(pp, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                              '/targetGenes_footprint_', motif, '_', species,  '_geneID_coordinate_4bed.rds'))
    
  }
  
  if(species == 'acoel'){
    
    species = 'acoel'
    cat(' star the analysis of footprinting target\n')
    
    gtf.file = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/acoel/Hmi_1.0/annotation/hmi_annotated_contig.CABFPE.gtf'
    
    # gene annotation
    gtf = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
    
    dir.list = list.dirs(path = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/",
                                       "other_species_atac_rnaseq/acoel_Gehrke2019/footprinting"), 
                         recursive = FALSE,  full.names = TRUE)
    
    dir.list = dir.list[grep('hmia_12h_tail|hmia_24h_tail|_48h_', dir.list)]
    
    ##### RUNX targets
    motif = 'RUNX' 
    for(m in 1:length(dir.list))
    {
      # m = 1
      subdir = list.dirs(dir.list[m], recursive = FALSE, full.names = TRUE)
      subdir = subdir[grep(motif, subdir)]
      #subdir = subdir[grep('RUNX3', subdir, invert = TRUE)]
      subdir = subdir[grep('MA0002.1', subdir, invert = TRUE)]
      
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
    
    saveFootprint = FALSE
    if(saveFootprint){
      saveFootprint = FALSE
      if(saveFootprint){
        saveDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/other_species_atac_rnaseq/CNEs_test'
        xx = resize(footprint, width = 300, fix = 'center');
        xx = xx[which(start(xx)>0)]
        export(xx,
               con = paste0(saveDir, '/', species, '_', motif, '.bed'), format = 'bed')
      }
      
    }
    
    
    pp.annots = annotatePeak(footprint, TxDb=gtf, tssRegion = c(-2000, 2000), level = 'transcript')
    pp.annots = as.data.frame(pp.annots)
    pp.annots = pp.annots[which(abs(pp.annots$distanceToTSS) < 5000), ]
    cat(length(unique(pp.annots$geneId)), 'targets by promoter footprinting \n')
    
    # DE genes from RNA-seq data
    #res = readRDS(file = paste0(RdataDir, '/RNAseq_fpm_DEgenes_lfcShrink_res.rds'))
    res = readRDS(file = paste0('../results/cross_species_20220621/acoel', 
                                '/Rdata/RNAseq_fpm_DEgenes_lfcShrink_res.rds'))
    fdr.cutoff = 0.1; logfc.cutoff = 0 # select only activated genes in regeneration
    
    # fdr.cutoff = 0.1; logfc.cutoff = 0 # select only activated genes in regeneration
    jj = which((res$padj_tail.1h.vs.0h < fdr.cutoff & res$log2FoldChange_tail.1h.vs.0h > logfc.cutoff) |
                 (res$padj_tail.1h.vs.0h < fdr.cutoff & (res$log2FoldChange_tail.3h.vs.0h) > logfc.cutoff)|
                 (res$padj_tail.6h.vs.0h < fdr.cutoff & (res$log2FoldChange_tail.6h.vs.0h) > logfc.cutoff) |
                 (res$padj_tail.12h.vs.0h < fdr.cutoff & (res$log2FoldChange_tail.12h.vs.0h) > logfc.cutoff)|
                 (res$padj_head.1h.vs.0h < fdr.cutoff & res$log2FoldChange_head.1h.vs.0h > logfc.cutoff) |
                 (res$padj_head.3h.vs.0h < fdr.cutoff & (res$log2FoldChange_head.3h.vs.0h) > logfc.cutoff)|
                 (res$padj_head.6h.vs.0h < fdr.cutoff & (res$log2FoldChange_head.6h.vs.0h) > logfc.cutoff) |
                 (res$padj_head.12h.vs.0h < fdr.cutoff & (res$log2FoldChange_head.12h.vs.0h) > logfc.cutoff)
               
    )
    cat(length(jj), '\n')
    
    DEgenes = res[jj, ]
    DEgenes = data.frame(DEgenes)
    DEgenes$gene = as.character(sapply(rownames(DEgenes), function(x) unlist(strsplit(x, '[|]'))[1]))
    DEgenes$geneID = DEgenes$gene  
    
    ggs = pp.annots$geneId
    ggs = as.character(sapply(ggs, function(x) unlist(strsplit(x, '[|]'))[1]))
    #ggs = gsub('_HUMAN', '', ggs)
    #ggs = gsub('[|]', '_', ggs)
    #ggs = get_geneName(ggs)
    ggs = unique(as.character(ggs))
    
    ggs = ggs[!is.na(match(ggs, DEgenes$gene))]
    cat(length(ggs), ' targets found for ', motif, ' \n')
    
    ggs = ggs[grep('_HUMAN', ggs)]
    
    ggs = as.character(sapply(ggs, function(x) unlist(strsplit(x, '[|]'))[1]))
    
    write.table(ggs, file = paste0("../results/Uniprot_EntryName.txt"), sep = '\t',
                col.names = FALSE, row.names = FALSE, quote =FALSE)
    
    ggs =read.delim('../data/uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.07.07-16.03.50.51.tsv', 
                    sep = '\t')
    ggs$geneSymbols = sapply(ggs$Gene.Names, function(x) unlist(strsplit(as.character(x), ' '))[1])
    #pp = pp.annots[!is.na(match(pp.annots$geneId, DEgenes$geneID)), ]
    #ggs = DEgenes$gene[match(pp$geneId, DEgenes$geneID)]
    #ggs = unique(as.character(ggs))
    #cat(length(ggs), ' targets \n')
    saveRDS(ggs$geneSymbols, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                               '/targetGenes_footprint_', motif, '_', species,  '.rds'))
    
    pp = pp.annots[!is.na(match(pp.annots$geneId, rownames(DEgenes))), ]
    
    pp$gene = DEgenes$geneID[match(pp$geneId, rownames(DEgenes))]
    pp$geneSymbols = ggs$geneSymbols[match(pp$gene, ggs$From)]
    
    
    saveRDS(pp, file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata',
                              '/targetGenes_footprint_', motif, '_', species,  '_geneID_coordinate_4bed.rds'))
    
  }
  
  return(ggs)
  
}  
  

extract_proteinSeq_acoel = function()
{
  library('rtracklayer')
  library(GenomicRanges)
  library('GenomicFeatures')
  gtf.file = paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/', 
                    'acoel/Hmi_1.0/annotation/hmi_annotated_contig.CABFPE_CDS.gtf')
  annot = import(gtf.file)
  annot = data.frame(annot)
  annot$name = paste0(annot$transcript_id, ';', annot$gene_id)
  annot = annot[, c(1, 2, 3, 12, 8, 5)]
  annot$score = 0
  
  export(annot, paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/', 
                       'acoel/Hmi_1.0/annotation/hmi_annotated_contig.CABFPE_CDS.bed'), format = 'bed')
}

find_RUNX_orthologues_orthofinder = function()
{
  outDir_orthofinder = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/', 
                              'Orthofinder/Results_Jul18/Orthologues/Orthologues_8296_Ambystoma_mexicanum')
  
  ax = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_axolotl_geneID_coordinate_4bed.rds'))
  
  targets = data.frame(ax$gene, ax$geneId, ax$transcriptId, stringsAsFactors = FALSE)
  targets$ortho_zebrafishFin = NA
  targets$ortho_zebrafishHeart = NA
  targets$ortho_acoel = NA
  
  ## start the orthologues in fish fin
  annotFish = 
  read.delim('/Volumes/groups/tanaka/People/current/jiwang/Genomes/zebrafish/GRCz11/annotation_ens_biomart_proteinID.txt')
  
  pp1 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_zebrafishFin_geneID_coordinate_4bed.rds'))
  pp2 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_zebrafishHeart_geneID_coordinate_4bed.rds'))
  ort1 = read.delim(paste0(outDir_orthofinder, '/8296_Ambystoma_mexicanum__v__7955_Danio_rerio.tsv'), 
                      header = TRUE)
  pp3 = readRDS(paste0(RdataDir, '/targetGenes_footprint_RUNX_acoel_geneID_coordinate_4bed.rds'))
  pp3$gene = sapply(pp3$geneId, function(x){
    x = unlist(strsplit(as.character(x), '[|]'))
    x[length(x)]
  })
  ort3 = read.delim(paste0(outDir_orthofinder, '/8296_Ambystoma_mexicanum__v__442651_Hofstenia_miamia.tsv'), 
                    header = TRUE)
  
  for(n in 1:nrow(targets))
  {
    # n = 647
    transcript = targets$ax.transcriptId[n]
    transcript = unlist(strsplit(as.character(transcript), '[|]'))
    transcript = transcript[length(transcript)]
    transcript = unlist(strsplit(as.character(transcript), '[.]'))[1]
    kk1 = grep(transcript, ort1$X8296_Ambystoma_mexicanum)
    kk3 = grep(transcript, ort3$X8296_Ambystoma_mexicanum)
    
    cat(transcript, " -- fish orth: ", kk1, " ; acoel orth:", kk3, '\n')
    if(length(kk1)>0){
      test1 = ort1$X7955_Danio_rerio[kk1]
      test1 = unlist(strsplit(as.character(test1), ', '))
      test1 = sapply(test1, function(x){
        x = unlist(strsplit(as.character(x), '[|]'))
        x = x[grep('tr', x, invert = TRUE)]
        x = unique(gsub('_DANRE', '', x))
        as.character(x)
      })
      
      test1 = annotFish$Gene.stable.ID[match(test1, annotFish$UniProtKB.Gene.Name.ID)]
      mm1 = match(test1, pp1$geneId)
      if(!all(is.na(mm1))) targets$ortho_zebrafishFin[n] = 1
      
      mm2 = match(test1, pp2$geneId)
      if(!all(is.na(mm2))) targets$ortho_zebrafishHeart[n] = 1
    }
    
    if(length(kk3)>0){
      test3 = ort3$X442651_Hofstenia_miamia[kk3]
      test3 = unlist(strsplit(as.character(test3), ', '))
      test3 = sapply(test3, function(x) {unlist(strsplit(as.character(x), '[.]'))[1]})
      #test3 = test3[grep('tr', test3, invert = TRUE)]
      #test3 = unique(gsub('_DANRE', '', test3))
      #test3 = annotFish$Gene.stable.ID[tch(test3, annotFish$UniProtKB.Gene.Name.ID)]
      mm3 = match(test3, pp3$gene)
      if(!all(is.na(mm3))) targets$ortho_acoel[n] = 1
      
    }
  }
  
  ss = apply(as.matrix(targets[, c(4:5)]), 1, function(x) sum(x, na.rm = TRUE))
  
  targets[ss>=2, c(1, 4:6)]
  
  
}

