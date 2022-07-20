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
source('Functions_rnaseq.R')
require(openxlsx)
require(ggplot2)
require(DESeq2)
require(dplyr)
require(gridExtra)

version.Data = 'RNAseq_data_used';
version.analysis = paste0("_", version.Data, "_20220408") # in the version 20220408, counts were done with update annot with Hox patch

## Directories to save results 
dataDir = "/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/rnaseq_using/"
design.file = paste0(dataDir, 'sampleInfos_parsed.txt')

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
# Section III: positonal genes from both microarray and RNA-seq
# mainly TFs that are changed in mature UA, LA and Hand (head as a control)
# the main analysis of segment-specific genes were using microarray data 
#
########################################################
########################################################

##########################################
# rerun the limma test for microarray
##########################################
Run.limma.test = FALSE
if(Run.limma.test){
  require(limma)
  
  Rdata.microarray = "../results/microarray/Rdata/"
  load(file = paste0(Rdata.microarray, 'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary.Rdata'))
  
  mm = match(rownames(res), annot$geneID)
  ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
  rownames(res)[!is.na(mm)] = ggs[!is.na(mm)]
  
  f <- factor(rep(c('mUA', 'mLA', 'mHand'), each = 3), levels=c("mUA","mLA","mHand")) 
  design <- model.matrix(~0+f)
  colnames(design) <- c("mUA","mLA","mHand")
  
  #To make all pair-wise comparisons between the three groups one could proceed
  fit <- lmFit(res, design)
  contrast.matrix <- makeContrasts(mLA-mUA, mHand-mUA, mHand-mLA, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  #A list of top genes for RNA2 versus RNA1 can be obtained from
  topTable(fit2, coef=1, adjust="BH")
  
  #The outcome of each hypothesis test can be assigned using
  results <- decideTests(fit2)
  
  xx = data.frame(fit2$p.value)
  colnames(xx) = c('mLA.vs.mUA', 'mHand.vs.mUA', 'mHand.vs.mLA')
  xx$pval.max = apply(as.matrix(-log10(xx)), 1, max)
  
  res = data.frame(res, xx)
  res = res[order(-res$pval.max), ]
  
  save(res, raw, fit2, file = paste0(RdataDir, 
                                     'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_DEpval.Rdata'))
  
}

require(limma)
# load microarray analysis results: log2 signal (res), comparison (fit2), and raw data (raw)
load(file = paste0("../results/microarray/Rdata/", 
                   'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_DEpval.Rdata'))

res = res[, c(1:9)] # drop the pval kept previously

tops = topTable(fit2, coef=1, adjust="BH", number = nrow(res), genelist = rownames(res))[, c(2, 5,6)]
colnames(tops) = paste0(colnames(tops), '_mLA.vs.mUA')
res = data.frame(res, tops[match(rownames(res), rownames(tops)), ])

tops = topTable(fit2, coef=2, adjust="BH", number = nrow(res), genelist = rownames(res))[, c(2, 5,6)]
colnames(tops) = paste0(colnames(tops), '_mHand.vs.mUA')
res = data.frame(res, tops[match(rownames(res), rownames(tops)), ])

tops = topTable(fit2, coef=3, adjust="BH", number = nrow(res), genelist = rownames(res))[, c(2, 5,6)]
colnames(tops) = paste0(colnames(tops), '_mHand.vs.mLA')
res = data.frame(res, tops[match(rownames(res), rownames(tops)), ])

# plot all position-dependent genes
library("pheatmap")

res$fdr.max = apply(-log10(res[, grep('adj.P.Val_', colnames(res))]), 1, max)
res$logFC.max = apply((res[, grep('logFC_', colnames(res))]), 1, function(x) return(x[which(abs(x)==max(abs(x)))][1]))

o1 = order(-res$fdr.max)
res = res[o1, ]

saveRDS(res, file = paste0("../results/microarray/Rdata/", 
                     'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_limma.DE.stats.rds'))


##########################################
# select DE genes and make plots 
##########################################
res = readRDS(file = paste0("../results/microarray/Rdata/", 
'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_limma.DE.stats.rds'))

ggs = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})

qv.cutoff = 0.05
logfc.cutoff = 1
select = which(res$fdr.max> -log10(qv.cutoff) & abs(res$logFC.max)> 0)

select = which(res$adj.P.Val_mHand.vs.mLA < qv.cutoff & abs(res$logFC_mHand.vs.mLA) > logfc.cutoff|
                 res$adj.P.Val_mHand.vs.mUA < qv.cutoff & abs(res$logFC_mHand.vs.mUA) > logfc.cutoff |
                 res$adj.P.Val_mHand.vs.mLA < qv.cutoff & abs(res$logFC_mHand.vs.mLA) > logfc.cutoff )
# cat(length(select), ' positional genes found \n')
cat(length(select), ' DE genes selected \n')

ggs = ggs[select]

cat(ggs[grep('MEIS', ggs)], '\n')
print(intersect(ggs, tfs))

saveRDS(intersect(ggs, tfs), file = paste0(RdataDir, '/postional_genes_tfs.rds'))

print(intersect(ggs, sps))
print(intersect(ggs, eps))
print(intersect(ggs, rbp))

yy = res[select, ]
saveRDS(yy, file = paste0(RdataDir, 'microarray_positionalGenes_data.rds'))

df <- data.frame(condition = rep(c('mUA', 'mLA', 'mHand'), each = 3))
rownames(df) = colnames(yy)
colnames(df) = 'segments'

annot_colors = c('springgreen4', 'steelblue2', 'gold2')
names(annot_colors) = c('mUA', 'mLA', 'mHand')
annot_colors = list(segments = annot_colors)

sample_colors = c('springgreen4', 'steelblue2', 'gold2')
names(sample_colors) = c('mUA', 'mLA', 'mHand')
annot_colors = list(segments = sample_colors)

pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(16), 
         show_colnames = FALSE,
         scale = 'row',
         cluster_cols=FALSE, annotation_col=df,
         annotation_colors = annot_colors,
         width = 6, height = 12, 
         filename = paste0(figureDir, '/Fig2A_heatmap_DEgenes_matureSample_fdr.0.05_log2fc.1_microarray.pdf')) 

pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(32), 
         show_colnames = FALSE,
         scale = 'none',
         cluster_cols=FALSE, annotation_col=df,
         annotation_colors = annot_colors,
         width = 6, height = 12, 
         filename = paste0(figureDir, '/Fig2A_heatmap_DEgenes_matureSample_fdr.0.05_log2fc.1_microarray_nonScaled.pdf')) 


if(saveTables){
  source('Functions_histM.R')
  
  yy = readRDS(file = paste0(RdataDir, 'microarray_positionalGenes_data.rds'))
  yy = data.frame(gene = get_geneName(rownames(yy)), geneID = rownames(yy), 
                  specificSegment = 'mHand', yy, stringsAsFactors = FALSE)
  yy$specificSegment[which(yy$logFC_mHand.vs.mUA<0)] = 'mUA'
  
  write.csv(yy[, c(1:3)],
            file = paste0(tableDir, '/position_geneList_microarray.csv'), 
            quote = FALSE, row.names = FALSE)
  
  #jj = grep('^CYP26|^RARA|^RARB|^RARG|SCUBE2|FGF9|FGF2_|FGFR2|FGFR4', rownames(res))
  jj = unique(grep('MEIS1_AMEX60DD024424|MEIS|LRAT|FGF9|FGF2_|FGFR2|FGFR4', rownames(res)))
  yy = res[jj, c(1:18)]
  
  yy = yy[, c(1:9, grep('_mHand.vs.mUA', colnames(yy)))]
  
  write.csv(yy,
           file = paste0(tableDir, '/MEIS_LRAT_FGF_geneList_microarray.csv'), 
           quote = FALSE, row.names = TRUE)
  
  write.csv(res[select, ],
            file = paste0(tableDir, '/position_dependent_genes_from_matureSamples_microarray_qv.0.05_log2FC.1.csv'), 
            quote = FALSE, row.names = TRUE)
  
  yy = cal_sample_means(res[, c(1:9)], conds = c('mUA', 'mLA', 'mHand'))
  
  yy = cbind(yy, res[,  grep('_mHand.vs.mUA', colnames(res))])
  
  write.csv(yy,
            file = paste0(tableDir, '/allGenes_microarray.csv'), 
            quote = FALSE, row.names = TRUE)
  
  
  
  
}

##########################################
# GO term analysis of DE genes 
##########################################
library(enrichplot)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(stringr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

clean_geneNames = function(gg.expressed)
{
  gg.expressed = unique(unlist(lapply(gg.expressed, 
                                      function(x) { x = unlist(strsplit(as.character(x), '_'));  return(x[length(x)])})))
  gg.expressed = unique(annot$gene.symbol.toUse[match(gg.expressed, annot$geneID)])
  gg.expressed = gg.expressed[which(gg.expressed != '' & gg.expressed != 'N/A' & !is.na(gg.expressed))]
  
  return(gg.expressed)
}

res = readRDS(file = paste0("../results/microarray/Rdata/", 
            'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_limma.DE.stats.rds'))
yy = readRDS(file = paste0(RdataDir, 'microarray_positionalGenes_data.rds'))

# background
bgs0 = unique(clean_geneNames(rownames(res)))

xx0 = unique(annot$gene.symbol.toUse)
xx0 = xx0[which(xx0 != '' & xx0 != 'N/A' & !is.na(xx0))]
bgs = unique(xx0)

gg.expressed = unique(clean_geneNames(rownames(yy)[which(yy$logFC_mHand.vs.mUA < 0)]))
cat('# of genes --', length(gg.expressed), '\n')

#gg.expressed = firstup(tolower(gg.expressed))
#bgs = firstup(tolower(bgs))
#bgs0 = firstup(tolower(bgs0))

gene.df <- bitr(gg.expressed, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

bgs.df <- bitr(bgs, fromType = "SYMBOL",
               toType = c("ENSEMBL", "ENTREZID"),
               OrgDb = org.Hs.eg.db)
head(bgs.df)

bgs0.df <- bitr(bgs0, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

#pval.cutoff = 0.05
ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                 universe     = bgs0.df$ENSEMBL,
                 #universe     = bgs.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 #OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.1,
                 qvalueCutoff  = 0.3, 
                 minGSSize = 3)

#head(ego)

barplot(ego) + ggtitle("Go term enrichment for positional genes")

kegg = enrichKEGG(gene = gene.df$ENTREZID,
                  organism = 'hsa', 
                  keyType = 'kegg', 
                  universe = bgs0.df$ENTREZID, 
                  minGSSize = 5)
barplot(kegg)

#edox <- setReadable(ego, 'org.Mm.eg.db', 'ENSEMBL')
pdfname = paste0(figureDir, 'GOterm_postionalGenes_mUASpecific.pdf')
pdf(pdfname, width = 12, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

dotplot(ego, showCategory=30) + ggtitle("positional genes")

dev.off()

write.csv(ego, file = paste0(tableDir, "GO_term_enrichmenet_for_positional_genes_Microarray.csv"), 
          row.names = TRUE)

SaveGoTermName = FALSE
if(SaveGoTermName){
  
  ego = read.csv(file = paste0(tableDir, "FigS1EGO_term_enrichmenet_for_positional_genes_Microarray.csv"))
  
  
}

##########################################
# volcano plot with highlighting genes
##########################################
library(ggrepel)
library(dplyr)
library(tibble)
res = readRDS(file = paste0("../results/microarray/Rdata/", 
      'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_limma.DE.stats.rds'))

res$gene = sapply(rownames(res), function(x) unlist(strsplit(as.character(x), '_'))[1])

for(comp in c('mHand.vs.mUA', 'mHand.vs.mLA', 'mLA.vs.mUA'))
{
  comp =  'mHand.vs.mUA'
  
  res$fdr = eval(parse(text = paste0('-log10(res$adj.P.Val_', comp, ')')))
  res$pval = eval(parse(text = paste0('-log10(res$P.Value_', comp, ')')))
  res$logfc = eval(parse(text = paste0('res$logFC_', comp)))
  
  #examples.sel = which(res$pval > 4 & abs(res$logfc) > 2)
  examples.sel = c()
  examples.sel = unique(c(examples.sel, grep('HOXA13|HOXD13', res$gene)))
  
  ggplot(data=res, aes(x=logfc, y=pval, label = gene)) +
    geom_point(size = 0.3) + 
    geom_point(data=res[which(res$logfc > 1 & res$pval > -log10(0.001)), ], aes(x=logfc, y=pval), colour="red", size=0.3) +
    geom_point(data=res[which(res$logfc < -1 & res$pval > -log10(0.001)), ], aes(x=logfc, y=pval), colour="blue", size=0.3) + 
    geom_point(data=res[c(examples.sel), ], aes(x=logfc, y=pval), colour="darkgreen", size=1.5) +
    theme_classic() + 
    geom_text_repel(data= res[c(examples.sel), ], 
                    aes(x=logfc, y=pval),
                    size = 5,
                    color = "darkgreen",
                    #family = 'Times',
                    fontface = 'bold',
                    # Add extra padding around each text label.
                    box.padding = unit(0.3, 'lines'),
                    # Add extra padding around each data point.
                    point.padding = unit(1.6, 'lines')) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14), 
          axis.text.y = element_text(angle = 0, size = 14), 
          axis.title =  element_text(size = 14),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14)
          #legend.position=c(0.2, 0.8),
          #plot.margin = margin()
          #legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    ) +
    geom_vline(xintercept=c(-1, 1), col='gray') +
    geom_hline(yintercept=-log10(0.001), col="gray") +
    labs(x = "log2FC (mHand/mUA)", y = "-log10(p-val)")
  
  ggsave(paste0(figureDir, "Fig2C_VolcanoPlot_log2FC_pval_microarray_noLabels_", comp, ".pdf"), width=6, height = 4)
  
  
}


##########################################
# highlight TFs and EPs in positional genes 
##########################################
# plot only TFs and SPs
ggs = rownames(yy)
ggs = sapply(ggs, function(x) unlist(strsplit(as.character(x), '_'))[1])

mm = match(ggs, unique(c(tfs, eps)))
yy1 = yy[unique(c(which(!is.na(mm)))), ]

pheatmap(yy1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(16), 
         scale = 'none',
         cluster_cols=FALSE, annotation_col=df, fontsize_row = 8, 
         width = 8, height = 10,
         annotation_colors = annot_colors,
         filename = paste0(figureDir, 'FigS2A_heatmap_DE.tfs_eps_mature_qv.0.05_log2fc.1_microarray_lgo2.geneExp_nonScaled.pdf')) 

pheatmap(yy1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(16), 
         scale = 'row',
         cluster_cols=FALSE, annotation_col=df, fontsize_row = 8, 
         width = 8, height = 10,
         annotation_colors = annot_colors,
         filename = paste0(figureDir, 'FigS2A_heatmap_DE.tfs_eps_mature_qv.0.05_log2fc.1_microarray_scaled.log2.geneExp_zscore.pdf')) 

mm = match(ggs, unique(c(toupper(sps))))
yy1 = yy[unique(c(which(!is.na(mm)), grep('CYP2', rownames(yy)))), ]

pheatmap(yy1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE,
         scale = 'none',
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(16), 
         cluster_cols=FALSE, annotation_col=df, fontsize_row = 8, 
         width = 8, height = 10,
         annotation_colors = annot_colors,
         filename = paste0(figureDir, 'FigS2A_heatmap_DE_sps_mature_qv.0.05_log2FC.1_microarray_log2.geneExp_nonScaled.pdf')) 

pheatmap(yy1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE,
         scale = 'row',
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(16), 
         cluster_cols=FALSE, annotation_col=df, fontsize_row = 8, 
         width = 8, height = 10,
         annotation_colors = annot_colors,
         filename = paste0(figureDir, 'FigS2A_heatmap_DE_sps_mature_qv.0.05_log2FC.1_microarray_log2.geneExp_zscore.pdf')) 

##########################################
# RNA-seq data of mature sample to confirm the postional gene candidates
# 
##########################################
Test.mature.smartseq2.using.HOX.patch.annot = FALSE
if(Test.mature.smartseq2){
  # # load dds normalized object and annotations
  # load(file = paste0(RdataDir, 'RNAseq_design_dds.object.Rdata'))
  # # select mature samples
  # sels = intersect(which(design.matrix$batch == 4), grep('Mature', design.matrix$condition))
  # dds = dds[, sels]
  # design.matrix = design.matrix[sels, ]
  # save(dds, design.matrix, file = paste0(RdataDir, '/Smartseq2_matureSamples_selected_forAnalysis.Rdata'))
  # saveRDS(design.matrix, file = paste0(RdataDir, '/design_smartseq2_matureSamples_selected.rds'))
  design = readRDS(file = paste0(RdataDir, '/design_smartseq2_matureSamples_selected.rds'))
  
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/rnaseq_using/'
  
  xlist = list.files(path=paste0(dataDir, 'featurecounts_Q10'),
                     pattern = "*featureCounts.txt$", full.names = TRUE) ## list of data set to merge
  all = cat.countTable(xlist, countsfrom = 'featureCounts')
  
  colnames(design)[1] = 'SampleID'
  counts = process.countTable(all=all, design = design[, c(1,3)], merge.technicalRep.sameID = FALSE)
  
  
  all = counts
  
  mm = match(all$gene, annot$geneID)
  ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
  all$gene[!is.na(mm)] = ggs[!is.na(mm)]
  
  raw = as.matrix(all[, -1])
  rownames(raw) = all$gene
  #design$batch = paste0(design$request, '_', design$protocol)
  
  dds <- DESeqDataSetFromMatrix(raw, DataFrame(design), design = ~ condition)
  
  save(design, dds, file = paste0(RdataDir, 'design_dds_all_matureSamples_8selectedSamples.batch4_v47.hox.patch.Rdata'))
  
  ##########################################
  # filering and normalization 
  ##########################################
  load(file = paste0(RdataDir, 'design_dds_all_matureSamples_8selectedSamples.batch4_v47.hox.patch.Rdata'))
  
  design$condition = droplevels(design$condition)
  design$batch = droplevels(design$batch)
  table(design$condition, design$batch)
  
  kk = which(design$condition != 'Mature_HEAD')
  design = design[kk, ]
  dds = dds[, kk]
  dds$condition = droplevels(dds$condition)
   
  ss = rowSums(counts(dds))
  
  hist(log10(ss), breaks = 100);abline(v = log10(20), lwd = 2.0, col = 'red')
  cat(length(which(ss>50)), ' gene selected \n')
  
  ddx = dds[which(ss>50), ]
  ddx = estimateSizeFactors(ddx)
  
  sizeFactors(dds) = sizeFactors(ddx)
  
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  pca=plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = FALSE)
  print(pca)
  
  pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE, ntop = 3000))
  ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 1, nudge_y = 1, size=2.5)
  
  ggsave(paste0(resDir, '/PCA_smartseq2_matureSamples.pdf'),  width=12, height = 8)
  
  #ggsave(paste0(resDir, '/PCA_smartseq2_regeneration.timepoints_filteredSamples_onlyR10724.pdf'),  width=12, height = 8)
  
  # cpm = fpm(dds)
  # cpm = log2(cpm + 2^-6)
  # 
  # plot.pair.comparison.plot(cpm[, c(1:6)], linear.scale = FALSE)
  # 
  cpm = fpm(dds)
  cpm = log2(cpm + 2^-4)
  
  colnames(cpm) = gsub('Mature_Hand', 'mHand', colnames(cpm))
  colnames(cpm) = gsub('Mature_LA', 'mLA', colnames(cpm))
  colnames(cpm) = gsub('Mature_UA', 'mUA', colnames(cpm))
  
  dds$condition = droplevels(dds$condition)
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
  
  saveRDS(data.frame(cpm, res, stringsAsFactors = FALSE), 
       file = paste0("../results/RNAseq_data_used/Rdata/matureSamples_cpm_DEgenes_8selectedSamples.batch4_v47.hox.patch.rds"))
  
  # load results from microarray with limma to compare with
  res = readRDS(file = paste0("../results/RNAseq_data_used/Rdata", 
                              "/matureSamples_cpm_DEgenes_8selectedSamples.batch4_v47.hox.patch.rds"))
  res0 = readRDS(file = paste0("../results/microarray/Rdata/", 
       'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_limma.DE.stats.rds'))
  
  res0$gene = sapply(rownames(res0), function(x) unlist(strsplit(as.character(x), '_'))[1])
  res0$geneID = sapply(rownames(res0), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
  
  mm = match(res0$geneID, res$geneID)
  missed = which(is.na(mm))
  mm_missed = match(res0$gene[missed], res$gene)
  mm[missed] = mm_missed 
  
  res = res[mm, ]
  
  select = which(res0$P.Value_mHand.vs.mUA<0.05)
  plot( res0$logFC_mHand.vs.mUA[select], res$log2FoldChange_mHand.vs.mUA[select], cex =0.4);
  abline(0, 1, lwd = 2.0, col = 'red')
  
  ## volcanoplot for smartseq2 data 
  for(comp in c('mHand.vs.mUA', 'mHand.vs.mLA', 'mLA.vs.mUA'))
  {
    # comp = 'mHand.vs.mUA'
    # comp = 'mHand.vs.mLA'
    # comp = 'mLA.vs.mUA'
    res$fdr = eval(parse(text = paste0('-log10(res$padj_', comp, ')')))
    res$pval = eval(parse(text = paste0('-log10(res$pvalue_', comp, ')')))
    res$logfc = eval(parse(text = paste0('res$log2FoldChange_', comp)))
    
    #examples.sel = which(res$pval > 4 & abs(res$logfc) > 2)
    examples.sel = c()
    examples.sel = unique(c(examples.sel, grep('HOXA13|HOXD13', res$gene)))
    
    ggplot(data=res, aes(x=logfc, y=pval, label = gene)) +
      geom_point(size = 0.5) + 
      geom_point(data=res[which(res$logfc > 1 & res$pval > -log10(0.001)), ], aes(x=logfc, y=pval), colour="red", size=1) +
      geom_point(data=res[which(res$logfc < -1 & res$pval > -log10(0.001)), ], aes(x=logfc, y=pval), colour="blue", size=1) +
      theme_classic() + 
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12)) + 
      #geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
      #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
      #scale_color_manual(values=c("blue", "black", "red")) +
      geom_vline(xintercept=c(-1, 1), col='gray') +
      geom_hline(yintercept=-log10(0.001), col="gray") +
      labs(x = "log2FC")
    
    ggsave(paste0(figureDir, "Fig2C_VolcanoPlot_log2FC_pval_smartseq2_noLabels_", comp, ".pdf"), width=12, height = 8)
    
    
    res0$fdr = eval(parse(text = paste0('-log10(res0$adj.P.Val_', comp, ')')))
    res0$pval = eval(parse(text = paste0('-log10(res0$P.Value_', comp, ')')))
    res0$logfc = eval(parse(text = paste0('res0$logFC_', comp)))
    
    qv.cutoff = 0.05; logfc.cutoff = 1
    select = which(res0$fdr > -log10(qv.cutoff) & abs(res0$logfc) > logfc.cutoff)
    
    # mm = match(res$geneID, res0$geneID[select])
    # yy = data.frame(log2fc.smartseq2 = res$logfc[!is.na(mm)],  log2fc.microarray = res0$logfc[mm[!is.na(mm)]])
    
    #plot(yy, xlim = c(-4, 4), ylim = c(-4, 4), cex = 0.3);
    #abline(0, 1, lwd = 2.0)
    cat(cor(yy[, 1], yy[, 2], use = "na.or.complete"), '\n')
    
  }
}

########################################################
########################################################
# Section V : dynamic gene and TFs, SPs in regeneration using smart-seq2 data
# first double check the regeneration sample quality
# because the quick technical replicate merging from different batches did not work well, e.g. dpa5_136150 from two R10724 and R11635  
# it turned out Akane's early samples were polyA and later samples were smartseq2
# the samples here were reprocessed because samples were not merged replicates.
########################################################
########################################################

##########################################
# Be carefull here !!!! 
# select the R10724 samples without sample 136150 with is names 13615x
# only keep the sample 136150 from the request R11635
##########################################
#load(file = paste0('/Users/jiwang/workspace/imp/positional_memory/results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/',
#                   'design_dds_all_regeneration.samples_allBatches.Rdata')) # import the old sample infos
#design$condition = sapply(design$condition, function(x) )
#load(file = paste0(RdataDir, 'RNAseq_design_dds.object.Rdata'))
#design = design.matrix

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/rnaseq_using/regeneration_RNAseq/'
design = read.csv(paste0(dataDir, 'regeneration_samples_toUse.csv'))

xlist = list.files(path=paste0(dataDir, 'featurecounts_Q10'),
                   pattern = "*featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

#colnames(all)[grep('136150s', colnames(all))] = '13615x.txt'
#colnames(all) =  gsub('[#]', '_', colnames(all))
colnames(design)[1] = 'SampleID'
counts = process.countTable(all=all, design = design[, c(1,2)], merge.technicalRep.sameID = FALSE)


all = counts

mm = match(all$gene, annot$geneID)
ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
all$gene[!is.na(mm)] = ggs[!is.na(mm)]

raw = as.matrix(all[, -1])
rownames(raw) = all$gene
#design$batch = paste0(design$request, '_', design$protocol)

dds <- DESeqDataSetFromMatrix(raw, DataFrame(design), design = ~ condition)

save(design, dds, file = paste0(RdataDir, 'design_dds_all_regeneration_12selectedSamples_v47.hox.patch.Rdata'))


##########################################
# filering and normalization 
##########################################
load(file = paste0(RdataDir, 'design_dds_all_regeneration_12selectedSamples_v47.hox.patch.Rdata'))

#design$protocol = gsub(' ', '', design$protocol)
#design$batch = paste0(design$request, '_', design$protocol)
table(design$condition, design$batch)

#sels = which(design$batch == 'R10724_smartseq2' & design$SampleID != '13615x')
#sels = which(design$batch == 'R11635_smartseq2' | (design$batch == 'R10724_smartseq2' & design$SampleID != '13615x'))
#design = design[sels, ]
#write.csv(design, 
#  file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/rnaseq_using/regeneration_RNAseq/', 
#                'regeneration_samples_toUse.csv'), row.names = FALSE, quote = FALSE)
# dds = dds[ ,sels]
#table(design$condition, design$batch)

#dds = DESeqDataSetFromMatrix(raw, DataFrame(design), design = ~ condition)
dds$condition = droplevels(dds$condition)

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100);abline(v = log10(20), lwd = 2.0, col = 'red')
cat(length(which(ss>50)), ' gene selected \n')

dds = dds[which(ss>50), ]

# annot = rbind(annot, c('AMEX60DD029624', rep(NA, 15), "POU6F1"))
# tfs.examples = c("POU6F1", # not expressed 
#                  "ZNF680", # not annotated
#                  "ISL1",  # not expressed
#                  "HNF1A", # not expressed
#                  "POU6F2", # not expressed
#                  "CTCFL", # not annotated
#                  "ZNF701", # found expressed
#                  "YY2", # no YY1, but found YY1
#                  "ZFP335", # not annotated
#                  "LHX1", # not expressed
#                  "NKX6-3", # not expressed not annotated 
#                  "ZNF263", # not annotated
#                  "SNAI3" ) # not annotated but two SNAIL2 found 
# gene2find = 'SNAI'
# counts(dds)[grep(gene2find, rownames(dds)),]
# annot[grep(gene2find, annot$hs.nr.new), ]
# 
# rownames(dds)[10783] = 'ZNF419.ZNF701_AMEX60DD023090'      

dds$condition = droplevels(dds$condition)

x <- DESeqDataSetFromMatrix(counts(dds), DataFrame(design), design = ~ batch + condition )
#dds$batch = droplevels(dds$batch)

dds = estimateSizeFactors(x)

#sels = unique(c(which((design.matrix$batch == 3 | design.matrix$condition == 'BL_UA_9days'))))
#sels = unique(c(which((design.matrix$batch == 3 | design.matrix$condition == 'BL_UA_9days') & 
#                        design.matrix$SampleID != '136150' & design.matrix$SampleID != '106351')))


vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE, ntop = 3000))
ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 1, nudge_y = 1, size=2.5)

ggsave(paste0(resDir, '/PCA_smartseq2_R10724_R11635_regeneration.pdf'),  width=12, height = 8)

#ggsave(paste0(resDir, '/PCA_smartseq2_regeneration.timepoints_filteredSamples_onlyR10724.pdf'),  width=12, height = 8)

cpm = fpm(dds)
cpm = log2(cpm + 2^-6)

plot.pair.comparison.plot(cpm[, c(1:7)], linear.scale = FALSE)

## batch correction with combat
require(sva)
tmm = cpm
bc = as.factor(design$batch)
mod = model.matrix(~ as.factor(condition), data = design)

# if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
# because there is no better batche other others 
#ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = 'R10724_smartseq2') 
#fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 

source('Functions_atac.R')
make.pca.plots(tmm, ntop = 3000, conds.plot = 'all')
ggsave(paste0(resDir, "/smartseq2_R10724_R11635_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)


make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
ggsave(paste0(resDir, "/smartseq2_R10724_R11635_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)

plot.pair.comparison.plot(fpm.bc[, c(1:6)], linear.scale = FALSE)

#dds$condition = droplevels(dds$condition)
dds$condition <- relevel(dds$condition, ref = "Mature_UA")

dds <- DESeq(dds, test="LRT", reduced = ~ batch, fitType = c("parametric"))
plotDispEsts(dds, ymin = 10^-3)

res <- results(dds)

resultsNames(dds)
#res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")

res = data.frame(res[, c(2, 5, 6)])
colnames(res) = paste0(colnames(res), '_LRT')

res.ii = results(dds, name='condition_BL_UA_5days_vs_Mature_UA', test = 'Wald')
res.ii <- lfcShrink(dds, coef="condition_BL_UA_5days_vs_Mature_UA")
colnames(res.ii) = paste0(colnames(res.ii), "_dpa5.vs.mUA")
res = data.frame(res, res.ii[, c(2, 5, 6)])
 
res.ii = results(dds, name='condition_BL_UA_9days_vs_Mature_UA', test = 'Wald')
res.ii <- lfcShrink(dds, coef="condition_BL_UA_9days_vs_Mature_UA")
colnames(res.ii) = paste0(colnames(res.ii), "_dpa9.vs.mUA")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, name='condition_BL_UA_13days_proximal_vs_Mature_UA', test = 'Wald')
res.ii <- lfcShrink(dds, coef="condition_BL_UA_13days_proximal_vs_Mature_UA")
colnames(res.ii) = paste0(colnames(res.ii), "_dpa13prox.vs.mUA")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, name='condition_BL_UA_13days_distal_vs_Mature_UA', test = 'Wald')
res.ii <- lfcShrink(dds, coef="condition_BL_UA_13days_distal_vs_Mature_UA")
colnames(res.ii) = paste0(colnames(res.ii), "_dpa13dist.vs.mUA")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res = data.frame(fpm.bc, res, stringsAsFactors = FALSE)
saveRDS(res, file = paste0(RdataDir, 'smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked',
                           version.analysis, '.rds'))


res = results(dds, contrast = c('condition', 'BL_UA_13days_distal', 'BL_UA_13days_proximal'), test = 'Wald')
res <- lfcShrink(dds, contrast = c('condition', 'BL_UA_13days_distal', 'BL_UA_13days_proximal'))
colnames(res) = paste0(colnames(res.ii), "_dpa13dist.vs.dpa13prox")
res = data.frame(res)

saveRDS(res, file = paste0(RdataDir, 'smartseq2_regeneration_dpa13dist.vs.dpa13prox.rds'))

##########################################
# visualize the DE genes 
##########################################
source('Functions_atac.R')
res = readRDS(file = paste0(RdataDir, 'smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked',
                            version.analysis, '.rds'))

cpm = res[, c(1,2, 5:12)]
res = res[, -c(1:12)]
#cpm = log2(fpm(dds) + 2^-7)

conds = c("Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal",  "BL_UA_13days_distal")
sample.sels = c();  
cc = c()
sample.means = c()
for(n in 1:length(conds)) 
{
  kk = grep(conds[n], colnames(cpm))
  sample.sels = c(sample.sels, kk)
  cc = c(cc, rep(conds[n], length(kk)))
  if(length(kk)>1) {
    sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
  }else{
    sample.means = cbind(sample.means, cpm[, kk])
  }
  
}  

colnames(sample.means) = conds
cpm = cpm[, sample.sels]

if(saveTable){
  jj = grep('^KDM6|^LRAT_', rownames(sample.means))
  
  write.csv(sample.means[jj, ], file = paste0(tableDir, 'KDM6_LRAT_smartseq2_regeneration.csv'), row.names = TRUE)
  
}
res$log2fc = apply(sample.means, 1, function(x) max(x) - min(x))
res$maxs = apply(sample.means, 1, max)
res$mins = apply(sample.means, 1, min)

### select the siganificant genes and visualize the restuls
require(corrplot)
require(pheatmap)
require(RColorBrewer)

#res = readRDS(file = paste0(RdataDir, 'regeneration_smartseq2_fpm_LRTtest_firstMergedSamples.rds'))
#rgs$gene = sapply(rownames(rgs), function(x){x = unlist(strsplit(as.character(x), '_')); return(x[1])})
#rgs$geneID = sapply(rownames(rgs), function(x){x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})

fdr.cutoff = 0.05
logfc.cutoff = 1

length(which(res$padj_LRT<fdr.cutoff))
length(which(res$padj_LRT<fdr.cutoff & res$log2fc>1))
length(which(res$padj_LRT<fdr.cutoff & res$log2fc>2))
length(which(res$padj_LRT<fdr.cutoff & res$log2fc>1.5))

select = which(
                 (res$padj_dpa5.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_dpa5.vs.mUA) > logfc.cutoff) | 
                 (res$padj_dpa9.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_dpa9.vs.mUA) > logfc.cutoff) |
                 (res$padj_dpa13prox.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_dpa13prox.vs.mUA) > logfc.cutoff) |
                 (res$padj_dpa13dist.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_dpa13dist.vs.mUA) > logfc.cutoff))

cat(length(select), ' DE genes \n')
#gg.select = rownames(res)[select]
#saveRDS(gg.select, file = paste0(RdataDir, 'RRGs_candidates_tempList.rds'))

# double check the developmental gene examples 
dev.example = c('HOXA13', 'HOXA11', 'HOXA9', 'HOXD13','HOXD11', 'HOXD9',
                'SHH', 'FGF8', 'FGF10', 'HAND2', 'BMP4', 'ALX1',
                'ALX4', 'PRRX1', 'GREM1', 'LHX2', 'LHX9', 
                'TBX2', 'TBX4', 'TBX5', 'LMX1B', 'MEIS1', 'MEIS2', 'SALL4')

examples.sel = unique(grep(paste0(dev.example, collapse = '|'), rownames(sample.means)))
test = data.frame(sample.means, DEs = !is.na(match(c(1:nrow(sample.means)), select)))


# yy = cpm[select, ]
yy = sample.means[select, ]

newcc = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
df = as.data.frame(newcc)
colnames(df) = 'condition'
rownames(df) = colnames(yy)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

yy <- t(apply(yy, 1, cal_z_score))


library(dendextend)
nb_clusters = 6
my_hclust_gene <- hclust(dist(yy), method = "complete")

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
my_gene_col <- paste0('G', my_gene_col)

## manual change the cluster names
my_gene_col[which(my_gene_col == 'G4')] = '1'
my_gene_col[which(my_gene_col == 'G6')] = '2'
my_gene_col[which(my_gene_col == 'G1')] = '3'
my_gene_col[which(my_gene_col == 'G2')] = '4'
my_gene_col[which(my_gene_col == 'G5')] = '5'
my_gene_col[which(my_gene_col == 'G3')] = '6'

my_gene_col <- data.frame(cluster =  paste0('G', my_gene_col))
rownames(my_gene_col) = rownames(yy)

#corrplot(cor(yy), method = 'number', type = 'upper', diag = TRUE)
#ggsave(filename = paste0(resDir, '/corrplot_smartseq2_regeneration.pdf'),  width = 10, height = 12)

sample_colors = c('springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2')
names(sample_colors) = newcc

col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
          "#33a02c", "#fb9a99", "#e31a1c",
          "#fdbf6f", "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99", "#b15928")
cluster_col = col3[1:nb_clusters]
names(cluster_col) = paste0('G', c(1:nb_clusters))
annot_colors = list(
  sample = sample_colors,
  cluster = cluster_col)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(yy, annotation_row = my_gene_col, 
         annotation_col = df, show_rownames = FALSE, scale = 'none', 
         color =  colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10), 
         show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = FALSE,  
         clustering_method = 'complete', cutree_rows = nb_clusters, 
         annotation_colors = annot_colors, 
         width = 5, height = 12, 
         clustering_callback = callback,
         treeheight_row = 30,
         filename = paste0(figureDir, '/heatmap_DEgenes_regeneration_fdr.0.01_log2fc.1_RNAseq_filtered.R10724.R11635.pdf')) 

##########################################
# save the gene expression clustesrs
##########################################
Save.regeneration.gene.clusters = FALSE
if(Save.regeneration.gene.clusters){
  geneClusters = my_gene_col
  geneClusters$groups = NA
  geneClusters$groups[which(geneClusters$cluster == 'G4')] = 'earlyTransient'
  geneClusters$groups[which(geneClusters$cluster == 'G5')] = 'earlyContinue'
  geneClusters$groups[which(geneClusters$cluster == 'G6')] = 'lateResp'
  
  mm = match(rownames(geneClusters), rownames(sample.means))
  
  geneClusters = data.frame(geneClusters, sample.means[mm, ], stringsAsFactors = FALSE)
  colnames(geneClusters)[-c(1:2)] = newcc
  
  mm = match(rownames(geneClusters), rownames(res))
  geneClusters = data.frame(geneClusters, res[mm, ], stringsAsFactors = FALSE)
  
  geneClusters = geneClusters[order(-geneClusters$log2fc), ]
  
  saveRDS(geneClusters, file = paste0(RdataDir, 'regeneration_geneClusters.rds'))
  
  mm = match(rownames(sample.means), rownames(res))
  colnames(sample.means) = newcc
  
  xx = data.frame(sample.means, res[mm, ], stringsAsFactors = FALSE)
  
  xx$cluster = NA
  mm = match(rownames(xx), rownames(geneClusters))
  xx$cluster[which(!is.na(mm))] = geneClusters$cluster[mm[!is.na(mm)]]
  
  saveRDS(xx, file = paste0(RdataDir, 'regeneration_dynamicGeneClusters_allGenes.rds'))
  
  geneClusters = readRDS(file = paste0(RdataDir, 'regeneration_geneClusters.rds'))
  
  xx = geneClusters[!is.na(geneClusters$groups), ]
  
  #gaps_row =  gaps.row, 
  #filename = paste0(saveDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_rmPeaks.head.pdf'), 
  #width = 6, height = 12)
  
  # gaps.col = c(3, 6)
  # pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
  #          color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
  #          show_colnames = FALSE,
  #          scale = 'row',
  #          cluster_cols=FALSE, annotation_col=df,
  #          annotation_colors = annot_colors,
  #          width = 6, height = 12, 
  #          filename = paste0(figureDir, '/heatmap_DEgenes_regeneration_fdr.0.01_log2fc.1_RNAseq_filtered.R10724.R11635.pdf')) 
  
}

##########################################
# highlight TF, eps and other 
##########################################
# yy = cpm[select, ]
ggs = rownames(yy)
ggs = sapply(ggs, function(x) unlist(strsplit(as.character(x), '_'))[1])
#rownames(yy) = ggs

print(intersect(ggs, tfs))
print(intersect(ggs, sps))
print(intersect(ggs, eps))
print(intersect(ggs, rbp))

for(subg in c('tfs', 'eps', 'sps', 'rbp'))
{
  
  # subg = 'rbp'
  mm = eval(parse(text = paste0('match(ggs, unique(', subg, '))')))
  yy1 = yy[unique(c(which(!is.na(mm)))), ]
  cat(nrow(yy1), ' ', subg, ' found \n')
  
  h = nrow(yy1)*0.1
  
  callback = function(hc, mat){
    sv = svd(t(mat))$v[,2]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  pheatmap(yy1, cluster_rows=TRUE, show_rownames=TRUE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
           show_colnames = FALSE,
           scale = 'row',
           cluster_cols=FALSE, annotation_col=df,
           annotation_colors = annot_colors,
           width = 5, height = h, 
           clustering_callback = callback,
           treeheight_row = 15,
           annotation_legend = FALSE,
           
           filename = paste0(figureDir, '/heatmap_DEgenes_regeneration_fdr.0.01_log2fc.2_smartseq2_', subg, '.pdf'))
  
  #write.table(yy, file = paste0(resDir, '/DEtfs_mUA_regeneration_dev.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
  
}


##########################################
# check the dev-specific genes and mature-specific genes during the regeneration 
##########################################
#rgs = res[select, ] 
dgs = readRDS(file = paste0("../results/Rxxxx_R10723_R11637_R12810_atac/Rdata", 
                            '/pseudoBulk_scRNAcellPooling_FluidigmC1_stage40.44.mUA_2500DEgenes.edgR.Rdata'))
dgs$geneID = sapply(rownames(dgs), function(x){x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})

select = match(dgs$geneID[which(dgs$DEgene == 'matureGene')], rownames(res))

yy = cpm[select[!is.na(select)], ]
df = as.data.frame(cc)
colnames(df) = 'condition'
rownames(df) = colnames(yy)


sample_colors = c('springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2')
names(sample_colors) = conds
annot_colors = list(segments = sample_colors)

pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(16), 
         show_colnames = FALSE,
         scale = 'none',
         cluster_cols=FALSE, annotation_col=df,
         annotation_colors = annot_colors,
         width = 6, height = 8, 
         filename = paste0(figureDir, '/heatmap_DEgenes_regeneration_fdr.0.01_log2fc.1_RNAseq_filtered.R10724_matureSpecificGenes.pdf')) 

select = match(dgs$geneID[which(dgs$DEgene == 'devGene')], rownames(res))

yy = cpm[select[!is.na(select)], ]
pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(16), 
         show_colnames = FALSE,
         scale = 'none',
         cluster_cols=FALSE, annotation_col=df,
         annotation_colors = annot_colors,
         width = 6, height = 8, 
         filename = paste0(figureDir, '/heatmap_DEgenes_regeneration_fdr.0.01_log2fc.1_RNAseq_filtered.R10724_devSpecificGenes.pdf')) 


########################################################
########################################################
# Section : regeneration samples by pooling the scRNA-seq data from Gerber et al. 2018
# 
########################################################
########################################################

##########################################
# reload the metadata and counts 
##########################################
load(file = paste0("../results/Rxxxx_R10723_R11637_R12810_atac/Rdata", 
                   '/Gerber_2018_Fluidigm_C1.Rdata'))

counts = scRNAseq.counts[, -1]
rownames(counts) = scRNAseq.counts$ENSEMBL_ID
rm(scRNAseq.counts)

# convert gene names to gene symbols
mm = match(rownames(counts), annot$geneID)
ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
rownames(counts)[!is.na(mm)] = ggs[!is.na(mm)]

accs = sapply(colnames(counts), function(x) unlist(strsplit(as.character(x), '_'))[1])

mm = match(metadata$metadata.run_accession, accs)
metadata = metadata[which(!is.na(mm)), ]
counts = counts[, mm[which(!is.na(mm))]]
colnames(counts) = as.character(metadata$metadata.run_accession)

metadata$condition = metadata$timepoint
metadata$sample = as.character(metadata$metadata.run_accession)

#raw = as.matrix(counts[ ,-1])
#rownames(raw) = counts$ENSEMBL_ID

metadata$batch = sapply(metadata$cell_id, function(x) unlist(strsplit(as.character(x), '[.]'))[1])
rownames(metadata) = metadata$cell_id

colnames(counts) = rownames(metadata)

## house-keeping genes annotated
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
hs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
ctl =  controls.tissue$geneIDs[which(controls.tissue$tissues  != 'housekeeping')]

##########################################
# analyze scRNA-seq data with Seurat 
##########################################

##########################################
# pool scRNA-seq counts to have pseudo-bulk and analyze with DESeq2 
##########################################
Pool.scRNAseq.pseudobulk = FALSE
if(Pool.scRNAseq.pseudobulk){
  
  require(DESeq2)
  raw = counts
  rm(counts)
  
  conds = c('0dpa', '3dpa', '5dpa', '8dpa', '11dpa', '18dpa', '1apa',  'Stage40', 'Stage44')
  pseudo = c()
  cc = c()
  plates = c()
  for(n in 1:length(conds))
  {
    cat(n, ' : ', conds[n], ' -- ')
    jj = which(metadata$condition == conds[n])
    
    mm = match(metadata$cell_id[jj], colnames(raw))
    if(length(which(is.na(mm)))>0) {
      cat('some cells lost \n')
    }else{
      cat(length(mm), ' cell found \n')
      
      bcs = unique(metadata$batch[jj])
      for(bc in bcs){
        cat('-- ', bc, '--')
        xx = as.matrix(raw[, jj[which(metadata$batch[jj] == bc)]])
        cat(ncol(xx), ' cells \n')
        xx[is.na(xx)] = 0
        pseudo = cbind(pseudo, apply(xx, 1, sum))
        cc = c(cc, conds[n])
        plates = c(plates, bc)
      }
    }
  }
  
  rownames(pseudo) = rownames(raw)
  design = data.frame(condition = cc, plates = plates, stringsAsFactors = FALSE)
  colnames(pseudo) = paste0(design$condition, '_', design$plates)
  
  dds <- DESeqDataSetFromMatrix(pseudo, DataFrame(design), design = ~ condition)
  
  # saveRDS(dds, file = paste0(Rdata.smartseq2, 'Geber_pooledscRNA_muA_regeneration_dev.rds'))
  
  # filtering with stage40/44 samples
  ss = rowMax(counts(dds))
  
  hist(log10(ss), breaks = 60)
  abline(v = log10(c(20, 50, 100)))
  length(which(ss>20))
  length(which(ss>50))
  length(which(ss>100))
  
  dd0 = dds[which(ss > 50), ]
  dd0 = estimateSizeFactors(dd0)
  
  sizeFactors(dds) = sizeFactors(dd0)
  
  ggs = sapply(rownames(dds), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
  
  mm = match(ggs, hs)
  length(which(!is.na(mm)))
  
  dd0 = dds[which(!is.na(mm)), ]
  dd0 = estimateSizeFactors(dd0)
  
  plot(sizeFactors(dds), sizeFactors(dd0))
  
  #dds = estimateSizeFactors(dds)
  #sizeFactors(dds) = sizeFactors(dd0)
  
  cpm = fpm(dds)
  cpm = data.frame(log2(cpm + 2^-7))
  
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  pca=plotPCA(vsd, intgroup = c('condition', 'plates'), returnData = FALSE)
  print(pca)
  
  pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'plates'), returnData = TRUE, ntop = 3000))
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = plates))  + 
    geom_point(size=3) + 
    geom_text(hjust = 1, nudge_y = 1, size=2.5)
  
  plot(ggp)
  
  ggsave(paste0(resDir, '/Geber_pooledscRNA_PCA_muA_regeneration_regeneration.pdf'),  width=12, height = 8)
  
  
  cc = c("0dpa", "3dpa",  "5dpa", "8dpa", "11dpa","18dpa")
  test = matrix(NA, nrow = nrow(cpm), ncol = length(cc))
  for(n in 1:length(cc))
  {
    jj = which(design$condition == cc[n])
    if(length(jj) == 1) test[,n] = cpm[,jj]
    if(length(jj) >1) test[,n] = apply(cpm[,jj], 1, mean)
  }
  rownames(test) = rownames(cpm)
  colnames(test) = cc
  
  saveRDS(test, file = paste0(RdataDir, 'pooled_scRNAseq_cpm.mean_4GRN.rds'))
  
  ##########################################
  # DE test
  ##########################################
  dds = dds[, which(dds$condition != '1apa')]
  
  #cpm = log2(fpm0[, sels] + 2^-6)
  dds$condition = droplevels(dds$condition)
  dds$condition <- relevel(dds$condition, ref = "0dpa")
  
  dds <- DESeq(dds, test="LRT", reduced=~1, fitType = c("local"))
  
  plotDispEsts(dds, ymin = 10^-3)
  
  res <- results(dds)
  
  resultsNames(dds)
    
  res = data.frame(res[, c(2, 5, 6)])
  colnames(res) = paste0(colnames(res), '_LRT')
  
  res.ii = results(dds, name='condition_3dpa_vs_0dpa', test = 'Wald')
  res.ii <- lfcShrink(dds, coef="condition_3dpa_vs_0dpa")
  colnames(res.ii) = paste0(colnames(res.ii), "_3dpa.vs.mUA")
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  res.ii = results(dds, name='condition_5dpa_vs_0dpa', test = 'Wald')
  res.ii <- lfcShrink(dds, coef="condition_5dpa_vs_0dpa")
  colnames(res.ii) = paste0(colnames(res.ii), "_5dpa.vs.mUA")
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  res.ii = results(dds, name='condition_8dpa_vs_0dpa', test = 'Wald')
  res.ii <- lfcShrink(dds, coef="condition_8dpa_vs_0dpa")
  colnames(res.ii) = paste0(colnames(res.ii), "_8dpa.vs.mUA")
  
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  res.ii = results(dds, name='condition_11dpa_vs_0dpa', test = 'Wald')
  res.ii <- lfcShrink(dds, coef="condition_11dpa_vs_0dpa")
  colnames(res.ii) = paste0(colnames(res.ii), "_11dpa.vs.mUA")
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  res.ii = results(dds, name='condition_Stage40_vs_0dpa', test = 'Wald')
  res.ii <- lfcShrink(dds, coef="condition_Stage40_vs_0dpa")
  colnames(res.ii) = paste0(colnames(res.ii), "_Stage40.vs.mUA")
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  res.ii = results(dds, name='condition_Stage44_vs_0dpa', test = 'Wald')
  res.ii <- lfcShrink(dds, coef="condition_Stage44_vs_0dpa")
  colnames(res.ii) = paste0(colnames(res.ii), "_Stage44.vs.mUA")
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  res.ii = results(dds, contrast = c('condition', '11dpa','Stage44'), test = 'Wald')
  res.ii <- lfcShrink(dds, contrast = c('condition', '11dpa','Stage44'))
  colnames(res.ii) = paste0(colnames(res.ii), "_11dpa.vs.Stage44")
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  res.ii = results(dds, contrast = c('condition', '11dpa','Stage40'), test = 'Wald')
  res.ii <- lfcShrink(dds, contrast = c('condition', '11dpa','Stage40'))
  colnames(res.ii) = paste0(colnames(res.ii), "_11dpa.vs.Stage40")
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  saveRDS(res, file = paste0(RdataDir, 'plated_based_pooling_LRT_Waldtest_lfcShrink_v2.rds'))
  save(res, dds, file = paste0(paste0(RdataDir, 'plated_based_pooling_LRT_Waldtest_lfcShrink_dds_res_v2.Rdata')))
  
  ##########################################
  # check the results
  ##########################################
  load(file = paste0(paste0(RdataDir, 'plated_based_pooling_LRT_Waldtest_lfcShrink_dds_res_v2.Rdata')))
  source('Functions_atac.R')
  
  cpm = fpm(dds)
  
  cpm = log2(fpm(dds) + 2^-4)
  #vsd = varianceStabilizingTransformation(dds, blind = FALSE)
  #cpm = assay(vsd)
  
  conds = c("Stage40", "Stage44", "0dpa", "3dpa",  "5dpa",  "8dpa",  "11dpa", '18dpa')
  sample.sels = c();  
  cc = c()
  sample.means = c()
  for(n in 1:length(conds)) 
  {
    kk = grep(conds[n], colnames(cpm))
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
    if(length(kk)>1) {
      sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
    }else{
      sample.means = cbind(sample.means, cpm[, kk])
    }
    
  }  
  
  colnames(sample.means) = conds
  
  cpm = cpm[, sample.sels]
  
  res$log2fc = apply(sample.means, 1, function(x) max(x) - min(x))
  res$maxs = apply(sample.means, 1, max)
  res$mins = apply(sample.means, 1, min)
  
  all = data.frame(sample.means, res, stringsAsFactors = FALSE)
  saveRDS(all, file = paste0(RdataDir, 'pooled_scRNAseq_mUA_regeneration_dev_LRTtest_DESeq2_sample.means.rds'))
  
  ##########################################
  # visualize the restuls
  ##########################################
  require(corrplot)
  require(pheatmap)
  require(RColorBrewer)
  
    
  ## select the significant changing genes
  fdr.cutoff = 0.1
  logfc.cutoff = 1
  
  length(which(res$padj_LRT<fdr.cutoff))
  length(which(res$padj_LRT<fdr.cutoff & res$log2fc>1))
  length(which(res$padj_LRT<fdr.cutoff & res$log2fc>2))
  length(which(res$padj_LRT<fdr.cutoff & res$log2fc>1.5))
  
  select = which( 
                 (res$padj_5dpa.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_5dpa.vs.mUA) > logfc.cutoff) | 
                 (res$padj_8dpa.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_8dpa.vs.mUA) > logfc.cutoff) |
                 (res$padj_11dpa.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_11dpa.vs.mUA) > logfc.cutoff) |
                 (res$padj_Stage40.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_Stage40.vs.mUA) > logfc.cutoff) | 
                 (res$padj_Stage44.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_Stage44.vs.mUA) > logfc.cutoff) )
  cat(length(select), 'DE genes found !\n')
  
  ggs = sapply(rownames(res), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
  
  length(intersect(hs, ggs[select]))
  
  yy = sample.means[select, ]
  
  
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
  library(dendextend)
  yy <- t(apply(yy, 1, cal_z_score))
  
  nb_clusters = 8
  my_hclust_gene <- hclust(dist(yy), method = "complete")
  my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
  res$clusters = NA; 
  res$clusters[match(names(my_gene_col), rownames(res))] = my_gene_col # save the cluster index in ther res
  
  my_gene_col <- data.frame(cluster =  paste0('cluster_', my_gene_col))
  rownames(my_gene_col) = rownames(yy)
  
  df = as.data.frame(conds)
  colnames(df) = 'condition'
  rownames(df) = colnames(yy)
  
  sample_colors = c('deepskyblue', 'deepskyblue3', 
                    'springgreen4', 
                    'springgreen', 'springgreen2', 'springgreen3', 'chartreuse', 'chartreuse4')[c(1:length(conds))]
  
  names(sample_colors) = conds
  #annot_colors = list(condition = sample_colors)
  
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
            "#33a02c", "#fb9a99", "#e31a1c",
            "#fdbf6f", "#ff7f00", "#cab2d6",
            "#6a3d9a", "#ffff99", "#b15928")
  cluster_col = col3[1:nb_clusters]
  names(cluster_col) = paste0('cluster_', c(1:nb_clusters))
  annot_colors = list(
    condition = sample_colors,
    cluster = cluster_col)
  
  gaps.col = c(2, 3)
  pheatmap(yy,  annotation_row = my_gene_col, 
           annotation_col=df,
           cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
           show_colnames = FALSE,
           scale = 'none', gaps_col = gaps.col,
           clustering_method = 'complete', cutree_rows = nb_clusters, 
           cluster_cols=FALSE, 
           annotation_colors = annot_colors,
           width = 6, height = 12, 
           filename = paste0(resDir, '/heatmap_DEgenes_regeneration_pooled.scRNAseq.pdf')) 
  
  pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           annotation_colors = annot_colors,
           width = 6, height = 12, 
           filename = paste0(resDir, '/heatmap_DEgenes_regeneration_pooled.scRNAseq_log2Exp.pdf')) 
  
  ##########################################
  # rough clustering for genes 
  ##########################################
  
  
  ##########################################
  # highlight TF, eps and other 
  ##########################################
  ggs = rownames(yy)
  ggs = sapply(ggs, function(x) unlist(strsplit(as.character(x), '_'))[1])
  #rownames(yy) = ggs
  
  print(intersect(ggs, tfs))
  print(intersect(ggs, sps))
  print(intersect(ggs, eps))
  print(intersect(ggs, rbp))
  
  for(subg in c('tfs', 'eps', 'sps', 'rbp'))
  {
    
    #subg = 'rbp'
    mm = eval(parse(text = paste0('match(ggs, unique(', subg, '))')))
    
    yy1 = yy[unique(c(which(!is.na(mm)))), ]
    
    h = nrow(yy1)*0.1 
    pheatmap(yy1, cluster_rows=TRUE, show_rownames=TRUE, fontsize_row = 5,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
             show_colnames = FALSE,
             scale = 'row',
             cluster_cols=FALSE, annotation_col=df,
             annotation_colors = annot_colors,
             width = 8, height = h, 
             filename = paste0(resDir, '/heatmap_DEgenes_regeneration_pooled.scRNAseq_', subg, '.pdf'))
    
    #write.table(yy, file = paste0(resDir, '/DEtfs_mUA_regeneration_dev.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
    
  }
  
}


########################################################
########################################################
# Section : highlight the comparison between Dev and Mature
# 
########################################################
########################################################
##########################################
# first import Tobias' scRNA-seq data 
##########################################
require(DESeq2)
library(ggrepel)
library(dplyr)
library(tibble)
library("cowplot")

# GO term annotation
limb.go = read.delim(file = '../data/limb_development_GO_term_summary_20220218_162116.txt', header = FALSE)
limb.go = unique(limb.go$V2[-1])
limb.go = toupper(limb.go)

##########################################
# use load DE genes from DESeq2 
##########################################
aa = readRDS(file = paste0(RdataDir, 'pooled_scRNAseq_mUA_regeneration_dev_LRTtest_DESeq2_sample.means.rds'))
res = aa[, -c(1:8)]
fpm = aa[, c(1:8)]

fpm = data.frame(fpm)
fpm$log2FoldChange_Stage40.vs.mUA = res$log2FoldChange_Stage40.vs.mUA
fpm$log2FoldChange_Stage44.vs.mUA = res$log2FoldChange_Stage44.vs.mUA
## house-keeping genes annotated
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
hs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
ctl =  controls.tissue$geneIDs[which(controls.tissue$tissues  != 'housekeeping')]

fdr.cutoff = 0.1
logfc.cutoff = 1

select = which( 
  #(res$padj_5dpa.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_5dpa.vs.mUA) > logfc.cutoff) | 
  #  (res$padj_8dpa.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_8dpa.vs.mUA) > logfc.cutoff) |
  # (res$padj_11dpa.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_11dpa.vs.mUA) > logfc.cutoff) |
    (res$padj_Stage40.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_Stage40.vs.mUA) > logfc.cutoff) | 
    (res$padj_Stage44.vs.mUA < fdr.cutoff & abs(res$log2FoldChange_Stage44.vs.mUA) > logfc.cutoff) )
cat(length(select), 'DE genes found !\n')


fpm$DEgene = NA
fpm$DEgene[select] = '1'

#fpm$genetype[which(fpm$genetype == 'dev.lowlyExp')] = 'Expr.E40.44'

fpm$gene =  sapply(rownames(fpm), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[1])})
fpm$geneID = sapply(rownames(fpm), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
fpm$genetype[!is.na(match(fpm$geneID, hs))] = 'house.keeping'

dev.example = c('HOXA13', 'HOXA11', 'HOXA9', 'HOXD13','HOXD11', 'HOXD9',
                'SHH', 'FGF8', 'FGF10', 'HAND2', 'BMP4', 'ALX1',
                'ALX4', 'PRRX1', 'GREM1', 'LHX2', 'LHX9', 
                'TBX2', 'TBX4', 'TBX5', 'LMX1', 'MEIS1', 'MEIS2', 'SALL4', 'IRX3', 'IRX5')
limb.go = unique(c(limb.go, dev.example))


fpm$limb.go = NA
fpm$limb.go[!is.na(match(fpm$gene, limb.go))] = '1'


fpm$mins = apply(fpm[, c(1, 3)], 1, min)
fpm$maxs = apply(fpm[, c(1, 3)], 1, max)

fpm = fpm[which(fpm$maxs > 0), ]

#fpm$genetype[!is.na(fpm$DEgene)] = 'DE.gene'
#fpm$genetype[which(fpm$genetype == 'devGene')] = 'Expr.E40.44'

fpm$dev = apply(fpm[, c(1, 2)], 1, mean)
fpm$mUA = fpm$X0dpa 

examples.sel = unique(grep(paste0(dev.example, collapse = '|') ,fpm$gene))

fpm$DEgene[which(!is.na(fpm$DEgene) & (fpm$log2FoldChange_Stage40.vs.mUA > 0 |fpm$log2FoldChange_Stage44.vs.mUA >0))] = 'DevGene'
fpm$DEgene[which(!is.na(fpm$DEgene) & (fpm$log2FoldChange_Stage40.vs.mUA < 0|fpm$log2FoldChange_Stage44.vs.mUA <0 ))] = 'MatureGene'

fpm$genetype[which(!is.na(fpm$DEgene))] = fpm$DEgene[which(!is.na(fpm$DEgene))]
fpm$genetype[is.na(fpm$genetype)] = 'others'

scatterplot = ggplot(fpm, aes(x = dev, y = mUA, color = genetype, label = gene)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values=c('red', "darkgray", 'forestgreen',  "orange",   'black')) + 
  #geom_point(data=fpm[which(!is.na(fpm$limb.go)), ], aes(x=Stage40, y=mUA), size=2) + 
  geom_text_repel(data= fpm[examples.sel, ], size = 4.0, color = 'darkblue') + 
  #geom_hline(yintercept=2.0, colour = "darkgray") + 
  #geom_vline(xintercept = 2.0, colour = "darkgray")
  geom_abline(slope = 1,  intercept = 0, colour = 'cyan3') +
  theme_classic() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position=c(0.8, 0.2),
        plot.margin = margin()
        #legend.key.size = unit(1, 'cm')
        #legend.key.width= unit(1, 'cm')
  ) + guides(colour = guide_legend(override.aes = list(size=4)))


# Define marginal histogram
marginal_distribution <- function(x, var, col = 'gray') {
  ggplot(x, aes_string(x = var)) +
    geom_histogram(bins = 40, alpha = 0.4, position = "identity") +
    geom_density(alpha = 0.4, size = 0.1) +
    scale_color_manual(values=col)+
    scale_fill_manual(values=col) +
    guides(fill = FALSE) +
    theme_void() +
    theme(plot.margin = margin())
}

# Set up marginal histograms
x_hist <- marginal_distribution(fpm, 'dev', 'red')
y_hist <- marginal_distribution(fpm, "mUA", 'forestgreen') +
  coord_flip()

# Align histograms with scatterplot
aligned_x_hist <- align_plots(x_hist, scatterplot, align = "v")[[1]]
aligned_y_hist <- align_plots(y_hist, scatterplot, align = "h")[[1]]

# Arrange plots
plot_grid(
  aligned_x_hist
  , NULL
  , scatterplot
  , aligned_y_hist
  , ncol = 2
  , nrow = 2
  , rel_heights = c(0.2, 1)
  , rel_widths = c(1, 0.2)
)

ggsave(paste0(figureDir, "Stage40.44_vs_mUA_devGenes_comparisons_Gerber2018.pseudobulk.pdf"), width=12, height = 10)

