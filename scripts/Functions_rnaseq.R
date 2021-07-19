##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: collection of functions used for RNA-seq analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Jul 19 10:42:46 2021
##########################################################################
##########################################################################
Refine.ax.gene.annot = function(annot)
{
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

Run.QC.for.RNA.replicates = function(design, raw)
{
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

save.scalingFactors.for.deeptools = function(dds)
{
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


########################################################
########################################################
# Section : try to solve the LA and Hand sample swapping
# 
########################################################
########################################################
PCA.and.Pairwise.Comparison.mature.LA.Hand = function(dds, design.matrix)
{
  # PCA 
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  pca=plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = FALSE)
  print(pca)
  ggp = ggplot(data=pca2save[which(pca2save$batch==4), ], aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 0.7, nudge_y = 1, size=2.5)
  
  plot(ggp) + ggsave(paste0(resDir, "/PCAplot_batch4.pdf"), width=12, height = 8)
  
  
  # pairwise comparisons
  fpm = fpm(dds, robust = TRUE)
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
  
  
 
}



