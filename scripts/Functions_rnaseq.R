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
  # PCA with all genes
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  pca=plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = FALSE)
  print(pca)
  pca2save = as.data.frame(plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = TRUE))
  
  jj = which(pca2save$batch == 4 & 
                     (pca2save$conditions == 'Mature_UA'| pca2save$conditions == 'Mature_LA'| pca2save$conditions == 'Mature_Hand'))
  ggp = ggplot(data=pca2save[jj, ], aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 0.7, nudge_y = 1, size=2.5)
  
  plot(ggp) + ggsave(paste0(resDir, "/PCAplot_batch4.pdf"), width=12, height = 8)
  
  
  ## PCA with potenital position-dependent genes from RNA-seq data
  load(file = '../results/Rdata/LA_Hand_potential.specific.genes_fromRNAseq.rds')
  
  xx = data.frame(rbind(xx1, xx2), gene = c(rownames(xx1), rownames(xx2)), stringsAsFactors = FALSE)
  xx = xx[match(unique(xx$gene), xx$gene), ]
  xx = xx[which(xx$mUA < 0 & xx$BL.day5<0 & (xx$BL.day9 > 2 | xx$BL.day13 >2)), ]
  
  
  #xx1 = xx1[which((xx1$BL.day9 - xx1$mUA) >2 & (xx1$BL.day9 - xx1$BL.day5) >2), ]
  mm = match(rownames(vsd), c(rownames(xx1), rownames(xx2)))
    
  pca2save = as.data.frame(plotPCA(vsd[which(!is.na(mm)), ], intgroup = c('conditions', 'batch'), returnData = TRUE))
  
  jj = which(pca2save$batch == 4 & 
               (pca2save$conditions == 'Mature_UA'| pca2save$conditions == 'Mature_LA'| 
                  pca2save$conditions == 'Mature_Hand'))
  ggp = ggplot(data=pca2save[jj, ], aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 0.7, nudge_y = 1, size=2.5)
  
  plot(ggp) + ggsave(paste0(resDir, "/PCA_top60genes_LA.Hand_RNAtimecourse.pdf"), width=12, height = 8)
  
  ## PCA with potenital position-dependent genes from ATAC-seq data
  aa = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                      'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID.rds'))
  fpm = fpm(dds, robust = TRUE)
  
  xx.m = readRDS(file = '../results/Rdata/LA_Hand_potential.specific.genes_fromATACseq.rds')
  mns = apply(as.matrix(xx.m[, c(2:3)]), 1, min)
  xx.m = xx.m[which(mns < 2), ]
  xx.m$vars = apply(as.matrix(xx.m[, c(2:3)]), 1, var)
  xx.m = xx.m[order(-xx.m$vars), ]
  
  #xx.m = xx.m[which(vars > 3), ]
  xx.m = xx.m[c(1:100), ] 
  
  ggs = xx.m$geneId
  ggs = sapply(ggs, function(x) {names = unlist(strsplit(as.character(x), '[|]')); return(names[length(names)])})
 
  mm = match(ggs, aa$transcritID)
  ggs = unique(aa$geneID[mm])
  
  axg = sapply(rownames(vsd), function(x) {names = unlist(strsplit(as.character(x), '_')); return(names[length(names)])})
  
  mm = which(!is.na(match(axg, ggs)))
  pca2save = as.data.frame(plotPCA(vsd[mm, ], intgroup = c('conditions', 'batch'), returnData = TRUE))
  jj = which(pca2save$batch == 4 & 
               (pca2save$conditions == 'Mature_LA'| 
                  pca2save$conditions == 'Mature_Hand'))
  
 
  ggp = ggplot(data=pca2save[jj, ], aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 0.2, nudge_y = 0.05, size=2.5)
  
  plot(ggp) + ggsave(paste0(resDir, "/PCA_top100genes_LA.Hand_ATACseq.pdf"), width=12, height = 8)
  
  cpm = fpm[mm, jj]
  
  plot.pair.comparison.plot(cpm)
  
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

Identify.LA.Hand.specific.genes.from.UA.regeneration.timeCourse = function(dds, design.matrix)
{
  jj = which(design.matrix$batch == 3 & design.matrix$conditions != 'Embryo_Stage40')
  
  ddx = dds[, jj]
  ddx$conditions = droplevels(ddx$conditions)
  
  cpm = fpm(ddx)
  
  xx = data.frame(mUA = apply(cpm[, grep('Mature_UA', colnames(cpm))], 1, mean), 
                  BL.day5 = cpm[, grep('BL_UA_5days', colnames(cpm))], 
                  BL.day9 = apply(cpm[, grep('BL_UA_9days', colnames(cpm))], 1, mean), 
                  BL.day13 = apply(cpm[, grep('BL_UA_13days', colnames(cpm))], 1, mean))
  xx = log2(xx + 2^-6)
  
  ss = apply(as.matrix(xx), 1, max)
  xx = xx[which(ss>1), ]
  
  kk1 = which((xx$BL.day9 - xx$BL.day5) >1 & (xx$BL.day9 - xx$mUA) >1) 
  kk2 = which((xx$BL.day13 - xx$BL.day5) >1 & (xx$BL.day13 - xx$mUA) >1 & (xx$BL.day13 - xx$BL.day9) > 1) 
  
  xx1 = xx[kk1, ]
  xx2 = xx[kk2, ]
  save(xx1, xx2, file = '../results/Rdata/LA_Hand_potential.specific.genes_fromRNAseq.rds')
  
}

