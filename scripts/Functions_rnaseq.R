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
  xx.m = xx.m[c(1:50), ] 
  
  #write.table(xx.m, file = paste0(resDir, '/top100_LA.Hand.specific.genes.fromATACseq.matureSamples_usedforHand.LA.sampleSwapping.txt'), 
  #             sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
  
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

##########################################
# Compare microarray data with Akane's RNA-seq data
##########################################
compare.mature.samples.RNAseq.with.microarray = function()
{
  # import mature samples in RNA-seq 
  load(file = paste0("../results/rnaseq_Rxxxx.old_R10724_R161513/Rdata/", 
                     'RNAseq_design_dds.object.Rdata'))
  cpm = fpm(dds, robust = TRUE)
  
  kk = which(design.matrix$batch == 4 & (
    design.matrix$conditions == 'Mature_UA'| design.matrix$conditions == 'Mature_LA'| 
      design.matrix$conditions == 'Mature_Hand'))
  #vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  cpm = cpm[, kk]
  dm = design.matrix[kk, ]
  
  # import microarray data and microarray data is already in log2 scale; position-dependent genes in the matrix
  load(file = paste0("../results/microarray/Rdata/", 
                     'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_DEpval.Rdata'))
  #annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
  #                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  #mm = match(rownames(res), annot$geneID)
  #mm = match(all$gene, annot$geneID)
  #ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
  #rownames(res)[!is.na(mm)] = ggs[!is.na(mm)]
  
  # select the same set of genes from RNAseq and microarray 
  ggs = unique(intersect(rownames(res), rownames(cpm)))
  
  cpm = cpm[match(ggs, rownames(cpm)), ]
  cpm = log2(cpm + 2^-6)
  res = res[match(ggs, rownames(res)), ]
  #res = res-(median(res) - median(cpm))
  
  ss = apply(as.matrix(cpm), 1, function(x) length(which(x > 0)))
  sels = which(ss>2)
  res = res[sels, ]
  cpm = cpm[sels, ]
  
  candidates = res[, c(10:12)]
  res = res[, c(1:9)]
  
  PCA.plot.Using.DEgenes.from.microarray = FALSE
  if(PCA.plot.Using.DEgenes.from.microarray)
  {
    ##########################################
    #  PCA plot using the identified postional genes from microarray
    ##########################################
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
    pca=plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = FALSE)
    print(pca)
    pca2save = as.data.frame(plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = TRUE))
    
    jj = which(pca2save$batch == 4 & 
                 (pca2save$conditions == 'Mature_UA'| pca2save$conditions == 'Mature_LA'| pca2save$conditions == 'Mature_Hand'))
    ggplot(data=pca2save[jj, ], aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
      geom_point(size=3) + 
      geom_text(hjust = 0.7, nudge_y = 1, size=2.5)
    
    
    sels = apply(as.matrix(candidates), 1, function(x) length(which(x < 0.0001)))
    cands = candidates[which(sels > 0), ]
    
    cat(nrow(cands), ' positional genes selected \n')
    
    mm = match(rownames(vsd), rownames(cands))
    pca2save = as.data.frame(plotPCA(vsd[which(!is.na(mm)), ], intgroup = c('conditions', 'batch'), returnData = TRUE))
    
    jj = which(pca2save$batch == 4 & 
                 (pca2save$conditions == 'Mature_UA'| pca2save$conditions == 'Mature_LA'| 
                    pca2save$conditions == 'Mature_Hand'))
    ggplot(data=pca2save[jj, ], aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
      geom_point(size=3) + 
      geom_text(hjust = 0.7, nudge_y = 1, size=2.5)
    
    ##########################################
    # calculate correlation between RNA-seq samples and microarray data using postional genes 
    ##########################################
    avg = data.frame(mUA = apply(res[, grep('mUA', colnames(res))], 1, mean), 
                     mLA = apply(res[, grep('mLA', colnames(res))], 1, mean),
                     mHand = apply(res[, grep('mHand', colnames(res))], 1, mean))
    
    pval.cutoff = 0.00001
    #sels = apply(as.matrix(candidates), 1, function(x) length(which(x < pval.cutoff)))
    #cands = candidates[which(sels > 0), ]
    cands = candidates[which(candidates$mLA.vs.mUA < pval.cutoff | candidates$mHand.vs.mLA < pval.cutoff),]
    cat(nrow(cands), ' positional genes selected \n')
    gene2use = rownames(cands)
    
    # top1 = topTable(fit2, coef = 1, number = 100)
    # top2 = topTable(fit2, coef = 3, number = 100)
    # gene2use = unique(c(rownames(top1)[which(top1$P.Value < pval.cutoff & top1$logFC > 0)], 
    #              rownames(top2)[which(top2$P.Value < pval.cutoff & top2$logFC > 0)]))
    
    mm = match(rownames(cpm), gene2use) 
    mm = which(!is.na(mm))
    xx = matrix(ncol = ncol(avg), nrow = ncol(cpm))
    colnames(xx) = colnames(avg)
    rownames(xx) = colnames(cpm)
    for(n in 1:ncol(xx))
    {
      for(m in 1:nrow(xx)){
        xx[m, n] = cor(cpm[mm, m], avg[mm, n], method = 'pearson')
      }
    }
    
    require(dplyr)
    require(ggplot2)
    library(tidyr)
    data.frame(xx, sample = rownames(xx))%>% as_tibble() %>%  gather(ref, correlation, 1:3) %>%
      ggplot(aes(x = sample, y = correlation, fill = ref)) + 
        geom_bar(stat = "identity", position="dodge")# +
     # theme(axis.text.x = element_text(angle = 90, size = 10))
    
  }
  
  # check PCA without batch correction
  Check.PCA.without.batch.correction = FALSE
  if(Check.PCA.without.batch.correction){
    #xx = cbind(cpm, res)
    xx = cpm
    library(factoextra)
    
    ntop = 1000
   
    xx = as.matrix(xx)
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
    
  }
 
  ## correct batch using UA samples
  require("sva")
  dm = rbind(dm[, c(2, 5)], data.frame(conditions = colnames(res), batch = as.factor(rep(0, ncol(res)))))
  dm$conditions[7:9] = 'Mature_UA'
  #dm$conditions[10:12] = 'Mature_LAx'
  #dm$conditions[13:15] = 'Mature_Handx'
  
  bc = droplevels(dm$batch)
  #bc = levelsdroplevels(bc)
  mod = model.matrix(~ as.factor(conditions), data = dm)
  
  yy = cbind(cpm, res)
  cpm.bc = ComBat(dat=yy, batch=bc, mod=mod, par.prior=TRUE, ref.batch = 4, mean.only = FALSE)    
  
  require(limma)
  cpm.bc = removeBatchEffect(yy, batch = bc)
  
  ntop = 1000
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
  
  
}




