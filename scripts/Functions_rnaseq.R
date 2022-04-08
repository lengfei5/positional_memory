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

##########################################
# find limb fibrobalst-expressed gene and non-expressed genes in dev, mature and regeneration
# mature and regeneration samples with smart-seq2 data
# dev samples will use the scRNA-seq data
##########################################
Define_limb_expressed_genes = function()
{
  # library('rtracklayer')
  # library(GenomicRanges)
  # library('GenomicFeatures')
  # import full gtf file with hox patch
  # gtf = import(paste0(annotDir, 'AmexT_v47_Hox.patch.gtf'))
  
  ## all transcripts and geneID from v47
  annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
  annot = readRDS(paste0(annotDir, 
                         'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
  
  load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
  hs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
  ctl =  controls.tissue$geneIDs[which(controls.tissue$tissues  != 'housekeeping')]
  
  ## genes with gene symbols 
  genes = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  
  ## define the matrix for the result
  
  ## import the feature count output file and extract the gene names in the hox patch and also the gene length information
  lls = read.table(paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/rnaseq_using/', 
                                 'featurecounts_Q10/Mature_UA_136149_featureCounts.txt'),
                   header = TRUE, sep = '\t')
  
  ggs = unique(lls$Geneid)
  ggs = data.frame(geneID = ggs, ORF.type = annot$ORF.type_gtf.transcript[match(ggs, annot$geneID)],
                   length = lls$Length[match(ggs, lls$Geneid)], stringsAsFactors = FALSE)
  
  ggs$geneSymbol = NA
  mm = match(ggs$geneID, genes$geneID)
  ggs$geneSymbol[!is.na(mm)] = genes$gene.symbol.toUse[mm[!is.na(mm)]]
  ggs$names = NA
  kk = which(is.na(ggs$geneSymbol))
  ggs$names[kk] = as.character(ggs$geneID[kk])
  ggs$names[-kk] = paste0(as.character(ggs$geneSymbol[-kk]), '_', as.character(ggs$geneID[-kk]))
  
  
  ### check first the microarray data for mature samples
  Define_limb.expressing.genes_microarray = FALSE
  if(Define_limb.expressing.genes_microarray){
    # load microarray analysis results: log2 signal (res), comparison (fit2), and raw data (raw)
    load(file = paste0("../results/microarray/Rdata/", 
                       'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_DEpval.Rdata'))
    cpm = as.matrix(res[, c(1:9)])
    
    conds = c("mUA", 'mLA', 'mHand') 
    sample.means = c()
    for(n in 1:length(conds)) 
    {
      kk = grep(conds[n], colnames(cpm))
      #sample.sels = c(sample.sels, kk)
      #cc = c(cc, rep(conds[n], length(kk)))
      if(length(kk)>1) {
        sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
      }else{
        sample.means = cbind(sample.means, cpm[, kk])
      }
      
    }  
    
    colnames(sample.means) = conds
    sample.means[grep('SHH|PAX6|CD68|VIM|FOXA2', rownames(sample.means)), ]
    
    ids = sapply(rownames(sample.means), function(x){x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
    
    
    for(n in 1:ncol(sample.means))
    {
      x = sample.means[,n]
      ii.control = which(!is.na(match(ids, ctl)))
      
      cutoff =  quantile(x[ii.control], 0.8)
      cat(n,  ' -- cutoff used ', cutoff, '--', length(which(x>cutoff)),  ' expressed gene \n')
      
      sample.means[,n] = sample.means[,n] > cutoff  
      
    }
    sample.means[grep('SHH|PAX6|CD68|VIM|FOXA2', rownames(sample.means)), ]
    
    nbs = apply(sample.means, 1, sum)
    
    nbs = nbs[which(nbs>0)]
    
    xx = sapply(names(nbs), function(x){x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
    
    ggs$expr.mature.microarray = NA
    ggs$expr.mature.microarray[!is.na(match(ggs$geneID, xx))] = 1
    
  }
  
  ###### check the smartseq2 mature samples
  # load dds normalized object and annotations
  Rdata.smartseq2 = "../results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/"
  #load(file = paste0(Rdata.smartseq2, 'dds_design.matrix_all29smartseq2_beforeFiltering.Rdata'))
  load(file = paste0(RdataDir, 'dds_design.matrix_all29smartseq2_beforeFiltering.Rdata'))

  # select mature samples
  sels = intersect(which(design.matrix$batch == 4), grep('Mature', design.matrix$condition))
  dds = dds[, sels]
  
  cpm = fpm(dds)
  cpm = cpm[, grep('HEAD', colnames(cpm), invert = TRUE)]
  
  conds = c("Mature_UA", 'Mature_LA', 'Mature_Hand') 
  sample.means = c()
  for(n in 1:length(conds)) 
  {
    kk = grep(conds[n], colnames(cpm))
    #sample.sels = c(sample.sels, kk)
    #cc = c(cc, rep(conds[n], length(kk)))
    if(length(kk)>1) {
      sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
    }else{
      sample.means = cbind(sample.means, cpm[, kk])
    }
    
  }  
  colnames(sample.means) = conds
  
  log2(sample.means[grep('SHH|PAX6|CD68|VIM|FOXA2', rownames(sample.means)), ])
  
  xx = sample.means
  for(n in 1:ncol(sample.means))
  {
    # n = 1
    x = sample.means[,n]
    x = x[which(x>0)]
    x = log2(x)
    ids = sapply(names(x), function(x){x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
    ii.control = which(!is.na(match(ids, ctl)))
    cutoff =  quantile(x[ii.control], 0.75)
    
    hist(x, breaks = 100);abline(v = log2(cutoff), lwd = 2.0, col = 'red')
    cat(n,  ' -- cutoff used ', log2(cutoff), '--', length(which(x>log2(cutoff))),  ' expressed gene \n')
    
    xx[,n] = sample.means[,n] > cutoff  
    
  }
  
  xx[grep('SHH|PAX6|CD68|VIM|FOXA2', rownames(xx)), ]
  nbs = apply(xx, 1, sum)
  nbs = nbs[which(nbs>0)]
  
  geneID.expr = sapply(names(nbs), function(x){x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
  
  ggs$expr.mature.smartseq = NA
  ggs$expr.mature.smartseq[!is.na(match(ggs$geneID, geneID.expr))] = 1
  
  
  
  
  Test.multiple.gaussian = FALSE
  if(Test.multiple.gaussian){
    for(n in 1:ncol(cpm)){
      # n = 1
      x = cpm[,n]
      x = x[x>0]
      #x = log2(x)
      #x = x[which(x>-2)]
      cat(length(x), '\n')
      ii.control = which(!is.na(match(ids, ctl)))
      
      gg.expr = unique(c(gg.expr, names(x)))
      
      PLOT.EM = FALSE
      if(PLOT.EM){
        require(mixtools)
        fit <- normalmixEM(x, lambda = c(0.4, 0.6), 
                           mu = c(5, 10), sigma = NULL, 
                           #mean.constr=m.constr, 
                           #sd.constr=sigma.constr, 
                           k=2, maxrestarts=20, maxit = 1500)
        #plot(fit, density=TRUE)
        lambda = fit$lambda;
        par1 = fit$mu
        par2 = fit$sigma
        
        xfit<-seq(min(x),max(x),length=100)
        yfit1<-dnorm(xfit, par1[1],par2[1])*lambda[1]
        yfit2<-dnorm(xfit, par1[2],par2[2])*lambda[2]
        
        #test = c(test, percent*lambda[c(1:2)])
        hist(x, breaks=60,  freq=FALSE)
        lines(xfit, yfit1, col='green', lwd=2.0)
        lines(xfit, yfit2, col='green', lwd=2.)
        abline(v = c(7:9), lwd = 2.0, col = 'red')
        
      }
      
    }
  }
  
  
  
}


##########################################
# process human chromatin remodelers and RBPs  
##########################################
process_database_huamn_chromatin.remodelers_RBP = function()
{
  eps = read.delim(file = 
      paste0('/Volumes/groups/tanaka/People/current/jiwang/annotations/mouse/mouse_chromatin.remodeler/Epifactors_database.csv'))
  eps = unique(as.character(eps$HGNC.approved.symbol))
  
  rbp = read.csv(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/annotations/mouse/', 
                               'mouse_RBPDB/RBPDB_v1.3.1_human_2012-11-21_CSV/RBPDB_v1.3.1_proteins_human_2012-11-21.csv'), 
                   header = FALSE)
  
  rbp = unique(as.character(rbp$V5))
  
  saveRDS(eps, file = paste0('../data/human_chromatin_remodelers_Epifactors.database.rds'))
  saveRDS(rbp, file = paste0('../data/human_RBPs_rbpdb.rds'))
  
}


##########################################
# some QC functions
##########################################
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


########################################################
# saturation curve from rseqc
# the r code from rseqc output
########################################################
RNAseq.sequence.saturation.test = function()
{
  rseqc.file = list.files('../Data/R10724_rnaseq/saturation_rseqc', pattern = 'junctionSaturation_plot.r', 
                          full.names = TRUE)
  library(stringr)
  
  yy = c()
  for(n in 1:length(rseqc.file))
  {
    cat(n, '\n')
    xx = read.delim(rseqc.file[n])
    xx = xx[grep('y=', xx[, 1 ]), ]
    #xx = gsub('y=c', '', xx)
    # Get the parenthesis and what is inside
    k <- str_extract_all(xx, "\\([^()]+\\)")[[1]]
    # Remove parenthesis
    k <- substring(k, 2, nchar(k)-1)
    #k = gsub('["]', '', k)
    k = as.numeric(unlist(strsplit(as.character(k), ',')))
    yy = rbind(yy, k)
  }
  
  rownames(yy) = gsub('_junction.junctionSaturation_plot.r', '', basename(rseqc.file))
  
  pdfname = paste0(resDir, '/saturation_curve_rseqc_knownJunctions.pdf')
  pdf(pdfname, width = 16, height = 8)
  par(cex = 1.0, las = 1,  mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0), tcl = -0.3)
  
  yy = yy[which(rownames(yy) != '136150_TTAACCTTCGAGGCCAGACA_HNF3KDSXY_3_20201223B_20201223'), ]
  span = 0.75
  # saturation curve with nb of peaks
  xlims = c(0, 120)
  ylims = range(yy/10^3)
  frac = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)/100
  
  library(RColorBrewer)
  cols = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(nrow(yy))
  plot(0, 0, xlim = xlims, ylim = ylims, type ='n', xlab = 'nb of TOTAL reads (Million)', 
       ylab = 'nb of known junctions (K)', main = paste0('saturation curve from rseqc'))
  abline(v = c(20, 30, 40, 50), col = 'blue', lwd = 1.0, lty =2)
  
  #legend('topleft', legend = sample.uniq, col = cols, bty = 'n', lwd = 2.0, cex = 0.7)
  
  for(n in 1:nrow(yy))
  {
    # n = 1
    cat(n, '\n')
    
    kk = which(design$sample == rownames(yy)[n])
    
    satt = data.frame(nb.reads = design$total.reads[kk]*frac/10^6, nb.junctions = yy[n, ]/10^3)
    
    points(satt[,1], satt[,2], type= 'p', col = cols[n])
    loessMod <- loess(nb.junctions ~ nb.reads, data=satt, span=span)
    smoothed <- predict(loessMod)
    lines(smoothed, x=satt$nb.reads, col=cols[n], lwd = 3.0)
    
    text(satt[nrow(satt), 1], smoothed[length(smoothed)], labels = paste0(design$fileName[kk], '_', design$sampleID[kk]), 
         cex = 0.7, pos = 4, offset = 0.2)
    
  }
  
  dev.off()
  
}

##########################################
# save scaling factors in DESeq2 for deeptools
##########################################
save.scalingFactors.for.deeptools = function(dds, saveDir)
{
  saveDir = dataDir
  #ss = colSums(counts(dds))
  #plot(ss[jj], (design.matrix$alignment.rate*design.matrix$trimmed.reads)[jj])
  
  #reads.mapped = design.matrix$trimmed.reads*design.matrix$alignment.rate/100
  reads.mapped = design.matrix$mapped
  xx = data.frame(sampleID = design.matrix$SampleID,  
                  scalingFactor = reads.mapped/(design.matrix$sizefactor*median(reads.mapped)),
                  stringsAsFactors = FALSE)
  xx[grep('161514|161517', xx[1, ]), ]
  
  write.table(xx, file = paste0(saveDir, '/DESeq2_scalingFactor_forDeeptools.txt'), sep = '\t',
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

########################################################
########################################################
# Section : Update or add newly sequenced samples
# 
########################################################
########################################################
Update.samples.160343.160344 = function()
{
  load(file=paste0(RdataDir, 'Design_stats_readCounts_', version.analysis, '.Rdata'))
  
  Correct_swapped.mature.samples = TRUE
  if(Correct_swapped.mature.samples){
    jj = which(design$SampleID == '161517')  
    design$conditions[jj] = 'Mature_Hand'
    design$conds[jj] = 'Mature_Hand_161517.batch4'
    
    jj = which(design$SampleID == '161518')  
    design$conditions[jj] = 'Mature_LA'
    design$conds[jj] = 'Mature_LA_161518.batch4'
    
    colnames(all)[-1] = design$conds
    
  }
  
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11635_rnaseq_redo.no.arab/'
  
  # update counts
  countFile = paste0(dataDir, 'nf_out/featureCounts/merged_gene_counts.txt')
  counts = read.delim(file = countFile, sep = '\t', header = TRUE)
  counts = counts[match(all$gene, counts$ENSEMBL_ID), ]
  
  par(mfrow = c(1, 2))
  plot(counts[, 2], all[, which(design$SampleID == '106343' & design$batch == '4')], cex = 0.6, log ='xy');
  abline(0, 1, lwd = 2.0, col = 'red')
  
  plot(counts[, 3], all[, which(design$SampleID == '106344' & design$batch == '4')], cex = 0.6, log ='xy');
  abline(0, 1, lwd = 2.0, col = 'red')
  
  all[, which(design$SampleID == '106343' & design$batch == '4')] = counts[, 2]
  all[, which(design$SampleID == '106344' & design$batch == '4')] = counts[, 3]
  
  # update the sequence statistics
  stat = read.delim(file = paste0(dataDir, '/nf_out/MultiQC/multiqc_data/multiqc_general_stats.txt'))
  
  #xx = data.frame(design, matrix(NA, nrow = nrow(design), ncol = 9), stringsAsFactors = FALSE)
  #colnames(xx)[7:15] = colnames(stat0)[-c(1:3)]
  #mm = match(stat0$sampleID, design$SampleID)
  #xx[mm, c(7:15)] = stat0[, -c(1:3)]
  
  kk = c(which(design$SampleID == '106343' & design$batch == '4'), which(design$SampleID == '106344' & design$batch == '4'))
  design$pct.duplication[kk] = stat$FastQC_mqc.generalstats.fastqc.percent_duplicates
  design$pct.GC[kk] = stat$FastQC_mqc.generalstats.fastqc.percent_gc
  design$readType[kk] = 'SE'
  design$avg.seq.length[kk] = stat$FastQC_mqc.generalstats.fastqc.avg_sequence_length
  design$total.reads[kk] = stat$FastQC_mqc.generalstats.fastqc.total_sequences
  design$trimmed.reads[kk] = NA
  design$assigned.reads[kk] = stat$featureCounts_mqc.generalstats.featurecounts.Assigned
  design$alignment.rate[kk] = stat$HISAT2_mqc.generalstats.hisat2.overall_alignment_rate
  design$pct.assigned.to.features[kk] = stat$featureCounts_mqc.generalstats.featurecounts.percent_assigned
 
  save(design, all, file=paste0(RdataDir, 'Design_stats_readCounts_updatedResequenced', version.analysis, '.Rdata'))
   
}


########################################################
########################################################
# Section : double check the replicates of regeneration samples
# 
########################################################
########################################################
Check.R10724.136150.sample.Quality = function()
{
  
  ##########################################
  # first check all regeneration samples from R10724
  # secondly compare the resequenced sample 136150; resequened 136150 is totally different from the its technical replicate 
  # and also biological replicate (what ???)
  # last add the polyA samples and collect all regeneration RNAseq samples
  
  ##########################################
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/rnaseq_using/regeneration_RNAseq/'
  
  design = read.xlsx(paste0(dataDir, 'sample_Infos.xlsx'))
  #colnames(design) = c('sampleID', 'fileName')
  #design = rbind(design, c('136150s', 'BL_UA_5days'))
  
  #design$SampleID[13] = '13615x'
  
  xlist = list.files(path=paste0(dataDir, 'gene_counts'),
                     pattern = "*featureCounts.txt$", full.names = TRUE) ## list of data set to merge
  
  all = cat.countTable(xlist, countsfrom = 'featureCounts')
  
  #colnames(all)[grep('136150s', colnames(all))] = '13615x.txt'
  colnames(all) =  gsub('[#]', '_', colnames(all))
  
  colnames(design)[1] = 'SampleID'
  
  
  counts = process.countTable(all=all, design = design[, c(1,2)], merge.technicalRep.sameID = TRUE)
  
  design$condition = design$fileName
  design$fileName = paste0(design$fileName, '_', design$SampleID)
  design$conds = design$condition
  
  design$batch = NA
  colnames(counts)[-1] = design$fileName
  
  design$batch = as.factor(design$batch)
  
  annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  all = counts
  
  mm = match(all$gene, annot$geneID)
  ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
  all$gene[!is.na(mm)] = ggs[!is.na(mm)]
  
  raw = as.matrix(all[, -1])
  rownames(raw) = all$gene
  
  design$batch = paste0(design$request, '_', design$protocol)
  
  
  dds <- DESeqDataSetFromMatrix(raw, DataFrame(design), design = ~ condition)
  
  save(design, dds, file = paste0(RdataDir, 'design_dds_all_regeneration.samples_allBatches.Rdata'))
  
  ##########################################
  # reload the data and select the samples needed 
  ##########################################
  load(file = paste0(RdataDir, 'design_dds_all_regeneration.samples_allBatches.Rdata'))
  
  sels = which(design$SampleID != '13615x')
  design = design[sels, ]
  raw = raw[, sels]
  
  design$protocol = gsub(' ', '', design$protocol)
  
  design$batch = paste0(design$request, '_', design$protocol)
  
  table(design$condition, design$batch)
  
  sels = which(design$batch == 'R10724_smartseq2' | design$batch == 'R9707_polyA')
  design = design[sels, ]
  
  dds = dds[,sels]
  
  dds$condition = droplevels(dds$condition)
  
  
  ss = rowSums(counts(dds))
  
  
  dds = dds[which(ss>20), ]
  
  dds = estimateSizeFactors(dds)
  
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  pca=plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = FALSE)
  print(pca)
  
  pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE, ntop = 1000))
  ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 1, nudge_y = 1, size=2.5)
  
  
  ##########################################
  # test batch correction 
  ##########################################
  source('Functions_atac.R')
  library(edgeR)
  require("sva")
  require(limma)
  
  d <- DGEList(counts=counts(dds), group=dds$condition)
  tmm <- calcNormFactors(d, method='TMM')
  tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  
  bc = as.factor(design$batch)
  mod = model.matrix(~ as.factor(condition), data = design)
  
  # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
  # because there is no better batche other others 
  #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)   
  
  table(design$conds, design$batch)
  
  fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
  
  #design.tokeep<-model.matrix(~ 0 + conds,  data = design.sels)
  #cpm.bc = limma::removeBatchEffect(tmm, batch = bc, design = design.tokeep)
  # plot(fpm.bc[,1], tmm[, 1]);abline(0, 1, lwd = 2.0, col = 'red')
  make.pca.plots(tmm, ntop = 3000, conds.plot = 'all')
  ggsave(paste0(resDir, "/matureSamples_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)
  
  
  make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
  ggsave(paste0(resDir, "/matureSamples_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)
  
  
  ggsave(paste0(resDir, '/PCA_smartseq2_regeneration.timepoints_R10724_R11635.pdf'),  width=12, height = 8)
  
  #cpm = fpm(dds)
  #cpm = log2(fpm(dds) + 2^-7)
  cpm = fpm.bc
  
  conds = c("Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal",  "BL_UA_13days_distal")
  
  sample.sels = c();  cc = c()
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
  
  require(corrplot)
  require(pheatmap)
  require(RColorBrewer)
  
  
  gg.select = readRDS(file = paste0(RdataDir, 'RRGs_candidates_tempList.rds'))
  select = match(gg.select, rownames(cpm))
  select = select[!is.na(select)]
  
  yy = cpm[select, sample.sels]
  df = as.data.frame(cc)
  colnames(df) = 'condition'
  rownames(df) = colnames(yy)
  
  corrplot(cor(yy), method = 'number', type = 'upper', diag = TRUE)
  
  plot.pair.comparison.plot(yy[, c(8:9)], linear.scale = FALSE)
  #ggsave(filename = paste0(resDir, '/corrplot_smartseq2_regeneration.pdf'),  width = 10, height = 12)
  
  sample_colors = c('springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2')
  names(sample_colors) = conds
  annot_colors = list(segments = sample_colors)
  
  pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
           show_colnames = FALSE,
           scale = 'row',
           cluster_cols=FALSE, annotation_col=df,
           annotation_colors = annot_colors,
           width = 8, height = 12, 
           filename = paste0(resDir, '/heatmap_DEgenes_regeneration_fdr.0.01_log2fc.2_RNAseq_polyA.smartseq2.pdf')) 
  
  ##########################################
  # sample means 
  ##########################################
  yy = sample.means[select, ]
  df = as.data.frame(conds)
  colnames(df) = 'condition'
  rownames(df) = colnames(yy)
  
  #corrplot(cor(yy), method = 'number', type = 'upper', diag = TRUE)
  #plot.pair.comparison.plot(yy[, c(8:9)], linear.scale = FALSE)
  #ggsave(filename = paste0(resDir, '/corrplot_smartseq2_regeneration.pdf'),  width = 10, height = 12)
  
  sample_colors = c('springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2')
  names(sample_colors) = conds
  annot_colors = list(segments = sample_colors)
  
  pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
           show_colnames = FALSE,
           scale = 'row',
           cluster_cols=FALSE, annotation_col=df,
           annotation_colors = annot_colors,
           width = 8, height = 12, 
           filename = paste0(resDir, '/heatmap_DEgenes_regeneration_fdr.0.01_log2fc.2_RNAseq_polyA.smartseq2_means.pdf'))
  
} 

edgeR.DE.test.without.replicates = function()
{
  edgeR.DE.test = FALSE
  if(edgeR.DE.test){
    ##########################################
    # DE genes using edgeR
    ##########################################
    require(edgeR)
    group = dds$condition
    y <- DGEList(counts=counts(dds), group=group)
    y <- calcNormFactors(y, method = 'TMM')
    
    cpm = cpm(y, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.5)
    
    plot.pair.comparison.plot(cpm, linear.scale = FALSE)
    
    y1 = y
    y1$samples$group <- 1
    
    index.hs = which(fpm$genetype == 'hs')
    y0 <- estimateDisp(y1[index.hs,], trend = 'none', tagwise=FALSE)
    
    y$common.dispersion <- y0$common.dispersion
    
    fit <- glmFit(y)
    lrt <- glmLRT(fit, coef=2:6)
    #lrt3 = glmLRT(fit, coef = 3)
    
    lrt = as.data.frame(topTags(lrt, adjust.method = "BH", n = nrow(y)))
    #lrt3 = as.data.frame(topTags(lrt3, adjust.method = "BH", n = nrow(y)))
    
    lrt$logfc = apply(lrt[, c(1:6)], 1, function(x) max(x) - min(x)) 
    #select = which(lrt$FDR< fdr.cutoff & lrt$logfc>5)
    
    select = c(1:3000)
    
    yy = cpm[match(rownames(lrt)[select], rownames(cpm)), ]
    df = as.data.frame(conds)
    colnames(df) = 'condition'
    rownames(df) = colnames(yy)
    
    
    #corrplot(cor(yy), method = 'number', type = 'upper', diag = TRUE)
    #ggsave(filename = paste0(resDir, '/corrplot_smartseq2_regeneration.pdf'),  width = 10, height = 12)
    
    sample_colors = c('springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2', 'darkgray', 'red')
    names(sample_colors) = conds
    annot_colors = list(samples = sample_colors)
    
    pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
             show_colnames = FALSE,
             scale = 'row',
             cluster_cols=FALSE, annotation_col=df,
             annotation_colors = annot_colors,
             width = 6, height = 12, 
             filename = paste0(resDir, '/heatmap_DEgenes_regeneration_pseudoBulk.pdf')) 
    
    ##########################################
    # highlight TF, eps and other 
    ##########################################
    ggs = rownames(yy)
    ggs = sapply(ggs, function(x) unlist(strsplit(as.character(x), '_'))[1])
    
    print(intersect(ggs, tfs))
    print(intersect(ggs, sps))
    print(intersect(ggs, eps))
    print(intersect(ggs, rbp))
    
    for(subg in c('tfs', 'eps', 'sps', 'rbp'))
    {
      
      mm = eval(parse(text = paste0('match(ggs, unique(', subg, '))')))
      
      yy1 = yy[unique(c(which(!is.na(mm)))), ]
      
      pheatmap(yy1, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
               show_colnames = FALSE,
               scale = 'row',
               cluster_cols=FALSE, annotation_col=df,
               annotation_colors = annot_colors,
               width = 8, height = 8, 
               filename = paste0(figureDir, '/heatmap_DEgenes_regeneration_fdr.0.01_log2fc.2_smartseq2_', subg, '.pdf'))
      
      #write.table(yy, file = paste0(resDir, '/DEtfs_mUA_regeneration_dev.txt'), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
      
    }
    
  }
  
}  
  
########################################################
########################################################
# a data exploration function
# compare embryo stage 44 proximal with distal and BL.UA.day13 proximal and distal 
# trying to check the what deposoits the positional information 
########################################################
########################################################
Compare_stage44.proximal.distal_BL.UA.day13.promximal.distal = function()
{
  
  # load dds normalized object and annotations
  load(file = paste0(RdataDir, 'RNAseq_design_dds.object.Rdata'))
  annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  tfs = readRDS(file = paste0('../results/motif_analysis/TFs_annot/curated_human_TFs_Lambert.rds'))
  sps = readRDS(file = '~/workspace/imp/organoid_patterning/results/Rdata/curated_signaling.pathways_gene.list_v2.rds')
  eps = readRDS(file = paste0('../data/human_chromatin_remodelers_Epifactors.database.rds'))
  rbp = readRDS(file = paste0('../data/human_RBPs_rbpdb.rds'))
  tfs = unique(tfs$`HGNC symbol`)
  sps = toupper(unique(sps$gene))
  
  saveTables = TRUE
  
  # select mature samples
  sels = grep('BL_UA_13days|Embryo_Stage46', design.matrix$condition)
  dds = dds[, sels]
  
  cpm = fpm(dds)
  cpm = log2(cpm + 2^-4)
  
  dds$condition = droplevels(dds$condition)
  dds <- estimateDispersions(dds, fitType = 'parametric')
  plotDispEsts(dds, ymin = 10^-3); abline(h = 0.1, col = 'blue', lwd = 2.0)
  
  dds = nbinomWaldTest(dds, betaPrior = TRUE)
  resultsNames(dds)
  
  ##########################################
  # heatmap to visualize all positional genes and regulators (TFs and SPs)
  ##########################################
  library("pheatmap")
  Comparison =  "BL.UA.distal.proximal" 
  
  if(Comparison == 'BL.UA.distal.proximal'){
    # BL_UA_day13 distal vs. proximal
    res.ii = results(dds, contrast=c("condition", 'BL_UA_13days_distal', 'BL_UA_13days_proximal'), alpha = 0.1)
    #colnames(res.ii) = paste0(colnames(res.ii), "_BL_UA_13days_distal.vs.proximal")
    res = data.frame(res.ii[, c(2, 5, 6)])
    
  }else{
    res.ii = results(dds, contrast=c("condition", 'Embryo_Stage46_distal', 'Embryo_Stage46_proximal'), alpha = 0.1)
    #colnames(res.ii) = paste0(colnames(res.ii), "_Embryo_Stage46_distal.vs.proximal")
    res = data.frame(res.ii[, c(2, 5, 6)])
    
    
    #pval.cutoff = 0.001
    #res$pval.max = apply(-log10(res[, grep('pvalue_', colnames(res))]), 1, max)
    #res$logFC.max = apply((res[, grep('log2FoldChange_', colnames(res))]), 1, function(x) return(x[which(abs(x)==max(abs(x)))][1]))
  }
  
  o1 = order(res$pvalue)
  res = res[o1, ]
  cpm = cpm[o1, ]
  
  ggs = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})
  select = which(res$padj < 0.1 & abs(res$log2FoldChange)> 1)
  ggs = ggs[select]
  
  print(intersect(ggs, tfs))
  print(intersect(ggs, sps))
  print(intersect(ggs, eps))
  print(intersect(ggs, rbp))
  
  #select = which(res$pvalue_Hand.vs.LA < pval.cutoff | res$pvalue_Hand.vs.UA < pval.cutoff | res$pvalue_LA.vs.UA < pval.cutoff)
  #cat(length(select), ' positional genes found \n')
  
  # fdr.cutoff = 0.2
  # select = which(res$padj_Hand.vs.LA < fdr.cutoff | 
  #                  res$padj_Hand.vs.UA < fdr.cutoff |
  #                  res$padj_LA.vs.UA < fdr.cutoff)
  # cat(length(select), ' positional genes found \n')
  
  # df <- as.data.frame(colData(dds)[,c("condition", 'batch')])
  # o1 = c(grep('UA', df$condition), grep('LA', df$condition), grep('Hand', df$condition))
  # 
  # yy = cpm[select, o1]
  # ss = apply(as.matrix(yy), 1, mean)
  # 
  # #yy = yy[which(ss>-2), ]
  # pheatmap(yy, cluster_rows=TRUE, show_rownames=TRUE, fontsize_row = 6,
  #          show_colnames = FALSE,
  #          scale = 'row',
  #          cluster_cols=FALSE, annotation_col=df[o1, ], 
  #          width = 8, height = 20, filename = paste0(figureDir, '/heatmap_DEgenes_matureSamples_pvalmax.3_log2FC.0_.pdf'))
  
  
  library(ggrepel)
  library(dplyr)
  library(tibble)
  res$gene = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})
  res$padj = -log10(res$padj)
  
  examples.sel = unique(res$gene[which(res$padj >10 | abs(res$log2FoldChange) > 5)])
  #examples.sel = unique(intersect(res$, tfs)
  
  
  ggplot(data=res, aes(x=log2FoldChange, y=padj, label = gene)) +
    geom_point(size = 1) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    #geom_text_repel(data=subset(res, log2FoldChange > 2), size = 4) +
    geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
                     size = 3) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(0), col="red") +
    geom_hline(yintercept=5, col="red") + 
    ggsave(paste0(figureDir, "VolcanoPlot_logFC_padj_EmbryoStage.proximal.vs.distal.pdf"), width=12, height = 8)
  
  # # narrow down to TFs and SPs
  # ggs = sapply(rownames(yy), function(x) unlist(strsplit(as.character(x), '_'))[1])
  # mm = match(ggs, unique(c(tfs[, 3], toupper(sps$gene))))
  # yy1 = yy[!is.na(mm), ]
  # 
  # pheatmap(yy1, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE,
  #          scale = 'row',
  #          cluster_cols=FALSE, annotation_col=df[o1, ], fontsize_row = 10, 
  #          width = 8, height = 4,
  #          filename = paste0(resDir, '/heatmap_DE.tfs.sps_mature_fdr.0.1.pdf')) 
  
  
  if(saveTables){
    xx = data.frame(gene = rownames(yy), yy, res[match(rownames(yy), rownames(res)), ], stringsAsFactors = FALSE)
    write.csv(xx, file = paste0(resDir, '/position_dependent_genes_from_matureSamples_RNAseq_fdr.0.1.csv'), 
              quote = FALSE, col.names = TRUE, row.names = FALSE)
    
  }
  
}

##########################################
# edgeR DE test for pooled Dev and Mature from scRNA-seq using house-keeping gene to test 
##########################################
Dev.vs.Mature_edgeR.DE.test.without.replicates = function()
{
  ## house-keeping genes annotated
  annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
  load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
  hs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
  ctl =  controls.tissue$geneIDs[which(controls.tissue$tissues  != 'housekeeping')]
  load(file = paste0(RdataDir, '/pseudoBulk_scRNAcellPooling_FluidigmC1_stage40.44.mUA_dev_geneSelection.Rdata'))
  require(edgeR)
  group = dds$condition
  y <- DGEList(counts=counts(dds), group=group)
  y <- calcNormFactors(y)
  
  y1 = y
  y1$samples$group <- 1
  #group1 = y1$samples$group
  #design <- model.matrix(~group1)
  
  index.hs = which(fpm$genetype == 'hs')
  y0 <- estimateDisp(y1[index.hs,], trend = 'none', tagwise=FALSE)
  
  y$common.dispersion <- y0$common.dispersion
  
  fit <- glmFit(y)
  lrt2 <- glmLRT(fit, coef = 2)
  lrt3 = glmLRT(fit, coef = 3)
  
  lrt2 = as.data.frame(topTags(lrt2, adjust.method = "BH", n = nrow(y)))
  lrt3 = as.data.frame(topTags(lrt3, adjust.method = "BH", n = nrow(y)))
  
  DE.genes = unique(c(rownames(lrt2)[which(lrt2$FDR<0.01 & abs(lrt2$logFC) > 1)],
                      rownames(lrt3)[which(lrt3$FDR <0.01 & abs(lrt3$logFC) >1)]))
  
  cat(length(DE.genes), ' DE genes found \n')
  
}





