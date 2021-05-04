##########################################################################
##########################################################################
# Project: positional memory 
# Script purpose: 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Mar 13 11:12:07 2021
##########################################################################
##########################################################################

########################################################
########################################################
# Section : identify group of house-keeping genes using 21 axolotl tissues 
# here the general idea is that there are three big groups of genes: 
# tissue-wide (expressed in all tissues), tissue-restricted (expressed in certain tissues not others), 
# tissue-specific (in only one tissue)
# the tissue-wide genes were commonly considered house-keeping genes, e.g. ribosome-related genes
# Regarding our limb system, there will be house-keeping genes, limb-related genes 
# (i.e. limb-specific genes and genes expressed in limb 
# some limited number of tissues)
########################################################
########################################################
find.houseKeepingGenes.limb.related.genes = function()
{
  require(ggplot2)
  require(DESeq2)
  library("dplyr")
  library("ggplot2")
  require("pheatmap")
  
  counts = read.delim(file = '../Data/axolotl_21tissues_rnaseq/nf_out/featureCounts/merged_gene_counts.txt',
                      sep = '\t', header = TRUE)
  
  rownames(counts) = as.character(counts[, 1])  
  counts = as.matrix(counts[, -1])
  design = data.frame(fileName = colnames(counts),  tissue = gsub('_left.sorted_gene.featureCounts.txt', '', colnames(counts)), 
                      stringsAsFactors = FALSE)
 
  dds <- DESeqDataSetFromMatrix(counts, DataFrame(design), design = ~ tissue)
  
  ss0 = colSums(counts(dds))
  ss = rowSums(counts(dds))
  
  hist(log10(ss), breaks = 200, main = 'log10(sum of reads for each gene)')
  
  ## here filter genes whose sums of all samples are < 20 reads  
  cutoff.peak = 20
  cat(length(which(ss > cutoff.peak)), 'genes selected \n')
  abline(v = log10(cutoff.peak), lwd = 2.0, col = 'red')
  
  dds <- dds[ss > cutoff.peak, ]
  
  # normalization and dimensionality reduction
  dds <- estimateSizeFactors(dds)
  fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  
  pca=plotPCA(vsd, intgroup = c('tissue'), returnData = FALSE)
  print(pca)
  
  pca2save = as.data.frame(plotPCA(vsd, intgroup = c('tissue'), returnData = TRUE))
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = tissue, color= tissue))  + 
    geom_point(size=4) + 
    geom_text(hjust = 1.2, nudge_y = 2, size=4)
  
  plot(ggp) + ggsave(paste0(resDir, "/PCAplot_axolotl_21tissues.pdf"), width=12, height = 8)
  
  ##########################################
  # for each tissue: group genes into not.detected, lowly.expressed, expressed 
  # those three group will be categrotelize as 0, 1, 2
  ##########################################
  plot(sizeFactors(dds), colSums(counts(dds))/median(colSums(counts(dds))), log = 'xy')
  
  fpm = fpm(dds, robust = TRUE)
  colnames(fpm) = design$tissue
  fpm = fpm[, grep('pancreas|nCC|oCC', colnames(fpm), invert = TRUE)]
  
  ggs = matrix(NA, nrow = nrow(fpm), ncol = ncol(fpm))
  rownames(ggs) = rownames(fpm)
  colnames(ggs) = colnames(fpm)
  
  # consider house-keeping genes should be on the top 10k at each condition and in each tissue
  rank.hsg = 10000
  
  for(n in 1:ncol(fpm))
  {
    # n = 1
    test = fpm[,n]
    kk = which(fpm[, n] <=10^-6)
    if(length(kk) > 0){
      ggs[kk, n] = 0;
      test = test[-kk]
    }
    
    test = log2(test)
    #hist(test, breaks = 100)
    
    o1 = order(-test)
    test = test[o1]
    hist(test[c(1:rank.hsg)], breaks = 100)
    
    ggs[match(names(test)[1:rank.hsg], rownames(ggs)), n] = 2
    ggs[match(names(test)[c((rank.hsg+1):length(test))], rownames(ggs)), n] = 1
    
  }
  
  # discard genes not detected across all samples  
  ss = apply(ggs, 1, function(x) return(length(which(x > 1)))) # count nb of samples in which genes were in the top 10k
  jj = which(ss>0)
  ss = ss[jj]
  fpm = fpm[jj, ]
  ggs = ggs[jj, ]
  
  hist(ss)
  
  ii.hkg = which(ss >= 16) ## house-keeping highly expressed in >= 16 tissue and conditions
  
  hkg = names(ss)[ii.hkg]
  
  annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
  annot = readRDS(file = paste0(annotDir, 
                'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
  
  hkg.hs = annot$gene.symbol.hs_gtf.gene[match(hkg, annot$geneID)]
  
  hkgs = data.frame(amexID = hkg, hs.symbol = hkg.hs, stringsAsFactors = FALSE)
  rm(list = c('hkg', 'hkg.hs'))
  
  saveRDS(hkgs, file = paste0(annotDir, 'axolotl_house.keep.genes_expressed.top10k.across18tissues.rds'))
  
  ##########################################
  # here we also save some mature tissue-specific genes
  # for the control of limb-related genes 
  ##########################################
  #ss = ss[-ii.hkg]
  #fpm = fpm[-ii.hkg, ]
  #ggs = ggs[-ii.hkg, ]
  ggs = data.frame(ggs)
  fpm = data.frame(fpm)
  
  rank.specific = 3000
  
  ii.liver = which(ss == 1 & ggs$mLiver == 2 & fpm$mLiver >= fpm$mLiver[order(-fpm$mLiver)][rank.specific])
  ii.islet = which(ss == 1 & ggs$Islets_L12905_Track_36524 == 2 & 
                     fpm$Islets_L12905_Track_36524 >= 
                     fpm$Islets_L12905_Track_36524[order(-fpm$Islets_L12905_Track_36524)][rank.specific])
  ii.testis = which(ss == 1 & ggs$mTestis == 2 & fpm$mTestis >= fpm$mTestis[order(-fpm$mTestis)][rank.specific])
  ii.heart = which(ss == 1 & ggs$mHeart == 2 & fpm$mHeart >= fpm$mHeart[order(-fpm$mHeart)][rank.specific])
  ii.spleen = which(ss == 1 & ggs$mSpleen == 2 & fpm$mHeart >= fpm$mHeart[order(-fpm$mHeart)][rank.specific])
  ii.lung = which(ss == 1 & ggs$mLung == 2 & fpm$mLung >= fpm$mLung[order(-fpm$mLung)][rank.specific])
  
  ii.controls = c(ii.hkg, ii.liver, ii.islet, ii.testis, ii.spleen, ii.heart, ii.lung)
  
  pheatmap(log2(fpm[ii.controls, ]+2^-6), cluster_rows=FALSE, show_rownames=FALSE, scale = 'none', show_colnames = TRUE,
           cluster_cols=FALSE)
  
  controls.tissue = data.frame(tissues = c(rep('housekeeping', length(ii.hkg)),
                                    rep('liver', length(ii.liver)), 
                                    rep('islet', length(ii.islet)), 
                                    rep('testis', length(ii.testis)), 
                                    rep('spleen', length(ii.spleen)), 
                                    rep('heart', length(ii.heart)), 
                                    rep('lung', length(ii.lung))),
                        geneIDs = c(names(ss)[c(ii.hkg, ii.liver, ii.islet, ii.testis, ii.spleen, 
                                                ii.heart, ii.lung)]),
                        stringsAsFactors = FALSE)
  
  controls.tissue$gene.hs = annot$gene.symbol.hs_gtf.gene[match(controls.tissue$geneIDs, annot$geneID)]
  
  controls.tissue = data.frame(controls.tissue, stringsAsFactors = FALSE)
  
  save(controls.tissue, 
       file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
  
  
}

processing.development.genes.from.Sergej.paper = function()
{
  ##########################################
  # focus on peaks around the developmental genes used in Sergej's paper 
  # not used here, because most of developmental genes were maintained open
  ##########################################
  Test.developmental.gene.from.Sergej = FALSE
  if(Test.developmental.gene.from.Sergej){
    require(openxlsx)
    library("readxl")
    annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
    gene.mapping = read.delim(paste0(annotDir, 'AmexT_v47.release_rm.contigs_geneID_geneSymbol.mapping.txt'), 
                              sep = '\t', header = TRUE)
    
    Use.Sergej.developmental.genes = TRUE
    if(Use.Sergej.developmental.genes){
      devg = read_excel(paste0(annotDir, 'developmental_geneList.xls'))
      devg = unique(devg$HGNC.symbol) 
    }
    
    mm = match(devg, gene.mapping$gene.symbol.hs)
    
    devg = gene.mapping[mm[which(!is.na(mm))], ]
    
    library("GenomicFeatures")
    library("GenomicRanges")
    
    devg.c = GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
                     ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)))
  }
  
}

########################################################
########################################################
# Section : axolotl promoter regions selection 
# and peak annotation
# 
########################################################
########################################################
select.promoters.regions = function(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', toSave = FALSE)
{
  require(GenomicRanges)
  cat('load gtf annotation \n')
  annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
  #gene.mapping = read.delim(paste0(annotDir, 'AmexT_v47.release_rm.contigs_geneID_geneSymbol.mapping.txt'), 
  #                          sep = '\t', header = TRUE)
  annot = readRDS(paste0(annotDir, 
        'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
  
  gs = readRDS(file = paste0(annotDir, 'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr.rds'))
  #gs$geneSymbol.synteny.ns.nr = gs$gene.evidence.synteny
  cat(length(unique(gs$geneSymbol.synteny.ns.nr)), 'gene symbols \n')
  
  #jj = which(is.na(gs$geneSymbol.synteny.ns.nr) & !is.na(gs$gene.evidence.same.gene.symbol.hs.nr))
  #gs$geneSymbol.synteny.ns.nr[jj] = gs$gene.evidence.same.gene.symbol.hs.nr[jj]
  #cat(length(unique(gs$geneSymbol.synteny.ns.nr)), 'gene symbols \n')
  #saveRDS(gs, file = paste0(annotDir, 'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr.rds'))
  #index.geneSymbols = which(!is.na(gene.mapping$gene.symbol.hs)|!is.na(gene.mapping$gene.symbol.nr))
  
  # select putative full length transcripts and C-terminal transcripts
  
  if(ORF.type.gtf == 'Putative'){
    cat('select putative full length transcripts \n')
    index.sel = which(annot$ORF.type_gtf.transcript == 'Putative')
  }else{
    cat('select putative full length transcripts and C-terminal transcripts \n')
    index.sel = which(annot$ORF.type_gtf.transcript == 'Putative' | annot$ORF.type_gtf.transcript == 'N-terminal')
  }
  index.sel = index.sel[grep('^chr', annot$chr_transcript[index.sel])] # rm genes in cotigs
  
  tss = data.frame(annot[index.sel, c(25, 26, 27, 1, 28, 7)], stringsAsFactors = FALSE)
  colnames(tss) = c('chr', 'start', 'end', 'transcriptID', 'strand', 'geneID')
  tss$chr = as.character(tss$chr)
  tss$start = as.integer(as.character(tss$start))
  tss$end = as.integer(as.character(tss$end))
  
  jj = which(tss$strand == '+')
  tss$end[jj] = tss$start[jj] +1
  tss$start[-jj] = tss$end[-jj] - 1
  
  tss$geneSymbol = tss$geneID
  mm = match(tss$geneSymbol, gs$geneID)
  tss$geneSymbol[which(!is.na(mm))] = gs$geneSymbol.synteny.ns.nr[mm[which(!is.na(mm))]]
  tss$geneSymbol[which(is.na(tss$geneSymbol))] = tss$geneID[which(is.na(tss$geneSymbol))]
  
  colnames(tss) = c('seqnames', 'start', 'end', 'transcriptID', 'strand', 'geneID', 'geneSymbol')
  
  tss = makeGRangesFromDataFrame(tss, seqnames.field=c("seqnames"),
                                 start.field="start", end.field=c("end"), strand.field="strand", keep.extra.columns = TRUE)
  
  promoters = promoters(tss, upstream=upstream, downstream=downstream, use.names=TRUE)
  
  #promoters.reduced = reduce(promoters, drop.empty.ranges = TRUE, ignore.strand = TRUE)
  if(toSave){
    saveRDS(promoters, 
            file = paste0(RdataDir, 
                          '/axolotl_promoters_Granges_AmexT_v47_putative.full.length_N.terminal_upstream2kb.downstream2kb.rds'))
  }
  
  return(promoters)  
  
}


annotatePeak.curateAxolotl = function(peaks) 
{
  # peaks = bb
  peaks$promoters = NA
  pp = makeGRangesFromDataFrame(peaks, seqnames.field = 'chr', 
                                start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  xx = overlapsAny(pp, promoter, ignore.strand = TRUE)
  
  peaks$promoters[xx] = 'promters'
  
  ##########################################
  # promoter annotation
  ##########################################
  Annotate.promoters = FALSE
  if(Annotate.promoters){
    select.promoters.regions();
  }
  
  ##########################################
  # try to customize annotation for ChIPpeakannot 
  ##########################################
  require(GenomicFeatures)
  require(ChIPpeakAnno)
  amex = makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))
  
  saveRDS(amex, file = paste0(annotDir, 'TxDb_ax6_UCSC_2021_01_26_genes.putative.full.length.rds'))

  peakAnnots = annotatePeak(p1, TxDb=amex, tssRegion = c(-2000, 2000))
  xx = data.frame(peakAnnots)
  
  ### plot the partitions of Peaks into different genomics features and distance to TSS
  pdf("PLOTS/Peaks_distributions_distance_TSS.pdf", width = 10, height = 6)
  #par(mfrow=c(1,2))
  print(plotAnnoBar(peakAnnots))
  print(plotDistToTSS(peakAnnots))
  dev.off()
  
  return(peaks)
  
}

########################################################
########################################################
# Section : promoter peak analysis
# 1) Compare the gene expression and promoter peaks
# 2) gene enrichment analysis to check if those promoters were limb-development and limb related 
########################################################
########################################################
DoubleCheck.promoter.peaks.enrichment = function(fpm)
{
  require(ChIPpeakAnno)
  require(ChIPseeker)
  
  ##########################################
  # only input is the normalized peak signals
  # first make granges for peaks
  ##########################################
  pp = data.frame(t(sapply(rownames(dds), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  rownames(pp) = gsub('_', '-', rownames(pp))
  
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
 
  pp.annots = as.data.frame(annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript'))
  
  rownames(fpm) = names(pp)
  
  ##########################################
  # select first the peaks at promoters in dev and mature samples 
  ##########################################
  promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative')
  
  sel.promoter = which(overlapsAny(pp, promoters, ignore.strand = TRUE) == TRUE)
  
  conds.sel = c('Embryo_Stage40_13', 'Embryo_Stage40_93', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
    'Mature_UA_13', 'Mature_UA_74938|Mature_UA_102655', 'Mature_LA', 'Mature_Hand')
  
  sample.sel = c()
  ii.gaps = c()
  for(n in 1:length(conds.sel)) {
    c = conds.sel[n]
    sample.sel = c(sample.sel, grep(c, colnames(fpm)))
    if(n == 1) {
      ii.gaps = c(ii.gaps, length(grep(c, colnames(fpm))))
    }else{
      if(n != length(conds.sel)) ii.gaps = c(ii.gaps, (ii.gaps[n-1] + length(grep(c, colnames(fpm)))))
    }
  }
  
  df <- data.frame(colData(dds)[,c("conds", 'batch')])[sample.sel, ]
  
  keep =fpm[sel.promoter, sample.sel]
  
  filtering.with.signal = FALSE
  if(filtering.with.signal){
    ss = apply((keep), 1, max)
    
    max.cutoff = 3
    hist(ss, breaks = 100)
    abline(v = max.cutoff, col='red', lwd =2.0)
    cat(length(which(ss>max.cutoff)), 'peak selected\n')
    
    keep = keep[which(ss > max.cutoff), ]
  }
  
  keep = as.matrix(keep)
  
  norm_minmax <- function(x, a = -1, b = 1){
    #a + (x- min(x)) /(max(x)-min(x))*(b-a)
    x - max(x)
  }
  
  #keep.xx = as.data.frame(t(apply(keep, 1, norm_minmax)))
  breaks = c(-6, -4, 0, 2, 4, 8, 12);
  cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks)-1)
  
  pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE, 
           color = cols,
           cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, 
           breaks = breaks)
  
  ##########################################
  # find the genes for those promoters
  ##########################################
  #pp = data.frame(t(sapply(rownames(keep), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  #rownames(pp) = gsub('_', '-', rownames(pp))
  #rownames(keep) = rownames(pp)
  
  #pp$strand = '*'
  #pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
  #                              start.field="X2", end.field="X3", strand.field="strand")
  
  #promoters = readRDS(file = 
  #paste0('../results/R10723_Rxxxx_atacseq_bowtie2.newParam_mtDNA_picardrmdup_20210208/Rdata', 
  # '/axolotl_promoters_Granges_AmexT_v47_putative.full.length_N.terminal_upstream2kb.downstream2kb.rds'))
  pp = pp[match(rownames(keep), names(pp))]
  
  pp.overlap = GenomicRanges::findOverlaps(pp, promoters)
  pp.overlap = data.frame(pp.overlap)
  
  mapping = data.frame(matrix(NA, ncol = 4, nrow = nrow(keep)))
  rownames(mapping) = rownames(keep)
  colnames(mapping) = c('peak.index', 'promoter.index', 'gene', 'nb.genes')
  mapping$peak.index = c(1:nrow(keep))
  
  for(n in 1:nrow(mapping))
  {
    # n = 2
    ii = pp.overlap$subjectHits[which(pp.overlap$queryHits == n)]
    if(length(ii)>0){
      mapping$promoter.index[n] = paste0(pp.overlap$subjectHits[ii], collapse = ';')
      gg = unique(promoters$geneID[ii])
      mapping$gene[n] = paste0(gg, collapse = ';')
      mapping$nb.genes[n] = length(gg)
    }
  }
  
  ##########################################
  # associate the promoter peaks with the RNA-seq data  
  ##########################################
  Integrate.gene.expression.data = FALSE
  if(Integrate.gene.expression.data){
    load.RNAseq.data = FALSE
    if(load.RNAseq.data){
      load(file=paste0('../results/rnaseq_RNAseqSamples_all/Rdata/Design_stats_readCounts__rnaseq_RNAseqSamples_all_20210301.Rdata'))
      design$batch = as.factor(design$batch)
      require(ggplot2)
      require(DESeq2)
      library("dplyr")
      library("ggplot2")
      
      raw = as.matrix(all[, -1])
      rownames(raw) = all$gene
      
      #sels = which(design$batch != 1)
      sels = which(design$assigned.reads > 10^6)
      
      design.matrix = design[sels, ]
      raw = raw[, sels]
      
      dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ conditions)
      
      ss = rowMax(counts(dds))
      
      hist(log10(ss), breaks = 200, main = 'log10(sum of reads for each gene)')
      cutoff.peak = 10
      abline(v = log10(cutoff.peak), lwd = 2.0, col = 'red')
      
      cat(length(which(ss > cutoff.peak)), 'selected genes with read cutoff', cutoff.peak,  '\n')
      
      dds <- dds[ss > cutoff.peak, ]
      
      # normalization and dimensionality reduction
      dds <- estimateSizeFactors(dds)
      
      plot(sizeFactors(dds), colSums(counts(dds))/median(colSums(counts(dds))), log = 'xy')
      
      cpm = fpm(dds, robust = TRUE)
      
      par(mfrow = c(1:2))
      plot(cpm[, c(1, 16)], log = 'xy', cex = 0.2)
      abline(0, 1, lwd = 2.0, col = 'red')
      plot(cpm[, c(11, 13)], log = 'xy', cex = 0.2)
      abline(0, 1, lwd = 2.0, col = 'red')
      
      
      save(cpm, design.matrix, file = paste0(RdataDir, '/RNAseq_allSamples_cpm_designMatrix.Rdata'))
      
    }
    
    load(file = paste0(RdataDir, '/RNAseq_allSamples_cpm_designMatrix.Rdata'))
    
    select.rnaSamples = FALSE
    if(select.rnaSamples){
      rna.sel = c(which(design.matrix$conditions == 'Embryo_Stage40' & design.matrix$batch == '3'), 
                  #which(design.matrix$conditions == 'Embryo_Stage40' & design.matrix$batch == '2'),
                  #which(design.matrix$conditions == 'Embryo_Stage46_proximal' & design.matrix$batch == '1'),
                  #which(design.matrix$conditions == 'Embryo_Stage46_distal' & design.matrix$batch == '1'),
                  which(design.matrix$conditions == 'Mature_UA' & design.matrix$batch == '3')
                  #which(design.matrix$conditions == 'Mature_UA' & design.matrix$batch == '2'),
                  #which(design.matrix$conditions == 'BL_UA_5days' & design.matrix$batch == '3' & design.matrix$SampleID != '136150'),
                  #which(design.matrix$conditions == 'BL_UA_9days' & design.matrix$batch == '3'),
                  #which(design.matrix$conditions == 'BL_UA_13days_proximal' & design.matrix$batch == '3'), 
                  #which(design.matrix$conditions == 'BL_UA_13days_distal' & design.matrix$batch == '3')
      )
      
      df.rna <- data.frame(conds = design.matrix$conditions[rna.sel], batch = design.matrix$batch[rna.sel])
      df.rna$batch = as.factor(df.rna$batch)
      df.rna$conds = as.factor(df.rna$conds)
      cpm = cpm[, rna.sel]
      rownames(df.rna) = colnames(cpm)
      
    }
    
    #cpm.xx = as.matrix(log2(cpm + 2^-6))
    mapping$rna.index = NA
    mapping$rna.index = match(mapping$gene, rownames(cpm))
    
    for(n in 1:nrow(mapping))
    {
      if(is.na(mapping$rna.index[n])){
        mapping$rna.index[n] = match(unlist(strsplit(as.character(mapping$gene[n]), ';'))[1], rownames(cpm))
      }
    }
    
    gg.expressed = unique(unlist(lapply(mapping$gene, function(x) unlist(strsplit(as.character(x), ';')))))
    
    cpm.bgs = cpm[which(is.na(match(rownames(cpm), gg.expressed))), ]
    
    cpm.xx = cpm[mapping$rna.index, ]
    cpm.xx[which(cpm.xx<10^-6)] = 2^-6
    cpm.xx[is.na(cpm.xx)] = 2^-6
    cpm.xx = log2(cpm.xx)
    
    #cpm.xx[which(is.na(cpm.xx))] = -6
    breaks = c(-6,-4, -2, 0, 2, 4, 6, 8, 20);
    cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks)-1)
    gap.rna = c(2)
    
    pheatmap(cpm.xx, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE, 
             color = cols, breaks = breaks,
             cluster_cols=FALSE, na_col = 'gray', annotation_col = df.rna, gaps_col = gap.rna)
    
    cpm.bgs[which(cpm.bgs < 10^-6)] = 2^-6
    cpm.bgs = log2(cpm.bgs)
    pheatmap(cpm.bgs, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE, 
             color = cols, breaks = breaks,
             cluster_cols=FALSE, na_col = 'gray', annotation_col = df.rna, gaps_col = gap.rna)
    
  }
  
  ##########################################
  # Quality controls again with negative controls for atac-seq peaks and RNA-seq data 
  ##########################################
  check.negative.controls.other.tissue.specific.genes = FALSE
  if(check.negative.controls.other.tissue.specific.genes){
    
    load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
    
    #xx = as.data.frame(promoters)
    #promoters.liver = promoters[which(!is.na(match(xx$geneID, controls.tissue$geneIDs[which(controls.tissue$tissues == 'liver')])))]
    #promoters.islet = promoters[which(!is.na(match(xx$geneID, controls.tissue$geneIDs[which(controls.tissue$tissues == 'islet')])))]
    #promoters.testis = promoters[which(!is.na(match(xx$geneID, controls.tissue$geneIDs[which(controls.tissue$tissues == 'testis')])))]
    
    #save(promoters.liver, promoters.islet, promoters.testis)
    #sel.promoter = which(overlapsAny(pp, promoters, ignore.strand = TRUE) == TRUE)
    #sel.promoter.dev = which(overlapsAny(pp, promoters.dev, ignore.strand = TRUE) == TRUE)
    
    # ii.promoter = grep('Promoter', peakAnnots$annotation)
    # sel.intron = grep('Intron ', peakAnnots$annotation)
    # sel.intergen = grep('Intergenic', peakAnnots$annotation)
    
    # cat(length(sel.promoter), ' peaks in promoters\n')
    # cat(length(sel.promoter.dev), ' peaks in promoters.dev \n')
    # cat(length(sel.promoter.liver), ' peaks in promoters.liver \n')
    # cat(length(sel.promoter.islet), ' peaks in promoters.islet \n')
    # cat(length(sel.promoter.testis), ' peaks in promoters.testis \n')
    gene.peaks = c()
    for(n in 1:nrow(mapping))
    {
      gx = unlist(strsplit(as.character(mapping$gene[n]), ';'))
      for(g in gx) gene.peaks = rbind(gene.peaks, c(g, mapping$peak.index[n]))
    }
    
    gene.peaks = data.frame(gene.peaks, stringsAsFactors = FALSE)
    colnames(gene.peaks) = c('gene', 'peak.index')
    gene.peaks$peak.index = as.integer(gene.peaks$peak.index)
    gene.peaks$peak = rownames(mapping)[gene.peaks$peak.index]
    
    library(ggplot2)
    # boxplot for peak signals
    for(n in 1:ncol(keep))
    {
      # n = 1
      cat(n, ' -- sample : ', colnames(keep)[n], '\n')
      xx = cbind(keep[,n], rep('promoter.peak', nrow(keep)))
      
      for(tissue.sel in c('liver', 'islet', 'spleen', 'testis'))
      {
        ctls = unique(controls.tissue$geneIDs[which(controls.tissue$tissues == tissue.sel)])
        
        keep.ctls = as.matrix(keep[gene.peaks$peak.index[match(ctls, gene.peaks$gene)], ])
        #cpm.ctls = as.matrix(cpm[match(ctls, rownames(cpm)), ])
        xx = rbind(xx, cbind(keep.ctls[,n], rep(tissue.sel, nrow(keep.ctls))))
      }
      
      xx = data.frame(xx)
      colnames(xx) = c('peak.signal', 'tissue')
      xx$peak.signal = as.numeric(as.character(xx$peak.signal))
      xx$peak.signal[is.na(xx$peak.signal)] = -6
      
      p <- ggplot(xx, aes(x=tissue, y=peak.signal, fill = tissue)) +
        geom_violin(trim=FALSE)
      p + geom_boxplot(width=0.1, fill="white")+
        labs(title=paste0('promoter peak signals - ', colnames(keep)[n]), x="tissue", y = "peak signals") + 
        ggsave(paste0(resDir, '/promoter.peaks.signals_vs_otherTissues.specific.genes.top3000_', colnames(keep)[n],  '.pdf'),
               width = 8, height = 4)
      
    }
    
    # boxplots for gene expression
    for(n in 1:ncol(cpm.xx))
    {
      # n = 1
      cat(n, ' -- sample : ', colnames(cpm.xx)[n], '\n')
      xx = cbind(cpm.xx[,n], rep('promoter.peak', nrow(cpm.xx)))
      
      for(tissue.sel in c('liver', 'islet', 'spleen', 'testis'))
      {
        ctls = unique(controls.tissue$geneIDs[which(controls.tissue$tissues == tissue.sel)])
        
        #cpm.ctls = as.matrix(keep[gene.peaks$peak.index[match(ctls, gene.peaks$gene)], ])
        cpm.ctls = as.matrix(cpm[match(ctls, rownames(cpm)), ])
        cpm.ctls[which(cpm.ctls <10^-6)] = 2^-6
        cpm.ctls[is.na(cpm.ctls)] = 2^-6
        cpm.ctls = log2(cpm.ctls)
        
        xx = rbind(xx, cbind(cpm.ctls[,n], rep(tissue.sel, nrow(cpm.ctls))))
      }
      
      xx = data.frame(xx)
      colnames(xx) = c('gene.expresion', 'tissue')
      xx$gene.expresion = as.numeric(as.character(xx$gene.expresion))
      xx$gene.expresion[is.na(xx$gene.expresion)] = -10
      
      p <- ggplot(xx, aes(x=tissue, y=gene.expresion, fill = tissue)) +
        geom_violin(trim=FALSE)
      p + geom_boxplot(width=0.1, fill="white")+
        labs(title=paste0('gene expression - ', colnames(cpm.xx)[n]), x="tissue", y = "gene expresion") + 
        ggsave(paste0(resDir, '/promoter.peaks.signals_vs_otherTissues.specific.genes.top3000_', colnames(cpm.xx)[n],  '.pdf'),
               width = 8, height = 4)
      
    }
    
    
    Visulize.with.heatmap = FALSE
    if(Visulize.with.heatmap){
      tissue.sel = 'liver'
      ctls = unique(controls.tissue$geneIDs[which(controls.tissue$tissues == tissue.sel)])
      
      keep.ctls = as.matrix(keep[gene.peaks$peak.index[match(ctls, gene.peaks$gene)], ])
      cpm.ctls = as.matrix(cpm[match(ctls, rownames(cpm)), ])
      
      xx = keep.ctls[which(!is.na(keep.ctls[,1])), ]
      
      xx = xx[order(-apply(xx, 1, max)), ]
      
      cpm.ctls[which(cpm.ctls < 10^-6)] = NA
      cpm.ctls = log2(cpm.ctls)
      
      breaks = c(-6,-4, -2, 0, 2, 4, 6, 8, 20);
      cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks)-1)
      gap.rna = c(2)
      
      cpm.ctls = as.matrix(cpm.ctls)
      rownames(cpm.ctls) = c(1:nrow(cpm.ctls))
      #colnames(cpm.ctls) = NULL
      pheatmap(cpm.ctls, cluster_rows=FALSE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE, 
               color = cols, breaks = breaks,
               cluster_cols=FALSE, na_col = 'gray', annotation_col = df.rna, gaps_col = gap.rna) 
      
      
      #keep.xx = as.data.frame(t(apply(keep, 1, norm_minmax)))
      breaks = c(-6,-4,-2,0, 2, 4, 8, 12);
      cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks)-1)
      pheatmap(keep.ctls, cluster_rows=FALSE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE, 
               color = cols,
               cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, na_col = 'darkgray',
               breaks = breaks)
    }
    
    
    
  }
  
  ##########################################
  # check the gene expression difference in dev and mature samples
  # separate house-keeping gene apart
  ##########################################
  Compare.Expression.Dev.vs.Mature = FALSE
  if(Compare.Expression.Dev.vs.Mature){
    
    cpm.peak = cbind(apply(cpm.xx[, c(1:2)], 1, mean), apply(cpm.xx[, c(3:4)], 1, mean))
    colnames(cpm.peak) = c('E40', 'mUA')
    
    rr = cpm.peak[, 2] - cpm.peak[, 1]
    
    #names.uniq = unique(rownames(cpm.peak))
    #cpm.peak = cpm.peak[match(names.uniq, rownames(cpm.peak)), ]
    gg.names = data.frame(geneID = rownames(cpm.peak), stringsAsFactors = FALSE)
    gg.names$hs.symbol = annot$gene.symbol.hs_gtf.gene[match(gg.names$geneID, annot$geneID)]
    
    hkgs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
    gg.names$hkg = NA
    gg.names$hkg[!is.na(match(gg.names$geneID, hkgs))] = 1
    gg.names$ratio.mUA.vs.E40 = rr
    gg.names$hkg[!is.na(match(gg.names$geneID, hkgs)) & abs(rr) < 1.5] = 1
    
    # quickly check promoter peaks after removing the house-keeping genes
    breaks = c(-6, -4, -2, 0, 2, 4, 8, 12);
    cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks)-1)
    
    pheatmap(keep[which(is.na(gg.names$hkg)), ], cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE, 
             color = cols,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps, 
             breaks = breaks)
    
    
    mature.genes = c( 'Twist2', 'Lect1',  'Tnmd', 'Sulf2', 'Col4a2', 'Otos')
    dev.genes = c('Ccnb1', 'Nrep', 'Prdx2', 'Hand1', 'Hand2', 'Meis1', 'Meis2', 'Hoxa11', 'Hoxa13', 
                  'Shh', 'Fgf8', 'Fgf10', 'Grem1', 'Gli1', 'Fgf9')
    mature.genes = toupper(mature.genes)
    dev.genes = toupper(dev.genes)
    
    Test.ggplot = FALSE
    if(Test.ggplot){
      library(ggplot2)
      cpm.peak$gene.symbols = annot$gene.symbol.hs_gtf.gene[match(rownames(cpm.peak), annot$geneID_gtf.gene)]
      cpm.peak$hkg = NA 
      cpm.peak$hkg[which(gg.names$hkg == 1)] = 'hkg'
      cpm.peak$dev.genes = NA; cpm.peak$mature.genes = NA
      cpm.peak$dev.genes[!is.na(match(cpm.peak$gene.symbols, dev.genes))] = 'dev.gene'
      cpm.peak$mature.genes[!is.na(match(cpm.peak$gene.symbols, mature.genes))] = 'mature.gene'
      
      ggplot(cpm.peak, aes(x=E40, y=mUA)) +
        geom_point() + # Show dots
        geom_text(
          label=, 
          nudge_x = 0.25, nudge_y = 0.25, 
          check_overlap = TRUE
        )
      
    }
    
    
    pdfname = paste0(resDir, "/Expression_Dev.vs.Mature_promoterPeakGenes.pdf")
    pdf(pdfname, width = 12, height = 10)
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
    
    
    # Keep 30 first rows in the mtcars natively available dataset
    #data=head(mtcars, 30)
    cpm.peak = data.frame(cpm.peak, stringsAsFactors = FALSE)
    
    plot(cpm.peak, cex = 0.2, main = 'gene expression - log2 (normalized counts)', type = 'n')
    points(cpm.peak[which(gg.names$hkg == 1), c(1:2)], cex = 0.3, col = 'darkgray', pch = 16)
    points(cpm.peak[which(is.na(gg.names$hkg)), c(1:2)], cex = 0.3, col = 'darkorange', pch =1)
    mm = match(dev.genes, gg.names$hs.symbol)
    points(cpm.peak[mm, c(1:2)], cex = 1.5, col = 'darkgreen', pch = 16)
    text(cpm.peak[mm, c(1:2)], gg.names$hs.symbol[mm], 
         cex = 1.0, pos = 4, offset = 0.8, col = 'darkblue')
    
    mm = match(mature.genes, gg.names$hs.symbol)
    points(cpm.peak[mm, c(1:2)], cex = 1.5, col = 'magenta', pch = 16)
    text(cpm.peak[mm, c(1:2)], gg.names$hs.symbol[mm], 
         cex = 1.0, pos = 2, offset = 0.8, col = 'darkblue')
    
    abline(0, 1, lwd = 2.0, col = 'red')
    abline(-1.5, 1, lwd = 2.0, col = 'red')
    abline(1.5, 1, lwd = 2.0, col = 'red')
    abline(v = 0, col = 'gray', lwd = 2.0)
    abline(h = 0, col = 'gray', lwd = 2.0)
    
    dev.off()
      
  }
  
  ##########################################
  # GO term enrichment analysis 
  ##########################################
  GO.term.analysis = FALSE
  if(GO.term.analysis){
    
    library(enrichplot)
    library(clusterProfiler)
    library(openxlsx)
    library(ggplot2)
    library(stringr)
    #library(org.Ce.eg.db)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
    annot = readRDS(file = paste0(annotDir, 'geneAnnotation_geneSymbols.synteny.evidence.hs.nr.same.evidence.rds'))
    annot = data.frame(annot, stringsAsFactors = FALSE)
    
    gg.expressed = unique(unlist(lapply(mapping$gene, function(x) unlist(strsplit(as.character(x), ';')))))
    
    Use.annot.synteny = TRUE
    if(Use.annot.synteny){
      #annot = readRDS(file = paste0(annotDir, 
      #        'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
      
      gg.expressed = unique(annot$gene.symbol.hs_gtf.gene[match(gg.expressed, annot$geneID)])
      gg.expressed = gg.expressed[which(gg.expressed != '' & gg.expressed != 'N/A')]
      
      xx0 = unique(annot$gene.symbol.hs_gtf.gene)
      xx0 = xx0[which(xx0 != '' & xx0 != 'N/A')]
      bgs = xx0
    }else{
      load(file =  paste0(annotDir, 
                          'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
      
      kk = which(!is.na(annot$gene.evidence.same.gene.symbol.hs.nr) & annot$gene.evidence.same.gene.symbol.hs.nr != '' &
                   annot$gene.evidence.same.gene.symbol.hs.nr != 'N/A')
      annot = annot[kk, ]
      
      bgs0 = unique(annot$gene.evidence.same.gene.symbol.hs.nr)
      
      # remove house-keeping gene 
      hkgs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
      mm = match(hkgs, annot$geneID)
      mm = mm[which(!is.na(mm))]
      annot = annot[-mm, ]
      
      gg.expressed = unique(annot$gene.evidence.same.gene.symbol.hs.nr[match(gg.expressed, annot$geneID)])
      bgs = unique(annot$gene.evidence.same.gene.symbol.hs.nr)
      
    }
    
    gg.expressed = firstup(tolower(gg.expressed))
    bgs = firstup(tolower(bgs))
    bgs0 = firstup(tolower(bgs0))
    
    gene.df <- bitr(gg.expressed, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)
    head(gene.df)
    
    bgs.df <- bitr(bgs, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)
    
    head(bgs.df)
    
    bgs0.df <- bitr(bgs0, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Mm.eg.db)
    
    #pval.cutoff = 0.05
    
    #pdfname = paste0(resDir, "/GO.Terms_enrichment_mir1.mutant_upregulated_downregulated_genes.pdf")
    #pdf(pdfname, width = 12, height = 10)
    #par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
    
    #kk.up = which(res$pvalue_mir1.mutant.vs.wt < pval.cutoff & res$log2FoldChange_mir1.mutant.vs.wt > 0)
    #kk.down = which(res$pvalue_mir1.mutant.vs.wt < pval.cutoff & res$log2FoldChange_mir1.mutant.vs.wt < 0)
    #cat('nb of upregulated genes : ', length(kk.up), '\n')
    #cat('nb of downregulated genes : ', length(kk.down), '\n')
    
    ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                        #universe     = bgs0.df$ENSEMBL,
                        #OrgDb         = org.Hs.eg.db,
                        OrgDb         = org.Mm.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
    
    #head(ego)
    barplot(ego, showCategory=20) + ggtitle("Go term enrichment for promoter peaks in Akane's ATAC-seq data")
    
    write.csv(ego, file = paste0(resDir, "GO_term_enrichmenet_for_upregulated_genes_pval_0.05.csv"), 
              row.names = TRUE)
    
    # kegg.up <- enrichKEGG(gene         = uniprotnames[kk.up],
    #                  organism     = 'cel',
    #                  universe     = uniprotnames,
    #                  keyType = 'uniprot',
    #                  pvalueCutoff = 0.05)
    # head(kegg.up)
    # barplot(kegg.up, showCategory=20) + ggtitle("kegg for upregulated genes")
    # 
    # 
    
    #dev.off()
    
  }
  
  ##########################################
  # After confirming that there is no contamination and the promoters will be grouped to
  # house-keeping genes, developmental genes, mature genes and both from the gene expression data
  ##########################################
  Grouping.promoter.peaks = FALSE
  if(Grouping.promoter.peaks){
    # promoters without house-keeping genes
    hkgs = readRDS(file = paste0(annotDir, 'axolotl_house.keep.genes_expressed.top10k.across18tissues.rds'))
    xx = as.data.frame(promoters)
    mm = match(xx$geneID, hkgs$amexID)
    ii = which(is.na(mm))
    
    promoters.dev = promoters[ii]
    promoters.dev.reduced = reduce(promoters.dev, drop.empty.ranges = TRUE, ignore.strand = TRUE)
    
    saveRDS(promoters.dev, 
            file = paste0(RdataDir, 
            '/axolotl_promoters_Granges_AmexT_v47_putative.full.length_N.terminal_upstream2kb.downstream2kb_rmHouseKeepingGenes.rds'))
    
    saveRDS(promoters.dev.reduced, 
            file = paste0(RdataDir, 
    '/axolotl_promoters_Granges_AmexT_v47_putative.full.length_N.terminal_upstream2kb.downstream2kb_reduced_reHouseKeepingGenes.rds'))
    
  }
  
  
}

########################################################
########################################################
# Section : atac-seq normalization and batch correction
# 
########################################################
########################################################
normalize.batch.correct = function(fpm, design, norm.batch.method ='TMM.combat')
{
  ##########################################
  # normalized by DESeq2 
  ##########################################
  if(norm.batch.method =='DESEq2.combat'){
    fpm = fpm(dds, robust = FALSE)
    fpm = log2(fpm + 2^-4)
    
    make.pca.plots(fpm, ntop = 1000, conds.plot = 'all')
    
    # combat batch correction
    require("sva")
    bc = as.factor(design$batch)
    mod = model.matrix(~ as.factor(conds), data = design)
    fpm.bc = ComBat(dat=fpm, batch=bc, mod=mod, par.prior=TRUE, ref.batch = '2021')    
    
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'all')
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'Dev.Mature')
    
    
    # limma batch correction
    design.tokeep<-model.matrix(~ 0 + conds,  data = design)
    cpm.bc = removeBatchEffect(fpm, batch = bc, design = design.tokeep)
    
    make.pca.plots(cpm.bc, ntop = 5000, conds.plot = 'all') 
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'Dev.Mature')
    
  }
  
  ##########################################
  # test TMM in edgeR 
  ##########################################
  if(norm.batch.method =='TMM.combat'){
    library(edgeR)
    
    d <- DGEList(counts=counts(dds), group=design$conds)
    tmm <- calcNormFactors(d, method='TMM')
    fpm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
    #fpm = log2(tmm + 1)
    
    require("sva")
    bc = as.factor(design$batch)
    mod = model.matrix(~ as.factor(conds), data = design)
    fpm.bc = ComBat(dat=fpm, batch=bc, mod=mod, par.prior=TRUE, ref.batch = '2021')    
    fpm = fpm.bc
    
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'all')
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'Dev.Mature')
    
  }
  
  ##########################################
  # test cpn, original code from 
  # https://github.com/koenvandenberge/bulkATACGC/blob/master/scripts/benchmark/20191209_calderon2019.Rmd 
  ##########################################
  if(norm.batch.method =='CPN.combat'){
    library(edgeR)
    library(cqn)
    library(GenomicAlignments)
    library(scone)
    #library(qsmooth)
    #library(rhdf5)
    
    # get GC content
    pp = data.frame(t(sapply(rownames(dds), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    rownames(pp) = gsub('_', '-', rownames(pp))
    pp$strand = '*'
    gr = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    rerun.get.GC.content = FALSE
    if(rerun.get.GC.content){
      #rn <- rownames(dds)
      #sn <- unlist(lapply(lapply(strsplit(rn,split=":"),"[",1),function(x) gsub(pattern="chr",x=x,replacement="chr")))
      #start <- as.numeric(unlist(lapply(strsplit(rn,split="_"),"[",2)))
      #end <- as.numeric(unlist(lapply(strsplit(rn,split="_"),"[",3)))
      #gr <- GRanges(seqnames=sn, ranges=IRanges(start, end), strand="*", mcols=data.frame(peakID=rn))
      ff <- FaFile("/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/AmexG_v6.DD.corrected.round2.chr.fa")
      peakSeqs <- getSeq(x=ff, gr)
      
      gcContentPeaks <- letterFrequency(peakSeqs, "GC", as.prob=TRUE)[,1]
      gcGroups <- Hmisc::cut2(gcContentPeaks, g=20)
      
      saveRDS(peakSeqs, file = paste0(RdataDir, '/peakSeqs.for.GCcontetent.estimation.rds' ))
      saveRDS(gcContentPeaks, file = paste0(RdataDir, '/peak.GCcontetent.estimation.rds' ))
    }
    
    gcContentPeaks = readRDS(file = paste0(RdataDir, '/peak.GCcontetent.estimation.rds' ))
    ss = colSums(counts(dds))
    ss = median(ss) * sizeFactors(dds)
    library(cqn)
    cqnModel <- cqn(counts = counts(dds), x=gcContentPeaks, lengths=1000, 
                    lengthMethod="fixed", 
                    sizeFactors = ss)
    
    fpm = (cqnModel$y + cqnModel$offset)
    #fpm = log2(countsCqn)
    #make.pca.plots(fpm, ntop = 10000, conds.plot = 'all') 
    
    # combat batch correction
    require("sva")
    bc = as.factor(design$batch)
    mod = model.matrix(~ as.factor(conds), data = design)
    fpm.bc = ComBat(dat=fpm, batch=bc, mod=mod, par.prior=TRUE, ref.batch = '2021')    
    
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'all')
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'Dev.Mature')
    
    # limma batch correction
    design.tokeep<-model.matrix(~ 0 + conds,  data = design)
    cpm.bc = removeBatchEffect(fpm, batch = bc, design = design.tokeep)
    
    make.pca.plots(cpm.bc, ntop = 10000, conds.plot = 'all') 
    
  }
  
  ##########################################
  # quantile normalization
  ##########################################
  if(norm.batch.method =='Quantile.combat'){
    fpm = fpm(dds, robust = TRUE)
    fpm = log2(fpm + 2^-4)
    library(preprocessCore)
    fpm.qn = normalize.quantiles(fpm)
    colnames(fpm.qn) = colnames(fpm)
    rownames(fpm.qn) = rownames(fpm)
    fpm = fpm.qn
    rm(fpm.qn)
    
    require("sva")
    bc = as.factor(design$batch)
    mod = model.matrix(~ as.factor(conds), data = design)
    fpm.bc = ComBat(dat=fpm, batch=bc, mod=mod, par.prior=TRUE, ref.batch = '2021')    
    
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'all')
    make.pca.plots(fpm.bc, ntop = 5000, conds.plot = 'Dev.Mature')
    
  }
  
  return(fpm)
  
}


make.pca.plots = function(fpm, ntop = 1000, conds.plot = 'Dev.Mature')
{
  library(factoextra)
  
  if(conds.plot == 'all'){
    conds.sel = c('Embryo_Stage40_93',  'Embryo_Stage40_13', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
                  'Mature_UA_13', 'Mature_UA_74938|Mature_UA_102655', 'Mature_LA', 'Mature_Hand', 
                  'BL_UA_5days_13', 'BL_UA_5days_89',  'BL_UA_9days_13', 'BL_UA_13days_proximal',  'BL_UA_13days_distal')
  }
  
  if(conds.plot == 'Dev.Mature')
  {
    conds.sel = c('Embryo_Stage40_93',  'Embryo_Stage40_13', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
                  'Mature_UA_13', 'Mature_UA_74938|Mature_UA_102655', 'Mature_LA', 'Mature_Hand')
  }
  
  sample.sel = c()
  ii.gaps = c()
  for(n in 1:length(conds.sel)) {
    c = conds.sel[n]
    sample.sel = c(sample.sel, grep(c, colnames(fpm)))
    if(n == 1) {
      ii.gaps = c(ii.gaps, length(grep(c, colnames(fpm))))
    }else{
      if(n != length(conds.sel)) ii.gaps = c(ii.gaps, (ii.gaps[n-1] + length(grep(c, colnames(fpm)))))
    }
  }
  
  #sels = sample(c(1:nrow(fpm)), size = 10000, replace = FALSE)
  
  xx = as.matrix(fpm[, sample.sel])
  vars = apply(xx, 1, var)
  xx = xx[order(-vars), ]
  xx = xx[1:ntop, ]
  #par(mfrow = c(2,2))
  #pairs(xx[, c(1:4)])
  #plot.pair.comparison.plot(xx[, c(1:4)], linear.scale = FALSE)
  #plot.pair.comparison.plot(xx[, c(9:16)], linear.scale = FALSE)
  
  res.pca <- prcomp(t(xx), scale = TRUE)
  #res.var <- get_pca_var(res.pca)
  
  fviz_pca_ind(res.pca,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  
}

##########################################
# make heatmap for atac-seq peaks
##########################################
make.heatmap.atacseq = function(cpm)
{
  # cpm = fpm
  # make Granges for non-background peaks 
  jj = grep('bg_', rownames(cpm), invert = TRUE)
  cpm = cpm[jj, ]
  rownames(cpm) = gsub('_', '-', rownames(cpm))
  
  pp = data.frame(t(sapply(rownames(cpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  #rownames(pp) = gsub('_', '-', rownames(pp))
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  # annotate
  require(ChIPpeakAnno)
  require(ChIPseeker)
  
  amex = GenomicFeatures::makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))
  #amex = makeTxDbFromGFF(file = paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf'))
  #amex = readRDS(file = paste0(annotDir, 'TxDb_ax6_UCSC_2021_01_26_genes.putative.full.length.rds')) # reimport object does not work
  pp.annots = as.data.frame(annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript'))
  
  
  # select peaks within regions
  HoxA = data.frame(chr = 'chr2p', start = 873085043, end = 884416919, strand = '*', stringsAsFactors = FALSE)
  HoxA = makeGRangesFromDataFrame(HoxA, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  HoxD1 = data.frame(chr = 'chr9q', start = 416423355, end = 427456848, strand = '*', stringsAsFactors = FALSE)
  HoxD1 = makeGRangesFromDataFrame(HoxD1, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  HoxD2 = data.frame(chr = 'chr9q', start = 426557960, end = 435711507, strand = '*', stringsAsFactors = FALSE)
  HoxD2 = makeGRangesFromDataFrame(HoxD2, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  ii.HoxA = which(overlapsAny(pp, HoxA, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))
  ii.HoxD1 = which(overlapsAny(pp, HoxD1, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))
  ii.HoxD2 = which(overlapsAny(pp, HoxD2, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))
  
  ii.HoxA = which(overlapsAny(pp, HoxA, ignore.strand = TRUE) == TRUE) 
  ii.HoxD1 = which(overlapsAny(pp, HoxD1, ignore.strand = TRUE) == TRUE)
  ii.HoxD2 = which(overlapsAny(pp, HoxD2, ignore.strand = TRUE) == TRUE)
  
  peak.sels = c(ii.HoxA, ii.HoxD1, ii.HoxD2)
  peak.sels = unique(peak.sels)
  
  peaks.promoters = which(overlapsAny(pp, promoters))
  
  ii.promoters = intersect(peak.sels, peaks.promoters)
  ii.nonpromoters = setdiff(peak.sels, ii.promoters)
  
  peak.sels = c(ii.promoters, ii.nonpromoters)
  
  # select conditions
  #conds.sel = c('Embryo_Stage40_93',  'Embryo_Stage40_13', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
  #              'Mature_UA_13',  'Mature_LA', 'Mature_Hand', 'BL_UA_5days_13', 'BL_UA_5days_89',  'BL_UA_9days_13', 
  #            'BL_UA_13days_proximal',  'BL_UA_13days_distal')
  conds.sel = c('Mature_UA_13', 'BL_UA_5days_13', 'BL_UA_9days_13', 'BL_UA_13days_proximal',  'BL_UA_13days_distal')
  #conds.sel = c('Embryo_Stage40_93',  'Embryo_Stage40_13', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
  #               'Mature_UA_13', 'Mature_UA_74938|Mature_UA_102655', 'Mature_LA', 'Mature_Hand')
  
  sample.sel = c()
  ii.gaps = c()
  for(n in 1:length(conds.sel)) {
    c = conds.sel[n]
    sample.sel = c(sample.sel, grep(c, colnames(cpm)))
    if(n == 1) {
      ii.gaps = c(ii.gaps, length(grep(c, colnames(cpm))))
    }else{
      if(n != length(conds.sel)) ii.gaps = c(ii.gaps, (ii.gaps[n-1] + length(grep(c, colnames(cpm)))))
    }
  }
  
  df <- data.frame(colData(dds)[,c("conds", 'batch')])[sample.sel, ]
  
   
  keep =cpm[peak.sels, sample.sel]
  cc = unique(df$conds)
  keep.m = matrix(NA, nrow = nrow(keep), ncol = length(cc))
  for(n in 1:length(cc))
  {
    keep.m[,n] = apply(keep[, which(df$conds == cc[n])], 1, mean)
    
  }
  
  colnames(keep.m) = cc
  rownames(keep.m) = rownames(keep)
  keep = keep.m
  df = df[match(cc, df$conds),]
  
  filtering.with.signal = FALSE
  if(filtering.with.signal){
    
    ss = apply((keep), 1, max)
    
    max.cutoff = 3
    
    hist(ss, breaks = 100)
    abline(v = max.cutoff, col='red', lwd =2.0)
    cat(length(which(ss>max.cutoff)), 'peak selected\n')
    
    keep = keep[which(ss > max.cutoff), ]
    
  }
  
  keep = as.matrix(keep)
  ddf = data.frame(conds = df[,1])
  rownames(ddf) = rownames(df)
  
  #breaks = c(0, 1, 2, 3, 4, 6);
  #cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaks)-1)
  pheatmap(keep, cluster_rows=FALSE, show_rownames=TRUE, fontsize_row = 7,
           scale = 'none', show_colnames = FALSE,
           cluster_cols=FALSE, annotation_col = ddf, 
           col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(5), 
           filename = paste0(resDir, '/HoxA_HoxD_clusters_selectedPeaks_promoters.nonpromoters_withCoordinates.pdf'), 
           width = 10, height = 12)
  
  
  #pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
  #         cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
  
}


########################################################
########################################################
# Section : check the stage40 sample quality by comparing 
# bulk-data with pseudo-bulk from scRNA-seq in Gerber et al., 2018
########################################################
########################################################
pseudo.bulk.by.pooling.scRNAseq = function()
{
  scRNAseq.counts = read.delim(file = 
                  '../../limbRegeneration_scRNA/raw_NGS/axolotl/Gerber_2018/Fluidigm_C1/nf_out/featureCounts/merged_gene_counts.txt',
                               header = TRUE)
  
  metadata = read.delim(file = '../../limbRegeneration_scRNA/raw_NGS/axolotl/Gerber_2018/Metadata.txt')
  metadata = data.frame(metadata[which(metadata$Technology == 'Fluidigm C1'), ])
  metadata$Sample = as.character(metadata$Sample)
  
  condition = metadata$Sample
  condition = gsub('forelimb upper arm ', '', condition)
  condition = gsub(' forelimb limb bud', '', condition)
  condition = gsub(' ', '_', condition)
  #condition = gsub('11_dpa', 'dpa11', condition)
  
  metadata$condition = condition
  
  raw = as.matrix(scRNAseq.counts[ ,-1])
  rownames(raw) = scRNAseq.counts$ENSEMBL_ID
  
  names = colnames(raw)
  names = gsub('_1.sorted_gene.featureCounts.txt', '', names)
  colnames(raw) = names
   
  conds = unique(metadata$condition)
  pseudo = matrix(NA, ncol = length(conds), nrow = nrow(raw))
  colnames(pseudo) = conds
  rownames(pseudo) = rownames(raw)
  
  for(n in 1:length(conds))
  {
    cat(n, ' : ', conds[n], ' -- ')
    jj = which(metadata$condition == conds[n])
    mm = match(metadata$Run[jj], colnames(raw))
    if(length(which(is.na(mm)))>0) {
      cat('some cells lost \n')
    }else{
      cat(length(mm), ' cell found \n')
      xx = as.matrix(raw[, mm])
      xx[is.na(xx)] = 0
      pseudo[,n] = apply(xx, 1, sum)
    }
  }
  
  saveRDS(pseudo, file = paste0(RdataDir, 'pseudoBulk_scRNAcellPooling_FluidigmC1.rds'))
  
}


compare.Akane.RNAseq.pseudobulk.scRNAseq = function()
{
  require(ggplot2)
  require(DESeq2)
  library("dplyr")
  library("ggplot2")
  
  pseudo = readRDS(file = paste0(RdataDir, 'pseudoBulk_scRNAcellPooling_FluidigmC1.rds'))
  
  xx = pseudo[, match(c('uninjured', 'Stage_40', 'Stage_44'), colnames(pseudo))]
  
  load(file=paste0(RdataDir, 'Design_stats_readCounts_', version.analysis, '.Rdata'))
  
  raw = as.matrix(all[, -1])
  rownames(raw) = all$gene
  jj = c(which(design$conditions == 'Mature_UA' & design$batch != '1'), 
         which(design$conditions == 'Embryo_Stage40'), 
         which(design$conditions == 'Embryo_Stage46_proximal'), 
         which(design$conditions == 'Embryo_Stage46_distal'))
  
  design = design[jj, ]
  raw = raw[,jj]
  
  design = design[, c(1, 2, 5, 14)]
  
  merge.stage46 = TRUE
  if(merge.stage46){
    kk = grep('Embryo_Stage46_', design$conditions)
    design = design[-kk, ]
    design = rbind(design, 
                   c('merged', 'Embryo_Stage46', '1'))
    merged.data = apply(raw[, kk], 1, sum)
    raw = raw[, -kk]
    raw = cbind(raw, merged.data)
    colnames(raw)[ncol(raw)] = 'merged.stage46.batch1'
    
    kk = which(design$conditions == 'Embryo_Stage40' & design$batch == '1')
    design = design[-kk, ]
    design = rbind(design, 
                   c('merged', 'Embryo_Stage40', '1'))
    merged.data = apply(raw[, kk], 1, sum)
    raw = raw[, -kk]
    raw = cbind(raw, merged.data)
    colnames(raw)[ncol(raw)] = 'merged.stage40.batch1'
    
     
  }
  
  ## merge samples
  yy = cbind(raw, xx[match(rownames(raw), rownames(xx)), ])
  design = rbind(design, 
                 c('xx0', 'Mature_UA', '0'), 
                 c('xx1', 'Embryo_Stage40', '0'), 
                 c('xx2', 'Embryo_Stage46', '0'))
  
  design = data.frame(design, stringsAsFactors = FALSE)
  rownames(design) = colnames(yy)
  
  
  
  # DESeq2 object and gene selection (probably important due to the variable sequencing depth)
  dds <- DESeqDataSetFromMatrix(yy, DataFrame(design), design = ~ conditions)
  counts = counts(dds)
  
  ss = rowSums(counts)
  hist(log10(ss), breaks = 200, main = 'log2(sum of reads for each gene)')
  
  quantile(ss, probs = c(0.75, 0.98))
  #cutoff.peak = 2^10
  #cat(length(which(ss > cutoff.peak)), 'genes selected \n')
  # dds <- dds[ss > cutoff.peak, ]
  kk = which(ss >= quantile(ss, probs = 0.75) & ss < quantile(ss, probs = 0.98) )
  dds = dds[kk, ]
  
  counts = counts(dds)
  
  jj = c(1:nrow(counts))
  for(bs in unique(design$batch)){
    ss = apply(counts[, which(design$batch == bs)], 1, mean)
    jj1 = which(ss > 50)
    jj = intersect(jj, jj1)
  }
  length(jj)
  #nb.expressed = apply(counts, 1, function(x) length(which(x > 50)))
  #hist(nb.expressed)
  #length(which(nb.expressed>=12))
  #dds = dds[which(nb.expressed >= 12), ]
  dds = dds[jj, ]
  
  dds <- estimateSizeFactors(dds)
  ss = colSums(counts(dds))
  
  plot(sizeFactors(dds), ss/10^6, log = 'xy', xlab = 'size factor DESeq2', ylab = 'library size (million reads)')
  text(sizeFactors(dds), ss/10^6, names(ss), cex = 0.7)
  
  nb.genes = apply(counts(dds), 2, function(x) length(which(x>10)))
  plot(ss/10^6, nb.genes)
  text(ss/10^6, nb.genes, names(ss), cex = 0.7)
  
  fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  ntop = 1000
  pca=plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = FALSE, ntop = ntop)
  print(pca)
  
  pca2save = as.data.frame(plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = TRUE, ntop = ntop))
  
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 0.7, nudge_y = 1, size=3)
  
  plot(ggp) 
  
  
  ##########################################
  # assume that all difference are from batch 
  # after batch correction, sample conditions should be close to each other
  ##########################################
  require('limma')
  logcpm = log2(fpm + 2^-6)
  
  #sels = which(design$batch != '1')
  sels = c(1:nrow(design))
  
  bc = as.factor(design$batch[sels])
  
  design.tokeep<-model.matrix(~ 0 + conditions,  data = design[sels, ])
  cpm.bc = removeBatchEffect(logcpm[, sels], batch = bc, design = design.tokeep)
  vv = apply(cpm.bc, 1, var)
  
  vv.cutoff = 3
  cat('nb of top genes selected : ', length(which(vv>vv.cutoff)), '\n')
  
  cpm.bc = cpm.bc[which(vv>vv.cutoff), ]
  pca = prcomp(t(cpm.bc), scale. = TRUE, center = TRUE)
  pca2save = data.frame(pca$x, conditions=design$conditions[sels], 
                        batch = bc, 
                        name=colnames(cpm.bc))
  
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
    geom_point(size=4) + 
    geom_text(hjust = 0.7, nudge_y = 1, size=4)
  plot(ggp) 
  
  require("sva")
  xx = as.matrix(logcpm[, sels])
  vv = apply(xx, 1, var)
  xx = xx[which(vv > 0.), ]
  
  mod = model.matrix(~ as.factor(conditions), data = design[sels, ])
  logcpm.bc = ComBat(dat=xx, batch=bc, mod=mod, par.prior=TRUE, ref.batch = '0')
  vv = apply(logcpm.bc, 1, var)
  
  vv.cutoff = 1
  cat('nb of top genes selected : ', length(which(vv>vv.cutoff)), '\n')
  
  logcpm.bc = logcpm.bc[which(vv>vv.cutoff), ]
  pca = prcomp(t(logcpm.bc), scale. = TRUE, center = TRUE)
  pca2save = data.frame(pca$x, conditions=design$conditions[sels], 
                        batch = bc, 
                        name=colnames(logcpm.bc))
  
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
    geom_point(size=4) + 
    geom_text(hjust = 0.7, nudge_y = 1, size=4)
  plot(ggp) 
  
  
  ##########################################
  # compare the logFC in different batch without any batch correction
  ##########################################
  logcpm = data.frame(log2(fpm))
  
  # batch 0
  ref1 = logcpm$Stage_40 - logcpm$uninjured
  ref2 = logcpm$Stage_44 - logcpm$uninjured
  
  # batch 3
  means = apply(cbind(logcpm$Embryo_Stage40_137331.batch3, logcpm$Embryo_Stage40_137332.batch3), 1, mean)
  xx1 = means - logcpm$Mature_UA_136143.batch3
  xx2 = means - logcpm$Mature_UA_136149.batch3
  # batch 2
  xx3 = apply(logcpm[, which(design$batch == '2' & design$conditions == 'Embryo_Stage40')], 1, mean) -
    apply(logcpm[, which(design$batch == '2' & design$conditions == 'Mature_UA')], 1, mean)
    
  xx4 = logcpm$merged.stage46.batch1 -
    apply(logcpm[, which(design$batch == '2' & design$conditions == 'Mature_UA')], 1, mean)  
  
    
  cal.correlation = function(x1, x2)
  {
    x1 = as.numeric(x1)
    x2 = as.numeric(x2)
    
    jj = which(!is.na(x1) & !is.na(x2) & abs(x1) != Inf & abs(x2) != Inf)
    cat('pearson.cor : ', cor(x1[jj], x2[jj], use = "na.or.complete", method = 'pearson'), ' -- spearma.cor : ', 
        cor(x1[jj], x2[jj], method = 'spearman'), '\n')
    
  }
  
  par(mfrow = c(2, 3))
  cex = 0.3
  cols = 'red'
  #jj = which(abs(ref1) < 3 & abs(xx1)<3)
  plot(ref1, xx1, cex = cex, xlab = 'scRNA pooling', ylab = 'batch 3 (mUA.136143)', main = 'stage40 vs. mUA (logFC)')
  abline(0, 1, lwd = 2.0, col = cols)
  
  #jj = which(abs(ref1) < 3 & abs(xx2)<3)
  plot(ref1, xx2, cex = cex, xlab = 'scRNA pooling', ylab = 'batch 3 (mUA.136149)', main = 'stage40 vs. mUA (logFC)')
  abline(0, 1, lwd = 2.0, col = cols)
  
  #jj =  which(abs(ref1) < 3 & abs(xx3)<3)
  plot(ref1, xx3, cex = cex, xlab = 'scRNA pooling', ylab = 'batch 2 ', main = 'stage40 vs. mUA (logFC)')
  abline(0, 1, lwd = 2.0, col = cols)
  
  
  plot(xx1, xx2, cex = cex, xlab = 'batch 3 (mUA.136143)', ylab = 'batch 3 (mUA.136149)', main = 'stage40 vs. mUA (logFC)')
  abline(0, 1, lwd = 2.0, col = cols)
  
  plot(xx1, xx3, cex = cex, xlab = 'batch 3 (mUA.136143)', ylab = 'batch 2', main = 'stage40 vs. mUA (logFC)')
  abline(0, 1, lwd = 2.0, col = cols)
  
  plot(ref2, xx4, cex = cex, xlab = 'scRNA pooling', ylab = 'batch 1', main = 'stage46 vs. mature UA')
  abline(0, 1, lwd = 2.0, col = cols)
  
  cal.correlation(ref1, xx1)
  cal.correlation(ref1, xx2)
  cal.correlation(ref1, xx3)
  cal.correlation(xx1, xx2)
  cal.correlation(xx1, xx3)
  cal.correlation(xx2, xx3)
  cal.correlation(ref2, xx4)
  
}


########################################################
########################################################
# Section : grouping atac-seq peaks or RNA-seq data
# 
########################################################
########################################################
plot.peak.profiles = function(peak.name, fpm = NULL, mains = NULL)
{
  # peak.name = names
  if(is.null(fpm)){
    
    load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.withinPeaks.pval6.Rdata'))
    sels = grep('Mature|Embryo|BL_UA', design$conds)
    
    dds <- DESeqDataSetFromMatrix(as.matrix(counts[, sels]), DataFrame(design[sels, ]), design = ~ conds)
    
    select.peaks.with.readThreshold = TRUE
    quantile.normalization = TRUE
    
    if(select.peaks.with.readThreshold){
      #ss = rowMax(counts(dds)[, grep('Embryo_', dds$conds)])
      ss = rowMax(counts(dds))
      
      #hist(log10(ss), breaks = 200, main = 'log2(sum of read within peaks) ')
      cutoff.peak = 10
      cat(length(which(ss >= cutoff.peak)), 'peaks selected with minimum read of the highest peak -- ', cutoff.peak,  '\n')
      #abline(v= log10(cutoff.peak), col = 'red', lwd = 2.0)
      
      #xx = dds[ss<= cutoff.peak, ]
      dds <- dds[ss >= cutoff.peak, ]
      
    }
    
    dds <- estimateSizeFactors(dds)
    
    fpm = fpm(dds, robust = TRUE)
  }
 
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  conds = as.character(unique(design$conds))
  conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", "Embryo_Stage44_distal",
            "Mature_UA", "Mature_LA", "Mature_Hand",
            "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", "BL_UA_13days_distal")
  cs = c('Em40', 'Em44.P', 'Em44.D', 'mUA', 'mLA', 'mHand', 'BL.UA.D5', 'BL.UA.D9', 'BL.UA.D13.P', 'BL.UA.D13.D')
  
  c = match(design$conds, conds)
  
  o1 = order(c)
  c = c[o1]
  cpm = log2(fpm[, o1] +1)
  
  for(n in 1:length(peak.name))
  {
    jj = which(rownames(cpm) == peak.name[n])
    if(length(jj) == 1){
      if(!is.null(mains)){
        plot(c, cpm[jj, ], xlab = NA, ylab = 'DESeq2 norm (log2)' ,cex = 1.2, pch = 16, 
             main = paste0(rownames(cpm)[jj], '--',  mains[n]))
        axis(1, at=c(1:length(cs)),labels=cs, col.axis="darkblue", las=2)
      }else{
        plot(c, cpm[jj, ], xlab = NA, ylab = 'DESeq2 norm (log2)' ,cex = 1.2, pch = 16, 
             main = paste0(rownames(cpm)[jj]))
        axis(1, at=c(1:length(cs)),labels=cs, col.axis="darkblue", las=2)
      }
     
    }
    
    
  }
  
}

spatial.peaks.test = function(x, c = c("Mature_UA", "Mature_UA", "Mature_LA", "Mature_LA"), 
                              test.Dev.Reg = FALSE, testPlot = FALSE)
{
  # x = fpkm[ii.test[1], sample.sels]; c = cc; bg.dist = fpkm.bg;
  library(qvalue)
  if(length(x) != length(c)){
    stop('nb of data is the same as nb of conditions')
  }else{
    # the main test of mature samples 
    ii1 = which(cc == 'Mature_UA')
    ii2 = which(cc == 'Mature_LA')
    ii3 = which(cc == 'Mature_Hand')
    y0 = as.numeric(x[c(ii1, ii2, ii3)])
    s0 = c(rep(1, length(ii1)), rep(2, length(ii2)), rep(3, length(ii3)))  
    fit0 = lm (y0 ~ poly(s0, degree = 2, raw = TRUE))
    fit01 = lm(y0 ~ s0)
    fit02 = lm(y0 ~ 1)
    
    bics = BIC(fit0, fit01, fit02)
    scores = bics$BIC
    scores.relavtive = scores-min(scores)
    prob.model = exp(-0.5*scores.relavtive)
    prob.model = prob.model/sum(prob.model)
    
    # test UA, LA and Hand fitting values are above backgrounds
    pred = predict(fit0)[match(unique(s0), s0)]
    #pvals = empPvals(pred, bg.dist, pool = TRUE)
    
    res = c(prob.model[3],  max(pred), min(pred), (max(pred) - min(pred)))
    names(res) = c('prob.M0', 'max', 'min', 'log2FC')
    
    if(test.Dev.Reg){
      # test distal vs promixal in embryo.stage 44
      ii1 = which(cc == 'Embryo_Stage44_proximal')
      ii2 = which(cc == 'Embryo_Stage44_distal')
      y1 = as.numeric(x[c(ii1, ii2)])
      s1 = c(rep(1, length(ii1)), rep(3, length(ii2)))  
      fit1 = lm (y1 ~ s1)
      pval1 = summary(fit1)$coefficients[2, 4]
      
      nb1.bg = sum(c(mean(x[ii1]), mean(x[ii2])) < cutoff.bg) 
      
      ii1 = which(cc == 'BL_UA_13days_proximal')
      ii2 = which(cc == 'BL_UA_13days_distal')
      y2 = as.numeric(x[c(ii1, ii2)])
      s2 = c(rep(1, length(ii1)), rep(3, length(ii2)))  
      fit2 = lm (y2 ~ s2)
      pval2 = summary(fit2)$coefficients[2, 4]
      #nb2.bg = sum(c(mean(x[ii1]), mean(x[ii2])) < cutoff.bg)
    }
    
    if(testPlot){
      plot(s0, y0, cex = 1, ylim = range(x))
      points(s0, predict(fit0), type = 'l', col = 'blue')
      points(s0, predict(fit01), type = 'l', col = 'orange', lty = 1)
      abline(h = mean(y0), lty = 1, col = 'red', lwd = 2.0)
      points(s1, y1, cex = 1, pch = 2)
      points(s2, y2, cex = 1, pch =3)
      
    }
    
    return(res)
    
  }
  
}
static.peaks.test = function(x, c = rep(c(1:10), each = 2), testPlot = FALSE)
{
  # x = fpm[ii.test[1], sample.sels]; c = match(design$conds, conds);
  
  if(length(x) != length(c)){
    stop('nb of data is the same as nb of conditions')
  }else{
    #library("mgcv")
    library('gam')
    o1 = order(c)
    c = c[o1]
    x = x[o1]
    #names = names(x)
    x = as.numeric(x)
    
    # model null static 
    fit0 = lm(x ~ 1)
    #rss0 = sum((x - mean(x))^2)
    # definition of BIC https://en.wikipedia.org/wiki/Bayesian_information_criterion
    # be careful !!! nb of parameters k including intercept, slope and constant variance
    #bic0 = length(x)*log(rss0/length(x)) + 2*log(length(x)) 
      
    # alternative model
    
    #dat = data.frame(x =x, c = c, stringsAsFactors = FALSE)
    fit <- gam(x ~ s(c, df=3), family = gaussian)
    # summary(fit)
    #plot(fit)
    #rss1 = sum(fit$residuals^2)
    #bic1 = length(x)*log(rss1/length(x)) +  (pen.edf(fit)+2)*log(length(x))
    
    bics = BIC(fit0, fit)
    scores = bics$BIC
    scores.relavtive = scores-min(scores)
    prob.model = exp(-0.5*scores.relavtive)
    prob.model = prob.model[1]/sum(prob.model)
    
    #prob.model = prob.model[1]
    #names(prob.model) = 'prob.m0'
    
    if(testPlot){
      plot(c, x, cex = 1)
      points(c, fit$fitted.values, lwd = 1.0, type = 'l')
      abline(h = mean(x))
    }
    
    return(prob.model)
    
  }
  
}
