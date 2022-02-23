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
# prepare.annotatioin.gtf.for.peak.annotation = function(annotDir)
# {
#    
# }

select.promoters.regions = function(upstream = 2000, downstream = 2000, 
                                    ORF.type.gtf = 'Putative', promoter.select = 'all',
                                    
                                    toSave = FALSE)
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
  tss$end[jj] = tss$start[jj] + 1
  tss$start[-jj] = tss$end[-jj] - 1
  
  tss$geneSymbol = tss$geneID
  mm = match(tss$geneSymbol, gs$geneID)
  tss$geneSymbol[which(!is.na(mm))] = gs$geneSymbol.synteny.ns.nr[mm[which(!is.na(mm))]]
  tss$geneSymbol[which(is.na(tss$geneSymbol))] = tss$geneID[which(is.na(tss$geneSymbol))]
  
  colnames(tss) = c('seqnames', 'start', 'end', 'transcriptID', 'strand', 'geneID', 'geneSymbol')
  
  if(promoter.select == 'all'){
    cat('all promoters selected \n')
  }else{
    load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
    hkgs = controls.tissue[which(controls.tissue$tissues == 'housekeeping'), ]
    
    if(promoter.select == 'nonhousekeeping'){
      
      tss = tss[is.na(match(tss$geneID, hkgs$geneIDs)), ]
      cat('non-housekeeping promoters selected :', nrow(tss),  '\n')
      
    }else{
      if(promoter.select == 'housekeeping'){
        
        tss = tss[!is.na(match(tss$geneID, hkgs$geneIDs)), ]
        cat('housekeeping promoters selected', nrow(tss), '\n') 
        
      }else{
        cat('unknown promoter selection \n')
      }
    }
  }
  
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
  
  peaks$promoters[xx] = 'promoters'
  
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
  
  #saveRDS(amex, file = paste0(annotDir, 'TxDb_ax6_UCSC_2021_01_26_genes.putative.full.length.rds'))
  
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
# Section : Normalization and batch correction
# 
########################################################
########################################################
Global.Normalization.BatchCorrect = function(design, dds)
{
  
  ##########################################
  # QCs again for embryo and mature samples
  # Control the peak quality by checking the house-keeping genes and other-tissue specific genes, 
  # GO-enrichment analysis
  ##########################################
  QC.PLOT = FALSE
  if(QC.PLOT){
    pdfname = paste0(resDir, "/atacseq_Embryo_Mature_QCs_pval6.pdf")
    pdf(pdfname, width = 12, height = 10)
    
    Check.RNAseq.Quality(read.count=counts(dds)[c(1:12000),], design.matrix = data.frame(design$SampleID, 
                                                                                         design$conds, design$batch))
    dev.off()
    
  }
  
  
  d <- DGEList(counts=counts(dds), group=design$conds)
  tmm <- calcNormFactors(d, method='TMM')
  tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  #rm(tmm)
  rm(d)
  #fpm = log2(tmm + 1)
  make.pca.plots(tmm[,which(design$batch == '2020')], ntop = 1000, conds.plot = 'all')
  
  design$batch[grep('749', design$SampleID)] = '2019'
  
  table(design$condition, design$batch)
  
  bc = as.factor(design$batch)
  mod = model.matrix(~ as.factor(conds), data = design)
  
  fpm.bc = ComBat(dat=tmm, batch=bc, mod=mod, par.prior=FALSE, ref.batch = NULL) # No reference is better for some reasons    
  
  
  xx = tmm[, grep('Mature', colnames(fpm.bc))]
  yy = fpm.bc[, grep('Mature', colnames(fpm.bc))]
  
  # fpm.bc = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
  make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
  make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'Dev.Mature')
  make.pca.plots(fpm.bc, ntop = 1000, conds.plot = 'Mature')
  
  fpm = fpm.bc
  
  rm(fpm.bc)
  saveRDS(fpm, file = paste0(RdataDir, '/fpm_TMM_combat.rds')) 
  
  
}

########################################################
########################################################
# Section : promoter peak analysis
# 1) Compare the gene expression and promoter peaks
# 2) gene enrichment analysis to check if those promoters were limb-development and limb related 
########################################################
########################################################
Test.promoter.openness.enrichment = function(res, pp, bg)
{
  ##########################################
  # first, annotate peaks as promoter, hs.promoter, enhancers
  ##########################################
  xx = res
  xx$annotation[grep('Promoter', xx$annotation)] = 'Promoter'
  xx$annotation[which(xx$annotation != 'Promoter')] = 'Enhancer'
  #xx$annotation[grep('Intron', xx$annotation)] = 'Intron'
  #xx$annotation[grep('Exon', xx$annotation)] = 'Exon'
  #xx$annotation[grep('Downstream', xx$annotation)] = 'Downstream'
  
  xx$gene = xx$geneId
  xx$gene = sapply(xx$gene, function(x) {x = unlist(strsplit(as.character(x), '[|]')); x[length(x)]})
  
  annot = readRDS(paste0(annotDir, 
            'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
  
  mm = match(xx$gene, annot$transcriptid_gtf.transcript)
  xx$gene = annot$geneID_gtf.gene[mm]
  
  
  # double check the promoter peak annotation 
  xx$annot.new = 'enhancer'
  promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
  
 
  xx$annot.new[overlapsAny(pp, promoters, ignore.strand = TRUE)] = 'promoter.nohkg'  
  kk = which(xx$annot.new == 'Promoter' & xx$annotation != 'Promoter' )
  
  promoters.hkg = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', 
                                           promoter.select = 'housekeeping')
  
  xx$annot.new[overlapsAny(pp, promoters.hkg, ignore.strand = TRUE)] = 'promoter.hkg'  
  
  
  #saveRDS(xx, file = paste0(RdataDir, '/peak_annotation_enhancers_promoters_nohkg_hkg.rds'))
  xx = readRDS(file = paste0(RdataDir, '/peak_annotation_enhancers_promoters_nohkg_hkg.rds'))
  
  thresholds = seq(2, 4, by = 0.1)
  opens = matrix(NA, ncol = 5, nrow = length(thresholds))
  colnames(opens) = c('threshold.bg', 'all', 'enhancer', 'promoter.hs', 'promoter.nonhs')
  opens[,1] = thresholds
  
  for(n in 1:length(thresholds))
  {
    cutoff = thresholds[n]
    opens[n, 2] = length(which(xx$min > cutoff))
    opens[n, 3] = length(which(xx$min[which(xx$annot.new == 'enhancer')] > cutoff))
    opens[n, 4] = length(which(xx$min[which(xx$annot.new == 'promoter.hkg')] > cutoff))
    opens[n, 5] = length(which(xx$min[which(xx$annot.new == 'promoter.nohkg')] > cutoff))
  }
  
  opens = data.frame(opens)
  
  library("tidyverse")
  df <- opens %>%
    select(threshold.bg, all, enhancer, promoter.hs, promoter.nonhs) %>%
    gather(key = "annotation", value = "nb.peaks", -threshold.bg)
  
  head(df)
  
  ggplot(df, aes(x = threshold.bg, y = nb.peaks)) + 
    geom_line(aes(color = annotation)) 
  
  # ggplot(opens, aes(x = thresholds, y = nb.open)) +
  #   geom_point(color = 'blue') + geom_line() +
  #   geom_vline(xintercept = 3.0, color = 'red') +
  #   ggtitle('nb of peaks above the background in function of thresholds') +
  
  aa = c(table(xx$annot.new), table(xx$annot.new[which(xx$min>3)]))
  aa = c(aa, sum(table(xx$annot.new)), sum(table(xx$annot.new[which(xx$min>3)])))
  names(aa)[7:8] = 'all'
  aa = data.frame(names(aa), aa, c(rep(c('total.nb.peaks', 'nb.peak.above.background.all.conditions'), each = 3),
                                   'total.nb.peaks', 'nb.peak.above.background.all.conditions'))
  colnames(aa) = c('peakAnnot', 'nb.peaks', 'groups')
  #df = df[which(df$threshold.bg == 3), ]
  #df$nb.peaks.total = c(nrow(xx), length(which(xx$annot.new == 'enhancer')), 
  #                      length(which(xx$annot.new == 'promoter.hkg')),
  #                      length(which(xx$annot.new == 'promoter.nohkg')))
  
  ggplot(data=aa, aes(x=peakAnnot, y=nb.peaks, fill=groups)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    theme_minimal() 
  
  
}

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
    conds.sel = c('Embryo_Stage40',  'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
                  'Mature_UA', 'Mature_LA', 'Mature_Hand', 'HEAD', 
                  'BL_UA_5days',  'BL_UA_9days', 'BL_UA_13days_proximal',  'BL_UA_13days_distal')
  }
  
  if(conds.plot == 'Dev.Regeneration')
  {
    conds.sel = c('Embryo_Stage40',  'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
                  'BL_UA_5days',  'BL_UA_9days', 'BL_UA_13days_proximal',  'BL_UA_13days_distal')
  }
  
  if(conds.plot == 'Mature'){
    conds.sel = c('HEAD', 'Mature_UA', 'Mature_LA', 'Mature_Hand')
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

# normalize the peak width to have fpkm
FPKM.normalization.for.batchCorrection = function(fpm)
{
  library(preprocessCore)
  
  peakNames = rownames(fpm)
  peakNames = gsub('bg_', '', peakNames)
  peakNames = gsub('_', '-', peakNames)
  
  pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  ll = width(pp)
  
  fpkm = fpm.bc
  for(n in 1:ncol(fpkm))
  {
    fpkm[,n] = fpkm[,n]/ll*10^3
  }
  
  fpkm.qn = normalize.quantiles(fpkm)
  colnames(fpkm.qn) = colnames(fpkm)
  rownames(fpkm.qn) = rownames(fpkm)
  fpkm = fpkm.qn
  #rm(fpm.qn)
  make.pca.plots(fpkm.qn, ntop = 5000, conds.plot = 'all')
  
  rm(fpkm.qn)
  
  save(fpm, fpkm, file = paste0(RdataDir, '/fpm_TMM_combat_fpkm_quantileNorm.Rdata'))
  
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
  scRNADir = '/Volumes/groups/tanaka/People/current/jiwang/projects/limbRegeneration_scRNA/raw_NGS/axolotl/Gerber_2018/'
  scRNAseq.counts = read.delim(file = paste0(scRNADir, 'Fluidigm_C1/nf_out/featureCounts/merged_gene_counts.txt'), header = TRUE)
  
  metadata = read.delim(file = paste0(scRNADir, 'Metadata_EBI.ENA_filereport_read_run_PRJNA416091_tsv.txt'), header = TRUE)
  metadata = data.frame(metadata$run_accession, metadata$sample_title)
  
  cellannot = read.csv(file = paste0(scRNADir, 'aaq0681_TableS7.csv'))
  cellannot = cellannot[, c(1:3)]
  
  mm = match(cellannot$cell_id, metadata$metadata.sample_title)
  metadata = data.frame(cellannot, metadata[mm, ], stringsAsFactors = FALSE)
  
  #metadata = data.frame(metadata[which(metadata$Technology == 'Fluidigm C1'), ])
  #metadata$Sample = as.character(metadata$Sample)
  
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


pseudo.bulk.by.pooling.scRNAseq_fibroblastCells.Dev = function()
{
  scRNADir = '/Volumes/groups/tanaka/People/current/jiwang/projects/limbRegeneration_scRNA/raw_NGS/axolotl/Gerber_2018/'
  scRNAseq.counts = read.delim(file = paste0(scRNADir, 'Fluidigm_C1/nf_out/featureCounts/merged_gene_counts.txt'), header = TRUE)
  
  metadata = read.delim(file = paste0(scRNADir, 'Metadata_EBI.ENA_filereport_read_run_PRJNA416091_tsv.txt'), header = TRUE)
  metadata = data.frame(metadata$run_accession, metadata$sample_title)
  
  cellannot = read.csv(file = paste0(scRNADir, 'aaq0681_TableS7.csv'))
  cellannot = cellannot[, c(1:3)]
  
  mm = match(cellannot$cell_id, metadata$metadata.sample_title)
  metadata = data.frame(cellannot, metadata[mm, ], stringsAsFactors = FALSE)
  
  save(metadata, scRNAseq.counts, file = paste0(RdataDir, '/Gerber_2018_Fluidigm_C1.Rdata'))
  
  ##########################################
  # reload the metadata and counts 
  ##########################################
  load(file = paste0(RdataDir, '/Gerber_2018_Fluidigm_C1.Rdata'))
  
  counts = scRNAseq.counts[, -1]
  rownames(counts) = scRNAseq.counts$ENSEMBL_ID
  
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
  
  ##########################################
  # check scRNA-seq data 
  ##########################################
  library(dplyr)
  library(Seurat)
  library(patchwork)
  
  aa = CreateSeuratObject(counts = counts, project = "limb_regeneration", assay = 'RNA', meta.data = as.data.frame(metadata),
                          min.cells = 10, min.features = 500)
  
  VlnPlot(aa, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  aa <- RunPCA(aa, features = VariableFeatures(object = aa))
  ElbowPlot(aa)
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, resolution = 0.5)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.2)
  p1 = DimPlot(aa, reduction = "umap", group.by = 'timepoint')
  p2 = DimPlot(aa, reduction = "umap", group.by = 'batch')
  
  p1 + p2
  
  ggsave(paste0(resDir, "/Overview_Gerber2018_Fluidigm.C1_batches.pdf"), width = 12, height = 8)
  
  
  ## find DE genes between mUA and stage 40 and stage 44
  Idents(aa) = aa$timepoint
  DE.genes1 = FindMarkers(aa, ident.1 = "Stage40", ident.2 = "0dpa", logfc.threshold = 0.5)
  DE.genes2 = FindMarkers(aa, ident.1 = "Stage44", ident.2 = "0dpa", logfc.threshold = 0.5)
  
  DE.genes1 = DE.genes1[which(DE.genes1$p_val_adj < 0.05), ]
  DE.genes2 = DE.genes2[which(DE.genes2$p_val_adj < 0.05), ]
  
  DE.genes = unique(c(rownames(DE.genes1), rownames(DE.genes2)))
  
  # check the difference between stage40 and stage44
  aa = subset(aa, cells = colnames(aa)[which(aa$condition == 'Stage40'|aa$condition == 'Stage44')])
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  aa <- RunPCA(aa, features = VariableFeatures(object = aa))
  ElbowPlot(aa)
  
  aa <- FindNeighbors(aa, dims = 1:15)
  aa <- FindClusters(aa, resolution = 0.5)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.2)
  p1 = DimPlot(aa, reduction = "umap", group.by = 'timepoint')
  p2 = DimPlot(aa, reduction = "umap", group.by = 'batch')
  
  p1 + p2
  
  
  ##########################################
  # pool scRNA-seq data to have pseudo-bulk 
  ##########################################
  #sels = which(metadata$timepoint == '0dpa'| metadata$timepoint == 'Stage40' | metadata$timepoint == 'Stage44')
  #metadata = metadata[sels, ]
  raw = counts
  conds = c('0dpa', 'Stage40', 'Stage44')
  pseudo = matrix(NA, ncol = length(conds), nrow = nrow(raw))
  colnames(pseudo) = conds
  rownames(pseudo) = rownames(raw)
  
  for(n in 1:length(conds))
  {
    cat(n, ' : ', conds[n], ' -- ')
    jj = which(metadata$condition == conds[n])
    
    mm = match(metadata$cell_id[jj], colnames(raw))
    if(length(which(is.na(mm)))>0) {
      cat('some cells lost \n')
    }else{
      cat(length(mm), ' cell found \n')
      xx = as.matrix(raw[, mm])
      xx[is.na(xx)] = 0
      pseudo[,n] = apply(xx, 1, sum)
    }
  }
  
  #pseudo[, 2] = pseudo[,2] + pseudo[, 3]
  #pseudo = pseudo[, c(1, 2)]
  #colnames(pseudo) = c('mUA', 'stage40.44')
  
  annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  
  mm = match(rownames(pseudo), annot$geneID)
  ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
  ggs[is.na(mm)] = rownames(pseudo)[is.na(mm)]
  rownames(pseudo) = ggs
  
  require(DESeq2)
  condition = factor(c('mUA', 'stage40', 'stage44'))
  dds <- DESeqDataSetFromMatrix(pseudo, DataFrame(condition), design = ~ condition)
  
  # filtering with stage40/44 samples
  ss = rowSums(counts(dds)[, c(2:3)])
  
  hist(log10(ss), breaks = 100)
  dds = dds[which(ss > 20), ]
  
  dds = estimateSizeFactors(dds)
  
  ss = rowSums(counts(dds))
  length(which(ss > quantile(ss, probs = 0.75)))
  
  saveRDS(dds, 
       file = paste0(RdataDir, '/pseudoBulk_scRNAcellPooling_FluidigmC1_mUA_stage40_44.rds'))
  
  ss = rowMeans(counts(dds))
  dd0 = dds[ss > quantile(ss, probs = 0.75), ]
  dd0 = estimateSizeFactors(dd0)
  
  fpm = fpm(dds, robust = TRUE)
  fpm = data.frame(log2(fpm + 2^-7))
  
  
  rr1 = fpm[, 2] - fpm[, 1]
  rr2 = fpm[, 3] - fpm[, 1]
  rr = cbind(rr1, rr2)
  rr = apply(as.matrix(rr), 1, function(x) {x[which(abs(x) == max(abs(x)))][1]})
  
  fpm$log2fc = rr
  
  load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
  hs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
  ctl =  controls.tissue$geneIDs[which(controls.tissue$tissues  != 'housekeeping')]
  ggs = sapply(rownames(fpm), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
   
 
  colnames(fpm)[1] = 'mUA'
  
  fpm$genetype = 'devGene'
  
  mm = match(ggs, hs)
  xx = fpm[!is.na(mm), ]
  fpm$genetype[!is.na(mm)] = 'hs'
  
  mm = match(ggs, ctl)
  fpm$genetype[!is.na(mm) & fpm$genetype == 'devGene'] = 'ctl'
  
  par(mfrow = c(2, 1))
  hist(fpm[, 2], breaks = 100, main = 'log2 expression of stage40.44');abline(v = 2)
  hist(fpm[, 3], breaks = 100, main = 'log2 expression of stage40.44');abline(v = 2)
  
  fpm$genetype[which(fpm[,2]<2 & fpm[, 3]<2 & fpm$genetype == 'devGene')] = 'dev.lowlyExp'
  fpm$DEgene = NA
  fpm$DEgene[match(DE.genes, ggs)] = 1
  
  colnames(fpm)[6] = 'DEgene_seuratFindMarker'
  
  save(dds, fpm,
       file = paste0(RdataDir, '/pseudoBulk_scRNAcellPooling_FluidigmC1_stage40.44.mUA_dev_geneSelection.Rdata'))
  
  
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

spatial.peaks.test = function(cpm, c = c("Mature_UA", "Mature_UA", "Mature_LA", "Mature_LA"), model.selection = TRUE,
                              test.Dev.Reg = FALSE, testPlot = FALSE)
{
  # cpm = fpm[, sample.sels];  c = cc; 
  library(edgeR)
  library(qvalue)
  
  if(ncol(cpm) != length(c)){
    stop('ncol of cpm is NOT the same as nb of conditions')
  }else{
    # the main test of mature samples
    logCPM = cpm
    
    f = factor(c, levels= c('Mature_UA', 'Mature_LA', 'Mature_Hand'))
    
    mod = model.matrix(~ 0 + f)
    colnames(mod) = c('Mature_UA', 'Mature_LA', 'Mature_Hand')
   
    #To make all pair-wise comparisons between the three groups one could proceed
    fit <- lmFit(logCPM, mod)
    contrast.matrix <- makeContrasts(Mature_LA - Mature_UA, 
                                     Mature_Hand - Mature_UA, 
                                     Mature_Hand - Mature_LA, levels=mod)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    #fit <- lmFit(logCPM, mod)
    #fit <- eBayes(fit, trend=TRUE)
    #topTable(fit, coef=ncol(mod))
    
    #f <- factor(rep(c('mUA', 'mLA', 'mHand'), each = 3), levels=c("mUA","mLA","mHand")) 
    #design <- model.matrix(~0+f)
    #colnames(design) <- c("mUA","mLA","mHand")
    
    #results <- decideTests(fit2)
    
    res = data.frame(fit2$p.value)
    colnames(res) = paste0(c('mLA.vs.mUA', 'mHand.vs.mUA', 'mHand.vs.mLA'), '.pval')
    
    xx = topTable(fit2, coef = 1, number = nrow(res))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '.mLA.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    xx = topTable(fit2, coef = 2, number = nrow(res))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '.mHand.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    xx = topTable(fit2, coef = 3, number = nrow(res))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '.mHand.vs.mLA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    res$fdr.mean = apply(as.matrix(res[, grep('adj.P.Val', colnames(res))]), 1, function(x) return(mean(-log10(x))))
    res$logFC.mean =  apply(as.matrix(res[, grep('logFC', colnames(res))]), 1, function(x) return(mean(abs(x))))
    
    if(model.selection){
      library(lmtest)
      ii1 = which(c == 'Mature_UA')
      ii2 = which(c == 'Mature_LA')
      ii3 = which(c == 'Mature_Hand')
      
      res0 = matrix(NA, ncol = 7, nrow = nrow(cpm))
      colnames(res0) = c('prob.M0', 'max', 'min', 'log2FC', 'pval.lrt', 'res2.mean', 'res2.max')
      for(n in 1:nrow(cpm))
      {
        # n = 1
        if(n%%1000 == 0) cat(n, '\n')
        y0 = as.numeric(cpm[n, c(ii1, ii2, ii3)])
        tt = c(rep(1, length(ii1)), rep(2, length(ii2)), rep(3, length(ii3)))
        fit0 = lm (y0 ~ poly(tt, degree = 2, raw = FALSE))
        # #fit01 = lm(y0 ~ tt)
        fit02 = lm(y0 ~ 1)
        #
        lrt = lrtest(fit0, fit02)
        resid2 = fit0$residuals^2
        
        bics = BIC(fit0, fit02)
        scores = bics$BIC
        scores.relavtive = scores-min(scores)
        prob.model = exp(-0.5*scores.relavtive)
        prob.model = prob.model/sum(prob.model)

        # test UA, LA and Hand fitting values are above backgrounds
        pred = predict(fit0)[match(unique(tt), tt)]
        # #pvals = empPvals(pred, bg.dist, pool = TRUE)
        #
        res0[n,] = c(prob.model[2],  max(pred), min(pred), (max(pred) - min(pred)), lrt$`Pr(>Chisq)`[2], mean(resid2), max(resid2))
       
      }
      
    }
    
    res = data.frame(res, res0, stringsAsFactors = FALSE)
    
    if(test.Dev.Reg){
      # test distal vs promixal in embryo.stage 44
      ii1 = which(cc == 'Embryo_Stage44_proximal')
      ii2 = which(cc == 'Embryo_Stage44_distal')
      y1 = as.numeric(x[c(ii1, ii2)])
      s1 = c(rep(1, length(ii1)), rep(3, length(ii2)))  
      fit1 = lm (y1 ~ s1)
      pval1 = summary(fit1)$coefficients[2, 4]
      pred1 = fit1$fitted.values[match(unique(s1), s1)]
      
      res1 = c(min(pred1), (max(pred1) - min(pred1)), pval1)
      names(res1) = paste0(c('min', 'log2FC', 'pval'), '.embryo')
      #nb1.bg = sum(c(mean(x[ii1]), mean(x[ii2])) < cutoff.bg) 
      
      ii1 = which(cc == 'BL_UA_13days_proximal')
      ii2 = which(cc == 'BL_UA_13days_distal')
      y2 = as.numeric(x[c(ii1, ii2)])
      s2 = c(rep(1, length(ii1)), rep(3, length(ii2)))  
      fit2 = lm (y2 ~ s2)
      pval2 = summary(fit2)$coefficients[2, 4]
      
      pred2 = fit2$fitted.values[match(unique(s2), s2)]
      
      res2 = c(min(pred2), (max(pred2) - min(pred2)), pval2)
      names(res2) = paste0(c('min', 'log2FC', 'pval'), '.BL')
      
      res = c(res, res1, res2)
      
    }
    
    if(testPlot){
      plot(tt, y0, cex = 1, ylim = range(x), xlab = 'time points', 
           main = paste0('prob.M1 (', signif(prob.model[1], d = 2), ') -- prob.M1 (', 
                         signif(prob.model[2], d = 2), ')'))
      #points(tt, predict(fit0), type = 'l', col = 'blue')
      newtt = seq(0, 7, by = 0.2) 
      points(newtt, predict(fit0, newdata = data.frame(tt = newtt)), type = 'l', col = 'orange', lty = 1, lwd = 2.0)
      abline(h = mean(y0), lty = 1, col = 'red', lwd = 2.0)
    }
    
    return(res)
    
  }
  
}

# compare the development peaks with mUA 
dev.mUA.peaks.test = function(cpm = cpm, c = cc, model.selection = FALSE)
{
  # cpm = fpm[, sample.sels];  c = cc; 
  library(edgeR)
  library(qvalue)
  
  if(ncol(cpm) != length(c)){
    stop('ncol of cpm is NOT the same as nb of conditions')
  }else{
    # the main test of mature samples
    logCPM = cpm
    
    f = factor(c, levels= c('Mature_UA', 'Mature_Hand', 'Embryo_Stage40', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal'))
    
    mod = model.matrix(~ 0 + f)
    colnames(mod) = c('Mature_UA', 'Mature_Hand', 'Embryo_Stage40', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal')
    
    #To make all pair-wise comparisons between the three groups one could proceed
    fit <- lmFit(logCPM, mod)
    contrast.matrix <- makeContrasts(Embryo_Stage40 - Mature_UA, 
                                     Embryo_Stage40 - Mature_Hand, 
                                     Embryo_Stage44_proximal - Mature_UA,
                                     Embryo_Stage44_distal - Mature_Hand, 
                                     levels=mod)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    res = data.frame(fit2$p.value)
    colnames(res) = paste0(c('E40.vs.mUA', 
                             'E40.vs.mHand',
                             'E44prox.vs.mUA', 
                             'E44dist.vs.mHand'), '.pval')
    
    xx = topTable(fit2, coef = 1, number = nrow(res))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_E40.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    xx = topTable(fit2, coef = 2, number = nrow(res))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_E40.vs.mHand')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    
    xx = topTable(fit2, coef = 3, number = nrow(res))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_E44prox.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    xx = topTable(fit2, coef = 4, number = nrow(res))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_E44dist.vs.mHand')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    res$fdr.mean = apply(as.matrix(res[, grep('adj.P.Val', colnames(res))]), 1, function(x) return(mean(-log10(x))))
    res$logFC.mean =  apply(as.matrix(res[, grep('logFC', colnames(res))]), 1, function(x) return(mean(abs(x))))
    
    if(model.selection){
      library(lmtest)
      ii1 = which(c == 'Mature_UA')
      ii2 = which(c == 'Mature_LA')
      ii3 = which(c == 'Mature_Hand')
      
      res0 = matrix(NA, ncol = 7, nrow = nrow(cpm))
      colnames(res0) = c('prob.M0', 'max', 'min', 'log2FC', 'pval.lrt', 'res2.mean', 'res2.max')
      for(n in 1:nrow(cpm))
      {
        # n = 1
        if(n%%1000 == 0) cat(n, '\n')
        y0 = as.numeric(cpm[n, c(ii1, ii2, ii3)])
        tt = c(rep(1, length(ii1)), rep(2, length(ii2)), rep(3, length(ii3)))
        fit0 = lm (y0 ~ poly(tt, degree = 2, raw = FALSE))
        # #fit01 = lm(y0 ~ tt)
        fit02 = lm(y0 ~ 1)
        #
        lrt = lrtest(fit0, fit02)
        resid2 = fit0$residuals^2
        
        bics = BIC(fit0, fit02)
        scores = bics$BIC
        scores.relavtive = scores-min(scores)
        prob.model = exp(-0.5*scores.relavtive)
        prob.model = prob.model/sum(prob.model)
        
        # test UA, LA and Hand fitting values are above backgrounds
        pred = predict(fit0)[match(unique(tt), tt)]
        # #pvals = empPvals(pred, bg.dist, pool = TRUE)
        #
        res0[n,] = c(prob.model[2],  max(pred), min(pred), (max(pred) - min(pred)), lrt$`Pr(>Chisq)`[2], mean(resid2), max(resid2))
        
      }
      
      res = data.frame(res, res0, stringsAsFactors = FALSE)
      
    }
    
    return(res)
    
  }
  
}

temporal.peaks.test = function(x, c = c("Mature_UA", "Mature_UA", "BL_UA_5days", "BL_UA_5days", 'BL_UA_9days', 'BL_UA_9days'), 
                               testPlot = FALSE)
{
  #  c = cc[sels]; x = cpm[1, sels]
  library(qvalue)
  library('gam')
  
  if(length(x) != length(c)){
    stop('nb of data is the same as nb of conditions')
  }else{
    
    # reorder the values with the time points
    c.uniq = unique(c)
    
    y0 = c()
    tt = c()
    for(n in 1:length(c.uniq))
    {
      # ii1 = which(cc == 'Embryo_Stage40')
      # ii2 = which(cc == 'Embryo_Stage44_proximal')
      # ii3 = which(cc == 'Mature_UA')
      # ii4 = which(cc == 'BL_UA_5days')
      # ii5 = which(cc == 'BL_UA_9days')
      # ii6 = which(cc == 'BL_UA_13days_proximal')
      jj = which(c == c.uniq[n])
      y0 = c(y0, x[jj])
      tt = c(tt, rep(n, length(jj)))
      
      # y0 = as.numeric(x[c(ii1, ii2, ii3, ii4, ii5, ii6)])
      # tt = c(rep(1, length(ii1)), rep(2, length(ii2)), rep(3, length(ii3)), 
      #        rep(4, length(ii4)), rep(5, length(ii5)), rep(6, length(ii6)))
    }
    y0 = as.numeric(y0)
    tt = as.numeric(tt)
    
    fit <- gam(y0 ~ s(tt, df=(length(c.uniq)-2)), family = gaussian)
    fit0 = lm(y0 ~ 1)
    # summary(fit)
    #plot(fit)
    #rss1 = sum(fit$residuals^2)
    #bic1 = length(x)*log(rss1/length(x)) +  (pen.edf(fit)+2)*log(length(x))
    
    bics = BIC(fit0, fit)
    scores = bics$BIC
    scores.relavtive = scores-min(scores)
    prob.model = exp(-0.5*scores.relavtive)
    prob.model = prob.model/sum(prob.model)
    
    #prob.model = prob.model[1]
    #names(prob.model) = 'prob.m0'
    
    # test UA, LA and Hand fitting values are above backgrounds
    pred = predict(fit)[match(unique(tt), tt)]
    #pvals = empPvals(pred, bg.dist, pool = TRUE)
    
    res = c(prob.model[1],  max(pred), min(pred), (max(pred) - min(pred)))
    names(res) = c('prob.M0', 'max', 'min', 'log2FC')
    
    if(testPlot){
      plot(tt, y0, cex = 1, ylim = range(x), xlab = 'time points', 
           main = paste0('prob.M1 (', signif(prob.model[1], d = 2), ') -- prob.M1 (', 
                                                           signif(prob.model[2], d = 2), ')'))
      points(tt, predict(fit0), type = 'l', col = 'blue')
      newtt = seq(0, 7, by = 0.2) 
      points(newtt, predict(fit, newdata = data.frame(tt = newtt)), type = 'l', col = 'orange', lty = 1, lwd = 2.0)
      abline(h = mean(y0), lty = 1, col = 'red', lwd = 2.0)
      
            
    }
    
    return(res)
    
  }
}



all.peaks.test = function(x, c = c("Embryo_Stage40", "Embryo_Stage40", "Mature_UA", "Mature_UA", "BL_UA_5days", "BL_UA_5days", 
                                   'BL_UA_9days', 'BL_UA_9days'), testPlot = FALSE)
{
  # x = fpm[ii.test[1], sample.sels]; c = cc; 
  if(length(x) != length(c)){
    stop('nb of data is the same as nb of conditions')
  }else{
    library(qvalue)
    # reorder the values with the time points
    ii.index = c()
    tt = c()
    cc.order = c('Embryo_Stage40', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
                 'Mature_UA', 'Mature_LA', 'Mature_Hand', 'BL_UA_5days', 'BL_UA_9days', 'BL_UA_13days_proximal', 'BL_UA_13days_distal')
    
    for(i in 1:length(cc.order))
    {
      ii.index = c(ii.index, which(c == cc.order[i]))
      tt = c(tt, rep(i, length(which(c == cc.order[i]))))
    }
    
    y0 = as.numeric(x[ii.index])
   
    library('gam')
    
    fit <- gam(y0 ~ s(tt, df=4), family = gaussian)
    fit0 = lm(y0 ~ 1)
    # summary(fit)
    #plot(fit)
    #rss1 = sum(fit$residuals^2)
    #bic1 = length(x)*log(rss1/length(x)) +  (pen.edf(fit)+2)*log(length(x))
    
    bics = BIC(fit0, fit)
    scores = bics$BIC
    scores.relavtive = scores-min(scores)
    prob.model = exp(-0.5*scores.relavtive)
    prob.model = prob.model/sum(prob.model)
    
    #prob.model = prob.model[1]
    #names(prob.model) = 'prob.m0'
    
    # test UA, LA and Hand fitting values are above backgrounds
    pred = predict(fit)[match(unique(tt), tt)]
    #pvals = empPvals(pred, bg.dist, pool = TRUE)
    
    res = c(prob.model[1],  max(pred), min(pred), (max(pred) - min(pred)))
    names(res) = c('prob.M0', 'max', 'min', 'log2FC')
    
    if(testPlot){
      plot(tt, y0, cex = 1, ylim = range(x), xlab = 'time points', 
           main = paste0('prob.M1 (', signif(prob.model[1], d = 2), ') -- prob.M1 (', 
                         signif(prob.model[2], d = 2), ')'))
      points(tt, predict(fit0), type = 'l', col = 'blue')
      newtt = seq(0, max(tt), by = 0.2) 
      points(newtt, predict(fit, newdata = data.frame(tt = newtt)), type = 'l', col = 'orange', lty = 1, lwd = 2.0)
      abline(h = mean(y0), lty = 1, col = 'red', lwd = 2.0)
    
    }
    
    return(res)
    
  }
}




# all.peaks.test = function(x, c = rep(c(1:10), each = 2), testPlot = FALSE)
# {
#   # x = fpm[ii.test[1], sample.sels]; c = match(design$conds, conds);
#   
#   if(length(x) != length(c)){
#     stop('nb of data is the same as nb of conditions')
#   }else{
#     #library("mgcv")
#     library('gam')
#     o1 = order(c)
#     c = c[o1]
#     x = x[o1]
#     #names = names(x)
#     x = as.numeric(x)
#     
#     # model null static 
#     fit0 = lm(x ~ 1)
#     #rss0 = sum((x - mean(x))^2)
#     # definition of BIC https://en.wikipedia.org/wiki/Bayesian_information_criterion
#     # be careful !!! nb of parameters k including intercept, slope and constant variance
#     #bic0 = length(x)*log(rss0/length(x)) + 2*log(length(x)) 
#       
#     # alternative model
#     
#     #dat = data.frame(x =x, c = c, stringsAsFactors = FALSE)
#     fit <- gam(x ~ s(c, df=3), family = gaussian)
#     # summary(fit)
#     #plot(fit)
#     #rss1 = sum(fit$residuals^2)
#     #bic1 = length(x)*log(rss1/length(x)) +  (pen.edf(fit)+2)*log(length(x))
#     
#     bics = BIC(fit0, fit)
#     scores = bics$BIC
#     scores.relavtive = scores-min(scores)
#     prob.model = exp(-0.5*scores.relavtive)
#     prob.model = prob.model[1]/sum(prob.model)
#     
#     #prob.model = prob.model[1]
#     #names(prob.model) = 'prob.m0'
#     
#     if(testPlot){
#       plot(c, x, cex = 1)
#       points(c, fit$fitted.values, lwd = 1.0, type = 'l')
#       abline(h = mean(x))
#     }
#     
#     return(prob.model)
#     
#   }
#   
# }

########################################################
########################################################
# big Section  : quatitative analysis of peaks: bineary 0 (no peak) and 1 (peak)  
# compare the peak numbers between development, mature and regeneration 
# QCs and first analysis based on binary data (peak, no peak)
########################################################
########################################################
ATACseq.peaks.binary.analysis = function()
{
  Compare.peaks.nb.overlapping = FALSE
  if(Compare.peaks.nb.overlapping){
    source('functions_chipSeq.R')
    peakDir = '../Data/R10723_atac/calledPeaks_pval.0.001/'
    
    peak.files = list.files(path = paste0(peakDir, 'macs2'), 
                            pattern = '*macs2_peaks.xls', full.names = TRUE)
    # filter low-quality regeneration samples
    peak.files = peak.files[grep('90392|90393|108072|108073|108070|108071', peak.files, invert = TRUE)]
    
    
    ##########################################
    # three big groups: dev, mature, regeneration
    ##########################################
    compare.dev.mature.reg = FALSE
    if(compare.dev.mature.reg){
      peaks.dev = peak.files[grep('_Stage', peak.files)]
      peaks.mature = peak.files[grep('Mature', peak.files)]
      peaks.reg = peak.files[grep('BL_UA', peak.files)]
      #bam.list = list.files(path = '../Data/R10723_atac/alignments/BAMs_uniq_rmdup', pattern = '*.bam$', full.names = TRUE)
      
      peaks.dev = merge.peaks.macs2(peaks.dev, pcutoff = 6)
      peaks.mature = merge.peaks.macs2(peaks.mature, pcutoff = 6)
      peaks.reg = merge.peaks.macs2(peaks.reg, pcutoff = 6)
      
      
      p1 = peaks.dev[overlapsAny(peaks.dev, promoter, ignore.strand = TRUE)]
      p2 = peaks.mature[overlapsAny(peaks.mature, promoter, ignore.strand = TRUE)]
      p3 = peaks.reg[overlapsAny(peaks.reg, promoter, ignore.strand = TRUE)]
      
      Save.peak.coordinates = FALSE
      if(Save.peak.coordinates){
        saveRDS(p1, file = paste0(RdataDir, '/ATACseq_embryo_peaks_within.axolotl.promoters_Granges_Amex.rds'))
        saveRDS(p2, file = paste0(RdataDir, '/ATACseq_mature_peaks_within.axolotl.promoters_Granges_Amex.rds'))
        
        xx = as.data.frame(p1)
        xx$peak.name = paste0(xx$seqnames, ":", xx$start, "_", xx$end)
        xx = xx[, c(1:3, 6, 5)]
        write.table(xx, file = paste0(peakDir, 'merge_peak_Embryo_pval.6_promoters.1000bpUp.200bpDown.bed'), 
                    sep = '\t', row.names = FALSE, 
                    col.names = FALSE, quote = FALSE)
      }
      
      # all peaks overlap between Dev and Mature 
      pp = list(peaks.dev = peaks.dev, peaks.mature= peaks.mature)
      ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='peak number and overlapping Dev vs Mature')
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
      # promoter peaks overlap between Dev  
      pp = list(peaks.dev = p1, peaks.mature=p2)
      ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='promoter peak overlapping Dev vs Mature')
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
      
      # all peaks overlap between Dev, Mature and Regeneration 
      pp = list(peaks.dev = peaks.dev, peaks.mature= peaks.mature, peaks.reg = peaks.reg)
      ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='peak number and overlapping Dev, Mature and Reg')
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
      # promoter peaks overlap between Dev  
      pp = list(peaks.dev = p1, peaks.mature=p2, peaks.reg = p3)
      ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='promoter peak overlapping Dev vs Mature')
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
    }
    
    ##########################################
    # compare mature UA, UA.day5, UA.day9   
    ##########################################
    compare.mature.UAregeneration = FALSE
    if(compare.mature.UAregeneration){
      p1 = merge.peaks.macs2(peak.files[grep('Mature_UA', peak.files)], pcutoff = 6)
      p2 = merge.peaks.macs2(peak.files[grep('BL_UA_5days', peak.files)], pcutoff = 6)
      p3 = merge.peaks.macs2(peak.files[grep('BL_UA_9days', peak.files)], pcutoff = 6)
      
      pp = list(Mature.UA = p1, BL.UA.D5= p2, BL.UA.D9 = p3)
      ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='peak number and overlapping Dev, Mature and Reg')
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
      ##########################################
      # compare stage 40, mature.UA, UA5, UA9, UA13P, UA13D (only samples from 2021) 
      ##########################################
      p0 = merge.peaks.macs2(peak.files[grep('Stage40_13', peak.files)], pcutoff = 6)
      p1 = merge.peaks.macs2(peak.files[grep('Mature_UA_13', peak.files)], pcutoff = 6)
      p10 = merge.peaks.macs2(peak.files[grep('Mature_UA_74938|Mature_UA_102655', peak.files)], pcutoff = 6)
      p11 = merge.peaks.macs2(peak.files[grep('Mature_LA', peak.files)], pcutoff = 6)
      p12 = merge.peaks.macs2(peak.files[grep('Mature_Hand', peak.files)], pcutoff = 6)
      
      p2 = merge.peaks.macs2(peak.files[grep('BL_UA_5days_13', peak.files)], pcutoff = 6)
      p3 = merge.peaks.macs2(peak.files[grep('BL_UA_9days_13', peak.files)], pcutoff = 6)
      p4 = merge.peaks.macs2(peak.files[grep('BL_UA_13days_proximal_13', peak.files)], pcutoff = 6)
      p5 = merge.peaks.macs2(peak.files[grep('BL_UA_13days_distal_13', peak.files)], pcutoff = 6)
      
      peaks.ref = list(stage40 = p0, mUA = p1, mUL = p11, mHand = p12,  UA5= p2, UA9 = p3, UA13P = p4, UA13D = p5)
      
      Save.peak.coordinates.HoxClusters = FALSE
      if(Save.peak.coordinates.HoxClusters){
        load(file = paste0(RdataDir, '/samplesDesign_readCounts.withinPeaks.pval6.Rdata'))
        
        pp = data.frame(t(sapply(counts$gene, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
        pp$strand = '*'
        pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                      start.field="X2", end.field="X3", strand.field="strand")
        
        Hoxs =GenomicRanges::union(GenomicRanges::union(HoxA, HoxD1), HoxD2)
        
        pxx = pp[overlapsAny(pp, Hoxs)]
        
        require(ChIPpeakAnno)
        require(ChIPseeker)
        
        peakAnnots = as.data.frame(annotatePeak(pxx, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript'))
        
        xx = matrix(0, ncol = length(peaks.ref), nrow = nrow(peakAnnots))
        colnames(xx) = names(peaks.ref)
        for(n in 1:ncol(xx)){
          xx[overlapsAny(pxx, peaks.ref[[n]]), n] = 1
        }
        
        xx = data.frame(peakAnnots, xx, stringsAsFactors = FALSE)
        
        cne.hoxa = data.frame(read.delim(file = '../AkaneToJingkuiShareFiles/HoxA_Candidates_axCNEs/Hoxa_neiboughr_CNEs.bed',
                                         header = FALSE))
        colnames(cne.hoxa) = c('chr', 'start', 'end', 'name')
        cne.hoxa$name = paste0(cne.hoxa$name, '.HoxA')
        
        cne.hoxd = data.frame(read.delim(file = '../AkaneToJingkuiShareFiles/HoxD_Candidates_axCNEs/Axolotl_HoxD_CNEs_NEW/combined.bed',
                                         header = FALSE))
        colnames(cne.hoxd) = c('chr', 'start', 'end', 'name')
        cne.hoxd$name = paste0(cne.hoxd$name, '.HoxD')
        
        cnes = data.frame(rbind(cne.hoxa, cne.hoxd), stringsAsFactors = FALSE)
        cnes$strand = '*'
        
        
        jj = which(cnes$end <= cnes$start)
        test = cnes$start[jj]
        cnes$start[jj] = cnes$end[jj]
        cnes$end[jj] = test
        cnes$coordinates = paste0(cnes$chr, ':', cnes$start, '-', cnes$end)
        cnes = cnes[match(unique(cnes$coordinates), cnes$coordinates), c(1:5)]
        
        cnes = makeGRangesFromDataFrame(cnes, seqnames.field=c("chr"),
                                        start.field="start", end.field="end", strand.field = 'strand',  keep.extra.columns = TRUE)
        
        kk = GenomicRanges::findOverlaps(pxx, cnes, ignore.strand = TRUE, type = 'any', select = 'all')
        kk = data.frame(kk)
        
        cnes = data.frame(cnes)
        
        xx$cnes = NA
        
        for(n in 1:nrow(kk))
        {
          # n = 1
          if(is.na(xx$cnes[kk[n, 1]])) {
            xx$cnes[kk[n, 1]] = cnes$name[kk[n, 2]]      
          }else{
            xx$cnes[kk[n, 1]] = paste0(xx$cnes[kk[n, 1]], ';', cnes$name[kk[n, 2]])      
          }
          
        }
        
        write.table(xx, file = paste0(resDir, '/HoxCluster.peaks_annotation_CNEs.txt'), 
                    sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
        
        write.table(xx[, c(1, 2, 3, 12, 5)], file = paste0(resDir, '/HoxCluster.peaks_annotation.bed'), sep = '\t', 
                    col.names = FALSE, rown)   
      }
      
      library(Vennerable)
      # did not work; there are too many groups
      ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='all peaks', fill=c(1,2,3,4,5))
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
      # annotate peaks and select intron and intergeic regions
      require(ChIPpeakAnno); require(ChIPseeker)
      xx = sapply(pp, annotatePeak, TxDb=amex, tssRegion = c(-2000, 2000))
      
      pxx = pp
      for(n in 1:length(xx)){
        # n = 1
        p0 = pp[[n]]
        x0 = xx[[n]]
        x0 = as.data.frame(xx[[n]])
        
        ii = which(overlapsAny(p0, HoxD2, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', x0$annotation))
        
        pxx[[n]] = p0[ii]
      }
      
      library(Vennerable)
      # did not work; there are too many groups
      ol.peaks <- makeVennDiagram(pxx, NameOfPeaks=names(pxx), connectedPeaks="keepAll", main='peak overlapping of HoxD late TAD ',
                                  fill=c(1,2,3,4,5))
      
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
    }
    
    ##########################################
    # genomic features distribution of peaks in all conditiosn 
    ##########################################
    Characterize.genomic.feature.distrubtion.all.samples = FALSE
    if(Characterize.genomic.feature.distrubtion.all.samples){
      p0 = merge.peaks.macs2(peak.files[grep('Stage40_13', peak.files)], pcutoff = 6)
      p00 = merge.peaks.macs2(peak.files[grep('Stage40_933', peak.files)], pcutoff = 6)
      p01 = merge.peaks.macs2(peak.files[grep('Stage44_proximal_933', peak.files)], pcutoff = 6)
      p02 = merge.peaks.macs2(peak.files[grep('Stage44_distal_933', peak.files)], pcutoff = 6)
      
      p1 = merge.peaks.macs2(peak.files[grep('Mature_UA_13', peak.files)], pcutoff = 6)
      p10 = merge.peaks.macs2(peak.files[grep('Mature_UA_74938|Mature_UA_102655', peak.files)], pcutoff = 6)
      p11 = merge.peaks.macs2(peak.files[grep('Mature_LA', peak.files)], pcutoff = 6)
      p12 = merge.peaks.macs2(peak.files[grep('Mature_Hand', peak.files)], pcutoff = 6)
      
      p2 = merge.peaks.macs2(peak.files[grep('BL_UA_5days_13', peak.files)], pcutoff = 6)
      p20 = merge.peaks.macs2(peak.files[grep('BL_UA_5days_89', peak.files)], pcutoff = 6)
      
      p3 = merge.peaks.macs2(peak.files[grep('BL_UA_9days_13', peak.files)], pcutoff = 6)
      p4 = merge.peaks.macs2(peak.files[grep('BL_UA_13days_proximal_13', peak.files)], pcutoff = 6)
      p5 = merge.peaks.macs2(peak.files[grep('BL_UA_13days_distal_13', peak.files)], pcutoff = 6)
      
      pp = list(stage40 = p0, stage40.old = p00,  stage44P = p01, stage44D = p02,
                Mature.UA = p1, Mature.UA.old = p10, Mature.LA = p11, Mature.Hand = p12,  
                UA5= p2, UA5.old = p20,  
                UA9 = p3,  UA13P = p4, UA13D = p5
      )
      
      xx = sapply(pp, annotatePeak, TxDb=amex, tssRegion = c(-2000, 2000))
      #peakAnnots = annotatePeak(p0, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
      peakAnnots = xx
      
      
      #xx = data.frame(peakAnnots)
      pdf(paste0(resDir, "/Peaks_genomicFeatures_distance_TSS.pdf"), width = 10, height = 6)
      #par(mfrow=c(1,2))
      print(plotAnnoBar(peakAnnots))
      print(plotDistToTSS(peakAnnots))
      dev.off()
      
      
    }
    
    ##########################################
    # qualitatively characterize position-restricted peaks in mature samples 
    ##########################################
    Characterize.position.restricted.peaks = FALSE
    if(Characterize.position.restricted.peaks){
      source('functions_chipSeq.R')
      
      p1 = merge.peaks.macs2(peak.files[grep('Mature_UA_13', peak.files)], pcutoff = 3, select.overlappingPeaks = TRUE)
      #xx = merge.peaks.macs2(peak.files[grep('Mature_UA_13', peak.files)], pcutoff = 6, select.overlappingPeaks = FALSE)
      p10 = merge.peaks.macs2(peak.files[grep('Mature_UA_74938|Mature_UA_102655', peak.files)], pcutoff = 3, 
                              select.overlappingPeaks = TRUE)
      p11 = merge.peaks.macs2(peak.files[grep('Mature_LA', peak.files)], pcutoff = 3, select.overlappingPeaks = TRUE)
      p12 = merge.peaks.macs2(peak.files[grep('Mature_Hand', peak.files)], pcutoff = 3, select.overlappingPeaks = TRUE)
      
      pp = list(mUA = p1, mUA.old = p10, mLA = p11, mHand = p12)
      
      library(Vennerable)
      ol.peaks <- makeVennDiagram(pp, NameOfPeaks=names(pp), connectedPeaks="keepAll", main='peak overlapping in Mature',
                                  fill=c(1:length(pp)))
      pxx = pp[c(1:2)]
      makeVennDiagram(pxx, NameOfPeaks=names(pxx), connectedPeaks="keepAll", main='peak overlapping in Mature',
                      fill=c(1:length(pxx)))
      
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
      # pool peaks in all mature samples and get a matrix indicating where they were detected
      pxx = union(union(union(p1, p10,  ignore.strand=TRUE), p11,  ignore.strand=TRUE), p12, ignore.strand=TRUE)
      pxx = reduce(pxx)
      
      peaks.ref = pp
      
      xx = matrix(0, ncol = length(peaks.ref), nrow = length(pxx))
      colnames(xx) = names(peaks.ref)
      rownames(xx) = paste0(seqnames(pxx), ':', start(pxx), '-', end(pxx))
      
      for(n in 1:ncol(xx)){
        xx[overlapsAny(pxx, peaks.ref[[n]]), n] = 1
      }
      ss = apply(xx, 1, sum)
      
      xx = xx[which(ss < 4), ]
      xx = xx[which(xx[,1] == xx[, 2]), ]
      
      ss = apply(xx, 1, sum)
      xx = as.matrix(xx)
      
      pheatmap(xx[c(1:10000), ], cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = TRUE,
               cluster_cols=FALSE)
      
    }
    
    xx = read.delim(file = 
      '../results/R10723_Rxxxx_atacseq_bowtie2.newParam_mtDNA_picardrmdup_20210208/HoxCluster.peaks_annotation_CNEs.txt',
                    sep = '\t', header = TRUE)
    
    #xx = xx[which(xx$seqnames == 'chr2p'), ]
    #ss = apply(as.matrix(xx[, c(16, 19:22)]), 1, sum)
    #jj = which(xx$UA13D == 1 & ss == 1)
    #xx = xx[jj, ]
    
    #write.table(xx, file = paste0(resDir, '/peakList_coordinates_BL.UA.D13.Distal.txt'), sep = '\t', row.names = FALSE,
    #            col.names = FALSE)
    
    names = paste0(xx$seqnames, ':', xx$start, '-', xx$end)
    
    source('Functions.R')
    
    pdfname = paste0(resDir, "/peak_profiles.pdf")
    pdf(pdfname, width = 10, height = 6)
    par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(6,3,2,0.2), tcl = -0.3)
    
    plot.peak.profiles(peak.name = names)
    
    dev.off()
    
  }
  
}

########################################################
########################################################
# Section : first analysis for subgrouped peaks, promoter, enhancers
# 
########################################################
########################################################
Grouping.peaks.for.promoters.enhancers = function()
{
  ##########################################
  # overview all peaks with heatmap for subgrouped peaks
  ##########################################
  source('Functions.R')
  
  make.heatmap.atacseq()
  
  sel.promoter = which(overlapsAny(pp, promoters, ignore.strand = TRUE) == TRUE)
  sel.promoter.dev = which(overlapsAny(pp, promoters.dev, ignore.strand = TRUE) == TRUE)
  
  ii.promoter = grep('Promoter', peakAnnots$annotation)
  sel.intron = grep('Intron ', peakAnnots$annotation)
  sel.intergen = grep('Intergenic', peakAnnots$annotation)
  
  cat(length(sel.promoter), ' peaks in promoters\n')
  cat(length(sel.promoter.dev), ' peaks in promoters.dev \n')
  
  cat(length(sel.intron), ' peaks in introns\n')
  cat(length(sel.intergen), ' peaks in intergenic regions\n')
  
  ii.HoxA = which(overlapsAny(pp, HoxA, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))
  ii.HoxD1 = which(overlapsAny(pp, HoxD1, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))
  ii.HoxD2 = which(overlapsAny(pp, HoxD2, ignore.strand = TRUE) == TRUE & grepl('Intron|Intergenic', peakAnnots$annotation))
  
  ##########################################
  #  select conditions of interest
  ##########################################
  #conds.sel = c('Embryo_Stage40_93',  'Embryo_Stage40_13', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal', 
  # 'Mature_UA_', 'Mature_LA', 
  #              'Mature_Hand', 'BL_UA_5days_13', 'BL_UA_5days_89',  'BL_UA_9days_13', 'BL_UA_13days_proximal', 
  # 'BL_UA_13days_distal')
  
  conds.sel = c('Embryo_Stage40_93',  'Embryo_Stage40_13', 'Embryo_Stage44_proximal', 'Embryo_Stage44_distal',
                'Mature_UA_13',  'Mature_LA', 'Mature_Hand', 'BL_UA_5days_13', 'BL_UA_5days_89',  'BL_UA_9days_13', 
                'BL_UA_13days_proximal',
                'BL_UA_13days_distal')
  #conds.sel = c('Mature_UA_13', 'Mature_UA_74938|Mature_UA_102655', 'Mature_LA', 'Mature_Hand')
  
  conds.sel = c('Mature_UA_13', 'BL_UA_5days_13', 'BL_UA_9days_13', 'BL_UA_13days_proximal',  'BL_UA_13days_distal')
  
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
  
  test.normalization.quantile = FALSE
  if(test.normalization.quantile){
    xx = fpm[, sample.sel]
    plot(xx[, c(1, 3)], cex = 0.1); abline(0, 1, lwd = 2.0, col ='red')
    plot(xx[, c(2, 3)], cex = 0.2); abline(0, 1, lwd = 2.0, col ='red')
    
    library(factoextra)
    res.pca <- prcomp(t(xx), scale = TRUE)
    #res.var <- get_pca_var(res.pca)
    
    fviz_pca_ind(res.pca,
                 col.ind = "cos2", # Color by the quality of representation
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
    )
    
  }
  
  Split.peaks.into.promoters.intron.integeic = FALSE
  if(Split.peaks.into.promoters.intron.integeic){
    
    ##########################################
    # promoter peaks
    ##########################################
    keep =fpm[sel.promoter.dev, sample.sel]
    
    filtering.with.signal = TRUE
    if(filtering.with.signal){
      
      ss = apply((keep), 1, max)
      
      max.cutoff = 3
      
      hist(ss, breaks = 100)
      abline(v = max.cutoff, col='red', lwd =2.0)
      cat(length(which(ss>max.cutoff)), 'peak selected\n')
      
      keep = keep[which(ss > max.cutoff), ]
      
    }
    
    keep = as.matrix(keep)
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    
    ##########################################
    # promoter introns
    ##########################################
    keep =fpm[sel.intron, sample.sel]
    
    filtering.with.signal = TRUE
    if(filtering.with.signal){
      
      ss = apply((keep), 1, max)
      
      max.cutoff = 3
      hist(ss, breaks = 100)
      abline(v = max.cutoff, col='red', lwd =2.0)
      cat(length(which(ss>max.cutoff)), '\n')
      
      keep = keep[which(ss>max.cutoff), ]
      
    }
    
    #subsample = sample(c(1:nrow(keep)), 15000)
    #keep = keep[subsample, ]
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    ##########################################
    # enhancer peaks too many and subseting is necessary
    ##########################################
    keep =fpm[sel.intergen, sample.sel]
    
    filtering.with.signal = TRUE
    if(filtering.with.signal){
      ss = apply((keep), 1, max)
      
      max.cutoff = 3.5
      hist(ss, breaks = 100)
      abline(v = max.cutoff, col='red', lwd =2.0)
      cat(length(which(ss>max.cutoff)), '\n')
      
      keep = keep[which(ss>max.cutoff), ]
      
    }
    
    #subsample = sample(c(1:nrow(keep)), 10000)
    #keep = keep[subsample, ]
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    
    ##########################################
    # peaks in introns and intergenic regions in HoxA and HoxD clusters 
    ##########################################
    keep =fpm[ii.HoxA, sample.sel]
    pheatmap(as.matrix(log2(keep + 1)), cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    keep =fpm[ii.HoxD1, sample.sel]
    pheatmap(as.matrix(log2(keep + 1)), cluster_rows=TRUE, show_rownames=FALSE, scale = 'none', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    
  }
  
  ##########################################
  # position-related mature peaks 
  ##########################################
  Test.position.related.peaks.in.mature = FALSE
  if(Test.position.related.peaks.in.mature){
    keep = fpm[, sample.sel]
    
    filtering.with.signal = TRUE
    if(filtering.with.signal){
      ss = apply((keep[, c(1:4)]), 1, mean)
      vars = apply(keep[, c(1:4)], 1, var)
      
      plot(ss, vars, cex = 0.15)
      #abline(h=c(1,2), col = 'red', lwd =2.0)
      abline(v = 1, col = 'red', lwd =2.0)
      jj = which(vars < 1.)
      
      keep = keep[jj, ]
      
      ss = apply(keep[, c(5:6)], 1, mean)
      vars = apply(keep[, c(5:6)], 1, var)
      plot(ss, vars, cex = 0.15)
      abline(h=1, col = 'red', lwd =2.0)
      
      jj = which(vars < 1)
      
      keep = keep[jj, ]
      
      ss = apply(keep[, c(7:8)], 1, mean)
      vars = apply(keep[, c(7:8)], 1, var)
      plot(ss, vars, cex = 0.15)
      
      jj = which(vars < 1)
      
      keep = keep[jj, ]
      
      ss = apply(keep, 1, mean)
      vars = apply(keep, 1, var)
      plot(ss, vars, cex = 0.15)
      
      jj = which(ss>0 & vars >1)
      
      keep = keep[jj, ]
      
      
      max.cutoff = 3.5
      hist(ss, breaks = 100)
      abline(v = max.cutoff, col='red', lwd =2.0)
      cat(length(which(ss>max.cutoff)), '\n')
      
      
    }
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
             cluster_cols=FALSE, annotation_col = df, gaps_col = ii.gaps)
    
    
  }
  
}

########################################################
########################################################
# Section : Identify potential LA- and Hand-specific genes
# to try to solve the LA and Hand sample swapping issue
########################################################
########################################################
Identify.LA.Hand.specific.genes.from.atacseq = function(xx)
{
  xx = fpm[, sample.sels]
  xx = data.frame(mUA =  apply(xx[, grep('Mature_UA', colnames(xx))], 1, mean), 
                  mLA =  apply(xx[, grep('Mature_LA', colnames(xx))], 1, mean), 
                  mHand = apply(xx[, grep('Mature_Hand', colnames(xx))], 1, mean))
  #lfc = apply(xx[, grep('Mature_LA', colnames(xx))], 1, mean) - apply(xx[, grep('Mature_Hand', colnames(xx))], 1, mean)
  xx  =data.frame(xx, geneId = res$geneId, stringsAsFactors = FALSE)
  
  #lfc = lfc[order(-abs(lfc$lfc)), ]
  
  saveRDS(xx, file = '../results/Rdata/LA_Hand_potential.specific.genes_fromATACseq.rds')
  
}





########################################################
########################################################
# Section : some handy utility functions
# 
########################################################
########################################################
check.atacseq.sample.ratios = function()
{
  aa = read.table(file = '../data/ed_AK_countStatTable_manualCleaned.txt', sep = '\t', header = TRUE)
  xx = aa[which(aa$Name != ''), ]
  colnames(xx)[2] = 'condition'
  colnames(xx)[8:10] = c('percent.required', 'percent.total.reads', 'percent.unique.reads')
  
  aa = xx
  aa$pct.total = aa$Total/sum(aa$Total)
  aa$pct.usable = aa$unique.rmdup/sum(aa$unique.rmdup)
  aa$pct.desired = as.numeric(as.character(aa$percent.required))
  plot(aa$pct.desired, aa$pct.total, log = '');
  text(aa$pct.desired, aa$pct.total, aa$condition)
  abline(0, 1, lwd = 2.0, col = 'red')
 
  aa$read.pred = aa$pct.total* 3.7*10^9/10^6
  
}

extract.stat.for.samples.manual = function()
{
  design = read.table(design_file, sep = '\t', header = TRUE)
  colnames(design)[2] = 'condition'
  design$fileName = paste0(design$condition, '_', design$sampleID)
  
  #stats = read.table(paste0(dataDir, 'nf_out/result/countStatTable.txt'), sep = '\t', header = TRUE)
  #colnames(stats)[c(1, 3)] = c('fileName', 'trimmed')
  
  #cnts = list.files(path = '../Data//R10723_atac/QCs/cnt_raw', pattern = '*.txt', full.names = TRUE)
  
  #design = design[order(design$fileName), ]
  stats = list.files(path = paste0(dataDir, '/merged_bam_QCs/BamStat'), pattern = '*_stat.txt', full.names = TRUE)
  
  #design = design[order(design$fileName), ]
  design$unique = NA
  design$usable = NA
  for(n in 1:nrow(design))
  {
    # n = 1;
    cat(n, '\n')
    #cc.files = cnts[grep(design$sampleID[n], cnts)]
    ii = grep(design$sampleID[n], stats)
    if(length(ii) == 1){
      #index = c(index, ii)
      ss = read.table(file = stats[ii], sep = '\t', header = TRUE)
      design$unique[n] = ss$uniq
      design$usable[n] = ss$rmdup_uniq
    }else{
      #index = c(index, NA)
      cat(length(ii), ' fileName Found for ', design$sampleID[n], '\n')
    }
    
    # total = 0
    # for(m in 1:length(cc.files))
    # {
    #   #cat(m, '\n')
    #   total = total + read.table(cc.files[m], sep = '\t', header = FALSE)
    # }
    # 
    # ss = read.table(files.stat[grep(design$sampleID[n], files.stat)], sep = '\t', header = TRUE)
    # samples = c(samples, as.character(ss[1, 1]))
    # stats = rbind(stats, c(total, ss[1, -1]))
    
  }
  
  design$usable = design$usable/10^6
  
  saveRDS(design, file = paste0(RdataDir, '/design_merged_technicalReplicates.rds'))
  
  write.table(design, file = paste0(resDir, '/QCs_nb.usableReads_mergedResequenced.samples.txt'), 
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
  
  #stats = data.frame(design, stats[index, ], stringsAsFactors = FALSE)
  #colnames(stats) = c('sampleID', 'samples', 'fileName', 'total',  'adapter.trimmed', 'mapped', 'chrM.rm', 'unique', 
  # 'unique.rmdup')
  
  #stats$trimming.pct = as.numeric(stats$adapter.trimmed)/as.numeric(stats$total)
  # stats$mapped.pct = as.numeric(stats$mapped)/as.numeric(stats$adapter.trimmed)
  # stats$mito.pct = as.numeric(stats$chrM.rm)/as.numeric(stats$adapter.trimmed)
  # stats$multimapper.pct = 1- as.numeric(stats$unique) / as.numeric(stats$mapped)
  # stats$dup.rate = 1.0 - as.numeric(stats$unique.rmdup)/as.numeric(stats$unique)
  # stats$pct.usable = stats$unique.rmdup / stats$total
  # 
  #stats$sample = gsub('_sorted', '', stats$sample)
  #stats = stats[, c(1, 2, 3, 6, 4, 7, 5, 8)]
  #colnames(stats)[c(5, 7)] = c('uniq.mapped', 'uniq.mapped.rmdup')
  # library(ggplot2)
  # library(dplyr)
  # 
  # xx = stats[, c(8:12)]
  # 
  # # Create data
  # data <- data.frame(
  #   name=c(rep(colnames(xx), each = nrow(xx))),
  #   value=as.numeric(unlist(xx)) %>% round(2)
  # )
  # 
  # 
  # ggplot(data, aes(x=name, y=value, fill=name)) + 
  #   geom_violin()
  # 
  # #boxplot(stats[, c(8:12)])
  # df <- apply(stats,2,as.character)
  #stats = as.data.frame(stats)
  
  write.csv(stats, file = paste0(resDir, '/R11876_CutTag_QCs_stats.csv'), row.names = FALSE)
  
  stats$usable = stats$unique.rmdup/10^6
  colnames(stats)[c(2,3)] = c('condition', 'samples')
  stats$samples = paste0(stats$condition, '_', stats$sampleID)
  
  plot(stats$usable, stats$mapped); text(stats$usable, stats$mapped, stats$samples)
  
  #save(stats, file = paste0(RdataDir, '/R11876_CutTag_samples_design_stats.Rdata'))
  save(stats, file = paste0(RdataDir, '/R11637_atacseq_samples_design_stats.Rdata'))
  
}

select.plot.top.positional.peaks = function()
{
  ##########################################
  ## select top peaks genome-wide
  ##########################################
  Select.top.peaks = FALSE
  if(Select.top.peaks){
    load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl', version.analysis, '.Rdata'))
    xx = xx[order(-xx$logFC.mean), ]
    
    yy = xx
    yy[grep('HOXA13', yy$geneId), ]
    
    # res2 here is residues in the fitting for each condition, UA, LA and Hand
    yy = yy[which(yy$max > 4 & yy$min < 3 & yy$res2.max< 1.0), ]
    
    #yy = yy[which(yy$max > 3 & yy$min < 1), ]
    
    yy = yy[order(yy$logFC.mean), ]
    
    keep = fpm[!is.na(match(rownames(fpm), rownames(yy))), sample.sels]
    gg = res$geneId[match(rownames(keep), rownames(res))]
    grep('HOXA13', gg)
    
    rownames(keep) = paste0(rownames(keep), '_', gg)
    #rownames(keep) = gg
    keep = as.matrix(keep)
    
    gg = rownames(keep)
    gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '_'))[2])
    gg = sapply(gg, function(x) unlist(strsplit(as.character(x), '[|]'))[1])
    
    pheatmap(keep, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = FALSE,
             cluster_cols=FALSE, annotation_col = df, fontsize_row = 8, gaps_col = ii.gaps,
             filename = paste0(resDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_top300_genomewide.pdf'), 
             width = 16, height = 40)
    
    if(saveTable){
      write.csv(data.frame(keep, yy, stringsAsFactors = FALSE), 
                file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_top300_genomewide.csv'), 
                quote = FALSE, row.names = TRUE)
      
    }
  }
  
}

########################################################
########################################################
# Section : plotting functions
# 
########################################################
########################################################
plotPeakAnnot_piechart = function(peakAnnots, ndigit = 2)
{
  features = c("Promoter (<=1kb)", "Promoter (1-2kb)",   
               "5' UTR","3' UTR", "1st Exon", "Other Exon", 
               "1st Intron", "Other Intron",
               "Downstream (<=300)",  "Distal Intergenic") 
  # Create Data
  data <- data.frame(peakAnnots@annoStat, stringsAsFactors = FALSE)
  data$Feature = as.character(data$Feature)
  data$Frequency = as.numeric(as.character(data$Frequency))
  data = data[match(features, data$Feature), ]
  data$Feature = features
  data$Frequency[which(is.na(data$Frequency))] = 0.0
  data$Frequency[which(data$Feature == 'Promoter (<=1kb)')] = sum(data$Frequency[grep('Promoter', data$Feature)])
  data = data[-which(data$Feature == 'Promoter (1-2kb)'), ]
  data$Feature[which(data$Feature == 'Promoter (<=1kb)')] = 'Promoter (<=2kb)'
  
  data$Feature = paste0(data$Feature, ' (', signif(data$Frequency, d =ndigit), "%)")
  
  p1 = ggplot(data, aes(x="", y=Frequency, fill = factor(Feature, levels = Feature))) +
    geom_bar(stat="identity", width=1,  color="white") +
    coord_polar("y", start=0) +
    theme_void() +
    scale_fill_brewer(palette="Set1", direction = 1)
  
  plot(p1)
  
  return(data)
  
}


########################################################
########################################################
# Section : old functions not Used 
# 
########################################################
########################################################
# promoter analysis
# grouping atac-seq peak profiles
# 1) first identify static peaks (probably most of them) with model selection 
# (linear regression with only intercept and nonlinear fitting with spline or GAM)
# 2) grouping dynamic peaks using DP-GP
# ##########################################
# all-peaks test M0 (static peaks), M1 (dynamic peaks above background), M2 (dyanmic peaks with some condtions below background)
# we will also loci that are always open across all conditions
promoter.peaks.analysis.across.mature.dev.regeneraiton = function()
{
  make.test.All.Peaks = FALSE
  if(make.test.All.Peaks){
    
    # to change for all peaks
    fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.rds'))
    design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_MatureSamples_batch2019.2020.2021.2021S.rds'))
    
    # prepare the background distribution
    fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
    fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
    rownames(fpm) = gsub('_', '-', rownames(fpm))
    
    hist(fpm.bg, breaks = 100, main = 'background distribution')
    abline(v = 1, col = 'red', lwd = 2.0)
    quantile(fpm.bg, c(0.95, 0.99))
    
    ##########################################
    ## make Granges and annotate peaks
    ##########################################
    Make.Granges.and.peakAnnotation = TRUE
    if(Make.Granges.and.peakAnnotation){
      require(ChIPpeakAnno)
      require(ChIPseeker)
      
      pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
      pp$strand = '*'
      
      save.peak.bed = FALSE
      if(save.peak.bed){
        bed = data.frame(pp[, c(1:3)], rownames(pp), 0, pp[, 4], stringsAsFactors = FALSE)
        write.table(bed, file = paste0(resDir, '/peakset_', version.analysis, '.bed'), col.names = FALSE, row.names = FALSE,
                    sep = '\t', quote = FALSE)
      }
      
      pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                    start.field="X2", end.field="X3", strand.field="strand")
      
      # annotation from ucsc browser ambMex60DD_genes_putative
      amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
      pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
      #plotAnnoBar(pp.annots)
      
      pp.annots = as.data.frame(pp.annots)
      rownames(pp.annots) = rownames(fpm)
      
      promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
      
    }
    
    
    Run.all.peaks.test = FALSE
    if(Run.all.peaks.test){
      
      # all conditions included
      conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", "Embryo_Stage44_distal",
                "Mature_UA",  "Mature_LA", "Mature_Hand", 
                "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", "BL_UA_13days_distal")
      
      # test examples
      test.examples = c('HAND2', 'FGF8', 'KLF4', 'Gli3', 'Grem1')
      #test.examples = c('Hoxa13')
      ii.test = which(overlapsAny(pp, promoters[which(!is.na(match(promoters$geneSymbol, test.examples)))]))
      #ii.Hox = which(overlapsAny(pp, Hoxs))
      #ii.test = unique(c(ii.test, ii.Hox))
      
      sample.sels = c()
      cc = c()
      for(n in 1:length(conds)) {
        kk = which(design$conds == conds[n] & design$SampleID != '136159')
        sample.sels = c(sample.sels, kk)
        cc = c(cc, rep(conds[n], length(kk)))
      }
      
      library(tictoc)
      ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
      
      source('Functions.R')
      tic() 
      res = t(apply(fpm[ii.test, sample.sels], 1, all.peaks.test, c = cc))
      res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
      toc()
      
      saveRDS(res, file = paste0(RdataDir, '/res_allPeaks_test.rds'))
      
    }else{
      res = readRDS(file = paste0(RdataDir, '/res_allPeaks_test.rds'))
      
    }
    
    jj = which(res$prob.M0 < 0.05 & res$log2FC >1)
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC), ]
    
    length(which(xx$min <3.0))
    length(which(xx$min <2.5))
    length(which(xx$min <2.))
    
    bgs = data.frame(bg = as.numeric(fpm.bg))
    
    ggplot(bgs, aes(x=bg)) + geom_histogram(binwidth = 0.25, color="darkblue", fill="lightblue") +
      theme(axis.text.x = element_text(angle = 0, size = 16)) + 
      xlab(" background signals (log2)") 
    
    Test.promoter.openness = FALSE
    if(Test.promoter.openness){
      
      source('Functions.R')
      Test.promoter.openness.enrichment(res, bg) # two different methods to test nonHS promoters have tendancy to be open across conditions
      
      #DoubleCheck.promoter.peaks.enrichment(fpm) # double check those open promoters are not due to sample contamination
      
    }
  }
  
}
