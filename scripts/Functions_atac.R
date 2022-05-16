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


make.pca.plots = function(fpm, ntop = 1000, conds.plot = 'Dev.Mature', conds.sel = NULL)
{
  library(factoextra)
  if(is.null(conds.sel)){
    
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


temporal.peaks.test = function(cpm, c = c("Mature_UA", "Mature_UA", "BL_UA_5days", "BL_UA_5days", 'BL_UA_9days', 'BL_UA_9days'), 
                               run.pairwise.comparison.edgeR = TRUE,  testPlot = FALSE)
{
  #  c = cc[sels]; x = cpm[1, sels]
  library(qvalue)
  library('gam')
  library('lmtest')
  
  if(ncol(cpm) != length(c)){
    stop('nb of data is the same as nb of conditions')
  }else{
    
    # reorder the values with the time points
    c.uniq = unique(c)
    
    res0 = matrix(NA, ncol = 7, nrow = nrow(cpm))
    colnames(res0) = c('prob.M0', 'max', 'min', 'log2FC', 'pval.lrt', 'res2.mean', 'res2.max')
    
    for(n in 1:nrow(cpm))
    {
      # n = which(rownames(cpm) == 'chr12p:6136134-6136818')
      if(n%%1000 == 0) cat(n, '\n')
      y0 = c()
      tt = c()
      for(ic in 1:length(c.uniq))
      {
        # ii1 = which(cc == 'Embryo_Stage40')
        # ii2 = which(cc == 'Embryo_Stage44_proximal')
        # ii3 = which(cc == 'Mature_UA')
        # ii4 = which(cc == 'BL_UA_5days')
        # ii5 = which(cc == 'BL_UA_9days')
        # ii6 = which(cc == 'BL_UA_13days_proximal')
        jj = which(c == c.uniq[ic])
        y0 = c(y0, cpm[n, jj])
        tt = c(tt, rep(ic, length(jj)))
        
      }
      
      y0 = as.numeric(y0)
      tt = as.numeric(tt)
      
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
      
      lrt = lrtest(fit, fit0)
      resid2 = fit0$residuals^2
      
      #prob.model = prob.model[1]
      #names(prob.model) = 'prob.m0'
      
      # test UA, LA and Hand fitting values are above backgrounds
      pred = predict(fit)[match(unique(tt), tt)]
      #pvals = empPvals(pred, bg.dist, pool = TRUE)
      
      res0[n,] = c(prob.model[2],  max(pred), min(pred), (max(pred) - min(pred)), lrt$`Pr(>Chisq)`[2], mean(resid2), max(resid2))
      
    }  
    
    
    if(testPlot){
      #y0[3] = -0.5
      # pdfname = paste0(resDir, "/peak_profiles.pdf")
      # pdf(pdfname, width = 10, height = 6)
      # par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(6,3,2,0.2), tcl = -0.3)
      
      plot(tt, y0, cex = 2.0, ylim = range(y0), xlab = '', pch = 16, col = 'blue',
           ylab = 'log2(peak signals)',  xaxt="n",
           main = paste0('prob.M0 (', signif(prob.model[1], d = 2), ') -- prob.M1 (', 
                                                           signif(prob.model[2], d = 2), ')'))
      points(tt, predict(fit0), type = 'l', col = 'darkgray', lwd = 6.0)
      newtt = seq(range(tt)[1], range(tt)[2], by = 0.1) 
      points(newtt, predict(fit, newdata = data.frame(tt = newtt)), type = 'l', col = 'darkgreen', lty = 1, lwd = 6.0)
      legend('topleft', lwd = c(6, 6), col = c('darkgray', 'darkgreen'), legend = c('M0', 'M1'), bty = 'n' )
      #abline(h = mean(y0), lty = 1, col = 'red', lwd = 2.0)
      axis(1, at=c(1:8), labels= c('s40', 's44_prox', 's44_dist', 'mUA', '5dpa', '9dpa', '13dpa_prox', '13dpa.dist'), 
           col.axis="blue", las=1, lty = 2.0, cex.axis = 1.2)
      
      #dev.off()
      
    }
    
  }
  
  res = data.frame(res0, stringsAsFactors = FALSE)
  rownames(res) = rownames(cpm)
  rm(res0)
  
  if(run.pairwise.comparison.edgeR){
    library(edgeR)
    logCPM = cpm
    
    clevels = c('Mature_UA', setdiff(unique(c), 'Mature_UA'))
    f = factor(c, levels= clevels)
    
    mod = model.matrix(~ 0 + f)
    colnames(mod) = clevels
    
    #To make all pair-wise comparisons between the three groups one could proceed
    fit <- lmFit(logCPM, mod)
    contrast.matrix <- makeContrasts(Embryo_Stage40 - Mature_UA, 
                                     Embryo_Stage44_proximal - Mature_UA, 
                                     Embryo_Stage44_distal - Mature_UA,
                                     BL_UA_5days - Mature_UA,
                                     BL_UA_9days - Mature_UA,
                                     BL_UA_13days_proximal - Mature_UA,
                                     BL_UA_13days_distal - Mature_UA,
                                     levels=mod)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    #res = data.frame(fit2$p.value)
    #colnames(res) = paste0(c('s40.vs.mUA', 
    # 's44p.vs.muA',
    # 's44d.vs.mUA', 
    # '5dpa.vs.mUA',
    # '9dpa.vs.mUA',
    # '13dpap.vs.mUA',
    # '13dpad.vs.mUA'
    # ), '.pval')
    
    xx = topTable(fit2, coef = 1, number = nrow(logCPM))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_s40.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    xx = topTable(fit2, coef = 2, number = nrow(logCPM))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_s44p.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    
    xx = topTable(fit2, coef = 3, number = nrow(logCPM))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_s44d.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    
    xx = topTable(fit2, coef = 4, number = nrow(logCPM))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_5dpa.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    xx = topTable(fit2, coef = 5, number = nrow(logCPM))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_9dpa.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    xx = topTable(fit2, coef = 6, number = nrow(logCPM))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_13dpap.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    xx = topTable(fit2, coef = 7, number = nrow(logCPM))
    xx = xx[, c(1, 4, 5)]
    colnames(xx) = paste0(colnames(xx), '_13dpad.vs.mUA')
    res = data.frame(res, xx[match(rownames(res), rownames(xx)), ])
    
    res$fdr.mean.edgeR = apply(as.matrix(res[, grep('adj.P.Val', colnames(res))]), 1, function(x) return(mean(-log10(x))))
    res$logFC.mean.edgeR =  apply(as.matrix(res[, grep('logFC', colnames(res))]), 1, function(x) return(mean(abs(x))))
    
  }
  
  return(res)
  
  
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


process.dynamic.peaks.clustering.GPDP = function(yy, res)
{
  require(data.table)
  dpgpDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/DPGP_clustering/out_10k_2d/'
  
  dpgp = read.table(file = paste0(dpgpDir, 
  'dpgp_clusters_optimal_clustering.txt'), sep = '\t', header = TRUE)
  #$'dpgp_clusters_posterior_similarity_matrix.txt'))
  dpgp$gene = sapply(dpgp$gene, function(x) {x = unlist(strsplit(as.character(x), '_')); paste0(x[1], ':', x[2], '_', x[3])})
  
  clusters = table(dpgp$cluster)
  print(clusters)
  clusters = clusters[which(clusters >= 0.02*nrow(dpgp))]
  
  # <- fread(input)
  #dpgp = read.table(file = )
  mm = match(dpgp$gene, rownames(yy))
  
  xx = yy[mm, ]
  
  cs = matrix(NA, nrow = length(clusters), ncol = ncol(xx))
  colnames(cs) = colnames(xx)
  rownames(cs) = names(clusters)
  
  for(n in 1:length(clusters))
  {
    # n = 1
    kk = match(dpgp$gene[dpgp$cluster == names(clusters)[n]], rownames(xx))
    cs[n,] = apply(xx[kk, ], 2, mean)
  }
  
  library(corrplot)
  M = cor(t(cs))
  corrplot.mixed(M, order = 'hclust', addrect = 10)
  
  cor_cutoff = 0.6
  dxx = dpgp
  missed = which(is.na(match(rownames(yy), dpgp$gene)))
  losts = c()
  
  for(n in missed)
  {
    cat(n, '\n')
    test = cor(as.numeric(yy[n, ]), as.matrix(t(cs)))
    ii_max = which.max(test)
    cluster.assigned = rownames(cs)[ii_max] 
    if(test[ii_max] >= cor_cutoff){
      dxx = rbind(dxx, c(cluster.assigned, rownames(yy)[n]))
    }else{
      losts = c(losts, n)
    }
  }
  
  cat(length(losts)/length(missed))
  
  saveRDS(dxx, file = paste0(RdataDir, '/DPGP_clusters_extended_allPeaks.rds'))
  
  ##########################################
  # manually merge DPGP clusters to have coarse groups 
  ##########################################
  dpgp = readRDS(file = paste0(RdataDir, '/DPGP_clusters_extended_allPeaks.rds'))
  
  cs = dpgp[!is.na(match(dpgp$cluster, names(clusters))), ]
  cs$mc = cs$cluster
  cs$mc[which(cs$cluster == '5' | 
                cs$cluster == '7'|
                cs$cluster == '8' |
                cs$cluster == '11')] = 'mc1'
  
  cs$mc[which(cs$cluster == '9')] = 'mc2'
  
  cs$mc[which(cs$cluster == '12' | 
                cs$cluster == '19'|
                cs$cluster == '10')] = 'mc3'
  
  
  cs$mc[which(cs$cluster == '4')] = 'mc4'
  cs$mc[which(cs$cluster == '15')] = 'mc5'
  
  cs$mc[which(cs$cluster == '18')] = 'mc6'
  cs$mc[which(cs$cluster == '1')] = 'mc7'
  cs$mc[which(cs$cluster == '13'| cs$cluster == '2')] = 'mc8'
  
  res$clusters = NA
  mm = match(cs$gene, rownames(res))
  
  res$clusters[mm] = cs$mc
  
  saveRDS(res, file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds'))
  
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
# Section : analysis of regeneration genes
# associate regeneration genes with chromatin features: atac, histone markers
# 
########################################################
########################################################
process.normalize.atac.histM.allTSS = function()
{
  ##########################################
  # process the sample info and count table
  ##########################################
  RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
  RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
  source(RNA.functions)
  source(RNA.QC.functions)
  
  ## import sample infos
  design = readRDS(file = paste0('../results/CT_merged_20220328/Rdata/histM_CT_design_info.rds'))
  #design = read.csv(file = paste0(dataDir, "R11876_R12810_R12965_CT_analysis_20220217_QCs_AK.csv"))
  #design = design[!is.na(design$sampleID), ]
  design$fileName = paste0(design$condition, '_', design$sampleID)
  design = design[, c(1:3, 15, 5, 16, 4, 6:14)]
  colnames(design)[5] = 'batch'
  
  design.histM = design[, c(1:5)]
  design.histM = design.histM[grep('rRep', design.histM$batch), ]
  
  design = readRDS(file = paste0('../results/Rxxxx_R10723_R11637_R12810_atac/Rdata', 
                                 '/design_sels_bc_TMM_combat_mUA_regeneration_dev_2Batches.R10723_R7977',
                                          'Rxxxx_R10723_R11637_R12810_atac', '.rds'))
  design = design[, c(1, 2, 6)]
  design$sample = design$condition
  design$marks = 'atac'
  design = design[, c(1, 4,5, 2:3)]
  
  colnames(design.histM) = colnames(design)
  
  design = rbind(design, design.histM)
  design$sample = gsub('BL_UA_13days_distal', 'BL13days.dist', design$sample)
  design$sample = gsub('BL_UA_13days_proximal', 'BL13days.prox', design$sample)
  design$sample = gsub('BL_UA_9days', 'BL9days', design$sample)
  design$sample = gsub('BL_UA_5days', 'BL5days', design$sample)
  design$sample = gsub('Mature_UA', 'mUA', design$sample)
  
  #design = design[grep('Embryo', design$sample, invert = TRUE), ]
  design = design[grep('IgG', design$marks, invert = TRUE), ]
  design$condition = paste0(design$marks, '_', design$sample)
  
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_histM_both/'
  xlist = list.files(path=paste0(dataDir, 'featurecounts_peaks.Q30'), pattern = 'featureCounts.txt.summary', 
                     full.names = TRUE)
  
  design$total.reads = NA
  for(n in 1:nrow(design))
  {
    kk = grep(design$SampleID[n], xlist)
    cat(length(kk), '\n')
    test = read.table(file = xlist[kk], sep = '\t', header = TRUE)
    design$total.reads[n] = sum(as.numeric(as.character(test[,2])))
    
  }
  
  saveRDS(design, file = paste0(RdataDir, '/atac_histM_sample_design_regeneration.rds'))
  
  
  ## process the counts
  design = readRDS(file = paste0(RdataDir, '/atac_histM_sample_design_regeneration.rds'))
    
  xlist<-list.files(path=paste0(dataDir, 'featurecounts_peaks.Q30'),
                    pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge
  
  sels = c()
  for(n in 1:nrow(design)) sels = c(sels, grep(design$SampleID[n], xlist))
  xlist = xlist[sels]
  
  all = cat.countTable(xlist, countsfrom = 'featureCounts')
  
  colnames(design)[1] = 'SampleID'
  
  counts = process.countTable(all=all, design = design[, c(1,4)])
  
  save(design, counts, file = paste0(RdataDir, '/atac_histM_sample_design_regeneration_counts_withEmbryo.Rdata'))
  
  ##########################################
  # clean TSS, for each gene, only keep one TSS if there are mulitple ones
  # the TSS with highest H3K4me3 were chosen
  ##########################################
  load(file = paste0(RdataDir, '/atac_histM_sample_design_regeneration_counts_withEmbryo.Rdata'))
  rownames(counts) = counts$gene
  
  # quick and temporal cpm normalization
  cpm = as.matrix(counts[, -1])
  for(n in 1:ncol(cpm))
  {
    cpm[,n] = (cpm[,n]+0.5)/design$total.reads[n]*10^6
  }
  cpmm = log2(cpm)
  
  ## use geneID 
  annot = readRDS(paste0(annotDir, 
                         'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
  
  mapping = data.frame(transcript = rownames(cpmm), stringsAsFactors = FALSE)
  mapping$gene = annot$geneID[match(rownames(cpmm), annot$transcriptID)]
  mapping$gene[which(is.na(mapping$gene))] = mapping$transcript[which(is.na(mapping$gene))]
  
  test = table(mapping$gene)
  test = test[which(test<=10)]
  
  gg.unique = unique(mapping$gene)
  index_tss = rep(NA, length(gg.unique))
  
  for(n in 1:length(gg.unique))
  {
    # n = 3
    kk = which(mapping$gene == gg.unique[n])
    
    if(length(kk) == 1) {
      index_tss[n] = kk
    }else{
      cat(n, '\n')
      #stop()
      pmarkers = apply(cpmm[kk, grep('H3K4me3', colnames(cpmm))], 1, median)
      index_tss[n] = kk[which.max(pmarkers)]
    }
  }
  
  tss = counts[index_tss, ]
  tss = data.frame(tss[,-1], geneID = gg.unique, transcript = tss[,1], stringsAsFactors = FALSE)
  rownames(tss) = gg.unique
  
  ## add coordinates of tss
  saf = read.table(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_histM_both/',
                'amex6_TSS_all_56670genes.saf'), sep = '\t',  header = TRUE)
  mm = match(tss$transcript, saf$GeneID)
  tss$coords = paste0(saf$Chr[mm], ':', saf$Start[mm], '-', saf$End[mm])
  
  save(tss, design, file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration_counts_withEmbryo.Rdata'))
  
  ##########################################
  # here we will consider the list of gene/tss to keep 
  ##########################################
  load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration_counts_withEmbryo.Rdata'))
  
  #res = readRDS(file = paste0(RdataDir,  '/atac_histoneMarkers_normSignals_axolotlAllTSS.2kb_TSS_genelevel.rds'))
  #res = data.frame(res)
  #res$geneID = rownames(res)
  
  geneClusters = readRDS(file = paste0("../results/RNAseq_data_used/Rdata/", 'regeneration_geneClusters.rds'))
  geneClusters$gene = rownames(geneClusters)
  geneClusters$geneID = sapply(rownames(geneClusters), function(x) {test = unlist(strsplit(as.character(x), '_')); return(test[length(test)])})
  
  mm = match(tss$geneID, geneClusters$geneID)
  tss = data.frame(tss, geneClusters[mm, c(1:2, 26, 3:7)], stringsAsFactors = FALSE)
  
  tss$groups[which(!is.na(tss$groups))] = 'reg_up'
  tss$groups[which(tss$cluster == 'G1'|tss$cluster == 'G2'|tss$cluster == 'G3')] = 'reg_down'
  
  # limb fibroblast expressing genes
  expressed = readRDS(file = '../data/expressedGenes_list_limb_fibroblast_using_smartseq2.mature.regeneration_pooledscRNAseq.dev.rds')
  tss$fibroblast.expressed = expressed$expressed[match(tss$geneID, expressed$geneID)]
  
  tss$groups[which(is.na(tss$groups) & tss$fibroblast.expressed == '1')] = 'expressed.stable'
  
  # house keeping and other non-expressing genes
  load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
  hkgs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
  nonexp = controls.tissue$geneIDs[which(controls.tissue$tissues != 'housekeeping')]
  
  mm = match(tss$geneID, hkgs)
  tss$groups[which(tss$groups == 'expressed.stable' & !is.na(mm))] = 'house_keep'
  
  mm = match(tss$geneID, nonexp)
  tss$groups[which(is.na(tss$groups) & !is.na(mm))] = 'non_expr'
  
  # random select 1000 non-expressed genes
  jj = which(is.na(tss$groups))
  jj = sample(jj, size = 2000, replace = FALSE)
  tss$groups[jj] = 'non_expr'
  
  tss = tss[which(!is.na(tss$groups)), ]
  
  save(tss, design, 
       file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration.embryo_counts_expressedGenes_controls.Rdata'))
  
  #saveRDS(cpm, file = paste0(RdataDir,  '/atac_histoneMarkers_normSignals_axolotlAllTSS.2kb_TSS_genelevel_beforeBatchCorrection.rds'))
  #save(tss, design, 
  #     file = paste0(RdataDir,  '/atac_histoneMarkers_normSignals_axolotlAllTSS.2kb_TSS_genelevel_beforeBatchCorrection.Rdata'))
  
  ##########################################
  # compare the tss and histM 
  ##########################################
  Check.overlaps.between.tss.with.histM = FALSE
  if(Check.overlaps.between.tss.with.histM){
    load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration.embryo_counts_expressedGenes_controls.Rdata'))
    
    tp = data.frame(t(sapply(tss$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    tp$strand = '*'
    tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    load(file = paste0('../results/CT_merged_20220328/Rdata', 
                       '/histoneMarkers_samplesDesign_ddsPeaksMatrix_filtered_incl.atacPeak140k.missedTSS.bgs_345k.Rdata'))
    
    peakNames = rownames(dds)
    peakNames = gsub('tss.', '', peakNames)
    pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    mapping = findOverlaps(tp, pp, ignore.strand=TRUE,  minoverlap=100L)
    jj = (unique(mapping@from))
    missed = setdiff(c(1:nrow(tss)), jj)
    
    ## we won't touch it, because most of tss were already there
  }
  
  Check.overlaps.between.tss.with.atacseq = FALSE
  if(Check.overlaps.between.tss.with.atacseq){
    load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration.embryo_counts_expressedGenes_controls.Rdata'))
    
    tp = data.frame(t(sapply(tss$coords, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    tp$strand = '*'
    tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    load(file = paste0(RdataDir, '/regeneration_samples_beforeBatchCorrection.Rdata'))
    
    pp = data.frame(t(sapply(rownames(ddx), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    
    mapping = findOverlaps(tp, pp, ignore.strand=TRUE,  minoverlap=100L)
    jj = (unique(mapping@from))
    missed = setdiff(c(1:nrow(tss)), jj)
    
    cat(length(missed), ' tss missed \n')
    
    ## add missing tss to regeneration atac data
    load(file = paste0(RdataDir, '/tss_perGene_atac_histM_sample_design_regeneration.embryo_counts_expressedGenes_controls.Rdata'))
    load(file = paste0(RdataDir, '/regeneration_samples_beforeBatchCorrection.Rdata'))
    
    counts = counts(ddx)
    #kk = grep('Embryo_', design.sels$condition, invert = TRUE)
    #design.sels = design.sels[kk, ]
    #counts = counts[, kk]
    
    kk2 = grep('atac', design$marks)
    raw = tss[missed, kk2]
    design.atac = design[kk2, ]
    
    rownames(raw) = tss$coords[missed]
    
    mm = match(design.sels$sampleID, design.atac$SampleID)
    design.atac = design.atac[mm, ]
    raw = raw[ ,mm]
    
    colnames(raw) = colnames(counts)
    counts = rbind(counts, raw)
    
    ## batch correction for those 69k loci 
    table(design.sels$conds, design.sels$batch)
    
    d <- DGEList(counts=counts, group=design.sels$conds)
    tmm <- calcNormFactors(d, method='TMM')
    tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
    
    bc = as.factor(design.sels$batch)
    mod = model.matrix(~ as.factor(conds), data = design.sels)
    
    # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
    # because there is no better batche other others 
    #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
    fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
    
    #design.tokeep<-model.matrix(~ 0 + conds,  data = design.sels)
    #cpm.bc = limma::removeBatchEffect(tmm, batch = bc, design = design.tokeep)
    # plot(fpm.bc[,1], tmm[, 1]);abline(0, 1, lwd = 2.0, col = 'red')
    
    make.pca.plots(tmm, ntop = 3000, conds.plot = 'all')
    ggsave(paste0(resDir, "/regeneration_embryo_TSS_Samples_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    make.pca.plots(fpm.bc, ntop = 3000, conds.plot = 'all')
    ggsave(paste0(resDir, "/regeneration_embryo_TSS_Samples_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    fpm = fpm.bc
    
    rm(fpm.bc)
    
    saveRDS(fpm, file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration.dev.TSS_2Batches.R10723_R7977_', version.analysis, '.rds'))
    saveRDS(design.sels, 
            file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration.dev.TSS_2Batches.R10723_R7977',
                                       version.analysis, '.rds'))
    
    ## dynamic peak analysis
    fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration.dev.TSS_2Batches.R10723_R7977_', version.analysis, '.rds'))
    design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration.dev.TSS_2Batches.R10723_R7977',
                                   version.analysis, '.rds'))
    
    # prepare the background distribution
    fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
    fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
    # rownames(fpm) = gsub('_', '-', rownames(fpm))
    
    
    conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", 'Embryo_Stage44_distal', 
              "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", 'BL_UA_13days_distal')
    
    sample.sels = c(); cc = c()
    sample.means = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n])
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
      sample.means = cbind(sample.means, apply(fpm[, kk], 1, mean))
    }
    colnames(sample.means) = conds
    
    source('Functions_atac.R')
    cpm = fpm[, sample.sels]
    
    source('Functions_atac.R')
    tic()
    ## define the dynamic enhancers with mature UA and BL.UA and check them if embryo samples
    #sels = grep('Embryo', cc, invert = TRUE) 
    res = temporal.peaks.test(cpm, c = cc)
    
    toc()
    
    res = data.frame(fpm, res, stringsAsFactors = FALSE)
    
    saveRDS(res, file = paste0(RdataDir, '/res_temporal_dynamicPeaks_regeneration.dev.TSS_2Batches.R10723_R7977_peakAnnot_v1.rds'))
    
    
  }
  
  ##########################################
  # batch correction of each chromatin features separately for TSS
  # not used here,
  # because all histM includes the majority of TSS and the batch correction is only redone for atac TSS (see above)
  ##########################################
  # Batch.correction.TSS = FALSE
  # if(Batch.correction.TSS){
  #   require(edgeR)
  #   require(sva)
  #   source('Functions_atac.R')
  #   
  #   load(file = paste0(RdataDir,  '/atac_histoneMarkers_normSignals_axolotlAllTSS.2kb_TSS_genelevel_beforeBatchCorrection.Rdata'))
  #   raw = data.frame(raw[,-1], transcript = raw[, 1], stringsAsFactors = FALSE)
  #   
  #   conds = c('atac','H3K4me3', 'H3K4me1', 'H3K27me3', 'H3K27ac')
  #   samples = c('mUA', 'BL5days', 'BL9days', 'BL13days.prox', 'BL13days.dist')
  #   
  #   
  #   tss_list = readRDS(file = paste0(RdataDir, '/list_TSS_considered.rds'))
  #   raw = raw[which(!is.na(match(rownames(raw), tss_list$geneID))), ]
  #   
  #   aa = c()
  #   
  #   for(n in 1:length(conds))
  #   {
  #     # n = 4
  #     sels = c()
  #     for(s in samples)
  #     {
  #       sels = c(sels, which(design$condition == paste0(conds[n], '_', s)))
  #     }
  #     
  #     cat(n, ' --', conds[n], '--', samples, '--', length(sels), 'samples\n')
  #     
  #     counts = raw[, sels]
  #     design.sel = design[sels, ]
  #     
  #     #ss = apply(counts, 1, sum)
  #     #ii_missed = which(ss<10)
  #     
  #     # normalization
  #     d <- DGEList(counts=counts, group=design.sel$condition)
  #     
  #     if(conds[n] == 'H3K27ac' | conds[n] == 'IgG'){
  #       tmm <- calcNormFactors(d, method='TMMwsp') # TMMwsp used for H3K27ac
  #     }else{
  #       tmm <- calcNormFactors(d, method='TMM')
  #     }
  #     
  #     tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0.5)
  #     
  #     ## batch correction
  #     #design.sel$condition = droplevels(design.sel$condition)
  #     #design.sel$batch = droplevels(design.sel$batch)
  #     bc = as.factor(design.sel$batch)
  #     mod = model.matrix(~ as.factor(condition), data = design.sel)
  #     
  #     make.pca.plots(tmm, ntop = 3000, conds.sel = as.character(unique(design.sel$condition)))
  #     ggsave(paste0(resDir, "/", conds[n], "_beforeBatchCorrect_",  version.analysis, ".pdf"), width = 16, height = 14)
  #     
  #     #tmm = as.matrix(tmm)
  #     #vars = apply(tmm, 1, var)
  #     #tmm = tmm[which(vars>0.5), ]
  #     # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
  #     # because there is no better batche other others 
  #     #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)
  #     
  #     fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
  #     
  #     make.pca.plots(fpm.bc, ntop = 3000, conds.sel = as.character(unique(design.sel$condition)))
  #     ggsave(paste0(resDir, "/", conds[n], "_afterBatchCorrect_",  version.analysis, ".pdf"), width = 10, height = 6)
  #     
  #     aa = cbind(aa, fpm.bc)
  #     
  #     rm(fpm.bc)
  #     
  #     #saveRDS(fpm, file = paste0(RdataDir, '/fpm_bc_TMM_combat_', conds[n], '_', version.analysis, '.rds'))
  #     #saveRDS(design.sel, file = paste0(RdataDir, '/design.sels_bc_TMM_combat_', conds[n], '_', version.analysis, '.rds'))
  #     
  #   }
  #   
  #   aa = data.frame(aa, raw$transcript, stringsAsFactors = FALSE)
  #   
  #   ## sample means regardless of batches
  #   #source('Functions_histM.R')
  #   #cpmm = cal_sample_means(cpm, conds = unique(design$condition))
  #   saveRDS(aa, file = paste0(RdataDir,  '/atac_histoneMarkers_normSignals_axolotlAllTSS.2kb_TSS_genelevel_batchCorrected.rds'))
  # }
  # 
}

##########################################
# finally realize that it is better to do batch correction altogether for each chromatin features
# instead of doing it for tss and enhancer
# so we are going to do:
# add missing tss of interest by checking if all tss to considered are already included in the 
# atac and all 4 histone markers
# not done yet
##########################################





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


##########################################
# compare smartseq2 mature with dev  
##########################################
Smartseq2.mature.vs.dev = function()
{
  
  ##########################################
  # similar comparison mUA vs stage samples for mRNAs 
  ##########################################
  RNAseqDir = "../results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/"
  load(file = paste0(RNAseqDir, 'RNAseq_design_dds.object.Rdata'))
  
  # select mature samples
  sels = grep('Mature_UA|Embryo', design.matrix$condition)
  dds = dds[, sels]
  
  dds$condition = droplevels(dds$condition)
  
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  pca=plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = FALSE)
  print(pca)
  
  pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE))
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 0.7, nudge_y = 1, size=2.5)
  
  plot(ggp)
  
  
  dds$condition = droplevels(dds$condition)
  dds <- estimateDispersions(dds, fitType = 'parametric')
  plotDispEsts(dds, ymin = 10^-3); abline(h = 0.1, col = 'blue', lwd = 2.0)
  # 
  dds = nbinomWaldTest(dds, betaPrior = TRUE)
  resultsNames(dds)
  
  # 
  res.ii = results(dds, contrast=c("condition", 'Embryo_Stage40', 'Mature_UA'), alpha = 0.1)
  colnames(res.ii) = paste0(colnames(res.ii), "_Hand.vs.UA")
  res = data.frame(res.ii[, c(2, 5, 6)])
  # 
  res.ii = results(dds, contrast=c("condition", 'Mature_Hand', 'Mature_LA'), alpha = 0.1)
  colnames(res.ii) = paste0(colnames(res.ii), "_Hand.vs.LA")
  res = data.frame(res, res.ii[, c(2, 5, 6)])
  
  
  cpm = fpm(dds)
  cpm = log2(cpm + 2^-7)
  
  xx = data.frame(mUA = apply(cpm[, grep('Mature_UA', colnames(cpm))], 1, median), 
                  stage40 = apply(cpm[, grep('Embryo_Stage40', colnames(cpm))], 1, mean), 
                  stage46_prox = apply(cpm[, grep('Embryo_Stage46_proximal', colnames(cpm))], 1, mean), 
                  stage46_dist = apply(cpm[, grep('Embryo_Stage46_distal', colnames(cpm))], 1, mean))
  
  maxs = apply(xx, 1, max)
  
  cpm = xx
  cpm$gene =  sapply(rownames(cpm), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[1])})
  
  jj = grep(paste0(dev.example, collapse = '|') ,cpm$gene)
  cpm = cpm[which(maxs> -1), ]
  
  load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
  hs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
  ctl =  controls.tissue$geneIDs[which(controls.tissue$tissues  != 'housekeeping')]
  ggs = sapply(rownames(cpm), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})
  
  cpm$genetype = NA
  
  mm = match(ggs, hs)
  cpm$genetype[!is.na(mm)] = 'hs'
  
  mm = match(ggs, ctl)
  cpm$genetype[!is.na(mm) & is.na(cpm$genetype)] = 'ctl'
  
  par(mfrow = c(1, 3))
  hist(cpm[, 2], breaks = 100, main = 'log2 expression of stage40');abline(v = 2)
  hist(cpm[, 3], breaks = 100, main = 'log2 expression of stage44_prox');abline(v = 2)
  hist(cpm[, 4], breaks = 100, main = 'log2 expression of stage44_dist');abline(v = 2)
  
  cpm$genetype[which(cpm[,2]< -1 & cpm[, 3]< -1 & cpm[, 4] < -1 & is.na(cpm$genetype))] = 'dev.lowlyExp'
  cpm$genetype[is.na(cpm$genetype)] = 'dev.Exp'
  
  cpm$genetype[!is.na(match(cpm$gene, limb.go))] = 'limb.dev.GO'
  
  examples.sel = unique(grep(paste0(dev.example, collapse = '|') ,cpm$gene))
  
  cpm = data.frame(cpm)
  cpm = cpm[which(cpm$genetype != 'ctl'), ]
  
  ggplot(cpm, aes(x = stage40, y = mUA, color = genetype, label = gene)) +
    geom_point(size = 0.4) + 
    scale_color_manual(values=c("#E69F00", "#56B4E9",  "#999999",  'darkblue')) + 
    geom_text_repel(data= cpm[examples.sel, ], size = 4.0, color = 'darkblue') +
    geom_point(data=cpm[which(cpm$genetype == 'limb.dev'), ], aes(x=stage40, y=mUA), colour="darkblue", size=1.5) +
    #geom_hline(yintercept=2.0, colour = "darkgray") + 
    #geom_vline(xintercept = 2.0, colour = "darkgray")
    geom_abline(slope = 1,  intercept = 0, colour = 'red') +
    theme_classic() +
    theme(legend.text = element_text(size=10), 
          #legend.title = element_text(color = "blue", size = 14),
          legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    )
  
  
  ggsave(paste0(figureDir, "Stage40_vs_mUA_devGenes_comparisons_smartseq2.pdf"), width=12, height = 10)
  
}

##########################################
# segment-specific chromatin 
# and subclustering atac-peak clusters
##########################################
aggregate_atacPeaks_histMpeaks = function()
{
  ### collect the four markers with means 
  histMs =  c('H3K4me3','H3K27me3', 'H3K4me1')
  shm = c()
  
  for(n in 1:length(histMs))
  {
    # n = 1
    cat(n, ' -- ', histMs[n], '\n')
    cpm = readRDS(file = paste0(RdataHistM, '/fpm_bc_TMM_combat_', histMs[n], '_', version.analysis.HistM, '.rds'))
    cpm = cpm[ , grep(histMs[n], colnames(cpm))]
    #design.sel = readRDS(file = paste0(RdataHistM, '/design.sels_bc_TMM_combat_DBedgeRtest_', histMs[n], '_', 
    #                                   version.analysis.HistM, '.rds'))
    
    ### select the samples and extract sample means
    conds_histM = c("mUA", "mLA", "mHand")
    sample.sels = c();  
    cc = c()
    
    sample.means = c()
    for(ii in 1:length(conds_histM)) 
    {
      kk = grep(conds_histM[ii], colnames(cpm))
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds_histM[ii], length(kk)))
      if(length(kk)>1) {
        sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
      }else{
        sample.means = cbind(sample.means, cpm[, kk])
      }
    }
    colnames(sample.means) = paste0(conds_histM, '_', histMs[n])
    
    if(n == 1){
      shm = sample.means
    }else{
      shm = cbind(shm, sample.means[match(rownames(shm), rownames(sample.means)), ])
    }
  }    
  
  ## quantile normalization for different histone markers
  library(preprocessCore)
  xx = normalize.quantiles(shm)
  colnames(xx) = colnames(shm)
  rownames(xx) = rownames(shm)
  shm = xx
  
  pp_histM = data.frame(t(sapply(rownames(shm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp_histM$strand = '*'
  pp_histM = makeGRangesFromDataFrame(pp_histM, seqnames.field=c("X1"),
                                      start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(pp, pp_histM, type = 'within')
  
  yy0 = shm[mapping@to, ]
  #rownames(yy0) = rownames(yy)
  
  # reorder the histM according to the atac-seq peak clusters
  yy0 = yy0[plt$tree_row$order, ]
  yy0 = t(apply(yy0, 1, cal_z_score))
  #test = yy[plt$tree_row$order, ]
  #rownames(yy0) = rownames(yy)
  #yy0 = cbind(yy, yy0)
  
  df_histM = data.frame(segments = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[1])
                        #markers = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[2]), stringsAsFactors = FALSE
  )
  colnames(df_histM) = c('seg')
  rownames(df_histM) = colnames(yy0)
  sample_colors_histM = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors_histM) = c('mUA', 'mLA', 'mHand')
  marker_colors_histM = c('blue', 'red', 'deepskyblue2', 'darkgreen')
  names(marker_colors_histM) = histMs
  annot_colors_histM = list(seg = sample_colors_histM)
  gaps.col_histM = seq(1, ncol(yy0), by = 1)
  #clusters <- rownames(subsetDat[heat$tree_row[["order"]],])
  clusters <- sort(cutree(plt$tree_row, k = nb_clusters))
  
  # cluster order : 6, 1, 5, 3, 4, 2
  gaps.row = c(32, 32+76, 31 + 32 + 76, 292 + 31 + 32 + 76, 292 + 31 + 32 + 76 + 103)
  
  #yy0 = yy0[plt$tree_row$order, ]
  p1 = pheatmap(test, 
                nnotation_row = my_gene_col, 
                annotation_col = df, show_rownames = FALSE, scale = 'none', 
                color = col, 
                show_colnames = FALSE,
                cluster_rows = FALSE, cluster_cols = FALSE,  
                #clustering_method = 'complete', cutree_rows = nb_clusters, 
                annotation_colors = annot_colors, 
                #clustering_callback = callback,
                legend = TRUE,
                gaps_col = gaps.col, 
                gaps_row = gaps.row)
  
  p2 = pheatmap(yy0, cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  #grid.arrange(p1, p2,  nrow = 1)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  #plot_list[['p3']]=p2[[4]]
  grid.arrange(grobs=plot_list, ncol=2)
  
  indexs = data.frame(peak = names(clusters), cluster= clusters, stringsAsFactors = FALSE)
  indexs = indexs[match(rownames(test), rownames(indexs)), ]
  c_order = unique(indexs$cluster)
  
  yy1 = yy0
  
  ### transform the histM to account for different backgrounds and scaling
  # for(m in 1:ncol(yy1))
  # {
  #   jj1 = which(yy1[ ,m] < 1.5)
  #   yy1[jj1, m] = 1.5
  #   jj2 = which(yy1[ ,m] > 4)
  #   yy1[jj2, m] = 4
  #    
  # }
  # 
  
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  
  for(ac in c_order)
  {
    # ac = 6
    kk = which(indexs$cluster == ac)
    cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
    
    hm_hclust <- hclust(dist(as.matrix(yy1[kk,])), method = "complete")
    #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
    peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
    
  }
  
  new_order = match(peakNm, rownames(yy1))
  
  p1 = pheatmap(test[new_order, ], 
                nnotation_row = my_gene_col, 
                annotation_col = df, show_rownames = FALSE, scale = 'none', 
                color = col, 
                show_colnames = FALSE,
                cluster_rows = FALSE, cluster_cols = FALSE,  
                #clustering_method = 'complete', cutree_rows = nb_clusters, 
                annotation_colors = annot_colors, 
                #clustering_callback = callback,
                legend = TRUE,
                annotation_legend = FALSE,
                gaps_col = gaps.col, 
                gaps_row = gaps.row)
  
  gaps.col_histM = c(1:2)
  df_histM_new = as.data.frame(df_histM[c(1:3),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(1:3)]
  
  p2 = pheatmap(yy1[new_order, c(1:3)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col= df_histM_new,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  df_histM_new = as.data.frame(df_histM[c(4:6),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(4:6)]
  p3 = pheatmap(yy1[new_order, c(4:6)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  df_histM_new = as.data.frame(df_histM[c(7:9),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(7:9)]
  p4 = pheatmap(yy1[new_order, c(7:9)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                annotation_legend = FALSE,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  plot_list[['p3']]=p3[[4]]
  plot_list[['p4']]=p4[[4]]
  
  pdf(paste0(figureDir, "/Segemet_specific_chromatin_landscape_atac_histM_firstTry.pdf"),
      width = 10, height = 12) # Open a new pdf file
  
  layout = matrix(c(1, 1, 2, 3, 4), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  dev.off()
  
}

##########################################
# combine atac-seq peaks and histone marker to make heatmap 
# here histone peaks were overlapping atac-seq peaks or the atac-seq peak regions
# the initial clusters were done based on the atac-seq peaks 
# subclusters were done with histone makrers
##########################################
Assembly_histMarkers_togetherWith_ATACseq_regeneration = function()
{
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  require(pheatmap)
  require(RColorBrewer)
  
  # import atac-seq peak
  # res = readRDS(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration.rds')) # stat of dynamic peaks
  # z-score of data and heatmap (variable: plt, yy and res)
  load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters.Rdata')) 
  #design_atac = design
  res_atac = res
  yy_atac = yy
  
  Test_cluster_order = FALSE
  if(Test_cluster_order){
    test = yy
    test = test[plt$tree_row$order, ]
    pheatmap(test, 
             nnotation_row = my_gene_col, 
             annotation_col = df, show_rownames = FALSE, scale = 'none', 
             color = col, 
             show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE,  
             #clustering_method = 'complete', cutree_rows = nb_clusters, 
             annotation_colors = annot_colors, 
             #clustering_callback = callback,
             gaps_col = ii.gaps)
    
    ## save group bed of each cluster for deeptools heatmap
    load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters.Rdata')) 
    
    res = res[match(rownames(yy), rownames(res)), ]
    mcs = unique(res$clusters)
    for(n in 1:length(mcs)){
      jj = which(res$clusters == mcs[n])
      pp_atac = data.frame(t(sapply(rownames(yy)[jj], function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
      pp_atac$name = rownames(yy)[jj]
      pp_atac$score = 0
      pp_atac$strand = '*'
      
      peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                            'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_groups')
      write.table(pp_atac, file = paste0(peakGroupDir, '/peak_group_', n, '.bed'), 
                  sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    }
    
    ## save bed files of promoters, enhancers of up- and down-regulated atac-seq peaks
    res = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds'))
    res = res[which(!is.na(res$clusters)), ]
    
    pp = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    pp.annots = as.data.frame(pp.annots)
    
    res$atacseq = 'up'
    res$atacseq[which(res$clusters == 'mc1'| res$clusters == 'mc6')] = 'down'
    
    res$annots = NA
    res$annots[grep('Promoter', pp.annots$annotation)] = 'promoter'
    res$annots[grep('Intergenic|Intron', pp.annots$annotation)] = 'enhancer'
    
    pp_atac = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    pp_atac$name = rownames(res)
    pp_atac$score = 0
    pp_atac$strand = '*'
    
    peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                          'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_promoter_enhancer')
    
    jj = which(res$atacseq == 'up' & res$annots == 'enhancer')
    write.table(pp_atac[jj, ], file = paste0(peakGroupDir, '/peak_enhancer_up.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    jj = which(res$atacseq == 'down' & res$annots == 'enhancer')
    write.table(pp_atac[jj, ], file = paste0(peakGroupDir, '/peak_enhancer_down.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    
    jj = which(res$atacseq == 'up' & res$annots == 'promoter')
    xx = pp_atac[jj, ]
    xx0 =  pp.annots[jj, ]
    xx$strand = xx0$geneStrand # keep tss strand
    kk = which(xx$strand == '1')
    xx$X2 = as.numeric(as.character(xx$X2))
    xx$X3 = as.numeric(as.character(xx$X3))
    xx$X2[kk] = xx0$geneStart[kk]; xx$X2[-kk] = xx0$geneStart[-kk] - 1 
    xx$X3[kk] = xx$X2[kk] + 1; xx$X3[-kk] = xx$X2[-kk]
    xx$strand[kk] = '+'
    xx$strand[-kk] = '-'
    
    write.table(xx, file = paste0(peakGroupDir, '/peak_tss_up.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    jj = which(res$atacseq == 'down' & res$annots == 'promoter')
    xx = pp_atac[jj, ]
    xx0 =  pp.annots[jj, ]
    xx$strand = xx0$geneStrand # keep tss strand
    kk = which(xx$strand == '1')
    xx$X2 = as.numeric(as.character(xx$X2))
    xx$X3 = as.numeric(as.character(xx$X3))
    xx$X2[kk] = xx0$geneStart[kk]; xx$X2[-kk] = xx0$geneStart[-kk] - 1 
    xx$X3[kk] = xx$X2[kk] + 1; xx$X3[-kk] = xx$X2[-kk]
    xx$strand[kk] = '+'
    xx$strand[-kk] = '-'
    
    
    write.table(xx, file = paste0(peakGroupDir, '/peak_tss_down.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    ## save promoters to test histone markers
    peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                          'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_tss')
    
    tss = promoters(amex, upstream=2000, downstream=2000, use.names=TRUE)
    ggs = genes(amex)
    ggs = as.data.frame(ggs)
    
    tss = ggs
    jj = which(tss$strand == '-')
    tss$start[jj] = tss$end[jj]
    tss$end = tss$start + 1
    #tss$strand = '*'
    tss$score = 0
    tss = tss[, c(1, 2, 3, 6, 7, 5)] 
    
    write.table(tss, file = paste0(peakGroupDir, '/tss_expressed.genes.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  # import histone marker 
  RdataHistM = '/Users/jiwang/workspace/imp/positional_memory/results/CT_merged_20220328/Rdata'
  load(file = paste0(RdataHistM, '/combined_4histMarkers_overlapped55kATACseq_DE_regeneration.Rdata')) # variables (keep and DE.locus)
  design = readRDS(file = paste0(RdataHistM, '/histM_CT_design_info.rds'))
  
  yy = keep[, grep('^H3K4me3|^H3K27me3|^H3K4me1|^H3K27ac', colnames(keep))]
  
  sampleID = sapply(colnames(yy), function(x) unlist(strsplit(as.character(x), '_'))[3])
  mm = match(sampleID, design$sampleID)
  colnames(yy) = paste0(design$condition[mm], '_', design$Batch[mm], '_', design$sampleID[mm])
  
  source('Functions_histM.R')
  conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
  conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')
  cc = paste(rep(conds_histM, each = length(conds)), conds, sep = "_")
  yy = cal_sample_means(yy, conds = cc)
  
  ## select the regions overlapping with atac-seq
  pp_atac = data.frame(t(sapply(rownames(yy_atac), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  pp_atac$strand = '*'
  pp_atac = makeGRangesFromDataFrame(pp_atac, seqnames.field=c("X1"),
                                     start.field="X2", end.field="X3", strand.field="strand")
  
  pp_histM = data.frame(t(sapply(rownames(yy), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp_histM$strand = '*'
  pp_histM = makeGRangesFromDataFrame(pp_histM, seqnames.field=c("X1"),
                                      start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(pp_atac, pp_histM)
  
  ## peaks mapped to postional atac-peaks 
  yy1 = yy[mapping@to, ]
  yy_atac = yy_atac[mapping@from,]
  res_atac = res_atac[mapping@from, ]
  
  ## peaks not mapped to positional atac-peaks
  yy0 = yy[-mapping@to, ]
  
  ss = apply(DE.locus[, c(1:3)], 1, sum)
  mm = match(names(ss[which(ss>0)]), rownames(yy0))
  length(which(!is.na(mm)))
  
  yy0 = yy0[mm[!is.na(mm)], ]
  
  #plot.pair.comparison.plot(yy1[, grep('H3K4me1_mUA_', colnames(yy1))], linear.scale = FALSE)
  #plot.pair.comparison.plot(yy1[, grep('H3K4me3_mUA_', colnames(yy1))], linear.scale = FALSE)
  #plot.pair.comparison.plot(yy1[, grep('H3K27me3_mUA_', colnames(yy1))], linear.scale = FALSE)
  saveRDS(yy1, file = paste0(RdataDir, '/peak_signals_atac_4histM_regenerationPeaks.rds'))
  saveRDS(yy0, file = paste0(RdataDir, '/peak_signals_atac_4histM_notOverlapped.regenerationPeaks.rds'))
  
  #yy1 = yy1[, grep('rRep', colnames(yy1))]
  #yy0 = yy0[, grep('rRep', colnames(yy0))]
  yy1 = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_regenerationPeaks.rds'))
  yy0 = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_notOverlapped.regenerationPeaks.rds'))
  
  # not consider H3K27ac in the main figure
  # yy1 = yy1[, grep('H3K27ac', colnames(yy1), invert = TRUE)]
  
  mm = match(rownames(yy1), rownames(DE.locus))
  DE.peaks = DE.locus[mm, ]
  ss = apply(DE.peaks, 1, sum)
  length(which(ss>0))
  
  source('Functions_histM.R')
  conds_histM = c('H3K4me1', 'H3K27me3', 'H3K4me3', 'H3K27ac')
  stable.multiplier = 2
  for(n in 1:length(conds_histM))
  {
    # n = 1
    jj = grep(conds_histM[n], colnames(yy1))
    
    # yy1[,jj] = t(apply(yy1[,jj], 1, cal_z_score))
    # dynamic.list = DE.locus[, which(colnames(DE.locus) == conds_histM[n])]
    # names(dynamic.list) = rownames(DE.locus)
    # dynamic.list = dynamic.list[which(dynamic.list>0)]
    # kk = which(is.na(match(rownames(yy1), names(dynamic.list))))
    # cat(length(kk), '-', conds_histM[n], ' stable \n')
    # yy1[kk, jj] = yy1[kk, jj] /stable.multiplier
    # yy1[,jj] = t(apply(yy1[,jj], 1, cal_centering))
    yy1[ ,jj] = t(apply(yy1[,jj], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
    
    jj0 = grep(conds_histM[n], colnames(yy0))
    # yy0[,jj] = t(apply(yy0[,jj], 1, cal_z_score))
    # kk0 = which(is.na(match(rownames(yy0), names(dynamic.list))))
    # #cat(length(kk), '-', conds_histM[n], ' stable \n')
    # yy0[kk0, jj] = yy0[kk0, jj] /stable.multiplier
    yy0[ ,jj0] = t(apply(yy0[,jj0], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
    
  }
  
  save(yy1, res_atac, yy_atac, file = paste0(RdataDir, '/atac_histM_data_regeneration.Rdata'))
  
  ## peaks overlapped with atac-seq 
  df = as.data.frame(sapply(colnames(yy1), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'time'
  rownames(df) = colnames(yy1)
  
  sample_colors = c('springgreen', 'springgreen2', 'springgreen3', 'gold2', 'red')
  annot_colors = list(time = sample_colors)
  
  gaps_col = c(5, 10, 15)
  pheatmap(yy1, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           gaps_col = gaps_col,
           legend = TRUE,
           #annotation_colors = annot_colors,
           width = 5, height = 10, 
           filename = paste0(resDir, '/heatmap_histoneMarker_16k.regenerationPeaks.pdf'))
  
  ## peaks not overlapped with atac
  df = as.data.frame(sapply(colnames(yy0), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'time'
  rownames(df) = colnames(yy1)
  sample_colors = c('springgreen', 'springgreen2', 'springgreen3', 'gold2', 'red')
  annot_colors = list(time = sample_colors)
  
  pheatmap(yy0, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           gaps_col = gaps_col,
           legend = TRUE,
           #annotation_colors = annot_colors,
           width = 5, height = 10, 
           filename = paste0(figureDir, '/heatmap_histoneMarker_DE_notoverlapped.with.atac.regenerationPeaks.pdf'))
  
  
  ##########################################
  #  # reorder the histM according to the atac-seq peak clusters
  # cluster order : 6, 1, 5, 3, 4, 2
  # gaps.row = c(32, 32+76, 31 + 32 + 76, 292 + 31 + 32 + 76, 292 + 31 + 32 + 76 + 103)
  ##########################################
  load(file = paste0(RdataDir, '/atac_histM_data_regeneration.Rdata'))
  peaks = res_atac
  table(peaks$clusters)
  cluster_order = paste0('mc', c(1:8))
  
  conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
  conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')
  
  ## specify row gaps as atac-seq peaks
  gaps.row = c()
  for(n in 1:(length(cluster_order)-1))
  {
    if(n == 1)  {gaps.row = c(gaps.row, length(which(peaks$clusters == cluster_order[n])))
    }else{
      gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(peaks$clusters == cluster_order[n])))
    }
  }
  
  ## refine histone subclusters for each ATACseq peak cluster
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  
  for(ac in cluster_order)
  {
    # ac = 6
    kk = which(peaks$cluster == ac)
    cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
    
    hm_hclust <- hclust(dist(as.matrix(yy1[kk,])), method = "complete")
    #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
    peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
  }
  
  new_order = match(peakNm, rownames(yy1))
  
  ## plot the histone marker side by side with replicates
  df_histM = data.frame(time = sapply(colnames(yy1), function(x) unlist(strsplit(as.character(x), '_'))[2])
                        #markers = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[2]), stringsAsFactors = FALSE
  )
  colnames(df_histM) = c('time')
  rownames(df_histM) = colnames(yy1)
  sample_colors_histM =  sample_colors = c('springgreen', 'springgreen2', 'springgreen3', 'gold2', 'red')
  names(sample_colors_histM) = conds
  #marker_colors_histM = c('blue', 'red', 'deepskyblue2', 'darkgreen')
  #names(marker_colors_histM) = histMs
  annot_colors_histM = list(time = sample_colors_histM)
  gaps.col_histM = c(1:4)
  
  kk = c(1:5)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = colorRampPalette(((brewer.pal(n = 7, name ="BrBG"))))(10)
  p1 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                #color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlGnBu")))(8), 
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col= df_histM_new,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  kk = c(6:10)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da"))
  p2 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  
  kk = c(11:15)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = rev(terrain.colors(10))
  p3 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                annotation_legend = FALSE,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  plot_list[['p3']]=p3[[4]]
  
  pdf(paste0(figureDir, "/regeneration_chromatin_landscape_atac_histM.pdf"),
      width = 6, height = 10) # Open a new pdf file
  
  layout = matrix(c(1, 2, 3), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  dev.off()
  
  ### plot the atac-seq peaks with the changed order
  Replot.atac.peaks.with.newOrder = FALSE
  if(Replot.atac.peaks.with.newOrder){
    conds.atac = c("Embryo_Stage40", "Embryo_Stage44_proximal", 'Embryo_Stage44_distal', 
                   "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", 'BL_UA_13days_distal'
    )
    df <- data.frame(conds.atac)
    rownames(df) = colnames(yy_atac)
    colnames(df) = 'time'
    
    sample_colors = c('magenta', 'darkblue', 'springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2',
                      'red')[c(1:length(conds.atac))]
    names(sample_colors) = conds.atac
    
    # col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
    #           "#33a02c", "#fb9a99", "#e31a1c",
    #           "#fdbf6f", "#ff7f00", "#cab2d6",
    #           "#6a3d9a", "#ffff99", "#b15928")
    # cluster_col = col3[1:nb_clusters]
    # names(cluster_col) = paste0('cluster_', c(1:nb_clusters))
    annot_colors = list(
      condition = sample_colors)
    
    ii.gaps = c(3, 4)
    col = colorRampPalette(c("navy", "white", "red3"))(16)
    
    pheatmap(yy_atac, 
             #annotation_row = my_gene_col, 
             annotation_col = df, show_rownames = FALSE, scale = 'none', 
             color = col, 
             show_colnames = FALSE,
             cluster_rows = FALSE, cluster_cols = FALSE,  
             #clustering_method = 'complete', cutree_rows = nb_clusters, 
             annotation_colors = annot_colors,
             gaps_row = gaps.row,
             #clustering_callback = callback,
             gaps_col = ii.gaps, 
             filename = paste0(figureDir, 'heatmap_regenerationPeaks_scaled.pdf'), 
             width = 6, height = 12)
  }
}

########################################################
########################################################
# Section : analyze TSS of regeneration response genes
# 
########################################################
########################################################
Update.TSS.for.regenerationData = function()
{
  ## first update: modifying gene names; define the lowly expressed genes
  tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM.rds'))
  tss$gene[which(tss$geneID == 'AMEX60DD018448'|tss$geneID == 'AMEX60DD018449')] = NA
  tss$gene[which(tss$geneID == 'AMEX60DD002780')] = NA
  tss$gene[which(tss$geneID == 'AMEX60DD003752'|tss$geneID=='AMEX60DD024053')] = NA
  
  tss$groups1 = tss$groups
  tss$groups = NA
  
  ## import the smart-seq2 data to further refine the gene groups
  load(file=paste0('../results/RNAseq_data_used/Rdata/', # design and dds for smartseq2 REGENERATION data
                   'design_dds_all_regeneration_12selectedSamples_v47.hox.patch.Rdata')) 
  rm(design); 
  dds$condition = droplevels(dds$condition)
  ss = rowSums(counts(dds))
  ids = get_geneID(rownames(dds))
  
  # double check the non-expressed genes
  ggs = tss$geneID[which(tss$groups == 'non_expr')]
  length(intersect(ggs, ids[which(ss==0)]))
  #tss$groups[which(tss$groups == 'non_expr')] = NA
  tss$groups[which(!is.na(match(tss$geneID, ids[which(ss==0)])))] = 'non_expr'
  
  tss$groups[which(!is.na(match(tss$geneID, ids[which(ss>0 & ss<100)])))] = 'lowlyExpr_stable' # by read counts
  tss$groups[which(!is.na(match(tss$geneID, ids[which(ss>100)])))] = 'highlyExpr_stable'
  
  ss = apply(tss[, grep('smartseq2_', colnames(tss))], 1, mean) # define house-keeping genes
  tss$groups[which(tss$groups == 'highlyExpr_stable' & ss>2.5 & tss$groups1 == 'house_keep')] = 'house_keep'
  
  jj = which(tss$groups1 == 'reg_down' & (tss$groups == "lowlyExpr_stable"|tss$groups == 'highlyExpr_stable'|tss$groups == 'house_keep'))
  tss$groups[jj] = 'DE_down'
  
  jj = which(tss$groups1 == 'reg_up' & (tss$groups == "lowlyExpr_stable"|tss$groups == 'highlyExpr_stable'|tss$groups == 'house_keep'))
  tss$groups[jj] = 'DE_up'
  
  
  saveRDS(tss, file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM_geneCorrection_v2.rds'))
  
  ## second update: refine lowly expressed genes and DE_up and DE_down using log2FC
  tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM_geneCorrection_v2.rds'))
  
  tss$groups2 = tss$groups # save the old groups
  
  ss = apply(tss[, grep('smartseq2_', colnames(tss))], 1, mean)
  jj = which(tss$groups == 'highlyExpr_stable')
  tss$groups[jj[which(ss[jj]<0)]] = 'lowlyExpr_stable'
  
  jj = which(tss$groups == 'non_expr')
  
  
  saveRDS(tss, file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM_geneCorrection_v3.rds'))
  
  Test.redefine.Up.down.groups = FALSE # it turn out to be not easy !!
  if(Test.redefine.Up.down.groups){
    rna = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/', 
                                'smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked',
                                '_RNAseq_data_used_20220408', '.rds'))
    rna$geneID = get_geneID(rownames(rna))
    rna = rna[match(tss$geneID, rna$geneID), ]
    rna = rna[, grep('log2FoldChange_dpa', colnames(rna))]
    
    jj = which(tss$groups == 'DE_up')
    tss$time_max = NA
    for(j in jj)
    {
      # j = 70
      test = max(rna[j,])
      if(test > 0.5){
        
        tt = colnames(rna)[which.max(rna[j, ])]
        tt = gsub('log2FoldChange_', '', tt)
        tt = gsub('.vs.mUA', '', tt)
        tss$time_max[j] = tt
        
      }else{
        cat('Error --', j, '\n')
        stop()
      }
    }
  }
  
}

##########################################
# add features for TSS to predict gene groups defined by RNA-seq data 
##########################################
prepare_for_promoterScanning = function()
{
  annot = data.frame(annot)
  
  yy$CpG = NA
  
  outDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/motif_analysis/FIMO_promoters/'
  
  ## file for TF motif analysis
  tp = data.frame(t(sapply(yy$promoters, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  tp$X3 = as.numeric(as.character(tp$X3)) - 1700
  rownames(tp) = rownames(yy)
  tp$strand = '*'
  tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  export(tp, format = 'bed', con = paste0(outDir, 'promoters_2kb.300bp_TFmotifs.bed'))
  
  ## file for TATA box searching
  tp = data.frame(t(sapply(yy$promoters, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  
  #tp$X2 = as.numeric(as.character(tp$X2)) + 1900
  tp$X3 = as.numeric(as.character(tp$X3)) - 2000 - 20
  tp$X2 = tp$X3 - 30
  rownames(tp) = rownames(yy)
  tp$strand = annot$strand[match(rownames(tp), annot$gene_id)]
  tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  export(tp, format = 'bed', con = paste0(outDir, 'promoters_50.20bp_TATAbox.bed'))
  
  ## file for CpG enrichment
  tp = data.frame(t(sapply(yy$promoters, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  tp$X3 = as.numeric(as.character(tp$X3)) - 2000 + 100
  tp$X2 = tp$X3 - 500
  
  rownames(tp) = rownames(yy)
  tp$strand = annot$strand[match(rownames(tp), annot$gene_id)]
  tp = makeGRangesFromDataFrame(tp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  export(tp, format = 'bed', con = paste0(outDir, 'promoters_400.100bp_CpG.bed'))
  
}

add_CpG_features = function(yy)
{
  library('rtracklayer')
  library(GenomicRanges)
  library('GenomicFeatures')
  require(Biostrings)
  annot = import(paste0(annotDir, 'AmexT_v47_Hox.patch.gtf'))
  
  annot = data.frame(annot)
  
  yy$CpG = NA
  
  seqs = read.table(paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/motif_analysis/FIMO_promoters/seq_CpG', 
                           '/seqs_CpG.txt'), sep = '\t', header = FALSE)
  seqs$geneID = gsub("\\s*\\([^\\)]+\\)","",as.character(seqs$V1))
  seqs = seqs[match(rownames(yy), seqs$geneID),]
  
  #test = DNAString(as.character(seqs$V2[1]))
  #countPattern("CG", test)
  test = DNAStringSet(seqs$V2)
  
  ## gc content and CpG observation/expectation (https://www.biostars.org/p/355876/)
  ## example code from https://web.stanford.edu/class/bios221/labs/biostrings/lab_1_biostrings.html
  gc = (letterFrequency(test, letters="G", OR=0) + letterFrequency(test, letters="C", OR=0))/width(test)
  oe = vcountPattern("CG", test)*width(test)/(letterFrequency(test, letters="G", OR=0)*letterFrequency(test, letters="C", OR=0))
  
  yy$gc_percent = gc
  yy$cpg_oe = oe
  
  yy = data.frame(yy, stringsAsFactors = FALSE)
  
  saveRDS(yy, file = paste0(RdataDir, '/chromatin_promoter_features_geneGroups.rds'))
  
  # similar example from 
  # https://github.com/mbh038/HDDA/blob/master/7x/DNA%20Methylation/1.%20Introduction/Introduction%20to%20DNA%20methylation.Rmd
  # res = alphabetFrequency(cgiseq)
  # L = rowSums(res)
  # cprop = res[,"C"]/L
  # gprop = res[,"G"]/L
  # expected=L*cprop*gprop
  # observed=vcountPattern("CG",cgiseq)
  # cpgoe=observed/expected
  
  # ggplot(data = yy, aes(x = groups, y = gc.percent, fill = groups)) +
  #   #geom_point(size = 0.1) +
  #   geom_boxplot(outlier.alpha = 0.1)
  # 
  ggplot(data = yy, aes(x = groups, y = cpg.oe, fill = groups)) +
    #geom_point(size = 0.1) +
    geom_boxplot(outlier.alpha = 0.1)
  
  #tp = data.frame(t(sapply(yy$promoters, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  
  # cpg = read.table(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/CpG/cpg.bed6'),
  #                   sep = '\t', header = FALSE)
  # cpg$strand = '*'
  # cpg = makeGRangesFromDataFrame(cpg, seqnames.field=c("V1"),
  #                                start.field="V2", end.field="V3", strand.field="strand")
  # 
  # # export(cpg, format = 'bed', con = paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/CpG/cpg_island.bed'))
  # mapping = findOverlaps(tp, cpg, ignore.strand=TRUE)
  test.Make.BSgenome = FALSE
  if(test.Make.BSgenome){
    cat('give up\n')
  }
   
}

add_TFmotifs_feature = function()
{
  ## process the fimo output to make a motif occurrency matrix
  Process.fimo.output = FALSE
  if(Process.fimo.output){
    library(data.table)
    
    fimo.out = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/motif_analysis/FIMO_promoters/promoter_tfs/'
    fimo = fread(paste0(fimo.out, 'fimo_out/fimo.tsv'), header = TRUE)
    
    motif.oc = table(fimo$motif_id, fimo$sequence_name, useNA = 'ifany')
    motif.oc = t(motif.oc)
    
    print(head(rownames(motif.oc), 20))
    
    saveRDS(motif.oc, file = '../results/motif_analysis/motif_oc_fimo_jaspar2022_pval.0.0001_regenerationTSS.rds')
  }
  
  motif.oc = readRDS(file = '../results/motif_analysis/motif_oc_fimo_jaspar2022_pval.0.0001_regenerationTSS.rds')
  yy = readRDS(file = paste0(RdataDir, '/chromatin_promoter_features_geneGroups.rds')) 
  
  tp = data.frame(t(sapply(yy$promoters, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  tp$X3 = as.numeric(as.character(tp$X3)) - 1700
  tp$X2 = as.numeric(as.character(tp$X2)) -1
  rownames(tp) = rownames(yy)
  tp$promoter = paste0(tp$X1, ':', tp$X2, '-', tp$X3)
  
  mm = match(tp$promoter, rownames(motif.oc))
  motif.oc = motif.oc[mm, ]
  rownames(motif.oc) = rownames(yy)
  motif.oc = as.matrix(motif.oc)
  ss = apply(motif.oc, 2, sum)
  
  xx = as.data.frame(unclass(motif.oc))
  yy = cbind(yy, xx)
  
  saveRDS(yy, file = paste0(RdataDir, '/chromatin_promoter_features_Tfmotifs_geneGroups.rds')) 
  
}

FeatureImportance.RF = function()
{
  library(randomForest)
  library("iml")
  require(ggplot2)
  require(ggrepel)
  
  #yy = readRDS(file = paste0(RdataDir, '/chromatin_promoter_features_geneGroups.rds'))
  yy = readRDS(file = paste0(RdataDir, '/chromatin_promoter_features_Tfmotifs_geneGroups.rds')) 
  
  X = yy[which(yy$groups == 'DE_up'|yy$groups == 'DE_down'), ]
  X$groups = gsub('DE_', 'DE', X$groups)
  
  ## check the relevant features
  Feature.checing = FALSE
  if(Feature.checing){
    examples.sel = unique(grep(paste0(dev.example, collapse = '|'), X$gene))
    
    ggplot(data = X, aes(x = rna_mUA, y = H3K27me3_mUA, color = groups, label = gene)) +   
      geom_point(size = 0.4) +
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12)) + 
      #geom_text_repel(data=subset(yy, pvalue_pos > 2), size = 4)
      geom_text_repel(data= X[c(examples.sel), ], size =5)
    
    
    ggplot(data = X, aes(x = H3K4me3_mUA, y = H3K27me3_mUA, color = groups, label = gene)) +   
      geom_point(size = 0.4) +
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12)) +
      geom_text_repel(data= X[c(examples.sel), ], size =5)
    
    ggplot(data = X, aes(x = H3K4me3_mUA, y = cpg.oe, color = groups, label = gene)) +   
      geom_point(size = 0.4) +
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12)) +
      geom_text_repel(data= X[c(examples.sel), ], size =5)
    
    ggplot(data = X, aes(x = H3K27me3_mUA, y = cpg.oe, color = groups)) +   
      geom_point(size = 0.4) +
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12))
    
    ggplot(data = X, aes(x = rna_mUA, y = cpg.oe, color = groups)) +   
      geom_point(size = 0.4) +
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12))
    #geom_text_repel(data=subset(yy, pvalue_pos > 2), size = 4)
    #geom_text_repel(size = 4)
  }
  
  ## test some interactions
  y = as.factor(X$groups)
  x = data.frame(X[, c(4:9, 11:ncol(X))])
  #x$K27me3.K4me3 = x
  
  tic()
  set.seed(42)
  rf <- randomForest::randomForest(x = x, y = y, ntree = 200, keep.forest = FALSE, 
                                   importance = TRUE)
  toc()
  
  varImpPlot(rf, type = 1, n.var = 10)
  
  # plot the feature importance to distinguish 
  imps = importance(rf, type = 1)
  imps = data.frame(imps)
  imps$names = rownames(imps)
  colnames(imps)[1] = 'scores'
  imps$rank = rank(imps$scores)
  imps = imps[order(-imps$scores), ]
  
  saveRDS(imps, file = paste0(RdataDir, '/RF_featuresImportance.rds'))
  
  ggplot(data = imps, aes(x = scores, y = rank, label = names)) +   
    geom_point(size = 3.0, color = 'blue') +
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) + 
    #geom_text_repel(data=subset(yy, pvalue_pos > 2), size = 4)
    geom_text_repel(size = 4)
  
  
  gp = ggplot(data = yy, aes(x = scores, y = rank, label = names)) +   
    geom_point(size = 2.0, color = 'blue') +
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) + 
    geom_text_repel(size = 4) + ggtitle('Importance score from Random Forest MARA')
  
  
  
  # predictor <- Predictor$new(rf, data = x, y = y)
  # imp <- FeatureImp$new(predictor, loss = "mae")
  # plot(imp)
  # 
  # interact <- Interaction$new(predictor)
  # plot(interact)
  # 
  # interact <- Interaction$new(predictor, feature = "crim")
  # plot(interact)
  # 
  
  
}

plot_rna_chromainFeatures_geneExamples = function(tss, 
                                                  geneList = c('FGF10', 'SALL4'), 
                                                  outDir = 'Gene_Examepls')
{
  source('Functions_histM.R')
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library("cowplot")
  require(gridExtra)
  library(tidyr)
  require(patchwork)
  
  # outDir = "/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/Gene_Examples"   
  # geneList = dev.genes
  
  if(!dir.exists(outDir)) dir.create(outDir)
  features = c('rna', 'atac', 'H3K4me3', 'H3K27me3', 'H3K4me1')
  samples = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
  
  for(n in 1:length(geneList))
  {
    # n = 1
    kk = which(tss$gene == geneList[n])
    if(length(kk) == 1){
      cat(n, ' -- ', geneList[n], '\n')
      test = matrix(NA, nrow = length(samples), ncol = length(features))
      rownames(test) = samples
      colnames(test) = features
      for(m in 1:ncol(test))
      {
        if(m == 1) {
          mm = grep('smartseq2_', colnames(tss))
          test[,m] = as.numeric(tss[kk, mm[order(mm)]])
        }
        if(m == 2){
          mm = grep(paste0(colnames(test)[m], '_mUA|', colnames(test)[m], '_X'), colnames(tss))
          test[,m] = as.numeric(tss[kk, mm[order(mm)]]) 
        }
        if(m > 2){
          mm = grep(paste0(colnames(test)[m], '_mUA|', colnames(test)[m], '_BL'), colnames(tss))
          test[ ,m] = apply(matrix(as.numeric(tss[kk, mm[order(mm)]]), nrow = 2), 2, mean)
        }
      }
      
      test = data.frame(test, sample = rownames(test))
      #test[which(test[,1]<(-2)),1] = -2
      test[,1] = (test[, 1] - min(test[, 1]))/(max(test[,1]) - min(test[, 1])) * 3
      ylims = range(test[, c(1:5)])
      as_tibble(test) %>%  gather(features, signals,  1:5) %>% 
        mutate(features = factor(features, levels=c('rna', 'atac', 'H3K4me3', 'H3K27me3', 'H3K4me1'))) %>%
        mutate(sample = factor(sample, levels = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d'))) %>%
        ggplot(aes(y=signals, x=sample, color = features, group = features)) + 
        geom_line(aes(linetype=features, color=features), size = 1.2) +
        geom_point(size = 4.0) +
        theme_classic() +
        geom_hline(yintercept=c(0), col="darkgray", size = 1.2) +
        theme(axis.text.x = element_text(angle = 0, size = 14), 
              axis.text.y = element_text(angle = 0, size = 14), 
              axis.title =  element_text(size = 14),
              legend.text = element_text(size=12),
              legend.title = element_text(size = 14)
              )+
        scale_color_manual(values=c('black', 'springgreen', 'blue', 'red','gold2')) +
        scale_linetype_manual(values=c("longdash", "solid", "solid", "solid", "longdash")) + 
        labs( x = '', y = 'feature signals') +
        ylim(ylims[1], ylims[2]) +
        ggtitle(geneList[n])
        
      ggsave(paste0(outDir, "/Regeneration_allFeatures_", geneList[n],  ".pdf"),  width = 7, height = 4)  
      
    }else{
      cat(n, ' -- ', geneList[n], 'FOUND in tss',  length(kk), '\n')  
    }
  }
   
}

########################################################
########################################################
# Section : another section
# 
########################################################
########################################################

##########################################
# test the comparisons of histmarker across group which are redefined RNA-seq data 
##########################################
Redefine.gene.groups.with.RNAseq = function()
{
  #aa = readRDS(file =  "../results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/TestStat_regeneration_RNAseq.rds") 
  aa[grep('DBP|SALL', rownames(aa)), grep('Mature_UA', colnames(aa))]
  
  aa$expr.mUA = apply(aa[, grep('Mature_UA', colnames(aa))], 1, mean)
  aa$expr.BLday5 = apply(aa[, grep('BL_UA_5days', colnames(aa))], 1, mean)
  aa$expr.BLday9 = apply(aa[, grep('BL_UA_9days', colnames(aa))], 1, mean)
  aa$expr.BLday13 = apply(aa[, grep('BL_UA_13days_proximal', colnames(aa))], 1, mean)
  aa$expr.BLday13.distal = apply(aa[, grep('BL_UA_13days_distal', colnames(aa))], 1, mean)
  aa = aa[, c(5:11, 23:29)]
  aa = aa[match(res$geneID, aa$geneID), ]
  res = data.frame(res, aa, stringsAsFactors = FALSE)
  
  res$x1 = res$UA_K4me3
  res$x2 = res$UA_K27me3
  res$ratio = res$x1 - res$x2
  
  examples.sel = c()
  examples.sel = unique(c(examples.sel, grep('HOXA13|HOXA11|HOXA9|HOXD13|HOXD11|HOXD9|MEIS|SALL|DBP', res$gene)))
  
  
  ggplot(data=res, aes(x=x1, y=x2, label = gene)) +
    geom_point(size = 0.25) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
    #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=4, col='darkgray') +
    geom_hline(yintercept=4, col="darkgray") +
    labs(x = "UA_H3K4me3", y= 'UA_H3K27me3')
  
  
  ggplot(data=res, aes(x=ratio, y=expr, label = gene)) +
    geom_point(size = 0.25) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
    #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=4, col='darkgray') +
    geom_hline(yintercept=4, col="darkgray") +
    labs(x = "UA_H3K4me3/UA_H3K27me3", y= 'Expr')
  
  load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
  hkgs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
  nonexp = controls.tissue$geneIDs[which(controls.tissue$tissues != 'housekeeping')]
  
  res$groups = 'limb'
  res$groups[!is.na(match(res$geneID, hkgs))] = 'house_keep'
  res$groups[!is.na(match(res$geneID, nonexp))] = 'other_tissues'
  
  xx = res[order(-res$expr), ]
  head(xx[which(xx$groups != 'house_keep'), ])
  
  ggplot(data=res, aes(x=expr.mUA, y=ratio, label = gene, color = groups)) +
    geom_point(size = 0.25) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
    #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=4, col='darkgray') +
    geom_hline(yintercept=4, col="darkgray") +
    labs(y = "UA_H3K4me3/UA_H3K27me3", x= 'Expr')
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'house_keep', 'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  
  table(res$groups)
  kk = which(res$groups == 'limb')
  
  ss = apply(res[kk, grep('expr', colnames(res))], 1, max)
  
  select = kk[which(ss< -1 | is.na(ss))]
  res$groups[select] = 'lowlyExpr_limb'
  table(res$groups)
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  #jj = which(res$x1>4 & res$x2 >4)
  #xx = res[jj, ]
  
  kk = which(res$groups == 'limb')
  
  diffs = res$expr.BLday5[kk] - res$expr.mUA[kk]
  select = kk[which(res$expr.mUA[kk] < 0 & diffs >2)]
  
  res$groups[select] = 'upregulated.d5'
  table(res$groups)
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) >2)]
  res$groups[select] = 'downregulated.d5'
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 'upregulated.d5', 'downregulated.d5', 'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  
  kk = which(res$groups == 'limb')
  
  select = kk[which(res$expr.mUA[kk] < 0 & (res$expr.BLday5[kk] - res$expr.mUA[kk]) < 2 &  (res$expr.BLday9[kk] - res$expr.mUA[kk]) > 2)]
  res$groups[select] = 'upregulated.d9'
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) < 2 & (res$expr.mUA[kk] - res$expr.BLday9[kk]) >2)]
  res$groups[select] = 'downregulated.d9'
  
  
  kk = which(res$groups == 'limb')
  
  select = kk[which(res$expr.mUA[kk] < 0 & (res$expr.BLday5[kk] - res$expr.mUA[kk]) < 2 &  (res$expr.BLday9[kk] - res$expr.mUA[kk]) < 2 &
                      (res$expr.BLday13[kk] - res$expr.mUA[kk]) >2)]
  res$groups[select] = 'upregulated.d13'
  
  select = kk[which(res$expr.mUA[kk] > 0 & (res$expr.mUA[kk] - res$expr.BLday5[kk]) < 2 & (res$expr.mUA[kk] - res$expr.BLday9[kk]) < 2 &
                      (res$expr.mUA[kk] - res$expr.BLday9[kk]) > 2) ]
  res$groups[select] = 'downregulated.d13'
  
  
  
  
  table(res$groups)
  
  res[,c(1, 4, 7:ncol(res))] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'lowlyExpr_limb',
                                             'house_keep', 
                                             'upregulated.d5', 'upregulated.d9','upregulated.d13', 
                                             'downregulated.d5',  'downregulated.d9',
                                             'downregulated.d13',
                                             'limb')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
  
  ggsave(paste0(figureDir, "histMarker_H3K27me3_H3K4me3_up.downregulatedGenes.UA.BL.days.pdf"), width=12, height = 8)
  
  fdr.cutoff = 0.05
  select = which(aa$padj < fdr.cutoff & aa$log2FC >1 & aa$log2FC.mUA.vs.others >0)
  res$groups[!is.na(match(res$geneID, aa$geneID[select]))] = 'mature_highlyExp'
  table(res$groups)
  
  select = which(aa$padj < fdr.cutoff & aa$log2FC >1 & aa$log2FC.mUA.vs.others < 0)
  res$groups[!is.na(match(res$geneID, aa$geneID[select]))] = 'regeneration'
  table(res$groups)
  
  select = which((aa$padj >= fdr.cutoff | aa$log2FC <=1) & log2(aa$baseMean) <4) 
  res$groups[!is.na(match(res$geneID, aa$geneID[select])) & res$groups == 'limb_static'] = 'limb_lowlyExp'
  
  
  ggplot(data=res, aes(x=x1, y=x2, label = gene, color = groups)) +
    geom_point(size = 0.4) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    geom_text_repel(data= res[examples.sel, ], size = 3.0, color = 'blue') +
    #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel), size = 2) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=4, col='darkgray') +
    geom_hline(yintercept=4, col="darkgray") +
    labs(x = "UA_H3K4me3", y= 'UA_H3K27me3')
  
  
  #xx = melt(res[], id.vars = c('UA_K4me3', 'UA_K27me3'), variable_name = 'markers')
  res[,c(1, 4, 7:13)] %>% 
    pivot_longer(cols = c('UA_K4me3', 'UA_K27me3'), names_to = 'markers') %>%
    ggplot(aes(x = factor(groups, levels = c('other_tissues', 'house_keep', 'limb_lowlyExp',
                                             'mature_highlyExp', 'regeneration', 'limb_static')), y=value, fill=markers)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    #geom_jitter(width = 0.1)+
    #geom_violin(width = 1.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14)) +
    labs(x = "", y= 'normalized data (log2)')
  
}

Test.other.mark.combinations.besides.bivalency = function()
{
  ##########################################
  # test other markers
  ##########################################
  res$x = res$H3K27me3_mUA
  res$y = res$H3K4me1_mUA
  ggplot(data=res, aes(x=x, y=y, label = gene)) +
    geom_point(size = 0.1, color = 'darkgray') + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    #scale_color_manual(values=c('black', "orange", 'darkgray',  "red",   'green')) + 
    #geom_point(data=res[res$groups == 'reg_up', ], aes(x=H3K4me3_mUA, y=H3K27me3_mUA),  size=0.7) +
    #geom_point(data=res[res$groups == 'reg_down', ], aes(x=H3K4me3_mUA, y=H3K27me3_mUA),  size=0.7) +
    geom_point(data=res[examples.sel, ], aes(x=x, y=y),  size=1.5, color = 'blue') +
    geom_text_repel(data= res[examples.sel, ], size = 4.0, color = 'blue') +
    geom_point(data=res[matures.sel, ], aes(x=x, y=y),  size=1.5, color = 'darkred') +
    geom_text_repel(data= res[matures.sel, ], size = 4.0, color = 'darkred') +
    
    #geom_text_repel(data= fpm[examples.sel, ], size = 4.0, color = 'darkblue') + 
    #geom_hline(yintercept=2.0, colour = "darkgray") + 
    #geom_vline(xintercept = 2.0, colour = "darkgray")
    #geom_abline(slope = 1,  intercept = 0, colour = 'cyan3') +
    theme_classic() +
    theme(legend.text = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.position=c(0.1, 0.8),
          plot.margin = margin()
          #legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    ) + 
    geom_vline(xintercept=1, col='orange') +
    geom_hline(yintercept=0, col="orange") +
    labs(x = "UA_H3K27me3", y= 'UA_H3K4me1') +
    guides(colour = guide_legend(override.aes = list(size=2)))
  
}

