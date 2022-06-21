library(ChIPseeker)
library(rtracklayer)
#library(UpSetR)
library("ChIPpeakAnno")
library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library("ggplot2")
library("GenomicFeatures")

########################################################
########################################################
# Section: legacy codes from Thomas Burkard
########################################################
########################################################
venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  
  Venn(SetNames=SetNames, Weight=Weight)
}

venn_cnt2barplot <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts")

  cnt <- venn_cnt[,(n+1):length(colnames(venn_cnt))]
  rownames(cnt) <- cbind(apply(venn_cnt[,1:(n-1)], 1, paste, collapse="|"))
  cnt <- cnt[rownames(cnt) != "0|0|0", ]
  cnt <- melt(cnt)
  colnames(cnt) <- c("intersections", "condition", "peaks")
  legend.title <- paste(colnames(venn_cnt[,1:(n-1)]), collapse="|")

  p <- ggplot(cnt, aes(x=condition, y=peaks, fill=intersections)) +  geom_bar(stat="identity", colour="grey50", alpha = 0.9) + scale_fill_brewer(palette="Set3") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank()) + labs(fill=legend.title)

  plot(p)  
}


filterPeaks <- function(files, pval = 10, qval=0, remove.blacklist=FALSE, barplot = FALSE, fileType="macs2", input1B=TRUE, addChr=FALSE)
{
  # files = macs2.files;remove.blacklist=FALSE; barplot = FALSE; fileType="macs2"; input1B=TRUE; addChr=FALSE; pval = 100; qval=0;
  bar <- NULL
  peaks <- NULL
  for (f in files)
  {
    p <- readPeakFile(f, as = "GRanges")
    if (addChr)
    {
      seqlevels(p) <- paste0("chr", seqlevels(p))
    }
    if (input1B)
    {
      start(p) <- start(p) - 1
    }
    all <- length(p)
    if (fileType == "macs14")
    {
      p <- p[mcols(p)[,"X.10.log10.pvalue."] > pval]
    } else if (fileType == "macs2") {
      if (qval == 0)
      {
        p <- p[mcols(p)[,"X.log10.pvalue."] > pval]
      } else {
        p <- p[mcols(p)[,"X.log10.qvalue."] > qval]
      }
    }
    p100 <- length(p)
    if(remove.blacklist)
    {
      p <- p[!overlapsAny(p, blacklist)]
      p.bl <- p100 - length(p)
    }
    all <- all - p100
    p100 <- length(p)
    
    if(remove.blacklist){
      bar <- cbind(bar, c(all, p100, p.bl))
    }else{
      bar <- cbind(bar, c(all, p100))
    }
    
    peaks <- c(peaks, p)
  }
  
  names(peaks) <- names(files)
  
  if (barplot)
  {
    if (qval !=0)
    {
      rownames(bar) <- c("all", paste0("q", pval), "blacklisted")
    } else {
      rownames(bar) <- c("all", paste0("p", pval), "blacklisted")
    }
    colnames(bar) <- names(peaks)
    p <- ggplot(melt(bar), aes(x=X2, y=value, fill=X1)) + geom_bar(stat="identity", colour="black") + scale_fill_brewer() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.title=element_blank())
    plot(p)
  }
  
  return(peaks)  
  
}


plotV <- function(peaks, main = "")
{
    if (length(peaks) > 5)
    {
        peaks <- peaks[1:5] #up to five!!!
    }
    ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=names(peaks), maxgap=0, minoverlap =1, main=main, connectedPeaks="keepAll")
 
    v <- venn_cnt2venn(ol.peaks$vennCounts)
    try(plot(v))
    #plot(v, doWeights=FALSE)
    #plot(v, doEuler = TRUE)
    venn_cnt2barplot(ol.peaks$vennCounts)
    return(ol.peaks)
}


doAllPeak <- function(set)
{
  if (set == "macs14")
  {
    allPeaks <- filterPeaks(macs14.files, barplot=TRUE)
    allPeaks.50 <- filterPeaks(macs14.files, pval=50, blacklist=blacklist)
  } else if (set == "macs2")
  {
    allPeaks <- filterPeaks(macs2.files, pval = 10, blacklist=blacklist, barplot=TRUE, fileType="macs2")
    allPeaks.50 <- filterPeaks(macs2.files, pval = 5, blacklist=blacklist, barplot=TRUE, fileType="macs2")
  } else if (set == "macs2qval")
  {
    allPeaks <- filterPeaks(macs2.files, qval = 10, blacklist=blacklist, barplot=TRUE, fileType="macs2")
    allPeaks.50 <- filterPeaks(macs2.files, qval = 5, blacklist=blacklist, barplot=TRUE, fileType="macs2")
  } else if (set == "rep")
  {
    allPeaks <- filterPeaks(rep.files, blacklist=blacklist, barplot=TRUE)
    allPeaks.50 <- filterPeaks(rep.files, pval=50, blacklist=blacklist)
  }
  
  promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000) #ChIPseeker method
  tagMatrixList <- lapply(allPeaks, getTagMatrix, windows=promoter)
  
  tagMatrixList <- tagMatrixList[!sapply(lapply(tagMatrixList, nrow), is.null)]
  tagMatrixList[sapply(tagMatrixList, nrow) == 0] <- NULL
  
  #promoter <- promoters(genes(txdb), upstream=5000, downstream=5000) #GRanges
  #tagMatrixList <- lapply(allPeaks, getTagMatrix, windows=promoter)
  #plotAvgProf(tagMatrixList, xlim=c(-5000, 5000))
  #tagHeatmap(tagMatrixList, xlim=c(-5000, 5000), color=NULL)
  
  peakAnnoList <- lapply(allPeaks, annotatePeak, TxDb=txdb, tssRegion=c(-2000, 200), verbose=FALSE)
  
  return(list(allPeaks, allPeaks.50, peakAnnoList, tagMatrixList))
  
}


doAllPlot <- function(set, allPeaksList)
{
    allPeaks <- allPeaksList[[1]]
    allPeaks.50 <- allPeaksList[[2]]
    peakAnnoList <- allPeaksList[[3]]
    tagMatrixList <- allPeaksList[[4]]
    
    if (set == "rep")
      {
        plotV(allPeaks, main = "Ints (< 1e-100)") #different order due to error otherwise ...
        plotV(allPeaks.50, main = "Ints (< 1e-50)")
      } else {
    
        plotV(allPeaks, main = "Ints (< 1e-100)") #different order due to error otherwise ...
        plotV(allPeaks.50, main = "Ints (< 1e-50)")
      }

#covplot(allPeaks[[9]], weightCol="V5", chrs = c(paste0("chr", c(1:19, "X", "Y"))))


#peakAnnoList <- list("Normal"=annotatePeak(allPeaks[[1]],  TxDb=txdb,tssRegion=c(-2000, 200), verbose=FALSE), 
    #"TSS"=annotatePeak(allPeaks[[1]],  TxDb=txdb,tssRegion=c(0, 0), verbose=FALSE), "5kb"=annotatePeak(allPeaks[[1]],  
    #TxDb=txdb,tssRegion=c(-5000,5000), verbose=FALSE), "3kb"=annotatePeak(allPeaks[[1]],  TxDb=txdb, verbose=FALSE))

    print(plotAnnoBar(peakAnnoList))
    print(plotDistToTSS(peakAnnoList))
    
    for (n in names(peakAnnoList))
      {
        par(mfrow=c(1,1))
        vennpie(peakAnnoList[[n]])
        text(x=0, y=-1, n)
        par(mfrow=c(1,1))
        plotAnnoPie(peakAnnoList[[n]])
#        par(mfrow=c(1,1))
#        upsetplot(peakAnnoList[[n]]) #only in TB-3.2.1-dev #  vennpie=TRUE,
      }
}

#for (set in c("macs14", "macs2", "macs2qval"))
#  {
#    pdf(paste0("peaksAnalysis.", set, ".pdf"), width = 10, height=10)
#    doAllPeak(set=set)
#    dev.off()
#  }

########################################################
########################################################
# Section: Quality control parts 
# including peak overlapping check
########################################################
########################################################
PLOT.Quality.Controls.Summary = function(stat, index)
{
  par(mfrow=c(2,2))
  ss = stat[index, ]
  samples = sapply(ss$filename, find.samples.conditions, ID='samples')
  conditions =  sapply(ss$filename, find.samples.conditions, ID='conditions')
  ss$conditions = conditions
  ss = data.frame(ss, stringsAsFactors = FALSE)
  ss = ss[with(ss, order(samples, conditions)), ]
  
  cols.conditions = data.frame(unique(conditions), c(1:length(unique(conditions))))
  cols = cols.conditions[match(ss$conditions, cols.conditions[,1]), 2]
  
  barplot(as.numeric(ss$nb.UniqueRead)/10^6, col=cols, main="nb of uniquely-mapped Reads", horiz=TRUE, names.arg = ss$filename, 
          las=1, xlim = c(0, 50))
  abline(v=c(10, 20), col='red', lty=1, lwd=2.0);abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
  
  pbc = ss$nb.UniqueRmdupRead/ss$nb.UniqueRead
  barplot(pbc, col=cols, main="NRF", horiz=TRUE, names.arg = ss$filename, 
          las=1, xlim = c(0, 1))
  abline(v=c(0.9), col='red', lty=1, lwd=2.0);
  
  
  plot(ss$NSC, ss$RSC, type='n', xlab='NSC', ylab='RSC', ylim=range(ss$RSC), xlim=range(ss$NSC))
  colfunc <- colorRampPalette(c("red", "green"))
  cols = colfunc(5)
  #jj = which(stat$samples !='Input')
  #stat = stat[jj,]
  for(n in 1:nrow(ss))
  {
    points(ss$NSC[n], ss$RSC[n], type='p', col=cols[(ss$Quality[n]+3)], pch=16, cex=2.0)
    #else  points(stat$NSC.1.05.[n], stat$RSC.0.8.[n], type='p', 
    #col=cols[(stat$QualityTag..2.verylow..1.low.0.medium.1.high.2.veryhigh.[n]+3)], pch=17, cex=2.0)
    text(ss$NSC[n], ss$RSC[n], n, cex=1.6, offset=0.5, pos = 1)
  }
  abline(h=c(0.8), col='red', lty=3, lwd=2.0);abline(v=c(1.05), col='red', lty=3, lwd=2.0)
  legend('topleft', legend=rev(c('very low', 'low', 'medium', 'high', 'very high')), col=rev(cols), pch=16, bty='n', cex=1.5)
  
  plot(ss$NSC, ss$RSC, type='n',ylim=c(0, 2.0), xlim=c(1, 1.5), xlab=NA, ylab=NA, axes=FALSE)
  
  if(nrow(ss)<=40){
    legend('topleft', legend=paste(c(1:nrow(ss)), ss$filename, sep='-'), col=cols[ss$Quality+3], pch=16, bty='n', cex=1.2)
  }else{
    legend('topleft', legend=paste(c(1:40), ss$filename[1:40], sep='-'), col=cols[ss$Quality+3][1:40], pch=16, bty='n', cex=1.)
    legend('topright', legend=paste(c(41:nrow(ss)), ss$filename[41:nrow(ss)], sep='-'), col=cols[ss$Quality+3][41:nrow(ss)], 
           pch=16, bty='n', cex=1.2)
  }
  
}

make.design.matrix.from.file.list = function(peak.files, varnames = c('condition','factor'),  reOrder = FALSE)
{
  cat("-- parsing design matrix\n")
  cat('-- ', length(peak.files), 'peak files \n')
  cat('-- experiment design conditions: \n')
  print(varnames)
  bname = basename(peak.files)
  
  design.matrix = data.frame(peak.files, bname, stringsAsFactors = FALSE)
  
  for(n in 1:length(varnames))
  {
    #test = bname[n];
    design.matrix = data.frame(design.matrix, sapply(bname, function(x) unlist(strsplit(as.character(x), "_"))[n]), 
                               stringsAsFactors = FALSE)
    #test = 
    #xx = c(xx, test[1])
    #test = unlist(strsplit(as.character(test), "[.]"))[-1]
    #test = paste0(test, collapse = ".")
    #yy = c(yy, test[2])
  }
  colnames(design.matrix) = c('file.path', 'file.name', varnames)
  
  #design = data.frame(sapply(bname, find.samples.conditions, ID='conditions'), sapply(bname, find.samples.conditions, ID='samples'),
  #                    stringsAsFactors = FALSE)
  
  
  #if(reOrder){
  #  design.matrix = design.matrix[with(design.matrix, order(factor, condition)), ];
  #}
  
  #factor.condition = paste0(design.matrix$factor, '_', design.matrix$condition)
  #design.matrix = data.frame(design.matrix, factor.condition, stringsAsFactors = FALSE)
  
  return(design.matrix)
  
}

Comparison.overlapping.peaks = function(design.matrix, peaks.list, toCompare="condition", pval=10, PLOT.p10 = FALSE, qval=0.1) 
{
  # design.matrix = design; peaks.list = peak.files; toCompare = 'condition'; pval = 10; PLOT.p10 = FALSE
  cat('-- Compare peak overlapping with ', toCompare, '--\n')
  index  = c()
  for(n in 1:nrow(design.matrix))
  {
    test = grep(design.matrix$sampleID[n], peaks.list)
    if(length(test) != 1) {
      cat(length(test), 'peak files Found \n')
    }else{
      index = c(index, test)
    }
  }
  peaks.list = peaks.list[index]
  
  if(toCompare == "factor"){
    cat("Warning ---- DB analysis for each factor across condition \n")
    cat("Warning ---- Should merge peaks for all replicates \n")
    
  }else{
    cat("checking peak overlapping for",  toCompare, "\n")
    if(toCompare=="factor.condition"){
     cat("--compare peak overlapping for replicates--\n") 
    }
    sels2compare = unique(design.matrix[, which(colnames(design.matrix)==toCompare)]);
    
  }
  
  for(nn in sels2compare)
  {
    # nn = "BL_UA_5days"
    
    cat(nn, '\n')
    kk = which(design.matrix[, which(colnames(design.matrix)==toCompare)]==nn)
    
    if(length(kk)>1)
    {
      peaks = c()
      peaks10 = c()
      peaknames = paste0(nn, '_', design.matrix$sampleID[kk])
      
      for(k in kk) 
      {
        p = readPeakFile(peaks.list[[k]], as = "GRanges");
        #eval(parse(text = paste0("p = pp.", k)));
        with.p.values = "X.log10.pvalue." %in% colnames(mcols(p))
        if(with.p.values) {
          p10 <- p[mcols(p)[,"X.log10.pvalue."] > pval];
          p10 = reduce(p10);
          peaks10= c(peaks10, p10);
        }else{ 
          cat("no p values conlumn found for -- ", design.matrix$file.name[k], "\n");
          PLOT.p10 = FALSE;
        }
        
        p = reduce(p)
        peaks= c(peaks, p)
      }
      
      ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=peaknames, connectedPeaks="keepAll", main=nn)
      v <- venn_cnt2venn(ol.peaks$vennCounts)
      try(plot(v))
      
      if(PLOT.p10){
        ol.peaks <- makeVennDiagram(peaks10, NameOfPeaks=paste(peaknames, '_p10', sep=''), connectedPeaks="keepAll", main=nn)
        v <- venn_cnt2venn(ol.peaks$vennCounts)
        try(plot(v))
      }
    }else{
      cat("Error --- less than 2 samples selected to compare for ", sels2compare[nn], "\n")
    }
  }
}



find.sample.names = function(x, n1=1, n2=3)
{
  temp = unlist(strsplit(as.character(x), '_'));
  kk = which(temp=='macs2')
  return(paste(temp[c(n1:n2)], sep='', collapse = '_'))
}


find.samples.conditions = function(x, ID='samples')
{
  if(ID=='samples') 
  { 
    j = 2
    temp = unlist(strsplit(as.character(x), '_'));
    return(temp[j])
  }
  if(ID=='conditions') 
  {
    j = 1;
    temp = unlist(strsplit(as.character(x), '_'));
    temp = temp[j]
    temp = unlist(strsplit(as.character(temp), '[.]'));
    return(temp[length(temp)])
    #kk = which(temp=='macs2')
    #return(paste(temp[c(n1:(kk-n2))], sep='', collapse = '_'))
  }
}

########################################################
########################################################
# Section: Differential Binding (DB) analysis after counting reads within peaks
########################################################
########################################################
merge.peaks.macs2 = function(peak.list, select.overlappingPeaks.acrossRep = FALSE,  cc = NULL, pcutoff = 4, 
                             select.overlappingPeaks = FALSE, rm.chrM = TRUE, merge.dist = NULL)
{
  require(GenomicRanges)
  # peak.list = peak.files; cc = design$condition;pcutoff = 6
  
  if(select.overlappingPeaks.acrossRep == TRUE & !is.null(cc)){
    cat('select overlapping peaks between replicates \n')
    cc.unique = unique(cc)
    for(n in 1:length(cc.unique))
    {
      jj = which(cc == cc.unique[n])
      cat(n, '--', cc.unique[n], ': ', length(jj),  'rep \n')
      
      ### find overlapping peaks between replicates
      px = readPeakFile(peak.list[jj[1]] , as = "GRanges")
      px <- px[mcols(px)[,"X.log10.pvalue."] > pcutoff]
      cat(length(px), 'peaks in one rep \n')
      
      if(length(jj) > 1){
        for(m in 2:length(jj))
        {
          pxx = readPeakFile(peak.list[jj[m]] , as = "GRanges")
          pxx <- pxx[mcols(pxx)[,"X.log10.pvalue."] > pcutoff]
          cat(length(pxx), 'peaks in one rep \n')
          #px = px[overlapsAny(px, pxx, type = 'any', ignore.strand=TRUE)]
          px = GenomicRanges::intersect(px, pxx, ignore.strand=TRUE)
        }
      }
      cat(length(px), ' peaks found across rep \n')
      
      if(rm.chrM){seqlevels(px, pruning.mode="coarse") <- seqlevels(px)[seqlevels(px) != "chrM"]; }
      
      # union of peaks from different conditions
      if(n == 1){
        peaks.merged <- px
      }else{
        peaks.merged = GenomicRanges::union(peaks.merged, px, ignore.strand=TRUE)
      }
    }
        
  }else{
    cat('merge all peaks \n')
    
    for(n in 1:length(peak.list)){
      cat(peak.list[n], "\n")
      px = readPeakFile(peak.list[n] , as = "GRanges")
      px <- px[mcols(px)[,"X.log10.pvalue."] > pcutoff]
      
      if(rm.chrM){seqlevels(px, pruning.mode="coarse") <- seqlevels(px)[seqlevels(px) != "chrM"]; }
      
      if(n == 1){
        peaks.merged <- px
      }else{
        if(!select.overlappingPeaks){
          peaks.merged = GenomicRanges::union(peaks.merged, px, ignore.strand=TRUE)
        }else{
          peaks.merged = peaks.merged[overlapsAny(peaks.merged, px, type= 'any', ignore.strand=TRUE)]
          px = px[overlapsAny(px, peaks.merged, type = 'any', ignore.strand=TRUE)]
          peaks.merged = GenomicRanges::union(peaks.merged, px, ignore.strand = TRUE)
          peaks.merged = reduce(peaks.merged)
          
        }
      }
    }
  }
  
  
  if(!is.null(merge.dist)){
    peaks.merged <- mergeWindows(peaks.merged, tol=merge.dist, ignore.strand = TRUE)
    peaks.merged = peaks.merged$region
    #peaks.merged = GenomicRanges::reduce(pps);
  }
  
  #peaks.merged = as.data.frame(peaks.merged)
  return(peaks.merged)
}


quantify.signals.forGenome.csaw = function(bam.list){
  library(csaw)
  library(edgeR)
  library(GenomicRanges)
  library(rtracklayer)
  
  #frag.len <- 110
  win.width <- 2000
  chr.selected = c(paste0("chr", c(1:19)))
  param <- readParam(pe="none", dedup = TRUE, minq=30, restrict = chr.selected)
  
  binned <- windowCounts(bam.files, bin=TRUE, width=win.width, param=param)
  
  ## save the big file in case it is lost 
  save(binned, baits, file = paste0(NormDir, "bam_unique_rmdup_readCounts_windows_csaw_forNormalization.Rdata"))
  
}

quantify.signals.within.peaks = function(peaks, bam.list, rpkm.normalization = FALSE, isPairedEnd = FALSE)
{
  require(Rsubread)
  
  ## peak regions (configure your peak regions into a data.frame)
  peaks = data.frame(peaks)
  colnames(peaks)[c(1:3)] = c("chr", "start", "end")
  peaks$peak.name = paste0(peaks$chr, "_", peaks$start, "_", peaks$end)
  jj = match(unique(peaks$peak.name), peaks$peak.name)
  df = peaks[jj, ];
  
  SAF = data.frame(GeneID=df$peak.name, Chr=df$chr, Start=df$start, End=df$end, Strand="+", stringsAsFactors = FALSE)
  
  ## count reads for those peak regions using 
  fc <- featureCounts(files=bam.list, annot.ext = SAF, countMultiMappingReads = FALSE, minMQS = 30, 
                      ignoreDup = TRUE, strandSpecific = 0, juncCounts = FALSE, nthreads = 6, isPairedEnd = isPairedEnd)
  
  if(!rpkm.normalization){
    
    res = fc;
    
  }else{
    stat = fc$stat;
    counts = fc$counts;
    counts.annot = fc$annotation
    
    rpkm = matrix(NA, ncol = ncol(counts), nrow = nrow(counts))
    colnames(rpkm) = colnames(counts)
    row.names(rpkm) = rownames(counts)
    kk = which(stat$Status=="Assigned" | stat$Status== "Unassigned_NoFeatures")
    
    for(n in 1:ncol(counts))
    {
      #n = 1 
      jj = which(colnames(stat) == colnames(counts)[n])
      rpkm[, n] = (counts[, n])/counts.annot$Length/sum(stat[kk, jj])*10^9
    }
    
    #rpkm = log2(rpkm)
    res = data.frame(SAF, Length=counts.annot$Length[match(SAF$GeneID, counts.annot$GeneID)], 
                     rpkm[match(SAF$GeneID, rownames(rpkm)), ], stringsAsFactors = FALSE)
  }
  
  return(res)
  
}

## inputs are counts, design.matrix
DB.analysis= function(counts, design.matrix, size.factors = NULL, batch = FALSE, Threshold.read.counts = 20, cex.pairwise = 0.01) 

{
  require(lattice);
  require(ggplot2)
  require('DESeq2');
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr"); 
      
  # design.matrix = design.matrix[, kk];
  # size.factors = NULL;
  # Threshold.read.counts = 50
  
  design.matrix = data.frame(design.matrix)
  if(ncol(design.matrix)==1) {
    colnames(design.matrix) = "conds";
    #newO = order(design.matrix$conds)
    #design.matrix = design.matrix[newO, ]
    #counts = counts[, newO]
    
    dds = DESeqDataSetFromMatrix(counts, DataFrame(design.matrix), design = ~ conds)
    conditions = design.matrix$conds; 
  }else{
    conds = factor(paste0(colnames(design.matrix), collapse = " + "))
    #newO = order(design.matrix$conds)
    #design.matrix = design.matrix[newO, ]
    #counts = counts[, newO]
    
    conditions = apply(design.matrix, 1, function(x) {paste0(x, collapse = "_")})
    eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(counts, DataFrame(design.matrix), design = ~ ", conds, ")")))
  }
  
  cc.uniq = unique(conditions);
  cols = match(conditions, cc.uniq)
  
  # show the raw counts
  par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(6,16,2,0.8)+0.1, tcl = -0.3)
  par(mfrow=c(1,1))
  
  total = colSums(counts(dds));
  barplot(total/10^6, horiz = TRUE, names.arg = colnames(raw), las=1, col = cols, main='Total nb of reads quantified for features', xlab='number of reads (Million)')
  abline(v=c(1, 2, seq(5, 20, by=5)), col='red', lty=1, lwd=2.0);#abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
    
  # filter the lowly signal peaks
  dds <- dds[ rowSums(counts(dds)) > Threshold.read.counts, ]
  if(!is.null(size.factors)) {
    sizeFactors(dds) = size.factors
  }else{
    dds <- estimateSizeFactors(dds)
  }
  fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
 
  xx = fpm;
  ### boxplot (distribution) of gene expression for each sample
  par(mfrow=c(1,1))
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(15,3,2,0.2), tcl = -0.3)
  for(n in 1:ncol(xx)){
    kk = which(xx[,n]>0);
    if(n==1) boxplot(log2(xx[kk, n]), horizontal = FALSE, las=2, at=(n), ylim=c(-6, 23), xlim=c(0, (ncol(xx)+1)), names=as.character(colnames(xx)[n]),
                     las=1, width=0.6, ylab='log2(norm.deseq2)', col=cols[n], main="Distribution of normalized signals (deseq2)")
    else boxplot(log2(xx[kk, n]), horizontal = FALSE, las=1, add=TRUE, at=(n), names=colnames(xx)[n], width=0.6, col=cols[n])
    mtext(colnames(xx)[n], side = 1, at=n, las=2)
  }
  
  ### pairwise correlations for fpm
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  library(corrplot)
  col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                             "cyan", "#007FFF", "blue","#00007F"))
  xx = as.matrix(fpm)
  xx[which(xx==0)] = NA
  M <- cor(xx, use = "na.or.complete")
  #corrplot(M, method="circle", type = 'upper', order="hclust")
  corrplot(M, method="ellipse", order="hclust", tl.cex=1.2, cl.cex=0.7, tl.col="black", addrect=ceiling(ncol(xx)/2), col=col1(100), rect.col=c('green'), rect.lwd=2.0)
  
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    #as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  colnames(df)[1:2] <- c("x", "y")
  vsd.transform=ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  #print(vsd.transform)
  sampleDists <- dist(t(assay(vsd)))
  #sampleDists
  #rld <- rlog(dds, blind=FALSE)
  sampleDistMatrix <- as.matrix( sampleDists)
  #rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  #plot(distClustering)
  if(ncol(design.matrix)>1){
    pca=plotPCA(vsd, intgroup = c(colnames(design.matrix)), returnData=FALSE)
  }else{
    pca=plotPCA(vsd, intgroup = c('conds'), returnData=FALSE)
  }
  plot(pca);
  
  Show.sample.names.PCA.Clusters = TRUE
  if(Show.sample.names.PCA.Clusters){
    if(ncol(design.matrix) == 1){
      pca2save = as.data.frame(plotPCA(vsd, intgroup =c('conds'), returnData = TRUE))
      #p = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=condition, shape=time)) + geom_point(size=3)
      #p + geom_text(hjust = 0.5, nudge_y = 0.1, size=2.5) 
      ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=conds)) + geom_point(size=3) +
        geom_text(hjust = 0.3, nudge_y = 0.5, size=2.5)
      plot(ggp);
    }else{
      if(!batch){
        pca2save = as.data.frame(plotPCA(vsd, intgroup =c(colnames(design.matrix)), returnData = TRUE))
        conditions = apply(design.matrix, 1, function(x) {paste0(x, collapse = "_")})
        pca2save = data.frame(pca2save, conditions, stringsAsFactors = FALSE)
        
        ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, col = conditions)) + geom_point(size=3) +
          geom_text(hjust = 0.3, nudge_y = 0.5, size=2.5)
        plot(ggp)
      }else{
        pca2save = as.data.frame(plotPCA(vsd, intgroup =c(colnames(design.matrix)), returnData = TRUE))
        
        jj = which(colnames(design.matrix)=='batch')
        if(length(jj)==0){
          cat('Error: no batch colunm found in design matrix \n')
        }else{
          if(ncol(design.matrix) == 2){
            conditions = factor(setdiff(colnames(design.matrix), 'batch'))
            pca2save = data.frame(pca2save, stringsAsFactors = TRUE)
          }else{
            conditions = apply(design.matrix[,-jj], 1, function(x) {paste0(x, collapse = "_")})
            pca2save = data.frame(pca2save, conditions, stringsAsFactors = TRUE)
          }
          pca2save$batch = factor(as.character(pca2save$batch))
          eval(parse(text = paste0("ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, col =", conditions, ", shape = batch)) + geom_point(size=3) +
            geom_text(hjust = 0.3, nudge_y = 0.5, size=2.5)")))
          plot(ggp)
        } 
      }
    }
  }
  
  ###  pairwise correlation and fitting (to chek if there is batch effect)
  yy = as.matrix(fpm)
  #yy = as.matrix(assay(vsd))
  #yy = as.matrix(xx[, c(1:3)])
  yy[which(yy==0)] = NA;
  yy = log2(yy)
  
  pairs(yy[, order(design.matrix$batch)], lower.panel=NULL, upper.panel=panel.fitting)
  #pairs(assay(vsd), lower.panel = NULL, upper.panel = panel.fitting)
  
  #cc.partA = c("AN312", "515D10", "515D10H3")
  #cc.partB = c("AN312", "515D10", "D10A8", "D10D8")
  #cc.partC = c("924E12", "515D10", "E12F01")
  
  #kk = match(dds$conds, cc.partA); ii = which(!is.na(kk));
  #if(length(unique(dds$conds[ii]))>2) pairs(yy[, ii], lower.panel=NULL, upper.panel=panel.fitting)
  
  #kk = match(dds$conds, cc.partB); ii = which(!is.na(kk));
  #if(length(unique(dds$conds[ii]))>2) pairs(yy[, ii], lower.panel=NULL, upper.panel=panel.fitting)
  
  #kk = match(dds$conds, cc.partC); ii = which(!is.na(kk));
  #if(length(unique(dds$conds[ii]))>2) pairs(yy[, ii], lower.panel=NULL, upper.panel=panel.fitting)
  
  
  return(dds)
  
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

panel.fitting = function (x, y, bg = NA, pch = par("pch"), cex = 0.05, col='black') 
{
  #x = yy[,1];y=yy[,2];
  #kk = which(x>0 & y>0); x=x[kk];y=y[kk]
  lims = range(c(x, y), na.rm = TRUE)
  points(x, y, pch = 1, col = col, cex = cex, xlim=lims, ylim=lims)
  abline(0, 1, lwd=1.0, col='red')
  R = cor(x, y, use="na.or.complete", method='pearson')
  text(lims[2]*0.2, lims[2]*0.9, paste0('R = ', signif(R, d=2)), cex=0.5, col='red')
  jj = which(!is.na(x) & !is.na(y))
  #fit = lm(y[jj] ~ x[jj])
  #slope=summary(fit)$coefficients[1]
  #slope = fit$coefficients[2]
  #intercept = fit$coefficients[1]
  #pval=summary(fit)$coefficients[4]
  #abline(intercept, slope, lwd=1.2, col='darkblue', lty=3)
  #text(lims[2]*0.1, lims[2]*0.7, paste0('slop = ', signif(slope, d=2)), cex=1., col='blue')
  #text(lims[2]*0.1, lims[2]*0.6, paste0('pval = ', signif(pval, d=2)), cex=1., col='blue')
  #ok <- is.finite(x) & is.finite(y)
  #if (any(ok)) 
  #lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
  #        col = col.smooth, ...)
}


##########################################
# remove batch and check batch removal  
##########################################
remove.batch.atacseq = function(cpm, design.matrix, method = "combat")
{
  # cpm = xx
  #logcpm = log2(cpm + 2^-6)
  
  if(method == 'limma'){
    cat('remove the batch effect using limma \n')
    require('limma')
    #design.tokeep = design.matrix
    #design.tokeep$tissue.cell[which(design.tokeep$treatment == "untreated")] = 'whole.body'
    #design.tokeep$tissue.cell[which(design.tokeep$genotype=="N2" & design.tokeep$treatment=="treated")] = "background"
    batch = design.matrix$batch
    design.tokeep<-model.matrix(~ 0 + factor.condition,  data = design.matrix)
    cpm.bc = removeBatchEffect(cpm, batch = design.matrix$batch, design = design.tokeep)
  }
  
  if(method == 'combat'){
    ## here we use the combat to remove the batch effect 
    ## the combat requires the log2cpm
    cat('remove the batch effect using ComBat \n')
    require("sva")
    # example from the ComBat function in the R package 'sva'
    TEST.example = FALSE
    if(TEST.example){
      library(bladderbatch)
      data(bladderdata)
      dat <- bladderEset[1:50,]
      pheno = pData(dat)
      edata = exprs(dat)
      batch = pheno$batch
      mod = model.matrix(~as.factor(cancer), data=pheno)
      combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
    }
    
    batch = design.matrix$batch;
    design.tokeep = design.matrix
    design.tokeep$tissue.cell[which(design.tokeep$treatment == "untreated")] = 'whole.body'
    design.tokeep$tissue.cell[which(design.tokeep$genotype=="N2" & design.tokeep$treatment=="treated")] = "background"
    mod = model.matrix(~ as.factor(tissue.cell), data = design.tokeep);
    #conds = data.frame(rep(c("untreated", "treated"), ncol(cpm)/2))
    #colnames(conds) = 'treatment'
    #mod = model.matrix(~ as.factor(treatment), conds)
    logcpm.bc = ComBat(dat=logcpm, batch=batch, mod=mod, par.prior=TRUE, ref.batch = NULL)
  }
  
  sds = apply(cpm.bc, 1, sd)
  ranks = order(-sds)
  jj = which(ranks<=2000)
  
  pca = prcomp(t((cpm.bc[jj,])), scale. = TRUE, center = TRUE)
  pca2save = data.frame(pca$x, condition=design.matrix$factor.condition, 
                        batch = design.matrix$batch, 
                        name=colnames(cpm))
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = batch, color = condition)) + geom_point(size=4) +
    geom_text(hjust = 0.1, nudge_y = 0.6, size=5) +
    ggtitle(paste0("PCA - "))
  plot(ggp);
  
  #pairs(cpm.bc[, order(design.matrix$batch)], lower.panel=NULL, upper.panel=panel.fitting)
  return(cpm.bc)
  
}


########################################################
########################################################
# Section : customized peak-to-gene assignment
# 
########################################################
########################################################
customizedAssignment.peak2gene = function(peak.coord, window.size=2000, annotation='wormbase')
{
  # a customized way of assigning peak to genes
  # Inputs: peak.coord -- data.frame for peak coordinates with at least 3 first colunms (chr, start, end)
  # 
  
  ##########################################
  # prepare the annotation
  ##########################################
  ### import Refseq annotation
  #annot = read.table('/Volumes/groups/cochella/jiwang/Projects/Julien/R4224/Refseq_annotation.txt', header = TRUE, sep='\t')
  #save(annot, file='Refseq_annotation.Rdata')
  if(annotation=='wormbase'){
    #annot = read.delim('/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.txt', sep='\t', header = TRUE)
    #colnames(annot)[c(3,)]
    #save(annot, file='/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata')
    load(file='/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata')
    ## filter the annotation
    sels = which(annot$Chromosome.scaffold.name != "MtDNA" & annot$Status..gene.=="KNOWN");
    annot = annot[sels, ]
    mm = match(unique(annot$Gene.stable.ID), annot$Gene.stable.ID)
    annot = annot[mm, ]
    aa = data.frame(annot$Chromosome.scaffold.name,
                    annot$Gene.Start..bp., annot$Gene.End..bp., 
                    annot$Strand, annot$Gene.stable.ID, annot$Gene.name, stringsAsFactors = FALSE);
    colnames(aa) = c('chr', 'start', 'end', 'strand', 'wormbase.ID', 'gene')
    aa$strand[which(aa$strand>0)] = '+';
    aa$strand[which(aa$strand<0)] = '-';
    aa$chr = paste0('chr', aa$chr)
  }else{
    if(annotation=='refseq'){
      load(file='Refseq_annotation.Rdata');
      aa = data.frame(annot$ce10.refGene.chrom, 
                      annot$ce10.refGene.txStart, annot$ce10.refGene.cdsEnd, 
                      annot$ce10.refGene.strand, annot$ce10.refGene.name, annot$ce10.refGene.name2, stringsAsFactors = FALSE);
      colnames(aa) = c('chr', 'start', 'end', 'strand', 'RefID', 'gene')
    }else{
      cat('Error :', annotation, ' not available')  
    }
  }
  aa = makeGRangesFromDataFrame(aa, keep.extra.columns=TRUE)
  
  ##########################################
  # assign peaks one by one
  ##########################################
  #load(file= paste0('/Volumes/groups/cochella/jiwang/Projects/Thomas/Thomas_ChIPseq/results/Rdata/Merged_Peaks_macs2_p_5_filtered_N2.Rdata'))
  #load(file='Rdata/Merged_Peaks_macs2_p_5_filtered_N2.Rdata')
  #pp = mergedpeaks;
  peaks = peak.coord
  
  chrom.size = read.delim('/Volumes/groups/cochella/jiwang/annotations/ce10_chrom.sizes', sep='\t', header = FALSE)
  
  peaks.assignment.window = function(px, window.size = 2000){
    #library("tictoc")
    #tic()
    ### px is a vector of chr, star, end ###
    # px = peaks[12282, ]; window.size=2000;
    px = data.frame(px)
    ll = chrom.size[which(chrom.size[,1]==px[,1]), 2]
    aa.sel = aa[seqnames(aa) == as.character(px[, 1])]
    wsize = window.size;
    find = 1
    
    while(find>0)
    {
      start = as.numeric(as.character(px[,2])) - wsize;
      end = as.numeric(as.character(px[,3])) + wsize;
      
      if(start<0) start=0;
      if(end>ll) end=ll;
      
      dd = data.frame(px[,1], start, end, stringsAsFactors = FALSE); 
      colnames(dd) = c('chr', 'start', 'end')
      newp = makeGRangesFromDataFrame(dd)
      #tic()
      #tt = aa.sel[overlapsAny(aa.sel, newp)]
      #toc()
      #tt = findOverlaps(aa.sel, newp, ignore.strand = TRUE)
      #tic()
      tt = subsetByOverlaps(aa.sel, newp, ignore.strand = TRUE)
      #toc()
      
      if(length(tt)>0){
        tt = data.frame(tt)
        targets = c(wsize, t(apply(tt, 2, function(x) gsub(' ', '',paste0(x, collapse = ",")))))
        find = -1;
      }else{
        wsize = wsize + window.size;
      }
      
    }
    #toc()
    return(targets)
  }
  
  targets = data.frame(matrix(NA, nrow = nrow(peaks), ncol = 8))
  colnames(targets) = c('window.size','chr.gene', 'start.gene', 'end.gene', 'width.gene', 'strand.gene', 
                     'WormBase.Gene.ID', 'gene')
  for(n in 1:nrow(peaks)){
    if(n%%1000 == 0) cat('n = ', n, "\n")
    if(as.character(peaks[n, 1]) != 'chrM') targets[n, ] = peaks.assignment.window(peaks[n, ], window.size = window.size)
  }
  
  return(targets)
}

Peak.GPS = function(pp)
{
  #pp = res[407, ]
  summit = floor(mean(as.numeric(pp[match(c('start.peak', 'end.peak'), names(pp))])));
  tss = as.numeric(pp[which(names(pp)=='TSS')]);
  strand = pp[which(names(pp)=='strand.gene')];
  start = as.numeric(pp[which(names(pp)=='start.gene')]);
  end = as.numeric(pp[which(names(pp)=='end.gene')]);
  if(is.na(tss)){
    if(strand=='+'){tss = start;
    }else{tss = end;}
  }
  dist.tss = summit - tss; if(strand =='-') dist.tss = - dist.tss;
  promoter.1kb = (summit< (tss+1000) & summit> (tss-1000))
  promoter.2kb = ((summit< (tss+2000) & summit > (tss+1000)) | (summit> (tss-2000) & summit<(tss-1000)))
  gene.body = (summit< end & summit> start)
  if(strand=='+')
  {
    dowstream.3kb = (summit > end & summit < (end+3000))
  }else{
    dowstream.3kb = (summit > (start -3000) & summit < (start))
  }
  test = c(promoter.1kb, promoter.2kb, gene.body, dowstream.3kb)
  if(sum(test)==0) {test = c(test, TRUE) 
  }else {test = c(test, FALSE)}
  
  return(c(dist.tss, test))
}

##########################################
# function to group peak signals 
##########################################
merge.biologicalReplicates = function(pks)
{
  sampleNames = colnames(pks)
  sampleNames = sapply(sampleNames, function(x) paste(unlist(strsplit(as.character(x), "_"))[1:2], collapse = "_"))
  names.uniq = unique(sampleNames)
  
  merged = matrix(NA, ncol = length(names.uniq), nrow = nrow(pks))
  rownames(merged) = rownames(pks)
  colnames(merged) = names.uniq
  for(n in 1:length(names.uniq))
  {
    kk = which(sampleNames == names.uniq[n])
    if(length(kk)==1){
      merged[,n] = as.numeric(pks[,kk])
    }else{
      merged[,n] = as.numeric(apply(pks[,kk], 1, mean))
    }
  }
  
  return(merged)
  
}

clustering.peak.signals = function(pks, design.matrix, sd.cutoff = 0.7, scale.data = TRUE, cor.cutoff = 0.6, plot.grouping.result = TRUE, 
                                   preprocessing.limma = TRUE, batch.correction = FALSE)
{
  ##########################################
  # this function is to specifically discover patterns for atac-seq peaks in lineage and time points
  ##########################################
  library("cluster")
  library("factoextra")
  library("magrittr")
  library("pheatmap")
  library("RColorBrewer")
  require('limma')
  
  grps.all = rep(NA, nrow(pks))
  names(grps.all) = rownames(pks)
  lsy6.peak = "chrV_10647106_10647681"
  
  ## batch correction with limma (does not change too much)
  if(batch.correction){ 
    #sds = apply(pks, 1, sd)
    bc = remove.batch.atacseq(pks, design.matrix, method = 'limma')
    pks = bc
  }
  pks[which(rownames(pks)== lsy6.peak), ]
  
  ## step 1: filter the peaks showing no changes for time points or lineages 
  if(preprocessing.limma){
    library(splines)
    times = as.numeric(gsub('min', '', design.matrix$condition))
    X <- ns(times, df=2)
    #Then fit separate curves for the control and treatment groups:
    Group <- factor(design.matrix$factor)
    design <- model.matrix(~Group*X)
    fit <- lmFit(pks, design)
    fit <- eBayes(fit)
    
    # define contrast and make selection
    
  }else{
    x = data.matrix(pks);
    
    x[which(rownames(x) == lsy6.peak), ]
    
    means = apply(x, 1, median)
    sds = apply(x, 1, sd)
    
    plot(means, sds, cex = 0.2); abline(h=sd.cutoff, col = 'red', lwd=1.5)
    
    jj = which(sds < sd.cutoff) # regions with low fold changes and low variance were considered no-changing
  }
  
  grps.all[jj] = -1
  
  ## step 2 filter peaks showing changes across time points but with the same pattern for two lineages 
  x = x[-jj, ]
  
  calculate.dist.for.ABa.ABp = function(x)
  {
    return(mean(abs(x[c(1:3)] - x[c(4:6)])))
  }
  
  dist.abs = apply(x, 1, calculate.dist.for.ABa.ABp)
  jj = which(dist.abs<0.25)
  mm = match(names(dist.abs[jj]), names(grps.all))
  grps.all[mm] = 0
  
  ## step 3: cluster the rest of peaks
  my_data = x[-jj, ]
  #my_data = t(apply(my_data, 1, function(x) x/max(x, na.rm = TRUE)))
  #my_data <- x[-jj, ] %>%
  #  na.omit() %>%          # Remove missing values (NA)
  #  scale()                # Scale variables
  
  # View the firt 3 rows
  head(my_data, n = 3)
  my_data[which(rownames(my_data)==lsy6.peak), ]
  my_data  = my_data - 4.0 # background signal of 4 in log2 scale
  my_data[which(my_data<0)] = 0
  #x = t(apply(x, 1, function(x) x/max(x, na.rm = TRUE)))
  my_data[which(rownames(my_data)==lsy6.peak), ]
  
  Kmeans.gasStat = TRUE
  if(Kmeans.gasStat){
    #library("factoextra")
    fviz_nbclust(my_data, kmeans, method = "gap_stat", k.max = 20,  diss = dist(x, method = "euclidean"), nboot = 100)
    
    set.seed(123)
    km.res <- kmeans(my_data, centers = 10, nstart = 25, iter.max = 50)
    groups = km.res$cluster
    
  }else{
    dissimilarity <- 1 - cor(t(my_data))
    distance <- as.dist(dissimilarity) 
    #res.dist <- get_dist(my_data, stand = TRUE, method = "euclidean")
    cluster.groups <- hclust(distance, method="complete") 
    groups <- cutree(cluster.groups, h = (1-cor.cutoff))
    nb.clusters = length(unique(groups))
    cat("after clustering there are : ", length(unique(groups)), " groups \n")
  }
  
  groups[which(names(groups)=="chrV_10647106_10647681")]
  my_data[which(rownames(my_data)=="chrV_10647106_10647681"), ]
  
  mm = match(names(groups), names(grps.all))
  grps.all[mm] = groups
  grps.all[which(names(grps.all)==lsy6.peak)]
  cat(length(which(grps.all == grps.all[which(names(grps.all)==lsy6.peak)])), 'of peaks in the same cluster of lsy-6 peak \n')
  
  if(plot.grouping.result){
    #groups.sels = groups
    nb.peaks.per.groups = table(groups)
    groups.sels = which(nb.peaks.per.groups > 30)
    newgroups = groups[!is.na(match(groups, groups.sels))]
    #o1.groups = order(newgroups)
    newgroups = newgroups[order(newgroups)]
    
    sels = match(names(newgroups), rownames(my_data))
    
    colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")) )(255)
    pheatmap(my_data[sels,],
           col = colors, cluster_rows = FALSE,
           cluster_cols = FALSE, show_rownames = FALSE,
           scale = 'row',
           main = paste0('lsy-6 peak in cluster ', 
                         groups[which(names(groups)=="chrV_10647106_10647681")], 
                         ' out of ', length(unique(groups)), ' clusters')
           )
    
    group2check = grps.all[which(names(grps.all)==lsy6.peak)]
    kk = which(grps.all==group2check)
    
    df = pks[kk, ]
    df = t(apply(df, 1, function(x) x/max(x, na.rm = TRUE)))
    df = data.frame(df)
    boxplot(df)
    #matplot(t(df))
    
    # p <- ggplot(df, aes(x=dose, y=len)) + 
    #     geom_dotplot(binaxis='y', stackdir='center')
    # p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
    #     geom="errorbar", color="red", width=0.2) +
    #     stat_summary(fun.y=mean, geom="point", color="red")
    
  }
  
  #sampleDistMatrix <- as.matrix( sampleDists)
  #rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
  #colnames(sampleDistMatrix) <- NULL
  
 
  
  # if(plot.grouping.result){
  #   #cutoff = 0.8;
  #   plot(cluster.groups, main="Dissimilarity = 1 - Correlation", xlab="") 
  #   rect.hclust(cluster.groups, k = length(unique(groups)), border="red") 
  # }
  
  # regroup proportion matrix
  # newdata = matrix(NA, ncol = nb.clusters, nrow = nrow(x))
  # rownames(newdata) = rownames(x)
  # colnames(newdata) = rep('X', ncol(newdata))
  # 
  # for(n in 1:nb.clusters)
  # {
  #   # n = 4
  #   kk = which(groups == n)
  #   cat("cluster", n, "-- nb of neuron groups --", length(kk), "--", names(groups)[kk],"\n")
  #   
  #   if(length(kk)==0){
  #     cat("Error-- no neuron classes found \n")
  #   }else{
  #     if(length(kk)>1){
  #       newdata[,n] = apply(x[, kk], 1, sum)
  #       colnames(newdata)[n] = paste0(names(groups)[kk], collapse = ".")
  #     }else{
  #       newdata[,n] = x[, kk]
  #       colnames(newdata)[n] = names(groups)[kk];
  #     }
  #   }
  # }
  
  return(grps.all)
  
}



