#############################
##################################################
## Project: ZFP445 binding 
## Script purpose: to correlat ZFP445 binding peaks with H3K27me3 and H3K9me3
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Jan 23 11:56:33 2018
##################################################
library(Rsubread)
library(csaw)
#library(edgeR)

##################################################
## Section: specify input, output and parameters
##################################################
DIR.bam = "Bams_4comparison";

DIR.res = "result_compare_ZFP462_otherChIPs"
if(!dir.exists(DIR.res)) system(paste0('mkdir -p ', DIR.res))

bam.list = list.files(path=DIR.bam, pattern = "*.bam$", full.names = TRUE, ignore.case = TRUE)
#bam.list = bam.list[grep("Zfp|309M|H3K27", bam.list)]

Count.Reads.genome.wide.CSAW = FALSE

####################
## quantify signals for the whole geneome and estimate normalization factors 
####################
if(Count.Reads.genome.wide.CSAW) # not used here
{
  ## quantify read counts from bam files
  frag.len <- 110
  win.width <- 500
  #chr.selected = c(paste0("chr", c(1:20)))
  param <- readParam(pe="none", dedup = TRUE, minq=30)
  #chr.selected = c(paste0("chr", c(1:20)))
  #param <- readParam(pe="none", dedup = TRUE, minq=30, restrict = chr.selected)
  data <- windowCounts(bam.list, ext=frag.len, width=win.width, param=param, bin = TRUE)
  
  save(data, bam.list, file = paste0(DIR.res, "/csaw_windowCounts_output.Rdata"))
  
  ## filtering out uninteresting regions
  load(file = paste0(DIR.res, "/csaw_windowCounts_output.Rdata"))
  library(edgeR)
  abundances = aveLogCPM(asDGEList(data))
  summary(abundances)
  
  # check the threshold for the logcpm
  hist(abundances, breaks = 100); abline(v=seq(-3, 2, by=1), col='darkred')
  length(which(abundances>-2))
  
  keep <- abundances >= 0
  data <- data[keep,]
  
  ## calculate normalization factors 
  #binned <- windowCounts(bam.files, bin=TRUE, width=50000, param=param)
  aa <- normOffsets(data, type="scaling", se.out=TRUE)
  cpm = cpm(asDGEList(aa), log=TRUE)
  
  res = cpm;
  
  ## table processing
  #coordiantes = as.data.frame(rowRanges(data)) 
  #ggs = apply(coordiantes[, c(1:3)], 1, function(x) gsub(" ", "",paste(as.character(unlist(x)), sep = "",  collapse = "_"), fixed = TRUE))
  #ggs = paste0(coordiantes[, c(1:3)], collapse = "_") 
  
  #aa = data.frame(assay(data), stringsAsFactors = FALSE)
  #colnames(aa) = basename(filtered.data$bam.files)
  
}else{ ## quantify signals in the pre-selected regions
  Quantify.signals.within.peaks = TRUE
  Count.Reads.for.Peak.Regions = TRUE
  Change.Sample.Names = FALSE
  Count.Reads.genome.wide.CSAW = FALSE
  Merge.bed.file = TRUE
  
  DIR.bed = "result_compare_ZFP462_otherChIPs";
  # import chipseq peqks 
  bed.list = list.files(path=DIR.bed, pattern = "*.bed$", full.names = TRUE, ignore.case = TRUE)
  
  
  if(Merge.bed.file){
    peaks = c();
    for(n in 1:length(bed.list))
    {
      bed = read.table(bed.list[n], sep = "\t", header = FALSE)
      
      if(length(unique(bed[,4])) != nrow(bed)) {
        cat("peak names are Redundant ! --- rename the peak namaes\n");
        newPeakNames = apply(bed[, c(1:3)], 1, function(x) gsub(" ", "", paste0(unlist(x), collapse = "_")))
        bed[,4] = newPeakNames;
      }
      if(length(grep('random', basename(bed.list[n])))>0) bed[, 4] = paste0("random_", seq(1:nrow(bed)));
      peaks = rbind(peaks, bed); 
    }
    peaks = data.frame(peaks, stringsAsFactors = FALSE)
    colnames(peaks) = c("chr.peak", "start.peak", "end.peak", "peak.name", "score.peak", "strand.peak")[c(1:ncol(peaks))]
  }
  
  ######################################################
  ## Quantify chipseq signals for all peaks regions
  ######################################################
  if(Quantify.signals.within.peaks)
  {
    if(Count.Reads.for.Peak.Regions)
    {
      ## peak regions (configure your peak regions into a data.frame)
      jj = match(unique(peaks$peak.name), peaks$peak.name)
      df = peaks[jj, ];
      SAF = data.frame(GeneID=df$peak.name, Chr=df$chr.peak, Start=df$start.peak, End=df$end.peak, Strand="+", stringsAsFactors = FALSE)
      
      ## count reads for those peak regions using 
      fc <- featureCounts(files=bam.list, annot.ext = SAF, countMultiMappingReads = FALSE, minMQS = 10, 
                          ignoreDup = TRUE, strandSpecific = 0, juncCounts = FALSE, nthreads = 20)
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
      
      ## change colnames
      if(Change.Sample.Names)
      {
        names = colnames(res)[-c(1:7)]
        for(n in 1:length(names))
        {
          #n = 1
          test = unlist(strsplit(as.character(names[n]), "[.]"))
          test = test[3]
          if(length(grep("Input", test))>0){
            ttest = unlist(strsplit(as.character(test), "_"))
            names[n] = paste0(ttest[c(1,3, 2)], collapse = "_")
          }
          if(length(grep("_H3K", test))>0){
            ttest = unlist(strsplit(as.character(test), "_"))
            names[n] = paste0(ttest[c(1, 4, 2)], collapse = "_")
          }
          if(length(grep("tbx3", test))>0) names[n] = test
        }
        colnames(res)[-c(1:7)] = names
        
        kk = grep("6291|6344|Input", colnames(res)) ## all inputs for histone modifications are the same one.
        if(length(kk)>0) res = res[, -kk]
      }
      
      
      ## save the peak information and quantified different rpkm signals within peaks
      write.table(res, file=paste0(DIR.res, "/rpkm_within_chipeakk.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      #save(peaks, res, file='Rdata/Peaks_merged_macs2_p_5_filtered_N2_gene_assignment_TSS_WBcel235_analysis_hisModif_tbx_signals_v1.Rdata')
    }
  }
}

###############################
# test correlation between ZFP462 and other chipseq data 
###############################
res = read.table(file=paste0(DIR.res, "/rpkm_within_chipeakk.txt"), sep = "\t", header = TRUE)

res = res[, -c(1:6)]
colnames(res) = basename(bam.list)
colnames(res) = sapply(colnames(res), function(x) paste0(unlist(strsplit(as.character(x), "_"))[1:2], collapse = "_"))
res = as.matrix(res)

jj = which(res[,19]>0.5 & res[,20]>0.5)
res =res[jj, ]

res = log2(res+2^-6)

corres = c()
pdfname = paste0(DIR.res, "/comparison_ZFP462_vs_others_withPeaks_rpkm.pdf")
pdf(pdfname, width = 12, height = 10)

cc = cbind(rep(20, 19), c(1:19))
cc = rbind(cc, c(5, 11), c(17,15), c(17, 16), 
           c(10, 11), c(10, 12), c(10, 13), c(10,14), c(10, 7))
cex = 0.1

for(n in 1:nrow(cc))
{
  # n = 1
  corres = c(corres, cor(res[ ,cc[n,1]], res[, cc[n, 2]]))
  #cc = cor(res[,20], res[,n])
  plot(res[,cc[n, 1]], res[, cc[n, 2]],
       main=paste0("R = ", signif(cor(res[ ,cc[n,1]], res[, cc[n, 2]]), d=2)), cex=cex,
       xlab = colnames(res)[cc[n, 1]], ylab=colnames(res)[cc[n, 2]], log='')
}

#plot(res[,5], res[,11], main=paste0("R = ", signif(cor(res[,5], res[,11]), d=2)), cex = cex)
#plot(res[,17], res[,15], main=paste0("R = ", signif(cor(res[,17], res[,15]), d=2)), cex =cex)
#plot(res[,17], res[,16], main=paste0("R = ", signif(cor(res[,17], res[,16]), d=2)), cex= cex)
dev.off()

corres = corres[1:19]
names(corres) = colnames(res)[1:19]


######################################################
### test the correlations
######################################################
Search4correlations = FALSE
if(Search4correlations)
{
  dd = as.matrix(res[, -c(1:6)])
  rownames(dd) = res$GeneID
  colnames(dd) = c("H3K27me3.Input", "H3K27me3", "H3K9me3_Input", "H3K9me3_309M3A", "H3K9me3_309M3B", "ZFP445")
  kk.random = grep("random", rownames(dd))
  rr = dd[kk.random,]
  #dd = dd[-kk.random, ]
  colnames(rr) = colnames(dd)
  dd = log2(dd+2^-10)
  rr = log2(rr+2^-10)
  groups = rep(NA, nrow(dd))
  groups[grep('MACS', rownames(dd))] = 1
  groups[grep('MACS', rownames(dd), invert = TRUE)] = 2
  
  pdfname = paste0(DIR.res, "/comparison_plots_zfp445_vs_H3K9_H3K27.pdf")
  pdf(pdfname, width = 12, height = 8)
  
  ## 
  ## boxplot for peaks and randome locations
  library(vioplot)
  par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(12,3,2,0.8)+0.1, tcl = -0.3)
  par(mfrow=c(1,1))
  plot(c(0, 12), c(-10, 6), type = 'n', axes = TRUE, ylab = "log2 (rpkm)", xlab=NA)
  for(n in 1:ncol(dd))
  {
    add = TRUE
    if(n>1) add = TRUE
    #boxplot(dd[,n] ~ groups, col= c("blue", "darkgray"), at=c((2*n-1), 2*n), las=2,
    #        names=paste0(colnames(dd)[n], c("_peak", "_random")), add =add)
    vioplot(dd[which(groups==1),n], names=paste0(colnames(dd)[n], c("_peak")), horizontal = FALSE,
            col=c("blue"), add = add, at=c((2*n-1)))
    mtext(paste0(colnames(dd)[n], "_peaks"), side=1, at = (2*n-1), las=3)
    vioplot(dd[which(groups==2), n], names=paste0(colnames(dd)[n], c("_random")), horizontal = FALSE,
            col=c("darkgray"), add = add, at=c(2*n))
    mtext(paste0(colnames(dd)[n], "_random"), side=1, at = (2*n), las=3)
  }
  
  kk.peaks = which(groups==1)
  kk.random = which(groups==2)
  
  counts = matrix(NA, ncol = 3, nrow = 2)
  colnames(counts) = colnames(dd)[1:3]
  for(n in 1:ncol(counts))
  {
    counts[1,n]
  }
  
  cex = 0.25
  ## scatterplot for all 
  par(pty="s")
  pairs(dd, cex=cex)
  ## scatterplot for peaks
  kk = grep("MACS", rownames(dd))
  pairs(dd[kk,], cex=cex)
  ## scatterplot for randome sequence
  kk = grep("MACS", rownames(dd), invert = TRUE)
  pairs(dd[kk,], cex=cex)
  
  dev.off()
  
}

##################################################
## Section: check the peak overlapping
##################################################
library("ChIPseeker");
library("rtracklayer")
library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library("ggplot2");
library("GenomicFeatures");
library("GenomicRanges")
#library("DiffBind");
#library("IRanges");
source("functions_chipSeq.R") ## these functions are from Thomas Burkard 

DIR.peaks = "Peaks/macs2_broad"
#patterns = "*.bed"
xlist<-list.files(path=DIR.peaks, pattern = "*.bed$", full.names = TRUE)
ylist<-list.files(path=DIR.peaks, pattern = "*.xls$", full.names = TRUE)
xlist = c(xlist, ylist)
#xlist = xlist[grep("GSM8550", xlist, invert = TRUE)]
#xlist = xlist[grep("Input", xlist, invert = TRUE)]
#xlist = xlist[grep("")]
#xlist = xlist[c(1,2, 4)]

peaks.list = c() 
for(n in 1:length(xlist)){
  cat(n, '\n')
  xx = readPeakFile(xlist[n], as = "GRanges"); 
  peaks.list = c(peaks.list, xx)
  #eval(parse(text = paste0("pp.", n, " = xx")))
}

pdf(paste0(DIR.res, "/Comparison_Peaks_overlapping_after_ZFN445_H3K9_Hek_309M3A_B.pdf"), width = 12, height = 8)

source("functions_chipSeq.R") #
#par(mfrow=c(1,2))
peaknames = c("ZNF445", "H3K9_309M3A", "H3K9_309M3B")
sels = c(1:3)

peaks = c()
peaks.extend = c()
for(k in sels) 
{
  p = peaks.list[[k]]
  pe = resize(p, width = 5000, fix = "center")
  #p = readPeakFile(ff$file.path[k], as = "GRanges");
  #eval(parse(text = paste0("p = pp.", k)));
  #p10 <- p[mcols(p)[,"X.log10.pvalue."] > pval];
  p = reduce(p); 
  pe = reduce(pe,  drop.empty.ranges=TRUE, min.gapwidth=3000L)
  #p10 = reduce(p10);
  peaks= c(peaks, p);
  peaks.extend = c(peaks.extend, pe)
  #peaks10= c(peaks10, p10);
}

ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=peaknames[sels], connectedPeaks="keepAll", main="peak overlapping")
v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

ol.peaks <- makeVennDiagram(peaks.extend, NameOfPeaks=peaknames[sels], connectedPeaks="keepAll", main="extended peaks overlapping")
v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))
#Comparison.overlapping.peaks (design.matrix, peaks.list, toCompare="factor.condition") 

dev.off()







