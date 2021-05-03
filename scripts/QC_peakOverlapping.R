##################################################
## Project: general purpused quality control specific for ChIP-seq data by comparing the peak overlapping between
## biological replicates (for the moment)
## Script purpose: 
## Usage example: Rscript ~/scripts/Chipseq/chipseqQC.R XXX (peak files)
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Jan 10 11:17:03 2018
## tested in cluster
##################################################
args = commandArgs(trailingOnly=TRUE)

# args is a list of files from peak calling
# test if there is at least one argument: if not, return an error
if (length(args)==0){
  cat("this R script is to compare peak overlapping for peak files \n")
  cat("a list of peak files as input required \n")
  cat("Example: \n")
  cat("Rscript ~/scripts/Chipseq/QC_peakOverlapping.R `ls Peaks/macs2/*macs2_pval_0.00001_peaks.xls` \n")
  cat("-----\n")
  stop("a list of peak files as input Missing", call.=FALSE)
  
  }else{
    library("ChIPseeker");
    library("rtracklayer")
    #library(UpSetR);library("ChIPpeakAnno")
    library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
    library("ggplot2");
    library("GenomicFeatures");
    library("GenomicRanges")
    #library("DiffBind");library("IRanges");
    source('/home/imp/jingkui.wang/scripts/Chipseq/functions_chipSeq.R')
    
    DIR.cwd = getwd();
    outDir = paste0(DIR.cwd, "/QCs/peakOverlapping")
    if(!dir.exists(outDir)) dir.create(outDir);
    
    cat("current directory --", DIR.cwd, "\n");
    cat("QC summary directory --", outDir, "\n");
    
    peak.files = args;
    #cat(args, "\n");
    #stop()

    ## parseing the design matrix
    cat("-- parsing design matrix\n")
    bname = basename(peak.files)
    xx = c()
    yy = c()
    for(n in 1:length(bname))
      {
        test = bname[n];
        test = unlist(strsplit(as.character(test), "_"))
        xx = c(xx, test[1])
        #test = unlist(strsplit(as.character(test), "[.]"))[-1]
        #test = paste0(test, collapse = ".")
        yy = c(yy, test[2])
      }
    design = data.frame(yy, xx)

    #design = data.frame(sapply(bname, find.samples.conditions, ID='conditions'), sapply(bname, find.samples.conditions, ID='samples'),
    #                    stringsAsFactors = FALSE)
    design.matrix = data.frame(peak.files, bname, design, stringsAsFactors = FALSE)
    colnames(design.matrix) = c('file.path', 'file.name', 'condition', 'factor')
    design.matrix = design.matrix[with(design.matrix, order(factor, condition)), ];
    #condition= sapply(bname, find.sample.names)
    factor.condition = paste0(design.matrix$factor, '_', design.matrix$condition)
    design.matrix = data.frame(design.matrix, factor.condition, stringsAsFactors = FALSE)

    #head(design.matrix)
    #design.matrix
    #print(design.matrix, "\n");
    #stop();
    ## import macs peaks as Genomic Range object
    cat("-- import peak files as GRanges objects \n")
    peaks.list = c()
    for(n in 1:nrow(design.matrix)){
      #cat(n, '\n')
      xx = readPeakFile(design.matrix$file.path[n], as = "GRanges");
      peaks.list = c(peaks.list, xx)
      #eval(parse(text = paste0("pp.", n, " = xx")))
    }
    
    cat("-- compare peak overlapping and make plots \n")
    pdf(paste0(outDir, "/Comparison_peaks_for_replicates_overlapping.pdf"),
        width = 12, height = 8)
    #source("functions_chipSeq.R") #
    #par(mfrow=c(1,2))
    Comparison.overlapping.peaks(design.matrix, peaks.list, toCompare="factor.condition")

    dev.off()
    
  }                                         

