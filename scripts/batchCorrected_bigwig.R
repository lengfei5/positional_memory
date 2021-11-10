##########################################################################
##########################################################################
# Project: positional memory project and general
# Script purpose: make batch corrected bigwig file
# Usage example: the script is design for one bigwig file and to run multiple files in paralell
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Nov  9 11:17:19 2021
##########################################################################
##########################################################################
require(GenomicFeatures)
require(rtracklayer)

bgDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/bigwigs_scalingFactor_consensusPeaks'

version.analysis = 'atac_rna_chipseq_analysis_20211007'
resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')

##########################################
# import batch corrected peak files
##########################################
load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical_v1.Rdata'))
fpm = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))

# prepare the background distribution
fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
rownames(fpm) = gsub('_', '-', rownames(fpm))

fpm = 2^fpm

##########################################
## make Granges for peak coordinates
##########################################
pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
pp$strand = '*'
pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")


##########################################
# read the orginial bw file
##########################################
bw.files = list.files(path = bgDir, pattern = '*.bw', full.names = TRUE)

bw = import(bw.files[30], format = 'bw')
bw = dropSeqlevels(bw, 'chrM', pruning.mode=c("coarse"))

# find the gap between peaks and make one Grange 
gp = GenomicRanges::gaps(pp)
pp$score = 1
gp$score = 0

xx  = c(pp, gp)
xx = sort(xx)

newb = bw
xx$scalingfactor = NA 

# loop over the peaks and correct the bins within peaks 
library(tictoc)
tic()
overlaps = findOverlaps(pp, bw)

for(n in 1:2)
{
  # n = 1
  bins = bw[overlaps@to[which(overlaps@from == n)]]
  if(length(bins) > 0){
    sf = fpm[m, n]/sum(width(bw.peak) * bw.peak$score) * 100 # 50 is constant to reflect the fragment size 
    newb$score[kk] = bw$score[kk] * sf 
    
  }
}

toc()

export.bw(xx, con = paste0(OutDir, bw.name))
