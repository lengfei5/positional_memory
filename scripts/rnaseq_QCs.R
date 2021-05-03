##########################################################################
##########################################################################
# Project: Positional memory
# Script purpose: Control quality of RNA-seq data and analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Feb 19 15:45:30 2021
##########################################################################
##########################################################################
version.analysis = 'Rxxxx.rnaseq.old_202102020'

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../Data/Rxxxx_rnaseq_old/'

design = read.table(paste0(dataDir, 'sampleInfos_parsed.txt'), sep = '\t', header = TRUE)

########################################################
########################################################
# Section : sequencing quality controls
# 
########################################################
########################################################
stats = read.delim(paste0(dataDir, 'nf_out_RNAseq/MultiQC/multiqc_data/multiqc_general_stats.txt'), sep = '\t', 
                   header = TRUE)
alignment = read.delim(paste0(dataDir, 'nf_out_RNAseq/MultiQC/multiqc_data/multiqc_hisat2.txt')
                       , sep = '\t')

stats = stats[, c(1, 2, 3, 4, 6, 8, 9, 10)]
colnames(stats) = c('sample', 'pct.duplication', 'pct.GC', 'avg.seq.length', 'total.reads', 
                    'pct.assign', 'assigned.reads', 'alignment.rate')

stats = data.frame(stats, alignment[match(stats$sample, alignment$Sample), c(2, 4, 5)], stringsAsFactors = FALSE)
colnames(stats)[c(9:11)] = c('trimmed.reads', 'unique.aligned', 'multimapper')
stats = stats[, c(1:5, 9, 8, 10, 11, 7)]

#design = design[order(design$fileName), ]
ii = c()
jj = c()
for(n in 1:nrow(design))
{
  # n = 1;
  cat(n, '\n')
  kk = grep(design$sampleID[n], stats$sample)
  ii = c(ii, rep(n, length(kk)))
  jj = c(jj, kk)
  #kk = c(kk, grep(design$sampleID[n], stats$sample))
  #kk = c(kk, which(design$sampleID == ))
}

xx = data.frame(design[ii, ], stats[jj, ], stringsAsFactors = FALSE)
xx = xx[order(xx$fileName), ]

write.csv(xx, file = paste0(resDir, '/QCs_stats.csv'), row.names = FALSE)

design = data.frame(design, stats[kk, ], stringsAsFactors = FALSE)

########################################################
########################################################
# Section : saturation curve from rseqc
# the r code from rseqc output
########################################################
########################################################
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
