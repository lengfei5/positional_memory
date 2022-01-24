##########################################################################
##########################################################################
# Project: positional memory 
# Script purpose: process and analyze Akane's histone modification of Cut & Tag
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jan 20 10:43:21 2022
##########################################################################
##########################################################################
rm(list=ls())

version.analysis = 'R11876_R12810_CT'
#peakDir = "Peaks/macs2_broad"

resDir = paste0("../results/", version.analysis)

RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/'

source('Functions_atac.R')
require(ggplot2)
require(GenomicFeatures)


########################################################
########################################################
# Section : sequencing quality controls
# fragment size distribution 
# sequence saturation analysis
########################################################
########################################################
design = read.table(file = paste0(dataDir, 'R11876_cut.run/sampleInfo_parsed.txt'), header = TRUE)
design = rbind(design, read.table(paste0(dataDir, 'R12810_cuttag/sampleInfos_conditions_parsed.txt'), header = TRUE))

stats = read.table(file = paste0(dataDir, 'R11876_cut.run/nf_out/result/countStatTable.txt'), header = TRUE)
stats = rbind(stats, read.table(file = paste0(dataDir, 'R12810_cuttag/nf_out/result/countStatTable.txt'), header = TRUE))

colnames(stats)[c(1, 3)] = c('fileName', 'trimmed')

index = c()
for(n in 1:nrow(design))
{
  # n = 1;
  cat(n, '\n')
  #cc.files = cnts[grep(design$sampleID[n], cnts)]
  ii = grep(design$sampleID[n], stats$fileName)
  if(length(ii) == 1){
    index = c(index, ii)
  }else{
    index = c(index, NA)
    cat(length(ii), 'fileName Found for ', design$sampleID[n], '\n')
  }
  
}

stats = data.frame(design, stats[index, ], stringsAsFactors = FALSE)

colnames(stats) = c('sampleID', 'samples', 'fileName', 'total',  'adapter.trimmed', 'mapped', 'chrM.rm', 'unique', 'unique.rmdup')

stats$trimming.pct = as.numeric(stats$adapter.trimmed)/as.numeric(stats$total)
stats$mapped.pct =  as.numeric(stats$mapped)/as.numeric(stats$adapter.trimmed)
stats$mito.pct = 1.0 - as.numeric(stats$chrM.rm)/as.numeric(stats$adapter.trimmed)
stats$multimapper.pct = 1- as.numeric(stats$unique) / as.numeric(stats$mapped)
stats$dup.rate = 1.0 - as.numeric(stats$unique.rmdup)/as.numeric(stats$unique)
stats$pct.usable = stats$unique.rmdup / stats$total

stats$usable = stats$unique.rmdup/10^6

sels = match(c('163648', '163651', '163656', '163659', 
               '182553', '182554', '182555', '182556'), design$sampleID)

write.csv(stats, file = paste0(resDir, '/', version.analysis,  '_QCs_stats.csv'), row.names = FALSE)



##########################################
# saturation analysis 
##########################################
source('Functions_utility.R')
Sequence.Saturation.Analysis(design)





