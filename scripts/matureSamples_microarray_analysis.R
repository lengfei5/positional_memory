##########################################################################
##########################################################################
# Project: Akane's positional memory project
# Script purpose: analyze microarray data of mature samples to double comfirm Akane's mature RNA-seq samples
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jul 30 13:36:06 2021
##########################################################################
##########################################################################
rm(list = ls())

source('functions_ExonsArray.R')
require(openxlsx)
require(ggplot2)
require(dplyr)
require(gridExtra)
library(data.table)

version.Data = 'microarray'
version.analysis = paste0("_", version.Data, "_20210730")

## Directories to save results 
design.file = "../exp_design/microarray_sample_overview_PM.xlsx"
#dataDir = "../data/quantseq/featurecounts_R7183"

resDir = paste0("../results/", version.Data)
tabDir =  paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}


########################################################
########################################################
# Section : import and process probe-to-transcript mapping and probe intensity
# 
########################################################
########################################################
cat("Loading sample info file \n")

dataDir = '/Volumes/groups/tanaka/Data/microarrays/Prayag/Proximal.Distal_Dunja.Prayag/ExpressionData'

design = read.xlsx(design.file)
design = design[, c(1:7)]
design = design[which(design$`Mat/Bl` == 'mat'), ]
design = design[, c(1, 7, 2:6)]
colnames(design)[2] = 'condition'
design$condition = gsub('Mat ', 'm', design$condition)
design$condition = gsub('mW', 'mHand', design$condition) # cleaned sample info and condition design

# probe intensity
files = list.files(path = dataDir, full.names = TRUE)

for(n in 1:nrow(design))
{
  # n = 1
  kk = grep(design$sample[n], files)
  
  if(length(kk) == 1){
    cat(n, ' : ', design$sample[n], ' ', design$condition[n], '\n')
    tmp = fread(files[kk], header='auto', stringsAsFactors=FALSE, skip = 9)
    # Sergej's script to read expression data
    # # tmp = read.table(fname, header=T, skip=9);
    # # ma_data = cbind(ma_data, tmp$gProcessedSignal); 
    # ma_data = NULL; # Microarray data
    # ma_ids = NULL;  # Microarray IDs, e.g. <Tissue>_<Replicate>
    # 
    # cat("Loading the data files\n");
    # nFiles = dim(overview)[1];
    # 
    # for (i in 1:nFiles)
    # {
    #   fname=paste(rootdir, '/', 'ExpressionData', '/', overview$Array.name[i], '.txt', sep=""); 
    #   cat(paste("   Loading the file ", i, " of ", nFiles, " ('", fname, "')\n", sep=""));
    #   
    #   tmp = read.table(fname, header=T, skip=9); 
    #   ma_data = cbind(ma_data, tmp$gProcessedSignal); 
    #   ma_ids = c(ma_ids, paste(overview$day.blastema[i], 'dpa', '_', 'R', overview$replicate[i], sep=""));
    # }
    # 
    # # Update the names of the rows and the columns
    # colnames(ma_data) = ma_ids;
    # rownames(ma_data) = tmp$ProbeName;
    if(n == 1){
      raw = data.frame(featureNum = tmp$FeatureNum,  probeName = tmp$ProbeName, tmp$gProcessedSignal, stringsAsFactors = FALSE)
    }else{
      raw = data.frame(raw, tmp$gProcessedSignal[match(raw$featureNum, tmp$FeatureNum)], stringsAsFactors = FALSE)
    }
    
  }else{
    cat(n, ' : ', design$sample[n], ' ', design$condition[n], ' Error NOT FOUDN \n')
  }
  
}

## clean the probe intensity table
raw = data.frame(raw, stringsAsFactors = FALSE)
colnames(raw)[c(3:ncol(raw))] = paste0(design$sample, '_', design$condition) 
#rownames(raw) = raw[, 1]
#raw = raw[, -1]

save(design, raw, file = paste0(RdataDir, 'design_probeIntensityMatrix.Rdata'))

##########################################
# process the probe-to-gene mapping  
##########################################
mapping = fread('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/microarray_probes/AmexMA_v2-AmexT_v47.bed', 
                header=FALSE, stringsAsFactors=FALSE)
colnames(mapping) = c('transcript', 'start', 'end', 'probeName', 'match', 'strand')

# filter the probes with >=5 mismatches
mapping = mapping[which(mapping$match > 55), ]





