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
library(Biobase)
library(oligoClasses)
library(oligo)
library(arrayQualityMetrics)

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

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
annot.genes = readRDS(paste0(annotDir, 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
annot.transcript = readRDS(file = paste0(annotDir, 
                'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))

mm = match(mapping$transcript, annot.transcript$transcriptID)
mapping$transcript.strand = annot.transcript$strand_transcript[mm]
mapping$geneID = annot.transcript$geneID[mm]

kk = match(mapping$geneID, annot$geneID)
mapping$geneSymbol = annot$gene.symbol.toUse[kk]

save(design, raw, mapping, file = paste0(RdataDir, 'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol.Rdata'))


########################################################
########################################################
# Section : probe filtering and normalization
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol.Rdata'))

# select only probes mapped to the positive strand of transcripts
cat(nrow(mapping), ' annotated mapping : keep only the positive strand probes \n')
mapping = mapping[which(mapping$strand == '+'), ]
cat(nrow(mapping), ' probes left \n')

# filter unannotated probes
cat(nrow(raw), ' probes in intensity matrix; now keep only those with annotated transcrips \n')
mm = match(raw$probeName, mapping$probeName)
raw = raw[!is.na(mm), ]
cat(nrow(raw), ' probes left in intensity matrix \n')

# combine the intensity matrix with transcript and genes
mm = match(raw$probeName, mapping$probeName)

raw = data.frame(mapping[mm, c(1,4, 5:6, 8:9)], raw, stringsAsFactors = FALSE)
raw = raw[, -8] 

mat = as.matrix(raw[, c(8:16)])
mat = log2(mat)

# background distribution; and keep the probes with >=2 probes above background with pvalue < 0.05
hist(log2(tmp$gBGMeanSignal), breaks = 100)
bg.cutoff = quantile(log2(tmp$gBGMeanSignal), 0.95)

nb.pass.cutoff = apply(mat, 1, function(x) length(which(x > bg.cutoff)))

sels = which(nb.pass.cutoff >= 2)

mat = mat[sels, ]
raw = raw[sels, ]

mat = mat[, c(grep('UA', colnames(mat)), grep('LA', colnames(mat)), grep('Hand', colnames(mat)))]

raw[, c(8:16)] = mat
colnames(raw)[c(8:16)] = colnames(mat)

##########################################
# Qunatile normalization
##########################################
library(preprocessCore)
mat = as.matrix(raw[, c(8:16)])

mat.norm = normalize.quantiles(mat)
colnames(mat.norm) = paste0(rep(c('mUA', 'mLA', 'mHand'), each = 3), '_', c(1:3))

ggs = unique(raw$geneID)
res = matrix(NA, ncol = ncol(mat.norm), nrow = length(ggs))
colnames(res) = colnames(mat.norm)
rownames(res) = ggs

for(n in 1:nrow(res))
{
  if(n%%500 == 0) cat(n, '\n')
  jj = which(raw$geneID == rownames(res)[n])
  if(length(jj) > 1) {
    res[n, ] = apply(mat.norm[jj, ], 2, median)
  }else{
    res[n, ] = mat.norm[jj, ]
  }
}

         save(res, raw, file = paste0(RdataDir, 'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary.Rdata'))

