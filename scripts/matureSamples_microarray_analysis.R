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
#library(oligoClasses)
#library(oligo)
#library(arrayQualityMetrics)
require(limma)
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
#design = readTargets(file = design.file)
design = design[, c(1:7)]
design = design[which(design$`Mat/Bl` == 'mat'), ]
design = design[, c(1, 7, 2:6)]
colnames(design)[2] = 'condition'
design$condition = gsub('Mat ', 'm', design$condition)
design$condition = gsub('mW', 'mHand', design$condition) # cleaned sample info and condition design
#design$sample = paste0(dataDir, '/', design$sample, '.txt')
#design = design[, c(1, 2)]
#colnames(design)[1] = 'FileName'

# probe intensity
#all <- read.maimages(files = design, source = 'agilent', )
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
      raw = data.frame(featureNum = tmp$FeatureNum, probeName = tmp$ProbeName, tmp$gProcessedSignal, 
                       tmp$gBGMeanSignal, stringsAsFactors = FALSE)
    }else{
      raw = data.frame(raw, tmp$gProcessedSignal[match(raw$featureNum, tmp$FeatureNum)], 
                       tmp$gBGMeanSignal[match(raw$featureNum, tmp$FeatureNum)], stringsAsFactors = FALSE)
    }
    
  }else{
    cat(n, ' : ', design$sample[n], ' ', design$condition[n], ' Error NOT FOUDN \n')
  }
  
}

## clean the probe intensity table
raw = data.frame(raw, stringsAsFactors = FALSE)

colnames(raw)[seq(3, ncol(raw), by=2)] = paste0(design$sample, '_', design$condition) 
colnames(raw)[seq(4, ncol(raw), by=2)] = paste0(design$sample, '_', design$condition, '_BG') 
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

kk = match(mapping$geneID, annot.genes$geneID)
mapping$geneSymbol = annot.genes$gene.symbol.toUse[kk]

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

mat = as.matrix(raw[, c(8:ncol(raw))])
mat = log2(mat)
bg = mat[, grep('_BG', colnames(mat))]
mat = mat[, grep('_BG', colnames(mat), invert = TRUE)]

# background distribution; and keep the probes with >=2 probes above background with pvalue < 0.05
par(mfrow = c(1, 2))
hist(bg[, 1], breaks = 100)
hist(bg[, 2], breaks = 100)

Correct.Background.limma = FALSE
if(Correct.Background.limma){
  mat.bc = limma::backgroundCorrect.matrix(E = 2^mat, Eb = 2^bg, method = 'normexp', offset = 0, normexp.method = 'rma')
  mat.bc = log2(mat.bc)
  bg.cutoff = quantile(log2(tmp$gBGMeanSignal), 0.95)
}

Correct.Background.manual = TRUE
mat.bc = mat
if(Correct.Background.manual){

  ss = apply(bg, 2, mean)
  ss = ss - mean(ss)
  for(n in 1:ncol(mat.bc)) mat.bc[,n] = mat.bc[,n] - ss[n]
}

bg.cutoff = quantile(bg, prob = 0.95)
nb.pass.cutoff = apply(mat.bc, 1, function(x) length(which(x > bg.cutoff)))

sels = which(nb.pass.cutoff >= 3)

mat.bc = mat.bc[sels, ]
raw = raw[sels, ]

## change the colomun order to put replicates together
mat.bc = mat.bc[, c(grep('UA', colnames(mat.bc)), grep('LA', colnames(mat.bc)), grep('Hand', colnames(mat.bc)))]
colnames(mat.bc) = paste0(rep(c('mUA', 'mLA', 'mHand'), each = 3), '_', c(1:3))

raw = data.frame(raw[, c(1:7)], mat.bc) 
#colnames(raw)[c(8:16)] = colnames(mat)

save(design, raw, mapping, file = paste0(RdataDir, 'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_BgCorrected.Rdata'))

##########################################
# Qunatile normalization
##########################################
load(file = paste0(RdataDir, 'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_BgCorrected.Rdata'))

library(preprocessCore)
mat = as.matrix(raw[, c(8:16)])

mat.norm = normalize.quantiles(mat)
#mat.norm = limma::normalizeBetweenArrays(mat, )
colnames(mat.norm) = colnames(mat)

mat.norm[grep('HOXA13', raw$geneSymbol),]
mat.norm[grep('RARRES1', raw$geneSymbol),]

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


##########################################
# position-dependent genes from microarray
# original code were from limma manual section 9.3
##########################################
require(limma)
load(file = paste0(RdataDir, 'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary.Rdata'))
annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

mm = match(rownames(res), annot$geneID)
ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
rownames(res)[!is.na(mm)] = ggs[!is.na(mm)]

f <- factor(rep(c('mUA', 'mLA', 'mHand'), each = 3), levels=c("mUA","mLA","mHand")) 
design <- model.matrix(~0+f)
colnames(design) <- c("mUA","mLA","mHand")

#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(res, design)
contrast.matrix <- makeContrasts(mLA-mUA, mHand-mUA, mHand-mLA, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#A list of top genes for RNA2 versus RNA1 can be obtained from
topTable(fit2, coef=1, adjust="BH")

#The outcome of each hypothesis test can be assigned using
results <- decideTests(fit2)

xx = data.frame(fit2$p.value)
colnames(xx) = c('mLA.vs.mUA', 'mHand.vs.mUA', 'mHand.vs.mLA')
xx$pval.max = apply(as.matrix(-log10(xx)), 1, max)

res = data.frame(res, xx)
res = res[order(-res$pval.max), ]

save(res, raw, fit2, file = paste0(RdataDir, 
                             'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_DEpval.Rdata'))

