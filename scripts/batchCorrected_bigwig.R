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
library(tools)

bw.file = commandArgs(trailingOnly=TRUE)

if (length(bw.file)==0 | length(bw.file) > 1){
  stop("one and only one bigwig file required", call.=FALSE)
}else{
  
  #filename = file_path_sans_ext(bw.file)
  extension = file_ext(bw.file)
  
  if(extension != 'bw') stop('one BIGWIG file required', call. = FALSE)  
}

outDir = './outs/'
if(!dir.exists(outDir)) dir.create(outDir)

# bw.file = list.files(path = bgDir, pattern = '*.bw', full.names = TRUE)[30]

##########################################
# import batch corrected peak files
##########################################
load(file = paste0('samplesDesign.cleaned_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical_v1.Rdata'))
fpm = readRDS(file = paste0('fpm_TMM_combat.rds'))

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
bw.name = basename(bw.file)
bw.name = gsub('_mq_30.bw', '', bw.name)

index.sample = which(colnames(fpm) == bw.name)

if(length(index.sample) != 1){
  print('Sample Not Found')
  stop()
}

cpm = fpm[, index.sample]

bw = import(bw.file, format = 'bw')
bw = dropSeqlevels(bw, 'chrM', pruning.mode=c("coarse"))

# find the gap between peaks and make one Grange 
gp = GenomicRanges::gaps(pp)
pp$score = 1
gp$score = 0

xx  = c(pp, gp)


xx$scalingfactor = NA 

seqlevels(pp) = sort(seqlevels(pp))

seqlevels(bw) = seqlevels(pp)
seqlevels(xx) = seqlevels(pp)
xx = sort(xx)
bw = sort(bw)


newb = bw

##########################################
# loop over the peaks and correct the bins within peaks
##########################################
overlaps = findOverlaps(bw, pp, type = 'within', select = 'all')
mappings = match(names(pp), names(xx))

library(tictoc)
tic()
for(n in 1:length(pp))
{
  # n = 7
  kk = overlaps@from[which(overlaps@to == n)]
  bins = bw[kk]
  
  if(length(bins) > 0){
    
    sf = cpm[which(names(cpm) == names(pp)[n])]/sum(as.numeric(width(bins)) * bins$score) * 50 # 50 is constant to reflect the fragment size
    
    #tic()
    # scale the bigwig scores
    if(!is.na(sf) & sf != Inf){
      mcols(newb)$score[kk] = mcols(bw)$score[kk] * sf
      # save the scaling factor
      xx$scalingfactor[mappings[n]] = sf
    }else{
      print(paste0('NA or Inf scaling factor found for ',  n, '-- peak name : ', names(pp)[n]))
    }
  }else{
    print(paste0('no bins within  ', n, '-- peak name : ', names(pp)[n]))
  }
}

toc()

##########################################
# impute the missing scaling factors from previous steps due to the sf = NA or sf = Inf
# impute them with the median of good sf values in a heuristic way so as to continue with the gap scaling
##########################################
index.missing.sf = mappings[which(is.na(xx$scalingfactor[mappings]) | xx$scalingfactor[mappings] == Inf)]
index.good.sf = setdiff(mappings, index.missing.sf)
xx$scalingfactor[index.missing.sf] = median(xx$scalingfactor[index.good.sf])

##########################################
# # loop over the gaps between peaks, chr by chr 
##########################################
for(chr in seqlevels(xx))
{
  # chr = 'chr10p'
  print(chr)
  xxchr = xx[which(seqnames(xx) == chr)]
  
  overlaps = findOverlaps(xxchr, bw)
  
  tic()
  for(n in 1:length(xxchr))
  #for(n in 1:20)
  {
    # n = length(xxchr) - 1
    # n = 33
    if(mcols(xxchr)$score[n] == 0){
      
      jj = overlaps@to[which(overlaps@from == n)]
      bins = bw[jj]
      
      if(n == 1){ # first gap in chr
        sf_final = mcols(xxchr)$scalingfactor[2]
        if(!is.na(sf_final) & sf_final > 0){
          mcols(newb)$score[jj] = mcols(bw)$score[jj] * sf_final
        }
        
      }else{
        if(n == length(xxchr)){  # last gap in chr
          sf_final = mcols(xxchr)$scalingfactor[(n-1)]
          if(!is.na(sf_final) & sf_final > 0){
            mcols(newb)$score[jj] = mcols(bw)$score[jj] * sf_final
          }
        }else{ # gaps with flanking peaks
          
          # default values 
          sf_prev = 1
          sf_next = 1
          sf_prev = mcols(xxchr)$scalingfactor[(n-1)]
          sf_next = mcols(xxchr)$scalingfactor[(n+1)]
          
          if(!is.na(sf_prev) & !is.na(sf_next)){
            #fractions = ((as.numeric(start(bins)) + as.numeric(end(bins)))/2.0 - start(xxchr)[n])/width(xxchr)[n]
            fractions = ((as.numeric(start(bins)) + as.numeric(end(bins)))/2.0 - as.numeric(start(xxchr)[n]))/as.numeric(width(xxchr)[n])
            
            if(any(fractions<0)) fractions[which(fractions<0)] = 0
            if(any(fractions>1)) fractions[which(fractions>1)] = 1
            if(!all(is.na(fractions))) fractions[is.na(fractions)] = median(fractions, na.rm = TRUE)
            
            sf_final = sf_prev*(1.0-fractions) + sf_next*fractions 
            
            if(all(!is.na(sf_final)) & all(sf_final != Inf)){
              mcols(newb)$score[jj] = mcols(bw)$score[jj] * sf_final
            }else{
              print(paste0('NA scaling factor found ',  chr, '-- index -- ', n,  '-- gap : ', 
                           seqnames(xxchr[n]), ':', start(xxchr[n]), '-', end(xxchr[n])))
            }
            
            
          }
        } # end of gaps with flanking peaks
      }
    }
    
  } # end of loop over gaps in chr
  toc()
  
}

# export the batch-corrected bigwig file
export.bw(newb, con = paste0(outDir, bw.name, '_bc.bw'))

##########################################
# debugging the script run in cluster  
##########################################
Debugging = FALSE
if(Debugging){
  # intermediate results from test example
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/bigwigs_bc/'
  bw.file = '/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/bigwigs_scalingFactor_consensusPeaks/Mature_LA_74939_mq_30.bw'
  
  load(file = paste0(dataDir, 'samplesDesign.cleaned_readCounts.within_manualConsensusPeaks.pval3_mergedTechnical_v1.Rdata'))
  fpm = readRDS(file = paste0(dataDir, 'fpm_TMM_combat.rds'))
  
  # prepare the background distribution
  fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
  fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  fpm = 2^fpm
  
  bw.name = basename(bw.file)
  bw.name = gsub('_mq_30.bw', '', bw.name)
  
  index.sample = which(colnames(fpm) == bw.name)
  
  if(length(index.sample) != 1){
    print('Sample Not Found')
    stop()
  }
  
  cpm = fpm[, index.sample]
  
  
  load(paste0(dataDir, 'outs/Mature_LA_74939_bc_peaks_gaps.Rdata'))
  
  jj = which(is.na(mcols(newb)$score))
  head(which(mcols(bw)$score[jj]>0))
  
  findOverlaps(bw[jj], pp)
  
}


Test.scalingFactor.track.improvement = FALSE
if(Test.scalingFactor.track.improvement){
  
  library(GenomicAlignments)
  library(rtracklayer)
  
  sfs =  readRDS(file = paste0(RdataDir, '/DESeq2_peaks.based_scalingFactors_forGenomicRanger.rds'))
  
  Normalized.coverage = TRUE
  Logtransform.coverage = FALSE
  
  baseDir = '/Volumes//groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/'
  
  InputDir = paste0(baseDir, '/bam_merged_using')
  OutDir = paste0(baseDir, "bigwigs_reNormed/")
  if(!dir.exists(OutDir)) dir.create(OutDir)
  
  bamlist = list.files(path = InputDir, pattern = "*.bam$", full.names = TRUE)
  bamlist = bamlist[grep('Mature_LA', bamlist)]
  
  Pairend = TRUE
  
  for(n in c(1:length(bamlist)))
  {
    # n = 3
    bam = bamlist[n]
    bw.name = basename(bam)
    bw.name = gsub(".bam", "", bw.name)
    bw.name = gsub("_uniq_rmdup", '', bw.name)
    bw.name = gsub('_mq_30.bw', '', bw.name)
    
    cat("bam file: ", bamlist[n], '-- ', "bw name: ", bw.name, "\n")
    
    if(!file.exists(paste0(OutDir, bw.name))){
      if(Pairend){
        ga = readGAlignmentPairs(bam)
        #ga = readGAlignmentPairs(bam,param = ScanBamParam(flag=scanBamFlag(isDuplicate =FALSE))
      }else{
        ga = readGAlignments(bam)
      }
      
      if(Normalized.coverage){
        scalingfactor = sfs$sf[grep(bw.name, sfs$sample)]
        xx = coverage(granges(ga))/scalingfactor*10^6
        
      }else{
        xx = ga
      }
      
      if(Logtransform.coverage) xx = log2(xx+2^-6)
      
      export.bw(xx, con = paste0(OutDir, bw.name, 'reNorm.bw'))
      
    }
  }
  
  
  
}




