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

# loop over the peaks and correct the bins within peaks
overlaps = findOverlaps(pp, bw)
mappings = match(names(pp), names(xx))

library(tictoc)
tic()
for(n in 1:length(pp))
{
  # n = 7
  kk = overlaps@to[which(overlaps@from == n)]
  bins = bw[kk]
  
  if(length(bins) > 0){
    sf = cpm[which(names(cpm) == names(pp)[n])]/sum(as.numeric(width(bins)) * bins$score) * 50 # 50 is constant to reflect the fragment size
    
    #tic()
    # scale the bigwig scores
    if(!is.na(sf)){
      mcols(newb)$score[kk] = mcols(bw)$score[kk] * sf
      # save the scaling factor
      xx$scalingfactor[mappings[n]] = sf
    }else{
      print(paste0('NA scaling factor found for ',  n, '-- peak name : ', names(pp)[n]))
    }
  }
  
}

toc()


# loop over the gaps between peaks, chr by chr 
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
            
            if(all(!is.na(fractions))){
              mcols(newb)$score[jj] = mcols(bw)$score[jj] * sf_final
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
  
  dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/bigwigs_bc/outs'
  load(paste0(dataDir, '/Mature_LA_74939_bc_peaks.Rdata'))
  
  kk = which(is.na(mcols(newb)$score))
  head(which(mcols(bw)$score[kk]>0))
  
  findOverlaps(bw[kk], pp)
  
  
}



