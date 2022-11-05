##########################################################################
##########################################################################
# Project: positional memeory
# Script purpose: prepare files for GEO submission
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jun 25 12:56:35 2020
##########################################################################
##########################################################################
library(openxlsx)
setwd('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/GEO_submission')

tableDir = paste0('~/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/',
                  'SupTables/metadata')

##########################################
# manually collect sample IDs and other information 
##########################################
file.list = list.files(path = tableDir, pattern = '*.csv', full.names = TRUE)
metadata = read.csv2(file = file.list[1])
metadata = metadata[, c(1, 3, 2, 6)]
metadata = rbind(metadata, metadata[c(13:14),])
metadata = metadata[c(1:13, 15, 14, 16),]
metadata$SampleID[c(13:16)] = c('177595', '177597', '177596', '177598')
metadata$samples[c(13:16)] = c('HEAD_177595', 'HEAD_177597', 'Mature_LA_177596', 'Mature_LA_177598')
metadata$protocol = 'atac'

xx = read.csv2(file = file.list[2])
xx = xx[, c(1, 3, 2, 6)]
xx$protocol = 'atac'
metadata = rbind(metadata, xx)
metadata = metadata[match(unique(metadata$SampleID), metadata$SampleID),]

xx = read.csv2(file = file.list[3])
xx = xx[, c(1, 7, 4, 5)]
xx$protocol = 'CT'
colnames(xx) = colnames(metadata)
metadata = rbind(metadata, xx)
#metadata = metadata[match(unique(metadata$SampleID), metadata$SampleID),]

xx = read.csv2(file = file.list[3])
xx = xx[, c(1, 7, 4, 5)]
xx$protocol = 'CT'
colnames(xx) = colnames(metadata)
metadata = rbind(metadata, xx)

xx = read.csv2(file = file.list[4])
xx = xx[, c(1, 2, 3,5)]
xx$protocol = 'smartseq2'
colnames(xx) = colnames(metadata)
metadata = rbind(metadata, xx)

xx = read.csv2(file = file.list[5])
xx = xx[, c(1, 2, 7,4)]
xx$protocol = 'smartseq2'
colnames(xx) = colnames(metadata)
metadata = rbind(metadata, xx)

metadata = metadata[match(unique(metadata$SampleID), metadata$SampleID),]

metadata$samples = paste0(metadata$protocol, '_', metadata$samples)

write.csv2(metadata, file = 'allSamples_merged.csv', row.names = FALSE)

##########################################
# from sample ID to find associated raw data and processed data 
##########################################
raw_list = list.files(path = 'geo_submission_positionMemory', pattern = '*.bam$', 
                      full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
raw_list = c(raw_list, 
            list.files(path = 'geo_submission_positionMemory', pattern = '*.fastq.gz$', 
                                 full.names = TRUE, recursive = TRUE, include.dirs = TRUE))

raw_list = c(raw_list, 
             list.files(path = 'geo_submission_positionMemory', pattern = '*.fq.gz$', 
                        full.names = TRUE, recursive = TRUE, include.dirs = TRUE))


processed_list = list.files(path = 'geo_submission_positionMemory', pattern = '*.bw$', 
                        full.names = TRUE, recursive = TRUE, include.dirs = TRUE)

raw_list = basename(raw_list)
processed_list = basename(processed_list)

metadata$processed = NA
metadata$raw = NA
metadata$raw_2 = ''
metadata$raw_3 = ''
metadata$raw_4 = ''
nb_raw = c()

pair_files = c()

for(n in 1:nrow(metadata))
{
  # n = 1
  cat(n, '--', metadata$samples[n], '\n')
  index_processed = grep(metadata$SampleID[n], processed_list)
  if(length(index_processed) != 1) {
    stop()
  }else{
    metadata$processed[n] = processed_list[index_processed]
  }
  
  index_raw = grep(metadata$SampleID[n], raw_list)
  if(length(index_raw)==0){
    stop()
  }else{
    nb_raw = c(nb_raw, length(index_raw))
    test = raw_list[index_raw]
    test = test[order(test)]
    
    if(length(test) == 1) {metadata$raw[n] = test}
    if(length(test) == 2) {
      metadata$raw[n] = test[1]
      metadata$raw_2[n] = test[2]
      
      pair_files = rbind(pair_files, c(test))
      
    }
    if(length(test) == 4) {
      metadata$raw[n] = test[1]
      metadata$raw_2[n] = test[2]
      metadata$raw_3[n] = test[3]
      metadata$raw_4[n] = test[4]
      
      pair_files = rbind(pair_files, test[c(1:2)])
      pair_files = rbind(pair_files, test[c(3:4)])
    }
    
    #metadata$raw[n] = paste0(raw_list[index_raw], collapse = ';')
    
  }
}
metadata$nb_raw = nb_raw



metadata$pairend = 'paired-end'
metadata$pairend[which(metadata$protocol == 'smartseq2')] = 'single-end'
metadata = metadata[, c(1:5, 8, 6,7,9:11)]


write.csv2(metadata, file = 'allSamples_merged_processed_raw.csv', row.names = FALSE)
saveRDS(metadata, file = 'allSamples_merged_processed.raw.files.rds')

write.csv2(pair_files,  file = 'paired_end_files.csv', 
           row.names = FALSE, quote = FALSE, col.names = FALSE)




##########################################
# get md5sum of processed files and raw data 
##########################################
raw.md5 = read.delim("md5sum_processed.file", header = FALSE, sep = ' ')

# here pair-end raw data and 2 fastq for one samples
samples = matrix(NA, ncol = 2, nrow = nrow(raw.md5))
samples = data.frame(samples)
colnames(samples) = c('file', 'md5sum')


for(n in 1:nrow(raw.md5)){
  cat(n, '--')
  cat(raw.md5$V3[n], '\n')
  samples$file[n] = basename(as.character(raw.md5$V3[n]))
  samples$md5sum[n] = as.character(raw.md5$V1[n])
  
}

write.csv2(samples,  file = 'md5sum_processed_files.csv', 
           row.names = FALSE, quote = FALSE)


raw.md5 = read.delim("md5sum_raw.file", header = FALSE, sep = ' ')

# here pair-end raw data and 2 fastq for one samples
samples = matrix(NA, ncol = 2, nrow = nrow(raw.md5))
samples = data.frame(samples)
colnames(samples) = c('file', 'md5sum')


for(n in 1:nrow(raw.md5)){
  cat(n, '--')
  cat(raw.md5$V3[n], '\n')
  samples$file[n] = basename(as.character(raw.md5$V3[n]))
  samples$md5sum[n] = as.character(raw.md5$V1[n])
  
}

write.csv2(samples,  file = 'md5sum_raw_files.csv', 
           row.names = FALSE, quote = FALSE)


##########################################
# collect raw data  
##########################################
Data.Collection = FALSE
if(Data.Collection){
  copyDir = 'geo_submission_25june'
  fromDir = '../R8024_R7846_R7708/ngs_raw'
  
  system(paste0("mkdir -p ", copyDir))
  file.format = '*.bam'
  
  files.all = list.files(path = fromDir, pattern = file.format, full.names = TRUE, recursive = TRUE)
  files.used = list.files(path = copyDir, pattern = file.format, full.names = TRUE)
  sId = samples$sampleID
  
  for(n in 1:length(sId))
  {
    kk = grep(sId[n], files.used)
    jj = grep(sId[n], files.all)
    
    cat(sId[n], " -- ")
    if(length(kk) != 1){
      
      if(length(kk) == 0){
        cat('Not Found -- ')
        if(length(jj)==1) {
          cat('one file to copy')
          system(paste0('cp ', files.all[jj], ' ', copyDir))
        }
        if(length(jj)==0) cat(" >> no file to move")
        if(length(jj) > 1) cat(" >>too many files to copy")
      }
      if(length(kk) > 1){
        cat("Too Many !!", files.used[kk])
      }
    }else{
      cat('copied already')
    }
    
    cat("\n")
  }
  
  system(paste0('cp ', files.all[jj], ' ', copyDir))
}


########################################################
########################################################
# Section : save count table for processed data
# 
########################################################
########################################################
load(file='/Volumes/groups/cochella/jiwang/Projects/Philipp/smallRNA_analysis_philipp/results/R8024_R7846_R7708_all/Design_Raw_readCounts_UMIfr_miRNAs_R8024_R7846_R7708_20190723.Rdata')

read.count = all[, -1];
rownames(read.count) = all$gene
sel.samples.with.spikeIns = match(samples$sampleID, design.matrix$SampleID)

read.count = read.count[, sel.samples.with.spikeIns]
read.count = read.count[-c(9:11), ]

colnames(read.count) = paste0(samples$sampleID, '.umifr.count')
write.table(read.count, file = 'processed/all_samples_umifr_count.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)

