##########################################################################
##########################################################################
# Project: positional memory 
# Script purpose: process the output after running chromHMM
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Sep 22 10:39:22 2023
##########################################################################
##########################################################################
rm(list = ls())

require(ggplot2)
library(tidyr)
require(dplyr)
require(gridExtra)
require(RColorBrewer)
require(tibble)
library(reshape2)
library(data.table)
library(tidyverse)
require(GenomicRanges)
library('rtracklayer')
library('GenomicFeatures')
library("ggalluvial")
library(tidyr)
library(dplyr)
require(pheatmap)

figureDir = paste0('/Users/jingkui.wang/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/',
                   'Jingkui/Hox Manuscript/DEVELOPMENTAL CELL/REVISION/analysis/plots_jiwang_regeneration') 
tableDir = paste0(figureDir, '/tables')

nb_chromStates = 8
states = paste0('E', c(1:nb_chromStates))

inputDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/',
                  'chromHMM/chromHMM_regeneration_mergedRep/clusters_', nb_chromStates)

resDir = '../results/chromHMM_test'

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##########################################
# functions 
##########################################
Read_chromHMM_segments_into_GenomicRanger = function(x)
{
  x = data.frame(read.table(x, header = FALSE), stringsAsFactors = FALSE)
  #x = data.frame(x, stringsAsFactors = FALSE)
  colnames(x)[4] = 'state'
  x$strand = '*'
  x = makeGRangesFromDataFrame(x, 
                               seqnames.field=c("V1"),
                               start.field="V2", end.field="V3", 
                               strand.field="strand", keep.extra.columns=TRUE)
  return(x)
  
}

intersect_GenomicRanges_keepMetadata = function(gr1, gr2)
{
  o = findOverlaps(gr1, gr2)
  grl1 = split(gr1[queryHits(o)], 1:length(o)) # You can't mendoapply on a GRanges object
  grl2 = split(gr2[subjectHits(o)], 1:length(o))
  foo = function(x, y) {
    rv = x
    start(rv) = max(start(x), start(y))
    end(rv) = min(end(x), end(y))
    return(rv)
  }
  
  return(unlist(mendoapply(foo, grl1, y=grl2)))
  
}

convert_fread_GenomicRange = function(x)
{
  x = data.frame(x, stringsAsFactors = FALSE)
  colnames(x)[4] = 'state'
  x$strand = '*'
  x = makeGRangesFromDataFrame(x, 
                               seqnames.field=c("V1"),
                               start.field="V2", end.field="V3", 
                               strand.field="strand", keep.extra.columns=TRUE)
  return(x)
  
}


########################################################
########################################################
# Section I : here we first show the emission probability of defined chromatin states
# and also the changes of states across segments mUA, mLA and mHand
########################################################
########################################################
emission = read.table(file = paste0(inputDir, '/emissions_', nb_chromStates,'.txt'), sep = '\t', header = TRUE)
emission = data.frame(emission)
colnames(emission)[1] = 'state'
emission1 = reshape2::melt(emission, id.vars = 'state', variable.name = "feature")

ggplot(emission1, aes(x = feature, y = state, fill = value)) +
  geom_tile() +
  labs(title = "emissions",
       x = "",
       y = "ChromHMM states") +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  #coord_fixed() +
  theme_classic() +
  scale_y_continuous(breaks = c(1:nb_chromStates)) +
  theme(axis.text.x = element_text(angle = 90, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
  ) + 
  #guides(fill=guide_legend(title="emission")) +
  ggtitle("")

ggsave(paste0(figureDir, '/Emission_chromStates_nb.clusters.', nb_chromStates, '.pdf'), 
       width=4, height = 5)

rm(emission)
rm(emission1)
#pheatmap(emission[, -1], cluster_rows=TRUE, show_rownames=TRUE, fontsize_row = 12,
# color = colorRampPalette((brewer.pal(n = 4, name ="Reds")))(4), 
# show_colnames = TRUE,
# scale = 'none',
# cluster_cols=FALSE)


##########################################
# transition between chromatin states and clustering
##########################################
library(pheatmap)
transition = read.table(file = paste0(inputDir, '/transitions_', nb_chromStates,'.txt'), sep = '\t', 
                        header = TRUE, row.names = c(1))

colnames(transition) = c(1:nb_chromStates)
rownames(transition) = c(1:nb_chromStates)

pheatmap(transition, cluster_rows=TRUE, show_rownames=TRUE, fontsize_row = 12,
         color = colorRampPalette((brewer.pal(n = 9, name ="OrRd")))(100), 
         show_colnames = TRUE,
         scale = 'none',
         cluster_cols=TRUE,
         #annotation_col=df,
         #annotation_colors = annot_colors,
         width = 5, height = 4, 
         treeheight_row = 20,
         treeheight_col = 20,
         filename = paste0(figureDir, '/chromState_transition_clustering_',
                           nb_chromStates, '.pdf'))

rm(transition)

##########################################
# genome-wide chromatin state coverages and changes in regeneration
# https://stackoverflow.com/questions/68487536/how-to-align-and-label-the-stratum-in-ggalluvial-using-ggrepel-or-otherwise
##########################################
## prepare gene regions to consider
Prepare_geneRegions = FALSE
if(Prepare_geneRegions){
  gtf.all = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/AmexT_v47_Hox.patch.gtf'
  
  amex.all = GenomicFeatures::makeTxDbFromGFF(file = gtf.all)
  # tss = GenomicFeatures::promoters(amex.all, upstream = 5000, downstream = 0, use.names = TRUE)
  gene = GenomicFeatures::genes(amex.all)
  #gene = resize(gene, , fix="end", use.names = TRUE)
  
  gene = as.data.frame(gene)
  gene$start[which(gene$strand == '+')] = gene$start[which(gene$strand == '+')] - 5000
  gene$end[which(gene$strand == '-')] = gene$end[which(gene$strand == '-')] + 5000
  
  #gene$tx_name = sapply(gene$tx_name, function(x){x = unlist(strsplit(as.character(x), '[|]')); x[length(x)]})
  gene$start[which(gene$start<=1)] = 1
  
  gtf.limb =  '../data/AmexT_v47_Hox.patch_limb.fibroblast.expressing.23585.genes.dev.mature.regeneration.gtf'
  ggs = GenomicFeatures::makeTxDbFromGFF(file = gtf.limb)
  ggs = GenomicFeatures::genes(ggs)
  
  mm = match(gene$gene_id, ggs$gene_id)
  gene = gene[which(!is.na(mm)), ]
  
  gene = makeGRangesFromDataFrame(gene, 
                                  seqnames.field=c("seqnames"),
                                  start.field="start", end.field="end", 
                                  strand.field="strand", keep.extra.columns=TRUE)
  
  saveRDS(gene, file = paste0("../data/genes_23585_used_for_histone_chromStates.rds"))
  strand(genes) = '*'
  export(object = genes,  con = paste0(inputDir, "/gene_to_use.bed"), format = 'bed')
  
  genes = data.frame(genes, stringsAsFactors = FALSE)
  write.table(genes[, c(1:3)], col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE,
              file = paste0(inputDir, '/genes_to_used.bed'))
  
  
  # consider only the regions of those significant genes
  limbgenes = readRDS(file = paste0("../data/genes_23585_used_for_histone_chromStates.rds"))
  ggs = readRDS(file = paste0(tableDir, '/significant_gene_list_matureSegment.rds'))
  ggs = sapply(ggs, function(x){test = unlist(strsplit(as.character(x), '_')); test[length(test)]})
  limbgenes = limbgenes[which(!is.na(match(limbgenes$gene_id, ggs)))]
  
  ggs =  data.frame(limbgenes, stringsAsFactors = FALSE)
  ggs$strand = '*'
  write.table(ggs[, c(1:3)], col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE,
              file = paste0(inputDir, '/significant_genes_chromswitch.bed'))
  
}

#overlap_models = bed_intersect(x = mLA[c(1:20), ], list(mUA = mUA[1:20, ], mHand = mHand[1:20, ]))%>%
#  group_by(state.col,state.ddm)%>%summarise(tot=sum(.overlap),.groups='drop_last')%>%
#  mutate(size_col=sum(tot))%>%group_by(state.ddm)%>%mutate(size_ddm=sum(tot))
files_bins = list.files(path = paste0(inputDir, '/binning_out_genes'),
                        pattern = '.bed', recursive = TRUE, full.names = TRUE)
conds = basename(files_bins)
conds = gsub('_8_segments_binned200bp_selectedGenes.bed', '', conds)

kk = c(1, 4, 5, 3, 2)
conds = conds[kk]
files_bins = files_bins[kk]


state_freq = matrix(NA, ncol = length(conds) , nrow = nb_chromStates^length(conds))
index = 1

for(n in 1:nb_chromStates)
{
  for(i in 1:nb_chromStates)
  {
    for(j in 1:nb_chromStates)
    {
      for(m in 1:nb_chromStates)
      {
        for(k in 1:nb_chromStates)
        {
          state_freq[index, 1] = n
          state_freq[index, 2] = i
          state_freq[index, 3] = j
          state_freq[index, 4] = m
          state_freq[index, 5] = k
          index = index + 1 
        }
      }
    }
  }
}

state_freq = data.frame(state_freq, stringsAsFactors = FALSE)
colnames(state_freq) = conds

state_freq$combineStates = apply(state_freq, 1, function(x) {paste0('E', x, collapse = '_')})

state_freq$freqs = 0

segments = vector("list", length = length(conds))
for(n in 1:length(segments))
{
  cat(n, ' -- ', conds[n], '\n')
  segments[[n]] = fread(files_bins[n])
}

library(tictoc)
tic()
transitions = paste0(segments[[1]]$V4, '_', 
                     segments[[2]]$V4, '_',
                     segments[[3]]$V4, '_',
                     segments[[4]]$V4, '_',
                     segments[[5]]$V4)

toc()

counts = table(transitions)
state_freq$freqs = counts[match(state_freq$combineStates, names(counts))]
state_freq = state_freq[which(!is.na(state_freq$freqs)), ]

# for(n in 1:nrow(segments[[1]]))
# #for(n in 1:10^5)
# {
#   # n = 1
#   if(n%%100000 == 0) cat(n, '\n')
#   
#   #xx = findOverlaps(mUA, genes[n])
#   kk = which(state_freq$combineStates == paste0(segments[[1]]$V4[n], "_", 
#                                                 segments, '_', mHand$V4[n]))
#   if(length(kk) == 1) state_freq$freqs[kk] =  state_freq$freqs[kk] + 1
# }
# 
# toc()

saveRDS(state_freq, file = paste0(tableDir, '/state_transition_frequency_nbStates_', nb_chromStates, '.rds'))

##########################################
# make summary of chromHMM and transitions  
##########################################
state_freq = readRDS(paste0(tableDir, '/state_transition_frequency_nbStates_', nb_chromStates, '.rds'))

# because E1_E1_E1 is empty state and not interesting
ss = apply(state_freq[, c(1:5)], 1, function(x){length(unique(x))}) # test if all sampels have the same state
state_freq = state_freq[which(ss>1), ]
#state_freq = state_freq[-which(state_freq$mUA == state_freq$mLA & state_freq$mUA == state_freq$mHand), ]

# count the distribution of each states
conds = colnames(state_freq)[1:5]

data = c()
for(n in 1:length(conds)){
  for(m in 1:nb_chromStates){
    data = rbind(data, c(colnames(state_freq)[n], m, sum(state_freq$freqs[which(state_freq[,n] == m)])))
  }
  #x = table(mUA$state)/nrow(mUA)
  #data = cbind(rep('mUA', nb_chromStates), names(x), x)
}

data = data.frame(data, stringsAsFactors = TRUE)
colnames(data) = c('condition', 'state', 'freq')
data$freq = as.numeric(as.character(data$freq))
data$counts = data$freq
data$freq = data$freq/sum(state_freq$freqs)

# Grouped
ggplot(data, aes(fill=factor(condition, levels = c('dpa0', 'dpa5', 'dpa9', 'dpa13_prox', 'dpa13_dist')), 
                 y=freq, 
                 x= state)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + 
  labs(x = 'ChromHMM states', y = '% of geomic coverage') + 
  guides(fill=guide_legend(title="regeneration")) + 
  theme_classic()  +
  scale_fill_manual(values = cbPalette[1:length(conds)]) + 
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
  ) + 
  #guides(fill=guide_legend(title="emission")) +
  ggtitle("")

ggsave(paste0(figureDir, '/distriubtions_chromStates_acrossSegments_nb.clusters.', 
              nb_chromStates, '.pdf'), width=8, height = 4)


states = data
states$log2fc = 0

for(n in 1:nb_chromStates){
  cat(n, ' --  chrom state \n')
  jj = which(states$state == n)
  
  states$log2fc[jj] = log2(states$freq[jj]/states$freq[which(states$condition == conds[1] & states$state == n)])
  #data = rbind(data, c(colnames(state_freq)[n], m, sum(state_freq$freqs[which(state_freq[,n] == m)])))
  #x = table(mUA$state)/nrow(mUA)
  #data = cbind(rep('mUA', nb_chromStates), names(x), x)
}

states = states[which(states$condition != 'dpa0'), ]

states$state = as.factor(states$state)
#as_tibble(stats) %>%  gather(group, freq,  2:ncol(stats)) %>% 

ggplot(data = states, aes(fill=factor(condition, levels = c('dpa5', 'dpa9', 'dpa13_prox', 'dpa13_dist')), 
                 y=log2fc, 
                 x= state)) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + 
  labs(x = 'ChromHMM states', y = 'log2fc of geomic coverage vs. dpa0') + 
  guides(fill=guide_legend(title="regeneration")) + 
  theme_classic()  +
  scale_fill_manual(values = cbPalette[2:length(conds)]) + 
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
  ) + 
  #guides(fill=guide_legend(title="emission")) +
  ggtitle("")

ggsave(paste0(figureDir, '/log2FC_coverage_chromStates_mHand.vs.mUA_nb.clusters.', 
              nb_chromStates, '.pdf'), width=8, height = 4)


########################################################
########################################################
# Section III: show the gene-chromatin states heatmap of genomic coverage changed between mHand/mUA
# 
########################################################
########################################################

##########################################
# read the replicates of chromHMM output
# not binning bed
##########################################
Import_chromHMM_state_perReplicate = FALSE
if(Import_chromHMM_state_perReplicate){
  nb_chromStates = 8
  bedDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/',
                    'chromHMM/chromHMM_regeneration_rep_makeSegementation/clusters_', nb_chromStates)
  
  bed_files = list.files(path = bedDir, 
                         pattern = '*.bed', full.names = TRUE)
  
  conds = basename(bed_files)
  conds = gsub('_8_segments.bed', '', conds)
  
  gr_list = vector("list", length = length(conds))
  
  for(n in 1:length(conds))
  {
    # n = 2
    cat(n, ' -- ', conds[n], '\n')
    gr_list[[n]] = Read_chromHMM_segments_into_GenomicRanger(bed_files[n])
    
  }
  
  names(gr_list) = conds
  
  gr_list = as(gr_list, "GRangesList")
  
  saveRDS(gr_list, file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/',
                                 'chromHMM/chromswitch',
                                 '/saved_UAregeneration_replicates_chromHMM.stated.rds'))
  
  rm(gr_list)
  
  
}

##########################################
## select genes otherwise it will take too long
## in addition refine the 4843 genes
##########################################
genes = readRDS(file = paste0("../data/genes_23585_used_for_histone_chromStates.rds"))

yy = readRDS(file = paste0('/Users/jingkui.wang/workspace/imp/positional_memory/results/RNAseq_data_used/Rdata/',
                                     'regeneration_geneClusters.rds'))

## manual change the gene annotation MEIS3
rownames(yy)[grep('AMEX60DD024424', rownames(yy))]
rownames(yy)[grep('AMEX60DD024424', rownames(yy))] = 'MEIS3_AMEX60DD024424' 

yy$gene = sapply(rownames(yy), function(x) unlist(strsplit(as.character(x), '_'))[1])
yy$geneID = sapply(rownames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); x[length(x)]})

Make_bedFile_regGenes = FALSE
if(Make_bedFile_regGenes){
  genes = readRDS(file = paste0("../data/genes_23585_used_for_histone_chromStates.rds"))
  
  strand(genes) = '*'
  
  gg.ids = sapply(rownames(yy), function(x){x = unlist(strsplit(as.character(x), '_')); x[length(x)]})
  
  genes = data.frame(genes, stringsAsFactors = FALSE)
  genes = genes[which(!is.na(match(genes$gene_id, gg.ids))), ]
  
  write.table(genes[, c(1:3)], col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE,
              file = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/chromHMM", 
                            "/genes_regeneration.bed"))
  
}

Further_filtering_regGenes = FALSE
if(Further_filtering_regGenes){
  #ggs = sapply(rownames(yy), function(x){unlist(strsplit(as.character(x), '_'))[1]})
  qv.cutoff = 0.01
  logfc.cutoff = 1
  #select = which(yy$fdr.max> -log10(qv.cutoff) & abs(yy$logFC.max)> 0)
  select = which(yy$padj_dpa5.vs.mUA < qv.cutoff & abs(yy$log2FoldChange_dpa5.vs.mUA) > logfc.cutoff|
                   yy$padj_dpa9.vs.mUA < qv.cutoff & abs(yy$log2FoldChange_dpa9.vs.mUA) > logfc.cutoff)
  cat(length(select), ' DE genes selected \n')
  
  ii_keep = grep('SHH|HOXA|HOXD|MEIS|SHOX|FGF', yy$gene)
  yy = yy[unique(c(select, ii_keep)), ]
  
  ggs = yy$gene
  
  ## add more HOXA and HOXD genes here
  sels = unique(c(grep('HOXA|HOXD', genes$gene_id), 
                  which(!is.na(match(genes$gene_id, yy$geneID)))))
  genes = genes[sels]
  
}

annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
            'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
ggs = genes$gene_id
mm = match(ggs, annot$geneID)
kk = which(!is.na(mm))
ggs[kk] = paste0(annot$gene.symbol.toUse[mm[kk]], '_',  annot$geneID[mm[kk]])

genes$gene_id = ggs

saveRDS(genes, file = paste0(tableDir, '/regGenes_chromswitch.rds'))

##########################################
# chromswitch : make a differential test using 
# because it takes long time >2 days for each time point 
# the actual script was running in CBE
##########################################
Test_chromswitch = FALSE
if(Test_chromswitch){
  library(chromswitch)
  #library(rtracklayer)
  #library(BiocParallel)
  #library(tictoc)
  
  states_regeneration = readRDS(file = paste0(tableDir, '/saved_UAregeneration_replicates_chromHMM.stated.rds'))
  genes = readRDS(file = paste0(tableDir, '/regGenes_chromswitch.rds'))
  
  query = genes
  strand(query) = '*'
  
  metadata_all <- data.frame(Sample = names(states_regeneration),
                             Condition = gsub('_rep1|_rep2', '', names(states_regeneration)), 
                             stringsAsFactors = TRUE)
  
  metadata_all
  
  for(cc in c('dpa5', 'dpa9', 'dpa13_prox'))
  {
    # cc = 'dpa5'
    #Paths to the BED files containing peak calls for each sample
    scores = matrix(NA, ncol = nb_chromStates, nrow = length(unique(genes$gene_id)))
    colnames(scores) = paste0('E', c(1:nb_chromStates))
    rownames(scores) = genes$gene_id
    
    jj = which(metadata_all$Condition == cc|metadata_all$Condition == 'dpa0')
    metadata = metadata_all[jj, ]
    
    states_sels = vector("list", length = length(jj))
    for(n in 1:length(jj))
    {
      states_sels[[n]] = states_regeneration[[jj[n]]]
    }
    names(states_sels) = names(states_regeneration)[jj]
    
    for(n in 1:nb_chromStates)
    {
      # n = 1
      peaks = states_sels;
      for(m in 1:length(peaks))
      {
        test =  states_sels[[m]]
        test = test[which(test$state == colnames(scores)[n])]
        peaks[[m]] = test
        rm(test)
      }
      
      #options(mc.cores = 1)
      #tic()
      #param <- SerialParam(bpnworkers = 1)
      #for(k in 1:nrow(scores))
      tic()
      out <- callBinary(query = query,       # Input 1: Query regions
                        metadata = metadata, # Input 2: Metadata dataframe
                        peaks = peaks, # Input 3: Peaks
                        #BPPARAM = multicoreParam,
                        #BPPARAM = SnowParam(workers = 2),
                        filter = FALSE)
      toc()
      
      scores[ ,n] = out$Consensus
      
    }
    
    saveRDS(scores, file = paste0(tableDir, '/savedConsensusScores_mUA_mHand_replicates_chromswitch_', cc,
                                  'vs.mUA.rds'))
    
  }
  
}

##########################################
# plot the chromswitch output run in CBE   
##########################################
file_scores = list.files(path = paste0('/Volumes/groups/tanaka/People/current/jiwang/',
                                       'projects/positional_memory/Data/chromHMM/chromswitch/perTimepoint'),
                         pattern = '*.vs.mUA.rds', full.names = TRUE)

scores_list = vector("list", length = length(file_scores))

for(n in 1:length(file_scores))
{
  # n = 1 
  scores = readRDS(file = file_scores[n])
  scores_list[[n]] = scores
  
  cc = basename(file_scores[n])
  cc = gsub('savedConsensusScores_chromswitch_regenerationReplicates_|.rds', '', cc)
  
  cat(n, ' -- ', cc, '\n')
  
  # thresold used in 
  # https://www.bioconductor.org/packages/release/bioc/vignettes/chromswitch/inst/doc/chromswitch_intro.html
  maxs = apply(scores, 1, max)
  scores = scores[which(maxs > 0.7), ]
  
  
  df <- data.frame(condition = colnames(scores))
  rownames(df) = colnames(scores)
  colnames(df) = 'chromState'
  
  #annotation_colors = ann_colors
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  chrom_cols = cbPalette[1:nrow(df)]
  names(chrom_cols) = rownames(df)
  annot_colors = list(chromSate = chrom_cols)
  
  #df$chromState = cbPalette[1:nrow(df)]
  
  pheatmap(scores, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 4,
           color = colorRampPalette((brewer.pal(n = 4, name ="Reds")))(4), 
           show_colnames = TRUE,
           scale = 'none',
           cluster_cols=FALSE,
           annotation_col=df,
           annotation_colors = annot_colors,
           width = 6, height = 12, 
           filename = paste0(resDir, '/heatmap_chromswitchSocres_', cc, '_nbState.',
                             nb_chromStates, '.pdf'))
  
}

saveRDS(scores_list, file = paste0(resDir, 'chromswitch_scores_significantThreshold.0.7_nbStates.', nb_chromStates, 
        '.rds'))

##########################################
# ## Construct cs*gene matrix for differential genes
##########################################
Import_chromHMM_state_mergedReplicate = FALSE
if(Import_chromHMM_state_mergedReplicate){
  nb_chromStates = 8
  
  bedDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/',
                  'chromHMM/chromHMM_regeneration_mergedRep/clusters_', nb_chromStates)
  
  bed_files = list.files(path = bedDir, 
                         pattern = '*segments.bed', full.names = TRUE)
  
  conds = basename(bed_files)
  conds = gsub('_8_segments.bed', '', conds)
  
  jj = c(1, 4, 5, 3, 2)
  conds = conds[jj]
  bed_files = bed_files[jj]
  
  gr_list = vector("list", length = length(conds))
  
  for(n in 1:length(conds))
  {
    # n = 2
    cat(n, ' -- ', conds[n], '\n')
    gr_list[[n]] = Read_chromHMM_segments_into_GenomicRanger(bed_files[n])
    
  }
  
  names(gr_list) = conds
  
  gr_list = as(gr_list, "GRangesList")
  
  saveRDS(gr_list, file = paste0(resDir, '/saved_UAregeneration_mergedRep_chromHMM.stated.rds'))
  
  rm(gr_list)
  
}

# original code from
# https://bioinformatics.stackexchange.com/questions/874/intersection-of-two-genomic-ranges-to-keep-metadata
Dir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/chromHMM/chromswitch/'
states_regeneration = readRDS(file = paste0(resDir, '/saved_UAregeneration_mergedRep_chromHMM.stated.rds'))
genes = readRDS(file = paste0(Dir, '/regGenes_chromswitch.rds'))

scores_list = readRDS(file = paste0(resDir, 'chromswitch_scores_significantThreshold.0.7_nbStates.',
                                    nb_chromStates, '.rds'))

res = matrix(NA, 
             ncol = nb_chromStates*length(states_regeneration), 
             nrow = length(unique(genes$gene_id)))
rownames(res) = genes$gene_id

colnames(res) = paste0(rep(paste0('E', c(1:nb_chromStates), '_'), each = length(states_regeneration)), 
                       rep(names(states_regeneration), times = nb_chromStates))

#pseudo_cutoff = 1000

library(tictoc)

Test_parallel = FALSE
if(Test_parallel){
  library(parallel)
  library(MASS)
  
  fx <- function(g, states_regeneration, states, res) 
  {
    # g = "HOXA10"
    #kmeans(Boston, 4, nstart=nstart)
    test = res[which(rownames(res) == g), ]
    for(i in 1:length(states_regeneration))
    {
      # i = 1
      overlap = intersect_GenomicRanges_keepMetadata(states_regeneration[[i]], genes[which(genes$gene_id == g)])
      #overlap2 = intersect_GenomicRanges_keepMetadata(mHand, genes[n])
      #xx = findOverlaps(mUA, genes[n])
      #x = GenomicRanges::intersect(mUA, genes[n], ignore.strand=TRUE)
      for(s in states)
      {
        x = sum(width(overlap)[which(overlap$state == s)])
        #cat(x, '\n')
        #x2 = sum(width(overlap2)[which(overlap2$state == colnames(res)[m])])
        test[which(colnames(res) == paste0(s,'_', names(states_regeneration)[i]))] = x
      }
    }
    
    return(test)
    
  }
  
  numCores <- detectCores()
  numCores
  
  ggs = rownames(res)
  
  system.time(
    results <- mclapply(ggs, fx, states_regeneration = states_regeneration, states = states, res = res,
                        mc.cores = numCores)
  )
  
  saveRDS(results, file = paste0(resDir, '/chromState_coverageChanges_mclaapy_regenerations_', 
                                 nb_chromStates, '_v1.rds'))
  
  
  for(n in 1:nrow(res))
  {
    res[n, ] = results[[n]]  
  }
  
  saveRDS(res, file = paste0(resDir, '/chromState_coverageChanges_regenerations_', nb_chromStates, '_v1.rds'))
  
  
}

##########################################
# make the heatmap 
##########################################
scores_list = readRDS(file = paste0(resDir, 'chromswitch_scores_significantThreshold.0.7_nbStates.',
                                    nb_chromStates, '.rds'))

res= readRDS(file = paste0(resDir, '/chromState_coverageChanges_regenerations_', nb_chromStates, '_v1.rds'))

gg_sels = c()
for(n in 1:length(scores_list))
{
  # n = 1
  scores = scores_list[[n]]
  maxs = apply(scores, 1, max)
  # thresold used in 
  # https://www.bioconductor.org/packages/release/bioc/vignettes/chromswitch/inst/doc/chromswitch_intro.html
  scores = scores[which(maxs > 0.7), ] 
  gg_sels = unique(c(gg_sels, rownames(scores)))
  
}

pseudo_cutoff = 1000
res = log2(res + pseudo_cutoff)
#res[which(abs(res) == Inf)] = 0
res = res[which(!is.na(match(rownames(res), gg_sels))), ]

test = res
for(s in states)
{
  # s = states[1]
  jj = grep(s, colnames(res))
  xx = res[, jj]
  range <- 4.0
  xx = t(apply(as.matrix(xx), 1, 
               function(x) {x = (x - x[1]); x[which(x >= range)] = range; 
                            x[which(x<= (-range))] = -range; x}
               )
         )
  
  test[, jj] = xx
  
}

maxs_fc = apply(abs(test), 1, max)

#res = res[order(-maxs),]
sels = which(maxs_fc > 1)
test = test[sels, ]


saveRDS(test, file = paste0(resDir, '/chromState_coverageChanges_regenerations_signifGenes', 
                            nb_chromStates, '_v1.rds'))
saveRDS(rownames(test), file = paste0(resDir, '/chromState_regeneration_signifiGene_list.rds'))


Add_smartseq2_regeneration = FALSE
if(Add_smartseq2_regeneration){
  res = readRDS(file = paste0(resDir, '/chromState_coverageChanges_regenerations_signifGenes', 
                              nb_chromStates, '_v1.rds'))
  
  kk = grep('_dpa0', colnames(res))
  res = res[ ,-kk]
  
  gg1 = sapply(rownames(res), function(x) {x = unlist(strsplit(as.character(x), '_')); x[length(x)]})
               
  ## import smartseq2 data
  yy = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/',
                             'smartseq2_regeneratoin_sampleMean.rds'))
  
  ## manual change the gene annotation MEIS3
  rownames(yy)[grep('AMEX60DD024424', rownames(yy))]
  rownames(yy)[grep('AMEX60DD024424', rownames(yy))] = 'MEIS3_AMEX60DD024424' 
  
  gg2 =  sapply(rownames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); x[length(x)]})
  
  mm = match(gg1, gg2)
  yy = yy[mm, ]
  
  yy = yy - yy[, 1]
  yy = yy[, -c(1)]
  
  range <- 4.0
  
  yy = t(apply(as.matrix(yy), 1, 
               function(x) {x[which(x >= range)] = range; 
               x[which(x<= (-range))] = -range; x})
         )
  
  colnames(yy) = paste0('smartseq2_', c('dpa5', 'dpa9', 'dpa13_prox', 'dpa13_dist'))
  
  res = cbind(res, yy)
  
  saveRDS(res, file = paste0(resDir, '/chromState_coverageChanges_smartseq2_regenerations_signifGenes', 
                              nb_chromStates, '_v1.rds'))
  
}


##########################################
# make heatmap for genes that changes chromatin states
# add also the smart-seq2 data
##########################################
res = readRDS(file = paste0(resDir, '/chromState_coverageChanges_smartseq2_regenerations_signifGenes', 
                            nb_chromStates, '_v1.rds'))

df <- data.frame(chromState = colnames(res),
                 time = colnames(res))
rownames(df) = colnames(res)

df$chromState = sapply(df$chromState, function(x) unlist(strsplit(as.character(x), '_'))[1])
df$time = sapply(df$time, function(x) {x = unlist(strsplit(as.character(x), '_'));
paste0(x[2:length(x)], collapse = '.')})

#annotation_colors = ann_colors
ss = unique(df$chromState)
cbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
chrom_cols = cbPalette[1:length(ss)]

names(chrom_cols) = ss
time_colors = c('springgreen', 'springgreen2', 'springgreen3', 'gold2')
names(time_colors) = unique(df$time)

library(dendextend)
nb_clusters = 12
my_hclust_gene <- hclust(dist(res), method = "complete")
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
#res$clusters = NA; 
#res$clusters[match(names(my_gene_col), rownames(res))] = my_gene_col # save the cluster index in ther res

my_gene_col <- data.frame(cluster =  paste0('cluster_', my_gene_col))
rownames(my_gene_col) = rownames(res)

col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
          "#33a02c", "#fb9a99", "#e31a1c",
          "#fdbf6f", "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99", "#b15928")

cluster_col = col3[1:nb_clusters]
names(cluster_col) = paste0('cluster_', c(1:nb_clusters))

annot_colors = list(time = time_colors,
  chromState = chrom_cols, 
                    cluster = cluster_col)

plt = pheatmap(res, cluster_rows=TRUE, 
         annotation_row = my_gene_col, 
         show_rownames=FALSE, fontsize_row = 4,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
         color =  colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
         show_colnames = TRUE,
         clustering_method = 'complete', cutree_rows = nb_clusters, 
         scale = 'none',
         cluster_cols=FALSE,
         annotation_col=df,
         gaps_col = c(1:8)*4, 
         #treeheight_row = 30,
         #annotation_colors = annot_colors,
         width = 18, height = 12, 
         filename = paste0(resDir, '/heatmap_chromState_coverageChange_regeneration.vs.mUA_chromswitchSignificant_',
                           'nbState.',  nb_chromStates,'.pdf'))

source('Functions_histM.R')
clusters = cutree(plt$tree_row, k = nb_clusters)
clusters = clusters[plt$tree_row$order]
cluster_order = unique(clusters)
clusters[which(clusters == 12| clusters == 9)] = 2
clusters[which(clusters == 10)] = 7
clusters[which(clusters == 11)] = 9

gaps.row = cal_clusterGaps(plt, nb_clusters = nb_clusters)
gaps.row = gaps.row[-1]
gaps.row = gaps.row[1:8]

clusters_old = clusters
cu = unique(clusters)

for(n in 1:length(cu))
{
   clusters[which(clusters_old == cu[n])] = n
}

my_gene_col <- data.frame(cluster =  paste0('cluster_', clusters))
rownames(my_gene_col) = names(clusters)

annot_colors = list(chromState = chrom_cols, time = time_colors,
                    
                    cluster = cluster_col[!is.na(match(names(cluster_col), 
                                                       paste0('cluster_', unique(clusters))))])

df = df[, c(2, 1)]
pheatmap(res[plt$tree_row$order, ], cluster_rows=FALSE, 
         annotation_row = my_gene_col, 
         show_rownames=FALSE, fontsize_row = 4,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
         color =  colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
         show_colnames = TRUE,
         #clustering_method = 'complete', cutree_rows = nb_clusters, 
         scale = 'none',
         cluster_cols=FALSE,
         annotation_col=df,
         gaps_col = c(1:8)*4, 
         gaps_row =  gaps.row, 
         #treeheight_row = 30,
         annotation_colors = annot_colors,
         width = 18, height = 12, 
         filename = paste0(resDir, '/heatmap_chromState_coverageChange_regeneration.vs.mUA_chromswitchSignificant_',
                           'nbState.',  nb_chromStates,'.pdf'))


pheatmap(res[plt$tree_row$order, ], cluster_rows=FALSE, 
         annotation_row = my_gene_col, 
         show_rownames=TRUE, fontsize_row = 0.5,
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
         color =  colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10),
         show_colnames = TRUE,
         #clustering_method = 'complete', cutree_rows = nb_clusters, 
         scale = 'none',
         cluster_cols=FALSE,
         annotation_col=df,
         gaps_col = c(1:8)*4, 
         gaps_row =  gaps.row, 
         #treeheight_row = 30,
         annotation_colors = annot_colors,
         width = 18, height = 80, 
         filename = paste0(resDir, '/heatmap_chromState_coverageChange_regeneration.vs.mUA_chromswitchSignificant_',
                           'nbState.',  nb_chromStates,'_withGeneSymbols.pdf'))


##########################################
# GO term enrichment for each cluster
##########################################
library(enrichplot)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(stringr)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

clean_geneNames = function(gg.expressed)
{
  gg.expressed = unique(unlist(lapply(gg.expressed, 
                                      function(x) { x = unlist(strsplit(as.character(x), '_'));  
                                      return(x[length(x)])})))
  gg.expressed = unique(annot$gene.symbol.toUse[match(gg.expressed, annot$geneID)])
  gg.expressed = gg.expressed[which(gg.expressed != '' & gg.expressed != 'N/A' & !is.na(gg.expressed))]
  
  return(gg.expressed)
}

annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
xx  = readRDS(file = paste0('../results/RNAseq_data_used/Rdata/',
                            'smartseq2_regeneratoin_sampleMean.rds'))
bgs0 = rownames(xx) 
bgs0 = sapply(bgs0, function(x) {x = unlist(strsplit(as.character(x), '_'))[1]})
bgs0 = unique(bgs0[grep('AMEX6', bgs0, invert = TRUE)])

for(c in unique(clusters))
{
  # c = 9
  jj = which(clusters == c)
 
  ggs  = names(clusters)[jj]
  ggs = sapply(ggs, function(x) {x = unlist(strsplit(as.character(x), '_'))[1]})
  ggs = unique(ggs[grep('AMEX6', ggs, invert = TRUE)])
  cat('cluster -- ', c, ': ', length(ggs), ' genes with symbols \n')
  
  # background
  #bgs0 = unique(clean_geneNames(rownames(xx)))
  
  gg.expressed = ggs
  gg.expressed = firstup(tolower(gg.expressed))
  bgs0 = firstup(tolower(bgs0))
  
  gene.df <- bitr(gg.expressed, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  head(gene.df)
  
  bgs0.df <- bitr(bgs0, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  
  #pval.cutoff = 0.05
  ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                   universe     = bgs0.df$ENSEMBL,
                   #universe     = bgs.df$ENSEMBL,
                   #OrgDb         = org.Hs.eg.db,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2, 
                   minGSSize = 10
                   )
  
  #barplot(ego) + ggtitle("Go term enrichment with BP")
  
  #edox <- setReadable(ego, 'org.Mm.eg.db', 'ENSEMBL')
  pdfname = paste0(resDir, '/GOterm_chromStates_cluster_', c, '.pdf')
  pdf(pdfname, width = 8, height = 6)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  dotplot(ego) + ggtitle(paste0('cluster_', c))
  
  dev.off()
  
  write.csv(ego, file = paste0(resDir, "/GOterm_chromStates_cluster_", c, ".csv"), 
            row.names = TRUE)
  
}

########################################################
########################################################
# Section : alluvial plots for state transcription across time points 
# 
########################################################
########################################################

##########################################
# test alluvial plots for significant genes from chromswitch 
##########################################
Import_chromHMM_state_mergedReplicate_binning.for.signifGenes = FALSE
if(Import_chromHMM_state_mergedReplicate_binning.for.signifGenes){
  genes = readRDS(file = paste0(Dir, '/regGenes_chromswitch.rds'))
  
  ggs = readRDS(file = paste0(resDir, '/chromState_regeneration_signifiGene_list.rds'))
  ggs =  data.frame(genes[which(!is.na(match(ggs, genes$gene_id)))], stringsAsFactors = FALSE)
  ggs$strand = '*'
  write.table(ggs[, c(1:3)], col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE,
              file = paste0(Dir, '/signifGenes_chromswitch_regeneration.bed'))
  
  ##########################################
  # import the binned segementation by significant genes
  ##########################################
  bedDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/',
                  'chromHMM/chromHMM_regeneration_mergedRep/clusters_', nb_chromStates, 
                  '/binning_out_chromswitchSigifGenes')
  
  bed_files = list.files(path = bedDir, 
                         pattern = '*.bed', full.names = TRUE)
  
  conds = basename(bed_files)
  conds = gsub('_8_segments_binned200bp_selectedGenes.bed', '', conds)
  
  jj = c(1, 4, 5, 3)
  conds = conds[jj]
  bed_files = bed_files[jj]
  
  gr_list = vector("list", length = length(conds))
  
  for(n in 1:length(conds))
  {
    # n = 2
    cat(n, ' -- ', conds[n], '\n')
    gr_list[[n]] = Read_chromHMM_segments_into_GenomicRanger(bed_files[n])
    
  }
  
  names(gr_list) = conds
  gr_list = as(gr_list, "GRangesList")
  
  saveRDS(gr_list, file = paste0(resDir, '/saved_UAregeneration_mergedRep_chromHMM.states_rebinning200bp_',
                                 'chromswitch.signif.genes.rds'))
  
  rm(gr_list)
  
}

states_regeneration = readRDS(file = paste0(resDir, 
                                            '/saved_UAregeneration_mergedRep_chromHMM.states_rebinning200bp_',
                                            'chromswitch.signif.genes.rds'))

state_freq = matrix(NA, ncol = 4, nrow = nb_chromStates^4)

index = 1
for(n in 1:nb_chromStates)
{
  for(i in 1:nb_chromStates)
  {
    for(j in 1:nb_chromStates)
    {
      for(m in 1:nb_chromStates)
      {
        state_freq[index, 1] = n
        state_freq[index, 2] = i
        state_freq[index, 3] = j
        state_freq[index, 4] = m
        index = index + 1
      }
    }
  }
}

state_freq = data.frame(state_freq, stringsAsFactors = FALSE)
colnames(state_freq) = conds

state_freq$combineStates = paste0('E', state_freq[,1], '_E', state_freq[, 2],
                                  '_E', state_freq[,3], '_E', state_freq[, 4])
state_freq$freqs = 0

library(tictoc)
tic()

test = paste0(states_regeneration[[1]]$state, "_", 
              states_regeneration[[2]]$state, '_', 
              states_regeneration[[3]]$state, '_',
              states_regeneration[[4]]$state)
counts = table(test)

mm = match(state_freq$combineStates, names(counts)) 
state_freq$freqs = counts[mm]
state_freq = state_freq[which(!is.na(state_freq$freqs)), ]

rm(list= c('mm', 'counts', 'test'))

Counting_by_loop = FALSE # too slow
if(Counting_by_loop){
  #for(n in 1:length(states_regeneration[[1]]))
  tic()
  for(n in 1:10^5)
  {
    # n = 1
    if(n%%100 == 0) cat(n, '\n')
    kk = which(state_freq$combineStates == paste0(states_regeneration[[1]]$state[n], "_", 
                                                  states_regeneration[[2]]$state[n], '_', 
                                                  states_regeneration[[3]]$state[n], '_',
                                                  states_regeneration[[4]]$state[n]))
    if(length(kk) == 1) state_freq$freqs[kk] =  state_freq$freqs[kk] + 1
  }
  
  toc()
  
}

saveRDS(state_freq, file = paste0(resDir, '/regeneration_state_transition_frequency_nbStates_', nb_chromStates, '
                                    chromswitch_signifGenes.rds'))


## carefully filtering for the alluvial plot 
state_freq = readRDS(paste0(resDir, '/regeneration_state_transition_frequency_nbStates_', nb_chromStates, '
                                    chromswitch_signifGenes.rds'))
dim(state_freq)

xx = as.matrix(state_freq[, c(1:4)])
xx[which(xx==4)] = 3
xx[which(xx==8)] = 7
xx[which(xx==5)] = 1
xx[which(xx==2)] = 1

xx = data.frame(xx)
kk = which(xx$dpa5 == xx$dpa0 & xx$dpa9 == xx$dpa0 & xx$dpa13_prox == xx$dpa0)

state_freq = state_freq[-kk, ]
dim(state_freq)

state_freq = state_freq[which(state_freq$freqs>0), ]
dim(state_freq)

state_freq = state_freq[order(-state_freq$freqs), ]

# sels = which((state_freq$mUA == 3 & state_freq$mHand == 4)|
#                (state_freq$mUA == 4 & state_freq$mHand == 3)|
#                (state_freq$mUA == 5 & state_freq$mHand == 6)|
#                (state_freq$mUA == 6 & state_freq$mHand == 5))
# state_freq = state_freq[-sels, ]

fractions = state_freq$freqs/sum(state_freq$freqs)
hist(log10(fractions), breaks = 100)

sels0 = which(fractions > 10^-6)
length(sels0)

sels0 = which(fractions > 10^-5)
length(sels0)

sels0 = which(fractions > 10^-4.5)
length(sels0)

sels0 = which(fractions > 10^-4)
length(sels0)

state_freq = state_freq[sels0, ]
dim(state_freq)

state_freq = state_freq[order(-state_freq$freqs), ]


par(mar = c(1,1,1,1)*12, cex = 0.6, xpd=NA)
df = as.data.frame(state_freq)
ll <- unique(c(as.character(df$dpa0)))
grid.col <- rainbow(length(ll))
grid.col <- setNames(grid.col, ll)
# set colours for alluvial plot (this is a little tricky as strata 
# colour ordering is required for ggalluvial strata)
#names(df) <- c("Condition1", "Condition2", "Condition2", "value")
levs1 <- levels(df$dpa5) 
levs2 <- levs1
res1 <- unique(df$dpa5)
res2 <- res1
cond1_cols <- grid.col[levs1[levs1 %in% res1]]
cond2_cols <- grid.col[levs2[levs2 %in% res2]]
columnCols <- c(cond1_cols, cond2_cols)
stratCols <- c(rev(cond1_cols), rev(cond2_cols))

df$dpa0 = as.factor(df$dpa0)
df$dpa5 = as.factor((df$dpa5))
df$dpa9 = as.factor(df$dpa9)
df$dpa13_prox = as.factor(df$dpa13_prox)

ggplot(as.data.frame(df),
       aes(y = freqs,
           axis1 = dpa0, axis2 = dpa5, axis3 = dpa9, axis4 = dpa13_prox)) +
  geom_alluvium(aes(fill = dpa0),
                width = 1/8, reverse = FALSE 
                #knot.pos = 0, 
  ) +
  scale_fill_manual(values = grid.col) +
  #scale_fill_manual(values = cbPalette[1:nb_chromStates]) + 
  #scale_fill_gradient(low = "white", high = "darkgreen") +
  #scale_fill_brewer(type = "qual", palette = "Set1")+
  guides(fill = "none") +
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE, size = 5) +
  scale_x_continuous(breaks = 1:4, labels = c("dpa0", "dpa5", "dpa9", "dpa13_prox")) +
  #coord_flip() +
  #theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, size = 16), 
        #axis.text.y = element_text(angle = 0, size = 16), 
        #axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14),
        axis.ticks.y = element_blank(), # remove tick marks on the vertical
        axis.text.y = element_blank(), # remove numbers from the vertical
        panel.background = element_blank(), # remove the gray background
        panel.grid.major = element_blank(), # remove the major grid lines
        panel.grid.minor = element_blank() # remove the minor grid lines
  ) + 
  ggtitle("ChromHMM state transitions for significant genes")

ggsave(paste0(figureDir, '/chromStates_transition_acrossSegments_nb.clusters.', 
              nb_chromStates, '_signifGenes.pdf'), 
       width=16, height = 10)



