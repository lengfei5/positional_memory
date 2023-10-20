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

figureDir = paste0('/Users/jingkui.wang/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/',
                   'Jingkui/Hox Manuscript/DEVELOPMENTAL CELL/REVISION/analysis/plots_jiwang') 
tableDir = paste0('/Users/jingkui.wang/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/',
                  'Jingkui/Hox Manuscript/DEVELOPMENTAL CELL/REVISION/analysis/tables_jiwang') 

nb_chromStates = 6
inputDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/',
                  'chromHMM/chromHMM_out_mergedRep_v2/clusters_', nb_chromStates)


Read_chromHMM_segments_into_GenomicRanger = function(segment = 'mUA', nb_chromStates,
                                                     inputDir)
{
  x = read.table(paste0(inputDir, '/', segment, '_', nb_chromStates, '_segments.bed'), header = FALSE)
  x = data.frame(x, stringsAsFactors = FALSE)
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

##########################################
# test alluvial plots  
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

#genes = readRDS(file = paste0("../data/genes_23585_used_for_histone_chromStates.rds"))
mUA = fread(paste0(inputDir, '/mUA_', nb_chromStates, '_segments_binned.200bp_genes.bed'))
mLA = fread(paste0(inputDir, '/mLA_', nb_chromStates, '_segments_binned.200bp_genes.bed'))
mHand = fread(paste0(inputDir, '/mHand_', nb_chromStates, '_segments_binned.200bp_genes.bed'))

state_freq = matrix(NA, ncol = 3, nrow = nb_chromStates^3)
index = 1
for(n in 1:nb_chromStates)
{
  for(i in 1:nb_chromStates)
  {
    for(j in 1:nb_chromStates)
    {
      state_freq[index, 1] = n
      state_freq[index, 2] = i
      state_freq[index, 3] = j
      index = index + 1
    }
  }
}

state_freq = data.frame(state_freq, stringsAsFactors = FALSE)
colnames(state_freq) = c('mUA', 'mLA', 'mHand')
state_freq$combineStates = paste0('E', state_freq$mUA, '_E', state_freq$mLA, '_E', state_freq$mHand)
state_freq$freqs = 0

library(tictoc)
tic()
for(n in 1:nrow(mUA))
#for(n in 1:10^5)
{
  # n = 1
  if(n%%100000 == 0) cat(n, '\n')
  #xx = findOverlaps(mUA, genes[n])
  kk = which(state_freq$combineStates == paste0(mUA$V4[n], "_", mLA$V4[n], '_', mHand$V4[n]))
  if(length(kk) == 1) state_freq$freqs[kk] =  state_freq$freqs[kk] + 1
  
}
toc()

saveRDS(state_freq, file = paste0(tableDir, '/state_transition_frequency_nbStates_', nb_chromStates, '.rds'))

##########################################
# make summary of chromHMM and transitions between mUA and mHand  
##########################################
library("ggalluvial")
library(tidyr)
library(dplyr)

state_freq = readRDS(paste0(tableDir, '/state_transition_frequency_nbStates_', nb_chromStates, '.rds'))

# because E1_E1_E1 is empty state and not interesting
state_freq = state_freq[-which(state_freq$mUA == state_freq$mLA & state_freq$mUA == state_freq$mHand), ]

# count the distribution of each states
data = c()
for(n in 1:3){
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
ggplot(data, aes(fill=factor(condition, levels = c('mUA', 'mLA', 'mHand')), y=freq, 
                 x= state)) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(x = 'ChromHMM states', y = '% of geomic coverage') + 
  guides(fill=guide_legend(title="segments")) + 
  theme_classic()  +
  scale_fill_manual(values = c('palegreen3', 'cyan4', 'goldenrod3')) + 
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
  ) + 
  #guides(fill=guide_legend(title="emission")) +
  ggtitle("")

ggsave(paste0(figureDir, '/distriubtions_chromStates_acrossSegments_nb.clusters.', 
              nb_chromStates, '.pdf'), width=6, height = 4)

states = data.frame(state = c(1:nb_chromStates), 
                    overlap = rep(NA, nb_chromStates),
                    log2fc = rep(NA, nb_chromStates))

for(n in 1:nrow(states)){
  jj = which(state_freq$mUA == states$state[n])
  jj1 = which(state_freq$mUA == states$state[n] & state_freq$mHand == states$state[n])
  jj2 = which(state_freq$mHand == states$state[n])
  states$overlap[n] = sum(state_freq$freqs[jj1])/sum(state_freq$freqs[jj])
  states$log2fc[n] = log2(sum(state_freq$freqs[jj2])/sum(state_freq$freqs[jj]))
}

states$state = as.factor(states$state)
#as_tibble(stats) %>%  gather(group, freq,  2:ncol(stats)) %>% 

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data = states, aes(fill=state, y=log2fc, x=state)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_flip() + 
  theme_classic() +
  labs(x = 'States', y = 'log2FC (mHand/mUA)') + 
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
  ) + 
  scale_fill_manual(values = cbPalette[1:nb_chromStates]) + 
  geom_hline(yintercept = c(-1, 1), colour = "gray", size = 1)

ggsave(paste0(figureDir, '/log2FC_coverage_chromStates_mHand.vs.mUA_nb.clusters.', 
              nb_chromStates, '.pdf'), width=5, height = 6)


########################################################
########################################################
# Section III: show the gene-chromatin states heatmap of genomic coverage changed between mHand/mUA
# 
########################################################
########################################################

##########################################
# read the replicates of chromHMM output
##########################################
Import_chromHMM_state_perReplicate = FALSE
if(Import_chromHMM_state_perReplicate){
  nb_chromStates = 6
  inputDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/',
                    'chromHMM/chromHMM_out_rep_v1/clusters_', nb_chromStates)
  
  mUA_rep1 = Read_chromHMM_segments_into_GenomicRanger(segment = 'mUA_rep1', nb_chromStates = nb_chromStates, 
                                                       inputDir = inputDir)
  mUA_rep2 = Read_chromHMM_segments_into_GenomicRanger(segment = 'mUA_rep2', nb_chromStates = nb_chromStates, 
                                                       inputDir = inputDir)
  mHand_rep1 = Read_chromHMM_segments_into_GenomicRanger(segment = 'mHand_rep1', nb_chromStates = nb_chromStates,
                                                         inputDir = inputDir)
  mHand_rep2 = Read_chromHMM_segments_into_GenomicRanger(segment = 'mHand_rep2', nb_chromStates = nb_chromStates,
                                                         inputDir = inputDir)
  states_mUA.mHand = GRangesList("mUA_rep1" = mUA_rep1, 
                                 "mUA_rep2" = mUA_rep2,
                                 "mHand_rep1" = mHand_rep1,
                                 "mHand_rep2" = mHand_rep2)
  
  #rm(c("mUA_rep1", "mUA_rep2")), mHand_rep1, mHand_rep2))
  
  saveRDS(states_mUA.mHand, file = paste0(tableDir, '/saved_mUA_mHand_replicates_chromHMM.stated.rds'))
  
}

##########################################
# # genes to consider
##########################################
states_mUA.mHand = readRDS(file = paste0(tableDir, '/saved_mUA_mHand_replicates_chromHMM.stated.rds'))

genes = readRDS(file = paste0("../data/genes_23585_used_for_histone_chromStates.rds"))
yy = readRDS(file = paste0("../results/microarray/Rdata/", 
                            'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_normalized_geneSummary_limma.DE.stats_RARB.rds'))

## manual change the gene annotation MEIS3
rownames(yy)[grep('AMEX60DD024424', rownames(yy))]
rownames(yy)[grep('AMEX60DD024424', rownames(yy))] = 'MEIS3_AMEX60DD024424' 
#ggs = sapply(rownames(yy), function(x){unlist(strsplit(as.character(x), '_'))[1]})
qv.cutoff = 0.05
logfc.cutoff = 1
#select = which(yy$fdr.max> -log10(qv.cutoff) & abs(yy$logFC.max)> 0)
select = which(yy$adj.P.Val_mHand.vs.mLA < qv.cutoff & abs(yy$logFC_mHand.vs.mLA) > logfc.cutoff|
                 yy$adj.P.Val_mHand.vs.mUA < qv.cutoff & abs(yy$logFC_mHand.vs.mUA) > logfc.cutoff |
                 yy$adj.P.Val_mLA.vs.mUA < qv.cutoff & abs(yy$logFC_mLA.vs.mUA) > logfc.cutoff )
# cat(length(select), ' positional genes found \n')
cat(length(select), ' DE genes selected \n')
#ggs = ggs[select]
yy = yy[select, ]
yy$gene = sapply(rownames(yy), function(x) unlist(strsplit(as.character(x), '_'))[1])
yy$geneID = sapply(rownames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); x[length(x)]})

sels = unique(c(grep('HOXA|HOXD', genes$gene_id), 
                which(!is.na(match(genes$gene_id, yy$geneID)))))
genes = genes[sels]

annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
            'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
ggs = genes$gene_id
mm = match(ggs, annot$geneID)
kk = which(!is.na(mm))
ggs[kk] = paste0(annot$gene.symbol.toUse[mm[kk]], '_',  annot$geneID[mm[kk]])

genes$gene_id = ggs

##########################################
# make a differential test using chromswitch 
##########################################
library(chromswitch)
#library(rtracklayer)
#library(BiocParallel)
#library(tictoc)

query = genes
strand(query) = '*'
metadata <- data.frame(Sample = names(states_mUA.mHand),
                       Condition = c('mUA', 'mUA', 'mHand', 'mHand'), 
                       stringsAsFactors = TRUE)
metadata

#Paths to the BED files containing peak calls for each sample
scores = matrix(NA, ncol = nb_chromStates, nrow = length(unique(genes$gene_id)))
colnames(scores) = paste0('E', c(1:nb_chromStates))
rownames(scores) = genes$gene_id

for(n in 1:nb_chromStates)
{
  # n = 1
  peaks = states_mUA.mHand;
  for(m in 1:length(peaks))
  {
    test =  states_mUA.mHand[[m]]
    test = test[which(test$state == colnames(scores)[n])]
    peaks[[m]] = test
    rm(test)
  }
  
  #options(mc.cores = 1)
  #tic()
  #param <- SerialParam(bpnworkers = 1)
  for(k in 1:nrow(scores)){
    cat(k, '--\n')
    out <- callBinary(query = query[k],       # Input 1: Query regions
                      metadata = metadata, # Input 2: Metadata dataframe
                      peaks = peaks, # Input 3: Peaks
                      #BiocParallel::multicoreWorkers(),
                      filter = FALSE)
    scores[k,n] = out$Consensus
  }
  #toc()
  
}

saveRDS(scores, file = paste0(tableDir, '/savedConsensusScores_mUA_mHand_replicates_chromswitch.rds'))


##########################################
# plot the chromswitch output run in CBE   
##########################################
scores = readRDS(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/',
                               'projects/positional_memory/Data/chromHMM/',
                               'savedConsensusScores_mUA_mHand_replicates_chromswitch_v1.rds'))
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
         filename = paste0(figureDir, '/heatmap_chromswitchSocres_mHand.vs.mUA_nbState.',
                           nb_chromStates, '.pdf')) 

##########################################
# ## Construct cs*gene matrix for differential genes
##########################################
# original code from
# https://bioinformatics.stackexchange.com/questions/874/intersection-of-two-genomic-ranges-to-keep-metadata
mUA = Read_chromHMM_segments_into_GenomicRanger(segment = 'mUA', nb_chromStates = nb_chromStates, 
                                                           inputDir = inputDir)
mHand = Read_chromHMM_segments_into_GenomicRanger(segment = 'mHand', nb_chromStates = nb_chromStates, 
                                                inputDir = inputDir)

res = matrix(NA, ncol = nb_chromStates, nrow = length(unique(genes$gene_id)))
colnames(res) = paste0('E', c(1:nb_chromStates))
rownames(res) = ggs

pseudo_cutoff = 1000

library(tictoc)
tic()

for(n in 1:nrow(res))
#for(n in 1:100)
{
  # n = 289
  #if(n%%100 == 0) 
  cat(n, '\n')
  
  overlap1 = intersect_GenomicRanges_keepMetadata(mUA, genes[n])
  overlap2 = intersect_GenomicRanges_keepMetadata(mHand, genes[n])
  #xx = findOverlaps(mUA, genes[n])
  #x = GenomicRanges::intersect(mUA, genes[n], ignore.strand=TRUE)
  
  for(m in 1:ncol(res))
  {
    x1 = sum(width(overlap1)[which(overlap1$state == colnames(res)[m])])
    x2 = sum(width(overlap2)[which(overlap2$state == colnames(res)[m])])
    
    if(x1 > pseudo_cutoff & x2 > pseudo_cutoff) {
      res[n, m] = x2/x1
    }else{
      if(x1 > pseudo_cutoff & x2 < pseudo_cutoff){
        res[n, m] = (x2 + pseudo_cutoff)/x1
      }
      if(x1 <= pseudo_cutoff & x2 > pseudo_cutoff){
        res[n, m] = (x2)/(x1 + pseudo_cutoff)
      }
    }
  }
}

toc()

saveRDS(res, file = paste0(tableDir, '/chromState_coverageChanges_mHand.vs.mUA_', nb_chromStates, '_v2.rds'))

##########################################
# make the heatmap 
##########################################
scores = readRDS(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/',
                               'projects/positional_memory/Data/chromHMM/',
                               'savedConsensusScores_mUA_mHand_replicates_chromswitch_v1.rds'))
# thresold used in 
# https://www.bioconductor.org/packages/release/bioc/vignettes/chromswitch/inst/doc/chromswitch_intro.html
maxs = apply(scores, 1, max)

scores = scores[which(maxs > 0.7), ] 

require(pheatmap)
res= readRDS(file = paste0(tableDir, '/chromState_coverageChanges_mHand.vs.mUA_', nb_chromStates, '_v2.rds'))

res = log2(res)
#res[which(abs(res) == Inf)] = 0
res[which(is.na(res))] = 0

maxs_fc = apply(abs(res), 1, max)
#res = res[order(-maxs),]

sels = which((maxs > 0.7 & maxs_fc >= 2) | (maxs_fc >= 2 & maxs > 0.2) ,)

saveRDS(names(sels), file = paste0(tableDir, '/significant_gene_list_matureSegment.rds'))
res = res[sels, ]
#res = res[which(maxs>1), ]

range <- 4.0
test = t(apply(res, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))

df <- data.frame(condition = colnames(res))
rownames(df) = colnames(test)
colnames(df) = 'chromState'

#annotation_colors = ann_colors
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
chrom_cols = cbPalette[1:nrow(df)]
names(chrom_cols) = rownames(df)
annot_colors = list(chromState = chrom_cols)

#annot_colors = c(1:ncol(res))
#names(annot_colors) = c('mUA', 'mLA', 'mHand')
#annot_colors = list(segments = annot_colors)
#sample_colors = c('springgreen4', 'steelblue2', 'gold2')
#names(sample_colors) = c('mUA', 'mLA', 'mHand')
#annot_colors = list(segments = sample_colors)

pheatmap(test, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 4,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
         show_colnames = TRUE,
         scale = 'none',
         cluster_cols=FALSE,
         annotation_col=df,
         #annotation_colors = annot_colors,
         width = 6, height = 12, 
         filename = paste0(figureDir, '/heatmap_chromState_coverageChange_mHand.vs.mUA_chromswitchSignificant_',
                           'nbState.',  nb_chromStates,'.pdf'))

pheatmap(test, cluster_rows=TRUE, show_rownames=TRUE, fontsize_row = 4,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
         show_colnames = TRUE,
         scale = 'none',
         cluster_cols=FALSE,
         annotation_col=df,
         annotation_colors = annot_colors,
         width = 6, height = 12, 
         filename = paste0(figureDir, '/heatmap_chromState_coverageChange_mHand.vs.mUA_chromswitchSignificant_',
                            'withGeneNames_nbState.',  nb_chromStates,'.pdf'))


##########################################
# test alluvial plots for significant genes from chromswitch 
##########################################
Make_chromState_transition_AlluvialPlot = FALSE
if(Make_chromState_transition_AlluvialPlot){
  
  #genes = readRDS(file = paste0("../data/genes_23585_used_for_histone_chromStates.rds"))
  mUA = fread(paste0(inputDir, '/mUA_', nb_chromStates, '_segments_binned.200bp_sigGenes.bed'))
  mLA = fread(paste0(inputDir, '/mLA_', nb_chromStates, '_segments_binned.200bp_sigGenes.bed'))
  mHand = fread(paste0(inputDir, '/mHand_', nb_chromStates, '_segments_binned.200bp_sigGenes.bed'))
  
  state_freq = matrix(NA, ncol = 3, nrow = nb_chromStates^3)
  index = 1
  for(n in 1:nb_chromStates)
  {
    for(i in 1:nb_chromStates)
    {
      for(j in 1:nb_chromStates)
      {
        state_freq[index, 1] = n
        state_freq[index, 2] = i
        state_freq[index, 3] = j
        index = index + 1
      }
    }
  }
  
  state_freq = data.frame(state_freq, stringsAsFactors = FALSE)
  colnames(state_freq) = c('mUA', 'mLA', 'mHand')
  state_freq$combineStates = paste0('E', state_freq$mUA, '_E', state_freq$mLA, '_E', state_freq$mHand)
  state_freq$freqs = 0
  
  library(tictoc)
  tic()
  for(n in 1:nrow(mUA))
  #for(n in 1:10^5)
  {
    # n = 1
    if(n%%10000 == 0) cat(n, '\n')
    #xx = findOverlaps(mUA, genes[n])
    kk = which(state_freq$combineStates == paste0(mUA$V4[n], "_", mLA$V4[n], '_', mHand$V4[n]))
    if(length(kk) == 1) state_freq$freqs[kk] =  state_freq$freqs[kk] + 1
    
  }
  toc()
  
  saveRDS(state_freq, file = paste0(tableDir, '/state_transition_frequency_nbStates_', nb_chromStates, '
                                    chromswitch_signifGenes.rds'))
  
  
  ## filtering 
  state_freq = readRDS(file = paste0(tableDir, '/state_transition_frequency_nbStates_', nb_chromStates, '
                                    chromswitch_signifGenes.rds'))
  
  
  state_freq = state_freq[-which(state_freq$mUA == state_freq$mLA & state_freq$mUA == state_freq$mHand), ]
  state_freq = state_freq[which(state_freq$freqs>0), ]
  state_freq = state_freq[order(state_freq$freqs), ]
  fractions = state_freq$freqs/sum(state_freq$freqs)
  
  state_freq = state_freq[which(fractions > 0.001),]
  
  sels = which(abs(state_freq[,3] - state_freq[,1])>1)
  state_freq = state_freq[sels, ]
  state_freq = state_freq[order(-state_freq$freqs), ]
  
  # count the distribution of each states
  data = c()
  for(n in 1:3){
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
  
  #state_freq =state_freq[grep('E3_E4_E6|E3_E6_E4|E4_E3_E6|E4_E6_E3|E6_E3_E4|E6_E4_E3', 
  #state_freq$combineStates), ]
  
  
  state_freq$mUA_nbStates = NA
  state_freq$mLA_nbStates = NA
  
  for(n in 1:nrow(state_freq))
  {
    state_freq$mUA_nbStates[n] = data$counts[which(data$condition == 'mUA' & data$state == state_freq$mUA[n])]
    state_freq$mLA_nbStates[n] = data$counts[which(data$condition == 'mLA' & data$state == state_freq$mLA[n])]
  }
  
  state_freq$freqs_norm = state_freq$freqs/state_freq$mUA_nbStates/state_freq$mLA_nbStates
  state_freq = state_freq[order(-state_freq$freqs_norm), ]
  
  par(mar = c(1,1,1,1)*12, cex = 0.6, xpd=NA)
  # generate some example data
  #somelongnames <- c("homo sapiens", "homo sapiens", letters[18],
  #                   "some other long name", letters[seq(4)])
  #df <- data.frame(x = factor(somelongnames),
  #y = factor(c("this label is long", "Golgi", 
  #              letters[13:18])),
  # count = c(2, 10, 4, 5, 5, 1, 9, 3))
  df = as.data.frame(state_freq)
  ll <- unique(c(as.character(df$mUA)))
  grid.col <- rainbow(length(ll))
  grid.col <- setNames(grid.col, ll)
  
  # set colours for alluvial plot (this is a little tricky as strata 
  # colour ordering is required for ggalluvial strata)
  #names(df) <- c("Condition1", "Condition2", "Condition2", "value")
  levs1 <- levels(df$mLA) 
  levs2 <- levs1
  res1 <- unique(df$mLA)
  res2 <- res1
  cond1_cols <- grid.col[levs1[levs1 %in% res1]]
  cond2_cols <- grid.col[levs2[levs2 %in% res2]]
  columnCols <- c(cond1_cols, cond2_cols)
  stratCols <- c(rev(cond1_cols), rev(cond2_cols))
  
  df$mUA = as.factor(df$mUA)
  df$mLA = as.factor((df$mLA))
  ggplot(as.data.frame(df),
         aes(y = freqs,
             axis1 = mUA, axis2 = mLA, axis3 = mHand)) +
    geom_alluvium(aes(fill = mUA),
                  width = 1/8, knot.pos = 0, reverse = FALSE) +
    scale_fill_manual(values = grid.col) +
    guides(fill = "none") +
    geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              reverse = FALSE) +
    scale_x_continuous(breaks = 1:3, labels = c("mUA", "mLA", "mHand")) +
    #coord_flip() +
    theme_classic() + 
    ggtitle("ChromHMM state transitions for significant genes")
  
  
  ggsave(paste0(figureDir, '/chromStates_transition_acrossSegments_nb.clusters.', 
                nb_chromStates, '_signifGenes.pdf'), 
         width=10, height = 8)
  
  
  ##########################################
  # plot the transition rate between states 
  ##########################################
  transition = c()
  for(n in 1:nb_chromStates)
  {
    for(m in 1:nb_chromStates)
    {
      jj = which(state_freq$mUA == n & state_freq$mHand == m)
      transition = rbind(transition, c(n, m, sum(state_freq$freqs[jj]), state_freq$mUA_nbStates[jj[1]]))
      
    }
  }
  
  transition = data.frame(transition)
  colnames(transition) = c('mUA', 'mHand', 'freqs', 'mUAcounts')
  transition$freqs_norm = (transition$freqs/transition$mUAcounts)
  
  transition$mUA = as.factor(transition$mUA)
  transition$mHand = as.factor(transition$mHand)
  
  ggplot(transition, aes(x = mHand, y = mUA, fill= freqs_norm)) + 
    geom_tile()
  
  
  # #df = as.data.frame(HairEyeColor)
  # ## move the lables
  # df = df[, c(1:3, 5)]
  # df_expanded <- df[rep(row.names(df), df$freqs), ]
  # df_expanded <- df_expanded %>%
  #   mutate(id = row_number()) %>%
  #   pivot_longer(-c(freqs, id), names_to = "Condition", values_to = "state")
  # 
  # # plot alluvial diagram
  # q <- ggplot(df_expanded, aes(x = Condition, stratum = state, alluvium = id, fill = state)) +
  #   geom_flow(width = 0) +
  #   scale_fill_manual(values = columnCols) +
  #   scale_color_manual(values = stratCols) +
  #   geom_stratum(width = 1 / 8, color = "white") +
  #   scale_x_discrete(
  #     expand = c(.25, .25)
  #   ) +
  #   scale_y_continuous(breaks = NULL) +
  #   theme_minimal() +
  #   theme(
  #     axis.ticks.y = element_blank(),
  #     axis.text.y = element_blank(),
  #     panel.grid.major.y = element_blank(),
  #     panel.grid.major.x = element_blank()
  #   ) +
  #   theme(legend.position = "none") +
  #   ylab(NULL)
  # 
  # q +
  #   geom_text(
  #     aes(
  #       label = after_stat(stratum),
  #       hjust = ifelse(Condition == "Condition1", 1, 0),
  #       x = as.numeric(factor(Condition)) + .075 * ifelse(Condition == "Condition1", -1, 1),
  #       color = after_stat(stratum)
  #     ),
  #     stat = "stratum", fontface = "bold", size = 3
  #   )
  # 
  # plot(q)
}


