##########################################################################
##########################################################################
# Project: postional memory project 
# Script purpose: save old version of code or  code not used anymore
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jun  2 14:52:35 2022
##########################################################################
##########################################################################

########################################################
########################################################
# Section : plot postional peaks and chromatin features together 
# 
########################################################
########################################################
##########################################
# combine atac-seq peaks and histone marker to make heatmap 
# here histone peaks were overlapping atac-seq peaks or the atac-seq peak regions
# the initial clusters were done based on the atac-seq peaks 
# subclusters were done with histone makrers
##########################################
Assembly_histMarkers_togetherWith_ATACseq_regeneration = function()
{
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  require(pheatmap)
  require(RColorBrewer)
  
  # import atac-seq peak
  # res = readRDS(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration.rds')) # stat of dynamic peaks
  # z-score of data and heatmap (variable: plt, yy and res)
  load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters.Rdata')) 
  #design_atac = design
  res_atac = res
  yy_atac = yy
  
  Test_cluster_order = FALSE
  if(Test_cluster_order){
    test = yy
    test = test[plt$tree_row$order, ]
    pheatmap(test, 
             nnotation_row = my_gene_col, 
             annotation_col = df, show_rownames = FALSE, scale = 'none', 
             color = col, 
             show_colnames = FALSE,
             cluster_rows = TRUE, cluster_cols = FALSE,  
             #clustering_method = 'complete', cutree_rows = nb_clusters, 
             annotation_colors = annot_colors, 
             #clustering_callback = callback,
             gaps_col = ii.gaps)
    
    ## save group bed of each cluster for deeptools heatmap
    load(file = paste0(RdataDir, '/dynamic_ATACpeaks_regeneration_data.heatmap_DPGPclusters.Rdata')) 
    
    res = res[match(rownames(yy), rownames(res)), ]
    mcs = unique(res$clusters)
    for(n in 1:length(mcs)){
      jj = which(res$clusters == mcs[n])
      pp_atac = data.frame(t(sapply(rownames(yy)[jj], function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
      pp_atac$name = rownames(yy)[jj]
      pp_atac$score = 0
      pp_atac$strand = '*'
      
      peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                            'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_groups')
      write.table(pp_atac, file = paste0(peakGroupDir, '/peak_group_', n, '.bed'), 
                  sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    }
    
    ## save bed files of promoters, enhancers of up- and down-regulated atac-seq peaks
    res = readRDS(file = paste0(RdataDir, '/renegeration_dynamicPeaks_GPDPclustering.merged.extended.rds'))
    res = res[which(!is.na(res$clusters)), ]
    
    pp = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    pp$strand = '*'
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    pp.annots = as.data.frame(pp.annots)
    
    res$atacseq = 'up'
    res$atacseq[which(res$clusters == 'mc1'| res$clusters == 'mc6')] = 'down'
    
    res$annots = NA
    res$annots[grep('Promoter', pp.annots$annotation)] = 'promoter'
    res$annots[grep('Intergenic|Intron', pp.annots$annotation)] = 'enhancer'
    
    pp_atac = data.frame(t(sapply(rownames(res), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
    pp_atac$name = rownames(res)
    pp_atac$score = 0
    pp_atac$strand = '*'
    
    peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                          'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_promoter_enhancer')
    
    jj = which(res$atacseq == 'up' & res$annots == 'enhancer')
    write.table(pp_atac[jj, ], file = paste0(peakGroupDir, '/peak_enhancer_up.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    jj = which(res$atacseq == 'down' & res$annots == 'enhancer')
    write.table(pp_atac[jj, ], file = paste0(peakGroupDir, '/peak_enhancer_down.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    
    jj = which(res$atacseq == 'up' & res$annots == 'promoter')
    xx = pp_atac[jj, ]
    xx0 =  pp.annots[jj, ]
    xx$strand = xx0$geneStrand # keep tss strand
    kk = which(xx$strand == '1')
    xx$X2 = as.numeric(as.character(xx$X2))
    xx$X3 = as.numeric(as.character(xx$X3))
    xx$X2[kk] = xx0$geneStart[kk]; xx$X2[-kk] = xx0$geneStart[-kk] - 1 
    xx$X3[kk] = xx$X2[kk] + 1; xx$X3[-kk] = xx$X2[-kk]
    xx$strand[kk] = '+'
    xx$strand[-kk] = '-'
    
    write.table(xx, file = paste0(peakGroupDir, '/peak_tss_up.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    jj = which(res$atacseq == 'down' & res$annots == 'promoter')
    xx = pp_atac[jj, ]
    xx0 =  pp.annots[jj, ]
    xx$strand = xx0$geneStrand # keep tss strand
    kk = which(xx$strand == '1')
    xx$X2 = as.numeric(as.character(xx$X2))
    xx$X3 = as.numeric(as.character(xx$X3))
    xx$X2[kk] = xx0$geneStart[kk]; xx$X2[-kk] = xx0$geneStart[-kk] - 1 
    xx$X3[kk] = xx$X2[kk] + 1; xx$X3[-kk] = xx$X2[-kk]
    xx$strand[kk] = '+'
    xx$strand[-kk] = '-'
    
    
    write.table(xx, file = paste0(peakGroupDir, '/peak_tss_down.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    ## save promoters to test histone markers
    peakGroupDir = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/',
                          'positional_memory/Data/histMod_CT_using/heatmaps_deeptools/peak_tss')
    
    tss = promoters(amex, upstream=2000, downstream=2000, use.names=TRUE)
    ggs = genes(amex)
    ggs = as.data.frame(ggs)
    
    tss = ggs
    jj = which(tss$strand == '-')
    tss$start[jj] = tss$end[jj]
    tss$end = tss$start + 1
    #tss$strand = '*'
    tss$score = 0
    tss = tss[, c(1, 2, 3, 6, 7, 5)] 
    
    write.table(tss, file = paste0(peakGroupDir, '/tss_expressed.genes.bed'), 
                sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  # import histone marker 
  RdataHistM = '/Users/jiwang/workspace/imp/positional_memory/results/CT_merged_20220328/Rdata'
  load(file = paste0(RdataHistM, '/combined_4histMarkers_overlapped55kATACseq_DE_regeneration.Rdata')) # variables (keep and DE.locus)
  design = readRDS(file = paste0(RdataHistM, '/histM_CT_design_info.rds'))
  
  yy = keep[, grep('^H3K4me3|^H3K27me3|^H3K4me1|^H3K27ac', colnames(keep))]
  
  sampleID = sapply(colnames(yy), function(x) unlist(strsplit(as.character(x), '_'))[3])
  mm = match(sampleID, design$sampleID)
  colnames(yy) = paste0(design$condition[mm], '_', design$Batch[mm], '_', design$sampleID[mm])
  
  source('Functions_histM.R')
  conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
  conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')
  cc = paste(rep(conds_histM, each = length(conds)), conds, sep = "_")
  yy = cal_sample_means(yy, conds = cc)
  
  ## select the regions overlapping with atac-seq
  pp_atac = data.frame(t(sapply(rownames(yy_atac), function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  pp_atac$strand = '*'
  pp_atac = makeGRangesFromDataFrame(pp_atac, seqnames.field=c("X1"),
                                     start.field="X2", end.field="X3", strand.field="strand")
  
  pp_histM = data.frame(t(sapply(rownames(yy), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp_histM$strand = '*'
  pp_histM = makeGRangesFromDataFrame(pp_histM, seqnames.field=c("X1"),
                                      start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(pp_atac, pp_histM)
  
  ## peaks mapped to postional atac-peaks 
  yy1 = yy[mapping@to, ]
  yy_atac = yy_atac[mapping@from,]
  res_atac = res_atac[mapping@from, ]
  
  ## peaks not mapped to positional atac-peaks
  yy0 = yy[-mapping@to, ]
  
  ss = apply(DE.locus[, c(1:3)], 1, sum)
  mm = match(names(ss[which(ss>0)]), rownames(yy0))
  length(which(!is.na(mm)))
  
  yy0 = yy0[mm[!is.na(mm)], ]
  
  #plot.pair.comparison.plot(yy1[, grep('H3K4me1_mUA_', colnames(yy1))], linear.scale = FALSE)
  #plot.pair.comparison.plot(yy1[, grep('H3K4me3_mUA_', colnames(yy1))], linear.scale = FALSE)
  #plot.pair.comparison.plot(yy1[, grep('H3K27me3_mUA_', colnames(yy1))], linear.scale = FALSE)
  saveRDS(yy1, file = paste0(RdataDir, '/peak_signals_atac_4histM_regenerationPeaks.rds'))
  saveRDS(yy0, file = paste0(RdataDir, '/peak_signals_atac_4histM_notOverlapped.regenerationPeaks.rds'))
  
  #yy1 = yy1[, grep('rRep', colnames(yy1))]
  #yy0 = yy0[, grep('rRep', colnames(yy0))]
  yy1 = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_regenerationPeaks.rds'))
  yy0 = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_notOverlapped.regenerationPeaks.rds'))
  
  # not consider H3K27ac in the main figure
  # yy1 = yy1[, grep('H3K27ac', colnames(yy1), invert = TRUE)]
  
  mm = match(rownames(yy1), rownames(DE.locus))
  DE.peaks = DE.locus[mm, ]
  ss = apply(DE.peaks, 1, sum)
  length(which(ss>0))
  
  source('Functions_histM.R')
  conds_histM = c('H3K4me1', 'H3K27me3', 'H3K4me3', 'H3K27ac')
  stable.multiplier = 2
  for(n in 1:length(conds_histM))
  {
    # n = 1
    jj = grep(conds_histM[n], colnames(yy1))
    
    # yy1[,jj] = t(apply(yy1[,jj], 1, cal_z_score))
    # dynamic.list = DE.locus[, which(colnames(DE.locus) == conds_histM[n])]
    # names(dynamic.list) = rownames(DE.locus)
    # dynamic.list = dynamic.list[which(dynamic.list>0)]
    # kk = which(is.na(match(rownames(yy1), names(dynamic.list))))
    # cat(length(kk), '-', conds_histM[n], ' stable \n')
    # yy1[kk, jj] = yy1[kk, jj] /stable.multiplier
    # yy1[,jj] = t(apply(yy1[,jj], 1, cal_centering))
    yy1[ ,jj] = t(apply(yy1[,jj], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
    
    jj0 = grep(conds_histM[n], colnames(yy0))
    # yy0[,jj] = t(apply(yy0[,jj], 1, cal_z_score))
    # kk0 = which(is.na(match(rownames(yy0), names(dynamic.list))))
    # #cat(length(kk), '-', conds_histM[n], ' stable \n')
    # yy0[kk0, jj] = yy0[kk0, jj] /stable.multiplier
    yy0[ ,jj0] = t(apply(yy0[,jj0], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
    
  }
  
  save(yy1, res_atac, yy_atac, file = paste0(RdataDir, '/atac_histM_data_regeneration.Rdata'))
  
  ## peaks overlapped with atac-seq 
  df = as.data.frame(sapply(colnames(yy1), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'time'
  rownames(df) = colnames(yy1)
  
  sample_colors = c('springgreen', 'springgreen2', 'springgreen3', 'gold2', 'red')
  annot_colors = list(time = sample_colors)
  
  gaps_col = c(5, 10, 15)
  pheatmap(yy1, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           gaps_col = gaps_col,
           legend = TRUE,
           #annotation_colors = annot_colors,
           width = 5, height = 10, 
           filename = paste0(resDir, '/heatmap_histoneMarker_16k.regenerationPeaks.pdf'))
  
  ## peaks not overlapped with atac
  df = as.data.frame(sapply(colnames(yy0), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'time'
  rownames(df) = colnames(yy1)
  sample_colors = c('springgreen', 'springgreen2', 'springgreen3', 'gold2', 'red')
  annot_colors = list(time = sample_colors)
  
  pheatmap(yy0, cluster_rows=TRUE, show_rownames=FALSE, fontsize_row = 5,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(12), 
           show_colnames = FALSE,
           scale = 'none',
           cluster_cols=FALSE, annotation_col=df,
           gaps_col = gaps_col,
           legend = TRUE,
           #annotation_colors = annot_colors,
           width = 5, height = 10, 
           filename = paste0(figureDir, '/heatmap_histoneMarker_DE_notoverlapped.with.atac.regenerationPeaks.pdf'))
  
  
  ##########################################
  #  # reorder the histM according to the atac-seq peak clusters
  # cluster order : 6, 1, 5, 3, 4, 2
  # gaps.row = c(32, 32+76, 31 + 32 + 76, 292 + 31 + 32 + 76, 292 + 31 + 32 + 76 + 103)
  ##########################################
  load(file = paste0(RdataDir, '/atac_histM_data_regeneration.Rdata'))
  peaks = res_atac
  table(peaks$clusters)
  cluster_order = paste0('mc', c(1:8))
  
  conds_histM = c('H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
  conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')
  
  ## specify row gaps as atac-seq peaks
  gaps.row = c()
  for(n in 1:(length(cluster_order)-1))
  {
    if(n == 1)  {gaps.row = c(gaps.row, length(which(peaks$clusters == cluster_order[n])))
    }else{
      gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(peaks$clusters == cluster_order[n])))
    }
  }
  
  ## refine histone subclusters for each ATACseq peak cluster
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  
  for(ac in cluster_order)
  {
    # ac = 6
    kk = which(peaks$cluster == ac)
    cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
    
    hm_hclust <- hclust(dist(as.matrix(yy1[kk,])), method = "complete")
    #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
    peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
  }
  
  new_order = match(peakNm, rownames(yy1))
  
  ## plot the histone marker side by side with replicates
  df_histM = data.frame(time = sapply(colnames(yy1), function(x) unlist(strsplit(as.character(x), '_'))[2])
                        #markers = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[2]), stringsAsFactors = FALSE
  )
  colnames(df_histM) = c('time')
  rownames(df_histM) = colnames(yy1)
  sample_colors_histM =  sample_colors = c('springgreen', 'springgreen2', 'springgreen3', 'gold2', 'red')
  names(sample_colors_histM) = conds
  #marker_colors_histM = c('blue', 'red', 'deepskyblue2', 'darkgreen')
  #names(marker_colors_histM) = histMs
  annot_colors_histM = list(time = sample_colors_histM)
  gaps.col_histM = c(1:4)
  
  kk = c(1:5)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = colorRampPalette(((brewer.pal(n = 7, name ="BrBG"))))(10)
  p1 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                #color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlGnBu")))(8), 
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col= df_histM_new,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  kk = c(6:10)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da"))
  p2 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  
  kk = c(11:15)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = rev(terrain.colors(10))
  p3 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                annotation_legend = FALSE,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  plot_list[['p3']]=p3[[4]]
  
  pdf(paste0(figureDir, "/regeneration_chromatin_landscape_atac_histM.pdf"),
      width = 6, height = 10) # Open a new pdf file
  
  layout = matrix(c(1, 2, 3), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  dev.off()
  
  ### plot the atac-seq peaks with the changed order
  Replot.atac.peaks.with.newOrder = FALSE
  if(Replot.atac.peaks.with.newOrder){
    conds.atac = c("Embryo_Stage40", "Embryo_Stage44_proximal", 'Embryo_Stage44_distal', 
                   "Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", 'BL_UA_13days_distal'
    )
    df <- data.frame(conds.atac)
    rownames(df) = colnames(yy_atac)
    colnames(df) = 'time'
    
    sample_colors = c('magenta', 'darkblue', 'springgreen4', 'springgreen', 'springgreen2', 'springgreen3', 'gold2',
                      'red')[c(1:length(conds.atac))]
    names(sample_colors) = conds.atac
    
    # col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
    #           "#33a02c", "#fb9a99", "#e31a1c",
    #           "#fdbf6f", "#ff7f00", "#cab2d6",
    #           "#6a3d9a", "#ffff99", "#b15928")
    # cluster_col = col3[1:nb_clusters]
    # names(cluster_col) = paste0('cluster_', c(1:nb_clusters))
    annot_colors = list(
      condition = sample_colors)
    
    ii.gaps = c(3, 4)
    col = colorRampPalette(c("navy", "white", "red3"))(16)
    
    pheatmap(yy_atac, 
             #annotation_row = my_gene_col, 
             annotation_col = df, show_rownames = FALSE, scale = 'none', 
             color = col, 
             show_colnames = FALSE,
             cluster_rows = FALSE, cluster_cols = FALSE,  
             #clustering_method = 'complete', cutree_rows = nb_clusters, 
             annotation_colors = annot_colors,
             gaps_row = gaps.row,
             #clustering_callback = callback,
             gaps_col = ii.gaps, 
             filename = paste0(figureDir, 'heatmap_regenerationPeaks_scaled.pdf'), 
             width = 6, height = 12)
  }
}

##########################################
# segment-specific chromatin 
# and subclustering atac-peak clusters
##########################################
aggregate_atacPeaks_histMpeaks = function()
{
  ### collect the four markers with means 
  histMs =  c('H3K4me3','H3K27me3', 'H3K4me1')
  shm = c()
  
  for(n in 1:length(histMs))
  {
    # n = 1
    cat(n, ' -- ', histMs[n], '\n')
    cpm = readRDS(file = paste0(RdataHistM, '/fpm_bc_TMM_combat_', histMs[n], '_', version.analysis.HistM, '.rds'))
    cpm = cpm[ , grep(histMs[n], colnames(cpm))]
    #design.sel = readRDS(file = paste0(RdataHistM, '/design.sels_bc_TMM_combat_DBedgeRtest_', histMs[n], '_', 
    #                                   version.analysis.HistM, '.rds'))
    
    ### select the samples and extract sample means
    conds_histM = c("mUA", "mLA", "mHand")
    sample.sels = c();  
    cc = c()
    
    sample.means = c()
    for(ii in 1:length(conds_histM)) 
    {
      kk = grep(conds_histM[ii], colnames(cpm))
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds_histM[ii], length(kk)))
      if(length(kk)>1) {
        sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
      }else{
        sample.means = cbind(sample.means, cpm[, kk])
      }
    }
    colnames(sample.means) = paste0(conds_histM, '_', histMs[n])
    
    if(n == 1){
      shm = sample.means
    }else{
      shm = cbind(shm, sample.means[match(rownames(shm), rownames(sample.means)), ])
    }
  }    
  
  ## quantile normalization for different histone markers
  library(preprocessCore)
  xx = normalize.quantiles(shm)
  colnames(xx) = colnames(shm)
  rownames(xx) = rownames(shm)
  shm = xx
  
  pp_histM = data.frame(t(sapply(rownames(shm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp_histM$strand = '*'
  pp_histM = makeGRangesFromDataFrame(pp_histM, seqnames.field=c("X1"),
                                      start.field="X2", end.field="X3", strand.field="strand")
  
  mapping = findOverlaps(pp, pp_histM, type = 'within')
  
  yy0 = shm[mapping@to, ]
  #rownames(yy0) = rownames(yy)
  
  # reorder the histM according to the atac-seq peak clusters
  yy0 = yy0[plt$tree_row$order, ]
  yy0 = t(apply(yy0, 1, cal_z_score))
  #test = yy[plt$tree_row$order, ]
  #rownames(yy0) = rownames(yy)
  #yy0 = cbind(yy, yy0)
  
  df_histM = data.frame(segments = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[1])
                        #markers = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[2]), stringsAsFactors = FALSE
  )
  colnames(df_histM) = c('seg')
  rownames(df_histM) = colnames(yy0)
  sample_colors_histM = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors_histM) = c('mUA', 'mLA', 'mHand')
  marker_colors_histM = c('blue', 'red', 'deepskyblue2', 'darkgreen')
  names(marker_colors_histM) = histMs
  annot_colors_histM = list(seg = sample_colors_histM)
  gaps.col_histM = seq(1, ncol(yy0), by = 1)
  #clusters <- rownames(subsetDat[heat$tree_row[["order"]],])
  clusters <- sort(cutree(plt$tree_row, k = nb_clusters))
  
  # cluster order : 6, 1, 5, 3, 4, 2
  gaps.row = c(32, 32+76, 31 + 32 + 76, 292 + 31 + 32 + 76, 292 + 31 + 32 + 76 + 103)
  
  #yy0 = yy0[plt$tree_row$order, ]
  p1 = pheatmap(test, 
                nnotation_row = my_gene_col, 
                annotation_col = df, show_rownames = FALSE, scale = 'none', 
                color = col, 
                show_colnames = FALSE,
                cluster_rows = FALSE, cluster_cols = FALSE,  
                #clustering_method = 'complete', cutree_rows = nb_clusters, 
                annotation_colors = annot_colors, 
                #clustering_callback = callback,
                legend = TRUE,
                gaps_col = gaps.col, 
                gaps_row = gaps.row)
  
  p2 = pheatmap(yy0, cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  #grid.arrange(p1, p2,  nrow = 1)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  #plot_list[['p3']]=p2[[4]]
  grid.arrange(grobs=plot_list, ncol=2)
  
  indexs = data.frame(peak = names(clusters), cluster= clusters, stringsAsFactors = FALSE)
  indexs = indexs[match(rownames(test), rownames(indexs)), ]
  c_order = unique(indexs$cluster)
  
  yy1 = yy0
  
  ### transform the histM to account for different backgrounds and scaling
  # for(m in 1:ncol(yy1))
  # {
  #   jj1 = which(yy1[ ,m] < 1.5)
  #   yy1[jj1, m] = 1.5
  #   jj2 = which(yy1[ ,m] > 4)
  #   yy1[jj2, m] = 4
  #    
  # }
  # 
  
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  
  for(ac in c_order)
  {
    # ac = 6
    kk = which(indexs$cluster == ac)
    cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
    
    hm_hclust <- hclust(dist(as.matrix(yy1[kk,])), method = "complete")
    #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
    peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
    
  }
  
  new_order = match(peakNm, rownames(yy1))
  
  p1 = pheatmap(test[new_order, ], 
                nnotation_row = my_gene_col, 
                annotation_col = df, show_rownames = FALSE, scale = 'none', 
                color = col, 
                show_colnames = FALSE,
                cluster_rows = FALSE, cluster_cols = FALSE,  
                #clustering_method = 'complete', cutree_rows = nb_clusters, 
                annotation_colors = annot_colors, 
                #clustering_callback = callback,
                legend = TRUE,
                annotation_legend = FALSE,
                gaps_col = gaps.col, 
                gaps_row = gaps.row)
  
  gaps.col_histM = c(1:2)
  df_histM_new = as.data.frame(df_histM[c(1:3),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(1:3)]
  
  p2 = pheatmap(yy1[new_order, c(1:3)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col= df_histM_new,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  df_histM_new = as.data.frame(df_histM[c(4:6),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(4:6)]
  p3 = pheatmap(yy1[new_order, c(4:6)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  df_histM_new = as.data.frame(df_histM[c(7:9),])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[c(7:9)]
  p4 = pheatmap(yy1[new_order, c(7:9)], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(8), 
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                annotation_legend = FALSE,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  plot_list[['p3']]=p3[[4]]
  plot_list[['p4']]=p4[[4]]
  
  pdf(paste0(figureDir, "/Segemet_specific_chromatin_landscape_atac_histM_firstTry.pdf"),
      width = 10, height = 12) # Open a new pdf file
  
  layout = matrix(c(1, 1, 2, 3, 4), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  dev.off()
  
}

##########################################
# reorder the positional histM for postional atac-seq peaks 
##########################################
subsampling.postional.histM.postioinalAtacPeaks = function()
{
  ##########################################
  #  # reorder the histM according to the atac-seq peak clusters
  # cluster order : 6, 1, 5, 3, 4, 2
  # gaps.row = c(32, 32+76, 31 + 32 + 76, 292 + 31 + 32 + 76, 292 + 31 + 32 + 76 + 103)
  ##########################################
  peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rdata/', 
                                'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds'))
  table(peaks$clusters)
  cluster_order = c(6, 1, 5, 3, 4, 2)
  
  ## specify row gaps as atac-seq peaks
  gaps.row = c()
  for(n in 1:(length(cluster_order)-1))
  {
    if(n == 1)  {gaps.row = c(gaps.row, length(which(peaks$clusters == cluster_order[n])))
    }else{
      gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(peaks$clusters == cluster_order[n])))
    }
  }
  
  ## refine histone subclusters for each ATACseq peak cluster
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  
  for(ac in cluster_order)
  {
    # ac = 6
    kk = which(peaks$cluster == ac)
    cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
    
    hm_hclust <- hclust(dist(as.matrix(yy1[kk,])), method = "complete")
    #hm_cluster <- cutree(tree = as.dendrogram(hm_hclust), h = 5)
    peakNm = c(peakNm, hm_hclust$labels[hm_hclust$order])
  }
  
  new_order = match(peakNm, rownames(yy1))
  
  
  ## plot the histone marker side by side with replicates
  df_histM = data.frame(segments = sapply(colnames(yy1), function(x) unlist(strsplit(as.character(x), '_'))[2])
                        #markers = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[2]), stringsAsFactors = FALSE
  )
  colnames(df_histM) = c('seg')
  rownames(df_histM) = colnames(yy1)
  sample_colors_histM = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors_histM) = c('mUA', 'mLA', 'mHand')
  #marker_colors_histM = c('blue', 'red', 'deepskyblue2', 'darkgreen')
  #names(marker_colors_histM) = histMs
  annot_colors_histM = list(seg = sample_colors_histM)
  gaps.col_histM = c(2, 4)
  
  kk = c(1:6)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = colorRampPalette((brewer.pal(n = 7, name ="BrBG")))(10)
  p1 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                #color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlGnBu")))(8), 
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col= df_histM_new,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  kk = c(7:12)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  
  cols = rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da"))
  p2 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  
  kk = c(13:18)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  cols = rev(terrain.colors(10))
  p3 = pheatmap(yy1[new_order, kk], cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                annotation_legend = FALSE,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  
  plot_list=list()
  plot_list[['p1']]=p1[[4]]
  plot_list[['p2']]=p2[[4]]
  plot_list[['p3']]=p3[[4]]
  
  pdf(paste0(figureDir, "/Segemet_specific_chromatin_landscape_atac_histM.pdf"),
      width = 8, height = 10) # Open a new pdf file
  
  layout = matrix(c(1, 2, 3), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  dev.off()
}

