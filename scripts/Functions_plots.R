##########################################################################
##########################################################################
# Project: positional memory 
# Script purpose: functions to make plots
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jun 28 10:48:37 2022
##########################################################################
##########################################################################

##########################################
# reorder the positional histM for postional atac-seq peaks 
##########################################
subclustering.postional.histM.postioinalAtacPeaks = function(test)
{
  ##########################################
  #  # reorder the histM according to the atac-seq peak clusters
  # cluster order : 6, 1, 5, 3, 4, 2
  # gaps.row = c(32, 32+76, 31 + 32 + 76, 292 + 31 + 32 + 76, 292 + 31 + 32 + 76 + 103)
  ##########################################
  peaks = readRDS(file = paste0('~/workspace/imp/positional_memory/results/Rdata/', 
                                'position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_with.clusters_6.rds')) # positional atac
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
  
  ## center the histone signals
  yy1 = readRDS(file = paste0(RdataDir, '/peak_signals_atac_4histM_positionalPeaks.rds'))
  yy1 = yy1[, grep('mRep', colnames(yy1))]
  test = yy1
  x = c()
  for(ht in c("H3K4me3", "H3K27me3", "H3K4me1"))
  {  # ht = "H3K4me3"
    cat(ht, '\n')
    ii.test = grep(ht, colnames(test))
    ttest = test[, ii.test]
    
    fc_sels = c('mUA', 'mLA', 'mHand')
    jj.test = c()
    for(fc in fc_sels) jj.test = c(jj.test, grep(fc, colnames(ttest)))
    ttest = ttest[, jj.test]
    #colnames(test) = fc_sels
    
    ttest =  cal_sample_means(ttest, conds = c('mUA', 'mLA', 'mHand'))
    ttest = t(apply(ttest, 1, cal_centering))
    colnames(ttest) = paste0(ht, '_', fc_sels)
    
    x = cbind(x, ttest)
  }
  test = x; rm(x)
  
  ## refine histone subclusters for each ATACseq peak cluster
  peakNm = c()
  library(dendextend)
  library(ggplot2)
  library(khroma)
  
  ## try dendsort methods
  library("dendsort")
  library("seriation")
  
  callback = function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  #new_order = c()
  for(ac in cluster_order)
  {
    # ac = 2
    kk = which(peaks$cluster == ac)
    cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
    
    if(ac == 6){
      ac = 6
      kk = which(peaks$cluster == ac)
      cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
      x <- as.matrix(test[kk,]) #drop the 5th colum
      nb_clusters = 3
      px = pheatmap(x[, c(4:9, 1:3)], #annotation_row = my_gene_col, 
                    #annotation_col = df, 
                    show_rownames = FALSE, scale = 'none', 
                    #color = col, 
                    show_colnames = FALSE,
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE,  
                    clustering_callback = callback,
                    clustering_method = 'complete', cutree_rows = nb_clusters, 
                    #annotation_colors = annot_colors, 
                    gaps_col = c(3, 5))
      peakNm = c(peakNm, px$tree_row$labels[px$tree_row$order])
      
    }
    if(ac == 1){
      ac = 1
      kk = which(peaks$cluster == ac)
      cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
      
      x <- as.matrix(test[kk,]) #drop the 5th colum
      nb_clusters = 4
      my_hclust_gene <- hclust(dist(x), method = "complete")
      my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
      
      px = pheatmap(x[, c(4:9, 1:3)], #annotation_row = my_gene_col, 
                    #annotation_col = df, 
                    show_rownames = FALSE, scale = 'none', 
                    #color = col, 
                    show_colnames = FALSE,
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE,  
                    clustering_callback = callback,
                    clustering_method = 'complete', cutree_rows = nb_clusters, 
                    #annotation_colors = annot_colors, 
                    gaps_col = c(3, 5))
      source('Functions_utility.R')
      gname = sort_geneNames_desiredClusterRank_c1(px, my_gene_col)
      peakNm = c(peakNm, gname)
      
      saveRDS(peakNm, file = paste0(RdataDir, '/peak_names_manual_order_c6.c1.rds'))
    }
    
    if(ac == 5){
      ac = 5
      kk = which(peaks$cluster == ac)
      cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
      x <- as.matrix(test[kk,]) #drop the 5th colum
      
      nb_clusters = 4
      my_hclust_gene <- hclust(dist(x), method = "complete")
      my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
      
      px = pheatmap(x[, c(4:9, 1:3)], #annotation_row = my_gene_col, 
                    #annotation_col = df, 
                    show_rownames = FALSE, scale = 'none', 
                    #color = col, 
                    show_colnames = FALSE,
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE,  
                    clustering_callback = callback,
                    clustering_method = 'complete', cutree_rows = nb_clusters, 
                    #annotation_colors = annot_colors, 
                    gaps_col = c(3, 5))
      
      peakNm = c(peakNm, px$tree_row$labels[px$tree_row$order])
      
      saveRDS(peakNm, file = paste0(RdataDir, '/peak_names_manual_order_c6.c1.c5.rds'))
    }
    
    if(ac == 3){
      ac = 3
      kk = which(peaks$cluster == ac)
      cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
      
      x <- as.matrix(test[kk,]) #drop the 5th colum
      nb_clusters = 6
      my_hclust_gene <- hclust(dist(x), method = "complete")
      my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
      
      px = pheatmap(x[, c(4:9, 1:3)], #annotation_row = my_gene_col, 
                    #annotation_col = df, 
                    show_rownames = FALSE, scale = 'none', 
                    #color = col, 
                    show_colnames = FALSE,
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE,  
                    clustering_callback = callback,
                    clustering_method = 'complete', cutree_rows = nb_clusters, 
                    #annotation_colors = annot_colors, 
                    gaps_col = c(3, 5))
      
      source('Functions_utility.R')
      gname = sort_geneNames_desiredClusterRank_c3(px, my_gene_col)
      peakNm = c(peakNm, gname)
      
      saveRDS(peakNm, file = paste0(RdataDir, '/peak_names_manual_order_c6.c1.c5.c3.rds'))
    }
    
    if(ac == 4){
      ac = 4
      kk = which(peaks$cluster == ac)
      cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
      
      x <- as.matrix(test[kk,]) #drop the 5th colum
      nb_clusters = 6
      my_hclust_gene <- hclust(dist(x), method = "complete")
      my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
      
      px = pheatmap(x[, c(4:9, 1:3)], #annotation_row = my_gene_col, 
                    #annotation_col = df, 
                    show_rownames = FALSE, scale = 'none', 
                    #color = col, 
                    show_colnames = FALSE,
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE,  
                    clustering_callback = callback,
                    clustering_method = 'complete', cutree_rows = nb_clusters, 
                    #annotation_colors = annot_colors, 
                    gaps_col = c(3, 5))
      
      source('Functions_utility.R')
      gname = sort_geneNames_desiredClusterRank_c4(px, my_gene_col)
      peakNm = c(peakNm, gname)
      
      
      saveRDS(peakNm, file = paste0(RdataDir, '/peak_names_manual_order_c6.c1.c5.c3.c4.rds'))
      
    }
    
    if(ac == 2){
      ac = 2
      kk = which(peaks$cluster == ac)
      cat('clsuter ', ac, ' -- ', length(kk), ' peaks \n')
      
      x <- as.matrix(test[kk,]) #drop the 5th colum
      nb_clusters = 6
      my_hclust_gene <- hclust(dist(x), method = "complete")
      my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
      
      px = pheatmap(x[, c(4:9, 1:3)], #annotation_row = my_gene_col, 
                    #annotation_col = df, 
                    show_rownames = FALSE, scale = 'none', 
                    #color = col, 
                    show_colnames = FALSE,
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE,  
                    clustering_callback = callback,
                    clustering_method = 'complete', cutree_rows = nb_clusters, 
                    #annotation_colors = annot_colors, 
                    gaps_col = c(3, 5))
      source('Functions_utility.R')
      gname = sort_geneNames_desiredClusterRank_c2(px, my_gene_col)
      peakNm = c(peakNm, gname)
      
      saveRDS(peakNm, file = paste0(RdataDir, '/peak_names_manual_order_c6.c1.c5.c3.c4.c2.rds'))
      
    }
  }
  
  new_order = match(peakNm, rownames(test))
  
  ## plot the histone marker side by side with replicates
  df_histM = data.frame(segments = sapply(colnames(test), function(x) unlist(strsplit(as.character(x), '_'))[2])
                        #markers = sapply(colnames(yy0), function(x) unlist(strsplit(as.character(x), '_'))[2]), stringsAsFactors = FALSE
  )
  colnames(df_histM) = c('seg')
  rownames(df_histM) = colnames(test)
  sample_colors_histM = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors_histM) = c('mUA', 'mLA', 'mHand')
  #marker_colors_histM = c('blue', 'red', 'deepskyblue2', 'darkgreen')
  #names(marker_colors_histM) = histMs
  annot_colors_histM = list(seg = sample_colors_histM)
  
  kk = c(1:3)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  nb_breaks = 8
  sunset <- colour("sunset")
  PRGn <- colour("PRGn")
  range <- 1.5
  xx = test[new_order, kk]
  xx = t(apply(xx, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
  cols = (sunset(nb_breaks))
  gaps.col_histM = c(1, 2)
  p1 = pheatmap(xx, cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                #color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlGnBu")))(8), 
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                cluster_cols=FALSE, annotation_col= df_histM_new,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM,
                breaks = seq(-range, range, length.out = nb_breaks), 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  kk = c(4:6)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  range <- 2.0
  xx = test[new_order, kk]
  xx = t(apply(xx, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
  
  cols = cols = rev(PRGn(nb_breaks-1));
  p2 = pheatmap(xx, cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                breaks = seq(-range, range, length.out = nb_breaks), 
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                gaps_col = gaps.col_histM, 
                annotation_legend = FALSE,
                gaps_row = gaps.row)
  
  
  kk = c(7:9)
  df_histM_new = as.data.frame(df_histM[kk,])
  colnames(df_histM_new) = colnames(df_histM)
  rownames(df_histM_new) = rownames(df_histM)[kk]
  
  xx = test[new_order, kk]
  range <- 2.0
  xx = t(apply(xx, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
  
  cols = cols = colorRampPalette(rev((brewer.pal(n = 8, name ="BrBG"))))(7)
  p3 = pheatmap(xx, cluster_rows=FALSE, show_rownames=FALSE, fontsize_row = 5,
                color = cols,
                show_colnames = FALSE,
                scale = 'none',
                breaks = seq(-range, range, length.out = nb_breaks), 
                cluster_cols=FALSE, annotation_col=df_histM,
                annotation_colors = annot_colors_histM,
                annotation_legend = FALSE,
                gaps_col = gaps.col_histM, 
                gaps_row = gaps.row)
  
  
  
  plot_list=list(); plot_list[['p1']]=p1[[4]]; plot_list[['p2']]=p2[[4]]; plot_list[['p3']]=p3[[4]]; 
  layout = matrix(c(1, 2, 3), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  pdf(paste0(figureDir, "/positional_chromatin_landscape_3histM_clusteringAllHistM_manual.pdf"),
      width = 6, height = 10) # Open a new pdf file
  
  layout = matrix(c(1, 2, 3), nrow = 1)
  grid.arrange(grobs=plot_list, nrow= 1,
               layout_matrix = layout)
  
  dev.off()
  
  ### reorder positional atacpeak acccordingly
  yy = readRDS(paste0(RdataDir, '/positional_atacPeaks_data_3reps_forHeatmap.rds'))
  df = data.frame(sapply(colnames(yy), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[2])}))
  colnames(df) = 'seg'
  rownames(df) = colnames(yy)
  sample_colors = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors) = c('UA', 'LA', 'Hand')
  annot_colors = list(
    seg = sample_colors)
  gaps.col = c(3, 6)
  col = colorRampPalette(c("navy", "white", "red3"))(8)
  
  p1 = pheatmap(yy[new_order, ], #annotation_row = my_gene_col, 
                annotation_col = df, show_rownames = FALSE, scale = 'none', 
                color = col, 
                show_colnames = FALSE,
                cluster_rows = FALSE, 
                cluster_cols = FALSE,  
                #clustering_method = 'complete', cutree_rows = nb_clusters, 
                annotation_colors = annot_colors, 
                #clustering_callback = callback,
                gaps_col = gaps.col, 
                treeheight_row = 20,
                gaps_row = gaps.row,
                annotation_legend = FALSE,
                filename = paste0(figureDir, '/positional_atacPeaks_fdr0.05_log2FC.1_rmHeadPeaks_reorder_final.pdf'), 
                width = 4, height = 12)
  
  ##########################################
  # old version of histone mark heatmap without subclustering
  ##########################################
  # ## heatmaps of histM for segment-specific atacPeaks (centered signals)
  # ## simiar plot to regeneration log2FC.sample.vs.mUA
  # table(peaks$clusters)
  # cluster_order = c(6, 1, 5, 3, 4, 2)
  # gaps.row = c() 
  # for(n in 1:(length(cluster_order)-1)) ## compute row gaps as atac-seq peaks
  # {
  #   if(n == 1)  {gaps.row = c(gaps.row, length(which(peaks$clusters == cluster_order[n])))
  #   }else{
  #     gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(peaks$clusters == cluster_order[n])))
  #   }
  # }
  # 
  # plt = readRDS(file = paste0(RdataDir, '/postional_atacPeaks_heatmap_orderSaved.rds'))
  # 
  # for(n in 1:3)
  # {
  #   # n = 2
  #   ii.test = grep(conds_histM[n], colnames(yy1))
  #   test = yy1[, ii.test]
  #   fc_sels = c('mUA', 'mLA', 'mHand')
  #   jj.test = c()
  #   for(fc in fc_sels) jj.test = c(jj.test, grep(fc, colnames(test)))
  #   test = test[, jj.test]
  #   #colnames(test) = fc_sels
  #   
  #   test =  cal_sample_means(test, conds = c('mUA', 'mLA', 'mHand'))
  #   test = t(apply(test, 1, cal_centering))
  #   #yy1[ ,jj] = t(apply(yy1[,jj], 1, cal_transform_histM, cutoff.min = 0, cutoff.max = 5, centering = FALSE, toScale = TRUE))
  #   
  #   range <- 2.0
  #   test = t(apply(test, 1, function(x) {x[which(x >= range)] = range; x[which(x<= (-range))] = -range; x}))
  #   
  #   nb_breaks = 8
  #   sunset <- colour("sunset")
  #   PRGn <- colour("PRGn")
  #   #highcontrast <- colour("high contrast")
  #   if(n == 1) cols = (sunset(nb_breaks))
  #   if(n == 2) cols = rev(PRGn(nb_breaks-1)); # paletteer_d("colorBlindness::Blue2DarkRed18Steps") 
  #   if(n == 3)   cols = colorRampPalette(rev((brewer.pal(n = 8, name ="BrBG"))))(7)
  #   #cols =  colorRampPalette(colour("high contrast"))(nb_breaks)
  #   #cols = rev(terrain.colors(10))
  #   
  #   pheatmap(test[plt$tree_row$order, ], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
  #            #color = c('darkgray', 'blue'), 
  #            #color = colorRampPalette((brewer.pal(n = 7, name ="PRGn")))(nb_breaks),
  #            color = cols, 
  #            breaks = seq(-range, range, length.out = nb_breaks), 
  #            gaps_row = gaps.row,
  #            filename = paste0(figureDir, '/positional_histM_centerend_1246positionalATACpeaks_', conds_histM[n], '.pdf'), 
  #            width = 3, height = 10)
  #}
  
  
}

##########################################
# process heatmap output from deeptools  
##########################################
Process.deeptools.heatmapTable = function()
{
  dpt = read.table(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/', 
                                 'Data/atacseq_using/heatmap_deeptools/heatmaps_all/heatmap_allchrFeatures_matrix4save.txt'),
                   comment.char = "#", sep = '\t', header = TRUE, skip = 2)
  
  peaks = c(rep(1, 76), 
            rep(2, 712), 
            rep(3, 292),
            rep(4, 103), 
            rep(5, 31), 
            rep(6, 32))
  
  samples = colnames(dpt)
  samples = sapply(samples, function(x){x = unlist(strsplit(as.character(x), '_')); paste0(x[1:3], collapse = '_')})          
  
  cluster_order = c(6, 1, 5, 3, 4, 2)
  
  index = c()
  gaps.row = c() 
  for(n in 1:length(cluster_order)) ## compute row gaps as atac-seq peaks
  {
    jj = which(peaks == cluster_order[n])
    index = c(index, jj)
    if(n == 1)  {gaps.row = c(gaps.row, length(jj))
    }else{
      if(n < length(cluster_order)){
        gaps.row = c(gaps.row,  gaps.row[n-1] + length(which(peaks == cluster_order[n])))
      }
    }
  }
  
  dpt = dpt[index, ]
  
  ##########################################
  # select representative samples 
  ##########################################
  library(ComplexHeatmap)
  library("RColorBrewer")
  library("circlize")
  
  ## select the sample to use for plotting
  idx =  c(grep('136164', samples), 
           grep('161521', samples), 
           grep('74940', samples), 
           intersect(c(grep('H3K27me3', samples), grep('H3K4me1', samples), grep('H3K4me3', samples)), grep('mRep2', samples))
  )
  samples = samples[idx]
  dpt = dpt[, idx]
  
  ## sort by H3K27me3 for each cluster
  index = grep('H3K27me3_mUA_mRep2', samples)
  ss = apply(as.matrix(dpt[, index]), 1, mean)
  new_order = c()
  for(c in cluster_order)
  {
    kk = which(peaks==c)
    new_order = c(new_order, kk[order(-ss[kk])])
  }
  
  features = c('atac','H3K4me3', 'H3K27me3', 'H3K4me1')
  
  for(n in 1:length(features))
  {
    # n = 4
    cat(n, ' -- ', features[n], '\n')
    
    index = grep(features[n], samples)
    test = as.matrix(dpt[new_order, index])*20
    
    splits = factor(peaks[new_order], levels = cluster_order)
    sample.sel = samples[index]
    sample.sel = sapply(sample.sel, function(x) unlist(strsplit(as.character(x), '_'))[2])
    ha <- HeatmapAnnotation(samples = sample.sel)
    
    unique(sample.sel)
    
    #col_fun = colorRamp2(c(0, 0.5, 1), c("#377EB8", "white", "#E41A1C"))
    quantile(test, c(0, 0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 0.8, 0.85,  0.90,  0.95, 0.99, 1))
    
    if(features[n] == 'atac') {range = 2.5; 
    cols = colorRamp2(c(0, 0.5, 1.5, range), rev(brewer.pal(n=4, name="RdBu")))
    #cols = colorRamp2(c(0, 2, 4), c("blue", "white", "red"));
    }
    if(features[n] == 'H3K27me3') {range = 5; cols = colorRamp2(c(0, range), colors = c('white', 'purple1'))}
    if(features[n] == 'H3K4me1'){range = 7; cols = colorRamp2(c(0, range), colors = c('white', 'orange2'))}
    if(features[n] == 'H3K4me3'){range = 5; cols = colorRamp2(c(0, range), c("white", "green"))
    #cols = colorRamp2(c(0.2, range), c("white", "#3794bf"))
    }
    
    pdf(paste0(figureDir, "/positional_peaks_intensity_heatmap_", features[n], "_v2.pdf"),
        width = 4, height = 10) # Open a new pdf file
    Heatmap(test, 
            cluster_rows = FALSE,
            cluster_columns = FALSE, 
            show_row_names = FALSE,
            show_column_names = FALSE,
            row_split = splits,
            cluster_column_slices = FALSE,
            column_split = factor(sample.sel, levels = c('mUA', 'mLA','mHand')),
            top_annotation = ha,
            show_heatmap_legend = TRUE,
            #col = colorRamp2(seq(0, range, length.out = breaks), rev(brewer.pal(n=breaks, name="RdBu")))
            col = cols
            #name = "mtcars", #title of legend
            #column_title = "Variables", row_title = "Samples",
            #row_names_gp = gpar(fontsize = 7) # Text size for row names
    )
    
    dev.off()
    
  }
  
}
