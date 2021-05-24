##########################################################################
##########################################################################
# Project: MARA motif activity analysis
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed May 19 10:47:23 2021
##########################################################################
##########################################################################
extract.TFs.annotation.from.TFClass = function()
{
  library(rdflib)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(jsonld)
  
  tfs = '/Users/jiwang/workspace/imp/positional_memory/results/motif_analysis/tfclass.ttl'
  
  x <- rdf()
  rdf <- rdf_parse(tfs)
  
  options(rdf_print_format = "turtle")
  rdf
  
  car_triples <- 
    rdf %>% 
    rownames_to_column("label") %>% 
    gather(attribute,measurement, -Model)
  
}


# manually add extra motifs for hnd-1, pha-4, unc-120 and nhr-67 from dm, mus and homo
convert.cisbp.format.to.meme = function()
{
  library(universalmotif)
  #pwmDIr = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/extra_pwm'
  pwm.old = '../results/motif_analysis/SwissRegulon_PWMs/hg19_weight_matrices_v2.txt'
  xx = read_matrix(file = pwm.old, skip = 0)
  yy = convert_type(xx, "PWM")
  
  xx <- query(MotifDb, andStrings=c("hsapiens"), orStrings=c("swissregulon"))
  yy = convert_type(xx, "PWM")
  
  write_meme(yy, file = '../results/motif_analysis/SwissRegulon_PWMs/hg19_weight_matrices_v2.meme', overwrite = TRUE)
  
}

extract.SwissRegulon.meme.from.MotifDb = function()
{
  library(universalmotif)
  #pwmDIr = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/extra_pwm'
  #pwm.old = '../results/motif_analysis/SwissRegulon_PWMs/hg19_weight_matrices_v2.txt'
  #xx = read_matrix(file = pwm.old, skip = 0)
  #yy = convert_type(xx, "PWM")
  
  library(MotifDb)
  xx <- query(MotifDb, andStrings=c("hsapiens"), orStrings=c("swissregulon"))
  yy = convert_type(xx, "PWM")
  
  write_meme(yy, file = '../results/motif_analysis/SwissRegulon_PWMs/hg19_weight_matrices_v2.meme', overwrite = TRUE)
  
  
  
  
}


generate.logos.for.motifs.pwm = function()
{
  library('universalmotif')
  pwm = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015_curated_extra.meme'
  #pwm = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.meme'
  
  meme= read_meme(file = pwm, skip = 0, readsites = FALSE, readsites.meta = FALSE)
  motifs = convert_motifs(meme, class = "universalmotif-universalmotif")
  # comparisons <- compare_motifs(motifs, method = "PCC", min.mean.ic = 0,
  #                               score.strat = "a.mean")
  # write.table(comparisons, file = '../data/motifs_tfs/pwm_similarity_correction_PCC.txt', sep = '\t', col.names = TRUE, 
  #             row.names = TRUE, quote = FALSE)
  
  p1 = view_motifs(motifs[[1]], use.type = 'ICM')
  plot(p1)
  
  for(n in 1:length(motifs))
  {
    cat(n, '\n')
    pdfname = paste0('../data/motifs_tfs/pwm_logos/', motifs[[n]]@name, '.pdf')
    pdf(pdfname, width=8, height = 6)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    p1 = view_motifs(motifs[[n]], use.type = 'ICM')
    plot(p1)
    
    dev.off()
  }
  
}


########################################################
########################################################
# Section : motif activity analysis 
# 
########################################################
########################################################
save.peak.bed.file.for.fimo = function()
{
  load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.withinPeaks.pval6.Rdata'))
  fpm = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
  
  # prepare the background distribution
  jj = grep('bg_', rownames(fpm), invert = TRUE)
  fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
  fpm = fpm[jj, ]
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  
  pp$name = rownames(fpm)
  pp$strand = '*'
  
  write.table(pp, file = '../results/motif_analysis/peaks/peaks_for_fimo.bed', row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = '\t')
  
}


##########################################
# after running FIMO, make motif occurrency matrix  
##########################################
make.motif.oc.matrix.from.fimo.output = function()
{
  library(data.table)
  #motif.tf = readRDS( '../data/motifs_tfs/motif_tf_mapping.rds')
  fimo.out = '../results/motif_analysis/FIMO/fimo_out/fimo.tsv'
  fimo = fread(fimo.out, header = TRUE)
  motif.oc = table(fimo$motif_id, fimo$sequence_name, useNA = 'ifany')
  motif.oc = t(motif.oc)
  
  print(head(rownames(motif.oc)))
  
  ##########################################
  # associate the scanned regions with gene
  # here we restrict the assignment to protein-coding genes using ChIPpeakAnno 
  # https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/pipeline.html
  ##########################################
  assign.regions.to.genes = FALSE
  if(assign.regions.to.genes){
    ## loading packages
    #library(ChIPseeker)
    library(ChIPpeakAnno)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(plyranges)
    
    annot = read.csv(file = '/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235_noFilters.csv') # all annotation included
    annot = annot[grep('protein_coding', annot$Gene.type), ] # keep only protein-coding genes
    #library(TxDb.Celegans.UCSC.ce11.ensGene)
    #txdb <- TxDb.Celegans.UCSC.ce11.ensGene
    annot = annot[match(unique(annot$Gene.name), annot$Gene.name), ]
    xx = data.frame(seqnames = paste0('chr', annot$Chromosome.scaffold.name), start = annot$Gene.start..bp., end= annot$Gene.end..bp., 
                    gene_id = as.character(annot$Gene.name),
                    strand=c("."), score=0, stringsAsFactors = FALSE)
    xx = makeGRangesFromDataFrame(xx, keep.extra.columns = TRUE)
    names(xx) = as.character(annot$Gene.name)
    
    get.peak.coord = function(x){
      x = unlist(strsplit(as.character(x), '[:]'))
      chr = x[1]
      x = unlist(strsplit(as.character(x[2]), '-'))
      return(c(chr, x))
    }
    peaks = rownames(motif.oc)
    peaks = t(sapply(peaks, get.peak.coord))
    colnames(peaks) = c('chr', 'start', 'end')
    peaks = data.frame(peaks, strand=c("."), score=0, stringsAsFactors = FALSE)
    peaks = makeGRangesFromDataFrame(peaks)
    #peaks = data.frame(, gsub('^*:')))
    #peak <- readPeakFile(peakfile = peak.file, as = 'GRanges')
    peakAnno = annotatePeakInBatch(peaks, 
                                   AnnotationData=xx, 
                                   output='nearestLocation', 
                                   bindingRegion=c(-2000, 500))
    #peakAnno <- annotatePeak(peak = peaks, tssRegion=c(-2000, 2000), level = 'gene', TxDb=txdb)
    assign = as.data.frame(peakAnno)
    
    #mm = match(assign$geneId, annot$wormbase.id)
    assign$genes = assign$feature
    rownames(assign) = assign$peak
    
    # keep scanned regions with mapped protein coding genes
    jj = which(!is.na(assign$genes))
    assign = assign[jj, ]
    motif.oc = motif.oc[jj, ]
    
    gene.uniq = unique(assign$genes)
    cat(nrow(motif.oc),  'scanned regions assigned to :', length(gene.uniq), ' genes \n')
    
    mocc = motif.oc[match(gene.uniq, assign$genes), ]
    rownames(mocc) = gene.uniq
    
    ## motif occurrency gene * motif
    for(n in 1:nrow(mocc))
    {
      kk = which(assign$genes == rownames(mocc)[n])
      if(length(kk)>1) {
        cat(n, '\n')
        mocc[n, ] = apply(motif.oc[kk, ], 2, sum) 
      }
    }
    
    motif.oc = mocc
    remove(mocc)
    
  }
  # else{
  #   load(file = '../data/Hashimsholy_et_al/annotMapping_ensID_Wormbase_GeneName.Rdata')
  #   
  #   mm = match(rownames(motif.oc), geneMapping$Wormbase)
  #   rownames(motif.oc) = geneMapping$Gene.name[mm]
  #   #kk = match(rownames(motif.oc, ))
  #   
  #   ss1 = apply(motif.oc, 1, sum)
  #   cat(length(which(ss1 == 0)), 'genes without scanned motifs \n')
  #   ss2 = apply(motif.oc, 2, sum)
  #   
  # }
  
  ##########################################
  # remove motif redundancy by merging occurrence for motifs in the same cluster
  ##########################################
  Remove.Motif.Redundancy = FALSE
  if(Remove.Motif.Redundancy){
    #mm = match(colnames(motif.oc), motif.tf$motifs)
    #colnames(motif.oc) = motif.tf$motifs.new[mm]
    names = unique(motif.tf$names[match(colnames(motif.oc), motif.tf$motifs)])
    xx = matrix(0, ncol = length(names), nrow = nrow(motif.oc))
    rownames(xx) = rownames(motif.oc)
    colnames(xx) = names
    
    ## to get motif occurrency gene * non-redundant-motif
    ## among the redundant motifs in the same cluster, the one with median of total occurrence is chosen.
    for(n in 1:ncol(xx))
    {
      # n = 3
      cat(n, ' -- ', colnames(xx)[n],  '\n')
      mtf = motif.tf$motifs[which(motif.tf$names == colnames(xx)[n])]
      
      kk = match(mtf, colnames(motif.oc))
      kk = kk[!is.na(kk)]
      
      if(length(kk) == 0){
        cat('Error : no motif found \n')
        
      }else{
        if(length(kk) == 1){
          
          xx[,n] = motif.oc[, kk]
          
        }else{
          cat('>>>>>>>>>',  length(kk), 'columns found \n')
          ss = apply(motif.oc[,kk], 2, sum)
          ss.o = ss[order(ss)]
          kk.sel = kk[which(ss == ss.o[ceiling(length(ss)/2)])]
          xx[,n] = motif.oc[, kk.sel[1]]
        }
      }
    }
    
    motif.oc = xx;
    remove(xx)
  }
  
  saveRDS(motif.oc, file = '../results/motif_analysis/motif_oc_fimo.rds')
  
}

##########################################
# run glmnet
##########################################
run.MARA.atac.temporal = function(keep, cc)
{
  library(pheatmap)
  library(RColorBrewer)
  require(glmnet)
  
  # prepare Y matrix 
  cc.uniq = unique(cc)
  Y = matrix(NA, ncol = length(cc.uniq), nrow = nrow(keep))
  
  for(n in 1:ncol(Y))
  {
    jj = which(cc == cc.uniq[n])
    if(length(jj) == 1) {
      Y[,n] = keep[,jj]
    }else{
      Y[,n] = apply(keep[ ,jj], 1, mean)
    }
  }
  colnames(Y) = cc.uniq
  rownames(Y) = rownames(keep)
  
  # X matrix
  motif.oc = readRDS(file = '../results/motif_analysis/motif_oc_fimo.rds')
  mm = match(rownames(motif.oc), rownames(Y))
  motif.oc = motif.oc[!is.na(mm), ]
  
  ss.m = apply(motif.oc, 2, sum)
  motif.oc = motif.oc[ , which(ss.m>0)]
  
  ss.p = apply(motif.oc, 1, sum)
  motif.oc = motif.oc[which(ss.p>0), ]
  
  kk = match(rownames(motif.oc), rownames(Y))
  Y = Y[kk, ]
  
  X = as.matrix(motif.oc)
  Y = as.matrix(Y)
  
  alpha = 0.5
  standardize = TRUE;
  use.lambda.min = TRUE;
  binarize.x = TRUE
  standardize.response=FALSE
  intercept=FALSE
  family = 'mgaussian'
  lambdas = 10^(seq(-4, 1, length.out = 100))
  x = X;
  y = scale(Y, center = TRUE, scale = FALSE); # center response and standardize X (cf. elastic-net paper)
  if(binarize.x) x = x > 0
  
  library(doMC) 
  registerDoMC(cores=6)
  
  library(tictoc)
  tic()
  cv.fit=cv.glmnet(x, y, family= family,
                   alpha=alpha, lambda = lambdas, standardize=standardize, 
                   standardize.response=standardize.response, parallel = TRUE)
  
  plot(cv.fit)
  toc()
  
  tic()
  fit=glmnet(x,y, alpha=alpha, lambda=cv.fit$lambda, family=family,
             standardize=standardize, standardize.response=standardize.response, intercept=intercept, 
             relax = FALSE)
  plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")
  toc()
    
  if(use.lambda.min){
    s.optimal = cv.fit$lambda.min
  }else{
    s.optimal = cv.fit$lambda.1se
  }
  
  xx = coef(fit, s = fit$lambda[100])
  aa = c()
  for(j in 1:length(xx))
  {
    aa = cbind(aa, as.numeric(xx[j]))
    
  }
  
  rownames(aa) = rownames(xx[[1]])
  #colnames(aa) = c('E40', 'E44.P', 'mUA', 'BL.UA.D5', 'BL.UA.D9', 'BL.UA.D13.P')
  colnames(aa) = colnames(y)
  aa = as.data.frame(aa[-1, ]) # ignore the intercept
  aa = apply(aa, 2, scale)
  #aa = scale(aa)
  rownames(aa) = rownames(xx[[1]])[-1]
  
  cutoff.activity = 0
  ss = apply(aa, 1, function(x) length(which(abs(x) > cutoff.activity)))
  print(aa[which(ss>0), ])
 
  
  pheatmap(aa[which(ss>0), ], cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
           scale = 'none', cluster_cols=FALSE, main = paste0("motif activity by MARA"), 
           na_col = "white", fontsize_col = 12) 
  
  
  # aa = as.data.frame(coef.glmnet(fit, s = s.optimal))
  # aa = aa[-1, ] # remove intecept
  # colnames(aa) = names(fit$beta)
  # if(alpha > 0.0){
  #   rownames(aa) = rownames(fit$beta[[2]])
  #   res = aa;
  #   ## collect result from the elastic-net
  #   #kk = apply(aa, 1, function(x) !all(x==0))
  #   #aa = aa[kk, ]
  #   #colnames(x)[which(fit$beta[[1]][,optimal]!=0)]
  #   #colnames(x)[which(fit$beta[[2]][,optimal]!=0)]
  #   
  #   # rerun lm with selected features
  #   relax.fitting.lm = FALSE
  #   if(relax.fitting.lm){
  #     fit.lm = lm(y ~ x[, match(rownames(aa), colnames(x))])
  #     res = data.frame(fit.lm$coefficients)
  #     res = res[-1, ] # remove intercept
  #     rownames(res) = rownames(aa)
  #     
  #     pheatmap(res, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
  #              scale = 'column', cluster_cols=FALSE, main = paste0("motif activity"), 
  #              na_col = "white", fontsize_col = 10)
  #   }
  #   
  # }else{
  #   if(zscore.output) aa = apply(aa, 2, scale)
  #   rownames(aa) = rownames(fit$beta[[2]])
  #   #ss = apply(aa, 1, function(x) !all(x==0))
  #   #aa = aa[ss, ]
  #   #head(rownames(aa)[order(-abs(aa$MSxp))], 10)
  #   #head(rownames(aa)[order(-abs(aa$MSxa))], 10)
  #   res = aa
  #   
  #   if(Test){
  #     ss = apply(aa, 1, function(x) length(which(abs(x) > Test.zscore.cutoff)))
  #     aa = aa[which(ss>0), ]
  #     print(aa)
  #   }
  # }
  
}


run.MARA.atac.spatial = function(keep, cc)
{
  library(pheatmap)
  library(RColorBrewer)
  require(glmnet)
  
  # prepare Y matrix 
  cc.uniq = unique(cc)
  Y = matrix(NA, ncol = length(cc.uniq), nrow = nrow(keep))
  
  for(n in 1:ncol(Y))
  {
    jj = which(cc == cc.uniq[n])
    if(length(jj) == 1) {
      Y[,n] = keep[,jj]
    }else{
      Y[,n] = apply(keep[ ,jj], 1, mean)
    }
  }
  colnames(Y) = cc.uniq
  rownames(Y) = rownames(keep)
  
  # X matrix
  motif.oc = readRDS(file = '../results/motif_analysis/motif_oc_fimo.rds')
  mm = match(rownames(motif.oc), rownames(Y))
  motif.oc = motif.oc[!is.na(mm), ]
  
  ss.m = apply(motif.oc, 2, sum)
  motif.oc = motif.oc[ , which(ss.m>0)]
  
  ss.p = apply(motif.oc, 1, sum)
  motif.oc = motif.oc[which(ss.p>0), ]
  
  kk = match(rownames(motif.oc), rownames(Y))
  Y = Y[kk, ]
  
  X = as.matrix(motif.oc)
  Y = as.matrix(Y)
  
  ### specify glment parameters
  alpha = 0
  standardize = TRUE;
  use.lambda.min = FALSE;
  binarize.x = TRUE
  standardize.response=FALSE
  intercept=FALSE
  family = 'mgaussian'
  
  if(binarize.x) x = x > 0
  library(doMC) 
  registerDoMC(cores=6)
  
  x = X;
  y = scale(Y, center = TRUE, scale = FALSE)
  
  #sels = c(1:5000)
  #x = x[sels, ]
  #y = y[sels, ]
  
  library(tictoc)
  tic()
  cv.fit=cv.glmnet(x, y, family= family, grouped=FALSE, 
                   alpha=alpha, nlambda=100, standardize=standardize, 
                   standardize.response=standardize.response, parallel = TRUE)
  
  plot(cv.fit)
  toc()
  
  
  if(use.lambda.min){
    s.optimal = cv.fit$lambda.min
  }else{
    s.optimal = cv.fit$lambda.1se
  }
  
  fit=glmnet(x,y,alpha=alpha, lambda=s.optimal, family=family, 
             standardize=standardize, standardize.response=standardize.response, intercept=intercept, 
             relax = FALSE)
  
  #fit=glmnet(x, y, family='mgaussian', standardize=standardize, standardize.response=standardize.response, intercept=TRUE)
  #plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")
  
  xx = coef(fit, s = s.optimal)
  aa = c()
  for(j in 1:length(xx))
  {
    aa = cbind(aa, as.numeric(xx[[j]]))
  }
  rownames(aa) = rownames(xx[[1]])
  #colnames(aa) = c('E40', 'E44.P', 'mUA', 'BL.UA.D5', 'BL.UA.D9', 'BL.UA.D13.P')
  colnames(aa) = colnames(y)
  aa = as.data.frame(aa[-1, ]) # ignore the intercept
  aa = apply(aa, 2, scale)
  #aa = scale(aa)
  rownames(aa) = rownames(xx[[1]])[-1] 
  
  #kk = apply(aa, 1, function(x) all(abs(x)>10^-6))
  #aa = aa[kk, ]
  
  #ss = apply(aa, 1, function(x) !all(x==0))
  #aa = aa[ss, ]
  #head(rownames(aa)[order(-abs(aa$MSxp))], 10)
  #head(rownames(aa)[order(-abs(aa$MSxa))], 10)
  Test.zscore.cutoff = 2.5
  ss = apply(aa, 1, function(x) length(which(abs(x) > Test.zscore.cutoff)))
  print(aa[which(ss>0), ])
  
  bb = aa[which(ss>0), ]
  bb[which(abs(bb)<Test.zscore.cutoff)] = 0
  pheatmap(bb, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, 
           scale = 'none', cluster_cols=FALSE, main = paste0("motif activity by MARA"), 
           na_col = "white", fontsize_col = 12) 
  
  
  # aa = as.data.frame(coef.glmnet(fit, s = s.optimal))
  # aa = aa[-1, ] # remove intecept
  # colnames(aa) = names(fit$beta)
  # if(alpha > 0.0){
  #   rownames(aa) = rownames(fit$beta[[2]])
  #   res = aa;
  #   ## collect result from the elastic-net
  #   #kk = apply(aa, 1, function(x) !all(x==0))
  #   #aa = aa[kk, ]
  #   #colnames(x)[which(fit$beta[[1]][,optimal]!=0)]
  #   #colnames(x)[which(fit$beta[[2]][,optimal]!=0)]
  #   
  #   # rerun lm with selected features
  #   relax.fitting.lm = FALSE
  #   if(relax.fitting.lm){
  #     fit.lm = lm(y ~ x[, match(rownames(aa), colnames(x))])
  #     res = data.frame(fit.lm$coefficients)
  #     res = res[-1, ] # remove intercept
  #     rownames(res) = rownames(aa)
  #     
  #     pheatmap(res, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
  #              scale = 'column', cluster_cols=FALSE, main = paste0("motif activity"), 
  #              na_col = "white", fontsize_col = 10)
  #   }
  #   
  # }else{
  #   if(zscore.output) aa = apply(aa, 2, scale)
  #   rownames(aa) = rownames(fit$beta[[2]])
  #   #ss = apply(aa, 1, function(x) !all(x==0))
  #   #aa = aa[ss, ]
  #   #head(rownames(aa)[order(-abs(aa$MSxp))], 10)
  #   #head(rownames(aa)[order(-abs(aa$MSxa))], 10)
  #   res = aa
  #   
  #   if(Test){
  #     ss = apply(aa, 1, function(x) length(which(abs(x) > Test.zscore.cutoff)))
  #     aa = aa[which(ss>0), ]
  #     print(aa)
  #   }
  # }
  
}




