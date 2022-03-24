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
  
  # x <- rdf()
  # rdf <- rdf_parse(tfs)
  # 
  # options(rdf_print_format = "turtle")
  # rdf
  # 
  # car_triples <- 
  #   rdf %>% 
  #   rownames_to_column("label") %>% 
  #   gather(attribute,measurement, -Model)
  # 
  
  ##########################################
  # try another code from 
  # https://bioconductor.org/packages/release/data/experiment/vignettes/ELMER.data/inst/doc/vignettes.html 
  ##########################################
  library(xml2)
  library(httr)
  library(dplyr)
  library(rvest)
  
  createMotifRelevantTfs <- function(classification = "family"){
    
    message("Accessing hocomoco to get last version of TFs ", classification)
    file <- paste0(classification,".motif.relevant.TFs.rda")
    
    # Download from http://hocomoco.autosome.ru/human/mono
    tf.family <- "http://hocomoco11.autosome.ru/human/mono?full=true" %>% read_html()  %>%  html_table()
    tf.family <- tf.family[[1]]
    # Split TF for each family, this will help us map for each motif which are the some ones in the family
    # basicaly: for a TF get its family then get all TF in that family
    col <- ifelse(classification == "family", "TF family","TF subfamily")
    family <- split(tf.family,f = tf.family[[col]])
    
    motif.relevant.TFs <- plyr::alply(tf.family,1, function(x){  
      f <- x[[col]]
      if(f == "") return(x$`Transcription factor`) # Case without family, we will get only the same object
      return(unique(family[as.character(f)][[1]]$`Transcription factor`))
    },.progress = "text")
    #names(motif.relevant.TFs) <- tf.family$`Transcription factor`
    names(motif.relevant.TFs) <- tf.family$Model
    # Cleaning object
    attr(motif.relevant.TFs,which="split_type") <- NULL
    attr(motif.relevant.TFs,which="split_labels") <- NULL
    
    return(motif.relevant.TFs)
  }
  
  updateTFClassList <- function(tf.list, classification = "family"){
    col <- ifelse(classification == "family","family.name","subfamily.name")
    TFclass <- getTFClass()
    # Hocomoco
    tf.family <- "http://hocomoco11.autosome.ru/human/mono?full=true" %>% read_html()  %>%  html_table()
    tf.family <- tf.family[[1]]
    
    tf.members <- plyr::alply(unique(TFclass %>% pull(col)),1, function(x){  
      TFclass$Gene[which(x == TFclass[,col])]
    },.progress = "text")
    names(tf.members) <- unique(TFclass %>% pull(col))
    attr(tf.members,which="split_type") <- NULL
    attr(tf.members,which="split_labels") <- NULL
    
    for(i in names(tf.list)){
      x <- tf.family[tf.family$Model == i,"Transcription factor"]
      idx <- which(sapply(lapply(tf.members, function(ch) grep(paste0("^",x,"$"), ch)), function(x) length(x) > 0))
      if(length(idx) == 0) next
      members <- tf.members[[idx]]
      tf.list[[i]] <- sort(unique(c(tf.list[[i]],members)))
    }
    return(tf.list)
  }
  
  getTFClass <- function(){
    # get TF classification
    file <- "TFClass.rda"
    if(file.exists(file)) {
      return(get(load(file)))
    }
    file <- "http://tfclass.bioinf.med.uni-goettingen.de/suppl/tfclass.ttl.gz"
    downloader::download(file,basename(file))
    char_vector <- readLines(basename(file))
    # Find TF idx
    idx <- grep("genus",char_vector,ignore.case = T)
    
    # get TF names
    TF <- char_vector[sort(c( idx +1, idx + 2, idx + 4))]
    TF <- TF[-grep("LOGO_|rdf:type",TF)]
    TF <- gsub("  rdfs:label | ;| rdfs:subClassOf <http://sybig.de/tfclass#|>","",TF)
    TF <- stringr::str_trim(gsub('"', '', TF))
    TF <- tibble::as.tibble(t(matrix(TF,nrow = 2)))
    colnames(TF) <- c("Gene", "class")
    
    # Get family and subfamily classification
    family.pattern <-  "^<http://sybig.de/tfclass#[0-9]+\\.[0-9]+\\.[0-9]+>"
    
    idx <- grep(family.pattern,char_vector)
    family.names <- char_vector[ sort(c(idx,idx+ 2))]
    family.names <- gsub("  rdfs:label | ;| rdfs:subClassOf <http://sybig.de/tfclass#|>|<http://sybig.de/tfclass#| rdf:type owl:Class","",family.names)
    family.names <- stringr::str_trim(gsub('"', '', family.names))
    family.names <- tibble::as.tibble(t(matrix(family.names,nrow = 2)))
    colnames(family.names) <- c("family", "family.name")
    
    
    subfamily.pattern <-  "^<http://sybig.de/tfclass#[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+>"
    
    idx <- grep(subfamily.pattern,char_vector)
    subfamily.names <- char_vector[ sort(c(idx,idx+ 2))]
    subfamily.names <- gsub("  rdfs:label | ;| rdfs:subClassOf <http://sybig.de/tfclass#|>|<http://sybig.de/tfclass#| rdf:type owl:Class","",subfamily.names)
    subfamily.names <- stringr::str_trim(gsub('"', '', subfamily.names))
    subfamily.names <- tibble::as.tibble(t(matrix(subfamily.names,nrow = 2)))
    colnames(subfamily.names) <- c("subfamily", "subfamily.name")
    subfamily.names$family <- stringr::str_sub(subfamily.names$subfamily,1,5)
    
    classification <- left_join(family.names,subfamily.names)
    classification$class <- ifelse(is.na(classification$subfamily),classification$family,classification$subfamily)
    
    # Add classification to TF list
    TFclass <- left_join(TF,classification)
    
    # Break ( into multiple cases)
    m <- grep("\\(|/",TFclass$Gene)
    df <- NULL
    for(i in m){
      gene <- sort(stringr::str_trim(unlist(stringr::str_split(TFclass$Gene[i],"\\(|,|\\)|/"))))
      gene <- gene[stringr::str_length(gene) > 0]
      aux <- TFclass[rep(i,length(gene)),]
      aux$Gene <- gene
      df <- rbind(df,aux)
    }
    TFclass <- rbind(TFclass,df)
    TFclass <- TFclass[!duplicated(TFclass),]
    
    # Break ( into multiple cases)
    m <- grep("-",TFclass$Gene)
    df <- NULL
    for(i in m){
      gene <- gsub("-","",sort(stringr::str_trim(unlist(stringr::str_split(TFclass$Gene[i],"\\(|,|\\)|/")))))
      gene <- gene[stringr::str_length(gene) > 0]
      aux <- TFclass[rep(i,length(gene)),]
      aux$Gene <- gene
      df <- rbind(df,aux)
    }
    TFclass <- rbind(TFclass,df)
    
    library(limma)
    require(org.Hs.eg.db)
    df <- NULL
    for(i in 1:length(TFclass$Gene)){
      m <- TFclass$Gene[i]
      gene <- unique(c(toupper(alias2Symbol(toupper(m))),toupper(m),toupper(alias2Symbol(m))))
      if(all(gene %in% TFclass$Gene)) next
      aux <- TFclass[rep(i,length(gene)),]
      aux$Gene <- gene
      df <- rbind(df,aux)
    }
    TFclass <- rbind(TFclass,df)
    TFclass <- TFclass[!duplicated(TFclass),]
    TFclass <- TFclass[TFclass$Gene %in% human.TF$external_gene_name,]
    save(TFclass,file = "TFClass.rda")
    return(TFclass)
  }
  TF.family <- createMotifRelevantTfs("family")
  TF.family <- updateTFClassList(TF.family,"family")
  TF.subfamily <- createMotifRelevantTfs("subfamily")
  TF.subfamily <- updateTFClassList(TF.subfamily,classification = "subfamily")
  save(TF.family,file = "~/ELMER.data/data/TF.family.rda", compress = "xz")
  save(TF.subfamily,file = "~/ELMER.data/data/TF.subfamily.rda", compress = "xz")
  
  ##########################################
  # human TF downalod and used from 
  # https://bioconductor.org/packages/release/data/experiment/vignettes/ELMER.data/inst/doc/vignettes.html
  # reference : A curated list of TF was retrieved from Lambert, Samuel A., et al. 
  # “The human transcription factors.” Cell 172.4 (2018): 650-665 (Lambert, Samuel A and Jolma, Arttu and 
  #Campitelli, Laura F and Das, Pratyush K and Yin, Yimeng and Albu, Mihai and Chen, #
  # Xiaoting and Taipale, Jussi and Hughes, Timothy R and Weirauch, Matthew T 2018) with the following code.
  ##########################################
  human.TF <- readr::read_csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")
  human.TF <- human.TF[human.TF$`Is TF?` == "Yes",]
  colnames(human.TF)[1:2] <- c("ensembl_gene_id","external_gene_name")
  human.TF = as.data.frame(human.TF)
  
  saveRDS(human.TF, 
          file = paste0('/Users/jiwang/workspace/imp/positional_memory/results/motif_analysis/TFs_annot/curated_human_TFs_Lambert.rds'))
  
  
}

collect.TFs.gene.expression = function()
{
  tfs = readRDS(file = paste0('../results/motif_analysis/TFs_annot/curated_human_TFs_Lambert.rds'))
  load(file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))
  annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr.rds'))
  
  #kk = match(annot$gene.symbol.hs, tfs)
  
  
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


# manually add extra motifs for hnd-1, pha-4, unc-120 and nhr-67 from dm, mus and homo
convert.cisbp.format.to.meme = function()
{
  library(universalmotif)
  #pwmDIr = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/extra_pwm'
  pwm.cisbp = '../data/test/PWM.txt'
  xx = read_cisbp(file = pwm.cisbp, skip = 0)
  yy = convert_type(xx, "PWM")
  
  write_meme(yy, file = '../data/test/PWM_converted.meme', overwrite = TRUE)
  
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

## modify the motif names with associated TFs
process.motif.tf.mapping = function()
{
  motif.tf = read.table('/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/motifs_tfs_mapping_curated_extra.txt', 
                        header = FALSE, sep = ' ')
  motif.tf = motif.tf[, c(2:3)]
  colnames(motif.tf) = c('motifs', 'tfs')
  
  # manually modify the motif names
  motif.tf = data.frame(motif.tf, stringsAsFactors = FALSE)
  motif.tf$motifs.new = motif.tf$motifs
  motif.tf$tfs.new = motif.tf$tfs
  
  xx = motif.tf
  #xx$motifs.new = gsub('1.02', '', xx$motifs.new)
  #xx$motifs.new = gsub('1.02', '', xx$motifs.new)
  #xx$tfs.new = paste0(xx$tfs.new, '_', xx$tfs)
  #xx$tfs.new = gsub('-', '', xx$tfs.new)
  xx$tfs.new = gsub(':', '.', xx$tfs.new)
  xx$tfs.new = gsub('/', '.', xx$tfs.new)
  xx$tfs.new = gsub("\\(","", xx$tfs.new)
  xx$tfs.new = gsub("\\)","", xx$tfs.new)
  xx$tfs.new = gsub("_Homo_sapiens_DBD*.*",".homo", xx$tfs.new)
  xx$tfs.new = gsub("_Caenorhabditis_briggsae_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Drosophila_melanogaster_DBD*.*",".dm", xx$tfs.new)
  xx$tfs.new = gsub("_Mus_musculus_DBD*.*",".mus", xx$tfs.new)
  xx$tfs.new = gsub("_Brugia_pahangi_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Wuchereria_bancrofti_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_PBM_CONSTRUCTS_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Tetraodon_nigroviridis_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub('Mus', '.mus', xx$tfs.new)
  xx$tfs.new = gsub('Dm', '.dm', xx$tfs.new)
  xx$tfs.new = gsub('Homo', '.homo', xx$tfs.new)
  xx$tfs.new = gsub('Hand1', 'hnd-1', xx$tfs.new)
  xx$tfs.new = gsub('Foxa1', 'pha-4', xx$tfs.new)
  xx$tfs.new = gsub('SRF', 'unc-120', xx$tfs.new)
  xx$tfs.new = gsub('Srf', 'unc-120', xx$tfs.new)
  xx$tfs.new = gsub('Nr2e1', 'nhr-67', xx$tfs.new)
  xx$tfs.new = gsub('NR2E1', 'nhr-67', xx$tfs.new)
  
  xx$motifs.new = paste0(xx$motifs.new, '_', xx$tfs.new)
  
  motif.tf = xx
  
  saveRDS(motif.tf, file = '../data/motifs_tfs/motif_tf_mapping.rds')
  
}

# remove motif redundancies
remove.motifs.redundancy.by.similarity.clustering = function()
{
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering visualization
  require(graphics)
  library(universalmotif)
  
  motif.tf = readRDS(file = '../data/motifs_tfs/motif_tf_mapping.rds')
  motif.tf$motifs.new = paste0(motif.tf$tfs.new, '_', motif.tf$motifs)
  
  pwm = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015_curated_extra.meme'
  #pwm = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.meme'
  
  meme= read_meme(file = pwm, skip = 0, readsites = FALSE, readsites.meta = FALSE)
  motifs = convert_motifs(meme, class = "universalmotif-universalmotif")
  cat(length(motifs), ' motifs \n')
  
  pwm.corr <- compare_motifs(motifs, method = "PCC", min.mean.ic = 0,
                             score.strat = "a.mean")
  
  newName = motif.tf$motifs.new[match(rownames(pwm.corr), motif.tf$motifs)]
  rownames(pwm.corr) = newName
  colnames(pwm.corr) = newName
  
  comparisons <- 1 - pwm.corr
  dd <- as.dist(comparisons)
  
  # Hierarchical clustering using Complete Linkage
  hc <- hclust(dd, method = "ward.D2" )
  
  # Plot the obtained dendrogram
  #plot(hc, cex = 0.6, hang = -1)
  #sub_grp <- cutree(hc, h = 0.1)
  pdfname = paste0(resDir, "/pwm_celegans_similarity_clustering.pdf")
  pdf(pdfname, width=20, height = 30)
  par(cex =0.5, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #plot(hc, cex = 0.5, hang = -1)
  plot(as.dendrogram(hc), cex=0.5, horiz=TRUE)
  abline(v = c(0.05, 0.1, 0.15), col = c('blue', 'red', 'green'))
  #rect.hclust(hc, h = hc.cutoff, border="darkred")
  #groups <- 
  length(unique(cutree(hc, h = 0.05)))
  length(unique(cutree(hc, h = 0.1)))
  length(unique(cutree(hc, h = 0.15)))
  length(unique(cutree(hc, h = 0.2)))
  length(unique(cutree(hc, h = 0.25)))
  
  dev.off()
  
  change.pwm.logo.names = FALSE
  if(change.pwm.logo.names){
    logoDir = '../data/motifs_tfs/pwm_logos'
    logo.file = list.files(path = logoDir, pattern = '*.pdf', full.names = TRUE)
    for(n in 1:nrow(motif.tf))
    {
      # n = 1 
      cmd = paste0('mv ', logo.file[grep(motif.tf$motifs[n], logo.file)], ' ', logoDir, '/',  motif.tf$motifs.new[n], '.pdf')
      system(cmd)
    }
    
  }
  #fviz_nbclust(diss = comparisons, FUN = hcut, method = "wss")
  #fviz_nbclust(df, FUN = hcut, method = "silhouette")
  
  ##########################################
  # merge motifs using height = 0.1 and change motif names
  ##########################################
  hc.cutoff = 0.1
  
  groups <- cutree(hc, h = hc.cutoff)
  motif.tf = data.frame(motif.tf, group = groups, stringsAsFactors = FALSE)
  motif.tf$names = NA
  for(nn in unique(motif.tf$group))
  {
    # nn = 5
    kk = which(motif.tf$group == nn)
    motif.tf$names[kk] = paste0(paste0(unique(motif.tf$tfs.new[kk]), collapse = '_'), '.M', nn)
    
  }
  
  # save motif-to-tf mapping with redundancy removal information
  saveRDS(motif.tf, file = '../data/motifs_tfs/motif_tf_mapping.rds') 
  
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
# Section II : all processing steps 
# after motif scanning 
# before motif activity analysis
########################################################
########################################################
save.peak.bed.file.for.fimo = function()
{
  #load(file = paste0(RdataDir, '/samplesDesign.cleaned_readCounts.withinPeaks.pval6.Rdata'))
  fpm = readRDS(file = paste0(RdataDir, '/fpm_TMM_combat.rds'))
  
  # prepare the background distribution
  jj = grep('bg_', rownames(fpm), invert = TRUE)
  fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
  fpm = fpm[jj, ]
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  
  pp$name = rownames(fpm)
  pp$strand = '*'
  
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  ll = width(pp)
  
  
  write.table(pp, file = '../results/motif_analysis/peaks/peaks_for_fimo.bed', row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = '\t')
  
  
}

##########################################
# process fimo output .tsv file to bed file 
##########################################
make.bed.file.from.fimo.out = function()
{
  # manually define big domain for HOXA cluster chr2p:870,541,422-875,470,873 
  HoxA = data.frame(chr = 'chr2p', start = 870541422, end = 875470873 , strand = '*', stringsAsFactors = FALSE)
  HoxA = makeGRangesFromDataFrame(HoxA, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'strand')
  
  fimo = read.delim(file = paste0('../results/motif_analysis/FIMO/fimo_out/fimo_chr2p.tsv'), header = FALSE)
  fimo = fimo[, -c(2, 10)]
  fimo = data.frame(fimo, stringsAsFactors = FALSE)
  
  x = data.frame(t(sapply(fimo$V3, function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  
  fimo =data.frame(x, fimo)
  colnames(fimo)[c(1:3)] = c('chr', 'start', 'end')
  
  
  x$strand = '*'
  x = makeGRangesFromDataFrame(x, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  head(fimo)
  
  jj = overlapsAny(x, HoxA)
  
  fimo = fimo[jj, ]
  colnames(fimo)[c(1:3)] = c('chr', 'start', 'end')
  colnames(fimo)[4:11] = c('motif_id', 'sequence_name', 'start_sequence',  'end_sequence', 'strand_sequence', 'score', 'pval', 'qv')
  
  x = fimo[which(fimo$pval < 10^-6), ]
  ss = table(x$motif_id)
  ss[ss>0]
  
  x$start = as.numeric(as.character(x$start))
  x$end = as.numeric(as.character(x$end))
  x$start_sequence = as.numeric(as.character(x$start_sequence))
  x$end_sequence = as.numeric(as.character(x$end_sequence))
  
  kk = which(x$strand_sequence == '+')
  x$start_sequence[kk] = x$start[kk]+ x$start_sequence[kk]
  x$end_sequence[kk] = x$start[kk] + x$end_sequence[kk]
  
  kk = which(x$strand_sequence == '-')
  lls = x$end_sequence[kk] - x$start_sequence[kk]
  x$start_sequence[kk] = x$end[kk] -  x$end_sequence[kk]
  x$end_sequence[kk] = x$start_sequence[kk] + lls
 
  y = x[, c(1, 6, 7, 4, 10, 8)]
    
  write.table(y, file = paste0(resDir, '/fimo_out_HoxACluster.bed'), sep = '\t', 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
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
  
  saveRDS(motif.oc, file = '../results/motif_analysis/motif_oc_fimo_v2.rds')
  
}

########################################################
########################################################
# Section III: run MARA 
# 
########################################################
########################################################
run.RF.otherMethods = function()
{
  library(randomForest)
  tic()
  rf <- randomForest::randomForest(x = x, y = y[, 2], ntree = 200, keep.forest = FALSE, 
                                         importance = TRUE)
  toc()
  
  varImpPlot(rf, n.var = 50)
  
  RF.parameter.tuing = FALSE
  if(RF.parameter.tuing)
  {
    library(randomForest)
    library(ranger)
    library(stats)
    library(mlbench)
    library(caret)
    library(tictoc)
    
    ## tune mtry parameter with tuneRF from randomForest package
    m2 <- tuneRF(
      x          = train,
      y          = y,
      ntreeTry   = ntree,
      mtryStart  = floor(sqrt(ncol(x))),
      stepFactor = 1.5,
      improve    = 0.01,
      trace      = FALSE      # to not show real-time progress 
    )
    print(m2)
    
    # hyperparameter grid search
    hyper_grid <- expand.grid(
      mtry       = seq(20, 60, by = 5),
      node_size  = seq(1, 10, by = 2),
      #sample_size = c(.55, .632, .70, .80),
      sample_size = c(.632),
      prediction_err   = 0
    )
    
    # total number of combinations
    nrow(hyper_grid)
    ## [1] 96
    
    for(i in 1:nrow(hyper_grid)) {
      # train model
      model <- ranger(
        x = train, 
        y = y,
        num.trees       = ntree,
        mtry            = hyper_grid$mtry[i],
        min.node.size   = hyper_grid$node_size[i],
        sample.fraction = hyper_grid$sample_size[i],
        seed            = 123
      )
      
      # add OOB error to grid
      hyper_grid$prediction_err[i] <- model$prediction.error
    }
    
    hyper_grid %>% 
      dplyr::arrange(prediction_err) %>%
      head(10)
    
    index.optimal = which.min(hyper_grid$prediction_err)
    tic()
    train.optimal = ranger::ranger(x = train, y = y, 
                                   num.trees = ntree, 
                                   mtry = hyper_grid$mtry[index.optimal], 
                                   min.node.size = hyper_grid$node_size[index.optimal],
                                   sample.fraction = hyper_grid$sample_size[index.optimal], 
                                   write.forest = TRUE, 
                                   probability = TRUE,
                                   classification = TRUE, 
                                   local.importance = TRUE)
    toc()
    
    saveRDS(train.optimal, file = paste0(RdataDir, 'MurrayData_classifier_test_RF_optimalParam.rds'))
    rf.fit = readRDS(file = paste0(RdataDir, 'MurrayData_classifier_test_RF_optimalParam.rds'))
    #pred_test <-stats::predict(rf.fit, test)
    
    err.test = mean(pred_test==factor(y.test, levels = levels(pred_test)))
  }
  
  #Prediction <- stats::predict(train_rf, test, type = "prob")
  #Prediction <- stats::predict(m1, test, type = "prob")
  #Prediction <- predict(train_rf, test, type = "response")$predictions
  #rf.res = data.frame(label = apply(Prediction, 1, function(x) colnames(Prediction)[which.max(x)]),
  #                    prob = apply(Prediction, 1, function(x) x[which.max(x)]), stringsAsFactors = FALSE)
  
}

########################################################
########################################################
# Section : MARA for position-dependent peaks 
# 
########################################################
########################################################
run.MARA.atac.spatial = function(keep, cc)
{
  require(glmnet)
  library(pheatmap)
  library(RColorBrewer)
  library(scchicFuncs)
  
  # prepare Y response matrix
  Prepare.Response.Matrix = FALSE
  if(Prepare.Response.Matrix){
    res = readRDS(file = paste0(RdataDir, '/res_position_dependant_test_', version.analysis, '_v10.rds'))
    
    # select the positional peaks with pairwise comparisions 
    # limma logFC is in log2 scale
    fdr.cutoff = 0.01; logfc.cutoff = 1
    jj = which((res$adj.P.Val.mLA.vs.mUA < fdr.cutoff & abs(res$logFC.mLA.vs.mUA) > logfc.cutoff) |
                 (res$adj.P.Val.mHand.vs.mUA < fdr.cutoff & abs(res$logFC.mHand.vs.mUA) > logfc.cutoff)|
                 (res$adj.P.Val.mHand.vs.mLA < fdr.cutoff & abs(res$logFC.mHand.vs.mLA) > logfc.cutoff)
    )
    cat(length(jj), '\n')
    
    jj1 = which(res$prob.M0<0.01 & res$log2FC>logfc.cutoff)
    cat(length(jj1), '\n')
    jj2 = which(res$pval.lrt < 0.001 & res$log2FC > logfc.cutoff)
    cat(length(jj2), '\n')
    
    xx = res[c(jj), ]
    #xx = xx[order(-xx$log2FC.mature), ]
    xx[grep('HOXA13|SHOX', xx$transcriptId), c(1:8)]
    fpm[which(rownames(fpm) == 'chr2p:873464923-873465440'), ]
    
    xx[grep('HOXA13|SHOX|MEIS', xx$transcriptId), ]
    
    cat(nrow(xx), ' peaks left\n')
    # sort positional peaks with logFC
    #xx = xx[order(-xx$logFC.mean), ]
    xx = xx[order(-xx$log2FC), ]
    
    ########
    ## asscociate the signifiant postional peaks with expression matrix
    ########
    conds = c("Mature_UA", "Mature_LA", "Mature_Hand", 'HEAD')
    
    sample.sels = c();  cc = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n]) 
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    keep = fpm[(match(rownames(xx), rownames(fpm))), sample.sels]
    keep = as.matrix(keep)
    
    ### filter the peaks from head control sample
    ## either filter soly based on the head signals or based on the logFC between max(mature samples/control) or both
    Filtering.peaks.in.Head.samples = FALSE
    if(Filtering.peaks.in.Head.samples){
      maxs = apply(keep[, grep('Mature_', colnames(keep))], 1, function(x) return(max(c(mean(x[1:5]), mean(x[6:9]), mean(x[10:12])))))
      
      ctl.mean = apply(keep[, grep('HEAD', colnames(keep))], 1, mean)
      
      rr = maxs - ctl.mean
      #p.ctl = pp[match(names(ctl.sels), names(pp))]
      #non.overlap = !overlapsAny(p0, p.ctl)
      plot(maxs, ctl.mean, cex = 0.2);
      abline(0, 1, lwd = 2.0, col = 'red')
      abline(v = 3, lwd = 2.0, col = 'red')
      abline(h = 3, lwd = 2.0, col = 'red')
      
      #xx = xx[non.overlap, ]
      plot(rr, ctl.mean, cex = 0.2);
      abline(h = 3, col = 'red', lwd = 2.0)
      abline(v = c(0.5, 1), col = 'red', lwd = 2.0)
      
      sels = which(rr>1 & ctl.mean<3 & maxs > 3)
      cat(length(sels), 'peaks selected \n')
      
      xx = xx[sels, ]
      keep = keep[sels, ]
      
    }
    
    #load(file = paste0(RdataDir, '/ATACseq_positionalPeaks_excluding.headControl', version.analysis, '.Rdata'))
    #conds = c("Mature_UA", "Mature_LA", "Mature_Hand", 'HEAD')
    
    sample.sels = c();  cc = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n]) 
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    
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
    
  }
  
  # prepare X matrix from motif occurrency matrix
  motif.oc = readRDS(file = '../results/motif_analysis/motif_oc_fimo_v2.rds')
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
  
  Run.Bayesian.ridge = FALSE
  if(Run.Bayesian.ridge){ 
    # the original code from Jake 
    # https://github.com/jakeyeung/scchic-functions/blob/master/scripts/motevo_scripts/lib/run_ridge_regression2.R
    
    E = as.matrix(Y) # exp: matrix of expression, row centered.
    N = as.matrix(X)
    
    E = t(apply(E, 1, scale, center = TRUE, scale = FALSE))
    colnames(E) = colnames(Y)
    
    opt =  scchicFuncs::optimize.lambda(N, E)
    
    r = ridge.regression(N, E, opt$lambda.opt)
    
    saveRDS(r, file = paste0(resDir, '/MARA_Bayesian_ridge.rds'))
    
    zz = r$Zscore
    zz = apply(as.matrix(zz[, c(1:3)]), 1, function(x){x.abs = abs(x); return(x[which(x.abs == max(x.abs))][1]); })
    r$max.Zscore = abs(zz)
    
    # = sort(r$combined.Zscore, decreasing=TRUE)[1:50]
    sort(r$combined.Zscore, decreasing=TRUE)[1:20]
    sort(r$max.Zscore, decreasing=TRUE)[1:20]
    
    sort(r$max.Zscore, decreasing=TRUE)[1:50]
    
    r$max.Zscore[grep('MEIS|RAR|RXR|SOX9', names(r$max.Zscore))]
    
    top20 = sort(r$max.Zscore, decreasing = TRUE)[1:20]
    motif.names = rownames(r$Zscore)
    bb = r$Zscore[match(names(top20), motif.names), ]
    
    df <- data.frame(colnames(bb))
    rownames(df) = colnames(bb)
    colnames(df) = 'segments'
    annot_colors = c('springgreen4', 'steelblue2', 'gold2', 'darkgray')
    names(annot_colors) = c('Mature_UA', 'Mature_LA', 'Mature_Hand', 'HEAD')
    annot_colors = list(segments = annot_colors)
    
    pheatmap(bb, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, 
             scale = 'none', cluster_cols=FALSE, main = paste0(" Motif activity (Z-scores) by MARA"), 
             na_col = "white", fontsize_row = 12, annotation_col = df, 
             annotation_colors = annot_colors,
             filename = paste0(resDir, '/MARA_bayesian_ridge_positional_peaks.pdf'), 
             width = 8, height = 6)
    
  }
  
  Run.RF = FALSE
  if(Run.RF){
    library(randomForest)
    library(tictoc)
    library(ggrepel)
    
    E = as.matrix(Y) # exp: matrix of expression, row centered.
    N = as.matrix(X)
    
    E = t(apply(E, 1, scale, center = TRUE, scale = FALSE))
    colnames(E) = colnames(Y)
    
    tic()
    
    aa = matrix(NA, ncol = ncol(E), nrow = ncol(N))
    rownames(aa) = colnames(N)
    bb = aa
    
    for(n in 1:ncol(E)){
      cat(n, '\n')
      rf <- randomForest::randomForest(x = N, y = E[, n], ntree = 200, keep.forest = FALSE, 
                                       importance = TRUE)
      #aa[, (2*n-1)] = rf$importance[ ,1]
      #aa[, (2*n)] = rf$importance[ ,2 ]
      aa[, n] = importance(rf, type = 1)
      bb[, n] = importance(rf, type = 2)
    }
    
    toc()
    
    #rownames(aa) = rownames(rf$importance)
    
    varImpPlot(rf, n.var = 20)
    imps = importance(rf, type = 1)
    
    zz1 = apply(aa[,c(1:3)], 1, max)
    aa = aa[order(-zz1), ]
    zz1 = zz1[order(-zz1)]
    
    zz2 = apply(bb[, c(1:3)], 1, max)
    bb = bb[order(-zz2), ]
    zz2 = zz2[order(-zz2)]
    
    tops = 20
    plot(zz2[1:tops], c(tops:1))
    
    
    xx = data.frame(names = rownames(aa)[1:tops], scores = zz1[1:tops], rank = c(tops:1))
    yy = data.frame(names = rownames(bb)[1:tops], scores = zz2[1:tops], rank = c(tops:1))
    
    gp = ggplot(data = xx, aes(x = scores, y = rank, label = names)) +   
      geom_point(size = 3.0, color = 'blue') +
      theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) + 
      #geom_text_repel(data=subset(yy, pvalue_pos > 2), size = 4)
      geom_label_repel(size = 4)
    
    
    gp = ggplot(data = yy, aes(x = scores, y = rank, label = names)) +   
      geom_point(size = 2.0, color = 'blue') +
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12)) + 
      geom_text_repel(size = 4) + ggtitle('Importance score from Random Forest MARA')
      #geom_label_repel(size = 4)
      #theme(axis.text.x = element_text(angle = 90, size = 8)) + 
      #geom_hline(yintercept = c(20, 50, 100)) + ylab("unique.rmdup (M)")
    
    plot(gp) + ggsave(paste0(resDir, "/positional_peaks_RF_importance.pdf"), width=12, height = 8)
    
    pheatmap(xx, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, 
             scale = 'none', cluster_cols=FALSE, main = paste0("RF importance scores "), 
             na_col = "white", fontsize_row = 12, 
             filename = paste0(resDir, '/positional_peaks_RF_importance.pdf'), 
             width = 8, height = 12) 
    
    
  }
  
  Run.glmnet = FALSE
  if(Run.glmnet){
    ### specify glment parameters
    alpha = 1
    standardize = TRUE;
    use.lambda.min  = FALSE
    binarize.x = FALSE
    standardize.response=FALSE
    intercept=TRUE
    family = 'mgaussian'
    
   
    library(doMC) 
    registerDoMC(cores=6)
    
    #x = t(scale(t(X), center = TRUE, scale = TRUE));
    x = X
    y = scale(Y[, c(1:4)], center = TRUE, scale = FALSE)
    if(binarize.x) x = x > 0
    #sels = c(1:5000)
    #x = x[sels, ]
    #y = y[sels, ]
    
    library(tictoc)
    tic()
    cv.fit=cv.glmnet(x, y, family= family, grouped=FALSE, 
                     alpha=alpha, nlambda=200, standardize=standardize, 
                     standardize.response=standardize.response, parallel = TRUE)
    
    plot(cv.fit)
    toc()
    
    if(use.lambda.min){
      s.optimal = cv.fit$lambda.min
    }else{
      s.optimal = cv.fit$lambda.1se
    }
    
    fit=glmnet(x,y,alpha=alpha, lambda=cv.fit$lambda, family=family, 
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
    
    aa = aa[apply(aa, 1, function(x) all(x !=0)), ]
    
    
    #aa = apply(aa, 2, scale)
    aa = scale(aa)
    rownames(aa) = rownames(xx[[1]])[-1] 
    
    #kk = apply(aa, 1, function(x) all(abs(x)>10^-6))
    #aa = aa[kk, ]
    
    #ss = apply(aa, 1, function(x) !all(x==0))
    #aa = aa[ss, ]
    #head(rownames(aa)[order(-abs(aa$MSxp))], 10)
    #head(rownames(aa)[order(-abs(aa$MSxa))], 10)
    Test.zscore.cutoff = 2.0
    ss = apply(aa, 1, function(x) length(which(abs(x) > Test.zscore.cutoff)))
    print(aa[which(ss>0), ])
    print(aa[grep('MEIS|RXR|HOXA13|HOXD13', rownames(aa)), ])
    
    bb = aa[which(ss>0), ]
    bb[which(abs(bb)<Test.zscore.cutoff)] = 0
    pheatmap(bb, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, 
             scale = 'none', cluster_cols=FALSE, main = paste0("motif activity by MARA"), 
             na_col = "white", fontsize_col = 12) 
    
    
    
    Test.ridge.package = FALSE
    if(Test.ridge.package){
      library(ridge)
      fit1 = linearRidge(y[, 1] ~ x)
      aa1 = pvals(fit1)
      aa1 = data.frame(aa1$pval)
      
      fit2 = linearRidge(y[, 2] ~ x)
      aa2 = pvals(fit2)
      aa2 = data.frame(aa2$pval)
      
      fit3 = linearRidge(y[, 3] ~ x)
      aa3 = pvals(fit3)
      aa3 = data.frame(aa3$pval)
      
      index = 20
      pvs = data.frame(aa1[, index], aa2[, index], aa3[, index])
      rownames(pvs) = gsub('x', '', rownames(aa1))
      
      print(rownames(pv)[which(pv[, index]<0.01)])
      
      print(pvs[grep('SOX9|MEIS|RXR|RAR|HOXA13', rownames(pvs)), ])
      
      pval.cutoff = 0.01
      ss = apply(pvs, 1, function(x) length(which(x < pval.cutoff)))
      bb = pvs[which(ss>=1), ]
      colnames(bb) = colnames(y)
      bb = -log10(bb)
      
      for(n in 1:ncol(bb)){bb[which(bb[,n] < (-log10(pval.cutoff))) ,n] = NA}
      #bb[which( bb< pval.cutoff)] = NA
      
      pheatmap(bb, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, 
               scale = 'none', cluster_cols=FALSE, main = paste0("motif significance (-log10 pval) by ridge regression"), 
               na_col = "gray", fontsize_col = 12)
    }
    
  }
  
}


run.MARA.atac.temporal = function(keep, cc)
{
  require(glmnet)
  library(pheatmap)
  library(RColorBrewer)
  library(scchicFuncs)
  
  # prepare Y response matrix
  Prepare.Response.Matrix = FALSE
  if(Prepare.Response.Matrix){
    fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_mUA_regeneration_embryoStages.rds'))
    design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_mUA_regeneration_embryoStages.rds'))
    
    # prepare the background distribution
    fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
    fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
    rownames(fpm) = gsub('_', '-', rownames(fpm))
    
    hist(fpm.bg, breaks = 100, main = 'background distribution')
    abline(v = 1, col = 'red', lwd = 2.0)
    quantile(fpm.bg, c(0.95, 0.99))
    
    res = readRDS(file = paste0(RdataDir, '/res_temporal_dynamicPeaks_test_v5.rds'))
    
    # select the temporal dynamic peaks
    length(which(res$prob.M0<0.05))
    length(which(res$prob.M0<0.05 & res$log2FC > 1))
    length(which(res$prob.M0<0.01 & res$log2FC > 1))
    length(which(res$prob.M0<0.01 & res$log2FC > 1.5))
    length(which(res$prob.M0<0.01 & res$log2FC > 2))
    
    jj = which(res$prob.M0 < 0.01 & res$log2FC > 1 )
    
    conds = c("Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal", "BL_UA_13days_distal")
    
    sample.sels = c(); cc = c()
    for(n in 1:length(conds)) {
      kk = which(design$conds == conds[n])
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
    }
    
    xx = res[c(jj), ]
    xx = xx[order(-xx$log2FC), ]
    
    keep = fpm[!is.na(match(rownames(fpm), rownames(xx))), sample.sels]
  }
  
  
  keep = as.matrix(keep)
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
  
  # prepare X matrix from motif occurrency matrix
  motif.oc = readRDS(file = '../results/motif_analysis/motif_oc_fimo_v2.rds')
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
  
  Run.Bayesian.ridge = FALSE
  if(Run.Bayesian.ridge){ 
    # the original code from Jake 
    # https://github.com/jakeyeung/scchic-functions/blob/master/scripts/motevo_scripts/lib/run_ridge_regression2.R
    
    E = as.matrix(Y) # exp: matrix of expression, row centered.
    N = as.matrix(X)
    
    E = t(apply(E, 1, scale, center = TRUE, scale = FALSE))
    colnames(E) = colnames(Y)
    
    opt =  scchicFuncs::optimize.lambda(N, E)
    
    r = ridge.regression(N, E, opt$lambda.opt)
    
    # saveRDS(r, file = paste0(RdataDir, '/MARA_Bayesian_ridge_regenerationPeaks.rds'))
    
    zz = r$Zscore
    zz = apply(as.matrix(zz[, c(1:3)]), 1, function(x){x.abs = abs(x); return(x[which(x.abs == max(x.abs))][1]); })
    r$max.Zscore = abs(zz)
    
    #zz = r$Zscore
    #zz = apply(as.matrix(zz), 1, max)
    #r$max.Zscore = zz
    
    # = sort(r$combined.Zscore, decreasing=TRUE)[1:50]
    sort(r$combined.Zscore, decreasing=TRUE)[1:20]
    
    sort(r$max.Zscore, decreasing=TRUE)[1:20]
    sort(r$max.Zscore, decreasing=TRUE)[1:30]
    sort(r$max.Zscore, decreasing=TRUE)[1:40]
    sort(r$max.Zscore, decreasing=TRUE)[1:50]
    
    topMotifs = sort(r$max.Zscore, decreasing = TRUE)[1:50]
    motif.names = rownames(r$Zscore)
    bb = r$Zscore[match(names(topMotifs), motif.names), ]
    
    df <- data.frame(colnames(bb))
    rownames(df) = colnames(bb)
    colnames(df) = 'regeneration time'
    
    pheatmap(bb, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE, 
             colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10),
             scale = 'row', cluster_cols=FALSE, main = '', 
             na_col = "white", fontsize_row = 12, 
             annotation_col = df, 
             filename = paste0(resDir, '/MARA_bayesianRidge_temporalpeaks.pdf'), 
             width = 10, height = 12) 
    
    ##########################################
    # compare with the smartseq2 data 
    ##########################################
    source('Functions_atac.R')
    rnaDir = "../results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/"
    smartseq2 = readRDS(file = paste0(rnaDir, 'smartseq2_R10724_R11635_cpm.batchCorrect_DESeq2.test.withbatch.log2FC.shrinked.rds'))
    cpm = smartseq2[, c(1,2, 5:12)]
    #res = res[, -c(1:12)]
        
    #conds = c("Mature_UA", "BL_UA_5days", "BL_UA_9days", "BL_UA_13days_proximal",  "BL_UA_13days_distal")
    sample.sels = c();  
    cc = c()
    sample.means = c()
    for(n in 1:length(conds)) 
    {
      kk = grep(conds[n], colnames(cpm))
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
      if(length(kk)>1) {
        sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
      }else{
        sample.means = cbind(sample.means, cpm[, kk])
      }
      
    }  
    colnames(sample.means) = conds
    
    ## calculate the correlation with TF expression 
    ggs = sapply(rownames(sample.means), function(x){unlist(strsplit(as.character(x), '_'))[1]})
    mm = match(rownames(bb), ggs)
    
    bb =data.frame(bb, stringsAsFactors = FALSE)
    bb$gene = rownames(bb)
    bb$combine.Zscore = r$combined.Zscore[match(rownames(bb), names(r$combined.Zscore))]
    bb$rank = c(1:nrow(bb))
    bb$cor = NA
    bb$tf.index = mm
    for(n in 1:nrow(bb))
    {
      if(!is.na(mm[n])) {
        bb$cor[n] = cor(as.numeric(bb[n, c(1:5)]), sample.means[mm[n], ]) 
      }
    }
    
    xx = bb[which(!is.na(bb$cor)), ]
    
    xx$activity = NA
    xx$activity[xx$cor>0] = 'activator'
    xx$activity[xx$cor<0] = 'repressor'
    
    
    require(ggplot2)
    require(tidyverse)
    library(ggrepel)
    
    ggplot(xx, aes(x=rank, y=cor, label = gene)) + 
      geom_point(aes(x=rank, y=cor, color = activity), size = 3.0) + 
      scale_color_manual(values=c('blue', 'red')) +
      geom_text_repel(data= xx , size = 4.0) +
      theme_classic() +
      geom_hline(yintercept=0.0, colour = "darkgray", size = 1.5) +
      xlab("motif activity importance rank") + ylab("Correlation to TF expression") +
      theme(axis.text.x = element_text(size = 14), 
            axis.text.y = element_text(size = 14), 
            axis.title.x = element_text(size=14, face="bold"),
            axis.title.y = element_text(size=14, face="bold")) 
     
    ggsave(paste0(resDir, '/MARA_bayesianRidge_temporalpeaks_motifActivity.rank_vs_TFexpression.pdf'), width=8, height = 6)
    
    ggplot(xx, aes(x=combine.Zscore, y=cor, label = gene)) + 
      geom_point(aes(x=combine.Zscore, y=cor, color = activity), size = 3.0) + 
      scale_color_manual(values=c('blue', 'red')) +
      geom_text_repel(data= xx , size = 4.0) +
      theme_classic() +
      geom_hline(yintercept=0.0, colour = "darkgray", size = 1.5) +
      xlab("Motif activity importance (z-score)") + ylab("Correlation to TF expression") +
      theme(axis.text.x = element_text(size = 14), 
            axis.text.y = element_text(size = 14), 
            axis.title.x = element_text(size=14, face="bold"),
            axis.title.y = element_text(size=14, face="bold")) 
    
    ggsave(paste0(resDir, '/MARA_bayesianRidge_temporalpeaks_motifActivity_vs_TFexpression.pdf'), width=8, height = 6)
    
    
    pdfname = paste0(resDir, "/motifActivity_TFexpr_comparison_2.pdf")
    pdf(pdfname, width = 5, height = 5)
    par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(6,3,2,0.2), tcl = -0.3)
    
    jj = 16
    y0 = scale(as.numeric(xx[jj, c(1:5)]))
    y1 = scale(as.numeric(sample.means[xx$tf.index[jj], ]))
    
    plot(c(1:5), y0, col = 'darkgreen', type = 'l', ylim = range(c(y0, y1)), main = xx$gene[jj], lwd = 2.5, xlab = '', ylab = '',
         xaxt="n")
    points(c(1:5), y0, col = 'darkgreen', type = 'p', cex = 2.0, pch = 16)
    points(c(1:5), y1, col = 'darkblue', type = 'l', lwd = 2.0)
    points(c(1:5), y1, col = 'darkblue', type = 'p', cex = 2.0, pch = 16)
    legend('topright', lwd = c(6, 6), col = c('darkgreen', 'darkblue'), legend = c('motif activity', 'TF expression'), bty = 'n' )
    abline(h = 0, lty = 1, col = 'darkgray', lwd = 2.0)
    axis(1, at=c(1:5), labels= c('mUA', '5dpa', '9dpa', '13dpa_prox', '13dpa.dist'), 
          las=1, lty = 2.0, cex.axis = 1.2)
    
    dev.off()
    
    
    
  }
  
  Run.glmnet = FALSE
  if(Run.glmnet){
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
    
    alpha = 0.2
    standardize = TRUE;
    use.lambda.min = FALSE;
    binarize.x = TRUE
    standardize.response=FALSE
    intercept=FALSE
    family = 'mgaussian'
    lambdas = 10^(seq(-4, 2, length.out = 100))
    
    x = X;
    y = scale(Y, center = TRUE, scale = FALSE); # center response and standardize X (cf. elastic-net paper)
    if(binarize.x) x = x > 0
    #y = y[, c(2, 3)]
    
    library(doMC) 
    registerDoMC(cores=6)
    
    library(tictoc)
    tic()
    cv.fit=cv.glmnet(x, y, family= family,
                     alpha=alpha, lambda = lambdas, standardize=standardize, 
                     standardize.response=standardize.response, parallel = TRUE)
    
    plot(cv.fit)
    toc()
    
    
    
    #library(ridge)
    #fit2 = linearRidge(y[,1] ~ x)
    
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
    
    #s.optimal = fit$lambda[16]
    
    xx = coef(fit, s = s.optimal)
    aa = c();  for(j in 1:length(xx)) { aa = cbind(aa, as.numeric(xx[[j]])); }
    
    rownames(aa) = rownames(xx[[1]])
    colnames(aa) = colnames(y)
    aa = as.data.frame(aa[-1, ]) # ignore the intercept
    aa = scale(aa) # calculate z score here is probably not a good idea to do feature selection
    #rownames(aa) = rownames(xx[[1]])[-1]
    
    cutoff.activity = 2.5
    ss = apply(aa, 1, function(x) length(which(abs(x) > cutoff.activity)))
    cat(length(which(ss>0)), ' nzero features \n')
    
    print(aa[which(ss>0), ])
    
    
    pheatmap(aa[which(ss>0), ], cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
             scale = 'none', cluster_cols=FALSE, main = paste0("motif activity by MARA"), 
             na_col = "white", fontsize_col = 12) 
    
  }
  
}


