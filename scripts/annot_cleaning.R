##########################################################################
##########################################################################
# Project:
# Script purpose: clean GTF file of axolotl
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Feb 24 15:30:59 2021
##########################################################################
##########################################################################
rm(list=ls())

##########################################
# functions used to clean the gtf 
##########################################
extract.transcript.gene.gtf.attribute = function(test)
{
  # test =keep$V9[24192]
  saved = rep(NA, 8)
  
  names(saved) = c('geneid', 'transcriptid', 'geneSymbol.hs', 'geneSymbol.nr', 'homolog.info', 'ORF.type', 
                   'attributes.gtf.transcript', 'CDS')
  test = as.character(test)
  
  test = unlist(strsplit(as.character(test), '; '))
  jj = which(test == ' '); if(length(jj) >0 ) test = test[-jj]
  
  test = test[grep('peptide ', test, invert = TRUE)]
  
  jj.cds = grep('CDS ', test)
  if(length(jj.cds) == 1){
    cds = test[jj.cds]
    test = test[-jj.cds]
    cds = gsub('CDS ', '', cds)
    saved[8] = cds
  }
  
  saved[7] = paste0(test, collapse = ';')
  
  gene.id = gsub('gene_id ', '', test[grep('gene_id', test)])
  if(length(gene.id) ==1) saved[1] = gene.id 
  
  hlg = gsub('homolog ', '', test[grep('homolog', test)])
  if(length(hlg) == 1) saved[5] = hlg
  
  orf = gsub('ORF_type ', '', test[grep('ORF_type', test)])
  if(length(orf) == 1) saved[6] = orf 
  
  transcripts = gsub('transcript_id ', '', test[grep('transcript_id', test)])
  transcripts = unlist(strsplit(as.character(transcripts), '[|]'))
  
  amextranscript = transcripts[grep('AMEX60DD', transcripts)]
  if(length(amextranscript) == 1) {
    transcripts = setdiff(transcripts, amextranscript)
    saved[2] = amextranscript
  }
  
  
  name.nr = transcripts[grep('nr', transcripts)]
  name.hs = transcripts[grep('hs', transcripts)]
  
  if(length(name.nr)>0) saved[4] = gsub('nr','', gsub("\\[|\\]", "", name.nr))
  if(length(name.hs)>0) saved[3] = gsub('hs','', gsub("\\[|\\]", "", name.hs))
  
  if(length(name.hs) == 0 & length(name.nr) == 0 & length(transcripts) == 1){
    if(transcripts != gene.id & transcripts != amextranscript & transcripts != ''){
      saved[3] = transcripts
      saved[4] = transcripts
    }
  }
  
  return(saved)
  
}

extract.gene.id.gene.names = function(test)
{
  # test = keep$V9[5]
  save = rep(NA, 4)
  names(save) = c('gene.id', 'gene.symbol.nr', 'gene.symbol.hs', 'others')
  
  test = as.character(test)
  test = unlist(strsplit(as.character(test), ';'))
  jj = which(test == ' ')
  if(length(jj) >0 ) test = test[-jj]
  
  ii.id = grep('gene_id', test)
  ii.name = grep('gene_name', test)
  ii.other = setdiff(c(1:length(test)), unique(c(ii.id, ii.name)))
  if(length(ii.other)>0) save[4] = paste0(test[ii.other], collapse = ';')
  
  gene.id = test[ii.id]
  gene.id = gsub('gene_id ', '', gene.id)
  save[1] = gene.id
  gene.names = test[ii.name]
  gene.names = gsub('gene_name ', '', gene.names)
  gene.names = gsub(' ', '', gene.names)
  gene.names = unlist(strsplit(as.character(gene.names), '[|]'))
  name.nr = gene.names[grep('nr', gene.names)]
  name.hs = gene.names[grep('hs', gene.names)]
  
  if(length(name.nr)>0) save[2] = gsub('nr','', gsub("\\[|\\]", "", name.nr))
  if(length(name.hs)>0) save[3] = gsub('hs','', gsub("\\[|\\]", "", name.hs))
  if(length(name.hs) == 0 & length(name.nr) == 0){
    if(gene.names != gene.id){
      save[2] = gene.names
      save[3] = gene.names
    }
  }
  
  return(save)
  
}


find.geneId.contigName.homolog.fromGeneSymbolFile = function(test)
{
  saved = rep(NA, 5)
  
  #test = ggs[n, ]
  test = as.character(test)
  test = gsub('>', '', test)
  test = unlist(strsplit(test, '[|]'))
  
  if(length(test) == 1){
    saved[1] = test   
  }
  if(length(test) == 2){
    saved[1] = test[1]
    saved[2] = test[2]
  }
  
  # three columns
  if(length(test) == 3){
    saved[1] = test[1]
    saved[2] = test[2]
    saved[5] = test[3]
    
    if(grepl('[hs;nr]', test[3], fixed = TRUE)){
      saved[3] = gsub(' hs;nr','', gsub("\\[|\\]", "", test[3]))
      saved[4] = gsub(' hs;nr','', gsub("\\[|\\]", "", test[3]))
    }else{
      if(grepl('[nr;hs]', test[3], fixed = TRUE)){
        #cat(n, '\n')
        saved[3] = gsub(' nr;hs','', gsub("\\[|\\]", "", test[3]))
        saved[4] = gsub(' nr;hs','', gsub("\\[|\\]", "", test[3]))
      }else{
        if(grepl('[nr]', test[3], fixed = TRUE)){
          #cat(n, '\n')
          saved[3] = gsub(' nr','', gsub("\\[|\\]", "", test[3]))
        }else{
          if(grepl('[hs]', test[3], fixed = TRUE)){
            #cat(n, '\n')
            saved[4] = gsub(' hs','', gsub("\\[|\\]", "", test[3]))
          }else{
            cat(n, 'Error \n')
            break
          }
        }
      }
    }
    #break
  }
  
  # four columns
  if(length(test) == 4) {
    saved[1] = test[1]
    saved[2] = test[2]
    saved[5] = paste0(test[3], '|', test[4])
    
    gene.names = test[c(3,4)]
    name.nr = gene.names[grep('[nr]', gene.names)]
    name.hs = gene.names[grep('[hs]', gene.names)]
    if(length(name.nr)>0) {
      saved[3] = gsub(' nr','', gsub("\\[|\\]", "", name.nr))
    }else{
      cat('Error \n')
    }
    if(length(name.hs)>0) {
      saved[4] = gsub(' hs','', gsub("\\[|\\]", "", name.hs))
    }else{
      cat('Error \n')
    }
    #break
  }
  
  # > 4 columns Error 
  if(length(test)>4){
    cat('> 4 column found !\n')
    break;
  }
  
  return(saved)
  
}

########################################################
########################################################
# Section : start with the mexT_v47.with_genesymbols.txt
# 
########################################################
########################################################
Save.annot = FALSE
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

##########################################
# process mapping of transcript.id - contig.id - gene.symbol.nr - gene.symbol.hs
##########################################
ggs = read.delim(file = paste0(annotDir, 'AmexT_v47.with_genesymbols.txt'), sep = '\t', header = FALSE)
ggs$V1 = as.character(ggs$V1)


xx = data.frame(t(sapply(ggs$V1, find.geneId.contigName.homolog.fromGeneSymbolFile)),
                stringsAsFactors = FALSE)
colnames(xx) = c('transcritID', 'transcript.contig', 'homolog.nr', 'homolog.hs', 'homolog.orig')

saved = xx
rm(xx)

length(saved$homolog.hs)
length(saved$homolog.hs[which(!is.na(saved$homolog.hs))])
length(saved$homolog.hs[which(!is.na(saved$homolog.hs) & saved$homolog.hs != 'N/A')])
length(which(!is.na(saved$homolog.nr) & saved$homolog.nr != 'N/A'))

length(which((!is.na(saved$homolog.nr) & saved$homolog.nr != 'N/A')| (!is.na(saved$homolog.hs) & saved$homolog.hs != 'N/A')))


hs.uniq = unique(saved$homolog.hs[which(!is.na(saved$homolog.hs) & saved$homolog.hs != 'N/A')])
length(unique(saved$homolog.hs[which(!is.na(saved$homolog.hs) & saved$homolog.hs != 'N/A')]))

saved$info.orig = rownames(saved)

if(Save.annot){
  saveRDS(saved, file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs.rds'))
  
  saved = readRDS(file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs.rds'))
  
  xx = saved
  xx$geneID = xx$transcritID
  xx$geneID = sapply(xx$geneID, function(x) unlist(strsplit(as.character(x), '[.]'))[1])
  
  # geneIDs.from.transcriptI
  jj = grep('AMEX60DDU', xx$geneID, invert = TRUE)
  xx$geneID[jj] = sapply(xx$geneID[jj], function(x) {paste0(unlist(strsplit(as.character(x),''))[-c(9:11)], collapse = '')})
  
  saved = xx
  
  saveRDS(saved, file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID.rds'))
  
}


########################################################
########################################################
# Section : collect mapping from GTF file
# 
########################################################
########################################################
annot = readRDS(file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID.rds'))

aa = read.delim(file = paste0(annotDir, 'AmexT_v47.release.gtf'), sep = '\t', header = FALSE)

# remove contigs 
if(Save.annot){
  keep = aa
  keep = keep[grep('chr', keep$V1), ]
  write.table(keep, file = paste0(annotDir, 'AmexT_v47.release_rm.contigs.gtf'), sep = '\t', quote = FALSE, col.names = FALSE, 
              row.names = FALSE)
}

##########################################
# # select genes and make genn id to gene symbol mapping
##########################################
keep = aa

keep = keep[which(keep$V3 == 'gene'), ]
keep$V9 = as.character(keep$V9)

mapping = t(sapply(keep$V9, extract.gene.id.gene.names))

mapping = data.frame(keep[, c(1:9)],  mapping, stringsAsFactors = FALSE)

length(unique(mapping$gene.symbol.hs))

xx = annot

mm = match(xx$geneID, mapping$gene.id)

cat(length(which(is.na(mm))), ' genes missing ')

mapping = mapping[, c(1, 4, 5, 7, 9:13)]
colnames(mapping)[1:5] = c('chr', 'start', 'end', 'strand', 'gtf.attributes')
colnames(mapping) = paste0(colnames(mapping), '_gene')

xx = data.frame(xx, mapping[mm, ], stringsAsFactors = FALSE)

annot = xx
rm(xx)

jj = which(!is.na(annot$homolog.hs) & annot$homolog.hs!= 'N/A' & is.na(annot$gene.symbol.hs_gene))

if(Save.annot){
  
  saveRDS(annot, file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo..rds'))
  
  write.table(mapping, file = paste0(annotDir, 'AmexT_v47.release_rm.contigs_geneID_geneSymbol.mapping.txt'), sep = '\t', 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
}

##########################################
# select transcripts and make transcript annotation 
##########################################
annot = readRDS(file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo..rds'))

keep = aa
keep$V9 = as.character(keep$V9)

#keep = keep[grep('chr', keep$V1), ]
keep = keep[which(keep$V3 == 'transcript'), ]

tgmap = data.frame(t(sapply(keep$V9, extract.transcript.gene.gtf.attribute)))
rownames(tgmap) = NULL

mm = match(annot$transcritID, tgmap$transcriptid)
cat(length(which(is.na(mm))), ' transcripts missing \n') # happy that 0 transcripts missing in the cross check
colnames(tgmap) = paste0(colnames(tgmap), '_gtf.transcript')

xx = data.frame(annot, tgmap[mm, ],  keep[mm, c(1, 4, 5, 7)],  stringsAsFactors = FALSE)
rownames(xx) = NULL
colnames(xx)[c(25:28)] = c('chr_transcript', 'start_transcript', 'end_transcript', 'strand_transcript') 

for(n in 1:ncol(xx))
{
  xx[,n] = as.character(xx[,n])
}

xx = data.frame(xx, stringsAsFactors = FALSE)
xx$geneSymbol.hs_gtf.transcript = gsub(' ', '', as.character(xx$geneSymbol.hs_gtf.transcript))

length(unique(xx$homolog.hs))
length(unique(xx$gene.symbol.hs_gene))
length(unique(xx$geneSymbol.hs_gtf.transcript))


symbols = data.frame(xx$homolog.hs, xx$gene.symbol.hs_gene, xx$geneSymbol.hs_gtf.transcript, xx$transcritID,  
                     stringsAsFactors = FALSE)

kk = which(symbols[,1] != symbols[,2] | symbols[,3] != symbols[,2] | symbols[,1] != symbols[,3] )

colnames(xx)[1:6] = c('transcriptID', 'transcriptID.contig', "homolog.nr.fa", 
                      "homolog.hs.fa",  "homolog.orig.fa",  "info.orig.fa") 
colnames(xx)[c(12:16)] = c('attributes_gtf.gene', 'geneID_gtf.gene', 
                           "gene.symbol.nr_gtf.gene", "gene.symbol.hs_gtf.gene",  "others_gtf.gene")

annot = xx

saveRDS(annot, 
        file = paste0(annotDir, 
                      'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))

write.table(annot, 
            file = paste0(annotDir, 'annot_AmexT_v47_transcriptID_transcriptCotig_geneID_geneSymbol.hs_nr_CDS.txt'),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


########################################################
########################################################
# Section : explore the gene synteny and improve the gene symbol
# 
########################################################
########################################################
Improve.hs.symbol.annotation.with.synteny.analysis = FALSE
if(Improve.hs.symbol.annotation.with.synteny.analysis){
  
  rm(list=ls())
  
  annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
  resDir = '../results/gene_annotations'
  if(!dir.exists(resDir)) dir.create(resDir)
  
  ##########################################
  # prepare and clean human protein coding genes
  ##########################################
  hs = read.delim(paste0(annotDir, 'BioMart_ens_GRCh38.p13.txt'), header = TRUE)
  hs = hs[which(hs$Gene.type == 'protein_coding'), ]
  
  name.uniq = unique(hs$Gene.stable.ID)
  hs = hs[match(name.uniq, hs$Gene.stable.ID), ]
  
  kk = grep('GL|KI|CHR_|MT', hs$Chromosome.scaffold.name, invert = TRUE)
  
  hs = hs[kk, ]
  hs = data.frame(hs$Gene.stable.ID, hs$Gene.start..bp., hs$Gene.end..bp., hs$Chromosome.scaffold.name, hs$Gene.name, 
                  stringsAsFactors = FALSE)
  colnames(hs) = c('geneID', 'start', 'end', 'chr', 'gene.name')
  
  hs = hs[with(hs, order(chr, start, end)), ] ## 19920 genes in chr1-22, X,Y
  
  hs$gene.name = as.character(hs$gene.name)
  
  ##########################################
  # prepare axolotl genes 
  ##########################################
  annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
  annot = readRDS(file = 
                    paste0(annotDir, 
          'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
  
  Summary.gene.annotation = FALSE
  if(Summary.gene.annotation){
    
    counts = table(annot$ORF.type_gtf.transcript)
    barplot(counts, las = 2, col = c(1:length(counts)), main = 'ORF.types of 181985 transcripts')
    
    for(n in 1:length(counts))
    {
      jj = which(annot$ORF.type_gtf.transcript == names(counts)[n])
      counts[n] = length(unique(annot$geneID[jj]))
    }
    
    barplot(counts, las = 2, col = c(1:length(counts)), main = 'ORF.types of 99088 genes')
    
    ccts = table(annot$ORF.type_gtf.transcript)
    counts = matrix(NA, nrow = 5, ncol = length(ccts))
    rownames(counts) = c('total', 'total.chr', 'hs', 'nr', 'hs.protein.coding.uniq')
    colnames(counts) = names(ccts)
    counts[1, ] = ccts
    
    for(n in 1:ncol(counts))
    {
      # n = 10
      jj = which(annot$ORF.type_gtf.transcript == colnames(counts)[n])
      counts[1, n] = length(unique(annot$geneID[jj]))
      
      jj = which(annot$ORF.type_gtf.transcript == colnames(counts)[n] & grepl('^chr', annot$chr_transcript))
      counts[2, n] = length(unique(annot$geneID[jj]))
      
      jj = which(annot$ORF.type_gtf.transcript == colnames(counts)[n] & grepl('^chr', annot$chr_transcript) 
                 & !is.na(annot$geneSymbol.hs_gtf.transcript) & annot$geneSymbol.hs_gtf.transcript != 'N/A' &
                   annot$geneSymbol.hs_gtf.transcript != ' ')
      counts[3, n] = length(unique(annot$geneID[jj]))
      
      mm = match(annot$geneSymbol.hs_gtf.transcript[jj], hs$gene.name)
      counts[5, n] = length(unique(hs$gene.name[mm[!is.na(mm)]]))
      
      jj = which(annot$ORF.type_gtf.transcript == colnames(counts)[n] & grepl('^chr', annot$chr_transcript) 
                 & !is.na(annot$geneSymbol.nr_gtf.transcript) & annot$geneSymbol.nr_gtf.transcript != 'N/A' &
                   annot$geneSymbol.nr_gtf.transcript != ' ')
      counts[4, n] = length(unique(annot$geneID[jj]))
      
     
      # jj = which(annot$ORF.type_gtf.transcript == colnames(counts)[n] & grepl('^chr', annot$chr_transcript) 
      #            & !is.na(annot$geneSymbol.hs_gtf.transcript) & annot$geneSymbol.hs_gtf.transcript != 'N/A' &
      #              annot$geneSymbol.hs_gtf.transcript != ' ')
      # counts[5, n] = length(unique(annot$geneSymbol.hs_gtf.transcript[jj]))
      # 
      # jj = which(annot$ORF.type_gtf.transcript == colnames(counts)[n] & grepl('^chr', annot$chr_transcript) 
      #            & !is.na(annot$geneSymbol.nr_gtf.transcript) & annot$geneSymbol.nr_gtf.transcript != 'N/A' &
      #              annot$geneSymbol.nr_gtf.transcript != ' ')
      # counts[6, n] = length(unique(annot$geneSymbol.nr_gtf.transcript[jj]))
      
    }
    
    barplot(counts, col=c(1:nrow(counts)), legend = rownames(counts), beside=TRUE, las = 2, 
            main = 'counts of genes with hs (14000) and nr gene symbols ')
    
  }
  
  ##########################################
  ## subset the annotation only for putative full lengnth genes and human protein coding genes
  ##########################################
  mm = match(annot$geneSymbol.hs_gtf.transcript, hs$gene.name)
  cat(length(unique(hs$gene.name[mm[!is.na(mm)]])), ' human protein coding genes found in all aloxlotl annotated genes \n')
  
  mm = match(annot$geneSymbol.hs_gtf.transcript[which(annot$ORF.type_gtf.transcript == 'Putative')], hs$gene.name)
  cat(length(unique(hs$gene.name[mm[!is.na(mm)]])), ' human protein coding genes found in axolotl genes annotated as putative \n')
  
  mm = match(annot$geneSymbol.hs_gtf.transcript[which(annot$ORF.type_gtf.transcript == 'Putative'| 
                                                        !is.na(annot$geneSymbol.hs_gtf.transcript))], hs$gene.name)
  cat(length(unique(hs$gene.name[mm[!is.na(mm)]])), 
      ' human protein coding genes found in axolotl genes annotated as putative and with hs \n')
  
  ax = annot[which(annot$ORF.type_gtf.transcript == 'Putative'| !is.na(annot$geneSymbol.hs_gtf.transcript)), ]
  
  
  ax = data.frame(ax$geneID, ax$gene.symbol.hs_gtf.gene, ax$gene.symbol.nr_gtf.gene,  ax$chr_gene,
                  ax$start_gene, ax$end_gene, ax$transcriptID, ax$geneSymbol.hs_gtf.transcript, ax$geneSymbol.nr_gtf.transcript,
                  ax$ORF.type_gtf.transcript,
                  stringsAsFactors = FALSE)
  colnames(ax) = c('geneID', 'hs.symbol.gene', 'nr.symbol.gene', 'gene.chr', 'gene.start', 'gene.end', 'transcriptID', 
                   'hs.symbol.transcript', 'nr.symbol.transcript',  'ORF.type')
  
  ax = ax[grep('chr', ax$gene.chr), ]
  
  ax = ax[with(ax, order(gene.chr, gene.start, gene.end)), ]
  
  ax = ax[, c(1, 7, 4:6, 2:3, 8:10)]
  #colnames(ax)[c(6,7)] = c('hs.symbol', 'nr.symbol')
  
  gs = data.frame(ax[match(unique(ax$geneID), ax$geneID), ], stringsAsFactors = FALSE)
  gs = gs[, c(1, 3:7, 10)]
  gs = gs[with(gs, order(gene.chr, gene.start, gene.end)), ]
  colnames(gs) = c('geneID', 'chr', 'start', 'end', 'gene.symbol.hs', 'gene.symbol.nr', 'ORF.type')
  
  # double check the gene symbols from hs and nr
  gs$hs.new = NA
  gs$nr.new = NA
  gs$hs.nr.new = NA
  
  for(n in 1:nrow(gs))
  {
    # n = 1
    if(n%%1000 == 0) cat(n, '\n')
    ii = which(annot$geneID == gs$geneID[n])
    
    if(length(ii) == 0){cat("Error ! -- ", n, '\n')}
    
    orftypes = annot$ORF.type_gtf.transcript[ii]
    orftypes = orftypes[!is.na(orftypes)]
    
    if(any(orftypes == 'Putative')) {
      gs$ORF.type[n] = 'Putative'
    }
    
    if(length(ii) == 1){
      gs$hs.new[n] = annot$geneSymbol.hs_gtf.transcript[ii];
      gs$nr.new[n] = annot$geneSymbol.nr_gtf.transcript[ii];
      if(grepl(' ', gs$hs.new[n])) gsub(' ', '', gs$hs.new[n])
      if(grepl(' ', gs$nr.new[n])) gs$nr.new[n] = gsub(' ', '', gs$nr.new[n])
    }
    
    if(length(ii)>1){
      #cat(n, '\n')
      hs.test = unique(annot$geneSymbol.hs_gtf.transcript[ii])
      nr.test = unique(annot$geneSymbol.nr_gtf.transcript[ii])
      hs.test = hs.test[!is.na(hs.test)]
      nr.test = nr.test[!is.na(nr.test)]
      
      if(length(hs.test) >= 1 ){
        gs$hs.new[n] = paste0(hs.test, collapse = ';')
        if(grepl(' ', gs$hs.new[n])) gsub(' ', '', gs$hs.new[n])
      }
      if(length(nr.test) >= 1 ){
        gs$nr.new[n] = paste0(nr.test, collapse = ';')
        if(grepl(' ', gs$nr.new[n])) gs$nr.new[n] = gsub(' ', '', gs$nr.new[n])
      }
    }
  }
  
  jj = which(gs$hs.new == gs$nr.new)
  gs$hs.nr.new[jj] = gs$hs.new[jj]
  gs$hs.nr.new[-jj] = paste0(gs$hs.new[-jj], ';', gs$nr.new[-jj])
  
  gs = data.frame(gs, stringsAsFactors = FALSE)
  gs$start = as.numeric(gs$start)
  gs$end = as.numeric(gs$end)
  gs = gs[with(gs, order(chr, start)), ]
  
  keep = gs; 
  
  nb.hs = 2
  nb.search = 4
  gs$hs.synteny = NA
  gs$synteny.left = NA
  gs$synteny.right = NA
  
  for(n in 10:(nrow(gs)-10))
  {
    #if(n%%1000 == 0) cat(n, '\n')
    # n  = 167
    ggs = gs$hs.new[n]
    ggs = unique(unlist(strsplit(as.character(ggs), ';')))
    ggs = ggs[!is.na(ggs) & ggs != 'N/A']
    
    if(length(ggs)>0){
      for(gg in ggs)
      {
        # gg = ggs[1]
        hs.kk = which(hs$gene.name == gg)
        
        if(length(hs.kk) ==1){
          hs.left = as.character(hs$gene.name[seq((hs.kk-nb.hs), (hs.kk-1))])
          hs.right = as.character(hs$gene.name[seq((hs.kk+1), (hs.kk+nb.hs))])
          #print(hs.left)
          #print(hs.right)
          
          query.left = gs$hs.nr.new[seq(n-nb.search, n-1)]
          query.right = gs$hs.nr.new[seq(n+1, n + nb.search)]
          
          query.left = as.character(unlist(sapply(query.left, function(x) unlist(strsplit(as.character(x), ';'))))) 
          query.left = unique(query.left[!is.na(query.left) & query.left != 'NA'])
          query.right = as.character(unlist(sapply(query.right, function(x) unlist(strsplit(as.character(x), ';'))))) 
          query.right = unique(query.right[!is.na(query.right) & query.right != 'NA'])
          
          ll = intersect(query.left, hs.left)
          rr = intersect(query.right, hs.right)
          
          lr = intersect(query.left, hs.right)
          rl = intersect(query.right, hs.left)
          
          if((length(ll) >=1 & length(rr) >=1) | (length(lr) >=1 & length(rl) >=1) | 
             length(ll) >=2 | length(rr) >=2 | length(rl) >= 2 | length(rl) >=2 ) {
            cat(n, ' --', gg, '\n')
            if(is.na(gs$hs.synteny[n])){
              gs$hs.synteny[n] = gg
              gs$synteny.left[n] = paste0(unique(c(ll, lr)), collapse = ';')
              gs$synteny.right[n] = paste0(unique(c(rr, rl)), collapse = ';')
              
            }else{
              gs$hs.synteny[n] = paste0(gs$hs.synteny[n], ';',  gg)
              gs$synteny.left[n] = paste0(gs$synteny.left[n], '|', paste0(unique(c(ll, lr)), collapse = ';'))
              gs$synteny.right[n] = paste0(gs$synteny.left[n], '|', paste0( unique(c(rr, rl)), collapse = ';'))
            }
          }
        }
      }  
    }
  }
  
  
 
  colnames(gs)[11] = 'gene.evidence.synteny'
  
  gs$gene.evidence.same.gene.symbol.hs.nr = NA
  for(n in 1:nrow(gs)){
    # n = 1
    xx = unlist(strsplit(as.character(gs$hs.new[n]), ';'))
    yy = unlist(strsplit(as.character(gs$nr.new[n]), ';'))
    xy = intersect(xx, yy)
    if(length(xy)>0) gs$gene.evidence.same.gene.symbol.hs.nr[n] = paste0(xy, collapse = ';')
      
  }
  
  
  length(unique(gs$gene.evidence.synteny))
  length(unique(gs$gene.evidence.same.gene.symbol.hs.nr))
  
  mm1 = match(unique(gs$gene.evidence.synteny), hs$gene.name)
  mm1 = mm1[!is.na(mm1)]
  length(unique(hs$gene.name[mm1]))
  
  mm2 = match(unique(gs$gene.evidence.same.gene.symbol.hs.nr), hs$gene.name)
  mm2 = mm2[!is.na(mm2)]
  length(unique(hs$gene.name[mm2]))
  
  length(unique(hs$gene.name[c(mm1, mm2)]))
  
  counts = c(length(unique(hs$gene.name[mm1])), length(unique(hs$gene.name[mm2])), length(unique(hs$gene.name[c(mm1, mm2)])))
  names(counts) = c('synteny', 'shared.gene.symbol.hs.nr', 'synteny.or.shared.geneSymbol')
  
  barplot(counts, col = c(1:length(counts)), main = '# of protein coding genes with extra evidence')
  
  saveRDS(gs, file = paste0(annotDir, 'geneAnnotation_geneSymbols.synteny.evidence.hs.nr.same.evidence.rds'))
  
  write.table(gs,
              file = paste0(resDir, 'annot_AmexT_v47_geneAnnot_extraEvidence_synteny.shared.geneSymbol.hs.nr.txt'),
              sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
    
}


Test.annotation.v1 = FALSE
if(Test.annotation.v1){
  annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                         'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
  xx = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                      'geneAnnotation_geneSymbols.synteny.evidence.hs.nr.same.evidence.rds'))
  mm = match(annot$geneID, xx$geneID)
  
  
}


########################################################
########################################################
# Section : some tests at the first place
#  not used anymore
########################################################
########################################################
Initial.test = FALSE
if(Initial.test){
  library("tictoc")
  tic()
  
  res = t(sapply(keep$V9, extract.gene.transcript.gtf.attribute))
  
  colnames(res) = c('gene.names', 'transcript.names', 'attributes.gtf', 'gene.id.Amex60DD', 'transcript.id.Amex60DD', 'other.infos')  
  
  toc()
  
  
  keep$V9 = res[,3]
  keep = keep[grep('chr', keep$V1), ]
  write.table(keep, file = paste0(annotDir, 'ax6_UCSC_2021_01_26_attributes.cleaned_rm.contigs.gtf'), sep = '\t', 
              quote = FALSE, col.names = FALSE, 
              row.names = FALSE)
  
  
  res = data.frame(keep[, c(1:7)], res, stringsAsFactors = FALSE)
  
  xx = res[grep('chr', res$V1),  ]
  
  Test.for.loop = FALSE
  if(Test.for.loop){
    
    library("tictoc")
    tic()
    for(n in 1:9000){
      #n = 1
      if(n%%1000 == 0) cat(n, '\n')
      
      test = keep$V9[n]
      test = unlist(strsplit(as.character(test), ';'))
      jj = which(test == ' ')
      if(length(jj) >0 ) test = test[-jj]
      ii.gene = grep('gene_id', test)
      ii.transcript = grep('transcript_id', test)
      ii.other = setdiff(c(1:length(test)), unique(c(ii.gene, ii.transcript)))
      
      gene.id = test[ii.gene]
      gene.id = unlist(strsplit(as.character(gene.id), '[|]'))
      keep$gene[n] = gsub('gene_id ', '', paste0(gene.id[grep('AMEX6', gene.id, invert = TRUE)], collapse = ';'))
      
      transcript.id = test[ii.transcript]
      transcript.id = unlist(strsplit(as.character(transcript.id), '[|]'))
      keep$transcript[n] = gsub('transcript_id ', '', paste0(gene.id[grep('AMEX6', transcript.id, invert = TRUE)], collapse = ';'))
      
      keep$new.attr[n] = paste0('gene_id "', gene.id[grep('AMEX6', gene.id)], '"; transcript_id "', 
                                transcript.id[grep('AMEX6', transcript.id)], '";')
      
      if(length(ii.other)>0)  keep$others[n] = paste0(test[ii.other], collapse = ';')
      
    }
    
    toc()
    
  }
  
}
