##########################################################################
##########################################################################
# Project: positional memory
# Script purpose: scripts for revision
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Mar  6 13:29:22 2023
##########################################################################
##########################################################################
rm(list = ls())

require(ggplot2)
require(DESeq2)
require(dplyr)
require(gridExtra)
require(RColorBrewer)

figureDir = paste0('~/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/',
                   'DEVELOPMENTAL CELL/Review/plots_jingkui/') 
tableDir = paste0('~/Dropbox (VBC)/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/SupTables/')

##########################################
# gene examples of smartseq2 and microarray data
##########################################
res = readRDS(file = paste0("../results/microarray/Rdata/", 
                            'design_probeIntensityMatrix_probeToTranscript.geneID.geneSymbol_',
                            'normalized_geneSummary_limma.DE.stats_RARB.rds'))

## manual change the gene annotation MEIS3
rownames(res)[grep('AMEX60DD024424', rownames(res))]
rownames(res)[grep('AMEX60DD024424', rownames(res))] = 'MEIS3_AMEX60DD024424' 

ggs = sapply(rownames(res), function(x){unlist(strsplit(as.character(x), '_'))[1]})

qv.cutoff = 0.05
logfc.cutoff = 1
select = which(res$fdr.max> -log10(qv.cutoff) & abs(res$logFC.max)> 0)

select = which(res$adj.P.Val_mHand.vs.mLA < qv.cutoff & abs(res$logFC_mHand.vs.mLA) > logfc.cutoff|
                 res$adj.P.Val_mHand.vs.mUA < qv.cutoff & abs(res$logFC_mHand.vs.mUA) > logfc.cutoff |
                 res$adj.P.Val_mLA.vs.mUA < qv.cutoff & abs(res$logFC_mLA.vs.mUA) > logfc.cutoff )
# cat(length(select), ' positional genes found \n')
cat(length(select), ' DE genes selected \n')
ggs.sel = ggs[select]

### plot gene expression of microarray data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

gene_examples = c('HOXA13', 'HOXA9', 'HOXD13', 'HOXD9', 'MEIS1', 'MEIS3', 'SHOX2')

for(n in 1:length(gene_examples)){
 
  # n = 1
  kk = which(ggs == gene_examples[n])
  if(length(kk) == 1){
    test = res[kk, c(1:9)]
    test = data.frame(segs = sapply(names(test), function(x){unlist(strsplit(as.character(x), '_'))[1]}), 
                      reps = sapply(names(test), function(x){unlist(strsplit(as.character(x), '_'))[2]}), 
                      signals = as.numeric(test))
    test$segs = factor(test$segs, levels = c('mUA', 'mLA', 'mHand'))
    
    test = data_summary(test, varname = 'signals', groupnames = 'segs')
    
    # Default bar plot
    p<- ggplot(test, aes(x=segs, y=signals, fill=segs)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=signals-sd, ymax=signals+sd), width=.2,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0, size = 14), 
            axis.text.y = element_text(angle = 0, size = 14), 
            axis.title =  element_text(size = 14),
            legend.text = element_text(size=12),
            legend.title = element_text(size = 14)
      )+
      labs(x = "") +
      ggtitle(gene_examples[n])
    
    print(p)
    
    ggsave(paste0(figureDir, 'microarray_boxplot.errorbar_geneExample_', gene_examples[n], '.pdf'), 
           width=8, height = 6)
    
    
  }else{
    cat(length(kk), ' row found for gene: ', gene_examples[n], '\n')
  }
}




