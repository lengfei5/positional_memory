##########################################################################
##########################################################################
# Project: Akane's postional memory 
# Script purpose: Compare the development and mature samples
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Feb 23 09:59:18 2022
##########################################################################
##########################################################################
rm(list=ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('functions_chipSeq.R')
source('Functions_atac.R')


version.analysis = 'Rxxxx_R10723_R11637_R12810_atac'
#peakDir = "Peaks/macs2_broad"

resDir = paste0("../results/", version.analysis)

RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using/'
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
gtf.file =  paste0(annotDir, 'ax6_UCSC_2021_01_26.gtf')


figureDir = '/Users/jiwang/Dropbox/Group Folder Tanaka/Collaborations/Akane/Jingkui/Hox Manuscript/figure/plots_4figures/' 
tableDir = paste0(figureDir, 'tables4plots/')


require(ggplot2)
require(DESeq2)
require(GenomicRanges)
require(pheatmap)
library(tictoc)


load(file = paste0(RdataDir, '/samplesDesign_readCounts.within_manualConsensusPeaks.pval6_mergedTechnical_', 
                   version.analysis, '.Rdata'))
design$sampleID = design$SampleID
design$usable = as.numeric(design$usable)
design$usable[c(31:32)] = design$usable[c(31:32)]/10^6

saveRDS(design, file = paste0('../data/design_sampleInfos_32atacSamplesUsed.rds'))


require(ChIPpeakAnno)
require(ChIPseeker)

pp = data.frame(t(sapply(counts$gene, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
pp$strand = '*'

pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                              start.field="X2", end.field="X3", strand.field="strand")


##########################################
# Select peak consensus across mature, regeneration and embryo 
# and choose the background 
##########################################
Peaks.Background.selection = TRUE

if(Peaks.Background.selection){
  
  design$conds = design$condition
  design$unique.rmdup = design$usable
  colnames(design)[which(colnames(design) == 'fileName')] = 'samples'
  #sels = grep('Mature|Embryo|BL_UA', design$conds)
  sels = c(1:nrow(design))
  rownames(counts) = counts$gene
  counts = as.matrix(counts[, -1])
  
  dds <- DESeqDataSetFromMatrix(counts, DataFrame(design[sels, ]), design = ~ conds)
  colnames(dds) = colnames(counts)
  
  # check the peak length
  peakNames = rownames(dds)
  pp = data.frame(t(sapply(peakNames, function(x) unlist(strsplit(gsub('_', ':', as.character(x)), ':')))))
  
  pp$strand = '*'
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  ll = width(pp)
  
  ##########################################
  # filter peaks below certain thrshold of read counts
  # And also consider those filtered peaks as background
  ##########################################
  select.peaks.with.readThreshold = TRUE
  select.background.for.peaks = TRUE
  
  if(select.peaks.with.readThreshold){
    #ss = rowMax(counts(dds)[, grep('Embryo_', dds$conds)])
    ss = rowMaxs(counts(dds))/ll*500
    hist(log10(ss), breaks = 200, main = 'log2(max of read counts within peaks) ')
    cutoff.peak = 50 # 30 as peak cutoff looks good
    cutoff.bg = 20
    cat(length(which(ss >= cutoff.peak)), 'peaks selected with minimum read of the highest peak -- ', cutoff.peak,  '\n')
    cat(length(which(ss < cutoff.bg)), 'peaks selected with minimum read of the highest peak -- ', cutoff.bg,  '\n')
    abline(v= log10(cutoff.peak), col = 'red', lwd = 2.0)
    abline(v= log10(cutoff.bg), col = 'blue', lwd = 2.0)
    
    nb.above.threshold = apply(counts(dds), 1, function(x) length(which(x>cutoff.peak)))
    ii = which(ss >= cutoff.peak)
    #ii = which(nb.above.threshold>=2)
    
    if(select.background.for.peaks){
      ii.bg = which(ss < cutoff.bg)
      ii.bg = sample(ii.bg, size = 1000, replace = FALSE)
      rownames(dds)[ii.bg] = paste0('bg_', rownames(dds)[ii.bg])
      dds = dds[c(ii, ii.bg), ]
      ll.sels = ll[c(ii, ii.bg)]
      
    }else{
      dds <- dds[ii, ]
      ll.sels = ll[ss >= cutoff.peak]
    }
    
  }
  
  dds <- estimateSizeFactors(dds)
  
  plot(sizeFactors(dds), colSums(counts(dds))/median(colSums(counts(dds))), log = 'xy')
  
  plot(sizeFactors(dds), design$usable, log = 'xy')
  text(sizeFactors(dds), design$usable, labels = design$samples, cex = 0.7)
  
  save.scalingFactors.for.deeptools = FALSE
  if(save.scalingFactors.for.deeptools){
    xx = data.frame(sampleID = design$SampleID,  
                    scalingFactor = design$unique.rmdup/(sizeFactors(dds)*median(design$unique.rmdup)),
                    stringsAsFactors = FALSE)
    
    write.table(xx, file = paste0(resDir, '/DESeq2_scalingFactor_forDeeptools.txt'), sep = '\t',
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    sfs = data.frame(sample = colnames(dds), sf = sizeFactors(dds)*median(colSums(counts(dds))), stringsAsFactors = FALSE)
    
    saveRDS(sfs, file = paste0(RdataDir, '/DESeq2_peaks.based_scalingFactors_forGenomicRanger.rds'))
    
  }
}

##########################################
# test normalization and batch correction of ATAC-seq data
# TMM and combat were selected for normalization and batch correction
##########################################
source('Functions_atac.R')
library(edgeR)
require("sva")
require(limma)

Split.Mature.Regeneration.samples = TRUE
if(Split.Mature.Regeneration.samples){
  ##########################################
  # batch correct samples separately for mature samples and regeneration samples  
  ##########################################
  table(design$condition, design$batch)
  
  # batch correct development samples and mUA
  Batch.Correct.mUA.embryoStage = FALSE
  if(Batch.Correct.mUA.embryoStage){
    
    sels = unique(c(which((design$condition == 'Mature_UA'| design$condition == 'Mature_Hand') & 
                            (design$batch == '2020'| design$batch == '2021')), 
                    grep('Embryo_Stage', design$condition)))
    
    design.sels = design[sels, ]
    #design.sels = design.sels[which(design.sels$SampleID != '74938'), ]
    
    design.sels$conds = droplevels(design.sels$conds)
    
    #design.sels$batch[grep('749', design.sels$SampleID)] = '2019'
    #design.sels$batch = droplevels(design.sels$batch)
    table(design.sels$conds, design.sels$batch)
    
    ddx = dds[, sels]
    ddx$conds = droplevels(ddx$conds)
    ss = rowSums(counts(ddx))
    
    d <- DGEList(counts=counts(ddx), group=design.sels$conds)
    tmm <- calcNormFactors(d, method='TMM')
    tmm = cpm(tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
    
    bc = as.factor(design.sels$batch)
    mod = model.matrix(~ as.factor(conds), data = design.sels)
    
    # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
    # because there is no better batche other others 
    #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
    fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
    
    #design.tokeep<-model.matrix(~ 0 + conds,  data = design.sels)
    #cpm.bc = limma::removeBatchEffect(tmm, batch = bc, design = design.tokeep)
    # plot(fpm.bc[,1], tmm[, 1]);abline(0, 1, lwd = 2.0, col = 'red')
    make.pca.plots(tmm, ntop = 1000, conds.plot = 'all')
    ggsave(paste0(resDir, "/Dev_mUA_Samples_batchCorrect_before_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    
    make.pca.plots(fpm.bc, ntop = 1000, conds.plot = 'all')
    ggsave(paste0(resDir, "/Dev_mUA_Samples_batchCorrect_after_",  version.analysis, ".pdf"), width = 16, height = 14)
    
    
    fpm = fpm.bc
    
    rm(fpm.bc)
    
    saveRDS(fpm, file = paste0(RdataDir, '/fpm.bc_TMM_combat_Dev.mUA.samples_batch2020.2021.rds'))
    saveRDS(design.sels, file = paste0(RdataDir, '/design_sels_bc_TMM_combat_Dev.mUA.samples_batch2020.2021.rds'))
    
  }
  
}


########################################################
########################################################
# Section  : compare embryo samples to mature UA at the mRNA levels
# key questions: how do developmental genes required in embryo stages 
# look like ? lowly expressed in mature sample ? total not expressed 
# what about their promoter and enhancers in ATAC-seq data ? 
########################################################
########################################################

##########################################
# first import Tobias' scRNA-seq data 
##########################################
require(DESeq2)
library(ggrepel)
library(dplyr)
library(tibble)

# GO term annotation
limb.go = read.delim(file = '../data/limb_development_GO_term_summary_20220218_162116.txt', header = FALSE)
limb.go = unique(limb.go$V2[-1])
limb.go = toupper(limb.go)


##########################################
# use edgR for DE without replicates but with housekeeping gene
##########################################
load(file = paste0(RdataDir, '/pseudoBulk_scRNAcellPooling_FluidigmC1_stage40.44.mUA_dev_geneSelection.Rdata'))
require(edgeR)
group = dds$condition
y <- DGEList(counts=counts(dds), group=group)
y <- calcNormFactors(y)

y1 = y
y1$samples$group <- 1
#group1 = y1$samples$group
#design <- model.matrix(~group1)

index.hs = which(fpm$genetype == 'hs')
y0 <- estimateDisp(y1[index.hs,], trend = 'none', tagwise=FALSE)

y$common.dispersion <- y0$common.dispersion

fit <- glmFit(y, design)
lrt2 <- glmLRT(fit, coef = 2)
lrt3 = glmLRT(fit, coef = 3)

lrt2 = as.data.frame(topTags(lrt2, adjust.method = "BH", n = nrow(y)))
lrt3 = as.data.frame(topTags(lrt3, adjust.method = "BH", n = nrow(y)))

DE.genes = unique(c(rownames(lrt2)[which(lrt2$FDR<0.01 & abs(lrt2$logFC) > 1)],
                  rownames(lrt3)[which(lrt3$FDR <0.01 & abs(lrt3$logFC) >1)]))

fpm$DEgene = NA
fpm$DEgene[!is.na(match(rownames(fpm), DE.genes))] = '1'

fpm = data.frame(fpm)
fpm = fpm[which(fpm$genetype != 'ctl'), ]

fpm$genetype[which(fpm$genetype == 'dev.lowlyExp')] = 'Expr.E40.44'
fpm$genetype[which(fpm$genetype == 'devGene')] = 'Expr.E40.44'

fpm$gene =  sapply(rownames(fpm), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[1])})


dev.example = c('HOXA13', 'HOXA11', 'HOXA9', 'HOXD13','HOXD11', 'HOXD9',
                'SHH', 'FGF8', 'FGF10', 'HAND2', 'BMP4', 'ALX1',
                'ALX4', 'PRRX1', 'GREM1', 'LHX2', 'LHX9', 
                'TBX2', 'TBX4', 'TBX5', 'LMX1', 'MEIS1', 'MEIS2', 'SALL4', 'IRX3', 'IRX5')
limb.go = unique(c(limb.go, dev.example))

fpm$limb.go = NA
fpm$limb.go[!is.na(match(fpm$gene, limb.go))] = '1'

#fpm$genetype[!is.na(match(fpm$gene, limb.go))] = 'limb.dev.GO'

#plot(apply(fpm[which(fpm$genetype == 'hs'), c(1, 2)], 1, mean), 
#     (fpm[which(fpm$genetype == 'hs'), 1] - fpm[which(fpm$genetype == 'hs'), 2]), cex = 0.2)
#ss = seq(-10, 15, length.out = 100)
# points(ss, ss*0 + exp(ss)/(1- exp(ss)), lwd = 2.0, col = 'red', type = 'l')

fpm$mins = apply(fpm[, c(1, 3)], 1, min)
fpm$maxs = apply(fpm[, c(1, 3)], 1, max)
fpm = fpm[which(fpm$maxs > 1), ]

fpm$genetype[!is.na(fpm$DEgene)] = 'DE.gene'

examples.sel = unique(grep(paste0(dev.example, collapse = '|') ,fpm$gene))

ggplot(fpm, aes(x = Stage40, y = mUA, color = genetype, label = gene)) +
  geom_point(size = 0.2) + 
  scale_color_manual(values=c('red', "#E69F00", "#999999", 'darkblue', "#56B4E9")) + 
  #geom_point(data=fpm[which(!is.na(fpm$limb.go)), ], aes(x=Stage40, y=mUA), size=2) + 
  geom_text_repel(data= fpm[examples.sel, ], size = 4.0, color = 'darkblue') + 
  #geom_hline(yintercept=2.0, colour = "darkgray") + 
  #geom_vline(xintercept = 2.0, colour = "darkgray")
  geom_abline(slope = 1,  intercept = 0, colour = 'green') +
  theme_classic() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
        #legend.key.size = unit(1, 'cm')
        #legend.key.width= unit(1, 'cm')
  ) + guides(colour = guide_legend(override.aes = list(size=4)))

# ggsave(paste0(figureDir, "Stage40_vs_mUA_devGenes_comparisons_Gerber2018.pseudobulk.pdf"), width=12, height = 10)

fpm$DEgene[which(!is.na(fpm$DEgene) & fpm$log2fc>0)] = 'devGene'
fpm$DEgene[which(!is.na(fpm$DEgene) & fpm$log2fc<0)] = 'matureGene'


saveRDS(fpm, file = paste0(RdataDir, '/pseudoBulk_scRNAcellPooling_FluidigmC1_stage40.44.mUA_2500DEgenes.edgR.Rdata'))


##########################################
# similar comparison mUA vs stage samples for mRNAs 
##########################################
RNAseqDir = "../results/rnaseq_Rxxxx.old_R10724_R161513_mergedTechRep/Rdata/"
load(file = paste0(RNAseqDir, 'RNAseq_design_dds.object.Rdata'))

# select mature samples
sels = grep('Mature_UA|Embryo', design.matrix$condition)
dds = dds[, sels]

dds$condition = droplevels(dds$condition)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('condition', 'batch'), returnData = TRUE))
ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5)

plot(ggp)


dds$condition = droplevels(dds$condition)
dds <- estimateDispersions(dds, fitType = 'parametric')
plotDispEsts(dds, ymin = 10^-3); abline(h = 0.1, col = 'blue', lwd = 2.0)
# 
dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

# 
res.ii = results(dds, contrast=c("condition", 'Embryo_Stage40', 'Mature_UA'), alpha = 0.1)
colnames(res.ii) = paste0(colnames(res.ii), "_Hand.vs.UA")
res = data.frame(res.ii[, c(2, 5, 6)])
# 
res.ii = results(dds, contrast=c("condition", 'Mature_Hand', 'Mature_LA'), alpha = 0.1)
colnames(res.ii) = paste0(colnames(res.ii), "_Hand.vs.LA")
res = data.frame(res, res.ii[, c(2, 5, 6)])


cpm = fpm(dds)
cpm = log2(cpm + 2^-7)

xx = data.frame(mUA = apply(cpm[, grep('Mature_UA', colnames(cpm))], 1, median), 
                stage40 = apply(cpm[, grep('Embryo_Stage40', colnames(cpm))], 1, mean), 
                stage46_prox = apply(cpm[, grep('Embryo_Stage46_proximal', colnames(cpm))], 1, mean), 
                stage46_dist = apply(cpm[, grep('Embryo_Stage46_distal', colnames(cpm))], 1, mean))

maxs = apply(xx, 1, max)

cpm = xx
cpm$gene =  sapply(rownames(cpm), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[1])})

jj = grep(paste0(dev.example, collapse = '|') ,cpm$gene)
cpm = cpm[which(maxs> -1), ]

load(file =  paste0(annotDir, 'axolotl_housekeepingGenes_controls.other.tissues.liver.islet.testis_expressedIn21tissues.Rdata'))
hs = controls.tissue$geneIDs[which(controls.tissue$tissues == 'housekeeping')]
ctl =  controls.tissue$geneIDs[which(controls.tissue$tissues  != 'housekeeping')]
ggs = sapply(rownames(cpm), function(x) {x = unlist(strsplit(as.character(x), '_')); return(x[length(x)])})

cpm$genetype = NA

mm = match(ggs, hs)
cpm$genetype[!is.na(mm)] = 'hs'

mm = match(ggs, ctl)
cpm$genetype[!is.na(mm) & is.na(cpm$genetype)] = 'ctl'

par(mfrow = c(1, 3))
hist(cpm[, 2], breaks = 100, main = 'log2 expression of stage40');abline(v = 2)
hist(cpm[, 3], breaks = 100, main = 'log2 expression of stage44_prox');abline(v = 2)
hist(cpm[, 4], breaks = 100, main = 'log2 expression of stage44_dist');abline(v = 2)

cpm$genetype[which(cpm[,2]< -1 & cpm[, 3]< -1 & cpm[, 4] < -1 & is.na(cpm$genetype))] = 'dev.lowlyExp'
cpm$genetype[is.na(cpm$genetype)] = 'dev.Exp'

cpm$genetype[!is.na(match(cpm$gene, limb.go))] = 'limb.dev.GO'

examples.sel = unique(grep(paste0(dev.example, collapse = '|') ,cpm$gene))

cpm = data.frame(cpm)
cpm = cpm[which(cpm$genetype != 'ctl'), ]

ggplot(cpm, aes(x = stage40, y = mUA, color = genetype, label = gene)) +
  geom_point(size = 0.4) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9",  "#999999",  'darkblue')) + 
  geom_text_repel(data= cpm[examples.sel, ], size = 4.0, color = 'darkblue') +
  geom_point(data=cpm[which(cpm$genetype == 'limb.dev'), ], aes(x=stage40, y=mUA), colour="darkblue", size=1.5) +
  #geom_hline(yintercept=2.0, colour = "darkgray") + 
  #geom_vline(xintercept = 2.0, colour = "darkgray")
  geom_abline(slope = 1,  intercept = 0, colour = 'red') +
  theme_classic() +
  theme(legend.text = element_text(size=10), 
        #legend.title = element_text(color = "blue", size = 14),
        legend.key.size = unit(1, 'cm')
        #legend.key.width= unit(1, 'cm')
  )


ggsave(paste0(figureDir, "Stage40_vs_mUA_devGenes_comparisons_smartseq2.pdf"), width=12, height = 10)

##########################################
# ATAC-seq peaks: development-specific peaks and mature peaks 
##########################################
Compare.dev.vs.mature.ATACseq.peaks = FALSE
if(Compare.dev.vs.mature.ATACseq.peaks){
  
  require(ChIPpeakAnno)
  require(ChIPseeker)
  
  ##########################################
  # import batch corrected gene expression and design 
  ##########################################
  fpm = readRDS(file = paste0(RdataDir, '/fpm.bc_TMM_combat_Dev.mUA.samples_batch2020.2021.rds'))
  design = readRDS(file = paste0(RdataDir, '/design_sels_bc_TMM_combat_Dev.mUA.samples_batch2020.2021.rds'))
  
  # prepare the background distribution
  fpm.bg = fpm[grep('bg_', rownames(fpm), invert = FALSE), ]
  fpm = fpm[grep('bg_', rownames(fpm), invert = TRUE), ]
  rownames(fpm) = gsub('_', '-', rownames(fpm))
  
  hist(fpm.bg, breaks = 100, main = 'background distribution')
  abline(v = 3, col = 'red', lwd = 2.0)
  quantile(fpm.bg, c(0.95, 0.99))
  
  ##########################################
  ## make Granges and annotate peaks
  ##########################################
  Make.Granges.and.peakAnnotation = TRUE
  if(Make.Granges.and.peakAnnotation){
    pp = data.frame(t(sapply(rownames(fpm), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
    pp$strand = '*'
    
    save.peak.bed = FALSE
    if(save.peak.bed){
      bed = data.frame(pp[, c(1:3)], rownames(pp), 0, pp[, 4], stringsAsFactors = FALSE)
      write.table(bed, file = paste0(resDir, '/peakset_', version.analysis, '.bed'), col.names = FALSE, row.names = FALSE,
                  sep = '\t', quote = FALSE)
    }
    
    pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                  start.field="X2", end.field="X3", strand.field="strand")
    
    # annotation from ucsc browser ambMex60DD_genes_putative
    amex = GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
    pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
    #plotAnnoBar(pp.annots)
    
    pp.annots.all = pp.annots
    
    pp.annots = as.data.frame(pp.annots)
    rownames(pp.annots) = rownames(fpm)
    
    promoters = select.promoters.regions(upstream = 2000, downstream = 2000, ORF.type.gtf = 'Putative', promoter.select = 'all')
    
  }
  
  ##########################################
  # run spatial test for mature samples
  ##########################################
  Run.test.Dev.mUA.peaks = FALSE
  if(Run.test.spatial.peaks){
    conds = c("Mature_UA",  "Mature_Hand", "Embryo_Stage40", "Embryo_Stage44_proximal", "Embryo_Stage44_distal")
    
    sample.sels = c();  cc = c()
    for(n in 1:length(conds)) {
      #kk = which(design$conds == conds[n] & design$SampleID != '136159')
      kk = which(design$conds == conds[n]) 
      sample.sels = c(sample.sels, kk)
      cc = c(cc, rep(conds[n], length(kk)))
      
    }
    
    library(tictoc)
    ii.test = c(1:nrow(fpm)) # takes about 2 mins for 40k peaks
    source('Functions_atac.R')
    
    # select the mature samples
    cpm = fpm[, sample.sels]
    
    # select the peaks that are above background with >3 samples
    quantile(fpm.bg, 0.99)
    
    hist(fpm.bg)
    
    tic() 
    
    res = dev.mUA.peaks.test(cpm = cpm, c = cc, test.Dev.Reg = FALSE)
    
    #res = data.frame(res, pp.annots[ii.test, ], stringsAsFactors = FALSE)
    toc()
    
    xx = data.frame(res, pp.annots[match(rownames(res), rownames(pp.annots)), ], stringsAsFactors = FALSE)
    
    res = xx
    saveRDS(res, file = paste0(RdataDir, '/results_ATACseq.peaks_Dev.vs.mUA.test_', version.analysis, '_v2.rds'))
    
    rm(xx)
    
  }
  
  ##########################################
  # select all positional-dependent loci with below threshold
  ##########################################
  res = readRDS(file = paste0(RdataDir, '/results_ATACseq.peaks_Dev.vs.mUA.test_', version.analysis, '_v2.rds'))
  
  # select the positional peaks with pairwise comparisions 
  # limma logFC is in log2 scale
  fdr.cutoff = 0.05; logfc.cutoff = 1;
  
  jj = which((res$adj.P.Val_E40.vs.mUA < fdr.cutoff & abs(res$logFC_E40.vs.mUA) > logfc.cutoff & 
                res$adj.P.Val_E40.vs.mHand < fdr.cutoff & abs(res$logFC_E40.vs.mHand) > logfc.cutoff) |
               (res$adj.P.Val_E44prox.vs.mUA < fdr.cutoff & abs(res$logFC_E44prox.vs.mUA) > logfc.cutoff &
                  res$adj.P.Val_E44_dist.vs.mHand < fdr.cutoff & abs(res$logFC_E44_dist.vs.mHand) > logfc.cutoff))
  
  cat(length(jj), '\n')
  
  xx = res[c(jj), ]
  
  xx[grep('HOXA13|SHOX', xx$transcriptId), ]
  # fpm[which(rownames(fpm) == 'chr2p:873464923-873465440'), ]
  
  xx[grep('HOXA13|SHOX|MEIS', xx$transcriptId), ]
  
  cat(nrow(xx), ' peaks left\n')
  # sort positional peaks with logFC
  #xx = xx[order(-xx$logFC.mean), ]
  xx = xx[order(-abs(xx$logFC.mean)), ]
  
  ########
  ## asscociate the signifiant postional peaks with expression matrix
  ########
  conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", "Embryo_Stage44_distal", "Mature_UA",  "Mature_Hand")
  
  sample.sels = c();  cc = c()
  sample.means = c()
  for(n in 1:length(conds)) {
    kk = which(design$conds == conds[n]) 
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
    sample.means = cbind(sample.means, apply(fpm[, kk], 1, mean))
    
  }
  
  keep = fpm[(match(rownames(xx), rownames(fpm))), sample.sels]
  keep = as.matrix(keep)
  sample.means = sample.means[(match(rownames(xx), rownames(fpm))), ]
  
  ### filter the peaks with low signals 
  Filtering.peaks.in.Head.samples = TRUE
  if(Filtering.peaks.in.Head.samples){
    
    maxs = apply(sample.means, 1, max)
    mins = apply(sample.means, 1, min)
    hist(fpm.bg)
    
    sels = which(maxs > 3)
    cat(length(sels), 'peaks selected \n')
    
    sels = which(maxs > 3 & mins < 3)
    cat(length(sels), 'peaks selected \n')
    
    #sels = which(rr >1 & maxs > 3.)
    #cat(length(sels), 'peaks selected \n')
    #nonsels = which(rr<=1 | maxs <=2)
    xx = xx[sels, ]
    keep = keep[sels, ]
    
  }
  
  dim(xx)
  
  save(xx, keep, file = paste0(RdataDir, '/ATACseq_dev.mautre.peaks_filteredLowSignal_', version.analysis, '.Rdata'))
  
  ##########################################
  # quick clustering of postional peaks using 2 replicates of each condition
  ##########################################
  library(dendextend)
  library(ggplot2)
  
  load(file = paste0(RdataDir, '/ATACseq_dev.mautre.peaks_filteredLowSignal_', version.analysis, '.Rdata'))
  
  pp = data.frame(t(sapply(rownames(xx), function(x) unlist(strsplit(gsub('-', ':', as.character(x)), ':')))))
  pp$strand = '*'
  
  pp = makeGRangesFromDataFrame(pp, seqnames.field=c("X1"),
                                start.field="X2", end.field="X3", strand.field="strand")
  
  # make.pca.plots(keep, ntop = 1246, conds.plot = 'Mature')
  rep.sels = grep('1373|1361', colnames(keep), invert = TRUE)
  
  
  yy = keep[, rep.sels]
  
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
  yy <- t(apply(yy, 1, cal_z_score))
  
  nb_clusters = 6
  
  saveDir = paste0(figureDir, 'dev.vs.mature_peaks_clusters_', nb_clusters)
  if(!dir.exists(saveDir)) dir.create(saveDir)
  
  my_hclust_gene <- hclust(dist(yy), method = "complete")
  
  my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
  xx$clusters = my_gene_col
  
  my_gene_col <- data.frame(cluster =  paste0('cluster_', my_gene_col))
  rownames(my_gene_col) = rownames(yy)
  
  df <- data.frame(rep(conds, each = 2))
  rownames(df) = colnames(yy)
  colnames(df) = 'samples'
  
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
            "#33a02c", "#fb9a99", "#e31a1c",
            "#fdbf6f", "#ff7f00", "#cab2d6",
            "#6a3d9a", "#ffff99", "#b15928")
  
  sample_colors = c('green', 'green2', 'yellow2', 'springgreen4', 'gold2')
  names(sample_colors) = conds
  cluster_col = col3[1:nb_clusters]
  names(cluster_col) = paste0('cluster_', c(1:nb_clusters))
  annot_colors = list(
    segments = sample_colors,
    cluster = cluster_col)
  
  gaps.col = c(2, 4, 6, 8)
  
  pheatmap(yy, annotation_row = my_gene_col, 
           annotation_col = df, show_rownames = FALSE, scale = 'none', 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlGn")))(8), 
           show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = FALSE,  
           clustering_method = 'complete', cutree_rows = nb_clusters, 
           annotation_colors = annot_colors, 
           gaps_col = gaps.col, 
           #gaps_row =  gaps.row, 
           filename = paste0(saveDir, '/heatmap_dev.vs.mature_peaks_clusters_', nb_clusters, '.pdf'), 
           width = 6, height = 12)
  
  write.csv(xx, file = paste0(saveDir, '/dev.vs.mature_peaks_with.clusters', 
                              nb_clusters, '.csv'), quote = FALSE, row.names = TRUE)
  
  
  ## test the distribution of features of different groups
  source('Functions_atac.R')
  
  pdfname = paste0(saveDir, "/Features_distribution_allPeaks.pdf")
  pdf(pdfname, width = 8, height = 6)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  pp.annots = annotatePeak(pp, TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  
  stats.all = plotPeakAnnot_piechart(pp.annots.all)
  
  dev.off()
  
  pdfname = paste0(saveDir, "/Features_distribution_DBpeak_promoterDepeletion_10-274.pdf")
  pdf(pdfname, width = 8, height = 6)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  pp.annots = annotatePeak(pp[match(rownames(yy), names(pp))], TxDb=amex, tssRegion = c(-2000, 2000), level = 'transcript')
  
  stats = plotPeakAnnot_piechart(pp.annots)
  
  dev.off()
  
  
  total = pp.annots.all@peakNum
  mm = round(total * stats.all[1, 2]/100)
  nn = total - mm
  qq = round(pp.annots@peakNum * stats[1, 2]/100)
  phyper(qq-1, mm, nn, pp.annots@peakNum, lower.tail = FALSE, log.p = FALSE)
  phyper(qq+1, mm, nn, pp.annots@peakNum, lower.tail = TRUE, log.p = FALSE)
  
  
  ##########################################
  # double check the promoter accessibility of dev genes and mature genes
  # for dev-/mature-specific peaks at promotes, check if there are another peaks with higher and static signals
  ##########################################
  annot = readRDS(paste0(annotDir, 
                         'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
  
  
  load(file = paste0(RdataDir, '/ATACseq_dev.mautre.peaks_filteredLowSignal_', version.analysis, '.Rdata'))
  
  promoter.sels = grep('Promoter', xx$annotation)
  yy = keep[promoter.sels, rep.sels]
  peaks = xx[promoter.sels, ]
  
  peaks$transcript = sapply(peaks$transcriptId, function(x) {x = unlist(strsplit(as.character(x), '[|]')); return(x[length(x)])})
  
  peaks$Gene = annot$geneID[match(peaks$transcript, annot$transcriptID)]
  
  refs = res[grep('Promoter', res$annotation), ]
  refs$transcript = sapply(refs$transcriptId, function(x) {x = unlist(strsplit(as.character(x), '[|]')); return(x[length(x)])})
  refs$Gene = annot$geneID[match(refs$transcript, annot$transcriptID)]
  
  gg.uniq = unique(peaks$Gene)
  
  index1 = rep(NA, length(gg.uniq))
  index2 = rep(NA, length(gg.uniq))
  
  for(n in 1:length(gg.uniq))
  {
    # n = 1230
    cat(n, '\n')
    
    ii.dynp = match(rownames(peaks)[which(peaks$Gene == gg.uniq[n])], rownames(cpm))
    
    if(length(ii.dynp)== 1) {
      mean.dynp = mean(cpm[ii.dynp, ])
      index1[n] = ii.dynp
    }
    if(length(ii.dynp) >1){
      mean.dynp = apply(cpm[ii.dynp, ], 1, mean)
      index1[n] = ii.dynp[which(mean.dynp == max(mean.dynp))]
      mean.dynp = max(mean.dynp)
    }
    
    ii.all = match(rownames(refs)[which(refs$Gene == gg.uniq[n])], rownames(cpm))
    ii.altern = setdiff(ii.all, ii.dynp)
    
    if(length(ii.altern) ==1 ){
      mean.altern = mean(cpm[ii.altern, ])
      if(mean.altern > mean.dynp){
        index2[n] = ii.altern
      }
    }
    
    if(length(ii.altern) >1){
      mean.altern = apply(cpm[ii.altern, ], 1, mean)
      max.mean.altern = max(mean.altern)
      if(max.mean.altern > mean.dynp){
        index2[n] = ii.altern[which(mean.altern == max.mean.altern)]
      }
    }
    
  }
  
  indexs = data.frame(gg.uniq, index1, index2, stringsAsFactors = FALSE)
  indexs$gene = annot$gene.symbol.hs_gtf.gene[match(indexs$gg.uniq, annot$geneID)]
  
  rm(index2)
  rm(index1)
  
  saveRDS(indexs, file = paste0(RdataDir, '/Dev_vs_Mature_dyanmic_promoters_alterantiveStaticPromoter.rds'))
  
  # reorder the samples 
  conds = c("Embryo_Stage40", "Embryo_Stage44_proximal", "Embryo_Stage44_distal", "Mature_UA",  "Mature_Hand")
  
  sample.sels = c();  cc = c()
  
  for(n in 1:length(conds)) {
    kk = grep(conds[n], colnames(cpm))
    sample.sels = c(sample.sels, kk)
    cc = c(cc, rep(conds[n], length(kk)))
    
  }
  
  yy = cpm[ ,sample.sels]
  yy = yy[, grep('1373|1361', colnames(yy), invert = TRUE)]
  yy = data.frame(yy[as.numeric(indexs$index1), ], yy[as.numeric(indexs$index2), ])
  rownames(yy) = paste0(indexs$gene, "_", indexs$gg.uniq)
  
  cal_z_score <- function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  yy <- t(apply(yy, 1, cal_z_score))
  
  
  nb_clusters = 6
  saveDir = paste0(figureDir, 'dev.vs.mature_promoterPeaks_clusters_', nb_clusters)
  if(!dir.exists(saveDir)) dir.create(saveDir)
  
  my_hclust_gene <- hclust(dist(yy), method = "complete")
  
  my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = nb_clusters)
  indexs$clusters = my_gene_col
  
  saveRDS(indexs, file = paste0(RdataDir, '/Dev_vs_Mature_dyanmic_promoters_alterantiveStaticPromoter_withClusters.rds'))
  
  my_gene_col <- data.frame(cluster =  paste0('cluster_', my_gene_col))
  rownames(my_gene_col) = rownames(yy)
  
  df <- data.frame(rep(rep(conds, each = 2), 2))
  rownames(df) = colnames(yy)
  colnames(df) = 'samples'
  
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
            "#33a02c", "#fb9a99", "#e31a1c",
            "#fdbf6f", "#ff7f00", "#cab2d6",
            "#6a3d9a", "#ffff99", "#b15928")
  
  sample_colors = c('green', 'green2', 'yellow2', 'springgreen4', 'gold2')
  names(sample_colors) = conds
  
  annot_colors = list(
    segments = sample_colors,
    cluster = cluster_col)
  
  gaps.col = c(2, 4, 6, 8, 10, 12, 14, 16)
  
  pheatmap(yy,  
           annotation_row = my_gene_col, clustering_method = 'complete', cutree_rows = nb_clusters, # cluster parameters
           annotation_col = df, show_rownames = FALSE, scale = 'none', 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlGn")))(8), 
           show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = FALSE, 
           annotation_colors = annot_colors, 
           gaps_col = gaps.col, 
           #gaps_row =  gaps.row, 
           filename = paste0(saveDir, '/heatmap_dev.vs.mature_peaks_clusters_alternativePromoter_625genes.pdf'), 
           width = 10, height = 12)
  
  ## further check the dyanmic promoters related to dev genes
  dm = readRDS(file = paste0(RdataDir, '/Dev_vs_Mature_dyanmic_promoters_alterantiveStaticPromoter_withClusters.rds'))
  rownames(dm) = rownames(cpm)[dm$index1]
  
  dm = dm[which(is.na(dm$index2 & dm$clusters == 2)), ]
  
  
  gg = xx$geneId
  grep('HOXA13', gg)
  rownames(yy) = paste0(rownames(yy), '_', gg)
  #rownames(keep) = gg
  yy = as.matrix(yy)
  
  gg = rownames(yy)
  gg = sapply(gg, function(x) {x = unlist(strsplit(as.character(x), '_')); return(paste0(x[2:length(x)], collapse = '_'))})
  gg = sapply(gg, function(x) {
    xx = unlist(strsplit(as.character(x), '[|]')); 
    if(xx[1] == 'N/A') 
    {
      return(x);
    }else{
      xx = xx[-length(xx)]
      xx = xx[length(xx)]
      return(xx);
    }})
  gg = as.character(gg)
  gg = gsub("\\[|\\]", "", gg)
  gg = gsub(' hs', '', gg)
  gg = gsub(' nr', '', gg)
  
  rownames(yy) = gg
  sample_colors = c('springgreen4', 'steelblue2', 'gold2')
  names(sample_colors) = c('Mature_UA', 'Mature_LA', 'Mature_Hand')
  annot_colors = list(segments = sample_colors)
  
  pheatmap(yy, 
           annotation_col = df, show_rownames = TRUE, scale = 'row', 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlGn")))(8), 
           show_colnames = FALSE,
           cluster_rows = TRUE, cluster_cols = FALSE, 
           annotation_colors = annot_colors, 
           gaps_col = gaps.col, 
           #gaps_row =  gaps.row, 
           filename = paste0(figureDir, '/heatmap_positionalPeaks_fdr0.01_log2FC.1_top.promoters.pdf'), 
           width = 10, height = 8)
  
  if(saveTable){
    write.csv(data.frame(keep, yy, stringsAsFactors = FALSE), 
              file = paste0(resDir, '/position_dependent_peaks_from_matureSamples_ATACseq_rmPeaks.head_top50_promoterPeaks.csv'), 
              quote = FALSE, row.names = TRUE)
    
  }
  
}



