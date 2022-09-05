##########################################################################
##########################################################################
# Project: positional memory 
# Script purpose: prepare the genome-wide target-to-regulator (TFs) matrix for GRN inference
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Aug 29 15:06:40 2022
##########################################################################
##########################################################################
E_all = readRDS(file = paste0(RdataDir, '/GRNinference_scExprMatrix_allGenes_290TFsymbols.rds'))

source("myGENIE3.R")

###  construct target-CRE matrix
target.ids_all = get_geneID(rownames(E_all))
mat_gcre_all = build_target_CRE_matrix(target.ids_all)
mm = match(rownames(mat_gcre_all), target.ids_all)
rownames(mat_gcre_all) = rownames(E_all)[mm]


### CRE-motif-occurence matrix (regulaory regions * motif occurrence)
mocs = readRDS(file = '../results/motif_analysis/motif_oc_fimo_atacPeaks.2kbTSS_jaspar2022.core.unvalided_pval.0.00001_v1.rds')
cres = unique(rownames(mocs))

## match the tf-CRE for all targets
jj = match(colnames(mat_gcre_all), rownames(mocs))
missed = which(is.na(jj))
cat(length(missed), ' CRE missed \n')

ii = which(!is.na(jj))
jj = jj[ii]

mat_gcre_sel = mat_gcre_all[,ii]
mocs_sel = mocs[jj, ]

rm(E_all); rm(mat_gcre_all); rm(mocs); rm(cres)

target.moc = mat_gcre_sel %*% mocs_sel 

ss = apply(target.moc, 2, sum)
length(which(ss==0))

saveRDS(target.moc, file = paste0(RdataDir, '/target_motifOccurrence_all.rds'))

### motif-TF-association matrix
mapping =readRDS(file = '../data/JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata_manual_rmRedundantUNVALIDED.rds')
mm = match(colnames(target.moc), mapping$name)

#target.moc_sel = target.moc[ ,which(!is.na(mm))] # filter the unvalided motifs if core motifs are present
#mat_gcre = mat_gcre[, which(!is.na(match(colnames(mat_gcre), mapping$name)))]

mapping$gene = as.character(mapping$gene)
atm = table(mapping$name, mapping$gene)

jj = match(colnames(target.moc), rownames(atm))
ii = which(!is.na(jj))

target.moc_sel = target.moc[, ii];
atm_sel = atm[jj[ii], ]

rm(target.moc);
target.tfs = target.moc_sel %*% atm_sel


x = target.tfs[grep('HOXA13', rownames(target.tfs)), ];
x[which(x>0)]

# check if some targets have no regulators
ss1 = apply(target.tfs, 1, sum)
cat(length(which(ss1==0)), 'targets with no regulators \n')
target.tfs = target.tfs[which(ss1>0), ]

# check if some regulators don't have any targets due to missing motifs 
ss2 = apply(target.tfs, 2, sum)
cat(length(which(ss2==0)), 'regulators with no targets found \n')

target.tfs = target.tfs[, which(ss2>0)]

cat(nrow(target.tfs), ' targets by ', ncol(target.tfs), ' TFs \n')

##########################################
# filter TFs with gene expression 
##########################################
E = readRDS(file = paste0(RdataDir, '/GRNinference_scExprMatrix_346geneID_289TFsymbols.rds'))
ggs = as.character(get_geneName(rownames(E)))
mm = match(colnames(target.tfs), ggs)

target.tfs_sel = target.tfs[, which(!is.na(mm))]

# check if some targets have no regulators
ss1 = apply(target.tfs_sel, 1, sum)
cat(length(which(ss1==0)), 'targets with no regulators \n')
target.tfs_sel = target.tfs_sel[which(ss1>0), ]

# check if some regulators don't have any targets due to missing motifs 
ss2 = apply(target.tfs_sel, 2, sum)
cat(length(which(ss2==0)), 'regulators with no targets found \n')

target.tfs_sel = target.tfs_sel[, which(ss2>0)]

cat(nrow(target.tfs_sel), ' targets by ', ncol(target.tfs_sel), ' TFs \n')

saveRDS(target.tfs_sel, file = paste0(RdataDir, '/GRN_priorNetwork_allTarget_expressedTFs.rds'))

ttf = readRDS(file = paste0(RdataDir, '/GRN_priorNetwork_allTarget_expressedTFs.rds'))

