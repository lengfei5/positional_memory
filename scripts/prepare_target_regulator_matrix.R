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


### motif-TF-association matrix
mapping =readRDS(file = '../data/JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata_manual_rmRedundantUNVALIDED.rds')
mm = match(colnames(target.moc), mapping$name)
target.moc = target.moc[ ,which(!is.na(mm))] # filter the unvalided motifs if core motifs are present
#mat_gcre = mat_gcre[, which(!is.na(match(colnames(mat_gcre), mapping$name)))]

#atm = table(mapping$name, mapping$gene)
atm = matrix(0, nrow = ncol(target.moc), ncol = nrow(target.moc))
colnames(atm) = rownames(target.moc)
rownames(atm) = colnames(target.moc)
genes = get_geneName(colnames(atm))

for(n in 1:nrow(atm))
{
  # n = 164
  tfs = unlist(strsplit(as.character(mapping$tfs[which(mapping$name == rownames(atm)[n])]), '_'))
  for(tf in tfs)
  {
    jj = which(genes == tf)
    if(length(jj)>0){
      cat(n, '--', rownames(atm)[n], ' -- tfs --', tf, '\n')
      atm[n, jj] = 1
    }
  }
}

target.tfs = target.moc %*% atm


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


saveRDS(target.tfs, file = paste0(RdataDir, '/GRN_priorNetwork_target_TFs.rds'))



