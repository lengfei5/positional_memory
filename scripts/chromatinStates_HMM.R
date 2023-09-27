##########################################################################
##########################################################################
# Project: positional memory 
# Script purpose: infer chromatin states with chromatin data 
# the analysis is inspired by https://www.biorxiv.org/content/10.1101/2020.11.20.391524v1.supplementary-material
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu May  5 11:22:31 2022
##########################################################################
##########################################################################
## test code for depmixS4, original code from https://stackoverflow.com/questions/66682218/multinomial-hmm-using-r
library(depmixS4)
set.seed(1)
df <- data.frame(id = rep(1, each = 40), 
                year = seq(1961,2000),
                X1 = rbinom(40, size = 1, prob = 1 - 0.6) * rpois(40, lambda = 4),
                X2 = rbinom(40, size = 1, prob = 1 - 0.7) * rpois(40, lambda = 4),
                X3 = rbinom(40, size = 1, prob = 1 - 0.6) * rpois(40, lambda = 5),
                X4 = rbinom(40, size = 1, prob = 1 - 0.7) * rpois(40, lambda = 6))
# matrix for single multinomial response variable
X <- as.matrix(df[,c("X1", "X2", "X3", "X4")])
 
# formulate model
mod<-depmix(X ~ 1, data=df, nstates=3,
           family=multinomial("identity"))

# fit model
fmod <- fit(mod)

# show results
summary(fmod)
posterior(fmod)

# predict the states by estimating the posterior
est.states <- posterior(fmod)
head(est.states)


########################################################
########################################################
# Section : prepare the input matrix 
# 
########################################################
########################################################
tss = readRDS(file = paste0(RdataDir, '/regeneration_tss_perGene_smartseq2_atac_histM.rds'))

histMs = c('atac', 'H3K4me3','H3K27me3', 'H3K4me1', 'H3K27ac')
newcc = c('mUA', '5dpa', '9dpa', '13dpa.p', '13dpa.d')
conds = c("mUA", "BL5days", "BL9days", 'BL13days.prox', 'BL13days.dist')

kk = grep('atac_mUA|atac_X|H3K4me3_mUA|H3K4me3_BL|H3K27me3_mUA|H3K27me3_BL|H3K4me1_mUA|H3K4me1_BL|H3K27ac_mUA|H3K27ac_BL', colnames(tss))
df = tss[, kk]
colnames(df)[1:5] = paste0('atac_', newcc)

newdf = matrix(NA, nrow = nrow(df), ncol = length(newcc) *5)
newedf = data.frame(newdf)

for(n in 1:length(histMs))
{
  index = seq(((n-1)*5+1), n*5, by = 1)
  cat(n, '--', histMs[n], ' --', index, '\n')
  if(n ==1){
    for(m in 1:length(index)) 
    {
      newdf[, index[m]] = df[, m]
      #colnames(newdf)[index[m]] = paste0(histMs[n], '_', newcc[m])
    }
  }else{
    test = df[, grep(histMs[n], colnames(df))]
    test = cal_sample_means(cpm = test, conds = c('mUA', 'BL5days', 'BL9days', 'BL13days.prox', 'BL13days.dist'))
    for(m in 1:length(index)) 
    {
      newdf[, index[m]] = test[, m]
      #colnames(newdf)[index[m]] = paste0(histMs[n], '_', newcc[m])
    }
    
    #newdf[, index] = test
  }
}

colnames(newdf) = paste0(rep(histMs, each = length(newcc)), '_', newcc)
rownames(newdf) = rownames(df)

df = matrix(NA, ncol = length(histMs), nrow = nrow(newdf)*length(conds))

for(n in 1:length(histMs))
{
  # n = 1
  test = newdf[, grep(histMs[n], colnames(newdf))]
  df[,n] = test
}
colnames(df) = histMs
rownames(df) = paste0(rep(newcc, each = nrow(newdf)), '_', rownames(newdf))

saveRDS(df, file = paste0(RdataDir, '/chromatin_features_inputMatrix_forHMM.rds'))

##########################################
# test the speed of HMM
##########################################
df = readRDS(file = paste0(RdataDir, '/chromatin_features_inputMatrix_forHMM.rds'))
df = data.frame(2^df)
ss = apply(df, 1, sum)
df = df[which((ss>0.1)), ]
#df[is.na(df)] = 0.01
# df = df + 0.01
df = df[sample(1:nrow(df), 10000), ]

X <- as.matrix(df[,c("atac", "H3K4me3", "H3K27me3", "H3K4me1", "H3K27ac")])

require(tictoc)
tic()
# formulate model
#set.seed(100)
nb.states = 5
mod<-depmix(X ~ 1, data=df, nstates=nb.states, 
            #instart =  rep(1/nb.states, nb.states),
            trstart = runif(nb.states^2),
            family=multinomial("identity"))
# fit model
set.seed(1)
fmod <- fit(mod, emcontrol=em.control(classification="soft", maxit = 5000, tol = 1e-10, rand=FALSE))

toc()


# show results
summary(fmod)
posterior(fmod)

# predict the states by estimating the posterior
est.states <- posterior(fmod)
head(est.states)
