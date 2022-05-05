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






