##########################################################################
##########################################################################
# Project: positional memory project
# Script purpose: process TSS and merge with regeneration TSS
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon May 16 10:38:58 2022
##########################################################################
##########################################################################

########################################################
########################################################
# Section : first process TSS from mature samples 
# batch correction
########################################################
########################################################
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
