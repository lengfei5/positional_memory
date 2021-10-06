Position memory project with Akane

# ATAC-seq data analysis:
## important files:

   /groups/tanaka/People/current/jiwang/projects/positional_memory/Data/atacseq_using
   technical replicates merged bam files:
   /groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11637_atac/bam_merged_uniq_rmdup

   called peaks with macs2 for all replicates:
   /groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11637_atac/calledPeaks

   feature counts of pooled above peak:
   /groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11637_atac/featurecounts_peaks.Q30

   
   bigwigs files using DESeq2 scaling factors using above peaks
   /groups/tanaka/People/current/jiwang/projects/positional_memory/Data/R11637_atac/bigwigs_deeptools.scalingFactor

## consensus peaks were manually idenitified by comparing the biological replicates and taking the intersect of peaks (p<0.001) from two replicates or two out of three replicates; and then the the union of consensus peaks were used to quantify the counts and to compute the scaling factors with DESeq2

# RNA-seq data
 /groups/tanaka/People/current/jiwang/projects/positional_memory/Data/rnaseq_using
 used raw data were collected in ngs_raw, trashed raw data were not saved here

# Cut&Tag ChIP-seq

