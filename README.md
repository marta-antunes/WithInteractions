# DiffExpWithInteractions
Differential gene expression analysis based on the negative binomial distribution with interaction terms.

DESeq2 allows to estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution but does not give p-values for interaction terms.

I wrote DiffExpWithInteractions to fill this gap.

input: file with normalized count data

usage: glmer_expressionData.R filename
