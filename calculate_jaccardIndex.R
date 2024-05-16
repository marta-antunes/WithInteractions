#read file PT
normalizedCounts <- read.table("Galaxy218-Normalized_counts.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]

#transpose
transposed_of_normalizedCounts <- t(dataInAtLeastXsamples)    #transpose the rows and columns

########## PT1 vs PT2 ##########

#subset to make 2 by 2 comparisons
PT1_PT2=transposed_of_normalizedCounts[c(1, 2), ]

#sums columns
numerator <- sum(apply(PT1_PT2, 2, min)) #2 indicates columns
denominator <- sum(apply(PT1_PT2, 2, max)) #2 indicates columns

JaccardIndex_PT1_PT2=numerator/denominator
JaccardIndex_PT1_PT2

########## PT2 vs PT3 ##########

#subset to make 2 by 2 comparisons
PT2_PT3=transposed_of_normalizedCounts[c(2, 3), ]

#sums columns
numerator <- sum(apply(PT2_PT3, 2, min)) #2 indicates columns
denominator <- sum(apply(PT2_PT3, 2, max)) #2 indicates columns

JaccardIndex_PT2_PT3=numerator/denominator
JaccardIndex_PT2_PT3

########## PT1 vs PT3 ##########

#subset to make 2 by 2 comparisons
PT1_PT3=transposed_of_normalizedCounts[c(1, 3), ]

#sums columns
numerator <- sum(apply(PT1_PT3, 2, min)) #2 indicates columns
denominator <- sum(apply(PT1_PT3, 2, max)) #2 indicates columns

JaccardIndex_PT1_PT3=numerator/denominator
JaccardIndex_PT1_PT3


############################################
#read file NL
############################################
normalizedCounts <- read.table("Galaxy225-Normalized_counts.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]

#transpose
transposed_of_normalizedCounts <- t(dataInAtLeastXsamples)    #transpose the rows and columns

########## NL1 vs NL2 ##########

#subset to make 2 by 2 comparisons
NL1_NL2=transposed_of_normalizedCounts[c(1, 2), ]

#sums columns
numerator <- sum(apply(NL1_NL2, 2, min)) #2 indicates columns
denominator <- sum(apply(NL1_NL2, 2, max)) #2 indicates columns

JaccardIndex_NL1_NL2=numerator/denominator
JaccardIndex_NL1_NL2

########## NL2 vs NL3 ##########

#subset to make 2 by 2 comparisons
NL2_NL3=transposed_of_normalizedCounts[c(2, 3), ]

#sums columns
numerator <- sum(apply(NL2_NL3, 2, min)) #2 indicates columns
denominator <- sum(apply(NL2_NL3, 2, max)) #2 indicates columns

JaccardIndex_NL2_NL3=numerator/denominator
JaccardIndex_NL2_NL3

########## NL1 vs NL3 ##########

#subset to make 2 by 2 comparisons
NL1_NL3=transposed_of_normalizedCounts[c(1, 3), ]

#sums columns
numerator <- sum(apply(NL1_NL3, 2, min)) #2 indicates columns
denominator <- sum(apply(NL1_NL3, 2, max)) #2 indicates columns

JaccardIndex_NL1_NL3=numerator/denominator
JaccardIndex_NL1_NL3


#####################
#read file PT again and then keep genes significant in NL that are not in PT (3974 genes)
###################
normalizedCounts <- read.table("Galaxy218-Normalized_counts.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]
row_numbers <- read.csv("interest.txt",sep=',',header=FALSE)

# Convert row_numbers to a vector
row_numbers <- as.vector(unlist(row_numbers))


#keep only 3974 genes
filtered_df <- dataInAtLeastXsamples[row_numbers, ]

write.table(filtered_df, file = "mtcars.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#transpose
transposed_of_normalizedCounts <- t(filtered_df)    #transpose the rows and columns

########## PT1 vs PT2 ##########

#subset to make 2 by 2 comparisons
PT1_PT2=transposed_of_normalizedCounts[c(1, 2), ]

#sums columns
numerator <- sum(apply(PT1_PT2, 2, min)) #2 indicates columns
denominator <- sum(apply(PT1_PT2, 2, max)) #2 indicates columns

JaccardIndex_PT1_PT2=numerator/denominator
JaccardIndex_PT1_PT2

########## PT2 vs PT3 ##########

#subset to make 2 by 2 comparisons
PT2_PT3=transposed_of_normalizedCounts[c(2, 3), ]

#sums columns
numerator <- sum(apply(PT2_PT3, 2, min)) #2 indicates columns
denominator <- sum(apply(PT2_PT3, 2, max)) #2 indicates columns

JaccardIndex_PT2_PT3=numerator/denominator
JaccardIndex_PT2_PT3

########## PT1 vs PT3 ##########

#subset to make 2 by 2 comparisons
PT1_PT3=transposed_of_normalizedCounts[c(1, 3), ]

#sums columns
numerator <- sum(apply(PT1_PT3, 2, min)) #2 indicates columns
denominator <- sum(apply(PT1_PT3, 2, max)) #2 indicates columns

JaccardIndex_PT1_PT3=numerator/denominator
JaccardIndex_PT1_PT3


############################################
#read file NL
############################################
normalizedCounts <- read.table("Galaxy225-Normalized_counts.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]
row_numbers <- read.csv("interest.txt",sep=',',header=FALSE)

# Convert row_numbers to a vector
row_numbers <- as.vector(unlist(row_numbers))


#keep only 3974 genes
filtered_df <- dataInAtLeastXsamples[row_numbers, ]


#transpose
transposed_of_normalizedCounts <- t(filtered_df)    #transpose the rows and columns
########## NL1 vs NL2 ##########

#subset to make 2 by 2 comparisons
NL1_NL2=transposed_of_normalizedCounts[c(1, 2), ]

#sums columns
numerator <- sum(apply(NL1_NL2, 2, min)) #2 indicates columns
denominator <- sum(apply(NL1_NL2, 2, max)) #2 indicates columns

JaccardIndex_NL1_NL2=numerator/denominator
JaccardIndex_NL1_NL2

########## NL2 vs NL3 ##########

#subset to make 2 by 2 comparisons
NL2_NL3=transposed_of_normalizedCounts[c(2, 3), ]

#sums columns
numerator <- sum(apply(NL2_NL3, 2, min)) #2 indicates columns
denominator <- sum(apply(NL2_NL3, 2, max)) #2 indicates columns

JaccardIndex_NL2_NL3=numerator/denominator
JaccardIndex_NL2_NL3

########## NL1 vs NL3 ##########

#subset to make 2 by 2 comparisons
NL1_NL3=transposed_of_normalizedCounts[c(1, 3), ]

#sums columns
numerator <- sum(apply(NL1_NL3, 2, min)) #2 indicates columns
denominator <- sum(apply(NL1_NL3, 2, max)) #2 indicates columns

JaccardIndex_NL1_NL3=numerator/denominator
JaccardIndex_NL1_NL3

JI<-c(JaccardIndex_NL1_NL2,JaccardIndex_NL2_NL3,JaccardIndex_NL1_NL3)
mean(JI)


