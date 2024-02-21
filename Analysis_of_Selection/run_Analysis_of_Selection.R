# Analysis of Selection

#necessary
library(lme4)
library(car)
library(glmmTMB)
library(bbmle) #necessary to AIC

#check if there is no space in the beginning of the header
#file "Galaxy218-Normalized_counts.tabular" was used for the analysis of southern populations
#file "Galaxy225-Normalized_counts.tabular" was used for the analysis of southern populations
normalizedCounts <- read.table("Galaxy218-Normalized_counts.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]

#transpose
transposed_of_normalizedCounts <- t(dataInAtLeastXsamples)    #transpose the rows and columns

#read samples
rownames <- rownames(transposed_of_normalizedCounts)
splited <- strsplit(rownames, split = "_")
AP <-lapply(splited, `[[`, 2)  #AP is the second element of list
AP<-as.character(AP) #transform list in character
Selection <-lapply(splited, `[[`, 1) #selection is the first element of list
Selection<-as.character(Selection) #transform list in character

#Creating new data frame - name: NewColumns.Df,   
#Contains one column History, Environment, AP and Selection
NewColumns.Df <- data.frame(AP,Selection)


data1 <- transposed_of_normalizedCounts               # Replicate example data
data2 <- NewColumns.Df


#iterate over each gene in file
for(i in 1:ncol(data1)) {       # for-loop over columns
  geneName <- colnames(data1)[i]      # gene name
  data3 <- cbind(data1[ , i],data2)
  meanTreatment <- (data3[4,1]+data3[5,1]+data3[6,1])/3    #TREATMENT_counts for selection are in line 4,5 and 6 of the first column
  meanControl <- (data3[1,1]+data3[2,1]+data3[3,1])/3       #CONTROL_counts for selection are in lines 1,2 and 3 of the first column
  upOrDown<-meanTreatment-meanControl
  Expression <- data3[,1]   #expression is always in the first column, this way we can vary the number of factors without needing to alter this line of code
  options(contrasts=c("contr.sum","contr.poly"))
  glmmTMB1_binN = try(glmmTMB(Expression ~ Selection+(1|AP),data=data3,family="nbinom2"))
  AnovaGlmmTMB1_binN = try(as.data.frame(Anova(glmmTMB1_binN, type='III')))
  if (class(AnovaGlmmTMB1_binN) == "try-error") {
    Chisq <- c("ERROR", "ERROR")
    Df <- c("ERROR", "ERROR")
    pvalue <- c("ERROR", "ERROR")
    
    AnovaGlmmTMB1_binN <- data.frame(Chisq, Df, pvalue)
  }
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDown)

  
  # changing row names of data frame
  rownames(AnovaGlmmTMB1_binN)[3] <- "upOrDown"

  
  #Sort Using Character row.names
  #subset
  subsetted <- as.data.frame(subset(AnovaGlmmTMB1_binN, select = -c(`Chisq`,Df)))
  #change column name
  colnames(subsetted) <- c(geneName)
  #transpose
  transposed_of_AnovaGlmmTMB1_binN <- as.data.frame(t(subsetted),stringsAsFactors = FALSE)
  #for high latitude populations change the name of the output file
  try(write.table(transposed_of_AnovaGlmmTMB1_binN, "PT_only_selection.csv" , append = TRUE, sep = '\t', col.names = FALSE, row.names = TRUE))
 
}

