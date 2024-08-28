# Analysis of ...

#necessary
library(lme4)
library(car)
library(glmmTMB)
library(bbmle) #necessary to AIC

#check if there is no space in the begining of the header
normalizedCounts <- read.table("Galaxy179-Normalized_counts.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]

#transpose
transposed_of_normalizedCounts <- t(dataInAtLeastXsamples)    #transpose the rows and columns

#read samples
rownames <- rownames(transposed_of_normalizedCounts)
splited <- strsplit(rownames, split = "_")
Environment<-lapply(splited, `[[`, 4) #Environment is the fourth element of list
Environment<-as.character(Environment) #transform list in character
AP <-lapply(splited, `[[`, 2)  #AP is the second element of list
AP<-as.character(AP) #transform list in character
History <- c(rep("PT",3), rep("NL",6), rep("PT",3))


#Creating new data frame - name: NewColumns.Df,   
#Contains four columns History, Environment, AP and Selection
NewColumns.Df <- data.frame(Environment, History)


data1 <- transposed_of_normalizedCounts               # Replicate example data
data2 <- NewColumns.Df


#iterate over each gene in file
for(i in 1:ncol(data1)) {       # for-loop over columns
  geneName <- colnames(data1)[i]      # gene name
  data3 <- cbind(data1[ , i],data2)
  
  # we will consider up-regulated gene when is more expressed in NL than in PT
  meanNLEnvC <- (data3[3,1]+data3[4,1]+data3[5,1])/3    #TREATMENT_counts are in lines 3,4 and 5 of the first column
  meanPTEnvC <- (data3[1,1]+data3[2,1]+data3[3,1])/3       #CONTROL_counts are in lines 1,2 and 3 of the first column
  upOrDownEnvC<-meanNLEnvC-meanPTEnvC
  meanNLEnvW <- (data3[7,1]+data3[8,1]+data3[9,1])/3    #TREATMENT_counts are in lines 7,8 and 9 of the first column
  meanPTEnvW <- (data3[10,1]+data3[11,1]+data3[12,1])/3       #CONTROL_counts are in lines 10,11 and 12 of the first column
  upOrDownEnvW<-meanNLEnvW-meanPTEnvW
  Expression <- data3[,1]   #expression is always in the first column, this way we can vary the number of factors without needing to alter this line of code
  options(contrasts=c("contr.sum","contr.poly"))
  glmmTMB1_binN = try(glmmTMB(Expression ~ Environment*History+(1|AP),data=data3,family="nbinom2"))
  AnovaGlmmTMB1_binN = try(as.data.frame(Anova(glmmTMB1_binN, type='III')))
  if (class(AnovaGlmmTMB1_binN) == "try-error") {
    Chisq <- c("ERROR", "ERROR", "ERROR","ERROR")
    Df <- c("ERROR", "ERROR","ERROR","ERROR")
    pvalue <- c("ERROR", "ERROR","ERROR","ERROR")
    
    AnovaGlmmTMB1_binN <- data.frame(Chisq, Df, pvalue)
  }
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownEnvC)
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownEnvW)

  
  # changing row names of data frame
  rownames(AnovaGlmmTMB1_binN)[5] <- "upOrDownEnvC"
  rownames(AnovaGlmmTMB1_binN)[6] <- "upOrDownEnvW"

  
  #Sort Using Character row.names
  #subset
  subsetted <- as.data.frame(subset(AnovaGlmmTMB1_binN, select = -c(`Chisq`,Df)))
  #change column name
  colnames(subsetted) <- c(geneName)
  #transpose
  transposed_of_AnovaGlmmTMB1_binN <- as.data.frame(t(subsetted),stringsAsFactors = FALSE)
  try(write.table(transposed_of_AnovaGlmmTMB1_binN, "results.csv" , append = TRUE, sep = '\t', col.names = FALSE, row.names = TRUE))

}

