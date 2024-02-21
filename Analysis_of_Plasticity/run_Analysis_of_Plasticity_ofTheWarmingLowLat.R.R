# Analysis of Plasticity of the evolved populations
#necessary
library(lme4)
library(car)
library(glmmTMB)
library(bbmle) #necessary to AIC

#check if there is no space in the begining of the header
normalizedCounts <- read.table("Galaxy359-Normalized_counts.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

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

#Creating new data frame - name: NewColumns.Df,   
#Contains four columns History, Environment, AP and Selection
NewColumns.Df <- data.frame(Environment)


data1 <- transposed_of_normalizedCounts               # Replicate example data
data2 <- NewColumns.Df


#iterate over each gene in file
for(i in 1:ncol(data1)) {       # for-loop over columns
  geneName <- colnames(data1)[i]      # gene name
  data3 <- cbind(data1[ , i],data2)
  meanPTW <- (data3[4,1]+data3[5,1]+data3[6,1])/3    #PT_TREATMENT_counts are in lines 4,5 and 6 of the first column
  meanPTC <- (data3[1,1]+data3[2,1]+data3[3,1])/3       #PT_CONTROL_counts are in lines 1,2 and 3 of the first column
  upOrDownPT<-meanPTW-meanPTC
  Expression <- data3[,1]   #expression is always in the first column, this way we can vary the number of factors without needing to alter this line of code
  options(contrasts=c("contr.sum","contr.poly"))
  glmmTMB1_binN = try(glmmTMB(Expression ~ Environment+(1|AP),data=data3,family="nbinom2"))
  AnovaGlmmTMB1_binN = try(as.data.frame(Anova(glmmTMB1_binN, type='III')))
  if (class(AnovaGlmmTMB1_binN) == "try-error") {
    Chisq <- c("ERROR", "ERROR")
    Df <- c("ERROR", "ERROR")
    pvalue <- c("ERROR", "ERROR")
    
    AnovaGlmmTMB1_binN <- data.frame(Chisq, Df, pvalue)
  }
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownWPT)

  
  # changing row names of data frame
  rownames(AnovaGlmmTMB1_binN)[3] <- "upOrDownPT"
  #rownames(AnovaGlmmTMB1_binN)[6] <- "upOrDownWPT"

  
  #Sort Using Character row.names
  #subset
  subsetted <- as.data.frame(subset(AnovaGlmmTMB1_binN, select = -c(`Chisq`,Df)))
  #change column name
  colnames(subsetted) <- c(geneName)
  #transpose
  transposed_of_AnovaGlmmTMB1_binN <- as.data.frame(t(subsetted),stringsAsFactors = FALSE)
  try(write.table(transposed_of_AnovaGlmmTMB1_binN, "PT_onlyPlasticity.csv" , append = TRUE, sep = '\t', col.names = FALSE, row.names = TRUE))

}

