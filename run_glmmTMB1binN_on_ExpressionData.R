# G23 dataset
#install lme4
#install.packages("lme4",
#                 repos=c("http://lme4.r-forge.r-project.org/repos",
#                         getOption("repos")[["CRAN"]]))

#necessary libraries
library(lme4)
library(car)
library(glmmTMB)
library(bbmle) #necessary to AIC


args <- commandArgs(TRUE)
##First remove space in the beginning of the file
normalizedCounts <- read.table(args[1], sep = '\t', header=TRUE, stringsAsFactors = TRUE)
#normalizedCounts <- read.table("Galaxy238.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)


#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]


#transpose
transposed_of_normalizedCounts <- t(dataInAtLeastXsamples)    #transpose the rows and columns


#Creating new data frame - name: NewColumns.Df,   
#Contains four columns: History, Environment, AP and Selection
NewColumns.Df <- data.frame(  
  History = c(rep("PT",3), rep("NL",6), rep("PT",3), rep("PT",6), rep("NL",6)),
  Environment = c(rep("C",4),rep("W",1),rep("C",1), rep("W",1), rep("C",1),rep("W",4),rep("C",3), rep("W",3),rep("C",1), rep("W",1),rep("C",1), rep("W",1),rep("C",1), rep("W",1)),
  AP = c("PT1", "PT2", "PT3", "NL1", "NL1", "NL2", "NL2", "NL3", "NL3", "PT1", "PT2", "PT3", "PT1", "PT2", "PT3", "PT1", "PT2", "PT3", "NL1", "NL1", "NL2","NL2", "NL3", "NL3"),
  Selection = c(rep("C",3), rep("C",6), rep("C",3), rep("W",6), rep("W",6)))


data1 <- transposed_of_normalizedCounts               # Replicate example data
data2 <- NewColumns.Df


#iterate over each gene in file
for(i in 1:ncol(data1)) {       # for-loop over columns
  geneName <- colnames(data1)[i]      # gene name
  data3 <- cbind(data1[ , i],data2)
  meanPTW <- (data3[10,1]+data3[11,1]+data3[12,1])/3    #PT_TREATMENT_counts are in lines 10,11 and 12 of the first column
  meanPTC <- (data3[1,1]+data3[2,1]+data3[3,1])/3       #PT_CONTROL_counts are in lines 1,2 and 3 of the first column
  upOrDownPT<-meanPTW-meanPTC
  #print(paste0("PT Up or Down? ", upOrDownPT))
  meanNLW <- (data3[5,1]+data3[7,1]+data3[9,1])/3    #NL_TREATMENT_counts are in lines 5,7 and 9 of the first column
  meanNLC <- (data3[4,1]+data3[6,1]+data3[8,1])/3    #NL_CONTROL_counts are in lines 4,6 and 8 of the first column
  upOrDownNL<-meanNLW-meanNLC
  #print(paste0("NL Up or Down?", upOrDownNL))
  # if you have a third factor...
  meanWPTW <- (data3[16,1]+data3[17,1]+data3[18,1])/3    #WPT_TREATMENT_counts are in lines 16,17 and 18 of the first column
  meanWPTC <- (data3[13,1]+data3[14,1]+data3[15,1])/3       #WPT_CONTROL_counts are in lines 13,14 and 15 of the first column
  upOrDownWPT<-meanWPTW-meanWPTC
  #print(paste0("WPT Up or Down? ", upOrDownWPT))
  meanWNLW <- (data3[20,1]+data3[22,1]+data3[24,1])/3    #WNL_TREATMENT_counts are in lines 20,22 and 24 of the first column
  meanWNLC <- (data3[19,1]+data3[21,1]+data3[23,1])/3    #WNL_CONTROL_counts are in lines 19,21 and 23 of the first column
  upOrDownWNL<-meanWNLW-meanWNLC
  #print(paste0("WNL Up or Down? ", upOrDownWNL))
  Expression <- data3[,1]   #expression is always in the first column, this way we can vary the number of factors without needing to alter this line of code
  options(contrasts=c("contr.sum","contr.poly"))
  glmmTMB1_binN = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History),data=data3,family="nbinom2"))
  AnovaGlmmTMB1_binN = try(as.data.frame(Anova(glmmTMB1_binN, type='III')))
  if (class(AnovaGlmmTMB1_binN) == "try-error") {
    Chisq <- c("ERROR", "ERROR", "ERROR","ERROR","ERROR", "ERROR", "ERROR","ERROR")
    Df <- c("ERROR", "ERROR","ERROR","ERROR","ERROR", "ERROR", "ERROR","ERROR")
    pvalue <- c("ERROR", "ERROR","ERROR","ERROR","ERROR", "ERROR", "ERROR","ERROR")
    
    AnovaGlmmTMB1_binN <- data.frame(Chisq, Df, pvalue)
  }
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownPT)
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownNL)
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownWPT)
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownWNL)
  
  # changing row names of data frame
  rownames(AnovaGlmmTMB1_binN)[9] <- "upOrDownPT"
  rownames(AnovaGlmmTMB1_binN)[10] <- "upOrDownNL"
  rownames(AnovaGlmmTMB1_binN)[11] <- "upOrDownWPT"
  rownames(AnovaGlmmTMB1_binN)[12] <- "upOrDownWNL"
  
  #Sort Using Character row.names
  #subset
  subsetted <- as.data.frame(subset(AnovaGlmmTMB1_binN, select = -c(`Chisq`,Df)))
  #change column name
  colnames(subsetted) <- c(geneName)
  #transpose
  transposed_of_AnovaGlmmTMB1_binN <- as.data.frame(t(subsetted),stringsAsFactors = FALSE)
  try(write.table(transposed_of_AnovaGlmmTMB1_binN, args[2] , append = TRUE, sep = '\t', col.names = FALSE, row.names = TRUE))

}

