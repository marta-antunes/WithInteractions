# G23 dataset

library(lme4)
library(car)
library(glmmTMB)
library(bbmle) #necessary to AIC

##verify if there is no the space in the beggining of the header
normalizedCounts <- read.table("Galaxy238.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)


#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]


#transpose
transposed_of_normalizedCounts <- t(dataInAtLeastXsamples)    #transpose the rows and columns


#Creating new data frame - name: NewColumns.Df,   
#Contains two columns History and Environment  
NewColumns.Df <- data.frame(  
  History = c(rep("PT",3), rep("NL",6), rep("PT",3), rep("PT",6), rep("NL",6)),
  Environment = c(rep("C",4),rep("W",1),rep("C",1), rep("W",1), rep("C",1),rep("W",4),rep("C",3), rep("W",3),rep("C",1), rep("W",1),rep("C",1), rep("W",1),rep("C",1), rep("W",1)),
  AP = c("PT1", "PT2", "PT3", "NL1", "NL1", "NL2", "NL2", "NL3", "NL3", "PT1", "PT2", "PT3", "PT1", "PT2", "PT3", "PT1", "PT2", "PT3", "NL1", "NL1", "NL2","NL2", "NL3", "NL3"),
  Selection = c(rep("C",3), rep("C",6), rep("C",3), rep("W",6), rep("W",6)))


data1 <- transposed_of_normalizedCounts               # Replicate example data
data2 <- NewColumns.Df

#sink("20220803_output.txt", append = TRUE)
for(i in 1:ncol(data1)) {       # for-loop over columns
  geneName <- colnames(data1)[i]
  #print(geneName)     # gene name
  data3 <- cbind(data1[ , i],data2)
  meanPTW <- (data3[10,1]+data3[11,1]+data3[12,1])/3    #PT_TREATMENT_counts are in lines 10,11 and 12 of the first column
  meanPTC <- (data3[1,1]+data3[2,1]+data3[3,1])/3       #PT_CONTROL_counts are in lines 1,2 and 3 of the first column
  upOrDownPT<-meanPTW-meanPTC
  #print(paste0("PT Up or Down? ", upOrDownPT))
  meanNLW <- (data3[5,1]+data3[7,1]+data3[9,1])/3    #NL_TREATMENT_counts are in lines 5,7 and 9 of the first column
  meanNLC <- (data3[4,1]+data3[6,1]+data3[8,1])/3    #NL_CONTROL_counts are in lines 4,6 and 8 of the first column
  upOrDownNL<-meanNLW-meanNLC
  # if you have a third factor...
  #print(paste0("NL Up or Down?", upOrDownNL))
  meanWPTW <- (data3[16,1]+data3[17,1]+data3[18,1])/3    #WPT_TREATMENT_counts are in lines 16,17 and 18 of the first column
  meanWPTC <- (data3[13,1]+data3[14,1]+data3[15,1])/3       #WPT_CONTROL_counts are in lines 13,14 and 15 of the first column
  upOrDownWPT<-meanWPTW-meanWPTC
  #print(paste0("WPT Up or Down? ", upOrDownWPT))
  meanWNLW <- (data3[20,1]+data3[22,1]+data3[24,1])/3    #WNL_TREATMENT_counts are in lines 20,22 and 24 of the first column
  meanWNLC <- (data3[19,1]+data3[21,1]+data3[23,1])/3    #WNL_CONTROL_counts are in lines 19,21 and 23 of the first column
  upOrDownWNL <- meanWNLW-meanWNLC
  #print(paste0("WNL Up or Down? ", upOrDownWNL))
  #print("NL Up or Down?", upOrDownNL)
  Expression <- data3[,1]   #expression is always in the first column, this way we can vary the number of factors without needing to alter this line of code
  options(contrasts=c("contr.sum","contr.poly"))
  glmmTMB1_binN = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History),data=data3,family="nbinom2"))
  #try(print(Anova(glmmTMB1_binN, type='III')))
  glmmTMB2_binN = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History)+(1|AP:History:Environment),data=data3,family="nbinom2"))
  #try(print(Anova(glmmTMB2_binN, type='III')))
  glmmTMB3_binN = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History)+(1|AP:History:Environment)+(1|AP:History:Selection),data=data3,family="nbinom2"))
  #try(print(Anova(glmmTMB3_binN, type='III')))
  glmmTMB4_binN = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History)+(1|AP:History:Environment)+(1|AP:History:Selection)+(1|AP:History:Selection:Environment),data=data3,family="nbinom2"))
  #try(print(Anova(glmmTMB4_binN, type='III')))
  glmmTMB1_poisson = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History),data=data3,family="poisson"))
  #try(print(Anova(glmmTMB1_poisson, type='III')))
  glmmTMB2_poisson = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History)+(1|AP:History:Environment),data=data3,family="poisson"))
  #try(print(Anova(glmmTMB2_poisson, type='III')))
  glmmTMB3_poisson = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History)+(1|AP:History:Environment)+(1|AP:History:Selection),data=data3,family="poisson"))
  #try(print(Anova(glmmTMB3_poisson, type='III')))
  glmmTMB4_poisson = try(glmmTMB(Expression ~ History*Environment*Selection+(1|AP:History)+(1|AP:History:Environment)+(1|AP:History:Selection)+(1|AP:History:Selection:Environment),data=data3,family="poisson"))
  #try(print(Anova(glmmTMB4_poisson, type='III')))
  #try(print(AICtab(glmmTMB1_binN,glmmTMB2_binN,glmmTMB3_binN,glmmTMB4_binN,glmmTMB1_poisson,glmmTMB2_poisson,glmmTMB3_poisson,glmmTMB4_poisson)))
  AIC <- try(as.data.frame(AICtab(glmmTMB1_binN,glmmTMB2_binN,glmmTMB3_binN,glmmTMB4_binN,glmmTMB1_poisson,glmmTMB2_poisson,glmmTMB3_poisson,glmmTMB4_poisson)))
  
  #dealing with errors
  if (class(AIC) == "try-error") {
    AIC <- data.frame(genes=geneName,glmmTMB1_binN="ERROR_in_AIC",glmmTMB2_binN="ERROR_in_AIC", glmmTMB3_binN="ERROR_in_AIC",glmmTMB4_binN="ERROR_in_AIC",glmmTMB1_poisson="ERROR_in_AIC",glmmTMB2_poisson="ERROR_in_AIC", glmmTMB3_poisson="ERROR_in_AIC",glmmTMB4_poisson="ERROR_in_AIC", bestModel="ERROR_in_AIC")
    write.table(AIC, args[2] , append = TRUE, sep = '\t', col.names = FALSE, row.names = FALSE)
  } else {
  #Sort Using Character row.names
  try(sorted <- AIC[order(row.names(AIC)), ])
  
  #print the best
  try(best <- rownames(sorted[which.min(sorted$dAIC),]))
  #add row to dataframe
  try(sorted[nrow(sorted) + 1,] <- c(best, 'NA'))
  # changing row names of data frame
  try(rownames(sorted)[9] <- "bestModel")
  
  #subset
  try(AICsubsetted <- subset(sorted, select = -c(df) ))
  
  
  #change column name
  try(colnames(AICsubsetted)[1] <- geneName)
  
  #transpose
  try(transposed_of_AICsubsetted <- as.data.frame(t(AICsubsetted),stringsAsFactors = FALSE, ))
  #try(print(transposed_of_AICsubsetted))
  try(write.table(transposed_of_AICsubsetted, args[2] , append = TRUE, sep = '\t', col.names = FALSE, row.names = TRUE))
}}

#sink()

