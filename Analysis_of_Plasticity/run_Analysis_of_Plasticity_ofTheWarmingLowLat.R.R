# Analysis of Plasticity
#necessary
library(lme4)
library(car)
library(glmmTMB)
library(bbmle) #necessary to AIC


#args <- commandArgs(TRUE)
##tive de remover o espaço no início do header
#normalizedCounts <- read.table(args[1], sep = '\t', header=TRUE, stringsAsFactors = TRUE)
normalizedCounts <- read.table("Galaxy359-Normalized_counts.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)


#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]


#transpose
transposed_of_normalizedCounts <- t(dataInAtLeastXsamples)    #transpose the rows and columns

#read samples
rownames <- rownames(transposed_of_normalizedCounts)
splited <- strsplit(rownames, split = "_")
#History<-substr(lapply(splited, `[[`, 2), start = 1, stop = 2) #list of lists. history is nested in the second element of list  #lapply gives the second element of each list #substr takes the first two characters of each string
Environment<-lapply(splited, `[[`, 4) #Environment is the fourth element of list
Environment<-as.character(Environment) #transform list in character
AP <-lapply(splited, `[[`, 2)  #AP is the second element of list
AP<-as.character(AP) #transform list in character
#Selection <-lapply(splited, `[[`, 1) #selection is the first element of list
#Selection<-as.character(Selection) #transform list in character

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
  #print(paste0("PT Up or Down? ", upOrDownPT))
  #meanNLW <- (data3[5,1]+data3[7,1]+data3[9,1])/3    #NL_TREATMENT_counts are in lines 5,7 and 9 of the first column
  #meanNLC <- (data3[4,1]+data3[6,1]+data3[8,1])/3    #NL_CONTROL_counts are in lines 4,6 and 8 of the first column
  #upOrDownNL<-meanNLW-meanNLC
  #print(paste0("NL Up or Down?", upOrDownNL))
  # if you have a third factor...
  #meanWPTW <- (data3[7,1]+data3[8,1]+data3[9,1])/3    #WPT_TREATMENT_counts are in lines 7,8 and 9 of the first column
  #meanWPTC <- (data3[4,1]+data3[5,1]+data3[6,1])/3       #WPT_CONTROL_counts are in lines 4,5 and 6 of the first column
  #upOrDownWPT<-meanWPTW-meanWPTC
  #print(paste0("WPT Up or Down? ", upOrDownWPT))
  #meanWNLW <- (data3[20,1]+data3[22,1]+data3[24,1])/3    #WNL_TREATMENT_counts are in lines 20,22 and 24 of the first column
  #meanWNLC <- (data3[19,1]+data3[21,1]+data3[23,1])/3    #WNL_CONTROL_counts are in lines 19,21 and 23 of the first column
  #upOrDownWNL<-meanWNLW-meanWNLC
  #print(paste0("WNL Up or Down? ", upOrDownWNL))
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
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownPT)
  #AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',upOrDownWPT)

  
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
  #try(write.table(transposed_of_AnovaGlmmTMB1_binN, args[2] , append = TRUE, sep = '\t', col.names = FALSE, row.names = TRUE))

}

