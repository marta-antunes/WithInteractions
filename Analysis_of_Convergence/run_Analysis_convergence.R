# G23 dataset

#necessary
library(lme4)
library(car)
library(glmmTMB)
library(bbmle) #necessary to AIC


# input_analysis4.csv was constructed based on "Galaxy232-Normalized_counts.tabular"
subtractions <- read.table("input_analysis4.csv", sep = ',', header=TRUE, stringsAsFactors = TRUE)

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- subtractions[ rowSums( subtractions > 0 ) >= 3, ]


genes <- read.csv("CandidatesLowLat.txt",sep='\n',header=FALSE)

#genes <- read.csv("CandidatesHighLat.txt",sep='\n',header=FALSE)

# Convert genes to a vector
genes <- as.vector(unlist(genes))

#subset dataInAtLeastXsamples
subsetData <- dataInAtLeastXsamples[rownames(dataInAtLeastXsamples) %in% genes, ]

#transpose
transposed_of_dataInAtLeastXsamples <- t(subsetData)    #transpose the rows and columns



#read samples
rownames <- rownames(transposed_of_dataInAtLeastXsamples)
Selection <- c(rep("C",3), rep("W",3))
Block <- c(rep("1",1), rep("2",1), rep("3",1), rep("1",1), rep("2",1), rep("3",1))

#Creating new data frame - name: NewColumns.Df,   
#Contains four columns History, Environment, AP and Selection
NewColumns.Df <- data.frame(Selection, Block)


data1 <- transposed_of_dataInAtLeastXsamples               # Replicate example data
data2 <- NewColumns.Df


#iterate over each gene in file
for(i in 1:ncol(data1)) {       # for-loop over columns
  geneName <- colnames(data1)[i]      # gene name
  data3 <- cbind(data1[ , i],data2)
  
  # we will consider up-regulated gene when is more expressed in NL than in PT
  meanControlPops <- (data3[1,1]+data3[2,1]+data3[3,1])/3    #TREATMENT_counts are in lines 3,4 and 5 of the first column
  meanWarmingPops <- (data3[4,1]+data3[5,1]+data3[6,1])/3    #TREATMENT_counts are in lines 7,8 and 9 of the first column
 
  
  Expression <- data3[,1]   #expression is always in the first column, this way we can vary the number of factors without needing to alter this line of code
  options(contrasts=c("contr.sum","contr.poly"))
  glmmTMB1_binN = try(glmmTMB(Expression ~ Selection+(1|Block),data=data3,family="nbinom2"))
  AnovaGlmmTMB1_binN = try(as.data.frame(Anova(glmmTMB1_binN, type='III')))
  if (class(AnovaGlmmTMB1_binN) == "try-error") {
    Chisq <- c("ERROR", "ERROR", "ERROR","ERROR")
    Df <- c("ERROR", "ERROR","ERROR","ERROR")
    pvalue <- c("ERROR", "ERROR","ERROR","ERROR")
    
    AnovaGlmmTMB1_binN <- data.frame(Chisq, Df, pvalue)
  }
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',meanControlPops)
  AnovaGlmmTMB1_binN[nrow(AnovaGlmmTMB1_binN) + 1,] = c('NA', 'NA',meanWarmingPops)
  
  # changing row names of data frame
  rownames(AnovaGlmmTMB1_binN)[3] <- "meanControlPops"
  rownames(AnovaGlmmTMB1_binN)[4] <- "meanWarmingPops"
  
  
  #Sort Using Character row.names
  #subset
  subsetted <- as.data.frame(subset(AnovaGlmmTMB1_binN, select = -c(`Chisq`,Df)))
  #change column name
  colnames(subsetted) <- c(geneName)
  #transpose
  transposed_of_AnovaGlmmTMB1_binN <- as.data.frame(t(subsetted),stringsAsFactors = FALSE)
  try(write.table(transposed_of_AnovaGlmmTMB1_binN, "diferenca_candidatesLowLat_ratio.csv" , append = TRUE, sep = '\t', col.names = FALSE, row.names = TRUE))
  #try(write.table(transposed_of_AnovaGlmmTMB1_binN, args[2] , append = TRUE, sep = '\t', col.names = FALSE, row.names = TRUE))

}

