# G23 dataset


#necessary
library(lme4)
library(car)
library(glmmTMB)
library(bbmle) #necessary to AIC

#read file with normalized counts
normalizedCounts <- read.table("Galaxy238_changedLabels.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- normalizedCounts[ rowSums( normalizedCounts > 0 ) >= 3, ]

#make PCA with the prcomp() function
#This function takes a matrix of data, where the columns are the variables that we want to use to transform our samples, which should be the rows of the matrix. So we need to transpose our matrix.
#We Scaled the expression matrix before PCA 
project.pca <- prcomp(t(dataInAtLeastXsamples), scale. = TRUE)
summary(project.pca)

#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100

par(mar=c(4,4,1,1), mfrow=c(1,1), cex=1.0, cex.main=0.8, cex.axis=0.8)

#Plots scatter plot for PC 1 and 2
plot(project.pca$x, type="n", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col=c('blue','blue','blue','blue','red','blue','red','blue','red','red','red','red','blue','blue','blue','red','red','red','blue','red','blue','red','blue','red'), pch=c(16,16,16,16,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15,15,15,15,15), cex=1,)
#initially we included the name of the samples but them we realize we dont need, if you want to rescue that uncoment the line below
#text(project.pca$x, labels=rownames(project.pca$x), cex=0.6, font=2, pos=2)
legend("topright",legend=c("Warming environment", "Control environment", "Control Regime", "Warming Regime"),col=c("red", "blue","grey38", "grey38"), lty = c(0, 0, 0), cex=0.8, pch = c(16, 16, 16, 15))

