# G23 dataset

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

par(mar=c(4,4,1,0), mfrow=c(1,1), cex=1.0, cex.main=0.8, cex.axis=0.8)

# Set up layout to have space for the legend
layout(matrix(c(1,2), nrow=1), widths=c(4,2))  # 4:1 ratio of plot to legend space
par(mar=c(5,4,4,0) + 0.1)  # Adjust margins

#Plots scatter plot for PC 1 and 2
#high latitude populations represented by empty center
plot(project.pca$x, type="n", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col=c('blue','blue','blue','blue','red','blue','red','blue','red','red','red','red','blue','blue','blue','red','red','red','blue','red','blue','red','blue','red'), pch=c(16,16,16,21,21,21,21,21,21,16,16,16,15,15,15,15,15,15,22,22,22,22,22,22), cex=1)
#initially we included the name of the samples but them we realize we dont need, if you want to rescue that uncoment the line below
#text(project.pca$x, labels=rownames(project.pca$x), cex=0.6, font=2, pos=2)
# Plot the legend in the separate space
par(mar=c(5,0,4,0))  # Adjust margins for the legend
plot.new()  # Create a new plot
legend("center", legend=c("Warming environment", "Control environment", "Control Regime", "Warming Regime", "Low latitude", "High latitude"), col=c("red", "blue", "grey38", "grey38", "grey38", "grey38"), lty = c(0, 0, 0), cex=0.6, pch = c(16, 16, 16, 15, 16, 21))


#code used for pvca
library(Biobase)
library(ExpressionNormalizationWorkflow)

env <- sapply(rownames(project.pca$x), function(x) substr(x, nchar(x), nchar(x)))
pop <- c(sapply(rownames(project.pca$x[1:12,]), function(x) substr(x, 8,9)),
         sapply(rownames(project.pca$x[13:24,]), function(x) substr(x, 9,10)))
evo <- c(rep("control",12),rep("evolved",12))
meta.table<-data.frame(evo=evo,pop=pop,env = env, row.names = rownames(project.pca$x))
annot<-data.frame(labelDescription=c("Factor levels","Factor levels","Factor levels"))
annot_factors<-AnnotatedDataFrame(data = meta.table,varMetadata = annot)
expr.set<-ExpressionSet(assayData = as.matrix(dataInAtLeastXsamples),phenoData = annot_factors)
pvca_res<-pvcAnaly(expr.set, 0.53,c("evo","pop","env"))
## I used 0.53 as the threshold since the first three PC explain about 52.6% of the variance.
