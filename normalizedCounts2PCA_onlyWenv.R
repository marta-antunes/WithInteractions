# G23 dataset

#read file with normalized counts
normalizedCounts <- read.table("Galaxy238.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

onlyWarmingEnv <- normalizedCounts[, c("counts_PT1_G23_W", "counts_PT2_G23_W", "counts_PT3_G23_W", "counts_NL1_G23_W", "counts_NL2_G23_W", "counts_NL3_G23_W", "counts_WPT1_G23_W", "counts_WPT2_G23_W", "counts_WPT3_G23_W", "counts_WNL1_G23_W", "counts_WNL2_G23_W", "counts_WNL3_G23_W")]

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- onlyWarmingEnv[ rowSums( onlyWarmingEnv > 0 ) >= 3, ]




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
points(project.pca$x, col=c('grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38'), pch=c(16,16,16,21,21,21,15,15,15,22,22,22), cex=1)
#initially we included the name of the samples but them we realize we dont need, if you want to rescue that uncoment the line below
#text(project.pca$x, labels=rownames(project.pca$x), cex=0.6, font=2, pos=2)
# Get coordinates of NL1_W and WNL1_W for the arrow
NL1_W_coords <- project.pca$x["counts_NL1_G23_W",]
WNL1_W_coords <- project.pca$x["counts_WNL1_G23_W",]
NL2_W_coords <- project.pca$x["counts_NL2_G23_W",]
WNL2_W_coords <- project.pca$x["counts_WNL2_G23_W",]
NL3_W_coords <- project.pca$x["counts_NL3_G23_W",]
WNL3_W_coords <- project.pca$x["counts_WNL3_G23_W",]

# Add the arrow linking NL1_W to WNL1_W
arrows(x0 = NL1_W_coords[1], y0 = NL1_W_coords[2], 
       x1 = WNL1_W_coords[1], y1 = WNL1_W_coords[2], 
       length = 0.1, angle = 30, code = 2, col = "blue", lwd = 2)
# Add the arrow linking NL2_W to WNL2_W
arrows(x0 = NL2_W_coords[1], y0 = NL2_W_coords[2], 
       x1 = WNL2_W_coords[1], y1 = WNL2_W_coords[2], 
       length = 0.1, angle = 30, code = 2, col = "blue", lwd = 2)
# Add the arrow linking NL3_W to WNL3_W
arrows(x0 = NL3_W_coords[1], y0 = NL3_W_coords[2], 
       x1 = WNL3_W_coords[1], y1 = WNL3_W_coords[2], 
       length = 0.1, angle = 30, code = 2, col = "blue", lwd = 2)

# Get coordinates of PT1_W and WPT1_W for the arrow
PT1_W_coords <- project.pca$x["counts_PT1_G23_W",]
WPT1_W_coords <- project.pca$x["counts_WPT1_G23_W",]
PT2_W_coords <- project.pca$x["counts_PT2_G23_W",]
WPT2_W_coords <- project.pca$x["counts_WPT2_G23_W",]
PT3_W_coords <- project.pca$x["counts_PT3_G23_W",]
WPT3_W_coords <- project.pca$x["counts_WPT3_G23_W",]

# Add the arrow linking PT1_W to WPT1_W
arrows(x0 = PT1_W_coords[1], y0 = PT1_W_coords[2], 
       x1 = WPT1_W_coords[1], y1 = WPT1_W_coords[2], 
       length = 0.1, angle = 30, code = 2, col = "orange", lwd = 2)
# Add the arrow linking PT2_W to WPT2_W
arrows(x0 = PT2_W_coords[1], y0 = PT2_W_coords[2], 
       x1 = WPT2_W_coords[1], y1 = WPT2_W_coords[2], 
       length = 0.1, angle = 30, code = 2, col = "orange", lwd = 2)
# Add the arrow linking PT3_W to WPT3_W
arrows(x0 = PT3_W_coords[1], y0 = PT3_W_coords[2], 
       x1 = WPT3_W_coords[1], y1 = WPT3_W_coords[2], 
       length = 0.1, angle = 30, code = 2, col = "orange", lwd = 2)

# Plot the legend in the separate space
par(mar=c(5,0,4,0))  # Adjust margins for the legend
plot.new()  # Create a new plot
legend("center", legend=c("Control Regime", "Warming Regime", "Low latitude", "High latitude"), col=c("grey38", "grey38", "grey38", "grey38"), lty = c(0, 0, 0), cex=0.6, pch = c(16, 15, 16, 21))

#code used for pvca
library(Biobase)
library(ExpressionNormalizationWorkflow)


pop <- c(sapply(rownames(project.pca$x[1:6,]), function(x) substr(x, 8,9)),
         sapply(rownames(project.pca$x[7:12,]), function(x) substr(x, 9,10)))
evo <- c(rep("control",6),rep("evolved",6))
meta.table<-data.frame(evo=evo,pop=pop, row.names = rownames(project.pca$x))
annot<-data.frame(labelDescription=c("Factor levels","Factor levels"))
annot_factors<-AnnotatedDataFrame(data = meta.table,varMetadata = annot)
expr.set<-ExpressionSet(assayData = as.matrix(dataInAtLeastXsamples),phenoData = annot_factors)
pvca_res<-pvcAnaly(expr.set, 0.64,c("evo","pop"))
## I used 0.64 as the threshold since the first three PC explain about 63.7% of the variance.

