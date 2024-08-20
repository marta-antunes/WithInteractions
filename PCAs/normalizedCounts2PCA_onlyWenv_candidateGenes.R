# G23 dataset

#read file with normalized counts
normalizedCounts <- read.table("Galaxy238_changedLabels.tabular", sep = '\t', header=TRUE, stringsAsFactors = TRUE)

onlyWarmingEnv <- normalizedCounts[, c("PT1_W", "PT2_W", "PT3_W", "NL1_W", "NL2_W", "NL3_W", "WPT1_W", "WPT2_W", "WPT3_W", "WNL1_W", "WNL2_W", "WNL3_W")]

#keep genes that have data in at least 3 of the samples (columns).
dataInAtLeastXsamples <- onlyWarmingEnv[ rowSums( onlyWarmingEnv > 0 ) >= 3, ]


#genes <- read.csv("CandidatesLowLat.txt",sep='\n',header=FALSE)

genes <- read.csv("CandidatesHighLat.txt",sep='\n',header=FALSE)

# Convert genes to a vector
genes <- as.vector(unlist(genes))

#subset dataInAtLeastXsamples
subsetData <- dataInAtLeastXsamples[rownames(dataInAtLeastXsamples) %in% genes, ]


#make PCA with the prcomp() function
#This function takes a matrix of data, where the columns are the variables that we want to use to transform our samples, which should be the rows of the matrix. So we need to transpose our matrix.
#We Scaled the expression matrix before PCA 
project.pca <- prcomp(t(subsetData), scale. = TRUE)
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
#plot(project.pca$x, type="n", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), main = "Candidate Genes in Low Latitude Populations")
plot(project.pca$x, type="n", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), main = "Candidate Genes in High Latitude Populations")
points(project.pca$x, col=c('grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38','grey38'), pch=c(16,16,16,21,21,21,15,15,15,22,22,22), cex=1)
#initially we included the name of the samples but them we realize we dont need, if you want to rescue that uncoment the line below
#text(project.pca$x, labels=rownames(project.pca$x), cex=0.6, font=2, pos=2)
# Get coordinates of NL1_W and WNL1_W for the arrow
NL1_W_coords <- project.pca$x["NL1_W",]
WNL1_W_coords <- project.pca$x["WNL1_W",]
NL2_W_coords <- project.pca$x["NL2_W",]
WNL2_W_coords <- project.pca$x["WNL2_W",]
NL3_W_coords <- project.pca$x["NL3_W",]
WNL3_W_coords <- project.pca$x["WNL3_W",]

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
PT1_W_coords <- project.pca$x["PT1_W",]
WPT1_W_coords <- project.pca$x["WPT1_W",]
PT2_W_coords <- project.pca$x["PT2_W",]
WPT2_W_coords <- project.pca$x["WPT2_W",]
PT3_W_coords <- project.pca$x["PT3_W",]
WPT3_W_coords <- project.pca$x["WPT3_W",]

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

