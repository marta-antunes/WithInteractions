library(qqman) # Load the library qqman
library(CMplot) # Load the library CMplot
library(ggplot2)

rawdata <- read.table("inputFileNL.csv", sep="\t", header=TRUE) #load data
# if you get an error:
#check if file has div by zero (replace by NA)
#check if error in log2FC (replace by NA)
#check if chr NA replace by 0"
#check if there are dots in the name of chromosomes. Remove dots and every following character 
#P-value cannot be writen with "-"

data <- rawdata[with(rawdata, order(chr_name,meanPos)),] ## Sort by chr_name, then by meanPos

data <- data[!is.na(data$Pvalue),] #remove lines with Pvalue equal to NaN

data <- data[!is.na(data$Log2FC),] #remove lines with Log2FC equal to NA

#creates a new column named positiveOrNegativeLog10p. if Log2FC is negative, write log10(pvalue) which is also a negative value
data$positiveOrNegativeLog10p <- ifelse(data$Log2FC < 0, log10(data$Pvalue), -log10(data$Pvalue))

data <- subset(data,chr_name!='NC_045530' ) #keep only rows different from chr mitochondrial

#all chromossomes
# create list with DEG
genesToHighlight <- data[(data$Pvalue < 0.001),] #select genes with Pvalue<0.001 and abs(|log2FC|)>1
genesToHighlight_id<-genesToHighlight$gene


data <- data[, c("gene","chr_name", "meanPos", "positiveOrNegativeLog10p")] #only make this subset after generating genesToHighlight_id



# Make the Manhattan plot on the RNAseq dataset, atention that mean position of genes is being used
# Let's highlight the genes of interest.
#band is the space between chromosomes, aqui temos Highlight a Group of genes

#all chr
#amplify=FALSE because we dont want to make some points bigger
CMplot(data, type="p",plot.type="m", band=0.5, LOG10=FALSE, ylab="-log10(p) or log10(p)",ylim=c(-50,50),col=c("grey30","grey60"), highlight=genesToHighlight_id,
highlight.col="lightblue2",highlight.cex=1,highlight.pch=c(16), threshold=c(-log10(0.001),log10(0.001)),
threshold.lty=1, threshold.lwd=2, threshold.col="red", amplify=FALSE, width=14,height=6,
signal.col=NULL, chr.den.col=NULL, file="jpg",dpi=300,file.output=TRUE,
verbose=TRUE,cex=0.8)
