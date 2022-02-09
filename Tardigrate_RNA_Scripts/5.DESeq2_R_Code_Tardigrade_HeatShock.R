#5.DESeq2_R_Code_Tardigrade_HeatShock.R

# remove anything before we start
rm(list = ls())

# Download and install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("DESeq2")
BiocManager::install("DEFormats")

install.packages("xlsx")
install.packages("statmod")
library(DESeq2)
library(data.table)
library(edgeR)
library(limma)
library(xlsx)
library(statmod)
library(DEFormats)

# set the working directory
setwd("directory")

############################
### load phenotypic data ###
############################

# Load phenotypic data: should contain sample name, phenotype, date of extraction, person who extracted the RNA, lane used for sequencing...etc
colData <- read.csv("Samplenames_Extractiondates_4.csv", header=TRUE)

# Format phenotypic data
colData$SampleID <- as.factor(colData$SampleID)
colData$Group <- as.factor(colData$Group) # control or experiment
colData$Experiment <- as.factor(colData$Experiment) # Heatshock 
colData$Duration <- as.factor(colData$Duration) # 2h or 24h
colData$Thawing_Date <- as.factor(colData$Thawing_Date)
colData$Extraction_Date <- as.factor(colData$Extraction_Date)
colData$Extractor <- as.factor(colData$Extractor)
colData$Nr_Specimens <- as.integer(colData$Nr_Specimens)
#colData$Lane <- as.factor(colData$Lane)

colData <- colData[,c(1:12)] # be careful and check when using this!


#################################################################
### load count data or depth data (with or without row.names) ###
#################################################################

# Load a count matrix
countData <- read.csv("GeneCount_Clean.csv", header=TRUE, row.names=1)

#########################
# Heatshock Subset data #
#########################

colData_HeatShock <- colData[grep("HeatShock", colData$Experiment), ]

ID_HeatShock <- as.character(colData_HeatShock$SampleID)

countData_HeatShock <- countData[,ID_HeatShock]

####################
### Start DESeq2 ###
####################

# rename objects HeatShock
colData <- colData_HeatShock
countData <- countData_HeatShock

# Change the count table to integer as DESeq2 doesn't allow decimals
countData[] <- lapply(countData, as.integer)

# Full data - make groups, and put that into "design"
# https://support.bioconductor.org/p/92941/
colData$group <- factor(paste0(colData$Group, "_", colData$Duration))

#this command allows linear regression, we now choose the group with small g
#we look at the countData and comes the counts to the group
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)
model.matrix(~ colData$group)
dds$group<-relevel(dds$group,ref="Ctl_24h") #HeatShock to set the control as reference

#The same as above but for HeatShock
dds$group <- factor(dds$group, levels = c("Ctl_24h","Exp_24h"))


### 4.2The variance stabilizing transformation and the rlog ###
# The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30).
#removes outliers 

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

#install.packages("dplyr")
#install.packages("ggplot2")
library("dplyr")
library("ggplot2")

### 4.3 Sample distances ###

# calculate the Euclidean distance between samples. To ensure we have a roughly equal contribution from all genes, we use it on the VST data.
sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$Group, vsd$Experiment, vsd$Duration, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

### 4.4PCA plot
#this allows to us to see variations 
# PC1 low --> high noise 
# plots can be customized using ggplot function too
plotPCA(vsd, intgroup=c("Group", "Experiment", "Duration"))
plotPCA(vsd, intgroup = "Thawing_Date")
plotPCA(vsd, intgroup = "Extraction_Date")
plotPCA(vsd, intgroup = "Extractor")

pcaData <- plotPCA(vsd, intgroup=c("Group", "Experiment", "Duration","Extraction_Date"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Group, shape=Duration)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#plotPCA(vsd, intgroup = "Lane")
#Look at the other variances and looks for the batch effect


### 5.1Running the differential expression pipeline
# https://support.bioconductor.org/p/63201/
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)


#this allows comparison between the groups we want
results(dds, contrast=c("group","Exp_24h","Ctl_24h"))


###

# check levels
levels(dds$group)

######

# gene with lowest p-value
plotCounts(dds, gene=which.min(res$padj), intgroup="group")

#######

# Ctl_24h vs. Exp_24h 
res_Ctl24h_Exp24h <- results(dds, contrast=c("group","Ctl_24h", "Exp_24h"), alpha=0.05) 
summary(res_Ctl24h_Exp24h)
sum(res_Ctl24h_Exp24h$padj < 0.05, na.rm=TRUE)
res_Ctl24h_Exp24h_ordered <- res_Ctl24h_Exp24h[order(res_Ctl24h_Exp24h$padj),]
write.csv(as.data.frame(res_Ctl24h_Exp24h_ordered), file="DESeq2_res_Ctl24h_Exp24h_HeatShock.csv") #Heatshock

plotCounts(dds, gene=which.min(res_Ctl24h_Exp24h$padj), intgroup="group")


#########################################################################
##### Creating log-fold change against log-counts per million plot ######
#########################################################################

#Plot log-fold change against log-counts per million, with DE genes highlighted:
dev.off()
DESeq2::plotMA(res, main = "MA-Plot Heating 24h", ylim = c(-4, 4))


############################# 
# Check Ctl_24h vs. Exp_24h #
# DE genes - 226            #
#############################

colData_24h <- colData[ which(colData$Duration =='24h'), ]
ID_24h <- as.character(colData_24h$SampleID)
countData_24h <- countData[,ID_24h]

# rename objects
colData <- colData_24h
countData <- countData_24h

# Change the count table to integer as DESeq2 doesn't allow decimals
countData[] <- lapply(countData, as.integer)

# Full data - make groups, and put that into "design"
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Group)
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
