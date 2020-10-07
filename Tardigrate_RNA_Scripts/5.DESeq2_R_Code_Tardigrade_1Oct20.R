# remove anything before we start
rm(list = ls())

# get the right path
.libPaths()
.libPaths(.libPaths()[2:2])
.libPaths()

# Download and install packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("DESeq2")
#BiocManager::install("DEFormats")

install.packages("xlsx")

library(DESeq2)
library(data.table)
library(edgeR)
library(limma)
library(xlsx)
library(statmod)
library(DEFormats)

# set the working directory
setwd("\\\\a00519.science.domain/mch497/Documents/Tardigrade")

############################
### load phenotypic data ###
############################

# Load phenotypic data: should contain sample name, phenotype, date of extraction, person who extracted the RNA, lane used for sequencing...etc
colData <- read.csv("Samplenames_Extractiondates_4.csv", header=TRUE)

# Format phenotypic data
colData$SampleID <- as.factor(colData$SampleID)
colData$Group <- as.factor(colData$Group) # control or experiment
colData$Experiment <- as.factor(colData$Experiment) # Heatshock or freeze
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

####################################
# Heatshock or Freeze? Subset data #
####################################

colData_HeatShock <- colData[grep("HeatShock", colData$Experiment), ]
colData_Freeze <- colData[grep("Freeze", colData$Experiment), ]

ID_HeatShock <- as.character(colData_HeatShock$SampleID)
ID_Freeze <- as.character(colData_Freeze$SampleID)

countData_HeatShock <- countData[,ID_HeatShock]
countData_Freeze <- countData[,ID_Freeze]

###################################################
# Further filter out the count data before DESeq2 #
# SKIP THIS STEP IF TAKEN CARE OF ALREADY BEFORE  #
###################################################

# check if any genes have a total count of zero
table(colSums(countData_missing_removed) == 0)

# check if any individuals that have lots of zero counts
which.max(colSums(countData_missing_removed) == 0)

# at least 50% of the samples with a count of 10 or higher
# this removes individuals with too many missing counts
keep <- rowSums((countData_missing_removed) >= 10) >= 8 # 19/2 = 8

countData_missing_removed <- countData_missing_removed[keep,]
nrow(countData_missing_removed)

# Change the count table to integer as DESeq2 doesn't allow decimals
countData_missing_removed[] <- lapply(countData_missing_removed, as.integer)

####################
### Start DESeq2 ###
####################

# rename objects
colData <- colData_Freeze
countData <- countData_Freeze

# Change the count table to integer as DESeq2 doesn't allow decimals
countData[] <- lapply(countData, as.integer)

# Full data - make groups, and put that into "design"
# https://support.bioconductor.org/p/92941/

colData$group <- factor(paste0(colData$Group, "_", colData$Duration))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 0 + group)
#dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group) # since there is no "reference" in this case, above with "+0" is more straightforward (results are the same)
#model.matrix(~ 0 + colData$group)
#dds$group<-relevel(dds$group,ref="Ctl_24h")
dds$group <- factor(dds$group, levels = c("Ctl_2h","Exp_2h","Ctl_24h","Exp_24h"))

### 4.2The variance stabilizing transformation and the rlog ###
# The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30).
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

#rld <- rlog(dds, blind = FALSE)  # took 50 minutes and didn't finish
#head(assay(rld), 3)

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
# plots can be customized using ggplot function too
plotPCA(vsd, intgroup=c("Group", "Experiment", "Duration"))
plotPCA(vsd, intgroup = "Thawing_Date")
plotPCA(vsd, intgroup = "Extraction_Date")
plotPCA(vsd, intgroup = "Extractor")
#plotPCA(vsd, intgroup = "Lane")

### 5.1Running the differential expression pipeline
# https://support.bioconductor.org/p/63201/
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
summary(res)

results(dds, contrast=c("group","Exp_2h","Ctl_2h"))
results(dds, contrast=c("group","Exp_24h","Ctl_24h"))
results(dds, contrast=c("group","Ctl_2h","Ctl_24h"))
results(dds, contrast=c("group","Exp_2h","Exp_24h"))

###

# check levels
levels(dds$group)

######

# plot
#plotMA(res, ylim=c(-2,2))

# gene with lowest p-value
plotCounts(dds, gene=which.min(res$padj), intgroup="group")

#######

# Ctl_2h vs. Ctl_24h - 20 DE genes
res_Ctl2h_Ctl24h <- results(dds, contrast=c("group","Ctl_2h","Ctl_24h"), alpha=0.05) 
summary(res_Ctl2h_Ctl24h)
sum(res_Ctl2h_Ctl24h$padj < 0.05, na.rm=TRUE)
res_Ctl2h_Ctl24h_ordered <- res_Ctl2h_Ctl24h[order(res_Ctl2h_Ctl24h$padj),]
head(res_Ctl2h_Ctl24h_ordered)
write.csv(as.data.frame(res_Ctl2h_Ctl24h_ordered), file="DESeq2_res_Ctl2h_Ctl24h.csv")

plotCounts(dds, gene=which.min(res_Ctl2h_Ctl24h$padj), intgroup="group")

# Exp_2h vs. Exp_24h - 47 DE genes
res_Exp2h_Exp24h <- results(dds, contrast=c("group","Exp_2h","Exp_24h"), alpha=0.05) 
summary(res_Exp2h_Exp24h)
sum(res_Exp2h_Exp24h$padj < 0.05, na.rm=TRUE)
res_Exp2h_Exp24h_ordered <- res_Exp2h_Exp24h[order(res_Exp2h_Exp24h$padj),]
head(res_Exp2h_Exp24h_ordered)
write.csv(as.data.frame(res_Exp2h_Exp24h_ordered), file="DESeq2_res_Exp2h_Exp24h.csv")

plotCounts(dds, gene=which.min(res_Exp2h_Exp24h$padj), intgroup="group")

# Ctl_2h vs. Exp_2h - 2 DE genes
res_Ctl2h_Exp2h <- results(dds, contrast=c("group","Ctl_2h","Exp_2h"), alpha=0.05) 
summary(res_Ctl2h_Exp2h)
sum(res_Ctl2h_Exp2h$padj < 0.05, na.rm=TRUE)
res_Ctl2h_Exp2h_ordered <- res_Ctl2h_Exp2h[order(res_Ctl2h_Exp2h$padj),]
write.csv(as.data.frame(res_Ctl2h_Exp2h_ordered), file="DESeq2_res_Ctl2h_Exp2h.csv")

plotCounts(dds, gene=which.min(res_Ctl2h_Exp2h$padj), intgroup="group")


# Ctl_24h vs. Exp_24h - 128 DE genes
res_Ctl24h_Exp24h <- results(dds, contrast=c("group","Ctl_24h", "Exp_24h"), alpha=0.05) 
summary(res_Ctl24h_Exp24h)
sum(res_Ctl24h_Exp24h$padj < 0.05, na.rm=TRUE)
res_Ctl24h_Exp24h_ordered <- res_Ctl24h_Exp24h[order(res_Ctl24h_Exp24h$padj),]
write.csv(as.data.frame(res_Ctl24h_Exp24h_ordered), file="DESeq2_res_Ctl24h_Exp24h.csv")

plotCounts(dds, gene=which.min(res_Ctl24h_Exp24h$padj), intgroup="group")

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
