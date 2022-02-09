#6.EdgeR_R_Code_Tardigrade_HeatShock.R

# remove anything before we start
rm(list = ls())

# Install packages
library(DESeq2)
library(data.table)
library(edgeR)
library(limma)
library(xlsx)
library(statmod)
library(DEFormats)

############################

setwd("Directory")

# Load phenotypic data: should contain sample name, phenotype, date of extraction, person who extracted the RNA, lane used for sequencing...etc
colData <- read.csv("Samplenames_Extractiondates_4.csv", header=TRUE)

# Format phenotypic data
colData$SampleID <- as.factor(colData$SampleID)
colData$Group <- as.factor(colData$Group) # control or experiment
colData$Experiment <- as.factor(colData$Experiment) 
colData$Duration <- as.factor(colData$Duration) 
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

# rename objects
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
dds$group<-relevel(dds$group,ref="Ctl_24h")

#here we insure that we have 2 groups in the DE data set and we set the references, it makes sure that R understands that we have 2 different groups
dds$group <- factor(dds$group, levels = c("Ctl_24h","Exp_24h"))

##############################################################################################
### Start EdgeR                                                                            ###
# NOTE: the following commands should be right after DESeq2, in order to use saved objects ###
##############################################################################################

### on the whole countData ###
x <- DGEList(counts=countData, group=colData$group)

# Normalize for RNA composition
x <- calcNormFactors(x)
x$samples

# Compute counts per million (CPM)
cpm <- cpm(x)

# MDS (Plot samples on a two-dimensional scatterplot)
#plotMDS(x, col=group)
plotMDS(x)

# Design matrix
design <- model.matrix(~ 0 + group, data=colData)
design
colnames(design) # to get the group names to use for contrasts

# Estimate the dispersions
x <- estimateDisp(x, design, robust=TRUE)

#Scatterplot of the biological coefficient of variation (BCV) against the average abundance of each gene.
#stabilisations of the variables 
plotBCV(x)


################################################
##### To perform quasi-likelihood F-tests ######
################################################

fit <- glmQLFit(x, design, robust=TRUE)
plotQLDisp(fit)

# Smaller prior df estimates indicate that the true unknown dispersions are highly variable, so weaker moderation towards the trend is appropriate.
summary(fit$df.prior)

##### Make contrasts #####

# Ctl_24h vs. Exp_24h 
Ctl24h_Exp24h <- makeContrasts(groupCtl_24h-groupExp_24h, levels=design)
qlf_Ctl24h_Exp24h  <- glmQLFTest(fit, contrast=Ctl24h_Exp24h)
is.de <- decideTestsDGE(qlf_Ctl24h_Exp24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl24h_Exp24h <- topTags(qlf_Ctl24h_Exp24h, n=Inf, adjust.method="BH")
keep_Ctl24h_Exp24h <- out_Ctl24h_Exp24h$table$FDR <= 0.05
output_Ctl24h_Exp24h <- out_Ctl24h_Exp24h[keep_Ctl24h_Exp24h,]
write.csv(as.data.frame(output_Ctl24h_Exp24h), file="EdgeR_glmQLF_res_Ctl24h_Exp24h_HeatShock.csv") #HeatShock


#########################################################################
##### Creating log-fold change against log-counts per million plot ######
#########################################################################

#Plot log-fold change against log-counts per million, with DE genes highlighted:
dev.off()
plotMD(qlf_Ctl24h_Exp24h, p.value = 0.05, main = "")


#############################################
##### To perform likelihood ratio tests #####
#############################################

fit_2 <- glmFit(x, design, robust=TRUE)
lrt_2 <- glmLRT(fit_2)


# Ctl_24h vs. Exp_24h 
Ctl24h_Exp24h <- makeContrasts(groupCtl_24h-groupExp_24h, levels=design)
lrt_Ctl24h_Exp24h  <- glmLRT(fit, contrast=Ctl24h_Exp24h)
is.de <- decideTestsDGE(lrt_Ctl24h_Exp24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl24h_Exp24h <- topTags(lrt_Ctl24h_Exp24h, n=Inf, adjust.method="BH")
keep_Ctl24h_Exp24h <- out_Ctl24h_Exp24h$table$FDR <= 0.05
output_Ctl24h_Exp24h <- out_Ctl24h_Exp24h[keep_Ctl24h_Exp24h,]
write.csv(as.data.frame(output_Ctl24h_Exp24h), file="EdgeR_glmLRT_res_Ctl24h_Exp24h_HeatShock.csv")


#########
# Limma #
#########

#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

dge <- DGEList(counts=countData)

design <- model.matrix(~ 0 + group, data=colData)

dge <- calcNormFactors(dge)

v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
head(coef(fit))

### Ctl_24h vs. Exp_24h 
contr <- makeContrasts(groupCtl_24h - groupExp_24h, levels = colnames(coef(fit)))
# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
# Empirical Bayes smoothing of standard errors
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
# How many DE genes are there?
length(which(top.table$adj.P.Val < 0.05))
summary(decideTests(tmp))
# Write top.table to a file
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.csv(top.table, file = "Limma_res_Ctl24h_Exp24h_HeatShock.csv")
