rm(list = ls())

###
.libPaths()
.libPaths(.libPaths()[2:2])
.libPaths()
###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("DESeq2")
#BiocManager::install("DEFormats")

library(DESeq2)
library(data.table)
library(edgeR)
library(limma)
library(xlsx)
library(statmod)
library(DEFormats)

############################

setwd("\\\\a00519.science.domain/mch497/Documents/Tardigrade")

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
plotMDS(x, col = group)

# Design matrix
design <- model.matrix(~ 0 + group, data=colData)
design
colnames(design) # to get the group names to use for contrasts

# Estimate the dispersions
x <- estimateDisp(x, design, robust=TRUE)

#Scatterplot of the biological coefficient of variation (BCV) against the average abundance of each gene.
plotBCV(x)

################################################
##### To perform quasi-likelihood F-tests ######
################################################

fit <- glmQLFit(x, design, robust=TRUE)
plotQLDisp(fit)

# Smaller prior df estimates indicate that the true unknown dispersions are highly variable, so weaker moderation towards the trend is appropriate.
summary(fit$df.prior)

##### Make contrasts #####

# Ctl_2h vs. Ctl_24h - 0 DE genes
Ctl_2h_Ctl_24h <- makeContrasts(groupCtl_2h-groupCtl_24h, levels=design)
qlf_Ctl_2h_Ctl_24h  <- glmQLFTest(fit, contrast=Ctl_2h_Ctl_24h)
is.de <- decideTestsDGE(qlf_Ctl_2h_Ctl_24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl_2h_Ctl_24h <- topTags(qlf_Ctl_2h_Ctl_24h, n=Inf, adjust.method="BH")
keep_Ctl_2h_Ctl_24h <- out_Ctl_2h_Ctl_24h$table$FDR <= 0.05
output_Ctl_2h_Ctl_24h <- out_Ctl_2h_Ctl_24h[keep_Ctl_2h_Ctl_24h,]
write.csv(as.data.frame(output_Ctl_2h_Ctl_24h), file="EdgeR_glmQLF_res_Ctl2h_Ctl24h.csv")

# Exp_2h vs. Exp_24h - 0 DE genes
Exp2h_Exp24h <- makeContrasts(groupExp_2h-groupExp_24h, levels=design)
qlf_Exp2h_Exp24h  <- glmQLFTest(fit, contrast=Exp2h_Exp24h)
is.de <- decideTestsDGE(qlf_Exp2h_Exp24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Exp2h_Exp24h <- topTags(qlf_Exp2h_Exp24h, n=Inf, adjust.method="BH")
keep_Exp2h_Exp24h <- out_Exp2h_Exp24h$table$FDR <= 0.05
output_Exp2h_Exp24h <- out_Exp2h_Exp24h[keep_Exp2h_Exp24h,]
write.csv(as.data.frame(output_Exp2h_Exp24h), file="EdgeR_glmQLF_res_Exp2h_Exp24h.csv")

# Ctl_2h vs. Exp_2h - 0 DE genes
Ctl2h_Exp2h <- makeContrasts(groupCtl_2h-groupExp_2h, levels=design)
qlf_Ctl2h_Exp2h  <- glmQLFTest(fit, contrast=Ctl2h_Exp2h)
is.de <- decideTestsDGE(qlf_Ctl2h_Exp2h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl2h_Exp2h <- topTags(qlf_Ctl2h_Exp2h, n=Inf, adjust.method="BH")
keep_Ctl2h_Exp2h <- out_Ctl2h_Exp2h$table$FDR <= 0.05
output_Ctl2h_Exp2h <- out_Ctl2h_Exp2h[keep_Ctl2h_Exp2h,]
write.csv(as.data.frame(output_Ctl2h_Exp2h), file="EdgeR_glmQLF_res_Ctl2h_Exp2h.csv")

# Ctl_24h vs. Exp_24h - 3 DE genes
Ctl24h_Exp24h <- makeContrasts(groupCtl_24h-groupExp_24h, levels=design)
qlf_Ctl24h_Exp24h  <- glmQLFTest(fit, contrast=Ctl24h_Exp24h)
is.de <- decideTestsDGE(qlf_Ctl24h_Exp24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl24h_Exp24h <- topTags(qlf_Ctl24h_Exp24h, n=Inf, adjust.method="BH")
keep_Ctl24h_Exp24h <- out_Ctl24h_Exp24h$table$FDR <= 0.05
output_Ctl24h_Exp24h <- out_Ctl24h_Exp24h[keep_Ctl24h_Exp24h,]
write.csv(as.data.frame(output_Ctl24h_Exp24h), file="EdgeR_glmQLF_res_Ctl24h_Exp24h.csv")

#############################################
##### To perform likelihood ratio tests #####
#############################################

fit_2 <- glmFit(x, design, robust=TRUE)
lrt_2 <- glmLRT(fit_2)

# Ctl_2h vs. Ctl_24h - 56 DE genes
Ctl_2h_Ctl_24h <- makeContrasts(groupCtl_2h-groupCtl_24h, levels=design)
lrt_Ctl_2h_Ctl_24h  <- glmLRT(fit, contrast=Ctl_2h_Ctl_24h)
is.de <- decideTestsDGE(lrt_Ctl_2h_Ctl_24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl_2h_Ctl_24h <- topTags(lrt_Ctl_2h_Ctl_24h, n=Inf, adjust.method="BH")
keep_Ctl_2h_Ctl_24h <- out_Ctl_2h_Ctl_24h$table$FDR <= 0.05
output_Ctl_2h_Ctl_24h <- out_Ctl_2h_Ctl_24h[keep_Ctl_2h_Ctl_24h,]
write.csv(as.data.frame(output_Ctl_2h_Ctl_24h), file="EdgeR_glmLRT_res_Ctl2h_Ctl24h.csv")

# Exp_2h vs. Exp_24h - 150 DE genes
Exp2h_Exp24h <- makeContrasts(groupExp_2h-groupExp_24h, levels=design)
lrt_Exp2h_Exp24h  <- glmLRT(fit, contrast=Exp2h_Exp24h)
is.de <- decideTestsDGE(lrt_Exp2h_Exp24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Exp2h_Exp24h <- topTags(lrt_Exp2h_Exp24h, n=Inf, adjust.method="BH")
keep_Exp2h_Exp24h <- out_Exp2h_Exp24h$table$FDR <= 0.05
output_Exp2h_Exp24h <- out_Exp2h_Exp24h[keep_Exp2h_Exp24h,]
write.csv(as.data.frame(output_Exp2h_Exp24h), file="EdgeR_glmLRT_res_Exp2h_Exp24h.csv")

# Ctl_2h vs. Exp_2h - 34 DE genes
Ctl2h_Exp2h <- makeContrasts(groupCtl_2h-groupExp_2h, levels=design)
lrt_Ctl2h_Exp2h  <- glmLRT(fit, contrast=Ctl2h_Exp2h)
is.de <- decideTestsDGE(lrt_Ctl2h_Exp2h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl2h_Exp2h <- topTags(lrt_Ctl2h_Exp2h, n=Inf, adjust.method="BH")
keep_Ctl2h_Exp2h <- out_Ctl2h_Exp2h$table$FDR <= 0.05
output_Ctl2h_Exp2h <- out_Ctl2h_Exp2h[keep_Ctl2h_Exp2h,]
write.csv(as.data.frame(output_Ctl2h_Exp2h), file="EdgeR_glmLRT_res_Ctl2h_Exp2h.csv")

# Ctl_24h vs. Exp_24h - 358 DE genes
Ctl24h_Exp24h <- makeContrasts(groupCtl_24h-groupExp_24h, levels=design)
lrt_Ctl24h_Exp24h  <- glmLRT(fit, contrast=Ctl24h_Exp24h)
is.de <- decideTestsDGE(lrt_Ctl24h_Exp24h, adjust.method="BH", p.value=0.05)
summary(is.de)
out_Ctl24h_Exp24h <- topTags(lrt_Ctl24h_Exp24h, n=Inf, adjust.method="BH")
keep_Ctl24h_Exp24h <- out_Ctl24h_Exp24h$table$FDR <= 0.05
output_Ctl24h_Exp24h <- out_Ctl24h_Exp24h[keep_Ctl24h_Exp24h,]
write.csv(as.data.frame(output_Ctl24h_Exp24h), file="EdgeR_glmLRT_res_Ctl24h_Exp24h.csv")

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

### Ctl_2h vs. Ctl_24h - 0 DE genes ###
contr <- makeContrasts(groupCtl_2h - groupCtl_24h, levels = colnames(coef(fit)))
# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
# Empirical Bayes smoothing of standard errors
tmp <- eBayes(tmp)
# What genes are most differentially expressed?
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
# How many DE genes are there?
length(which(top.table$adj.P.Val < 0.05))
# Write top.table to a file
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.csv(top.table, file = "Limma_res_Ctl2h_Ctl24h.csv")

### Exp_2h vs. Exp_24h - 0 DE genes ###
contr <- makeContrasts(groupExp_2h - groupExp_24h, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.csv(top.table, file = "Limma_res_Exp2h_Exp24h.csv")

### Ctl_2h vs. Exp_2h - 0 DE genes ###
contr <- makeContrasts(groupCtl_2h - groupExp_2h, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.csv(top.table, file = "Limma_res_Ctl2h_Exp2h.csv")

### Ctl_24h vs. Exp_24h - 2 DE genes ###
contr <- makeContrasts(groupCtl_24h - groupExp_24h, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.csv(top.table, file = "Limma_res_Ctl24h_Exp24h.csv")
