# R
# set a working directory
setwd("/groups/hologenomics/hmd579/data/Tardigrade_Analysis/Salmon/")

# remove anything before we start
rm(list = ls())

#################################################################
### load count data or depth data (with or without row.names) ###
#################################################################

# load a count data
# genes in a row, samples in a column
countData <- read.csv("genecount.csv", header=TRUE, row.names=1)

### remove rows with  ###
# https://stackoverflow.com/questions/32618583/counting-number-of-instances-of-a-condition-per-row-r

countData$missing_number <- rowSums(countData == 0)
countData_missing_removed <- subset(countData, missing_number <= 9) # genes typed at least by 50% of the samples; relax the threshold if needed

# delete the column
countData_missing_removed$missing_number <- NULL

# change the counts to numeric; maybe this needs to be changed to integers
all(is.numeric(countData_missing_removed))

write.csv(countData_missing_removed, file=".csv", row.names = TRUE)
