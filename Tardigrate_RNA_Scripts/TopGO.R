############################################
# Enrichment analysis using topGO library
############################################ 

# Prepare annotation data
# save BGI_result/3.Annotation/Ro-001-Unigene.fa.Gene2GO.xls as a csv file
# paste -d "\t" <(cut -f 1 -d ";" Ro001UnigeneGene2GO.csv |tail -n +2) <(cut -f 2- -d ";" Ro001UnigeneGene2GO.csv |tail -n +2|perl -pe 's/;/, /g;' |perl -pe 's/ , ,//g;' |perl -pe 's/, ,//g;' |perl -pe 's/, \n/\n/g;' |perl -pe 's/,\n/\n/g;' ) > Ro001UnigeneGene2GOedited


# Set directory #
setwd("C:/Users/mch497/Dropbox/Carrot/DE/Blast2Go/TopGO")

library("topGO")
library("Rgraphviz")

# read annotation data 
geneID2GO <- readMappings("Ro001UnigeneGene2GOedited", sep = "\t", IDsep = ",")
str(head(geneID2GO))

# Predefine a list of DE genes #
geneNames <- names(geneID2GO)
head(geneNames)
genesOfInterest<-readLines("ListofInterestingGenes.txt")

geneList <- factor(as.integer(geneNames %in% genesOfInterest)) #binary list indicating the DE genes
names(geneList) <- geneNames
str(geneList)

# create an topGO data object #
myGOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

# run the Fisher's exact tests
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.05:
mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.05

# create and write a table summarising the top 'numsignif' results:
allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)

write.csv(allRes, "TopGO_Results_DE.csv")

# print a graph (to a pdf file) with the top 'numsignif' results:
printGraph(myGOdata, resultTopgo, firstSigNodes = numsignif, fn.prefix = "topGO", useInfo = "all", pdfSW = TRUE)
dev.off()

# check topGO manual

browseVignettes("topGO")



