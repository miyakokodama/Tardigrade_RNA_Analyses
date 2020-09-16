# Tardigrade_RNA_Analyses
This is a collection of scripts that are used to perform differential expression analyses in Tardigrade RNA data

##  Tardigrade_RNA_Analyses
In the Tardigrade_RNA_Analyses depository we have scripts for

• Perfome FastQC on "CLEAN" reads, trimming using AdapterRemoval, FastQC on "CLEAN" and trimmed reads, run MultiQC on before and after trimmed reads

• Run Salmon to estimate gene counts

• Make an input .csv file from individual gene count file

• Preprocess the .csv input file by removing the genes with too many missing counts

• Perform DE analyses using DESeq2

• Perform DE analyses using EdgeR and Limma
