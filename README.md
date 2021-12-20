# Tardigrade_RNA_Analyses
This is a collection of scripts that are used to perform differential expression analyses in Tardigrade RNA data.

##  Tardigrade_RNA_Analyses
In the Tardigrade_RNA_Analyses depository we have scripts to:

• Perfome FastQC on "CLEAN" reads, 1) trimming using AdapterRemoval, 2) running FastQC on "CLEAN" and trimmed reads, 3) running MultiQC on before and after trimmed reads. 

• Run Salmon to estimate gene counts.

• Make an input .csv file from individual gene count file. 

• Preprocess the .csv input file by removing the genes with too many missing counts. 

• Perform DE analyses using DESeq2. 

• Perform DE analyses using EdgeR and Limma. 

• Summarize the results from DESeq2, EdgeR and Limma. 

• Perform the GO enrichment analyses on the differentially expressed genes. 
