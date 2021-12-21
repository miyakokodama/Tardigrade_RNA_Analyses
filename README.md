# Tardigrade_RNA_Analyses
This is a collection of scripts that are used to perform differential expression analyses in Tardigrade RNA data.

## Pipeline
Script used for the analyses can be found in the the `Tardigrate_RNA_Scripts` folder in this repository.

Here is a brief description of the steps followed:

- Initial quality assessment of the sequencing data and adapter removal (`1.Tardigrade_script_FastQC_AdapterRemoval_FastQC_20Aug20.sh`)
    - Run FastQC on the "CLEAN" reads
    - Trim the sequencing adapter from the reads using AdapterRemoval2.0
    - Run FastQC on "CLEAN" and trimmed reads
    - Run MultiQC on before and after trimmed reads
- Run Salmon to estimate gene counts (`2.Tardigrade_Salmon.sh`)
- Make an input ´CSV´ file from individual gene count file (`3.Make_Genecount_File_From_Salmon_16Sep20.sh`)
- Preprocess the ´CSV´ input file by removing the genes with too many missing counts (`4.Data_Filtering.R`)
- Perform differential expression (DE) analyses using DESeq2 (`5.DESeq2_R_Code_Tardigrade_1Oct20.R`)
- Perform DE analyses using EdgeR and Limma (`6.EdgeR_R_Code_Tardigrade_7Oct20.R`)
- Summarize the results from DESeq2, EdgeR and Limma (`Summarise_DE_Results.R`)
- Perform the GO enrichment analyses on the differentially expressed genes (`TopGO.R`)
