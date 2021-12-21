# documentation
#https://combine-lab.github.io/salmon/getting_started/

# unload salmon
module load salmon/v1.1.0

# to get help
salmon quant --help-reads

###################
# Build index
###################

# go to the directory containing the reference transcriptome
cd /groups/hologenomics/hmd579/data/BGI_2019/Rv/Transcriptome_Denovo_Report.tar/Transcriptome_Denovo_Report/BGI_result/2.Assembly/Ro-001

# keep only the longest isoform for each of the genes in BGI's transcriptome assembly
Rscript KeepUniqCLs.R Ro-001-Unigene.fa Ro-001-Unigene_uniqCLs.fa

# index the reference, get it ready for salmon
salmon index -t Ro-001-Unigene_uniqCLs.fa -i tardigrade_bgi_de_novo_index

###################
# Run Salmon
###################

# Run salmon
# -l A tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.).
for SAMPLE in T_1 T_2  T_3  T_4  T_5  T_6  T_7  T_8  T_9  T_10  T_11  T_12  T_13  T_14  T_15  T_16  T_17  T_18
do
    cd /groups/hologenomics/hmd579/data/Tardigrade_Analysis/Salmon # make this directory (or use name of your choice)
    salmon quant -i /groups/hologenomics/hmd579/data/BGI_2019/Rv/Transcriptome_Denovo_Report.tar/Transcriptome_Denovo_Report/BGI_result/2.Assembly/Ro-001/tardigrade_bgi_de_novo_index -l A -1 /groups/hologenomics/hmd579/data/Tardigrade_Analysis/adapterremoval/$SAMPLE/*.pair1.truncated.gz -2 /groups/hologenomics/hmd579/data/Tardigrade_Analysis/adapterremoval/$SAMPLE/*.pair2.truncated.gz -p 2 -o $SAMPLE"_quants"
done
