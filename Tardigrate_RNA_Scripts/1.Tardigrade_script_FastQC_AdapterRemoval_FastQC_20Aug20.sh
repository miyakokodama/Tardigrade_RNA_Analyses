module load java/v1.8.0_202-jdk
#module load python/v3.5.2
module load python/v2.7.17
module load cutadapt/v2.6
module load fastqc/v0.11.8
module load bowtie2/v2.3.4.3
module load htslib/v1.9
module load samtools/v1.9
module load AdapterRemoval/v2.3.1
module load MultiQC/v1.8

###################
# Analyses
###################

#ref_genome=/groups/hologenomics/miyako/data
raw_data=/groups/hologenomics/miyako/data/Tardigrade_BGI/CLEAN
out_direc=/groups/hologenomics/miyako/data/Tardigrade

cat Tardigrade_Sample_List | while read LINES
do
	PrefixName=`echo $LINES`
	echo $PrefixName

# run FastQC on the CLEAN data
	mkdir -p $out_direc/fastqc_before/$PrefixName
	cd $raw_data/$PrefixName
	fastqc -t 2 $PrefixName"_1".fq.gz $PrefixName"_2".fq.gz -o $out_direc/2_fastqc_before/$PrefixName

# run AdapterRemoval
	mkdir -p $out_direc/adapterremoval/$PrefixName
	cd $out_direc/adapterremoval/$PrefixName
	AdapterRemoval --threads 2 --file1 $raw_data/$PrefixName/*_1.fq.gz --file2 $raw_data/$PrefixName/*_2.fq.gz --basename $PrefixName --trimns --trimqualities --minlength 25 --gzip --settings $PrefixName.settings

# run FastQC on the CLEAN data
	mkdir -p $out_direc/fastqc_after/$PrefixName
	cd $out_direc/adapterremoval/$PrefixName
	file=`ls -a $PrefixName.pair*.gz | sort`
	fastqc $file -t 2 -o $out_direc/4_fastqc_after/$PrefixName

done

# run MultiQC on them
	#multiqc $out_direc/2_fastqc_before/*/* --filename Tardigrade_BGI_20samples_multiqc_2_fastqc_before

# run MultiQC on them
	#multiqc $out_direc/4_fastqc_after/*/* --filename Tardigrade_BGI_20samples_multiqc_4_fastqc_after
