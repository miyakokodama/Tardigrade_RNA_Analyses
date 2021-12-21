####################################
# Summarise the results from the three DE methods
####################################

# R

# Read files
# Here add the names of the cvs files where you have the results from each of the three methods used:
m1<-read.csv("results_method1.csv")
m2<-read.csv("results_method2.csv")
m3<-read.csv("results_method3.csv")


# read annotation files:
annotation_info1<-read.table("Ro-001-Unigene.fa.blast.swissprot.txt", as.is=T, sep="\t", h=T)
annotation_info2<-read.table("Ro001UnigeneGene2GOedited", as.is=T, sep="\t")

# the 1 inside the brackets correspond to the column that contain the gene id (unigene or contig ID), so if it is not column 1 then change accordingly
DE_genes<-unique(m1[,1], m2[,1], m3[,1])

m1_p<-NULL
m2_p<-NULL
m3_p<-NULL
gene_id<-NULL
GO_id<-NULL
descrip<-NULL
m1_FDR<-NULL
m2_FDR<-NULL
m3_FDR<-NULL
for(i in 1:length(DE_genes)){
	# change here the 1 if that column doesn't corresponds to the gene id
	if(sum(m1[,1]==DE_genes[i])>0){
		# Here you could change "method 1" for the name of the program you used (e.g. limma, DESeq)
		m1_p<-c(m1_p, "Method 1")
		# Here the 1 in "m1[,1]" should be the column that contains the gene id, and the 2 in "[i],2]" should be the column that contains the FDR
		m1_FDR<-c(m1_FDR, m1[m1[,1]==DE_genes[i],2])
	}else{
		m1_p<-c(m1_p, NA)
		m1_FDR<-c(m1_FDR, NA)
	}

	# change here the 1 if that column doesn't corresponds to the gene id
	if(sum(m2[,1]==DE_genes[i])>0){
		# Here you could change "method 2" for the name of the program too
		m2_p<-c(m2_p, "Method 2")
		# Here the 1 in "m2[,1]" should be the column that contains the gene id, and the 2 in "[i],2]" should be the column that contains the FDR
		m2_FDR<-c(m2_FDR, m2[m2[,1]==DE_genes[i],2])

	}else{
		m2_p<-c(m2_p, NA)
		m2_FDR<-c(m2_FDR, NA)
	}

	# change here the 1 if that column doesn't corresponds to the gene id
	if(sum(m3[,1]==DE_genes[i])>0){
		# Here you could change "method 1" for the name of the program too
		m3_p<-c(m3_p, "Method 3")
		# Here the 1 in "m3[,1]" should be the column that contains the gene id, and the 2 in "[i],2]" should be the column that contains the FDR
		m3_FDR<-c(m3_FDR, m3[m3[,1]==DE_genes[i],2])

	}else{
		m3_p<-c(m3_p, NA)
		m3_FDR<-c(m3_FDR, NA)
	}

	if(sum(annotation_info1[,1]==DE_genes[i])>0){
		temp<-annotation_info1[annotation_info1[,1]==DE_genes[i],2]
		gene_id<-c(gene_id, temp[1])
		temp<-annotation_info1[annotation_info1[,1]==DE_genes[i],13]
		descrip<-c(descrip, temp[1])
	}else{
		gene_id<-c(gene_id, NA)
		descrip<-c(descrip, NA)
	}
	if(sum(annotation_info2[,1]==DE_genes[i])>0){
		GO_id<-c(GO_id, annotation_info2[annotation_info2[,1]==DE_genes[i],2])
	}
}

tab<-cbind(DE_genes, m1_p, m1_FDR, m2_p, m2_FDR, m3_p, m3_FDR, gene_id, GO_id, descrip)

colnames(tab)<-c("Seq_id", "m1", "m1_FDR", "m2", "m2_FDR", "m3", "m3_FDR", "gene_id", "GO_kegg", "Description")

write.table(tab, quote=F, sep="\t", col.names=T, row.names=F, file="Summary_table.txt")

####################################
