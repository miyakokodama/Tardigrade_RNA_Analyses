###############################################################
# Filter Trinity assembly and keep longest isoform from the CLs 
###############################################################


library(ShortRead)

# Assembly with the longest CLs and Unigenes
faFile<-readFasta("Ro-001-Unigene.fa")

faUniGene<-faFile[grep("Unigene", as.character(id(faFile)))]
faCL<-faFile[grep("CL", as.character(id(faFile)))]

CLnames<-sapply(strsplit(as.character(id(faCL)), "\\."), "[[", 1)
uCLnames<-unique(CLnames)

keepseqs<-NULL
keepids<-NULL
for(i in 1:length(uCLnames)){
	temp<-faCL[CLnames==uCLnames[i]]
	keepseqs<-c(keepseqs, as.character(sread(temp[which.max(width(temp))])))
	keepids<-c(keepids, as.character(id(temp[which.max(width(temp))])))
}

keepseqs<-c(keepseqs, as.character(sread(faUniGene)))
keepids<-c(keepids, as.character(id(faUniGene)))

newfa<-ShortRead(sread=DNAStringSet(keepseqs), id=BStringSet(keepids))

writeFasta(newfa, "Ro-001-Unigene_uniqCLs.fa")

###############################################################