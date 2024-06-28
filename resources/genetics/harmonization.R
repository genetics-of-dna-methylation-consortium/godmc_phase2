library(data.table)

#Thanks to EasyQC harmonization

args <- (commandArgs(TRUE));
bim_file <- as.character(args[1]);
SNPfail.out <- as.character(args[2]);

#bim_file="./processed_data/genetic_data/data.bim"

bim <- as.data.frame(fread(bim_file))

#bim<-data.frame(f$V4,f[,1],f$V2,f$V3,f$V2,f$V3)
bim[,2]<-as.character(bim[,2])
bim[,5]<-as.character(bim[,5])
bim[,6]<-as.character(bim[,6])


harmonize.alleles <- function(bim=bim,SNPfail=SNPfail) {

message("Checking allele coding")
a1<-bim[,5]
a2<-bim[,6]

SNPfail<-data.frame()

### A1=NA & A2=NA
isBothMissing <- which(is.na(a1) & is.na(a2))
if(length(isBothMissing)>0) {SNPfail<-rbind(SNPfail,bim[isBothMissing,2])}

message("Allele harmonization:",length(isBothMissing)," alleles with  NA are going to be removed")
	
## Recode missing single allele or '<DEL>' or '-' to 'D' and set the other to 'I'

    a1<-bim[,5]
	a2<-bim[,6]
    isRecodea1 <- which(is.na(a1)|a1=='<DEL>'|a1=='-')
	#bim[isRecodea1,5] <- "D"
	#bim[isRecodea1,6] <- "I"
	if(length(isRecodea1)>0) {SNPfail<-rbind(SNPfail,bim[isRecodea1,2])}
	message("Allele harmonization:",length(isRecodea1)," SNPs with A1 alleles with  NA or <DEL> or - are removed")

    a1<-bim[,5]
	a2<-bim[,6]
	isRecodea2 <- which(is.na(a2)|a2=='<DEL>'|a2=='-')
	#bim[isRecodea2,5] <- "I"
	#bim[isRecodea2,6] <- "D"
	if(length(isRecodea2)>0) {SNPfail<-rbind(SNPfail,bim[isRecodea2,2])}
	message("Allele harmonization:",length(isRecodea2)," SNPs with A2 alleles with  NA or <DEL> or - are removed")
	
## Recode MACH R/I -> D/I and R/D -> I/D
	a1<-bim[,5]
	a2<-bim[,6]
    isRecodea1_mach1 <- which(a1=="R" & a2=="I")
	#bim[isRecodea1_mach1,5] <- "D"
	if(length(isRecodea1_mach1)>0) {SNPfail<-rbind(SNPfail,bim[isRecodea1_mach1 ,2])}
	message("Allele harmonization:",length(isRecodea1_mach1)," SNPs with A1 alleles with R/I are removed")
	
	a1<-bim[,5]
	a2<-bim[,6]
	isRecodea1_mach2 <- which(a1=="R" & a2=="D")
	#bim[isRecodea1_mach2,5] <- "I"
	if(length(isRecodea1_mach2)>0) {SNPfail<-rbind(SNPfail,bim[isRecodea1_mach2 ,2])}
	
	message("Allele harmonization:",length(isRecodea1_mach2)," SNPs with A1 alleles with R/D are removed")
	
	a1<-bim[,5]
	a2<-bim[,6]
	isRecodea2_mach1 <- which(a2=="R" & a1=="I")
	#bim[isRecodea2_mach1,6] <- "D"
	if(length(isRecodea2_mach1)>0) {SNPfail<-rbind(SNPfail,bim[isRecodea2_mach1 ,2])}
	
	message("Allele harmonization:",length(isRecodea2_mach1)," SNPs with A2 alleles with I/R are removed")
	
	a1<-bim[,5]
    a2<-bim[,6]
	isRecodea2_mach2 <- which(a2=="R" & a1=="D")
	#bim[isRecodea2_mach2,6] <- "I"
	if(length(isRecodea2_mach2)>0) {SNPfail<-rbind(SNPfail,bim[isRecodea2_mach2 ,2])}
	message("Allele harmonization:",length(isRecodea2_mach2)," SNPs with A2 alleles with D/R are removed")
	
	a1<-bim[,5]
	a2<-bim[,6]
	isRecodea1_Y_Z <- which(a1=="Y" & a2=="Z")
	if(length(isRecodea1_Y_Z)>0) {SNPfail<-rbind(SNPfail,bim[isRecodea1_Y_Z ,2])}
	message("Allele harmonization:",length(isRecodea1_Y_Z)," SNPs with A1 alleles with Y/Z are removed")
	
	a1<-bim[,5]
    a2<-bim[,6]
	isRecodea2_Y_Z <- which(a2=="Y" & a1=="Z")
	if(length(isRecodea2_Y_Z)>0) {SNPfail<-rbind(SNPfail,bim[isRecodea2_Y_Z ,2])}
	message("Allele harmonization:",length(isRecodea2_Y_Z)," SNPs with A2 alleles with D/R are removed")

## Recode Sequence coding to D/I
	a1<-bim[,5]
	a2<-bim[,6]
	isRecode_seq1 <- which(nchar(a1)>nchar(a2))
	#bim[isRecode_seq1,5] <- "I"
	#bim[isRecode_seq1,6] <- "D"
	if(length(isRecode_seq1)>0) {SNPfail<-rbind(SNPfail,bim[isRecode_seq1 ,2])}
	message("Allele harmonization:",length(isRecode_seq1)," alleles with sequence coding are removed")
	
	a1<-bim[,5]
	a2<-bim[,6]
	isRecode_seq2 <- which(nchar(a1)<nchar(a2))
	#bim[isRecode_seq2,5] <- "D"
	#bim[isRecode_seq2,6] <- "I"
	if(length(isRecode_seq2)>0) {SNPfail<-rbind(SNPfail,bim[isRecode_seq2 ,2])}
	message("Allele harmonization:",length(isRecode_seq2)," alleles with sequence coding are removed")

	a1<-bim[,5]
	a2<-bim[,6]
	isRecode_seq3 <- which(nchar(a1)>1&nchar(a1)==nchar(a2))
	if(length(isRecode_seq3)>0) {SNPfail<-rbind(SNPfail,bim[isRecode_seq3 ,2])}
	message("Allele harmonization:",length(isRecode_seq3)," alleles with sequence coding are removed")

	a1<-bim[,5]
    a2<-bim[,6]
    #isInvalid <- !(a1%in%c("A","C","G","T","I","D")&a2%in%c("A","C","G","T","I","D")&a1!=a2)
	isInvalid <- !(a1%in%c("A","C","G","T")&a2%in%c("A","C","G","T")&a1!=a2)
	if(length(which(isInvalid))>0) {
	SNPfail<-rbind(SNPfail,bim[which(isInvalid),2])}
    message("Allele harmonization:",length(which(isInvalid))," alleles with coding other than A,C,T,G are going to be removed")
	

	rm(a1,a2)
	rm(isRecode_seq1,isRecode_seq2)
	rm(isRecodea1_mach1,isRecodea1_mach2,isRecodea2_mach1,isRecodea2_mach2)
	rm(isInvalid)
	
	recoded.bim <- list(bim, SNPfail)
	
	return(recoded.bim)
}

recoded.bim<-harmonize.alleles(bim,SNPfail)
if(length(recoded.bim[[2]])>0){
SNPfailures<-as.character(t(recoded.bim[[2]]))
write.table(SNPfailures,SNPfail.out,sep="\t",quote=F,col.names=F,row.names=F)}
#write.table(recoded.bim[[1]],bim_file,sep="\t",quote=F,col.names=F,row.names=F)



