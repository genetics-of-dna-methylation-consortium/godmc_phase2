args <- (commandArgs(TRUE));
frq_file <- as.character(args[1]);
out_file <- as.character(args[2]);
easyQC_file <- as.character(args[3]);
easyQC_script <- as.character(args[4]);

library(data.table)
frq<-as.data.frame(fread(frq_file,header=T))
message("read frq file")

spl<-strsplit(frq$ID,split=":")
spl<-do.call("rbind",spl)
spl2<-strsplit(spl[,2],split="_")
spl2<-do.call("rbind",spl2)
out<-data.frame(CHR=spl[,1],SNP=frq$ID, POS=spl2[,1],EFFECT_ALLELE=frq$ALT,OTHER_ALLELE=frq$REF,EAF=frq$ALT_FREQS,N=frq$OBS_CT/2,BETA=0,SE=0,PVAL=1,IMPUTATION=1,STRAND="+")

write.table(out,paste(out_file,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
message("created easyQC input file")

library(EasyQC)
scripts_dir<-getwd()
setwd(dirname(frq_file))
EasyQC(easyQC_script)
setwd(scripts_dir)
message("easyQC strand check finished successfully")
message("read easyQC_file")
r<-read.table(easyQC_file,header=T,sep="\t")
cat(easyQC_file,"\n")
head(easyQC_file)
message("read easyQC file successfully")

if(file.exists(paste(out_file,".mismatch.txt",sep=""))){
mm<-read.table(paste(out_file,".mismatch.txt",sep=""),header=T)
message("read mismatch file successfully")
}

if(!file.exists(paste(out_file,".mismatch.txt",sep=""))){
mm<-data.frame()
mm$SNP<-NULL
message("WARNING: couldn't find file data.easyqc.mismatch.txt")
}

if(file.exists(paste(out_file,".notinref.txt",sep=""))){
notinref<-read.table(paste(out_file,".notinref.txt",sep=""),header=T)
message("read not in ref file successfully")
}

if(!file.exists(paste(out_file,".notinref.txt",sep=""))){
not_in_ref<-data.frame()
not_in_ref$SNP<-NULL
message("WARNING: couldn't find file data.easyqc.notinref.txt")
}

if(file.exists(paste(out_file,".AFCHECK.outlier.txt",sep=""))){
afoutlier<-read.table(paste(out_file,".AFCHECK.outlier.txt",sep=""),sep="\t",header=T)
message("read AF outlier file successfully")
}

if(!file.exists(paste(out_file,".AFCHECK.outlier.txt",sep=""))){
afoutlier<-data.frame()
afoutlier$SNP<-NULL
message("WARNING: couldn't find file data.easyqc.AFCHECK.outlier.txt")
}

#strand.rm<-unique(c(as.character(mm$SNP),as.character(notinref$SNP),as.character(afoutlier$SNP)))
#w<-which(bim[,2]%in%strand.rm) #4494

SNP.rm<-unique(c(as.character(mm$SNP),as.character(afoutlier$SNP)))
write.table(SNP.rm,paste(out_file,".mismatch_afcheck.failed.SNPs.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
message("extracted mismatched and misaligned SNPs")
