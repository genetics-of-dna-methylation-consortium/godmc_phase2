#suppressPackageStartupMessages(library(meffil))
library(data.table)

arguments <- commandArgs(T)

methylationfile_untransformed <- as.character(arguments[1])
methylationfile_transformed <- as.character(arguments[2])
cov_file <- as.character(arguments[3])
hase_cov <- as.character(arguments[4])
hase_cov_males <- as.character(arguments[5])
hase_cov_females <- as.character(arguments[6])


message("Loading covariate data")
cat(cov_file,"\n")
covdat <- read.table(cov_file, header=TRUE, stringsAsFactors=FALSE)
male_ids <- subset(covdat, Sex_factor == "M")$IID
female_ids<- subset(covdat, Sex_factor == "F")$IID

message("Loading untransformed methylation data")
load(paste0(methylationfile_untransformed, ".RData"))
ids<-colnames(norm.beta)
#norm.beta.t<-data.frame(IID=colnames(norm.beta),t(norm.beta))
#norm.beta.t[1:5,1:5]
norm.beta<-data.frame(cpg=rownames(norm.beta),norm.beta,check.names=F)
message("Writing HASE format (untransformed data): ",ncol(norm.beta)," colums and ",nrow(norm.beta)," rows")
fwrite(norm.beta, file=paste0(methylationfile_untransformed, ".csv"), row=F, col=T, qu=F, sep=",")

message("Loading transformed methylation data")
load(paste0(methylationfile_transformed, ".RData"))
#norm.beta.t<-data.frame(IID=colnames(norm.beta),t(norm.beta))
#norm.beta.t[1:5,1:5]
norm.beta<-data.frame(cpg=rownames(norm.beta),norm.beta,check.names=F)
message("Writing HASE format (transformed data): ",ncol(norm.beta)," colums and ",nrow(norm.beta)," rows")
fwrite(norm.beta, file=paste0(methylationfile_transformed, ".csv"), row=F, col=T, qu=F, sep=",")

c<-data.frame(IID=colnames(norm.beta)[-1],cov=rnorm(n=ncol(norm.beta)-1, mean = 0, sd = 1))
fwrite(c,paste0(hase_cov,"/covariates.txt"),sep="\t",quote=F,col.names=T,row.names=F)

if(sum(ids %in% male_ids > 0)){
message("Loading X chromosome untransformed methylation data for males")
load(paste0(methylationfile_untransformed, ".Male.chrX.RData"))
#norm.beta.t<-data.frame(IID=colnames(norm.beta),t(norm.beta))
#norm.beta.t[1:5,1:5]
norm.beta<-data.frame(cpg=rownames(norm.beta),norm.beta,check.names=F)
message("Writing HASE format (X chromosome untransformed methylation data for males): ",ncol(norm.beta)," colums and ",nrow(norm.beta)," rows")
fwrite(norm.beta, file=paste0(methylationfile_untransformed, ".Male.chrX.csv"), row=F, col=T, qu=F, sep=",")

message("Loading Y chromosome untransformed methylation data for males")
load(paste0(methylationfile_untransformed, ".Male.chrY.RData"))
#norm.beta.t<-data.frame(IID=colnames(norm.beta),t(norm.beta))
#norm.beta.t[1:5,1:5]
norm.beta<-data.frame(cpg=rownames(norm.beta),norm.beta,check.names=F)
message("Writing HASE format (Y chromosome untransformed methylation data for males): ",ncol(norm.beta)," colums and ",nrow(norm.beta)," rows")
fwrite(norm.beta, file=paste0(methylationfile_untransformed, ".Male.chrY.csv"), row=F, col=T, qu=F, sep=",")

message("Loading X chromosome transformed methylation data for males")
load(paste0(methylationfile_transformed, ".Male.chrX.RData"))
#norm.beta.t<-data.frame(IID=colnames(norm.beta),t(norm.beta))
#norm.beta.t[1:5,1:5]
message("Writing HASE format (X chromosome transformed methylation data for males): ",ncol(norm.beta)," colums and ",nrow(norm.beta)," rows")
norm.beta<-data.frame(cpg=rownames(norm.beta),norm.beta,check.names=F)
fwrite(norm.beta, file=paste0(methylationfile_transformed, ".Male.chrX.csv"), row=F, col=T, qu=F, sep=",")

message("Loading Y chromosome transformed methylation data for males")
load(paste0(methylationfile_transformed, ".Male.chrY.RData"))
#norm.beta.t<-data.frame(IID=colnames(norm.beta),t(norm.beta))
#norm.beta.t[1:5,1:5]
message("Writing HASE format (Y chromosome transformed methylation data for males): ",ncol(norm.beta)," colums and ",nrow(norm.beta)," rows")
norm.beta<-data.frame(cpg=rownames(norm.beta),norm.beta,check.names=F)
fwrite(norm.beta, file=paste0(methylationfile_transformed, ".Male.chrY.csv"), row=F, col=T, qu=F, sep=",")

c<-data.frame(IID=colnames(norm.beta)[-1],cov=rnorm(n=ncol(norm.beta)-1, mean = 0, sd = 1))
fwrite(c,paste0(hase_cov_males,"/covariates_males.txt"),sep="\t",quote=F,col.names=T,row.names=F)


}


if(sum(ids %in% female_ids > 0))

{
message("Loading X chromosome untransformed methylation data for females")
load(paste0(methylationfile_untransformed, ".Female.chrX.RData"))
#norm.beta.t<-data.frame(IID=colnames(norm.beta),t(norm.beta))
#norm.beta.t[1:5,1:5]
message("Writing HASE format (X chromosome untransformed methylation data for females): ",ncol(norm.beta)," colums and ",nrow(norm.beta)," rows")
norm.beta<-data.frame(cpg=rownames(norm.beta),norm.beta,check.names=F)
fwrite(norm.beta, file=paste0(methylationfile_untransformed, ".Female.chrX.csv"), row=F, col=T, qu=F, sep=",")

message("Loading X chromosome transformed methylation data for females")
load(paste0(methylationfile_transformed, ".Female.chrX.RData"))
#norm.beta.t<-data.frame(IID=colnames(norm.beta),t(norm.beta))
#norm.beta.t[1:5,1:5]
message("Writing HASE format (X chromosome transformed methylation data for females): ",ncol(norm.beta)," colums and ",nrow(norm.beta)," rows")
norm.beta<-data.frame(cpg=rownames(norm.beta),norm.beta,check.names=F)
fwrite(norm.beta, file=paste0(methylationfile_transformed, ".Female.chrX.csv"), row=F, col=T, qu=F, sep=",")

c<-data.frame(IID=colnames(norm.beta)[-1],cov=rnorm(n=ncol(norm.beta)-1, mean = 0, sd = 1))
fwrite(c,paste0(hase_cov_females,"/covariates_females.txt"),sep="\t",quote=F,col.names=T,row.names=F)

}




