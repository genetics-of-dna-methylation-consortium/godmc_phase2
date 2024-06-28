# Create covariates file

# Original covariates
# Genetic PCs
# Cell counts


arguments <- commandArgs(T)

cov_file <- arguments[1]
pca_file <- arguments[2]
cellcount_file <- arguments[3]
smoking_file <- arguments[4]
fam_file <- arguments[5]
out_file <- arguments[6]

pca <- read.table(pca_file)[,-1]
names(pca) <- c("IID", paste("genetic_pc", 1:(ncol(pca)-1), sep=""))

if(cov_file=="NULL")
{
	covs <- data.frame(IID=pca$IID)
} else {
	covs <- read.table(cov_file, header=TRUE, stringsAsFactors=TRUE, colClass=c("Sex_factor"="character"))
}

# Remove FID column if it is present
if("FID" %in% names(covs)) covs <- subset(covs, select=-c(FID))

smok <- read.table(smoking_file, he=T)[,c('IID', 'Smoking')]

if(cellcount_file=="NULL")
{
  allcovs <- merge(pca, covs, by="IID", all=TRUE)
  allcovs <- merge(allcovs, smok, by="IID", all=TRUE) 

} else {
  cellcount <- read.table(cellcount_file, header=TRUE)

#cellcount <- cellcount[,-ncol(cellcount)]

covs <- subset(covs, IID %in% cellcount$IID)

allcovs <- merge(cellcount, pca, by="IID", all=TRUE)
allcovs <- merge(allcovs, covs, by="IID", all=TRUE)
allcovs <- merge(allcovs, smok, by="IID", all=TRUE)
}

fam<-read.table(fam_file)
m<-match(fam[,2],allcovs[,1])
allcovs<-data.frame(FID=fam[,2],allcovs[m,])
write.table(allcovs,file=paste(out_file,".plink",sep=""),row.names=F,col.names=F,quote=F)
write.table(allcovs[,-1], file=paste0(out_file, ".txt"), row=F, col=T, qu=F)

mat <- t(as.matrix(allcovs[,-1]))
write.table(mat, file=paste0(out_file, ".matrixeqtl"), row=TRUE, col=FALSE, qu=FALSE, sep="\t")

