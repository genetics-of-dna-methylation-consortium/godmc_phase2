# make grm txt file for glint and convert pheno and meth to text files
# to go in godmc2 pipeline

# Args
arguments <- commandArgs(T)

grm_name <- arguments[1] # [base name of input files]
beta_file <- arguments[2] # methylation Robj - named norm.beta
phen_name <- arguments[3] # ADHD
covs_file <- arguments[4] # methylation data already adjusted for covs
output_path <- arguments[5] 
output_extension <- arguments[6]
home_dir<-arguments[7]


library(genio)


# impute.matrix borrowed from meffil
impute.matrix <- function(x, margin=1, fun=function(x) mean(x, na.rm=T)) {
  if (margin == 2) x <- t(x)
  
  idx <- which(is.na(x) | !is.finite(x), arr.ind=T)
  if (length(idx) > 0) {
    na.idx <- unique(idx[,"row"])
    v <- apply(x[na.idx,,drop=F],1,fun) ## v = summary for each row
    v[which(is.na(v))] <- fun(v)      ## if v[i] is NA, v[i] = fun(v)
    x[idx] <- v[match(idx[,"row"],na.idx)] ##
    stopifnot(all(!is.na(x)))
  }
  
  if (margin == 2) x <- t(x)
  x
}

# read in GRM
grm <- read_grm(
  grm_name,
  n_ind = NA,
  verbose = TRUE,
  ext = "grm",
  shape = "triangle",
  size_bytes = 4,
  comment = "#"
)

grm_mat <- as.data.frame(grm$kinship)

# load in DNAm, pheno, cellcount data
# in beta_file cols should be samples and rows should be cpgs
load(beta_file)
# impute missing DNAm values with row means
# as glint can't handle missing values
norm.beta <- impute.matrix(norm.beta,1)

# organise pheno data
phen_file <- paste(home_dir,"/processed_data/genetic_data/PRS_",phen_name,".sscore",sep="")
message(paste("Loading PRS data for",phen_name))

pheno <- read.table(phen_file, header=T, stringsAsFactors=FALSE)
rownames(pheno) <- pheno$IID
pheno<-pheno[,which(names(pheno)%in%c("IID","SCORE"))]
pheno <- subset(pheno, IID %in% colnames(norm.beta), select="SCORE")

if(sum(!is.na(pheno[["SCORE"]])) < 10)
{
  message("There are fewer than 10 individuals remaining. Stopping the analysis.")
  q()
}
pheno<-na.omit(pheno)
lab_name<-which(names(pheno)%in%c("SCORE"))
names(pheno)[lab_name]<-paste0("PRS_",phen_name)

# load covariates
message("Loading covariates")
covs<-read.table(covs_file,sep=" ",header=T)
rownames(covs) <- covs$IID 
w<-which(names(covs)%in%c("IID"))
covs<-covs[,-w]

# match up DNAm and grm IDs and make sure they're in the same order
participants <- as.character(intersect(colnames(norm.beta),colnames(grm_mat)))
norm.beta <- norm.beta[,participants]
grm_mat <- grm_mat[,participants]
grm_mat <- grm_mat[participants,]
stopifnot(identical(colnames(norm.beta),colnames(grm_mat)))

idx<-match(participants,row.names(pheno))
pheno<-pheno[idx, , drop = FALSE]
idx<-match(participants,row.names(covs))
covs<-covs[idx, , drop = FALSE]

stopifnot(identical(rownames(pheno),rownames(covs)))
stopifnot(identical(rownames(pheno),colnames(norm.beta)))

# save out DNAm
write.table(norm.beta, file=paste0(output_path,"/dnam_for_glint_",output_extension,".txt"),na = "NaN", sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)

# save out grm removing col and row names 
write.table(grm_mat, file=paste0(output_path,"/grm_for_glint_",output_extension,".txt"),sep = "\t", quote = FALSE, col.names = F, row.names = F)

# Write out
write.table(
  pheno,
  file = paste0(output_path,"/phenotypes_for_glint_",output_extension,".txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

if(length(which(names(covs)%in%c("Sex_factor")))>0){
covs$Sex_factor <- ifelse(covs$Sex_factor == "M", 1,
                                 ifelse(covs$Sex_factor == "F", 0, NA))}

# Write out
write.table(
  covs,
  file = paste0(output_path,"/covariates_for_glint_",output_extension,".txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

