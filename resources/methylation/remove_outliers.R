library(parallel)
suppressMessages(library(matrixStats))
suppressMessages(library(ewaff))
library(data.table)

arguments <- commandArgs(T)
methylationfile <- arguments[1]
out_file <- arguments[2]
cohort_descriptives_commonids_file <- arguments[3]
methylation_summary_file <- arguments[4]
commonids <- arguments[5]
covar_file<-arguments[6]
covs_commonids_file<-arguments[7]
fam_file <- arguments[8]
bim_file <- arguments[9]
	
message("Reading methylation data...")
load(methylationfile)
message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")
ids<-read.table(commonids,he=F)
message("Data size: ", nrow(ids), " individuals with qc-ed genetic data")
shared_samples <- intersect(colnames(norm.beta), ids[,1])
norm.beta <- norm.beta[,shared_samples]
message("Data size after removal of samples without genetic data: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")

message("Identifying methylation outliers")

summariseMeth <- function(X)
{
	message("Counting outliers in methylation matrix")
	
	norm.beta.copy<-ewaff.handle.outliers(X, method=c("iqr"), iqr.limit=3)
	outliers <-norm.beta.copy[[2]]
	keep_cpgs <- rownames(outliers)[which(outliers$n > 0.9*ncol(X))]
	norm.beta <- norm.beta.copy[[1]][rownames(norm.beta.copy[[1]]) %in% keep_cpgs,]
	save(norm.beta, file=out_file)
	message("Data size after removal CpGs with >10% missing values across all samples: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")
	
	message("Estimating means")
	means <- rowMeans(norm.beta.copy[[1]], na.rm=T)

	message("Estimating SDs")
	sds <- rowSds(norm.beta.copy[[1]], na.rm=T)

	message("Estimating medians")
	medians <- rowMedians(norm.beta.copy[[1]], na.rm=T)
	
	dat <- data.frame(cpg=rownames(norm.beta.copy[[1]]), mean=means, median=medians, sd=sds, outlier=outliers)
	return(dat)
}

message("Generating summary stats of methylation")

meth_summary <- summariseMeth(norm.beta)
save(meth_summary, file=methylation_summary_file)

message("Reading covariate data and remove outliers")
covs <- read.table(covar_file, header = T, colClasses=c('Sex_factor'='factor'))
m<-match(ids[,1],covs$IID)
covs<-covs[m,]
write.table(covs,covs_commonids_file,sep="\t",quote=F,row.names=F,col.names=T)

bim <- fread(bim_file)
fam <- read.table(fam_file,header=F,stringsAsFactors=F)

cohort_summary <- list()
cohort_summary$mqtl_sample_size <- nrow(ids)
cohort_summary$methylation_sample_size <- ncol(norm.beta)
cohort_summary$n_CpGs <- nrow(norm.beta)
cohort_summary$geno_sample_size <- nrow(fam)
cohort_summary$n_snp <- nrow(bim)
cohort_summary$mqtl_n_males <- sum(covs$Sex_factor == "M",na.rm=T)
cohort_summary$mqtl_n_females <- sum(covs$Sex_factor == "F",na.rm=T)
cohort_summary$mqtl_mean_age <- mean(covs$Age_numeric,na.rm=T)
cohort_summary$mqtl_median_age <- median(covs$Age_numeric,na.rm=T)
cohort_summary$mqtl_sd_age <- sd(covs$Age_numeric,na.rm=T)
cohort_summary$mqtl_max_age <- max(covs$Age_numeric,na.rm=T)
cohort_summary$mqtl_min_age <- min(covs$Age_numeric,na.rm=T)
cohort_summary$covariates <- names(covs)[-1]

save(cohort_summary, file=cohort_descriptives_commonids_file)


