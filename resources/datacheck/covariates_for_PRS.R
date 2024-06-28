errorlist <- list()
warninglist <- list()

library(data.table)
suppressMessages(library(matrixStats))

args <- (commandArgs(TRUE));
covariates_file <- as.character(args[1]);
fam_file <- as.character(args[2]);
meth_ids_file <- as.character(args[3])
cohort_descriptives_file <- as.character(args[4])

if(covariates_file == "NULL") 
{
	msg <- "No covariates present"
	warninglist <- c(warninglist, msg)
	message(msg)
	q()
}

message("Checking covariates file: ", covariates_file)
covar <- read.table(covariates_file,header=T)
cov1 <- dim(covar)[1]
cov2 <- dim(covar)[2]

meth_ids <- scan(meth_ids_file, what="character")
fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)

commonids_mgc <- Reduce(intersect, list(meth_ids, covar$IID, fam[,2]))
message("Number of samples with covariate, methylation and genetic data: ", length(commonids_mgc))

if(length(commonids_mgc) < 200)
{
	msg <- paste0("must have at least 200 individuals with covariate, methylation and genetic data.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}


w <- which(names(covar)[1] %in% c("IID"))
if(w!=1)
{
	msg <- paste0("first column from cellcounts file should be the sample identifier with the name IID")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

g1<-grep("_factor",names(covar))
g2<-grep("_numeric",names(covar))
g<-unique(c(g1,g2))

if(length(g)!=(cov2-1))
{
	msg <- paste0("have you specified whether your covariates are factors or numeric in the header of the covariates file? Please make sure your column headers are e.g. 'IID Sex_factor Age_numeric' etc")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

for (i in 1:length(g1))
{
	if(length(table(na.omit(covar[,g1[i]])))==(cov1))
	{
		msg <- paste0(g1[i], " is specified as a factor but has the same number of levels as individuals")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}
}

a <- apply(covar,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*length(meth_ids))))
{
	msg <- paste0("more than 10% of missingness in the following covariates:\n", 
		paste(names(covar)[a>0.1*length(meth_ids)], collapse="\n")
	)
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

covar <- subset(covar, IID %in% commonids_mgc)

cohort_summary <- list()

for (var in names(covar[g1]))
{
  cohort_summary[[var]] <- table(covar[[var]])
}

for (var in names(covar[g2]))
{
  sum <- as.list(summary(covar[[var]]))
  sum$sd <- sd(covar[[var]])
  cohort_summary[[var]] <- as.data.frame(sum)
}

save(cohort_summary, file=cohort_descriptives_file)


message("\n\nCompleted checks\n")

message("Summary of data:\n")
for(i in 1:length(cohort_summary))
{
	a <- cohort_summary[[i]]
	if(is.numeric(a)) a <- round(a, 2)
	message(names(cohort_summary)[i], ": ", paste(a, collapse=", "))
}


if(length(warninglist) > 0)
{
	message("\n\nPlease take note of the following warnings, and fix and re-run the data check if you see fit:")
	null <- sapply(warninglist, function(x)
	{
		message("- ", x)
	})
}


if(length(errorlist) > 0)
{
	message("\n\nThe following errors were encountered, and must be addressed before continuing:")
	null <- sapply(errorlist, function(x)
	{
		message("- ", x)
	})
#	q(status=1)
}
message("\n\n")
