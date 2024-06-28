errorlist <- list()
warninglist <- list()

library(data.table)
suppressMessages(library(matrixStats))

args <- (commandArgs(TRUE));
phenotypes_for_PRS_file <- as.character(args[1]);
fam_file <- as.character(args[2]);
cohort_descriptives_file <- as.character(args[3])

message("Checking phenotype for PRS file: ", phenotypes_for_PRS_file)

if(phenotypes_for_PRS_file == "NULL") 
{
	msg <- paste0("No phenotypes for PRS have been provided.")
	warninglist <- c(warninglist, msg)
	message( msg)
	q()
}

phenoPRS <- read.table(phenotypes_for_PRS_file,header=T)
phenoPRS1 <- dim(phenoPRS)[1]
phenoPRS2 <- dim(phenoPRS)[2]

fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)

commonids_mgc <- Reduce(intersect, list(phenoPRS$IID, fam[,2]))
message("Number of samples with pheno_for_PRS and genetic data: ", length(commonids_mgc))

w <- which(names(phenoPRS)[1] %in% c("IID"))
if(w!=1)
{
	msg <- paste0("first column from phenotypes for PRS file should be the sample identifier with the name IID")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

g1<-grep("_factor",names(phenoPRS))
g2<-grep("_numeric",names(phenoPRS))
g<-unique(c(g1,g2))

if(length(g)!=(phenoPRS2-1))
{
	msg <- paste0("have you specified whether your phenotypes are factors or numeric in the header of the phenotypes for PRS file? Please make sure your column headers are e.g. 'IID ADHD_factor ADHD_numeric' etc")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}


for (i in 1:length(g1))
{
	if(length(table(na.omit(phenoPRS[,g1[i]])))==(phenoPRS1))
	{
		msg <- paste0(g1[i], " is specified as a factor but has the same number of levels as individuals")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}
}


phenoPRS <- subset(phenoPRS, IID %in% commonids_mgc)

cohort_summary <- list()

for (var in names(phenoPRS[g1]))
{
  cohort_summary[[var]] <- table(phenoPRS[[var]])
}

for (var in names(phenoPRS[g2]))
{
  sum <- as.list(summary(phenoPRS[[var]]))
  sum$sd <- sd(phenoPRS[[var]])
  cohort_summary[[var]] <- as.data.frame(sum)
}

save(cohort_summary, file=cohort_descriptives_file)

message("\n\nCompleted checks\n")

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
