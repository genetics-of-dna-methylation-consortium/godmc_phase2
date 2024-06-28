errorlist <- list()
warninglist <- list()


library(data.table)
suppressMessages(library(matrixStats))
suppressMessages(library(meffil))
suppressMessages(library(ggplot2))
library(ewaff)

args <- (commandArgs(TRUE));
betas_file <- as.character(args[1]);
fam_file <- as.character(args[2]);
sorted_methylation <- as.character(args[3]);
methylation_array <- as.character(args[4]);
cellcounts_file_measured <- as.character(args[5]);
meth_ids_file <- as.character(args[6]);
cohort_descriptives_file <- as.character(args[7]);
methylation_summary_file <- as.character(args[8]);
covar_file <- as.character(args[9]);
sex_pred_plot_path <- as.character(args[10]);
ids <- as.character(args[11]);
ids_plink <- as.character(args[12]);

#predicted_cellcounts_type <- "houseman"
#cellcounts_file_measured <- NULL
#meth_ids_file <- "/user/home/epzjlm/scratch/repo/EPIC_mQTL/processed_data/ids/meth_ids.txt"
#cohort_descriptives_file <- "/user/home/epzjlm/scratch/repo/EPIC_mQTL/processed_data/ids/methylation_descriptives.RData"
#methylation_summary_file <-"/user/home/epzjlm/scratch/repo/EPIC_mQTL/results/01/methylation_summary.RData"
#betas_file<-"/user/home/epzjlm/scratch/repo/EPIC_mQTL/input_data/betas_ARIES_F17_EPIC.Robj"
#fam_file<-"/user/home/epzjlm/scratch/repo/EPIC_mQTL/input_data/combined.fam"
#cellcounts_file<-"/user/home/epzjlm/scratch/repo/EPIC_mQTL/input_data/cellcounts.txt"
#ids<-"/user/home/epzjlm/scratch/repo/EPIC_mQTL/processed_data/ids/intersect_ids.txt"
#ids_plink<-"/user/home/epzjlm/scratch/repo/EPIC_mQTL/processed_data/ids/intersect_ids_plink.txt"

message("Checking methylation data: ", betas_file)
load(betas_file)


fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)

if(! "norm.beta" %in% ls())
{
	msg <- paste0("please save methylation matrix with the object name norm.beta")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("norm.beta object found")
# is methylation data a matrix

if(!is.matrix(norm.beta))
{
	msg <- paste0("please transform methylation norm.beta to a matrix")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
# are CpGs in rows?
d1 <- nrow(norm.beta)
nid_meth <- ncol(norm.beta)
message("Number of individuals with methylation data: ", nid_meth)
message("Number of CpGs: ", d1)
if(d1 < nid_meth)
{
	msg <- paste0("please transpose methylation matrix (CpGs in rows; samples in columns)")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("Data is a correctly oriented matrix")

# are individuals unique
if(sorted_methylation == "no") {
if(nid_meth < 200)

{
	msg <- paste0("Number of individuals with methylation is less than ", nid_meth," please contact developers")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("More than 200 individuals with methylation data")
}

#number of individuals
if(any(duplicated(colnames(norm.beta))))
{
	msg <- paste0("please remove duplicate samples from methylation data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

# check for NAs in beta matrix
if(any(is.na(norm.beta)))
{
	msg <- paste0("please remove NAs from methylation matrix")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("No NAs in data")

# check for negative values in beta matrix
if(any(norm.beta < 0))
#if(any(norm.beta < 0,na.rm=T))
{
	msg <- paste0("please remove negative values from methylation matrix. Are these beta values?")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

# check for values above 1 in beta matrix
#if(any(norm.beta > 1,na.rm=T)) 
if(any(norm.beta > 1)) 
{
	msg <- paste0("please remove values > 1 from methylation matrix. Are these beta values?")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("All values are within 0-1")

if(any(grepl("rs", rownames(norm.beta))))
{
	msg <- paste0("there are SNPs in the methylation data. Please remove all rows with rs IDs.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

feat <- meffil.get.features(methylation_array)
#866553
xy<-which(feat$chromosome%in%c("chrX","chrY"))
probes_xy<-as.character(feat[xy,"name"])
xy_overlap<-intersect(row.names(norm.beta),probes_xy)
n.xy_overlap <- length(xy_overlap)
no.overlap<-0.80*length(probes_xy)
if(n.xy_overlap < no.overlap)
{
	msg <- paste0("fewer than 80% chrx and chry probes. Please include chrx and y probes")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

#
probes<-unique(feat$name)
p.overlap<-intersect(row.names(norm.beta),probes)
no.overlap<-0.80*length(probes)
if(length(p.overlap) < no.overlap)
{
	msg <- paste0("fewer than 80% of ",methylation_array," probes. Please include more probes.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}



#
# Predict sex using DNA-methylation, and remove subjects with sex discrepancies.
#

# In preparation of the sex prediction, do a quick check of the covariate data (full check will be done later).
covar <- read.table(covar_file, header = T)

if(! "Sex_factor" %in% names(covar))
{
  msg <- paste0("There is no Sex_factor variable in the covariate file. Please provide M/F values, even if they are all the same sex.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}else{
  covar <- read.table(covar_file,header=T,colClasses=c('Sex_factor'='factor'))
}

if(any(is.na(covar$Sex_factor)))
{
  msg <- paste0("There are some missing values in the Sex_factor column. Please make sure all individuals have data for this column.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

index <- covar$Sex_factor %in% c("M", "F")
if(any(!index))
{
  msg <- paste0("There are some values in the Sex_factor column that are neither M nor F. Please make sure all individuals have data for this column.")
  errorlist <- c(errorlist, msg)
  warning("ERROR: ", msg)
}

#Check methylation array.
if(methylation_array == "450k"){
  message("According to the config-file, the DNA-methylation data were obtained using the 450K array.")
} else if(methylation_array == "epic"){
  message("According to the config-file, the DNA-methylation data were obtained using the EPIC array.")
} else if (methylation_array == "epic2"){
  message("According to the config-file, the DNA-methylation data were obtained using the EPIC v2 array.")
}

#Select the samples that are present in both the methylation data and the covariate file.
shared_samples <- intersect(colnames(norm.beta), covar$IID)
norm.beta <- norm.beta[,shared_samples]
covar <- covar[match(shared_samples, covar$IID),]

#Sex prediction won't work if the distribution of sex is too skewed. If this is the case (less than 10% for one sex), skip sex prediction.
if((length(covar$Sex_factor[covar$Sex_factor == "M"]) < 0.05 * length(covar$Sex_factor) |
    length(covar$Sex_factor[covar$Sex_factor == "F"]) < 0.05 * length(covar$Sex_factor))){
  
  message("Sexes are extremely skewed in this cohort (<5% males or females). Skipping sex-prediction.")
  
} else{
  
  feat_noNA <- feat[!is.na(feat$chromosome),]
  
  x_probes <- feat_noNA$name[feat_noNA$chromosome == "chrX"]
  y_probes <- feat_noNA$name[feat_noNA$chromosome == "chrY"]
  
  # Select the X and Y chromosome probes.
  beta_x <- norm.beta[rownames(norm.beta) %in% x_probes,]
  beta_y <- norm.beta[rownames(norm.beta) %in% y_probes,]
  
  #Remove NAs (these will throw off the PCA).
  beta_x <- na.omit(beta_x)
  beta_y <- na.omit(beta_y)
  
  message(paste0(nrow(beta_x), " X-chromosome probes were selected."))
  message(paste0(nrow(beta_y), " Y-chromosome probes were selected."))
  
  #If the dataset contains no sex-chromosome probes, skip sex-prediction.
  if(nrow(beta_x) == 0 | nrow(beta_y) == 0){
    message("No sex-chromosomal probes detected for the X and/or Y chromosomes. Skipping sex-prediction.")
    
    
  } else{
    
    message("Predicting sex from DNA-methylation on sex-chromosomes.")
    
    #Perform a PCA on the X and Y chromosomes.
    pca_x <- prcomp(t(beta_x))
    pca_y <- prcomp(t(beta_y))
    
    # Make a dataframe containing the first PC of X and Y methylation per person.
    pred_sex <- data.frame(
      IID = colnames(norm.beta),
      x_PC1 = pca_x$x[,1],
      y_PC1 = pca_y$x[,1],
      assumed_sex = factor(covar$Sex_factor)
    )
    
    #The sign of principal components is arbitrary. To enable the use of PC1 as a cutoff, define it as follows:
    # - Females higher than males for PC1 of the X-chromosome. 
    # - Males higher than females for PC1 of the Y-chromosome.
    if(mean(pred_sex$x_PC1[pred_sex$assumed_sex == "F"], na.rm = T) < mean(pred_sex$x_PC1[pred_sex$assumed_sex == "M"], na.rm = T)){
      pred_sex$x_PC1 <- -(pred_sex$x_PC1)
    }
    
    if(mean(pred_sex$y_PC1[pred_sex$assumed_sex == "F"], na.rm = T) > mean(pred_sex$y_PC1[pred_sex$assumed_sex == "M"], na.rm = T)){
      pred_sex$y_PC1 <- -(pred_sex$y_PC1)
    }
    
    #The first PC of the X/Y methylation should be strongly correlated to sex. If this is not the case, sex may be mislabeled. 
    x_cor <- cor(as.numeric(pred_sex$assumed_sex), pred_sex$x_PC1, use = "complete.obs")
    y_cor <- cor(as.numeric(pred_sex$assumed_sex), pred_sex$y_PC1, use = "complete.obs")
    
    if(abs(x_cor) < 0.7 | abs(y_cor) < 0.7){
      warning("The first PC of X/Y methylation does not correlate sufficiently (R > 0.7) with sex. Sex may be mislabeled!")
      
    }
    
    # Based on the PC of the X and Y-chromosomes, divide the samples into males and females.
    x_cutoff <- mean(c(min(pred_sex$x_PC1), max(pred_sex$x_PC1)))
    y_cutoff <- mean(c(min(pred_sex$y_PC1), max(pred_sex$y_PC1)))
    
    pred_sex$X <- ""
    pred_sex$Y <- ""
    
    pred_sex$X[pred_sex$x_PC1 > x_cutoff] <- "XX"
    pred_sex$X[pred_sex$x_PC1 <= x_cutoff] <- "X"
    
    pred_sex$Y[pred_sex$y_PC1 <= y_cutoff] <- ""
    pred_sex$Y[pred_sex$y_PC1 > y_cutoff] <- "Y"
    
    pred_sex$sex_chromosomes <- paste0(pred_sex$X, pred_sex$Y)
    
    # Add the predicted sex according to the X and Y chromosomes.
    # If the two match, add their shared result as the predicted sex.
    # If the X-prediction says female (marking two X-chromsomes) while the Y-prediction says male (marking a Y-chromosome), the subject is XXY (Klinefelter).
    # If the X-prediction says male (no two X-chromosomes) while the Y-prediction says male (no Y-chromosome), the subject is single X (Turner Syndrome).
    pred_sex$predicted_sex <- factor(NA, levels = c("F", "M", "Klinefelter", "Turner"))
    
    pred_sex$predicted_sex[pred_sex$sex_chromosomes == "XX"] <- "F" # XX = female
    pred_sex$predicted_sex[pred_sex$sex_chromosomes == "XY"] <- "M" # XY = male
    pred_sex$predicted_sex[pred_sex$sex_chromosomes == "XXY"] <- "Klinefelter" # XXY = Klinefelter
    pred_sex$predicted_sex[pred_sex$sex_chromosomes == "X"] <- "Turner" # single X = Turner
    
    table(pred_sex$assumed_sex, pred_sex$predicted_sex, useNA = "always")
    
    # Plot X/Y methylation PCs against assumed and predicted sex.
    plot_theme <- list(
      geom_text(size = 1, hjust = 0.5, vjust = 1.5),
      scale_color_manual(values = c("red", "blue", "darkgreen", "purple")),
      theme_bw(),
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank()),
      labs(x = "Subject ID", y = "PC1")
    )
    
    pdf(file = sex_pred_plot_path, width = 6, height = 6)
    
    #Plot the first PC of the X/Y methylation against the assumed sex supplied in the covariate file.
    #X PC1
    print(ggplot(data = pred_sex, aes(x = IID, y = x_PC1, color = assumed_sex, label = IID)) + 
            geom_hline(yintercept = x_cutoff, linetype = "dashed", alpha = 0.5) + 
            plot_theme + 
            labs(title = "PC1 of X methylation vs assumed sex"))
    
    #Y PC1
    print(ggplot(data = pred_sex, aes(x = IID, y = y_PC1, color = assumed_sex, label = IID)) + 
            geom_hline(yintercept = y_cutoff, linetype = "dashed", alpha = 0.5) + 
            plot_theme + 
            labs(title = "PC1 of Y methylation vs assumed sex"))
    
    
    #Plot the first PC of the X/Y methylation against the methylation-predicted sex.
    #X PC1
    print(ggplot(data = pred_sex, aes(x = IID, y = x_PC1, color = predicted_sex, label = IID)) + 
            geom_hline(yintercept = x_cutoff, linetype = "dashed", alpha = 0.5) + 
            plot_theme + 
            labs(title = "PC1 of X methylation vs predicted sex"))
    
    #Y PC1
    print(ggplot(data = pred_sex, aes(x = IID, y = y_PC1, color = predicted_sex, label = IID)) + 
            geom_hline(yintercept = y_cutoff, linetype = "dashed", alpha = 0.5) + 
            plot_theme + 
            labs(title = "PC1 of Y methylation vs predicted sex"))
    
    dev.off()
    
    # Print how many samples were sex-mismatched.
    total_samples <- pred_sex$IID
    
    matched_idx <- as.character(pred_sex$assumed_sex) == as.character(pred_sex$predicted_sex)
    matched_idx[is.na(matched_idx)] <- FALSE
    matched_samples <- total_samples[matched_idx]
    
    mismatched_samples <- total_samples[!(matched_idx)]
    
    total_length <- length(total_samples)
    matched_length <- length(matched_samples)
    mismatched_length <- length(mismatched_samples)
    
    message(sprintf("For %s out of %s samples, predicted sex matched the sex supplied in the covariate file.", matched_length, total_length))
    message(sprintf("For %s out of %s samples, a sex discrepancy was detected. Removing these samples.", mismatched_length, total_length))
    
    # Remove samples with a sex discrepancy from norm.beta.
    # sex_discrepancies <- mismatched_samples
    # save(sex_discrepancies, file=paste0(out_file))
    idx <- colnames(norm.beta) %in% mismatched_samples
    norm.beta <- norm.beta[, !idx]
    
    # If the mismatched samples make up at least 10% of all samples, stop execution and return an error.
    mismatched_percentage <- mismatched_length / total_length
    if (mismatched_percentage > 0.10) {
      stop("Percentage of sex discrepancies is too high (> 10%). Please check your input data!")
    }
  }
}

#
#Sex prediction done
#

# extract list of individuals with geno+methylation data
overlap <- intersect(colnames(norm.beta),fam[,2])
n.overlap <- length(overlap)
if(sorted_methylation == "no"){
if(n.overlap < 200)
{
	msg <- paste0("fewer than 200 subjects with methylation and genotype data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
}

w <- which(fam[,2] %in% overlap)
fam2 <- fam[w,1:2]
message(nrow(fam2), " individuals present in both genetic and methylation datasets")

message("Writing ID lists to ", ids, " and ", ids_plink)
if(file.exists(ids)) {
	message("Deleting old files")
	unlink(ids)
}
if(file.exists(ids_plink)) {
	message("Deleting old files")
	unlink(ids_plink)
}
write.table(fam2[,2],ids,sep="\t",quote=F,row.names=F,col.names=F)
write.table(fam2[,1:2],ids_plink,sep="\t",quote=F,row.names=F,col.names=F)


# Measured cell counts

message("Checking measured cell counts data: ", cellcounts_file_measured)

ccm.name<-"NULL"
if(cellcounts_file_measured != "NULL")
{
	ccm <- read.table(cellcounts_file_measured,header=T)
	ccm.name<-names(ccm)[-1]
	c1 <- dim(ccm)[1]
	c2 <- dim(ccm)[2]

	if(c1!=nid_meth)
	{
		msg <- paste0("number of samples in measured cell counts file is not the same as in beta matrix")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)   
	}

	w <- which(names(ccm)[1] %in% c("IID"))
	if(w!=1)
	{
		msg <- paste0("first column from measured cellcounts file should be the sample identifier with the name IID")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}

    
    w <- match(names(ccm),c("IID","Bcells","Tcells","Eos","Mono","Neu","Baso"))
	if(length(na.omit(w))!=length(names(ccm)))
	{
		msg <- paste0("the names in the measured cellcounts file do not match Bcells, Tcells, Eos, Mono, Neu, Baso")
		#errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}

	if(c2<3)
	{
		msg <- paste0("are there any columns with cell counts missing in the measured cell counts file?")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}

	a <- apply(ccm,2,function(x) y<-length(which(is.na(x))))
	if(length(which(a > 0.1*nid_meth)))
	{
		msg <- paste0("more than 10% of missingness in one of the measured cellcounts")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}
	message("Number of measured cell types: ", c2)
	message("Cell types:\n", paste(names(ccm)[-1], collapse="\n"))
}

if(cellcounts_file_measured == "NULL")
{	
	message("No cell counts are provided")
}



cohort_summary <- list()
cohort_summary$n_CpGs <- nrow(norm.beta)
cohort_summary$methylation_sample_size <- ncol(norm.beta)
cohort_summary$geno_meth_common_ids <- length(overlap)

#if(length(cc.name)>0){cohort_summary$predicted_cellcounts <- as.character(cc.name)}
#if(length(cc.name)==0){cohort_summary$measured_cellcounts <- "not_provided"}
#cohort_summary$predicted_cellcounts_type <- predicted_cellcounts_type
if(length(ccm.name)>0){cohort_summary$measured_cellcounts <- as.character(ccm.name)}
if(length(ccm.name)==0){cohort_summary$measured_cellcounts <- "not_provided"}

save(cohort_summary, file=cohort_descriptives_file)



summariseMeth <- function(X)
{
	message("Counting outliers in methylation matrix")
	
    norm.beta.copy <- X
	#for(i in 1:niter)
	#{
	#	sds <- rowSds(norm.beta.copy, na.rm=T)
	#	means <- rowMeans(norm.beta.copy, na.rm=T)
	#	norm.beta.copy[norm.beta.copy > means + sds*outlier_threshold | norm.beta.copy < means - sds*outlier_threshold] <- NA
	#}
	#outliers <- apply(norm.beta.copy, 1, function(x) sum(is.na(x)))
	
	norm.beta.copy<-ewaff.handle.outliers(norm.beta, method=c("iqr"), iqr.limit=3, winsorize.pct=0.05)
	outliers <-norm.beta.copy[[2]]

	message("Estimating means")
	means <- rowMeans(norm.beta.copy[[1]], na.rm=T)

	message("Estimating SDs")
	sds <- rowSds(norm.beta.copy[[1]], na.rm=T)

	message("Estimating medians")
	medians <- rowMedians(norm.beta.copy[[1]], na.rm=T)

	dat <- data.frame(cpg=rownames(norm.beta.copy[[1]]), mean=means, median=medians, sd=sds, outlier=outliers)
	return(dat)
}

#message("Generating summary stats of methylation")

#meth_summary <- summariseMeth(norm.beta)

#save(cohort_summary, file=cohort_descriptives_file)
#save(meth_summary, file=methylation_summary_file)
write.table(colnames(norm.beta), file=meth_ids_file, row=F, col=F, qu=F)


message("\n\nCompleted checks\n")

message("Summary of data:")
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
	q(status=1)
}
message("\n\n")
