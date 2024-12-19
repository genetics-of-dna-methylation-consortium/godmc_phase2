###############################discription###################################

# Thanks to Hannah Elliott for sending the original version of this script
# This script is dedicated to prediction on smoking score based on DNAm matrix

# If any errors occur while running this script, 
# please don't hesitate to reach out to the developer for assistance at s.w.wang@exeter.ac.uk.


###############################packages###################################

suppressPackageStartupMessages(library(impute))
suppressPackageStartupMessages(library(diptest))
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RPMM))

#################################Sub-functions################################

predict.smoking <- function(Illig_data, mbeta)
{
	# Remove CpGs missing from data
	Illig_data <- subset(Illig_data, cpgs %in% rownames(mbeta))

	# Split data into CPG sites that increase or decrease in smokers
	Illig_data_up <- subset(Illig_data, Illig_data$all_effect >=0)
	Illig_data_down <- subset(Illig_data, Illig_data$all_effect <0)

	# subset SABRE data and order
	mbeta_down <- mbeta[rownames(mbeta) %in% Illig_data_down$cpgs,]
	mbeta_up <- mbeta[rownames(mbeta) %in% Illig_data_up$cpgs, ]

  mbeta_down <- imputation(mbeta_down)
  mbeta_up <- imputation(mbeta_up)
	
	# sort Illig data by Cpg name
	Illig_data_up <- Illig_data_up[order(Illig_data_up$cpgs),]
	Illig_data_down <- Illig_data_down[order(Illig_data_down$cpgs),]

	# sort SABRE data by Cpg name
	mbeta_up <- mbeta_up[order(rownames(mbeta_up)),]
	mbeta_down <- mbeta_down[order(rownames(mbeta_down)),]

	# Calculate diffs between SABRE beta values and the reference for each CpG site
	matrix_up <- mbeta_up - Illig_data_up$reference_never_median_beta_all
	matrix_down <- Illig_data_down$reference_never_median_beta_all - mbeta_down

	# Calculate scores
	scores_up <- c(t(matrix_up) %*% Illig_data_up$weights)
	scores_down <- c(t(matrix_down) %*% Illig_data_down$weights)

	# Combine scores
	scores_combined <- scores_up + scores_down

	dat <- data.frame(IID = colnames(mbeta), Smoking = scores_combined)
	return(dat)
}


imputation <- function(x){
	# assumes you are only going to impute a few hundred sites. 
  nSamples = ncol(x)
  datMethUsed = t(x)
  dimnames1 = dimnames(datMethUsed)
  noMissingPerSample = rowSums(is.na(datMethUsed))
  
  if(nSamples > 1 & max(noMissingPerSample) > 0){
	  message("Imputing missing data")
	  datMethUsed = data.frame(t(impute.knn(datMethUsed)$data))
    colnames(datMethUsed) = dimnames1[[1]]
	  return(datMethUsed)
  } else {
	  message("No missing data to impute")
	  return(x)
  }

}

intersect <- function(m1, sd1, m2, sd2, p1, p2){
  
  B <- (m1/sd1^2 - m2/sd2^2)
  A <- 0.5*(1/sd2^2 - 1/sd1^2)
  C <- 0.5*(m2^2/sd2^2 - m1^2/sd1^2) - log((sd1/sd2)*(p2/p1))
  
  if (A!=0){
    (-B + c(1,-1)*sqrt(B^2 - 4*A*C))/(2*A)
  } else {-C/B}
} 

###############################main function###################################
main <- function()
{
	arguments <- commandArgs(T)
	methylationfile <- arguments[1]
	fam_file <- arguments[2]
	out_file <- arguments[3]
	smok_plot <- arguments[4]
	SD <- as.numeric(arguments[5])
	cov_file <- arguments[6]
	

	# Dataframe of Illig data (supplementary table 2) http://www.ncbi.nlm.nih.gov/pubmed/23691101
	load("resources/smoking/illig.RData")
	fam <- read.table(fam_file, stringsAsFactors=FALSE)[,1:2]
	colnames(fam) <- c("FID", "IID")
	covs <- read.table(cov_file, header=T, stringsAsFactors=FALSE)
	m <- match(fam[,2], covs[,"IID"])
	covs <- covs[m,]
	
	# Load methylation data - this contains the object "norm.beta"
	# This is an ncpg (rows) x nid (cols) matrix
	# Columns and rows labelled with ID and CPG respectively

	message("Loading methylation data")
	load(methylationfile)
	m <- match(fam[,2], colnames(norm.beta))
	norm.beta <- norm.beta[,m]
	
	message("Predicting smoking levels")
	smok <- predict.smoking(Illig_data, norm.beta)
	nom <- names(smok)
	smok <- merge(smok, fam, by.x="IID", by.y="IID")
	smok <- subset(smok, select=c("FID", nom))

	message("Fitting bimodality on distribution of smoking score.")
	biomod <- dip.test(smok$Smoking)
	if (biomod$p.value < 0.05){
	  message("The distribution of smoking score is bimodal.")
	}else{
	  message("The distribution of smoking score is unimodal.")
	}
	
	mixmdl <- normalmixEM(smok$Smoking)
	
	# find the intersection
	m1 <- mixmdl$mu[1]
	sd1 <- mixmdl$sigma[1]
	m2 <- mixmdl$mu[2]
	sd2 <- mixmdl$sigma[2]
	p1 <- mixmdl$lambda[1]
	p2 <- mixmdl$lambda[2]
	is <- intersect(m1,sd1,m2,sd2,p1,p2)
	d <- density(smok$Smoking)
	mini_d <- optimize(approxfun(d$x,d$y), interval=c(min(m1, m2),max(m1,m2)))$minimum
	

	if ( dim(smok)[1] > 5000 ){
		shapiroscore = signif(as.numeric(shapiro.test(sample(smok$Smoking, 3000))[2]),2)
	}else{
		shapiroscore = signif(as.numeric(shapiro.test(smok$Smoking)[2]),2)
	}
	
	pdf(smok_plot, width=12, height=10)
	par(mfrow=c(3,2))
	
	plot(smok$Smoking, xlab="", main=paste("Smoking prediction (N=", length(which(!is.na(smok$Smoking))),")",sep=""),cex.main=0.7)
	hist(smok$Smoking, xlab="", main=paste("Smoking prediction (N=", length(which(!is.na(smok$Smoking))),")",sep=""),cex.main=0.7)
	abline(v=mean(smok$Smoking,na.rm=T)-SD*sd(smok$Smoking,na.rm=T),lty=2)
	abline(v=mean(smok$Smoking,na.rm=T)+SD*sd(smok$Smoking,na.rm=T),lty=2)
	qqnorm(smok$Smoking, main=paste("Smoking prediction (N=", length(which(!is.na(smok$Smoking))),"; shapiroP=",shapiroscore,")",sep=""),cex.main=0.7)
	qqline(smok$Smoking)

	plot(covs$Age_numeric, smok$Smoking, xlab="Age", ylab="predicted smoking",main="Age vs Smoking prediction",cex.main=0.7)
	
	plot(density(smok$Smoking), xlab="", main=paste("Smoking prediction (N=", length(which(!is.na(smok$Smoking))),")",sep=""),cex.main=0.7)
	rug(smok$Smoking)
	abline(h=0, col = "grey")
	
	plot(mixmdl, which = 2)
	lines(density(smok$Smoking), lty=2, lwd=2)
	abline(v=is[2], col = "purple", lwd = 2, lty = 2)
	abline(v=mini_d, col = "orange", lwd = 2, lty = 2)
	
	write.table(smok, file=paste0(out_file, ".txt"), row=F, col=T, qu=F)
 
	dev.off()

}

main()

	

