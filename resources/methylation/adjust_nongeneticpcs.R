library(parallel)
library(tidyr)
suppressMessages(library(meffil))
suppressMessages(library(matrixStats))

main <- function()
{
	arguments <- commandArgs(T)

	methylationfile <- arguments[1]
	cov_file <- arguments[2]
	out_file <- arguments[3]
	transform <- arguments[4]
        nthreads <- as.numeric(arguments[5])
	chunks <- as.numeric(arguments[6])
	jid <- as.numeric(arguments[7])
       
	message("Reading methylation data...")
	load(methylationfile)

	if(!is.na(jid))
	{
		chunksize <- ceiling(nrow(norm.beta) / chunks)
		i1 <- chunksize * (jid-1) + 1
		i2 <- min(nrow(norm.beta), chunksize * jid)
		norm.beta <- norm.beta[i1:i2,]
		out_file <- paste0(out_file, ".", jid, ".RData")
	} else {
		out_file <- paste0(out_file, ".RData")
	}

	message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")

	covs <- read.table(cov_file, he=T)
	index <- apply(covs, 1, function(x) any(is.na(x) | is.nan(x) | is.infinite(x)))
	covs <- covs[!index, ]
	rownames(covs) <- covs$IID
	covs <- subset(covs, IID %in% colnames(norm.beta), select=-c(IID))

	norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(covs)]
	index_beta <- rowMeans(is.nan(norm.beta)) == 1
	norm.beta <- norm.beta[!index_beta, ]
	covs <- covs[match(colnames(norm.beta), rownames(covs)), ]
	stopifnot(all(rownames(covs) == colnames(norm.beta)))

	g <- grep("_factor",names(covs))

        if (length(g)>0){
		for (i in 1:length(g))
		{
			covs[,g[i]]<-as.factor(covs[,g[i]])
        	}
        }

	norm.beta.copy <- is.na(norm.beta)

	if(is.na(nthreads) | nthreads == 1)
	{
		norm.beta <- adjust.covs.serial(norm.beta, covs, transform)
	} else {
		message("Running with ", nthreads, " threads")
		norm.beta <- adjust.covs(norm.beta, covs, nthreads, transform)	
	}
	
	norm.beta[norm.beta.copy] <- NA

	index <- which(is.na(norm.beta), arr.ind = TRUE)
	if (length(index)>0){
    	message("Replace ",length(index)," missing values with rowmeans")
    	norm.beta[index] <- rowMeans(norm.beta, na.rm = TRUE)[index[, "row"]] }

	save(norm.beta, file=out_file)
	message("Successfully completed adjustments")
}


adjust.covs.1 <- function(x, covs, transform)
{	
    if(transform == "transformed"){
	    d <- data.frame(X=x, covs)
	    form <- as.formula(paste0("X ~ ", paste(names(d)[-1], collapse=" + ")))
	    as.numeric(rntransform(residuals(lm(form, data=d, na.action=na.exclude))))
    } else{
		d <- data.frame(X=x, covs)
		form <- as.formula(paste0("X ~ ", paste(names(d)[-1], collapse=" + ")))
		as.numeric(residuals(lm(form, data=d, na.action=na.exclude)))
	}
}

adjust.covs <- function(B, covs, mc.cores=mc.cores, transform)
{
	l1 <- get.index.list(nrow(B), mc.cores)
	l <- lapply(l1, function(ii)
	{
		res <- mclapply(ii, function(i)
		{
			message("Probe ", i, " of ", nrow(B))
			adjust.covs.1(B[i,], covs, transform)
		}, mc.cores=mc.cores)
		return(do.call(rbind, res))
	})
	l <- do.call(rbind, l)
	rownames(l) <- rownames(B)
	colnames(l) <- colnames(B)
	return(l)
}

get.index.list <- function(n, mc.cores)
{
	mc.cores <- ifelse(mc.cores < 1, 1, mc.cores)
	div <- floor(n / mc.cores)
	rem <- n %% mc.cores
	l1 <- lapply(1:div, function(x) (x-1) * mc.cores + 1:mc.cores)
	if(rem != 0) l1[[div+1]] <- l1[[div]][mc.cores] + 1:rem
	return(l1)
}

adjust.covs.serial <- function(B, covs, transform)
{
	for(i in 1:nrow(B))
	{
		cat(i, "\n")
		B[i, ] <- adjust.covs.1(B[i,], covs, transform)
	}
	return(B)
}

rntransform <- function(x)
{
	out <- rank(x) - 0.5
	out[is.na(x)] <- NA
	mP <- 0.5/max(out, na.rm = T)
	out <- out/(max(out, na.rm = T) + 0.5)
	out <- scale(qnorm(out))
	out
}

main()
