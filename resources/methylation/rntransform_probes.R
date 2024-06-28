library(parallel)
suppressMessages(library(matrixStats))


main <- function()
{
	arguments <- commandArgs(T)

	methylationfile <- arguments[1]
	out_file <- arguments[2]
	nthreads <- as.numeric(arguments[3])
	chunks <- as.numeric(arguments[4])
	jid <- as.numeric(arguments[5])


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

	message("Identifying methylation outliers")

	niter <- 3
	outlier_threshold <- 10
	norm.beta.copy <- norm.beta
	for(i in 1:niter)
	{
		sds <- rowSds(norm.beta.copy, na.rm=T)
		means <- rowMeans(norm.beta.copy, na.rm=T)
		norm.beta.copy[norm.beta.copy > means + sds*outlier_threshold | norm.beta.copy < means - sds*outlier_threshold] <- NA
	}
	outlier_count <- apply(norm.beta.copy, 1, function(x) sum(is.na(x)))
	norm.beta.copy <- is.na(norm.beta.copy)

	if(is.na(nthreads) | nthreads == 1)
	{
		norm.beta <- rntransform_betas.serial(norm.beta)
	} else {
		message("Running with ", nthreads, " threads")
		norm.beta <- rntransform_betas(norm.beta, nthreads)	
	}
	
	norm.beta[norm.beta.copy] <- NA
	save(norm.beta, file=out_file)
	message("Successfully completed adjustments")
}


rntransform_betas <- function(B, mc.cores=mc.cores)
{
	l1 <- get.index.list(nrow(B), mc.cores)
	l <- lapply(l1, function(ii)
	{
		res <- mclapply(ii, function(i)
		{
			message("Probe ", i, " of ", nrow(B))
			rntransform(B[i,])
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

adjust.covs.serial <- function(B, covs)
{
	for(i in 1:nrow(B))
	{
		cat(i, "\n")
		B[i, ] <- adjust.covs.1(B[i,], covs)
	}
	return(B)

	# apply(B, 1, function(x) adjust.relatedness.1(x, kin))
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

removeOutliers <- function(x, thresh, remove=FALSE)
{
	m <- mean(x, na.rm=T)
	s <- sd(x, na.rm=T)
	index <- x > m + thresh*s | x < m - thresh*s

	if(remove)
	{
		x <- x[!index]
	} else {
		x[index] <- NA
	}
	return(x)
}

remRec <- function(x, thresh, iterations)
{
	d <- array(0, iterations+1)
	d[1] <- sum(!is.na(x))
	for(i in 1:iterations)
	{
		x <- removeOutliers(x, thresh)
		d[i+1] <- sum(!is.na(x))
	}
	return(list(x=x, its=d))
}



main()
