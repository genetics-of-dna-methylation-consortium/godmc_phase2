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
	meth_array <- arguments[8]
    
	message("Reading methylation data...")
	load(methylationfile)

	if(!is.na(jid))
	{
		chunksize <- ceiling(nrow(norm.beta) / chunks)
		i1 <- chunksize * (jid-1) + 1
		i2 <- min(nrow(norm.beta), chunksize * jid)
		norm.beta <- norm.beta[i1:i2,]
	}

	message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")

	covs <- read.table(cov_file, he=T, stringsAsFactors=F,colClasse=c("Sex_factor"="character"))
	covs <- covs[,!colnames(covs)=="Treg"]

	rownames(covs) <- covs$IID
	covs <- subset(covs, IID %in% colnames(norm.beta), select=-c(IID))

	norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(covs)]
	covs <- covs[match(colnames(norm.beta), rownames(covs)), ]
	stopifnot(all(rownames(covs) == colnames(norm.beta)))

	## process CpGs on sex chromosomes
	annots <- meffil.get.features(meth_array)
   	annots <- annots[!is.na(annots$chromosome),]
   	x_probes <- annots$name[annots$chromosome == "chrX"]
   	y_probes <- annots$name[annots$chromosome == "chrY"]
  
    	beta.x.female <- norm.beta[rownames(norm.beta) %in% x_probes, colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="F"]]
   	beta.x.male <- norm.beta[rownames(norm.beta) %in% x_probes, colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="M"]]
   	beta.y.male <- norm.beta[rownames(norm.beta) %in% y_probes, colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="M"]]
    	
    	if(is.matrix(beta.x.female) == F){
   		beta.x.female <- matrix(beta.x.female, nrow=1)
   	}
	rownames(beta.x.female) <- rownames(norm.beta)[rownames(norm.beta) %in% x_probes]
   	colnames(beta.x.female) <- colnames(norm.beta)[colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="F"]]
   	

    	if(is.matrix(beta.x.male) == F){
   		beta.x.male <- matrix(beta.x.male, nrow=1)
	}
	rownames(beta.x.male) <- rownames(norm.beta)[rownames(norm.beta) %in% x_probes]
   	colnames(beta.x.male) <- colnames(norm.beta)[colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="M"]]
   	

    	if(is.matrix(beta.y.male) == F){
   		beta.y.male <- matrix(beta.y.male, nrow=1)
   	}
	rownames(beta.y.male) <- rownames(norm.beta)[rownames(norm.beta) %in% y_probes]
  	colnames(beta.y.male) <- colnames(norm.beta)[colnames(norm.beta) %in% rownames(covs)[covs$Sex_factor=="M"]]
    ## process CpGs on autosomes
        autosomal_probes <- annots$name[!annots$chromosome %in% c("chrX","chrY")]
        norm.beta <- norm.beta[rownames(norm.beta) %in% autosomal_probes,]
  	
    	g <- grep("_factor",names(covs))
    	if (length(g)>0){
        	for (i in 1:length(g))
        	{
            	covs[,g[i]]<-as.factor(covs[,g[i]])
        	}
    	}

	covs.female <- covs[covs$Sex_factor=="F",]
	covs.male <- covs[covs$Sex_factor=="M",]

    	if(nrow(norm.beta)>0 & ncol(norm.beta)>0){
		index<-sapply(covs,function(.col){all(is.na(.col) | .col[1L] == .col)})
    		index[is.na(index)] <- FALSE
		covs <- covs[match(colnames(norm.beta), rownames(covs)),!index]
        	run.adjust.cov(norm.beta, covs, nthreads, transform, paste0(out_file,".",jid,".RData"))
    	}

    	if(nrow(beta.x.female)>0 & ncol(beta.x.female)>0){
		index_female <-sapply(covs.female,function(.col){all(is.na(.col) | .col[1L] == .col)})
    		index_female[is.na(index_female)] <- FALSE
		covs.female <- covs.female[match(colnames(beta.x.female), rownames(covs.female)),!index_female]
        	run.adjust.cov(beta.x.female, covs.female, nthreads=1, transform, paste0(out_file, ".Female.chrX.", jid, ".RData"))
    	}

    	if(nrow(beta.x.male)>0 & ncol(beta.x.male)>0){
		index_male <-sapply(covs.male,function(.col){all(is.na(.col) | .col[1L] == .col)})
    		index_male[is.na(index_male)] <- FALSE
		covs.male <- covs.male[match(colnames(beta.x.male), rownames(covs.male)),!index_male]
        	run.adjust.cov(beta.x.male, covs.male, nthreads=1, transform, paste0(out_file,".Male.chrX.", jid, ".RData"))
    	}

    	if(nrow(beta.y.male)>0 & ncol(beta.y.male)>0){
        	index_male <-sapply(covs.male,function(.col){all(is.na(.col) | .col[1L] == .col)})
    		index_male[is.na(index_male)] <- FALSE
		covs.male <- covs.male[match(colnames(beta.y.male), rownames(covs.male)),!index_male]
		run.adjust.cov(beta.y.male, covs.male, nthreads=1, transform, paste0(out_file,".Male.chrY.", jid, ".RData"))
    	}
}

run.adjust.cov <- function(betas, covs, nthreads, transform, out_file)
{
    if(nrow(betas)>0){
        betas.copy <- is.na(betas)

        if(is.na(nthreads) | nthreads == 1){
            betas <- adjust.covs.serial(betas, covs, transform)
        } else {
            message("Running with ", nthreads, " threads")
            betas <- adjust.covs(betas, covs, nthreads, transform)
        }

        betas[betas.copy] <- NA

        index <- which(is.na(betas), arr.ind = TRUE)
        if (length(index)>0){
            message("Replace ",length(index)," missing values with rowmeans")
            betas[index] <- rowMeans(betas, na.rm = TRUE)[index[, "row"]] 
        }
        save(betas, file=out_file)
    }
}

adjust.covs.1 <- function(x, covs, transform)
{	
	if(transform == "transformed"){
		d <- data.frame(X=rntransform(x), covs)
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
		}, mc.cores=mc.cores, mc.preschedule=FALSE)
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
		message("Probe ",i, " of ", nrow(B))
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
