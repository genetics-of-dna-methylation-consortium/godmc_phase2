library(MatrixEQTL)
library(parallel)

main <- function()
{
	arguments <- commandArgs(T)
	geno_root <- arguments[1]
	phen_file <- arguments[2]
	out_file <- arguments[3]
	threshold <- as.numeric(arguments[4])
	mc.cores <- as.numeric(arguments[5])

	message("Loading genetic data")
	slicesize <- 1000

	gene <- SlicedData$new()
	gene$fileDelimiter = "\t"
	gene$fileOmitCharacters = "NA"
	gene$fileSkipRows = 1
	gene$fileSkipColumns = 1
	gene$fileSliceSize = slicesize
	gene$LoadFile( paste0(phen_file, ".txt"))

	bn <- basename(geno_root)
	dn <- dirname(geno_root)
	geno_file <- list.files(dn, pattern=paste0("^", bn))
	
	l <- run_all_chunks(dn, geno_file, gene, threshold, slicesize, mc.cores)

	load(paste0(phen_file, ".RData")) #loads pc
	message("\n\nIdentified ", length(l), " out of ", ncol(pc), " PCs with significant genetic component")

	if(length(l) == ncol(pc))
	{
		message("It appears that all the PCs have a genetic component\n",
			"This is a little worrying because it suggests family effects or stratification\n",
			"Please check that 03b is adjusting for these factors\n",
			"You could also try increasing the 'n_meth_pcs and/or meth_pc_cutoff' value in the config file.\n"
		)
		q()
	}

	pc <- pc[,!colnames(pc) %in% l, drop=FALSE]
	
	pc <- data.frame(IID=rownames(pc), pc)
	write.table(pc, file=paste0(out_file, ".txt"), row=F, col=T, qu=F)

}

run_all_chunks_serial <- function(dn, geno_file, gene, threshold, slicesize)
{
	l <- list()
	for(i in 1:length(geno_file))
	{
		message(i, " of ", length(geno_file))
			snps = SlicedData$new()
			snps$fileDelimiter = "\t"
			snps$fileOmitCharacters = "NA"
			snps$fileSkipRows = 1
			snps$fileSkipColumns = 1
			snps$fileSliceSize = slicesize
			snps$LoadFile( file.path(dn, geno_file[i]) )
			ids <- Reduce(intersect, list(snps$columnNames, gene$columnNames))
		l[[i]] <- run_chunk(snps, gene, threshold, slicesize)
	}
	return(unique(unlist(l)))
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

run_all_chunks <- function(dn, geno_file, gene, threshold, slicesize, mc.cores)
{
	message("Performing analysis using ", mc.cores, " threads")
	l1 <- get.index.list(length(geno_file), mc.cores)
	l <- lapply(l1, function(ii)
	{
		res <- mclapply(ii, function(i)
		{
			message("Chunk ", i, " of ", length(geno_file))
			snps = SlicedData$new()
			snps$fileDelimiter = "\t"
			snps$fileOmitCharacters = "NA"
			snps$fileSkipRows = 1
			snps$fileSkipColumns = 1
			snps$fileSliceSize = slicesize
			snps$LoadFile( file.path(dn, geno_file[i]) )
			run_chunk(snps, gene, threshold, slicesize)
		}, mc.cores=mc.cores, mc.preschedule=FALSE)
		return(res)
	})
	l <- unique(unlist(l))
	return(l)
}

run_chunk <- function(snps, gene, threshold, slicesize)
{

	useModel = modelLINEAR
	errorCovariance = numeric()

	ids <- Reduce(intersect, list(snps$columnNames, gene$columnNames))

	snps$ColumnSubsample(match(ids, snps$columnNames))
	gene$ColumnSubsample(match(ids, gene$columnNames))

	stopifnot(all(snps$columnNames==gene$columnNames))

	me <- Matrix_eQTL_engine(
		snps = snps,
		gene = gene,
		cvrt = SlicedData$new(),
		output_file_name = NULL,
		pvOutputThreshold = threshold,
		useModel = useModel, 
		errorCovariance = errorCovariance, 
		verbose = FALSE,
		pvalue.hist = FALSE,
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE
	)
	snp<-as.character(me$all$eqtls$snps)
    	spl<-do.call("rbind",strsplit(snp,":"))
    
    	if(!is.null(spl)){
        	spl2<-do.call("rbind",strsplit(spl[,2],"_"))
        	chrpos<-data.frame(chr=spl[,1],pos=spl2[,1])
        	w<-which(chrpos$chr=="6"&chrpos$pos>25000000&chrpos$pos<35000000)
	    	if(length(w)>0){me$all$eqtls<-me$all$eqtls[-w,]}
	}
    
    return(as.character(unique(me$all$eqtls$gene)))
}


main()


