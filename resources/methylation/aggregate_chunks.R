suppressMessages(library(dplyr))
suppressMessages(library(meffil))

main <- function(){
	arguments <- commandArgs(T)
	rootname <- arguments[1]
	n_chunks <- as.numeric(arguments[2])
	meth_array <- arguments[3]
	methylation_no_outlier <- arguments[4]
        rootname1 <- arguments[5]
	
        message("Checking all files are present")

	if(is.na(rootname1)){
    		fs <- paste0(rootname, ".", 1:n_chunks, ".RData")
    	}else{
		fs <- paste0(rootname, ".", rootname1, ".", 1:n_chunks, ".RData")
	}

	index <- file.exists(fs)
		if(any(!index) & is.na(rootname1))
	{
		message("Warning: The following files containing autosome CpGs don't exist, please check them before continuing.\n", paste(fs[!index], collapse="\n"))
		message(paste0("If all the CpGs in that methylation chunk ", fs[!index], " are on sex chromosomes, please ignore this warning .\n", collapse="\n"))
	}

	fs <- fs[index]
	load(fs[1])

    	if(grepl("pc", rootname)){
        	betas <- norm.beta
    	}

    	if(is.na(rootname1)){
	    message("Reading in data")
	    
	    cpg <- list()
	    for(i in 1:length(fs)){
		    load(fs[i])
		    if(grepl("pc", rootname)){
        		betas <- norm.beta
		    }
            cpg[[i]] <- data.frame(betas)
            message(i, " of ", length(fs))
	    }
    
            beta.temp <- as.matrix(bind_rows(cpg))
            colnames(beta.temp) <- colnames(betas)
            beta.sub <- check_cpg_number(rootname1, methylation_no_outlier, meth_array)

            if(all(rownames(beta.sub) == rownames(beta.temp))){
		norm.beta <- beta.temp
	        save(norm.beta, file=paste0(rootname, ".RData"))
	        message("The number of CpGs before and after adjustment are consistent")
                if(is.na(rootname1)){
                message(paste0("Successfully aggregated all RData chunks of ",rootname, " autosomal CpGs"))
                }else{message(paste0("Successfully aggregated all RData chunks of ",rootname, ".", rootname1))}
            }else{
                message("The number of CpGs before and after adjustment are inconsistent")
                q()
                }

        }else{
            f <- \(x) {
            fenv <- new.env()
            load(x, envir=fenv)
            get(ls(fenv), envir=fenv)
            }
        alldat<-do.call(rbind, lapply(fs, f))
        beta.sub <- check_cpg_number(rootname1, methylation_no_outlier, meth_array)

        if(all(rownames(beta.sub) == rownames(alldat))){
            norm.beta <- alldat
            save(norm.beta, file=paste0(rootname,".", rootname1, ".RData"))
            message("The number of CpGs before and after adjustment are consistent")
            message(paste0("Successfully aggregated all RData chunks of ",rootname, ".", rootname1))
        }else{
            message("The number of CpGs before and after adjustment are inconsistent")
            q()
        }
    }
}

check_cpg_number <- function(rootname1, methylation_no_outlier, meth_array){
    # check if the number of CpGs in aggregated data is consistent with the number of CpGs before adjustment
    load(methylation_no_outlier)
    annots <- meffil.get.features(meth_array)
    annots <- annots[!is.na(annots$chromosome),]
    x_probes <- annots$name[annots$chromosome == "chrX"]
    y_probes <- annots$name[annots$chromosome == "chrY"]
    autosomal_probes <- annots$name[!annots$chromosome %in% c("chrX","chrY")]
    
    if(is.na(rootname1)){
        beta.sub <- norm.beta[rownames(norm.beta) %in% autosomal_probes,]
    }else{
        if(rootname1 %in% c("Female.chrX","Male.chrX")){
            beta.sub <- norm.beta[rownames(norm.beta) %in% x_probes,] 
        }
        if(rootname1 %in% c("Male.chrY")){
            beta.sub <- norm.beta[rownames(norm.beta) %in% y_probes,]
        }
    }
    return(beta.sub)
}

main()
