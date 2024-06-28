library(plyr)

arguments <- commandArgs(T)

rootname <- arguments[1]
n_chunks <- as.numeric(arguments[2])
out <- arguments[3]
rootname1 <- arguments[4]

main <- function(){

    if(is.na(rootname1)){
   	fs <- paste0(rootname, ".", 1:n_chunks, ".RData_classes")
    }else{
        fs <- paste0(rootname, ".", rootname1, ".", 1:n_chunks, ".RData_classes")
    }

    index <- file.exists(fs)
    fs <- fs[index]

    message("Reading in data")
    l <- list()
    for(i in c(1:length(fs))){
		load(fs[i])
		message(i, " of ", length(fs))
		l[[i]] <- classes
    }

    classes <- rbind.fill(l)
    save(classes, file=out)
    message(sum(classes$cl == "try-error"), " out of ", nrow(classes), " probes couldn't be adjusted for polygenic effects")
    message(paste0("Successfully aggregated all ", rootname, " chunks class RData"))
}

main()

