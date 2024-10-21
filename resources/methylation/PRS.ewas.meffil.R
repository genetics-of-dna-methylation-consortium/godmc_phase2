suppressMessages(library(matrixStats))
suppressMessages(library(meffil))

# Perform EWAS

main <- function()
{

  arguments <- commandArgs(T)

  beta_file <- arguments[1]
  phen_name <- arguments[2]
  covs_file <- arguments[3]
  pc_file <- arguments[4]
  home_dir <- arguments[5]
  res_dir <- arguments[6] 
  study_name <- arguments[7]

  phen_file <- paste(home_dir,"/processed_data/genetic_data/PRS_",phen_name,".sscore",sep="")
  out_file <- paste(res_dir,"/",study_name,"_PRS_",phen_name,"_EWAS_results.RData",sep="")
  qqplot_file <- paste(res_dir,"/",study_name,"_PRS_",phen_name,"_EWAS_qqplot",sep="")
  report_file <- paste(res_dir,"/",study_name,"_PRS_",phen_name,"_EWAS",sep="")

	message("Loading non genetic methylation PCs")
  pcs<-read.table(paste(pc_file,".txt",sep=""),sep=" ",header=T)
  rownames(pcs) <- pcs$IID

  if(covs_file != "NULL") 
  {
    msg <- paste("Loading covariates for", phen_name)
 	  message(msg)

    covs <- read.table(covs_file, he=T, stringsAsFactors=FALSE)
    g<-grep("factor",names(covs))
    if(length(g)>1){ 
      for (i in 1:length(g)){
        covs[,g[i]]<-as.factor(covs[,g[i]])
        if(length(levels(covs[,g[i]]))==1)
        covs<-covs[,-g[i]]
      }
    }

    g<-grep("numeric",names(covs))
    if(length(g)>1){ 
      for (i in 1:length(g)){
        covs[,g[i]]<-as.numeric(as.character(covs[,g[i]]))
      }
    }

    rownames(covs) <- covs$IID
	
    m<-match(covs$IID,pcs$IID)
    covs<-data.frame(covs,pcs[m,-1])
  }
    
  if(covs_file == "NULL") 
  {
    msg <- paste("No covariates file present for", phen_name)
    message(msg)

    covs <- pcs
  }

	if(nrow(covs) < 10)
	{
		message("There are fewer than 10 individuals remaining. Stopping the analysis.")
		q()
	}

	message("Loading methylation data")
 	load(beta_file)

	message(paste("Loading PRS data for",phen_name))
	phen <- read.table(phen_file, header=T, stringsAsFactors=FALSE)
  rownames(phen) <- phen$IID

  phen<-phen[,which(names(phen)%in%c("IID","SCORE"))]

	phen <- subset(phen, IID %in% colnames(norm.beta), select="SCORE")
	if(sum(!is.na(phen[["SCORE"]])) < 10)
	{
		message("There are fewer than 10 individuals remaining. Stopping the analysis.")
		q()
	}
    
  phen<-na.omit(phen)

	norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(phen)]
	phen <- phen[match(colnames(norm.beta), rownames(phen)), , drop=FALSE]
	stopifnot(all(rownames(phen) == colnames(norm.beta)))
    
  covs <- covs[match(colnames(norm.beta), rownames(covs)), , drop=FALSE]
  covs <- subset(covs, , select=colnames(covs)[!colnames(covs)%in%c("IID")])
  nr<-nrow(norm.beta)
  nc<-ncol(norm.beta)

  message("\nPerforming EWAS for ", phen_name, " on ",nc," individuals and ",nr," CpGs")    

  if(ncol(covs) > 0){ 
    ewas.ret <- meffil.ewas(norm.beta, variable=phen[,1], covariates=covs,winsorize.pct = NA, most.variable = min(nrow(norm.beta), 20000), sva=F, isva=T)       
  }
  
  if(ncol(covs) == 0){ 
    ewas.ret <- meffil.ewas(norm.beta, variable=phen[,1], winsorize.pct = NA, most.variable = min(nrow(norm.beta), 20000), sva=F, isva=T)       
  }
 
  res<-as.data.frame(ewas.ret$p.value)
  
	message("Generating Q-Q plot")
  qqplot_pval(res$none, file=paste(qqplot_file,"nocovs.pdf",sep="."))
  if(ncol(covs) > 0){ qqplot_pval(res$all, file=paste(qqplot_file,"allcovs.pdf",sep=".")) }
  if(ncol(covs) == 0){ message("no QQ plot provided with covariate adjustment: no covariates provided") }  
  qqplot_pval(res$isva, file=paste(qqplot_file,"isvacovs.pdf",sep="."))

  main_model <- "none" #model with no additional covariates for the EWAS report
    
  ewas.parameters <- meffil.ewas.parameters(sig.threshold=1e-7, max.plots=100,qq.inflation.method="regression", model=main_model)
  ewas.summary<-meffil.ewas.summary(ewas.ret,norm.beta,parameters=ewas.parameters)                              
  meffil.ewas.report(ewas.summary, output.file=paste(report_file,".ewas.report.html",sep=""))

  #remove items not needed and save
  ewas.ret$samples <- NULL
  ewas.ret$variable <- NULL
  if(ncol(covs) > 0){
    ewas.ret$covariates <- NULL
    ewas.ret$analyses$all$design <- NULL
  }  
  ewas.ret$analyses$none$design <- NULL  
  ewas.ret$analyses$isva$design <- NULL
  ewas.ret$sva.ret$isv <- NULL
    
	message("Saving results")
	save(ewas.ret, file=out_file)



}


est_lambda <- function (data, plot = FALSE, proportion = 1, method = "regression", filter = TRUE, df = 1, ...)
{
	data <- data[which(!is.na(data))]
	if (proportion > 1 || proportion <= 0) 
		stop("proportion argument should be greater then zero and less than or equal to one")
	ntp <- round(proportion * length(data))
	if (ntp < 1) 
		stop("no valid measurements")
	if (ntp == 1) {
		warning(paste("One measurement, lambda = 1 returned"))
		return(list(estimate = 1, se = 999.99))
	}
	if (ntp < 10) 
		warning(paste("number of points is too small:", ntp))
	if (min(data) < 0) 
		stop("data argument has values <0")
	if (max(data) <= 1) {
		data <- qchisq(data, 1, lower.tail = FALSE)
	}
	if (filter) {
		data[which(abs(data) < 1e-08)] <- NA
	}
	data <- sort(data)
	ppoi <- ppoints(data)
	ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
	data <- data[1:ntp]
	ppoi <- ppoi[1:ntp]
	out <- list()
	if (method == "regression") {
		s <- summary(lm(data ~ 0 + ppoi))$coeff
		out$estimate <- s[1, 1]
		out$se <- s[1, 2]
	}
	else if (method == "median") {
		out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 
			df)
		out$se <- NA
	}
	else {
		stop("'method' should be either 'regression' or 'median'!")
	}
	if (plot) {
		lim <- c(0, max(data, ppoi, na.rm = TRUE))
		oldmargins <- par()$mar
		par(mar = oldmargins + 0.2)
		plot(ppoi, data, xlab = expression("Expected " ~ chi^2), 
			ylab = expression("Observed " ~ chi^2), pch=".",...)
		abline(a = 0, b = 1)
		abline(a = 0, b = out$estimate, col = "red")
		par(mar = oldmargins)
	}
	out
}

qqplot_pval <- function(P, filename=NULL)
{
	l <- est_lambda(P, method="median")
	nom <- paste("lambda = ", round(l$estimate, 3), sep="")
	if(!is.null(filename))
	{
		pdf(filename)
	}
	est_lambda(P, method="median", plot=TRUE, main=nom)
	if(!is.null(filename))
	{
		dev.off()
	}
}


main()
