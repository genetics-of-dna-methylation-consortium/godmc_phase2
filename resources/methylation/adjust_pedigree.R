library(parallel)
library(MASS)
library(tidyr)
library(dplyr)
suppressMessages(library(meffil))
suppressMessages(library(matrixStats))

main <- function()
{
  arguments <- commandArgs(T)
  
  filepath <- arguments[1]
  methylationfile <- arguments[2]
  grmfile <- arguments[3]
  cov_file <- arguments[4]
  out_file <- arguments[5]
  transform <- arguments[6]
  nthreads <- as.numeric(arguments[7])
  chunks <- as.numeric(arguments[8])
  jid <- as.numeric(arguments[9])
  meth_array <- arguments[10]
  
  source(paste0(filepath, "/resources/methylation/polygenic_genable.R"))
  source(paste0(filepath, "/resources/methylation/polylik_genable.R"))
  
  message("Reading methylation data...")
  load(methylationfile)
  
  if(!is.na(jid))
  {
    chunksize <- ceiling(nrow(norm.beta) / chunks)
    i1 <- chunksize * (jid-1) + 1
    i2 <- min(nrow(norm.beta), chunksize * jid)
    norm.beta <- norm.beta[i1:i2,]
  } 
  
  # Remove all IDs that have any NAs in the covariate file
  covs <- read.table(cov_file, he=T, stringsAsFactors=F, colClasses=c("Sex_factor"="character"))
  covs <- covs[,!colnames(covs)=="Treg"]
  covs <- subset(covs, IID %in% colnames(norm.beta))
  g <- grep("_factor",names(covs))
  if(length(g) > 0)
  {
    for (i in 1:length(g))
    {
      covs[,g[i]]<-as.factor(covs[,g[i]])
    }
  }
  
  index<-sapply(covs,function(.col){all(is.na(.col) | .col[1L] == .col)})
  index[is.na(index)] <- FALSE
  covs <- covs[,!index]
  
  grm <- readGRM(grmfile)
  kin <- makeGRMmatrix(grm)
  
  ## observe the intersect IDs across grm, covs, betas to ensure the same order of the three datasets
  common_ids <- Reduce(intersect, list(colnames(norm.beta), covs$IID, colnames(kin)))
  female_ids <- covs$IID[covs$Sex_factor=="F"]
  male_ids <- covs$IID[covs$Sex_factor=="M"]
  
  covs.all <- covs[match(common_ids, covs$IID),]
  kin.all <- kin[match(common_ids, colnames(kin)), match(common_ids, rownames(kin))]

  covs.female <- covs[match(female_ids, covs$IID),]
  covs.male <- covs[match(male_ids, covs$IID),]
  rownames(covs.female) <- covs.female$IID
  rownames(covs.male) <- covs.male$IID
  
  kin.female <- kin[match(female_ids, colnames(kin)), match(female_ids, rownames(kin))]
  kin.male <- kin[match(male_ids, colnames(kin)), match(male_ids, rownames(kin))]
  
  ## process CpGs on sex chromosomes
  annots <- meffil.get.features(meth_array)
  iannots <- annots[!is.na(annots$chromosome),]
  x_probes <- annots$name[annots$chromosome == "chrX"]
  y_probes <- annots$name[annots$chromosome == "chrY"]
  
  beta.x.female <- norm.beta[rownames(norm.beta) %in% x_probes, match(female_ids, colnames(norm.beta))] %>% as.matrix()
  beta.x.male <- norm.beta[rownames(norm.beta) %in% x_probes, match(male_ids, colnames(norm.beta))] %>% as.matrix()
  beta.y.male <- norm.beta[rownames(norm.beta) %in% y_probes, match(male_ids, colnames(norm.beta))] %>% as.matrix()
  beta.autosal.all <- norm.beta[(rownames(norm.beta) %in% c(x_probes, y_probes) == F),
                                match(common_ids, colnames(norm.beta))] %>% as.matrix()
  
  message("Kinship matrix of all samples ", nrow(kin.all), " by ", nrow(kin.all))
  message("Kinship matrix of female samples ", nrow(kin.female), " by ", nrow(kin.female))
  message("Kinship matrix of all samples ", nrow(kin.male), " by ", nrow(kin.male))
  
  if(nrow(norm.beta)>0 & ncol(norm.beta)>0 & !file.exists(paste0(out_file,".",jid,".RData"))){
  eig.all <- cal.eig(kin.all)
  beta.autosal.all <- transpose_check(beta.autosal.all)  
  run.adjust.cov(beta.autosal.all, covs.all %>% select(-IID), nthreads, kin.all, eig.all, transform, paste0(out_file,".",jid,".RData"))
  }
  
  if(nrow(beta.x.female)>0 & ncol(beta.x.female)>0 & !file.exists(paste0(out_file,".Female.chrX.",jid,".RData"))){
  eig.female <- cal.eig(kin.female)
  probename <- rownames(norm.beta)[rownames(norm.beta) %in% x_probes]
  beta.x.female <- transpose_check(beta.x.female, probename)
  run.adjust.cov(beta.x.female, covs.female %>% select(-IID, -Sex_factor), nthreads=1, kin.female, eig.female, transform, paste0(out_file, ".Female.chrX.", jid, ".RData"))
  }
  
  if(nrow(beta.x.male)>0 & ncol(beta.x.male)>0 & !file.exists(paste0(out_file,".Male.chrX.",jid,".RData"))){
  eig.male <- cal.eig(kin.male)
  probename <- rownames(norm.beta)[rownames(norm.beta) %in% x_probes]
  beta.x.male <- transpose_check(beta.x.male, probename)
  run.adjust.cov(beta.x.male, covs.male %>% select(-IID, -Sex_factor), nthreads=1, kin.male, eig.male, transform, paste0(out_file,".Male.chrX.", jid, ".RData"))
  }
  
  if(nrow(beta.y.male)>0 & ncol(beta.y.male)>0 & !file.exists(paste0(out_file,".Male.chrY.",jid,".RData"))){
  eig.male <- cal.eig(kin.male)
  probename <- rownames(norm.beta)[rownames(norm.beta) %in% y_probes]
  beta.y.male <- transpose_check(beta.y.male, probename)
  run.adjust.cov(beta.y.male, covs.male %>% select(-IID, -Sex_factor), nthreads=1, kin.male, eig.male, transform, paste0(out_file,".Male.chrY.", jid, ".RData"))
  }
}

cal.eig <- function(kin){
  relmat <- kin * 2
  tmp <- t(relmat)
  relmat[upper.tri(relmat)] <- tmp[upper.tri(tmp)]
  eig <- eigen(relmat, symmetric=TRUE)
  class(eig) <- "list"
  message(class(eig))
  print(str(eig))
  return(eig)
}

transpose_check <- function(mat, probename){
  if(ncol(mat) == 1){
    mat_tmp <- t(mat)
    rownames(mat_tmp) <- probename
} else{
    mat_tmp <- mat}
  return(mat_tmp)
}

run.adjust.cov <- function(betas, covs, nthreads, kin, eig, transform,  out_file)
{
  if(nrow(betas) > 0){
    betas.copy <- is.na(betas)
    
    if(is.na(nthreads) | nthreads == 1){
      out <- adjust.relatedness.serial(betas, covs, kin, eig, transform)
    } else {
      message("Running with ", nthreads, " threads")
      out <- adjust.relatedness(betas, covs, kin, eig, nthreads, transform)
  c  }
    
    betas <- out$x
    classes <- data.frame(cpg=rownames(betas), cl=out$cl)
    betas[betas.copy] <- NA
    
    index <- which(is.na(betas), arr.ind = TRUE) 
    if (length(index)>0){
      message("Replace ",length(index)," missing values with rowmeans")
      betas[index] <- rowMeans(betas, na.rm = TRUE)[index[, "row"]] 
    }
    save(betas, file=out_file)
    save(classes, file=paste0(out_file, "_classes"))
  }
}

readGRM <- function(rootname)
{
  bin.file.name <- paste(rootname, ".grm.bin", sep="")
  n.file.name <- paste(rootname, ".grm.N.bin", sep="")
  id.file.name <- paste(rootname, ".grm.id", sep="")
  
  cat("Reading IDs\n")
  id <- read.table(id.file.name)
  n <- dim(id)[1]
  cat("Reading GRM\n")
  bin.file <- file(bin.file.name, "rb")
  grm <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(bin.file)
  cat("Reading N\n")
  n.file <- file(n.file.name, "rb")
  N <- readBin(n.file, n=n*(n+1)/2, what=numeric(0), size=4)
  close(n.file)
  
  cat("Creating data frame\n")
  l <- list()
  for(i in 1:n)
  {
    l[[i]] <- 1:i
  }
  col1 <- rep(1:n, 1:n)
  col2 <- unlist(l)
  grm <- data.frame(id1=col1, id2=col2, N=N, grm=grm)	
  
  ret <- list()
  ret$grm <- grm
  ret$id <- id
  return(ret)
}

makeGRMmatrix <- function(grm)
{
  mat <- diag(nrow(grm$id))
  mat[upper.tri(mat, diag=TRUE)] <- grm$grm$grm
  mat <- t(mat)
  nsnpvec <- subset(grm$grm, id1 != id2)$N
  mat[upper.tri(mat, diag=FALSE)] <- nsnpvec
  rownames(mat) <- grm$id$V2
  colnames(mat) <- grm$id$V2
  return(mat)
}

adjust.relatedness.fast.1 <- function(x, covs, kin, eig, transform, quiet=TRUE)
{
  x[!is.finite(x)] <- mean(x, na.rm=T)
  
  if(transform == "transformed"){
    d <- data.frame(X=rntransform(x), covs)
  } else {
    d <- data.frame(X=x, covs)
  }
  rownames(d) <- colnames(kin)
  form <- as.formula(paste0("X ~ ", paste(names(d)[-1], collapse=" + ")))
  d$X <- residuals(lm(form, d, na.action = na.exclude))
  p_out <- try(polygenic(X, data=d, kinship.matrix=kin, eigenOfRel=eig, quiet=quiet))
  
  iter <- 1
  while(class(p_out) == "try-error" & iter < 5)
  {
    message("trying again...")
    iter <- iter + 1
    p_out <- try(polygenic(X, data=d, kinship.matrix=kin, eigenOfRel=eig, quiet=quiet))
  }
  if(class(p_out) == "try-error")
  {
    message("giving up, just using fixed effects model")
    a <- d$X
  } else {
    if(transform == "Yes"){
      a <- as.numeric(rntransform(p_out$grresidualY))
    }
    else{
      a <- as.numeric(p_out$grresidualY)
    }
  }
  return(list(x=a, cl=class(p_out)))
}

adjust.relatedness.serial <- function(B, covs, kin, eig, transform)
{
  cl <- array(0, nrow(B))
  for(i in 1:nrow(B))
  {
    message("Probe ",i, " of ", nrow(B))
    out <- adjust.relatedness.fast.1(B[i,], covs, kin, eig, transform)
    B[i, ] <- out$x
    cl[i] <- out$cl
  }
  return(list(x=B, cl=cl))
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

rntransform <- function(x)
{
  out <- rank(x) - 0.5
  out[is.na(x)] <- NA
  mP <- 0.5/max(out, na.rm = T)
  out <- out/(max(out, na.rm = T) + 0.5)
  out <- scale(qnorm(out))
  out
}

adjust.relatedness <- function(B, covs, kin, eig, mc.cores=mc.cores, transform)
{
  
  l1 <- get.index.list(nrow(B), mc.cores)
  l <- lapply(l1, function(ii)
  {
    res <- mclapply(ii, function(i)
    {
      message("Probe ", i, " of ", nrow(B))
      out <- adjust.relatedness.fast.1(B[i,], covs, kin, eig, transform)
    }, mc.cores=mc.cores, mc.preschedule=FALSE)
    a <- do.call(rbind, lapply(res, function(x) x$x))
    b <- sapply(res, function(x) x$cl)
    return(list(x=a, cl=b))
  })
  x <- do.call(rbind, lapply(l, function(x) x$x))
  cl <- unlist(lapply(l, function(x) x$cl))
  rownames(x) <- rownames(B)
  colnames(x) <- colnames(B)
  return(list(x=x, cl=cl))
}

main()
