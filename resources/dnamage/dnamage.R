###############################discription###################################

#This script, dedicated to age acceleration prediction using DunedinPACE, 
#Horvath's clock, and Levine's clock, is neatly organized into four key sections:

#1.Resources and Materials: Details the necessary resources for script execution.
#2.Packages: Lists the essential packages to be utilized.
#3.Sub-functions: Contains several sub-functions that will be employed within the main function.
#4.Main Function: Focuses on the primary function, generating output files from the DNAm matrix.

# If any errors occur while running this script, 
# please don't hesitate to reach out to the developer for assistance at s.w.wang@exeter.ac.uk.

###############################resources###################################

## USE:
## http://labs.genetics.ucla.edu/horvath/dnamage/TUTORIAL1.pdf
## http://labs.genetics.ucla.edu/horvath/dnamage/probeAnnotation21kdatMethUsed.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/datMiniAnnotation27k.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/AdditionalFile3.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/StepwiseAnalysis.txt
## http://labs.genetics.ucla.edu/horvath/dnamage/NORMALIZATION.R
## DunedinPACE calculation
## https://github.com/danbelsky/DunedinPACE

#################################packages####################################
suppressPackageStartupMessages(library(DunedinPACE))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RPMM))
suppressPackageStartupMessages(library(impute))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(vioplot))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))

#################################Sub-functions################################


# Calculate DNAmAge
dnam.age <- function(x, DatClock) {
  set.seed(1)
  # function to model age non-linearly
  anti.trafo <- function(x, adult.age=20) {
    ifelse(x<0,
           (1+adult.age)*exp(x)-1,
           (1+adult.age)*x+adult.age)
  }
  
  nsamples <- dim(x)[[2]]
  datMethUsedNormalized = imputation.knn(x = x)
  gc()
  selectCpGsClock=is.element(dimnames(datMethUsedNormalized)[[1]], as.character(DatClock$CpGmarker[-1]))
  if ( nsamples>1 ) {
    datMethClock=datMethUsedNormalized[selectCpGsClock,]
    predictedAge=as.numeric(anti.trafo(DatClock$CoefficientTraining[1]+t(as.matrix(datMethClock)) %*% as.numeric(DatClock$CoefficientTraining[match(rownames(datMethClock), DatClock$CpGmarker)])))
  }else if ( nsamples==1 ) {
    predictedAge=as.numeric(anti.trafo(DatClock$CoefficientTraining[1]+sum(datMethUsedNormalized *as.numeric(DatClock$CoefficientTraining[match(colnames(datMethClock), DatClock$CpGmarker)]))))
  }
  
  datout = data.frame(IID = colnames(x), PredAge = predictedAge)
  
  return(datout)
}

# Calculate PhenoAge
pheno.age <- function(beta, PhenoAgeCoeff){
  intercept <- PhenoAgeCoeff$Weight[1]
  coeffs <- PhenoAgeCoeff$Weight[-1]
  probeOverlap <- match(PhenoAgeCoeff$CpG[-1], rownames(beta))
  coeffs <- coeffs[!is.na(probeOverlap)]
  probeOverlap <- probeOverlap[!is.na(probeOverlap)]
  beta <- beta[probeOverlap,]
  imputed_beta <- imputation.knn(beta)
  imputed_dat <- as.matrix(imputed_beta)
  predAge <- intercept+coeffs %*% imputed_dat
  return(t(predAge))
}

# Calculate the summary statistics
cal.stats <- function(AgePredTable, PhenVal, AgeValid, SexValid, SD){
  # pre-phase
  temp = merge(AgePredTable, PhenVal, by.x = 'IID', by.y = 'IID')
  temp = subset(temp, select = -c(IID))
  cor <- NA
  subcor <- NA
  sumvar <- c('Smoking')
  # basic statistic values and identify outliers
  mu <- mean(temp$PredAge, na.rm=T)
  sigma <- sd(temp$PredAge, na.rm=T)
  keepindex <- which(temp$PredAge < mu + SD*sigma & temp$PredAge > mu - SD*sigma) 
  subtemp <- temp[keepindex,]
  # ckeck age and sex
  if (AgeValid){
    sumvar <- c(sumvar, 'Age_numeric') 
    cor <- cor(temp$PredAge, temp$Age_numeric, use = "p") 
    subcor <- cor(subtemp$PredAge, subtemp$Age_numeric, use = "p")
  } 
  if (SexValid) sumvar <- c(sumvar, 'Sex_factorM') 
  nrep <- 3 - length(sumvar)
  outtemp <- c(sum(!is.na(temp$PredAge)), mu, sigma, cor)
  # module for all the sample
  moduleall <- lm(PredAge ~ ., data=temp, na.action=na.exclude)
  outtemp <- c(outtemp, c(summary(moduleall)$coefficients[sumvar,c(1,2,4)], rep(NA, 3*nrep)))
  
  outtemp <- c(outtemp, sum(!is.na(subtemp$PredAge)), mean(subtemp$PredAge, na.rm=T), sd(subtemp$PredAge, na.rm=T), subcor)
  # module for filtered sample
  modulesub <- lm(PredAge ~ ., data=subtemp, na.action=na.exclude)
  outtemp <- c(outtemp, c(summary(modulesub)$coefficients[sumvar,c(1,2,4)], rep(NA, 3*nrep)))
  
  return(outtemp)
}


# Calculate residual divided by SD
generate.aar <- function(AgePredTable, PhenVal, ClockName){
  
  phen = subset(PhenVal, select = -c(IID))
  temp = merge(AgePredTable,PhenVal, by.x = 'IID', by.y = 'IID')
  if(ncol(phen) == 1 && colnames(phen) == "Smoking"){
    AgePredTable$PredAgeSD = AgePredTable$PredAge/sd(AgePredTable$PredAge)
  } else {
    module = residuals(lm(temp$PredAge ~ ., subset(temp, select = -c(Smoking, IID)), na.action=na.exclude))
    AgePredTable$PredAgeSD = module/sd(module)
  }
  
  moduless = residuals(lm(temp$PredAge ~ ., subset(temp, select = -c(IID)), na.action=na.exclude))
  AgePredTable$PredAgessSD = moduless / sd(moduless)
  colnames(AgePredTable) = c('IID', ClockName, paste0(ClockName, 'SD'), paste0(ClockName, 'ssSD'))
  
  return(AgePredTable)
  
}

# Plot
age.plot = function(AgePredTable, PhenVal, AgeValid, ClockNames, SD){

  if (AgeValid == T){
    par(mfrow=c(3,3))
    temp = merge(AgePredTable, PhenVal[,c('IID', 'Age_numeric')], by.x = 'IID', by.y = 'IID')
    cAge=temp$Age_numeric
  } else {
    par(mfrow=c(3,2))
  }
  
  for (cname in ClockNames) {
    if (AgeValid == T){
      pAge=temp[,which(colnames(temp) == cname)]
      message("The correlation between predicted ", cname, " and actual age is ", cor(cAge, pAge,use="pair"))
      plot(cAge, pAge, cex.main=1, cex=0.7, 
           main=paste("correlation between ", cname, "\nand actual age=",signif(cor(cAge, pAge,use="pair"),2), sep=""), 
           xlab = "Chronological Age (years)", ylab = cname)
      abline(lm(pAge ~ cAge))
    }else{
      pAge = AgePredTable[,which(colnames(AgePredTable)==cname)]
      message("No correlation between predicted ", cname, " and actual age.")
    }
    
    hist(pAge, xlab="", main=paste(cname,"\n(N=", length(which(!is.na(pAge))),"; Mean=", round(mean(pAge, na.rm=T),3),")",sep=""), cex.main=1)
    abline(v=mean(pAge, na.rm=T)-SD*sd(pAge, na.rm=T), lty=2)
    abline(v=mean(pAge,na.rm=T)+SD*sd(pAge, na.rm=T), lty=2)
    
    qqnorm(pAge, main=paste(cname, "\n(N=", length(which(!is.na(pAge))),"; shapiroP=",signif(as.numeric(shapiro.test(pAge)[2]),2),")",sep=""),cex.main=1)
    qqline(pAge)
}

}

# Imputation
imputation.knn <- function(x){
  # assumes you are only going to impute a few hundred sites. 
  nsamples = ncol(x)
  datMethUsed = t(x)
  dimnames1 = dimnames(datMethUsed)
  noMissingPerSample = rowSums(is.na(datMethUsed))
  
  if(nsamples > 1 & max(noMissingPerSample) > 0){
    datMethUsed = as.data.frame(t(impute.knn(datMethUsed)$data))
    colnames(datMethUsed) = dimnames1[[1]]
    return(datMethUsed)
  } else {
    return(x)
  }
}

# Matrix Scatter plot
# lower panel - correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=3)
  txt <- paste0("cor: ", r)
  text(0.5, 0.5, txt, cex=1.5)
}

# upper panel - scatter plots
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = "grey50")
}
  
# diagnal panel - sd calculation
std <- function(x) round(sd(x)/sqrt(length(x)), digits=4)
panel.se <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  se = round(std(x), digits = 3)
  txt <- paste0("se:", se)
  mtext(txt, side = 1, line = -1.5)
}


main <- function()
{
  # Args
  arguments <- commandArgs(T)
  beta_file <- arguments[1]
  cov_file <- arguments[2]
  fam_file <- arguments[3]
  out_file <- arguments[4]
  age_plot <-arguments[5]
  SD <- as.numeric(arguments[6])
  age_stats <- arguments[7]
  smoking_file <- arguments[8]
  cellcount_file <- arguments[9]
  
  
  message("Getting probe parameters")#######################################
  #dnamage_probeAnnotation=read.csv("resources/dnamage/probeAnnotation21kdatMethUsed.csv.gz")
  dnamage_datclock=read.csv("resources/dnamage/AdditionalFile3.csv.gz")
  # this is a table of coefficients taken from the orginial PhenoAge manuscript
  phenoage_coeff=read.csv("resources/dnamage/PhenoAgeCoeff.csv", stringsAsFactors = FALSE, fill = TRUE) 
  
  
  message("Reading in data and matching up samples across files")#######################################
  covs <- read.table(cov_file, header=T)
  fam <- read.table(fam_file)[,c(1,2)]
  
  # assumes that everyone in the genetic data has methylation & covariate data
  covs <- covs[match(fam[,2], covs[,"IID"]),]
  smoking <- read.table(smoking_file, header = T)[,c('IID', 'Smoking')]
  smoking <- smoking[match(fam[,2], smoking[,"IID"]),]
  load(beta_file)
  norm.beta <- norm.beta[, match(fam[,2], colnames(norm.beta))]
  message(paste(nrow(fam), "samples with genetic data matched to methylation data"))
  
  
  message("Checking covariates")#######################################
  # Check if the age and sex are valid for predicting biological age
  age_index <- grep("^age_numeric$", names(covs), ignore.case=TRUE)
  count_age <- length(unique(covs[, age_index]))
  sex_index <- grep("^Sex_factor$", names(covs), ignore.case = TRUE)
  count_sex <- length(unique(covs[, sex_index]))
  
  age_valid = FALSE
  sex_valid = FALSE
  if(length(age_index) != 1 | length(sex_index) != 1){
    message("There should be only one column in the covariate file called 
         'Age_numeric' and 'Sex_factor'. Neither variable will be adjusted for.")
    valid_vector <- c('IID')
  }else if (count_age == 1 & count_sex == 1){
    message("Age variable is constant, so will not be adjusted for.")
    message("Sex variable only has 1 level, so will not be adjusted for.")
    valid_vector <- c('IID')
  }else if (count_age == 1 & count_sex != 1){
    message("Age variable is constant, so will not be adjusted for.")
    message("Sex variable has ", count_sex, " levels.")
    sex_valid = TRUE
    valid_vector <- c('IID', 'Sex_factor')
  }else if (count_age != 1 & count_sex == 1){
    message("Sex variable has one level, so will not be adjusted for.")
    age_valid = TRUE
    valid_vector <- c('IID', 'Age_numeric')
  }else {
    age_valid = TRUE
    sex_valid = TRUE
    valid_vector <- c('IID', 'Sex_factor', 'Age_numeric')
  }
  
  if(length(age_index) != 1){
    names(covs)[age_index] <- "Age_numeric"
  }  
  
  if(length(sex_index) != 1){
    names(covs)[sex_index] <- "Sex_factor"
  }
  
  
  phen_value <- subset(covs, select = valid_vector) 
  phen_value <- merge(phen_value, smoking, by.x="IID", by.y="IID", all.x=TRUE)
  if (sex_valid == T){
    cortable <- subset(phen_value, select = -c(Sex_factor))
  }else{
    cortable <- phen_value
  }
  
  if (cellcount_file == 'NULL') {
    message("No predicted cell count matrix provided.")
  } else {
    cellcount <- read.table(cellcount_file, header = T, stringsAsFactors=FALSE)
    cellcount <- subset(cellcount, select=-c(Treg))
    m <- match(fam[,2], cellcount[,"IID"])
    cellcount <- cellcount[m,]
    phen_value <- merge(phen_value, cellcount,  by.x = 'IID', by.y = 'IID', all.x = TRUE)
  }
  
  
  pdf(paste0(age_plot, '.pdf'), width=12, height=12)
  name_sumstats <- c()
  sumstats <- c()
  
  
  message("Predicting DNAmAge")#############################################
  # filter to clock probes
  DNA_overlap <- intersect(dnamage_datclock$CpGmarker[-1], rownames(norm.beta))
  DNA_beta <- norm.beta[match(DNA_overlap, rownames(norm.beta)),]
  nsamples <- dim(DNA_beta)[[2]]
  nprobes <- dim(DNA_beta)[[1]]
  cat(paste( "The methylation dataset used to calculate DNAmAge contains", nsamples, "samples (e.g. arrays) and ", nprobes, " probes.\n"))
  
  if (length(DNA_overlap) == 0 ) {
    message("ERROR: No overlapping CpGs: Can't proceed with age prediction for DNAmAge")
  } else if (nsamples == 0 | nprobes == 0){
    message("ERROR: There must be a data input error since there seem to be no either samples or zero probes.")
  } else{
    dnampred <- dnam.age(x = DNA_beta, DatClock = dnamage_datclock)
    # summary on module statistic
    name_sumstats <- c(name_sumstats, 'DNAmAge')
    sumstats<- rbind(sumstats, cal.stats(AgePredTable=dnampred, PhenVal=phen_value, AgeValid=age_valid, SexValid=sex_valid, SD))
    # age acceleration residual prediction
    dnampred <- generate.aar(AgePredTable=dnampred, PhenVal=phen_value, ClockName = 'DNAmAge')
    # plot
    age.plot(AgePredTable=dnampred, PhenVal=phen_value, AgeValid=age_valid, ClockNames=colnames(dnampred)[-1], SD=SD)
    dnadensity <- density(dnampred[,2])
  }
  
  
  message("Predicting Phenoage")############################################
  Pheno_overlap <- match(phenoage_coeff$CpG[-1], rownames(norm.beta))
  Pheno_overlap <- Pheno_overlap[!is.na(Pheno_overlap)]
  
  if (length(!is.na(Pheno_overlap)) < 50) {
    message("ERROR: Less overlapping CpGs: Can't proceed with age prediction for PhenoAge")
  } else {
    pre_phen <- pheno.age(beta = norm.beta, PhenoAgeCoeff = phenoage_coeff)
    phenpred <- data.frame(IID = row.names(pre_phen), PredAge = pre_phen[,1])
    # summary on module statistic
    name_sumstats <- c(name_sumstats, 'PhenoAge')
    sumstats<- rbind(sumstats, cal.stats(AgePredTable=phenpred, PhenVal=phen_value, AgeValid=age_valid, SexValid=sex_valid, SD))
    # age acceleration residual prediction
    phenpred <- generate.aar(AgePredTable=phenpred, PhenVal=phen_value, ClockName='PhenoAge')
    # plot
    age.plot(AgePredTable=phenpred, PhenVal=phen_value, AgeValid=age_valid, ClockNames=colnames(phenpred)[-1], SD)
    phendensity <- density(phenpred$PhenoAge)
  }
  
  
  message("Predicting, adjusting and standardizing DunedinPACE")########################################
  PredAgeDun <- try(PACEProjector(norm.beta, proportionOfProbesRequired=0.7), silent = TRUE)
  paceresult <- data.frame("IID" = names(PredAgeDun[['DunedinPACE']]), "PredAge" = PredAgeDun[['DunedinPACE']])
  
  if (inherits(paceresult, 'data.frame') & inherits(paceresult$PredAge, 'numeric')) {
    pacepred <- paceresult
    # summary on module statistic
    name_sumstats <- c(name_sumstats, 'DunedinPACE')
    sumstats <- rbind(sumstats, cal.stats(AgePredTable=pacepred, PhenVal=phen_value, AgeValid=age_valid, SexValid=sex_valid, SD))
    # age acceleration residual prediction
    pacepred <-  generate.aar(AgePredTable=pacepred, PhenVal=phen_value, ClockName='DunedinPACE')
    # plot
    age.plot(AgePredTable=pacepred, PhenVal=phen_value, AgeValid=age_valid, ClockName=colnames(pacepred)[-1], SD)
    pacedensity <- density(pacepred$DunedinPACE)
  } else {
    message("ERROR: Failure on prediction on DunedinPACE by using PACEProjector function")
  }
  
  
  message("Outputing age plots, statistic table and age acceleration files.")########################################
  if (age_valid){
    par(mfrow=c(2,1))
    legendname = c("Chronological Age")
    densitycolor = c("#ff595e")
    agedensity = density(covs$Age_numeric)
    plot(agedensity, xlab = "Age", col = "white", cex.main=1, cex=0.7,
         xlim = c(0, max(agedensity$x + 10)), 
         ylim = c(0, max(agedensity$y + 0.05)), 
         main = "Density plot of chronological age and predicted age")
    polygon(agedensity, col = alpha("#ff595e", 0.6))
    abline(v =mean(covs$Age_numeric), lty=2, col="#ff595e")
    if (exists('dandensity')) {
      polygon(dandensity, col = alpha("#ffca3a", 0.6))
      abline(v=mean(dnampred[,2]), lty=2, col= "#ffca3a")
      legendname = c(legendname, "DNAmAge")
      densitycolor = c(densitycolor, "#ffca3a")
    }
    if (exists('phendensity')){
      polygon(phendensity, col = alpha( "#8ac926", 0.6))
      abline(v=mean(phenpred[,2]), lty=2, col= "#8ac926")
      legendname = c(legendname, "PhenoAge")
      densitycolor = c(densitycolor, "#8ac926")
    }
    legend("topleft", legend = legendname, pch = 19, col = densitycolor, inset = 0.01)
    if (exists('pacedensity')){
      plot(pacedensity, xlab = "Pace of Aging", col = "white", cex.main=1, cex=0.7,
           main = "Density plot of DunedinPCAE")
      polygon(pacedensity, col = alpha("#1982c4", 0.6))
      abline(v=mean(pacepred[,2]), lty=2, col="#1982c4")
    }
  }
  dev.off()
  
  # statistic table of models
  rownames(sumstats) <- name_sumstats
  colnames(sumstats) <- c("SampleSize", "mean", "SD", "cor", "SmokeEst", "SmokeSE", "SmokeP",
                          "AgeEst", "AgeSE", "AgeP", "SexEst", "SexSE", "SexP", 
                          "NPostFilter", "meanPostFilter", "sdPostFilter", "corPostFilter",
                          "SmokeEstPostFilter", "SmokeSEPostFilter", "SmokePPostFilter",
                          "AgeEstPostFilter", "AgeSEPostFilter", "AgePPostFilter", 
                          "SexEstPostFilter", "SexSEPostFilter", "SexPPostFilter")
  write.csv(sumstats, file = paste0(age_stats, ".csv"))	
  
  # output age acceleration files and correlation matrix files
  colnames(fam) <- c('FID', 'IID')
  if (exists('dnampred')) {
    fam <- merge(fam, dnampred[,c(1,3,4)], by.x = 'IID', by.y = 'IID', all.x=TRUE)
    cortable <- merge(cortable, dnampred, by.x = 'IID', by.y = 'IID', all.x=TRUE)
  }
  
  if (exists('phenpred')) {
    fam <- merge(fam, phenpred[,c(1,3,4)], by.x = 'IID', by.y = 'IID', all.x=TRUE)
    cortable <- merge(cortable, phenpred, by.x = 'IID', by.y = 'IID', all.x=TRUE)
  }
  if (exists('pacepred')) {
    fam <- merge(fam, pacepred[,c(1,3,4)], by.x = 'IID', by.y = 'IID', all.x=TRUE)
    cortable <- merge(cortable, pacepred, by.x = 'IID', by.y = 'IID', all.x=TRUE)
  }
  fam <- fam[,c(2,1,seq(3,ncol(fam)))]
  write.table(fam, file=paste0(out_file, ".txt"), row=F, col=T, qu=F)
  
  # scatter plot of correlation matrix
  png(file = paste0(age_plot, "_correlation.png"), width=1400, height=800)
  
  cortable=subset(cortable, select=-c(IID))
  pairs(cortable, 
        lower.panel = panel.cor,
        upper.panel = upper.panel,
        diag.panel = panel.se,
        cex.labels= 1.2, gap = 0.3,
        main = "Scatterplot Matrix of Age Accelerations")
  dev.off()
  
  # correlation and sd table for phenotypes matrix 
  corstat = as.data.frame(cor(cortable))
  corstat$SD = apply(cortable,2,std)
  corstat$mean = apply(cortable,2,mean)
  write.csv(corstat, file = paste0(age_stats, "_corrsd.csv"))

  
  
} 

main()


