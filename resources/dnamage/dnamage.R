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

# Calculate the summary statistics
cal.M.stats <- function(AgePredTable, PhenVal, AgeValid, SexValid, ClockName){
  # pre-phase
  temp = merge(AgePredTable, PhenVal, by.x = 'IID', by.y = 'IID')

  M1 = lm(PredAge ~ ., subset(temp, select = -c(Smoking, IID)), na.action=na.exclude)
  M2 = lm(PredAge ~ ., subset(temp, select = -c(IID)), na.action=na.exclude)

  coefM1 = as.data.frame(summary(M1)$coefficients)
  coefM2 = as.data.frame(summary(M2)$coefficients)

  stats_name = c("Min","Mean", "Median", "Max", "SD", "Rsquared", "AdjRsquared", 
               "Intercept", "SEIntercept","PvalIntercept")

  statsM0 = c(min(temp$PredAge), mean(temp$PredAge), median(temp$PredAge), max(temp$PredAge), sd(temp$PredAge),rep(NA,5))
  statsM1 = c(min(residuals(M1)), mean(residuals(M1)),median(residuals(M1)),max(residuals(M1)), summary(M1)$sigma,
              summary(M1)$r.squared, summary(M1)$adj.r.squared,
              coefM1["(Intercept)","Estimate"], coefM1["(Intercept)","Std. Error"],
              coefM1["(Intercept)","Pr(>|t|)"])
  statsM2 = c(min(residuals(M2)), mean(residuals(M2)),median(residuals(M2)),max(residuals(M2)), summary(M2)$sigma,
              summary(M2)$r.squared, summary(M2)$adj.r.squared,
              coefM2["(Intercept)","Estimate"], coefM2["(Intercept)","Std. Error"],
              coefM2["(Intercept)","Pr(>|t|)"])

  if (AgeValid){
    stats_name = c(stats_name,"CorAge_numeric","SEAge_numeric","PvalAge_numeric")
    
    statsM0Age = cor.test(temp$PredAge,temp$Age, method = "p", na.action=na.exclude)
    statsM0 = c(statsM0, statsM0Age$estimate, unname(sqrt((1 - statsM0Age$estimate^2)/statsM0Age$parameter)),
                statsM0Age$p.value)
    statsM1 = c(statsM1, coefM1["Age_numeric","Estimate"], coefM1["Age_numeric","Std. Error"],
                coefM1["Age_numeric","Pr(>|t|)"])
    statsM2 = c(statsM2, coefM2["Age_numeric","Estimate"], coefM2["Age_numeric","Std. Error"],
                coefM2["Age_numeric","Pr(>|t|)"])
  } 


  if (SexValid){
    stats_name = c(stats_name,"CorSex_factor","SESex_factor","PvalSex_factor")
    
    statsM0Sex = t.test(PredAge ~ Sex_factor, temp, na.action=na.exclude)
    statsM0 = c(statsM0, statsM0Sex$statistic, statsM0Sex$stderr, statsM0Sex$p.value)
    statsM1 = c(statsM1, coefM1["Sex_factorM","Estimate"], coefM1["Sex_factorM","Std. Error"],
                coefM1["Sex_factorM","Pr(>|t|)"])
    statsM2 = c(statsM2, coefM2["Sex_factorM","Estimate"], coefM2["Sex_factorM","Std. Error"],
                coefM2["Sex_factorM","Pr(>|t|)"])

 }

  stats_name = c(stats_name,"CorSmoking","SESmoking","PvalSmoking")
  statsM0Smok = cor.test(temp$PredAge,temp$Smoking, method = "p", na.action=na.exclude)
  statsM0 = c(statsM0, statsM0Smok$estimate, unname(sqrt((1 - statsM0Smok$estimate^2)/statsM0Smok$parameter)),
              statsM0Smok$p.value)
  statsM1 = c(statsM1, rep(NA,3))
  statsM2 = c(statsM2, coefM2["Smoking","Estimate"], coefM2["Smoking","Std. Error"],
              coefM2["Smoking","Pr(>|t|)"])

  sigM1 <- coefM1[coefM1$`Pr(>|t|)` < 0.05, ]
  sigM2 <- coefM2[coefM2$`Pr(>|t|)` < 0.05, ] 
  
  if (nrow(sigM1) > 0 ){
    print(paste0(ClockName,":Coefficients table of Module1 with Sig variables========="))
    print(sigM1)
  }

  if (nrow(sigM1) > 0 ){
    print(paste0(ClockName,":Coefficients table of Module2 with Sig variables========="))
    print(sigM1)
  }
  
  out = data.frame(StatsValue = stats_name,
                  StatsPredAge = statsM0,
                  StatsM1 = statsM1,
                  StatsM2 = statsM2)
  colnames(out) = c("StatsValue", ClockName, paste0(ClockName,"M1"),paste0(ClockName,"M2"))
  rownames(out) = out$StatsValue
  out = out[,-1]

  return(out)
}

# Calculate residual divided by SD
generate.aar <- function(AgePredTable, PhenVal, ClockName){
  
  phen = subset(PhenVal, select = -c(IID))
  temp = merge(AgePredTable,PhenVal, by.x = 'IID', by.y = 'IID')
  temp=temp[match(AgePredTable$IID, temp$IID),]

  if(ncol(phen) == 1 && colnames(phen) == "Smoking"){
    AgePredTable$PredAgeSD = AgePredTable$PredAge/sd(AgePredTable$PredAge)
  } else {
    module = residuals(lm(PredAge ~ ., data = subset(temp, select = -c(Smoking, IID)), na.action=na.exclude))
    AgePredTable$PredAgeSD = module/sd(module)
  }
  
  moduless = residuals(lm(PredAge ~ ., data = subset(temp, select = -c(IID)), na.action=na.exclude))
  AgePredTable$PredAgessSD = moduless / sd(moduless)
  AgePredTable = subset(AgePredTable, select = c(IID, PredAge, PredAgeSD, PredAgessSD))
  colnames(AgePredTable) = c('IID', ClockName, paste0(ClockName, 'SD'), paste0(ClockName, 'ssSD'))
  
  return(AgePredTable)
  
}

# Plot
age.plot = function(AgePredTable, PhenVal, AgeValid, ClockNames, SD){

  if (AgeValid == T){
    par(mfrow=c(3,3))
    temp = merge(AgePredTable, PhenVal[,c('IID', 'Age_numeric')], by.x = 'IID', by.y = 'IID')
    temp = temp[match(AgePredTable$IID, temp$IID),]
    cAge=temp$Age_numeric
  } else {par(mfrow=c(3,2))}
  
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

    if ( length(pAge) > 5000 ){
		  shapiroscore = signif(as.numeric(shapiro.test(sample(pAge, 3000))[2]),2)
    }else{
      shapiroscore = signif(as.numeric(shapiro.test(pAge)[2]),2)
    }
	
    
    hist(pAge, xlab="", main=paste(cname,"\n(N=", length(which(!is.na(pAge))),"; Mean=", round(mean(pAge, na.rm=T),3),")",sep=""), cex.main=1)
    abline(v=mean(pAge, na.rm=T)-SD*sd(pAge, na.rm=T), lty=2)
    abline(v=mean(pAge,na.rm=T)+SD*sd(pAge, na.rm=T), lty=2)
    
    qqnorm(pAge, main=paste(cname, "\n(N=", length(which(!is.na(pAge))),"; shapiroP=",shapiroscore,")",sep=""),cex.main=1)
    qqline(pAge)
  }
}

# Collect stistics baed on the sex
cal.Sex.stats <- function(CorTable, ClockName){
    age_stats = c(min(CorTable[,c(ClockName)], na.rm=T), 
                  mean(CorTable[,c(ClockName)], na.rm=T), 
                  median(CorTable[,c(ClockName)], na.rm=T), 
                  max(CorTable[,c(ClockName)], na.rm=T), 
                  sd(CorTable[,c(ClockName)], na.rm=T))
    age_stats_F = c(min(CorTable[CorTable$Sex_factor == "F", ClockName], na.rm=T), 
                    mean(CorTable[CorTable$Sex_factor == "F", ClockName], na.rm=T), 
                    median(CorTable[CorTable$Sex_factor == "F", ClockName], na.rm=T), 
                    max(CorTable[CorTable$Sex_factor == "F", ClockName], na.rm=T), 
                    sd(CorTable[CorTable$Sex_factor == "F", ClockName], na.rm=T))
    age_stats_M = c(min(CorTable[CorTable$Sex_factor == "M", ClockName], na.rm=T), 
                    mean(CorTable[CorTable$Sex_factor == "M", ClockName], na.rm=T), 
                    median(CorTable[CorTable$Sex_factor == "M", ClockName], na.rm=T), 
                    max(CorTable[CorTable$Sex_factor == "M", ClockName], na.rm=T), 
                    sd(CorTable[CorTable$Sex_factor == "M", ClockName], na.rm=T))
              
    outlist = data.frame(StatsValue = c("min", "mean", "median", "max", "sd"),
                         All = age_stats,
                         F = age_stats_F,
                         M = age_stats_M)
    colnames(outlist) = c("StatsValue", ClockName, paste0(ClockName,"_F"), paste0(ClockName,"_M"))
    return(outlist)
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
  
  if (cellcount_file == 'NULL') {
    message("No predicted cell count matrix provided.")
  } else {
    cellcount <- read.table(cellcount_file, header = T, stringsAsFactors=FALSE)
    cellcount <- subset(cellcount, select=-c(Treg))
    m <- match(fam[,2], cellcount[,"IID"])
    cellcount <- cellcount[m,]
    phen_value <- merge(phen_value, cellcount,  by.x = 'IID', by.y = 'IID', all.x = TRUE)
  }
  phen_value <- phen_value[match(fam[,2], phen_value$IID),]
  
  pdf(paste0(age_plot, '.pdf'), width=12, height=12)
  name_sumstats <- c()
  sumstats <- c()
  
  message(" ")
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
    modulestatsDNA <- cal.M.stats(AgePredTable=dnampred, PhenVal=phen_value, AgeValid=age_valid, SexValid=sex_valid, ClockName = 'DNAmAge')
    # age acceleration residual prediction
    dnampred <- generate.aar(AgePredTable=dnampred, PhenVal=phen_value, ClockName = 'DNAmAge')
    # plot
    age.plot(AgePredTable=dnampred, PhenVal=phen_value, AgeValid=age_valid, ClockNames=colnames(dnampred)[-1], SD=SD)
    dnadensity <- density(dnampred[,2])
    dna_valid <- TRUE
  }
  
  message(" ")
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
    modulestatsPheno <- cal.M.stats(AgePredTable=phenpred, PhenVal=phen_value, AgeValid=age_valid, SexValid=sex_valid, ClockName = 'PhenoAge')
    # age acceleration residual prediction
    phenpred <- generate.aar(AgePredTable=phenpred, PhenVal=phen_value, ClockName='PhenoAge')
    # plot
    age.plot(AgePredTable=phenpred, PhenVal=phen_value, AgeValid=age_valid, ClockNames=colnames(phenpred)[-1], SD)
    phendensity <- density(phenpred$PhenoAge)
    phen_valid <- TRUE
  }
  
  message(" ")
  message("Predicting, adjusting and standardizing DunedinPACE")########################################
  suppressWarnings(PredAgeDun <- PACEProjector(norm.beta, proportionOfProbesRequired=0.7))
  paceresult <- data.frame("IID" = names(PredAgeDun[['DunedinPACE']]), "PredAge" = PredAgeDun[['DunedinPACE']])
  
  if (inherits(paceresult, 'data.frame') & inherits(paceresult$PredAge, 'numeric')) {
    pacepred <- paceresult
    # summary on module statistic
    name_sumstats <- c(name_sumstats, 'DunedinPACE')
    sumstats <- rbind(sumstats, cal.stats(AgePredTable=pacepred, PhenVal=phen_value, AgeValid=age_valid, SexValid=sex_valid, SD))
    modulestatsDun <- cal.M.stats(AgePredTable=pacepred, PhenVal=phen_value, AgeValid=age_valid, SexValid=sex_valid, ClockName = 'DunedinPACE')
    # age acceleration residual prediction
    pacepred <-  generate.aar(AgePredTable=pacepred, PhenVal=phen_value, ClockName='DunedinPACE')
    # plot
    age.plot(AgePredTable=pacepred, PhenVal=phen_value, AgeValid=age_valid, ClockName=colnames(pacepred)[-1], SD)
    pacedensity <- density(pacepred$DunedinPACE)
    pace_valid <- TRUE
  } else {
    message("ERROR: Failure on prediction on DunedinPACE by using PACEProjector function")
  }

  message(" ")
  message("Checking the output files.")########################################
  message("Check if Age_numeric is valid: ", age_valid)
  message("Check if Sex_factor is valid: ", sex_valid)
  message("Check if DNAmAge, DNAmAgeSD and DNAmAgessSD have been generated: ", dna_valid)
  message("Check if PhenoAge, PhenoAgeSD and PhenoAgessSD have been generated: ", phen_valid)
  message("Check if DunedinPACE, DunedinPACESD and DunedinPACEssSD have been generated: ", pace_valid)
  
  message(" ")
  message("Outputing statistic table and age acceleration files.")########################################
  # close the pdf file
  dev.off()
  message(paste0("Scatterplot plots saved to ", age_plot, ".pdf"))

  # statistic table of models
  rownames(sumstats) <- name_sumstats
  colnames(sumstats) <- c("SampleSize", "mean", "SD", "cor", "SmokeEst", "AgeEst", "SexEst",
                          "SmokeSE", "AgeSE", "SexSE", "SmokePr(>|t|)", "AgePr(>|t|)", "Pr(>|t|)", 
                          "NPostFilter", "meanPostFilter", "sdPostFilter", "corPostFilter",
                          "SmokeEstPostFilter", "SmokeSEPostFilter", "SmokePPostFilter",
                          "AgeEstPostFilter", "AgeSEPostFilter", "AgePPostFilter", 
                          "SexEstPostFilter", "SexSEPostFilter", "SexPPostFilter")
  write.csv(sumstats, file = paste0(age_stats, ".csv"))	
  message(paste0("Basic statistic table of vairables saved to ", age_stats, ".csv"))
  
  modulestats <- cbind(modulestatsDNA, modulestatsPheno, modulestatsDun)
  write.csv(modulestats, file = paste0(age_stats, "_models.csv"), row.names = TRUE, quote = FALSE, na = "NA")	
  message(paste0("Statistic table of modules saved to ", age_stats, "_models.csv"))

  # output age acceleration files and correlation matrix files
  colnames(fam) <- c('FID', 'IID')
  cortable <- phen_value
  if (age_valid){sex_colnames <- c('Age_numeric')}else{sex_colnames <- c()}

  if (exists('dnampred')) {
    sex_colnames <- c(sex_colnames, 'DNAmAge', 'DNAmAgeSD', 'DNAmAgessSD')
    fam <- merge(fam, dnampred[,c(1,3,4)], by.x = 'IID', by.y = 'IID', all.x=TRUE)
    cortable <- merge(cortable, dnampred, by.x = 'IID', by.y = 'IID', all.x=TRUE)
  }
  
  if (exists('phenpred')) {
    sex_colnames <- c(sex_colnames, 'PhenoAge', 'PhenoAgeSD', 'PhenoAgessSD')
    fam <- merge(fam, phenpred[,c(1,3,4)], by.x = 'IID', by.y = 'IID', all.x=TRUE)
    cortable <- merge(cortable, phenpred, by.x = 'IID', by.y = 'IID', all.x=TRUE)
  }

  if (exists('pacepred')) {
    sex_colnames <- c(sex_colnames, 'DunedinPACE', 'DunedinPACESD', 'DunedinPACEssSD')
    fam <- merge(fam, pacepred[,c(1,3,4)], by.x = 'IID', by.y = 'IID', all.x=TRUE)
    cortable <- merge(cortable, pacepred, by.x = 'IID', by.y = 'IID', all.x=TRUE)
  }
  
  fam <- fam[,c(2,1,seq(3,ncol(fam)))]
  write.table(fam, file=paste0(out_file, ".txt"), row=F, col=T, qu=F)
  message(paste0("Age prediction table save to", out_file, ".txt"))
  save(age_valid, sex_valid, dna_valid, phen_valid, pace_valid, cortable, fam, file=paste0(out_file, ".RData"))
  message(paste0("Variables and age prediction for plotting save to", out_file, ".RData"))

  # correlation and sd table for phenotypes matrix 
  if (sex_valid) {
    corstat <- subset(cortable, select = -c(IID, Sex_factor))
  } else {
    corstat <- subset(cortable, select = -c(IID))
  }
  std <- function(x) round(sd(x)/sqrt(length(x)), digits=4)
  corstats = as.data.frame(cor(corstat))
  corstats$SD = apply(corstat,2,std)
  corstats$mean = apply(corstat,2,mean)
  write.csv(corstats, file = paste0(age_stats, "_corrsd.csv"))
  message(paste0("Correlation table among epi age preidcitons save to", age_stats, "_corrsd.csv"))

  # sex information
  if (sex_valid){
    sex_out = data.frame(StatsValue = c("min", "mean", "median", "max", "sd"))
    for (cname in sex_colnames){
      sex_stats_table = cal.Sex.stats(cortable, cname)
      sex_out = merge(sex_out, sex_stats_table, by = "StatsValue", all.x = TRUE)
    }
    rownames(sex_out) = sex_out$StatsValue
    write.csv(sex_out, file = paste0(age_stats, "_sex.csv"))
    message(paste0("Statistic table of epi age preidction based on sex save to ", age_stats, "_sex.csv"))
  }
}
main()


