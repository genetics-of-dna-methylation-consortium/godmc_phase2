
###############################discription###################################

# This script is dedicated to generating smoking residuals for performing GWAS on smoking.
# The phenotypes used to generate smoking residuals are age, sex, and cell count.

# If any errors occur while running this script, 
# please don't hesitate to reach out to the developer for assistance at s.w.wang@exeter.ac.uk.

###############################mian function###################################


main = function () {
  library(tibble)
  arguments <- commandArgs(T)
  
  cellcount_file <- arguments[1]
  covs_file <- arguments[2]
  fam_file <- arguments[3]
  out <- arguments[4]
  stats_out <- arguments[5]



  message("Generate phenotype for GWAS on smoking.")
  
  # Read in files
  covs <- read.table(covs_file, header=T, stringsAsFactors=FALSE)
  fam <- read.table(fam_file, stringsAsFactors=FALSE)[,1:2]
  smok <- read.table(paste0(out, ".txt"), header = T, stringsAsFactors = FALSE)
  covs <- covs[match(fam[,2], covs$IID),]
  smok <- smok[match(fam[,2], smok$IID),]

  # Check if the age and sex are validate
  age_index <- grep("^age_numeric$", names(covs), ignore.case=TRUE)
  count_age <- length(unique(covs[,age_index]))
  sex_index <- grep("^Sex_factor$", names(covs), ignore.case = TRUE)
  count_sex <- length(unique(covs[, sex_index]))

  
  age_vaild = FALSE
  sex_vaild = FALSE
  if(length(age_index) != 1 | length(sex_index) != 1){
    message("There should be only one column in the covariate file called 
         'Age_numeric' and 'Sex_factor', regardless of case.")
    valid_vector <- c('IID')
  }else if (count_age == 1 & count_sex == 1){
    message("Age variable is contant, which should not be considered.")
    message("Sex variable has 1 level.")
    message("There is no sex factor but age numeric should be considered.")
    valid_vector <- c('IID')
  }else if (count_age == 1 & count_sex != 1){
    message("Age variable is contant, which should not be considered.")
    message("Sex variable has ", count_sex, " levels.")
    sex_vaild = TRUE
    valid_vector <- c('IID', 'Sex_factor')
  }else if (count_age != 1 & count_sex == 1){
    message("Sex variable has one level, which should not be considered.")
    age_vaild = TRUE
    valid_vector <- c('IID', 'Age_numeric')
  }else {
    age_vaild = TRUE
    sex_vaild = TRUE
    valid_vector <- c('IID', 'Sex_factor', 'Age_numeric')
  }
  
  phen_value <- subset(covs, select = valid_vector) 
  phen_value <- merge(phen_value, smok[,c('IID', 'Smoking')], by.x="IID", by.y="IID", all.x=TRUE)

  # Check if the cellcount is available
  if (cellcount_file == 'NULL') {
    message("No predicted cell count matrix provided.")
  } else {
    cellcount <- read.table(cellcount_file, header = T, stringsAsFactors=FALSE)
    m <- match(fam[,2], cellcount[,"IID"])
    cellcount <- cellcount[m,]
    # exclude one of the cell types with highest number of zero or lowest median value
    nomissing_cell = as.integer(which(colSums(cellcount==0) == 0))
    missing_cell = order(colSums(cellcount==0), decreasing = TRUE)
    if (length(setdiff(missing_cell, nomissing_cell)) > 1){
      cellcount = subset(cellcount, select = -c(setdiff(missing_cell, nomissing_cell)[1]))
    } else {
      lowest_median_cell = which(colMedians(as.matrix(cellcount[,-1])) == min(colMedians(as.matrix(cellcount[,-1]))))
      cellcount = cellcount[, -(lowest_median_cell+1)]
    }
    phen_value <- merge(phen_value, cellcount, by.x = 'IID', by.y = 'IID', all.x = TRUE)
  }
  

  # calculate the residual model of smoking
  phen <- phen_value[match(fam[,2], phen_value$IID),]
  phen <- subset(phen, select = -c(IID))
  stats_name <- c("Min","Mean", "Median", "Max", "SD")
  smok_M0 <- c(min(phen$Smoking), mean(phen$Smoking), median(phen$Smoking), max(phen$Smoking), sd(phen$Smoking))
  if (ncol(phen) == 1 && colnames(phen) == 'Smoking'){
    message("There is no covariate to adjust for smoking residuals.")
    write.table(smok, file=paste0(out, ".smok.plink"), row=F, col=F, qu=F)
    smok_sum <- data.frame( stats_name, smok_M0)
  }else{
    SmokModel <- residuals(lm(phen$Smoking ~ ., phen, na.action=na.exclude))
    fam$residual <- SmokModel / sd(SmokModel)
    M1 <- lm(phen$Smoking ~ ., phen, na.action=na.exclude)
    coefM1 <- summary(M1)$coefficients
    stats_name <- c(stats_name, "Rsquared", "AdjRsquared", "Intercept", "SEIntercept","PvalIntercept")
    smok_sum <- data.frame( stats_name, smok_M0 = c(smok_M0, rep(NA, 5)), 
                            smok_M1 = c(min(SmokModel), mean(SmokModel), 
                            median(SmokModel), max(SmokModel), sd(SmokModel),
                            summary(M1)$r.squared, summary(M1)$adj.r.squared,
                            coefM1["(Intercept)","Estimate"], coefM1["(Intercept)","Std. Error"],
                            coefM1["(Intercept)","Pr(>|t|)"]))
    print("Summary table for smoking model==========================")
    print(summary(M1))
    write.table(fam, file=paste0(out, ".smok.plink"), row=F, col=F, qu=F)
  }

  write.table(smok_sum, file=paste0(stats_out, ".csv"), row.names=F, col.names=T, sep=",", quote=F)

  
}

main()

