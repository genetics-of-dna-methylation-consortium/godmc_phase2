#############################################################################################################################################################################################################################
# Script to calculate MZ-EpiScores (epigenetic signature of MZ twins) # https://www.nature.com/articles/s41467-021-25583-7
# Developer contact: Jenny van Dongen, j.van.dongen@vu.nl
#############################################################################################################################################################################################################################

rm(list = ls(all = TRUE))
options(stringsAsFactors=FALSE)


##### R packages  #####
if(!require(glmnet)){
        install.packages("glmnet")
        library(glmnet)
}


if(!require(ROCR)){
  install.packages("ROCR")
  library(ROCR)
}

##### command arguments  #####
arguments <- commandArgs(T)
methylation=arguments[1]
fam_file <- arguments[2]
pheno_file = arguments[3]
out_file = arguments[4]
covariates_dir =arguments[5]
pc_file =arguments[6]

#############################################################################################################################################################################################################################
# Prepare covariates files for GWAS 
#############################################################################################################################################################################################################################
covariates <- read.table(paste0(covariates_dir,"/covariates_intersectids.txt"),header=T)
cols <- grep("numeric", names(covariates), value = TRUE)
sel <- c("IID",cols)
sel
covariates_numeric <- covariates[,sel]

pc <- read.table(pc_file)
colnames(pc) <- c("FID","IID","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")
numcov <- merge(pc,covariates_numeric,by.x="IID",by.y="IID",all.x=T)
numcov <- numcov[,c("FID","IID","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",cols)]
write.table(numcov, paste0(covariates_dir,"covariates_intersectids.numeric"),col.names=F, row.names=F,quote=F)

colsf <- grep("factor", names(covariates), value = TRUE)
sel2 <- c("IID",colsf)
sel2
covariates_fac <- covariates[,sel2]
faccov <- merge(pc,covariates_fac,by.x="IID",by.y="IID",all.x=T)
faccov <- faccov[,c("FID","IID",colsf)]

write.table(faccov, paste0(covariates_dir,"covariates_intersectids.factor"),col.names=F, row.names=F,quote=F)


#############################################################################################################################################################################################################################
# LOAD DATA
#############################################################################################################################################################################################################################

##### 1 fam file  #####
fam <- read.table(fam_file, stringsAsFactors=FALSE)[,1:2]  
colnames(fam) <- c("FID", "IID")

##### 2 Pheno file ##### 
# Twin Cohorts should supply a pheno_file that contains information on zygosity of the twins.
# The script below creates a phenotype file for cohorts that do not include twins, and simply assigns a "non-twin" label to all samples. 
if(pheno_file=="NULL") {
message("Creating phenotype file (note: only twin cohorts should supply a phenotype file with zygosity information)")
pheno=data.frame(fam[,"IID"],rep("non-twin",nrow(fam)))
colnames(pheno) <- c("IID", "Twinzygosity")
} else {pheno <- read.table(pheno_file, header=T, stringsAsFactors=FALSE)}


##### 3 Methylation beta values #####
message("Loading methylation data")
# Load your DNA methylation data object "beta" (rows=samples, columns=CpGs). Values=methylation beta-values.
load(methylation)
print(dim(norm.beta))
message("Checking if beta value object contains any missing values")
print(any(is.na(norm.beta)))
m <- match(fam[,"IID"], colnames(norm.beta))
beta <- norm.beta[,m]
message("Checking if IIDs match, fam file vs methylation beta-values")
print(table(fam[,"IID"]==colnames(beta)))
m2 <- match(fam[,"IID"],pheno[,"IID"])
pheno <- pheno[m2,]
message("Checking if IIDs match, phenotype file vs methylation beta-values")
print(table(pheno[,"IID"]==colnames(beta)))
rm(norm.beta)
gc()


#############################################################################################################################################################################################################################
# Calculate MZ-EpiScore
#############################################################################################################################################################################################################################

message("Predicting MZ-EpiScore")


# Load elastic net prediction model
# This model is based on 352 methylation sites that are present on the Illumina 450k and EPIC array, and was trained to distinguish MZ twins from DZ twins and non-twins
#load("/data/jvandongen/2024_GoDMC/godmc_phase2/resources/MZtwin/MZEpiScore.RData") #CHECK
load("resources/MZtwin/MZEpiScore.RData")
CpGlist <- rownames(as.matrix(coef(cv.glmmod,s="lambda.min")))[-1]
message("N CpGs elasticnet:")
length(CpGlist)
message("Number of CpGs present in this dataset:")
length(intersect(CpGlist,rownames(beta))) # 756
missingCpGs <- CpGlist[which(!CpGlist %in% rownames(beta))]
# Add missing CpGs to the input dataset
# For these CpGs, all samples will receive a value of 0 (equivalent to the mean standardized beta-value - these CpGs will not contribute to prediction).
zeros<- matrix(NA,ncol=ncol(beta),nrow=length(missingCpGs))
rownames(zeros) <- missingCpGs    # CpGnames
colnames(zeros) <- colnames(beta) # sample names
beta_imp <- rbind(beta,zeros)
beta_imp <- as.matrix(beta_imp[CpGlist,])
beta_imp <- t(beta_imp)
message("Checking if cpgids match")
table(colnames(beta_imp)==CpGlist)
rm(beta)
gc()



coeffs <- data.frame(CpGlist,coef(cv.glmmod,s="lambda.min")[-1,1]) # [-1] # - intercept
non0CpGs <- CpGlist[which(coeffs[,2]!=0)]
message("Number of CpGs utilized by MZ-EpiScore predictor")
print(length(non0CpGs)) 
message ("Number of CpGs used by MZ-EpiScore predictor missing in this dataset:")
print(length(intersect(non0CpGs,missingCpGs))) 
percentagemissingCpGs <- 100*(length(intersect(non0CpGs,missingCpGs)) /length(non0CpGs))
message("Percentage of missing CpGs used by MZ-Epi predictor:")
print(percentagemissingCpGs) 


# Standardize DNA methylation beta-values. Apply to columns (CpGs) 
IIDs <- rownames(beta_imp)
beta_imp <- apply(beta_imp,2,scale)
beta_imp[,missingCpGs] <- 0
rownames(beta_imp) <- IIDs

# Classification: Predicted MZ twin status
message("Running Classification")
predicted <- as.matrix(predict(cv.glmmod, newx =beta_imp, s = "lambda.min", type = "class"))
predicted[which(predicted==1)] <- 'Predicted MZ'
predicted[which(predicted==0)] <- 'Predicted non-MZ'
message("Frequency of predicted")
print(table(predicted))

#AUC
message("Frequency of observed")
obs_zyg <- rep(0,nrow(pheno))
obs_zyg[which(pheno$Twinzygosity=="MZ")] <- 1
obs_zyg <- obs_zyg[which(!pheno$Twinzygosity=="UZ")]
beta_imp_tmp <- beta_imp[which(!pheno$Twinzygosity=="UZ"),]
print(table(obs_zyg))

message("computing AUC")
   
if (length(unique(obs_zyg)) < 2) {
    auc <- NA
} else {
    prob <-  predict(cv.glmmod,type="response", newx =beta_imp_tmp, s = "lambda.min")
    pred <- prediction(prob,obs_zyg)
    auc <- performance(pred, "auc")@y.values[[1]]
}
message("auc (only computed in twin cohorts)")
print(auc)

# Obtain continuous MZ-EpiScores
continousscore <- as.matrix(predict(cv.glmmod, newx =beta_imp, s = "lambda.min", type = "link"))

#collect in one object
EpiMZ <- data.frame(rownames(predicted),predicted, continousscore)
colnames(EpiMZ) <- c("IID","EpiMZClassifier","MZEpiscore")
EpiMZ <- merge(EpiMZ, fam, by.x="IID", by.y="IID")
m <- match(fam[,"IID"], EpiMZ[,"IID"])
EpiMZ <- EpiMZ[m,]
message("Checking if IIDs match")
print(table(fam[,"IID"]==EpiMZ[,"IID"]))


m <- match(fam[,"IID"], pheno[,"IID"])
pheno <- pheno[m,]
message("Checking if IIDs match")
print(table(pheno[,"IID"]==EpiMZ[,"IID"]))


observed_frequency            <- table(pheno[,"Twinzygosity"])
predicted_frequency           <- table(EpiMZ[,"EpiMZClassifier"])
observedvspredicted_frequency <- table(EpiMZ[,"EpiMZClassifier"],pheno[,"Twinzygosity"])


#############################################################################################################################################################################################################################
# SAVE OUTPUT
#############################################################################################################################################################################################################################
message("observed frequency")
print(observed_frequency) 
message("predicted frequency")
print(predicted_frequency)
message("observed versus predicted frequency")
print(observedvspredicted_frequency)

# 1 Distribution figure
MZEpiscore <- EpiMZ[,"MZEpiscore"]
Zygosity   <- pheno[,"Twinzygosity"]
pdf(paste0(out_file,"MZEpiscore_distribution.pdf"))
boxplot(MZEpiscore~Zygosity,xlab="Zygosity",ylab="MZEpiscore", main="MZ-Episcore distribution")
stripchart(MZEpiscore~Zygosity, vertical = TRUE,  method = "jitter", add = TRUE, pch = 20, col = "lightblue", cex=0.5)
dev.off()

# 2 Cohort descriptives
save(observed_frequency, predicted_frequency, observedvspredicted_frequency,missingCpGs,percentagemissingCpGs,auc, file=paste0(out_file,"MZEpiscore_Frequencies.RData"))

#3 GWAS pheno files
write.table(EpiMZ[,c("FID","IID","MZEpiscore")], file=paste0(out_file, "MZEpi_all.pheno"), row=F, col=T, qu=F)

if("MZ" %in% names(observed_frequency) & observed_frequency["MZ"] > 100)
{write.table(EpiMZ[which(pheno$Twinzygosity=="MZ"),c("FID","IID","MZEpiscore")], file=paste0(out_file, "MZEpi_MZtwins.pheno"), row=F, col=T, qu=F) 
message("MZEpi_MZtwins.pheno has been created") 
} else { message("There are not enough MZ twins --> GWAS will not be run separately for MZ twin pairs") }

if ("non-twin" %in% names(observed_frequency) & observed_frequency["non-twin"] > 100) 
{write.table(EpiMZ[which(pheno$Twinzygosity=="non-twin"),c("FID","IID","MZEpiscore")], file=paste0(out_file, "MZEpi_nontwins.pheno"), row=F, col=T, qu=F) 
message("MZEpi_nontwins.pheno has been created") 
} else { message("There are not enough non-twins --> GWAS will not be run separately for non-twins") }


rm(list = ls(all = TRUE))

