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

##### command arguments  #####
arguments <- commandArgs(T)
methylation=arguments[1]
fam_file <- arguments[2]
pheno_file = arguments[3]
out_file = arguments[4]



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
message("Creating phenotype file (note: twin cohorts only should supply a phenotype file with zygosity information)")
pheno=data.frame(fam[,"IID"],rep("non-twin",nrow(fam)))
colnames(pheno) <- c("IID", "Twinzygosity")
} else {pheno <- read.table(pheno_file, header=T, stringsAsFactors=FALSE)}


##### 3 Methylation beta values #####
message("Loading methylation data")
# Load your DNA methylation data object "beta" (rows=samples, columns=CpGs). Values=methylation beta-values.
load(methylation)
print(dim(norm.beta))
m <- match(fam[,"IID"], colnames(norm.beta))
beta <- norm.beta[,m]
print(table(fam[,"IID"]==colnames(beta)))
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
message("N CpGs elasticnet")
length(CpGlist)
message("Number of CpGs present in this dataset")
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
table(colnames(beta_imp)==CpGlist)
rm(beta)
gc()



coeffs <- data.frame(CpGlist,coef(cv.glmmod,s="lambda.min")[-1,1]) # [-1] # - intercept
non0CpGs <- CpGlist[which(coeffs[,2]!=0)]
message("Number of CpGs utilized by MZ-EpiScore predictor")
print(length(non0CpGs)) 
message ("Number of CpGs used by MZ-EpiScore predictor missing in this dataset")
print(length(intersect(non0CpGs,missingCpGs))) 
percentagemissingCpGs <- 100*(length(intersect(non0CpGs,missingCpGs)) /length(non0CpGs))
message("Percentage of missing CpGs used by MZ-Epi predictor")
print(percentagemissingCpGs) 


# Standardize DNA methylation beta-values. Apply to columns (CpGs) 
IIDs <- rownames(beta_imp)
beta_imp <- apply(beta_imp,2,scale)
beta_imp[,missingCpGs] <- 0
rownames(beta_imp) <- IIDs

# Classification: Predicted MZ twin status
predicted <- as.matrix(predict(cv.glmmod, newx =beta_imp, s = "lambda.min", type = "class"))
predicted[which(predicted==1)] <- 'Predicted MZ'
predicted[which(predicted==0)] <- 'Predicted non-MZ'

# Obtain continuous MZ-EpiScores
continousscore <- as.matrix(predict(cv.glmmod, newx =beta_imp, s = "lambda.min", type = "link"))

#collect in one object
EpiMZ <- data.frame(rownames(predicted),predicted, continousscore)
colnames(EpiMZ) <- c("IID","EpiMZClassifier","MZEpiscore")
EpiMZ <- merge(EpiMZ, fam, by.x="IID", by.y="IID")
m <- match(fam[,"IID"], EpiMZ[,"IID"])
EpiMZ <- EpiMZ[m,]
print(table(fam[,"IID"]==EpiMZ[,"IID"]))


m <- match(fam[,"IID"], pheno[,"IID"])
pheno <- pheno[m,]
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
save(observed_frequency, predicted_frequency, observedvspredicted_frequency,missingCpGs,percentagemissingCpGs, file=paste0(out_file,"MZEpiscore_Frequencies.RData"))

#3 GWAS pheno file
write.table(EpiMZ[,c("FID","IID","MZEpiscore")], file=paste0(out_file, "MZEpi.pheno"), row=F, col=T, qu=F)
