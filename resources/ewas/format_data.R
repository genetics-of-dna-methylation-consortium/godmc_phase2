# make grm txt file for glint and convert pheno and meth to text files
# to go in godmc2 pipeline

# Args
arguments <- commandArgs(T)

grm_name <- arguments[1] # [base name of input files]
beta_file <- arguments[2] # methylation Robj - named norm.beta
phen_name <- arguments[3] # ADHD
#covs_file <- arguments[4] # methylation data already adjusted for covs
#cellcounts_cov <- arguments[4] # cell counts already adjusted for at this point in godmc
meth_pcs_file <- arguments[4] # meth PCs will be covariates in the EWAS
DEEP_scripts_directory <- arguments[5] # DEEP_mqtls github repository
# study_specific_vars <- arguments[6] # not used in this case - EWAS covariates excluding smoking, cell counts, sex (which are already in this script)
output_path <- arguments[6] 
output_extension <- arguments[7] 


library(genio)

source(paste0(scripts_directory,"/resources/datacheck/fn_rm_constant_col.R"))
source(paste0(scripts_directory,"/resources/datacheck/fn_rm_highlycor.R"))

# impute.matrix borrowed from meffil
impute.matrix <- function(x, margin=1, fun=function(x) mean(x, na.rm=T)) {
  if (margin == 2) x <- t(x)
  
  idx <- which(is.na(x) | !is.finite(x), arr.ind=T)
  if (length(idx) > 0) {
    na.idx <- unique(idx[,"row"])
    v <- apply(x[na.idx,,drop=F],1,fun) ## v = summary for each row
    v[which(is.na(v))] <- fun(v)      ## if v[i] is NA, v[i] = fun(v)
    x[idx] <- v[match(idx[,"row"],na.idx)] ##
    stopifnot(all(!is.na(x)))
  }
  
  if (margin == 2) x <- t(x)
  x
}

# read in GRM
grm <- read_grm(
  grm_name,
  n_ind = NA,
  verbose = TRUE,
  ext = "grm",
  shape = "triangle",
  size_bytes = 4,
  comment = "#"
)

grm_mat <- as.data.frame(grm$kinship)

# load in DNAm, pheno, cellcount data
# in beta_file cols should be samples and rows should be cpgs
load(beta_file)
# impute missing DNAm values with row means
# as glint can't handle missing values
norm.beta <- impute.matrix(norm.beta,1)

# organise pheno data
phen_file <- paste(home_dir,"/processed_data/genetic_data/PRS_",phen_name,".sscore",sep="")
message(paste("Loading PRS data for",phen_name))
pheno <- read.table(phen_file, header=T, stringsAsFactors=FALSE)
rownames(pheno) <- pheno$IID
pheno<-pheno[,which(names(pheno)%in%c("IID","SCORE"))]
pheno <- subset(pheno, IID %in% colnames(norm.beta), select="SCORE")
if(sum(!is.na(phen[["SCORE"]])) < 10)
{
  message("There are fewer than 10 individuals remaining. Stopping the analysis.")
  q()
}
pheno<-na.omit(pheno)

# load meth PCs
message("Loading non genetic methylation PCs")
pcs<-read.table(paste(meth_pcs_file,".txt",sep=""),sep=" ",header=T)
rownames(pcs) <- pcs$IID  

#cell_counts <- read.table(cellcounts_cov, header=T)
#rownames(cell_counts) <- cell_counts$IID

### no covs being used
### #logic variables of whether we have cov file; then load in covs
### l_cov <- ifelse(covs_file != "NULL",TRUE,FALSE)
### if(l_cov) 
### {
###   msg <- paste("Loading covariates for", phen_name)
###   message(msg)
###   
###   covs <- read.table(covs_file, he=T, stringsAsFactors=FALSE)
###   g<-grep("factor",names(covs))
###   if(length(g)>1){ 
###     for (i in 1:length(g)){
###       covs[,g[i]]<-as.factor(covs[,g[i]])
###       if(length(levels(covs[,g[i]]))==1)
###         covs<-covs[,-g[i]]
###     }
###   }
###   
###   g<-grep("numeric",names(covs))
###   if(length(g)>1){ 
###     for (i in 1:length(g)){
###       covs[,g[i]]<-as.numeric(as.character(covs[,g[i]]))
###     }
###   }
###   
###   rownames(covs) <- covs$IID
###   
### }



# match up DNAm and grm IDs and make sure they're in the same order
participants <- as.character(intersect(colnames(norm.beta),colnames(grm_mat)))
norm.beta <- norm.beta[,participants]
grm_mat <- grm_mat[,participants]
grm_mat <- grm_mat[participants,]
stopifnot(identical(colnames(norm.beta),colnames(grm_mat)))
#cell_counts <- cell_counts[participants,]
pheno <- pheno[participants,]
pcs <- pcs[participants,]
stopifnot(identical(rownames(pheno),rownames(pcs)))
stopifnot(identical(rownames(pheno),colnames(norm.beta)))


### cell counts already adjusted for
### # Cell counts file organisation
### # 1. detect cell count panel prefixes
### # This will not be needed if you have only generated cell counts with one panel
### cell_count_cols <- setdiff(colnames(cell_counts), c("FID","IID"))
### cellcount_panel_prefixes <- unique(sub("\\..*", "", cell_count_cols))
### message("Detected cell count panel prefixes: ", paste(cellcount_panel_prefixes, collapse = ", "))
### # 2. del treg cell types if present
### celltypes <- grep(paste0("^(", paste0(cellcount_panel_prefixes, collapse="|"), ")"), colnames(cell_counts), value = TRUE)
### celltypes <- celltypes[!grepl("treg", celltypes, ignore.case = TRUE)]
### cell_counts <- cell_counts[,c("IID",celltypes)]
### # 3. del nRBC if mean age > 1
### if(mean(pheno$Age_numeric) < 1){
###   message("Keeping nRBC as mean age is less than 1")
### } else {
###   message("Removing nRBC as mean age is greater than 1")
###   celltypes <- celltypes[!grepl("nRBC", celltypes, ignore.case = TRUE)]
### }
### cell_counts <- cell_counts[,c("IID",celltypes)]
### # 4. del columns with no variation
### cell_counts <- remove_constant_cols(cell_counts, "cell_counts")

### no covs
#### rm Age_numeric and Sex_factor if no variation
#### change these col names as appropriate
###if (length(unique(na.omit(pheno$Age_numeric))) < 3) {
###  message("Age_numeric has no variation (only one value). Removing Age_numeric column from pheno.")
###  pheno$Age_numeric <- NULL
###}
###if (length(unique(na.omit(pheno$Sex_factor))) < 2) {
###  message("Sex_factor has no variation (only one value). Removing Sex_factor column from pheno.")
###  pheno$Sex_factor <- NULL
###}
###
###has_study_vars <- !(length(study_specific_vars) == 1 && is.na(study_specific_vars))
###make_covs <- function(base) {
###  if (has_study_vars) {
###    unique(c(base, celltypes, study_specific_vars))
###  } else {
###    unique(c(base, celltypes))
###  }
###}

# remove if only one cell count panel removed:
#cellcount_panel <- cellcount_panel_prefixes[1]
pheno_panel <- pheno

# add cell counts to pheno_panel
##cell_counts_colnames <- grep(paste0("^", cellcount_panel, "\\."), colnames(cell_counts), value = TRUE)
##cellcounts_temp <- cell_counts[,c("IID",cell_counts_colnames)]
##cellcounts_temp <- filter_correlated_cols(cellcounts_temp, setdiff(colnames(cellcounts_temp), "IID"), thresh = 0.9, method ="pearson")
##pheno_panel <- merge(pheno_panel, cellcounts_temp, by="IID")
rownames(pheno_panel) <- pheno_panel$IID
pheno_panel <- pheno_panel[participants,]

message("Final phenotype dataframe columns: ", paste(colnames(pheno_panel), collapse = ", "))

# save out DNAm
write.table(norm.beta, file=paste0(output_path,"dnam_for_glint_",output_extension,".txt"),na = "NaN", sep = "\t", quote=FALSE, col.names = NA, row.names = TRUE)

# save out grm removing col and row names 
write.table(grm_mat, file=paste0(output_path,"grm_for_glint_",output_extension,".txt"),sep = "\t", quote = FALSE, col.names = F, row.names = F)

# write out phenotype file (age is phenotype - this may need to be changed)
phenofile <- pheno_panel[, c("IID", phen_name)]
# Write out
write.table(
  phenofile,
  file = paste0(output_path,"phenotypes_for_glint_",output_extension,".txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

pheno_panel$Sex_factor <- ifelse(pheno_panel$Sex_factor == "M", 1,
                                 ifelse(pheno_panel$Sex_factor == "F", 0, NA))

# write out covariates file 
#needs to be meth PCs
#ewas_covars_age <- make_covs(c("Sex_factor", "p_smoking_mcigarette"))
covarfile <- pcs
print(head(covarfile))
# Write out
write.table(
  covarfile,
  file = paste0(output_path,"covariates_for_glint_",output_extension,".txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
