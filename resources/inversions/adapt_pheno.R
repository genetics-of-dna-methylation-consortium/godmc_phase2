args <- commandArgs(trailingOnly=TRUE)
phenotype_file <- args[1]
fam_file <- args[2]
pheno <- args[3]
out <- args[4]

phenodf <- read.table(phenotype_file, header = TRUE)
fam <- read.table(fam_file, header = FALSE)
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")

if (!pheno %in% colnames(phenodf)){
    stop(paste(pheno, "not present in phenotypes file. Check the variable name passed to ./08d-IWAS_phenotype.sh."))
}

combdf <- merge(phenodf, fam, by = "IID")

outdf <- combdf[, c("FID", "IID", pheno)]
write.table(outdf, col.names = FALSE, row.names = FALSE, file = paste0(out, pheno, ".plink"),
            quote = FALSE)