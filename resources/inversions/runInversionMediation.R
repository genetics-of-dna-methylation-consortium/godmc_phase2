args <- (commandArgs(TRUE))
betas <- as.character(args[1])
inversions <- as.character(args[2])
phenofile <- as.character(args[3])
phenotype <- as.character(args[4])
covfile <- as.character(args[5])
res_folder <- args[6]
threads <- as.numeric(args[7])

out_file <- paste0(res_folder, phenotype, ".mediation_invs.txt")
library(minfi)
library(mediation)


load(betas)
invs <- read.table(inversions, header = TRUE)
inv_ranges <- read.table("resources/inversions/inversion_ranges.txt")
phenos <- read.table(phenofile, header = TRUE)
covars <- read.table(covfile, header = TRUE)

if (!phenotype %in% colnames(phenos)){
    stop(paste(phenotype, "not present in phenotypes file. Check the variable name passed to ./08e-mediation_inversion_CpG_phenotype.sh."))
}


## Convert betas to GenomicRatioSet
grset <- makeGenomicRatioSetFromMatrix(norm.beta)

## Convert ranges to GRanges
colnames(inv_ranges) <- c("chromosome", "start", "end", "inversion")
inv_GR <- makeGRangesFromDataFrame(inv_ranges)
names(inv_GR) <- inv_ranges$inversion

## Ensure ranges are expressed in the same coordinates
seqlevelsStyle(inv_GR) <- "NCBI"
cpg_ranges <- rowRanges(grset)
seqlevelsStyle(cpg_ranges) <- "NCBI"

invs_df <- t(invs[, -1])
colnames(invs_df) <- invs$id
rownames(phenos) <- phenos$IID
rownames(covars) <- covars$IID
covars <- covars[, -1]

# extra_geneticPCs <- colnames(covars)[grep("genetic", colnames(covars))]
# extra_geneticPCs <- extra_geneticPCs[!extra_geneticPCs %in% paste0("genetic_pc", 1:10)]

# covars <- covars[, !colnames(covars) %in% extra_geneticPCs]

## Select samples in common between phenotypes, covariates, inversions and methylation
com_samps <- Reduce(intersect, list(rownames(phenos), rownames(invs_df), colnames(grset), rownames(covars)))

## Select inversions present in the cohort
all_invs <- colnames(invs_df)

covars <- covars[, colnames(covars) != "Treg"]

computeMediation <- function(cpg, inv, pheno){
  df <- data.frame(cpg = as.numeric(getBeta(grset[cpg, com_samps])),
                   inv = as.numeric(invs_df[com_samps, inv]),
                   pheno = phenos[com_samps, pheno])
  df <- cbind(df, covars[com_samps, ])
  mod.med <- lm(formula(paste("cpg ~ inv + ", paste(colnames(covars), collapse = "+"))), df)
  if (grepl("numeric", pheno)){
    mod.out <- lm(formula(paste("pheno ~ cpg + inv + ", paste(colnames(covars), collapse = "+"))), df)
  } else if (grepl("factor", pheno)){
    df$pheno <- as.numeric(factor(df$pheno)) - 1
    mod.out <- glm(formula(paste("pheno ~ cpg + inv + ", paste(colnames(covars), collapse = "+"))), family = "binomial", df)
  }

  med <- mediate(mod.med, mod.out, treat = "inv", mediator = "cpg", long = FALSE, sims = 400)
  out <- data.frame(cpg = cpg, inv = inv, pheno = pheno, prop = med$n.avg, 
    ci_low = med$n.avg.ci[1], ci_high = med$n.avg.ci[2], pval = med$n.avg.p)
  out
}
#computeMediation(names(cpg_ranges)[1], inv = all_invs[1], pheno = phenotype)

inv_res <- lapply(all_invs, function(inv){

    message("Computing mediation for inversion ", inv)
    inv_grset <- subsetByOverlaps(cpg_ranges, inv_GR[inv] + 150e3) ## Include CpGs closer than 150Kb to breakpoints
    med_res <- mclapply(names(inv_grset), computeMediation, inv = inv, pheno = phenotype, mc.cores = threads)
    med_res <- Reduce(rbind, med_res)

})
out_df <- Reduce(rbind, inv_res)
write.table(out_df, file = out_file, col.names = TRUE, row.names = FALSE, quote = FALSE)




# computeMediation2 <- function(cpg, inv, pheno, sims = 1000){
#   df <- data.frame(cpg = as.numeric(getBeta(grset[cpg, com_samps])),
#                    inv = as.numeric(invs_df[com_samps, inv]),
#                    pheno = phenos[com_samps, pheno])
#   df <- cbind(df, covars[com_samps, ])
#   mod.med <- lm(formula(paste("cpg ~ inv + ", paste(colnames(covars), collapse = "+"))), df)
#   mod.out <- lm(formula(paste("pheno ~ cpg + inv + ", paste(colnames(covars), collapse = "+"))), df)
#   med <- mediate(mod.med, mod.out, treat = "inv", mediator = "cpg", long = FALSE, sims = sims)
#   med
# }

# cot <- microbenchmark(s100 =  computeMediation2("cg08670715", "inv17_007", "Height_numeric", sims = 100),
# s200 = computeMediation2("cg08670715", "inv17_007", "Height_numeric", sims = 200),
# s1000 = computeMediation2("cg08670715", "inv17_007", "Height_numeric", sims = 1000), times = 20L)

# s100_l <- mclapply(1:100, function(i) computeMediation2("cg01445052", "inv1_008", "Height_numeric", sims = 100), mc.cores = 20 )
# s200_l <- mclapply(1:100, function(i) computeMediation2("cg08670715", "inv17_007", "Height_numeric", sims = 200), mc.cores = 20 )
# s300_l <- mclapply(1:100, function(i) computeMediation2("cg08670715", "inv17_007", "Height_numeric", sims = 300), mc.cores = 20 )
# s400_l <- mclapply(1:100, function(i) computeMediation2("cg08670715", "inv17_007", "Height_numeric", sims = 400), mc.cores = 20 )
# s500_l <- mclapply(1:100, function(i) computeMediation2("cg08670715", "inv17_007", "Height_numeric", sims = 500), mc.cores = 20 )
# s1000_l <- mclapply(1:100, function(i) computeMediation2("cg08670715", "inv17_007", "Height_numeric", sims = 1000), mc.cores = 20 )

# a <- computeMediation("cg08670715", "inv17_007", "Height_numeric")

# getTab <- function(med){

#     c(med$n.avg, med$n.avg.ci)
# }

# combTab <- Reduce(rbind, lapply(list(s100_l, s200_l, s300_l, s400_l, s500_l, s1000_l), function(x) t(sapply(x, getTab)))) %>%
#     as_tibble()

# colnames(combTab) <- c("Prop", "CI_low", "CI_high")
# combTab <- combTab %>%
#     mutate(index = seq_len(nrow(combTab)),
#            Dataset = rep(c(100, 200, 300, 400, 500, 1000), each = 100),
#            Dataset = factor(Dataset, levels = c(100, 200, 300, 400, 500, 1000)))

# ggplot(filter(combTab, !Dataset %in% c(100, 200)), aes(x = Prop, xmin = CI_low, xmax = CI_high, y = index, color = Dataset)) +
#     geom_point() +
#     geom_errorbar()
# dev.off()