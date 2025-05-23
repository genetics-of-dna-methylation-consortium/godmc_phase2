arguments <- commandArgs(T)
genetic_descriptives <- arguments[1]
methylation_descriptives <- arguments[2]
covariate_descriptives <- arguments[3]
phenotype_descriptives <- arguments[4]
cohort_descriptives <- arguments[5]

l <- list()
load(genetic_descriptives)
l <- c(l, cohort_summary)

load(methylation_descriptives)
l <- c(l, cohort_summary)

load(covariate_descriptives)
l <- c(l, cohort_summary)

if(file.exists(phenotype_descriptives))
{
	load(phenotype_descriptives)
	l <- c(l, cohort_summary)
}

cohort_summary <- l
save(cohort_summary, file=cohort_descriptives)
message("All descriptives are successfully saved to ./results/01/cohort_descriptives.RData")
