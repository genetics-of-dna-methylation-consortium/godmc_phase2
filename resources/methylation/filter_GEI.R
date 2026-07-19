library(tidyr)
library(dplyr)
library(arrow)
library(future.apply)

args <- commandArgs(T)
GEIpath <- args[1]
files <- list.files(path = GEIpath, pattern = "genetic_pc")
message(paste0(length(files), " G-genetic PC interaction files available"))

setwd(GEIpath)
df_list <- future_lapply(files, function(file) {
  dat <- read_parquet(file)
  pc <- regmatches(file, regexpr("(?<=genetic_pc)\\d+", file, perl = TRUE))
  dat$pc <- pc 
  return(dat[dat$pval_gi < 5e-8,])
})

df <- rbindlist(df_list)
write.table(df, file="GEI_geneticPC_interaction_5e-8.csv", col=T, row=F, sep="\t", quote=F)
