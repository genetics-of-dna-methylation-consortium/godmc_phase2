library(tidyr)
library(readr)
library(dplyr)

args <- commandArgs(TRUE)

inputfile <- args[1]
cpglist <- args[2]
outfile <- args[3]

load(inputfile)

cpglist <- read_delim(cpglist, col_names="cpg", delim="\t")
sub <- as.data.frame(norm.beta[rownames(norm.beta) %in% cpglist$cpg,])
sub$CpG <- rownames(sub)
sub1 <- sub[,c(ncol(sub),1:(ncol(sub)-1))]

write.table(sub1, outfile, col=T, row=F, sep="\t", quote=F)
