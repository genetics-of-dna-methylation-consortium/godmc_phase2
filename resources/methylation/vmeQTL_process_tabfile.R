library(tidyr)
library(dplyr)
library(readr)

args <- commandArgs(TRUE)
inputfile <- args[1]
outfile <- args[2]

load(inputfile)
annotation <- read_delim("./resources/methylation/vmeQTL/EPIC_annotation_autochr.opi", col_names=c("chr","CpG","pos","gene","strand"))

runtab <- function(i){
  tabfile <- annotation %>% filter(chr==i)
  temp <- norm.beta[rownames(norm.beta) %in% tabfile$CpG,]
  temp1 <- temp %>% as_tibble() %>% mutate(CpG=rownames(temp))
  temp1 <- temp1[,c(ncol(temp1), 1:ncol(temp1)-1)]
  write.table(temp1, file=paste0(outfile,i), col=T, row=F, sep="\t", quote=F)
}

lapply(c(1:22),function(x) runtab(x))
