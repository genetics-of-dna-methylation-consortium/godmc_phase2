library(readr)
library(tidyr)
library(dplyr)

args <- commandArgs(T)

chr <- args[1]
vmeQTL_list2 <- args[2]
meth_input <- args[3]
out_file <- args[4]

list2 <- read_delim(vmeQTL_list2, col_names="cpg", delim="\t")
DNAm <- read_delim(meth_input)
DNAm_sub <- DNAm %>% filter(CpG %in% list2$cpg)
write.table(DNAm_sub, out_file, col=T, row=F, sep="\t", quote=F)
