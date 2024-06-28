library(tidyr)
library(dplyr)
library(readr)

arguments <- commandArgs(T)

methylationfile <- arguments[1]
meth_array <- arguments[2]
out_file <- arguments[3]

load(methylationfile)
norm.beta1 <- as.data.frame(norm.beta)
rownames(norm.beta1) <- rownames(norm.beta)

annotation <- read_delim(paste0("./resources/methylation/", meth_array, "_annotation.opi")) %>%
                mutate("#chr"=CHR, start = MAPINFO, end = MAPINFO+1, gene_id = ID) %>%
                select(`#chr`, start, end, gene_id)

norm.beta1$gene_id <- rownames(norm.beta1)
meth <- merge(annotation, norm.beta1, by.x="gene_id")
meth$gene_id <- paste0("ENSG",meth$gene_id)

meth <- meth[,c(2:4,1,5:ncol(meth))]
meth <- meth[order(meth$`#chr`, meth$start),] 

for (chr in c(1:22)){
    temp <- meth %>% filter(`#chr`==chr)
    write.table(temp, paste0(out_file, "_chr",chr,".bed"), col=T, row=F, sep="\t", quote=F)
    system(paste0("gzip ",out_file, "_chr",chr,".bed"))
}
