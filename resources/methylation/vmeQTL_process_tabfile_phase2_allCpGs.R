library(tidyr)
library(dplyr)
library(readr)
library(meffil)

args <- commandArgs(TRUE)

inputfile <- args[1]
meth_chunk <- as.numeric(args[2])
outfile <- args[3]

load(inputfile)

annots1 <- meffil.get.features('epic') %>%
                mutate(`#chr` = chromosome, start = position, end = position+1, 
                        gene_id = name, CpG = name) %>%
                dplyr::select(`#chr`, start, end, gene_id, CpG) %>% na.omit() %>% filter(!(`#chr` %in% c("chrX","chrY")))
annots1 <- annots1[order(annots1$`#chr`),]
chunksize <- round(nrow(norm.beta)/meth_chunk)

DNAm_for_GE_input <- function(i){
    message(paste0("generating methylation data chunk ",i))
    i1 <- (i-1)*chunksize + 1
    i2 <- i*chunksize
    annots1_tmp <- annots1[c(i1:i2), ]
    temp <- norm.beta[rownames(norm.beta) %in% annots1_tmp$CpG,] %>% as.data.frame()
    temp$CpG <- rownames(temp)
    temp1 <- merge(annots1, temp, by.x="CpG") %>% dplyr::select(-CpG)
    temp1_sorted <- temp1[order(temp1$start),]
    write_tsv(temp1_sorted, file=paste0(outfile,i,".bed.gz"))
    return(invisible())
}

lapply(c(1:meth_chunk),function(x) DNAm_for_GE_input(x))
