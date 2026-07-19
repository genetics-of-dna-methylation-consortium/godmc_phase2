suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
library(readr)
suppressMessages(library(meffil))

args <- commandArgs(TRUE)

inputfile <- args[1]
outfile <- args[2]
list1 <- args[3]
bychr <- args[4]

load(inputfile)

if (bychr==TRUE){
    combined_df <- read_delim(list1, delim="\t", col_names=c("SNP","CpG","pair"))
}else{
    combined_df <- read_delim(list1, delim="\t", col_names=c("CpG"))
}

annots1 <- meffil.get.features('epic') %>%
                mutate(`#chr` = chromosome, start = position, end = position+1, 
                        gene_id = name, CpG = name) %>%
                dplyr::select(`#chr`, start, end, gene_id, CpG)

DNAm_for_GE_input_by_chr <- function(i){
    annots1_tmp <- annots1 %>% filter(`#chr`==paste0("chr",i)) %>% filter(CpG %in% combined_df$CpG)
    temp <- norm.beta[rownames(norm.beta) %in% annots1_tmp$CpG,] %>% as.data.frame()
    temp$CpG <- rownames(temp)
    temp1 <- merge(annots1, temp, by.x="CpG") %>% dplyr::select(-CpG)
    temp1_sorted <- temp1[order(temp1$start),]
    message(paste0("The number of CpGs on chr", i, " for GE interaction analysis is: ", nrow(temp1_sorted)))
    write_tsv(temp1_sorted, file=paste0(outfile,i,".bed.gz"))
    return(invisible())
}

DNAm_for_GE_input_once <- function(){
    annots1_tmp <- annots1 %>% filter(CpG %in% combined_df$CpG)
    temp <- norm.beta[rownames(norm.beta) %in% annots1_tmp$CpG,] %>% as.data.frame()
    temp$CpG <- rownames(temp)
    temp1 <- merge(annots1, temp, by.x="CpG") %>% dplyr::select(-CpG)
    temp1_sorted <- temp1[order(temp1$start),]
    write_tsv(temp1_sorted, file=paste0(outfile,".bed.gz"))
    return(invisible())
}

if (bychr==TRUE){lapply(c(1:22),function(x) DNAm_for_GE_input_by_chr(x))}else{DNAm_for_GE_input_once()}
