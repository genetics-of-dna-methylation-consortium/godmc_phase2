library(readr)

args <- commandArgs(T)
oii_file <- args[1]
fam_file <- args[2]

oii <- read_delim(oii_file, col_names=c("FID","IID","col1","col2","col3"))
fam <- read_delim(fam_file, col_names=c("FID","IID","col1","col2","col3","col4"))

fam1 <- fam[match(oii$IID, fam$IID),]

if (all(fam1$IID == oii$IID) == TRUE){
    oii$FID <- fam1$FID
    write.table(oii, oii_file, col=F, row=F, sep="\t", quote=F)
}else{
    message("Unmatched FID between oii file and fam file!")}
