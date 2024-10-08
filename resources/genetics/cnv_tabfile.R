arguments <- commandArgs(T)

cnvfile <- arguments[1]
tabfile <- arguments[2]
id_file <- arguments[3]
nchunks <- as.numeric(arguments[4])

message("Loading CNV data")
load(cnvfile)

message("Extracting QC'd individuals")
ids <- scan(id_file, what="character")
cnv <- cnv[,colnames(cnv) %in% ids]
message("Keeping ", ncol(cnv), " individuals")

message("Splitting ", nrow(cnv), " CNV variables into ", nchunks, " chunks.")
chunks <- split(1:nrow(cnv), 1:nchunks)
for(i in 1:length(chunks))
{
	message(i, " of ", length(chunks))
	write.table(cnv[chunks[[i]], ], file=paste0(tabfile, ".tab.", i), row=T, col=T, qu=F, sep="\t")
}
