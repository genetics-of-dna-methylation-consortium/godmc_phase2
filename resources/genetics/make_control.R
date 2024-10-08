arguments <- commandArgs(T)

in_file <- arguments[1]
ids_file <- arguments[2]
out_file <- arguments[3]

a <- scan(in_file, nlines=1, what=character(),sep=",")
b <- scan(in_file, skip=1, nlines=1, what=character(),sep=",")
b <- as.numeric(b[-1])
a <- as.character(a[-1])

#a <- read.table(in_file,he=T)
ids <- read.table(ids_file, stringsAsFactors=FALSE)
#ids <- ids[match(a$IID, ids$V2), ]
ids <- ids[match(a, ids$V2), ]
#stopifnot(all(ids$V2 == a$IID))
stopifnot(all(ids$V2 == a))
#ids$phen <- a[,2]
ids$phen <- b

write.table(ids, file=out_file, row=F, col=F, qu=F)
