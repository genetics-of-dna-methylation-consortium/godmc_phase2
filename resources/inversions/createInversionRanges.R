library(scoreInvHap)

data(inversionGR)

inv_df <- data.frame(inversionGR)
inv_df$seqnames <- gsub("chr", "", inv_df$seqnames)
write.table(inv_df[, c(1:3, 6)], file = "resources/inversions/inversion_ranges.txt", sep = "\t",
    col.names = FALSE, row.names = FALSE, quote = FALSE)