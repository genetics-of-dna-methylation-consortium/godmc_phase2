
arguments <- commandArgs(T)

methy_in <- arguments[1]
methy_out <- arguments[2]

ori <- read.table(methy_in)
new <- t(ori)
df <- data.frame(id = rownames(new), new)
write.table(df, methy_out, sep = "\t", row.names = FALSE, quote = FALSE)