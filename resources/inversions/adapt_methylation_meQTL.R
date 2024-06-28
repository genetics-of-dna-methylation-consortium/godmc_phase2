
arguments <- commandArgs(T)

methy_Rdata <- arguments[1]
methy_out <- arguments[2]

load(methy_Rdata)
df <- data.frame(id = rownames(norm.beta), norm.beta)
colnames(df) <- c("id", colnames(norm.beta))
write.table(df, methy_out, sep = "\t", row.names = FALSE, quote = FALSE)