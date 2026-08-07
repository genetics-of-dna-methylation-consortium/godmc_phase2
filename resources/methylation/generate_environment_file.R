library(tidyr)
library(readr)
library(dplyr)

args <- commandArgs(T)
cov_file <- args[1]
BMI_file <- args[2]
out_file <- args[3]
E_plots <- args[4]
E_summary <- args[5]

cov <- read.table(cov_file, he=T, stringsAsFactors=F, colClass=c("Sex_factor"="character"))
BMI <- read.table(BMI_file, he=T, stringsAsFactors=F)

df <- merge(cov, BMI, by.x="IID", all=T)
df$Sex_numeric[df$Sex_factor=="F"] <- 1
df$Sex_numeric[df$Sex_factor=="M"] <- 0

if (length(unique(df$Age_numeric))>1 & length(unique(df$Sex_numeric))>1){
    df$Age_sex <- df$Age_numeric*df$Sex_numeric
    df$Age2 <- (as.numeric(df$Age_numeric))^2
    df$Age2_sex <- df$Age2*df$Sex_numeric
}

if (length(unique(df$Age_numeric))>1 & length(unique(df$Sex_numeric))==1){
    df$Age2 <- (as.numeric(df$Age_numeric))^2
}

inc_names <- c("IID","Baso","Bmem","Bnv","CD4Tmem","CD4Tnv","CD8Tmem","CD8Tnv","Eos","Neu","NK","Mono","Treg","Smoking","Sex_numeric","Age_numeric","Age2","Age_sex","Age2_sex","BMI", paste0("genetic_pc",c(1:20)))
removed <- colnames(df)[(colnames(df) %in% inc_names)==F]
removed_names <- paste(removed, collapse = ", ")
message(paste0(length(removed)," factors were excluded in the GxE interaction analysis. They are ", removed_names))

df_tmp <- df[,colnames(df) %in% inc_names]
df_cleaned <- df_tmp[sapply(df_tmp, function(x) length(unique(x))) > 1]
df_cleaned[,2:ncol(df_cleaned)] <- df_cleaned[,2:ncol(df_cleaned)] %>% mutate_if(is.character, as.numeric)
factor_names <- paste(colnames(df_cleaned)[2:ncol(df_cleaned)], collapse = ", ")
message(paste0(ncol(df_cleaned)-1," factors were included in the GxE interaction analysis. They are ", factor_names))

write.table(df_cleaned, out_file, col=T, row=F, sep="\t", quote=F)

message("Calculating effective number of Es")
effectE <- function(df){
    pca <- prcomp(x = df %>% as.matrix(), scale = TRUE)
    eigenvalues <- (pca$sdev)^2
    var1 <- (sum(eigenvalues))^2
    var2 <- sum(eigenvalues^2)
    effectEnvNum <- var1/var2
    return(effectEnvNum)
}

env1 <- df_cleaned[,(colnames(df_cleaned) %in% c("IID",paste0("genetic_pc", c(1:20)))) == F] %>% na.omit()
env1_names <- paste(colnames(env1), collapse = ", ")
n_env1 <- effectE(env1)
message(paste0(ncol(env1)," factors to calculate effective numbers of E. They are ", env1_names, ". \nThe effective number is ", n_env1))

env2 <- df_cleaned[,(colnames(df_cleaned) %in% c("IID", "BMI", paste0("genetic_pc", c(1:20)))) == F] %>% na.omit()
env2_names <- paste(colnames(env2), collapse = ", ")
n_env2 <- effectE(env2)
message(paste0(ncol(env2)," factors to calculate effective numbers of E. They are ", env2_names, ". \nThe effective number is ", n_env2))


pdf(E_plots, width=8, height=8)
par(mfrow = c(2,2))
for (i in 2:ncol(df_cleaned)) {
  # Get the cell type for the current iteration
  E_type <- colnames(df_cleaned)[i] 
  # Generate plots
  plot(df_cleaned[,i], main = E_type, xlab = "Sample", ylab = "Value") # Cell count per sample.
  hist(df_cleaned[,i], main = E_type, xlab = "Value", ylab = "Distribution", # Histogram of cell counts. 
          col = "lightgrey", border = "black")
}

suppressMessages(dev.off())

# Save the distribution of cell counts for this cohort.
library(matrixStats)
E_mat <- as.matrix(df_cleaned[-1])
E_sum <- data.frame(
  mean = colMeans(E_mat, na.rm = T),
  sd = colSds(E_mat, na.rm = T),
  min = colMins(E_mat, na.rm = T),
  perc_0.25 = colQuantiles(E_mat, probs = 0.25, na.rm = T),
  median = colMedians(E_mat, na.rm = T),
  perc_0.75 = colQuantiles(E_mat, probs = 0.75, na.rm = T),
  max = colMaxs(E_mat, na.rm = T),
  NAs = colSums(is.na(E_mat))
)

write.table(E_sum, file = E_summary, quote = FALSE, row.names = TRUE)


