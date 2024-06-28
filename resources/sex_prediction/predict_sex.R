suppressMessages(library(meffil))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

arguments <- commandArgs(T)

methylation_file <- arguments[1]
meth_array <- arguments[2]
covar_file <- arguments[3]
out_file <- arguments[4]
plot_path <- arguments[5]

#The expected format of the methylation array name (as given in the config-file) is one of "450k", "epic", or "epic2". Return an error if it's anything else.
if(!(meth_array %in% c("450k", "epic", "epic2"))){
  stop("The name of the methylation array should be one set to one of the following: 450k, epic, or epic2 (case sensitive). Please check the config file!")
}

if(meth_array == "450k"){
  message("According to the config-file, the DNA-methylation data were obtained using the 450K array.")
} else if(meth_array == "epic"){
  message("According to the config-file, the DNA-methylation data were obtained using the EPIC array.")
} else if (meth_array == "epic2"){
  message("According to the config-file, the DNA-methylation data were obtained using the EPIC v2 array.")
}

#Load methylation data.
load(methylation_file)

#Load covariate file.
covar <- read.table(covar_file, header = T)

#Sex prediction won't work if the distribution of sex is too skewed. If this is the case (less than 10% for one sex), skip sex prediction.
if((length(covar$Sex_factor[covar$Sex_factor == "M"]) < 0.05 * length(covar$Sex_factor) |
    length(covar$Sex_factor[covar$Sex_factor == "F"]) < 0.05 * length(covar$Sex_factor))){
  
  message("Sexes are extremely skewed in this cohort (<5% males or females). Skipping sex-prediction.")
  
  sex_discrepancies <- character(0)
  save(sex_discrepancies, file=paste0(out_file))
  
} else{
  
  annots <- meffil.get.features(meth_array)
  annots <- annots[!is.na(annots$chromosome),]
  
  x_probes <- annots$name[annots$chromosome == "chrX"]
  y_probes <- annots$name[annots$chromosome == "chrY"]
  
  # Select the X and Y chromosome probes.
  beta_x <- norm.beta[rownames(norm.beta) %in% x_probes,]
  beta_y <- norm.beta[rownames(norm.beta) %in% y_probes,]
  
  #Remove NAs (these will throw off the PCA).
  beta_x <- na.omit(beta_x)
  beta_y <- na.omit(beta_y)
  
  message(paste0(nrow(beta_x), " X-chromosome probes were selected."))
  message(paste0(nrow(beta_y), " Y-chromosome probes were selected."))
  
  #If the dataset contains no sex-chromosomal probes, skip sex-prediction.
  if(nrow(beta_x) == 0 | nrow(beta_y) == 0){
    message("No sex-chromosomal probes detected for the X and/or Y chromosomes. Skipping sex-prediction.")
    
    sex_discrepancies <- character(0)
    save(sex_discrepancies, file=paste0(out_file))
    
    
  } else{
    
    message("Predicting sex from DNA-methylation on sex-chromosomes.")
    
    #Perform a PCA on the X and Y chromosomes.
    pca_x <- prcomp(t(beta_x))
    pca_y <- prcomp(t(beta_y))
    
    # Make a dataframe containing the mean X and Y methylation per person.
    pred_sex_temp <- data.frame(
      IID = colnames(norm.beta),
      x_PC1 = pca_x$x[,1],
      y_PC1 = pca_y$x[,1]
      )
#      assumed_sex = factor(covar$Sex_factor)
 
    pred_sex <- merge(pred_sex_temp, covar[,c("IID","Sex_factor")], by.x="IID") %>% 
        mutate(assumed_sex = factor(Sex_factor)) %>% select(-Sex_factor)

    #The sign of principal components is arbitrary. To enable the use of PC1 as a cutoff, define it as follows:
    # - Females higher than males for PC1 of the X-chromosome. 
    # - Males higher than females for PC1 of the Y-chromosome.
    if(mean(pred_sex$x_PC1[pred_sex$assumed_sex == "F"], na.rm = T) < mean(pred_sex$x_PC1[pred_sex$assumed_sex == "M"], na.rm = T)){
      pred_sex$x_PC1 <- -(pred_sex$x_PC1)
    }
    
    if(mean(pred_sex$y_PC1[pred_sex$assumed_sex == "F"], na.rm = T) > mean(pred_sex$y_PC1[pred_sex$assumed_sex == "M"], na.rm = T)){
      pred_sex$y_PC1 <- -(pred_sex$y_PC1)
    }
    
    #The first PC of the X/Y methylation should be strongly correlated to sex. If this is not the case, skip sex prediction. 
    x_cor <- cor(as.numeric(pred_sex$assumed_sex), pred_sex$x_PC1, use = "complete.obs")
    y_cor <- cor(as.numeric(pred_sex$assumed_sex), pred_sex$y_PC1, use = "complete.obs")
    
    if(abs(x_cor) < 0.7 | abs(y_cor) < 0.7){
      warning("The first PC of X/Y methylation does not correlate sufficiently (R > 0.7) with sex. Sex may be mislabeled!")
      
    }
    
    # Based on the PC of the X and Y-chromosomes, divide the samples into males and females.
    x_cutoff <- mean(c(min(pred_sex$x_PC1), max(pred_sex$x_PC1)))
    y_cutoff <- mean(c(min(pred_sex$y_PC1), max(pred_sex$y_PC1)))
    
    pred_sex$X <- ""
    pred_sex$Y <- ""
    
    pred_sex$X[pred_sex$x_PC1 > x_cutoff] <- "XX"
    pred_sex$X[pred_sex$x_PC1 <= x_cutoff] <- "X"
    
    pred_sex$Y[pred_sex$y_PC1 <= y_cutoff] <- ""
    pred_sex$Y[pred_sex$y_PC1 > y_cutoff] <- "Y"
    
    pred_sex$sex_chromosomes <- paste0(pred_sex$X, pred_sex$Y)
    
    # Add the predicted sex according to the X and Y chromosomes.
    # If the two match, add their shared result as the predicted sex.
    # If the X-prediction says female (marking two X-chromsomes) while the Y-prediction says male (marking a Y-chromosome), the subject is XXY (Klinefelter).
    # If the X-prediction says male (no two X-chromosomes) while the Y-prediction says male (no Y-chromosome), the subject is X (Turner Syndrome).
    pred_sex$predicted_sex <- factor(NA, levels = c("F", "M", "Klinefelter", "Turner"))
    
    pred_sex$predicted_sex[pred_sex$sex_chromosomes == "XX"] <- "F" # XX = female
    pred_sex$predicted_sex[pred_sex$sex_chromosomes == "XY"] <- "M" # XY = male
    pred_sex$predicted_sex[pred_sex$sex_chromosomes == "XXY"] <- "Klinefelter" # XXY = Klinefelter
    pred_sex$predicted_sex[pred_sex$sex_chromosomes == "X"] <- "Turner" # X = Turner
    
    table(pred_sex$assumed_sex, pred_sex$predicted_sex, useNA = "always")
    
    # Plot the mean X/Y methylation for assumed males and females.
    plot_theme <- list(
      geom_text(size = 1, hjust = 0.5, vjust = 1.5),
      scale_color_manual(values = c("red", "blue", "darkgreen", "purple")),
      theme_bw(),
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank()),
      labs(x = "Subject ID", y = "PC1")
    )
    
    pdf(file = plot_path, width = 6, height = 6)
    
    #Plot the first PC of the X/Y methylation against the assumed sex supplied in the covariate file.
    print(ggplot(data = pred_sex, aes(x = IID, y = x_PC1, color = assumed_sex, label = IID)) + 
            geom_hline(yintercept = x_cutoff, linetype = "dashed", alpha = 0.5) + 
            plot_theme + 
            labs(title = "PC1 of X methylation vs assumed sex"))
    
    print(ggplot(data = pred_sex, aes(x = IID, y = y_PC1, color = assumed_sex, label = IID)) + 
            geom_hline(yintercept = y_cutoff, linetype = "dashed", alpha = 0.5) + 
            plot_theme + 
            labs(title = "PC1 of Y methylation vs assumed sex"))
    
    
    #Plot the first PC of the X/Y methylation against the methylation-predicted sex.
    print(ggplot(data = pred_sex, aes(x = IID, y = x_PC1, color = predicted_sex, label = IID)) + 
            geom_hline(yintercept = x_cutoff, linetype = "dashed", alpha = 0.5) + 
            plot_theme + 
            labs(title = "PC1 of X methylation vs predicted sex"))
    
    print(ggplot(data = pred_sex, aes(x = IID, y = y_PC1, color = predicted_sex, label = IID)) + 
            geom_hline(yintercept = y_cutoff, linetype = "dashed", alpha = 0.5) + 
            plot_theme + 
            labs(title = "PC1 of Y methylation vs predicted sex"))
    
    dev.off()
    
    # Print how many samples were sex-mismatched.
    total_samples <- pred_sex$IID
    
    matched_idx <- as.character(pred_sex$assumed_sex) == as.character(pred_sex$predicted_sex)
    matched_idx[is.na(matched_idx)] <- FALSE
    matched_samples <- total_samples[matched_idx]
    
    mismatched_samples <- total_samples[!(matched_idx)]
    
    total_length <- length(total_samples)
    matched_length <- length(matched_samples)
    mismatched_length <- length(mismatched_samples)
    
    message(sprintf("For %s out of %s samples, predicted sex matched the sex supplied in the covariate file.", matched_length, total_length))
    message(sprintf("For %s out of %s samples, a sex discrepancy was detected. These samples will be removed.", mismatched_length, total_length))
    
    # Save samples with a sex discrepancy. These samples will be removed later.
    sex_discrepancies <- mismatched_samples
    save(sex_discrepancies, file=paste0(out_file))
    
    # If the mismatched samples make up at least 10% of all samples, stop execution and return an error.
    mismatched_percentage <- mismatched_length / total_length
    if (mismatched_percentage > 0.10) {
      stop("Percentage of sex discrepancies is too high (> 10%). Please check your input data!")
    }
  }
}
