###############################description###################################

#This script is dedicated to visualizing GWAS results, specifically for QQ plots, Manhattan plots, 
#and tables with significant SNPs for each phenotype. It comprises three sections:

#1.Package Listing
#2.Functions: qqplot for QQ plot generation to identify significant SNPs.
#3.Main Function: The main function that reads GWAS data and generates the specified outputs.

# If any errors occur while running this script, 
# please don't hesitate to reach out to the developer for assistance at s.w.wang@exeter.ac.uk.

#################################packages####################################
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggtext))
suppressPackageStartupMessages(library(qqman))


#################################Sub-functions################################
qqplot = function(data, filename, lambda) {
  png(file=paste0(filename, '_qqplot.png'))
  exp_pvalues = (rank(data, ties.method="first")+.5)/(length(data)+1)
  plot(-log10(exp_pvalues), -log10(data), 
       xlim = c(0, round(max(-log10(exp_pvalues)))),
       ylim = c(0, round(max(-log10(data), -log10(exp_pvalues)))),
       main = paste0("lambda:", round(lambda, 3)))
  abline(0,1, col = 'red')
  dev.off()
}


#################################Main#########################################

main = function(){
  arguments <- commandArgs(T)

  filenames <- arguments[1]
  pval_column <- as.numeric(arguments[2])
  chr_column <- as.numeric(arguments[3])
  pos_column <- as.numeric(arguments[4])
  snp_column <- as.numeric(arguments[5])
  header <- as.logical(arguments[6])
  control_chr <- as.numeric(arguments[7])
  control_pos <- as.numeric(arguments[8])
  control_window <- as.numeric(arguments[9])
  control_threshold <- as.numeric(arguments[10])

  filenames = read.table(filenames, header = F, sep = "\t")[,1]
  for (filename in filenames) {
    message("Reading in ", filename ," GWAS results")
    GWAS_result = fread(filename, header = T, data.table=F)
    outname = unlist(strsplit(filename, "[.]"))[[1]]
    
    if(length(unique(GWAS_result[,chr_column])) > 30){
      stop("Wrong chromosome column specified")
    }

    if(min(GWAS_result[,pval_column], na.rm=TRUE) < 0 | max(GWAS_result[,pval_column], na.rm=TRUE) > 1){
      stop("Wrong column specified for p-values")
    }

    if(any(GWAS_result[, pos_column] < 0)){
      stop("Negative values in position column")
    }

    if(any(GWAS_result[, pval_column] == 0, na.rm=TRUE)){
      w <- which(GWAS_result[,pval_column]==0)
      GWAS_result[w,pval_column] <- as.numeric(.Machine$double.xmin)
    }
    w <- which(GWAS_result[,chr_column]!=control_chr)
    a_minuschr <- GWAS_result[w,]

    
    if(control_chr != 0){
      index <- GWAS_result[,pos_column] > (control_pos - control_window) & GWAS_result[,pos_column] < (control_pos + control_window)
      
      GWAS_result_filter <- GWAS_result[index, ]
      min_pval <- min(GWAS_result_filter[,pval_column], na.rm=TRUE)
      
      message("\n\nExpecting a large meQTL near ", control_chr, ":", control_pos)
      message("Lowest p-value within ", control_window, " base pairs: ", min_pval)
      
	  if(min_pval > control_threshold) {
        message("WARNING!")
        message("There doesn't appear to be a QTL for this positive control")
        message("Please upload this section and contact GoDMC analysts before continuing.\n\n")
      	}
      
    chisq = qchisq(a_minuschr[,pval_column],1,lower.tail=FALSE)
    lambda = median(chisq) / qchisq(0.5,1)
    qqplot(data=a_minuschr[,pval_column], filename=paste0(outname, "_nocisChr"),lambda=lambda)
    message("Generating QQ-plot without cis chromosome for", outname, " with lambda ", lambda)

    message(paste0('Generating manhantten plot without cis chromosome ', outname))
    man_data = a_minuschr[order(a_minuschr[,pos_column], decreasing = F),]
    man_data = subset(man_data, -log10(man_data[,pval_column]) > 2)
      
    pdf(file=paste0(outname, '_nocisChr_manhattan.pdf'), width=50, height=10)
    manhattan(man_data, bp=names(man_data)[pos_column], 
            chr=names(man_data)[chr_column], 
            snp=names(man_data)[snp_column],
            ylim=c(2,max(-log10(man_data[,pval_column])+1)))
    dev.off()
      
    message("The following plots have been generated, please check!\n",
          paste0(outname , "_nocisChr_qqplot.png\n"),
          paste0(outname ,"_nocisChr_manhattan.png"))
    }
    
    
    chisq = qchisq(GWAS_result[,pval_column],1,lower.tail=FALSE)
    lambda = median(chisq) / qchisq(0.5,1)
    qqplot(data=GWAS_result[,pval_column], filename=outname, lambda=lambda)
    message("Generating QQ-plot for", outname, " with lambda ", lambda)

    
    message(paste0('Generating manhattan plot for ', outname))
    #man_data = GWAS_result[order(GWAS_result[,pos_column], decreasing = F),]
    man_data = subset(GWAS_result, -log10(GWAS_result[,pval_column]) > 2)
    
    pdf(file=paste0(outname, '_manhattan.pdf'), width=50, height=10)
    manhattan(man_data, bp=names(man_data)[pos_column], 
              chr=names(man_data)[chr_column], 
              snp=names(man_data)[snp_column],
              ylim=c(2,max(-log10(man_data[,pval_column])+1)))
    dev.off()

    message("The following plots have been generated, please check!\n",
            paste0(outname , "_qqplot.png\n"),
            paste0(outname ,"_manhattan.pdf"))
  }


}
  

main()
