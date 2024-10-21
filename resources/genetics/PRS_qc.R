#load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))

#load arguments
arguments <- commandArgs(T)

PRS <- arguments[1]
PRS_file <- arguments[2]
res_dir <- arguments[3]
pheno_file <- arguments[4]
cell_counts_file <- arguments[5]
pc_file <- arguments[6]
study_name <- arguments[7]

file_PRS_hist <- paste(res_dir,"/",study_name,"_PRS_",PRS,"_hist.pdf",sep="")
file_PRS_pheno_qc_plots <- paste(res_dir,"/",study_name,"_PRS_",PRS,"_pheno_qc_plots.pdf",sep="")
file_PRS_pheno_cell_counts_plots <- paste(res_dir,"/",study_name,"_PRS_",PRS,"_cell_counts_plots.pdf",sep="")
file_PRS_pheno_pc_plots <- paste(res_dir,"/",study_name,"_PRS_",PRS,"_meth_pcs_nongenetic_untransformed_plots.pdf",sep="")

#read PRS data
df_PRS <- read.table(PRS_file,header=T,comment.char="")

#standarise and save PRS
df_PRS$SCORE <- scale(df_PRS$SCORE1_AVG)
write.table(df_PRS,file=PRS_file,quote=F,row.names=F)

#generate PRS density plot 
hist <- ggplot(df_PRS, aes(x=SCORE)) +
    geom_histogram()

pdf(file_PRS_hist)
hist
dev.off()

#generate correlation plots for cell counts
message("")
message("Correlation plots for PRS and cell counts being generated")
message("")

df_cc <- fread(paste(cell_counts_file))
df_merge_cc <- inner_join(df_PRS, df_cc, by = "IID")

list_cc <- list()
for(cc in names(df_cc)[-1]){

  x <- df_merge_cc$SCORE
  y <- df_merge_cc[[cc]]
  cor <- cor(x,y)

  p <- ggplot(df_merge_cc, aes(x=SCORE, y=.data[[cc]])) + 
    geom_point()+
    geom_smooth(method=lm) +
    ylim(min(y),max(y)+(max(y)-min(y))/10)

  p <-annotate_figure(p,
    fig.lab = paste("r =",round(cor,2)),
    fig.lab.pos = "top.right")

  list_cc[[cc]] <- p

}

#save correlation plots for cell counts
pdf(file_PRS_pheno_cell_counts_plots)
print(ggarrange(plotlist=list_cc, ncol = 3, nrow = 4))
dev.off()


#generate correlation plots for non genetic methylation PCs

df_pc <- fread(paste(pc_file,".txt",sep=""))

if(ncol(df_pc) == 1){
  message("")
  message("No correlation plots generated: no non genetic methylation PCs available")
  message("")
}

if(ncol(df_pc) > 1){
  message("")
  message("Correlation plots for PRS and non genetic methylation PCs being generated")
  message("")

  df_merge_pc <- inner_join(df_PRS, df_pc, by = "IID")

  list_pc <- list()
  for(pc in names(df_pc)[-1]){

    x <- df_merge_pc$SCORE
    y <- df_merge_pc[[pc]]
    cor <- cor(x,y)

    p <- ggplot(df_merge_pc, aes(x=SCORE, y=.data[[pc]])) +
      geom_point()+
      geom_smooth(method=lm) +
      ylim(min(y),max(y)+(max(y)-min(y))/10)

    p <-annotate_figure(p,
      fig.lab = paste("r =",round(cor,2)),
      fig.lab.pos = "top.right")

  list_pc[[pc]] <- p

  }

  #save correlation plots for PCs
  ncol <- 3
  nrow <- ceiling(length(names(df_pc)[-1])/ncol)

  pdf(file_PRS_pheno_pc_plots)
  print(ggarrange(plotlist=list_pc, ncol = ncol, nrow = nrow))
  dev.off()

}

#check whether pheno file available
if(pheno_file == "NULL"){q()}

message("")
message("Phenotypes for PRS provided: QC plots being generated")
message("")

#read pheno data
df_pheno <- read.table(pheno_file,header=T)

#format pheno data
vect_factor <- grep("factor",names(df_pheno),value=T)
vect_numeric <- grep("numeric",names(df_pheno),value=T)

for(pheno in vect_factor){df_pheno[[pheno]] <- as.factor(df_pheno[[pheno]])}
for(pheno in vect_numeric){df_pheno[[pheno]] <- as.numeric(as.character(df_pheno[[pheno]]))}

#merge
df_merge <- inner_join(df_PRS, df_pheno, by = "IID")

#generate pheno QC plots 
list_num <- list()
for(pheno in vect_numeric){

  x <- df_merge$SCORE
  y <- df_merge[[pheno]]
  cor <- cor(x,y)

 p <- ggplot(df_merge, aes(x=SCORE, y=.data[[pheno]])) + 
    geom_point() +
    geom_smooth(method=lm) +
    ylim(min(y),max(y)+(max(y)-min(y))/10)

 p <-annotate_figure(p,
    fig.lab = paste("r =",round(cor,2)),
    fig.lab.pos = "top.right")

 list_num[[pheno]] <- p

}

list_factor<- list()
for(pheno in vect_factor){

  list_factor[[pheno]] <- ggplot(df_merge, aes(x=SCORE, color=.data[[pheno]])) +
    geom_density()

}
 
#save QC plots
pdf(file_PRS_pheno_qc_plots)
print(list_num)
print(list_factor)
dev.off()



