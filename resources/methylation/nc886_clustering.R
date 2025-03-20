
arguments <- commandArgs(T)

methylation=arguments[1]
out_file = arguments[2]

if(!require(matrixStats)){
#        install.packages("matrixStats")
        library(matrixStats)
}

if(!require(janitor)){
#        install.packages("janitor")
        library(janitor)
}

if(!require(ggplot2)){
#        install.packages("ggplot2")
        library(ggplot2)
}

#The methylation data is an R matrix object where CpGs should be in rows and sample IDs should be columns.
#The rownames must be the cg identifiers and the column names the IDs that correspond to the samples, and that correspond to sample IDs in the other datasets. 
load(methylation)

methylation = norm.beta

#these are the CpG sites of the nc886 locus
cpgs=c("cg07158503", "cg11608150", "cg06478886", "cg04481923", "cg18678645", "cg06536614", "cg25340688", "cg26896946", "cg00124993", "cg08745965", "cg18797653")

#now we subset all the probes to only the 11 that are of interest to us and the IDs
probes <- methylation[rownames(methylation) %in% cpgs,]

probes <-as.data.frame(t(probes))

probes$medians=rowMedians(as.matrix(probes)) #calculate medians and add them to the dataset


########## Making a scatter plot ########## 
#Based on this we will decide where exactly the cut-off points for methylation groups

ggplot(probes, aes(x=1:nrow(probes), y=medians)) + geom_point() #visualize the plot

ggsave(paste0(out_file, "nc886_scatter.jpeg"), plot = last_plot(),
       device = NULL,
       path = NULL,
       scale = 1,
       width = NA,
       height = NA,
       units = c("in", "cm", "mm", "px"),
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)



########## Clustering ########## 
#Clustering by setting a threshold for different methylation groups

probes$methylation_group=NA

probes$methylation_group[probes$medians<0.2]=1 #here we take all the individuals with the median value below 0.2 as non-methylated
probes$methylation_group[probes$medians>=0.2]=2 #individuals between 0.2 and 0.4 are intermediates
probes$methylation_group[probes$medians>0.4]=3 #individuals between 0.4 and 0.7 are imprinted
probes$methylation_group[probes$medians>0.7]=4 #these will get discarded for all further analyses

frequency_table <- probes %>% tabyl("methylation_group", show_na=TRUE)

write.table(frequency_table, file=paste0(out_file, "nc886_frequency.txt"), row.names=F, sep="\t", quote=F) 

if(sum(probes$methylation_group==2)<10) stop('There are not enough intermediatelly methylated individuals to continue the analysis')

########## GWAS variables ##########  
#Creating the variables for the GWAS analysis

probes$intermediate_non[probes$methylation_group==1]=1
probes$intermediate_non[probes$methylation_group==2]=2
probes$intermediate_non[probes$methylation_group>=3]=NA

probes$intermediate_imprinted[probes$methylation_group==1]=NA
probes$intermediate_imprinted[probes$methylation_group==2]=2
probes$intermediate_imprinted[probes$methylation_group==3]=1
probes$intermediate_imprinted[probes$methylation_group==4]=NA

probes$non_imprinted[probes$methylation_group==1]=2
probes$non_imprinted[probes$methylation_group==2]=NA
probes$non_imprinted[probes$methylation_group==3]=1
probes$non_imprinted[probes$methylation_group==4]=NA

nc886_groups=probes[,13:16]
nc886_groups=cbind(IID=rownames(nc886_groups), nc886_groups) #these IDs will need to match the IDs in genetic data. Moreover, the column name needs to be the same, too

write.table(nc886_groups, file="${home_directory}/processed_data/methylation_data/nc886_groups.txt", row.names=F, sep="\t", quote=F)
