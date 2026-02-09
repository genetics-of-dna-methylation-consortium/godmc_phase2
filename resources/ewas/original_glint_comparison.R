arguments <- commandArgs(T)

glint_ewas <- arguments[1] 
#original_ewas <- arguments[2]
output_path <- arguments[2]
res_dir <- arguments[3] 
study_name <- arguments[4]
phen_name <- arguments[5]
output_extension <- arguments[6] 


library(ggplot2)

glint_ewas <- read.table(
  glint_ewas,
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)


load(paste0(res_dir,"/",study_name,"_PRS_",phen_name,"_EWAS_results.RData")) # ewas.ret

# you may need to change these two lines depending on the format of the original EWAS
original_ewas <- ewas.ret$all$table
# this is creating a column of cpg IDs
original_ewas$name <- as.character(rownames(original_ewas))

joint_df <- merge(
  original_ewas,
  glint_ewas,
  by.x="name",by.y = "LMM.ID",
  suffixes = c("_original", "_glint")
)

# you may need to change 'coefficient' to the beta/coefficient from the original EWAS
print(paste("correlation between original and glint EWAS coefficients is:",round(cor(joint_df$coefficient, joint_df$beta), 2)))

plot_out <- ggplot() +
  geom_point(data=joint_df, aes(x=coefficient,y=beta), colour="#1F968BFF")+
  labs(title=paste0("Correlation original vs glint EWAS\ncor=",round(cor(joint_df$coefficient, joint_df$beta), 2)))+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

jpeg(filename = paste0(output_path,"original_vs_glint_correlation_plot_",output_extension,".jpg"),width = 12, height = row_dimensions, units = "in", res = 600)
print(makeplots)
dev.off()

