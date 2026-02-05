arguments <- commandArgs(T)

glint_ewas_all <- arguments[1] 
glint_ewas_unrelated <- arguments[2]
output_path <- arguments[3]
res_dir <- arguments[4] 
study_name <- arguments[5]
phen_name <- arguments[6]


library(ggplot2)

glint_ewas_all <- read.table(
  glint_ewas_all,
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

glint_ewas_all$se <- glint_ewas_all$beta/glint_ewas_all$statistic 

glint_ewas_unrelated <- read.table(
  glint_ewas_unrelated,
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

glint_ewas_unrelated$se <- glint_ewas_unrelated$beta/glint_ewas_unrelated$statistic 


joint_df <- merge(
  glint_ewas_all,
  glint_ewas_unrelated,
  by = "LMM.ID",
  suffixes = c("_all", "_unrelated")
)

# you may need to change 'coefficient' to the beta/coefficient from the original EWAS
print(paste("correlation between all and unrelated glint EWAS coefficients is:",round(cor(joint_df$beta_all, joint_df$beta_unrelated), 2)))

plot_out <- ggplot() +
  geom_point(data=joint_df, aes(x=beta_all,y=beta_unrelated), colour="#1F968BFF")+
  labs(title=paste0("Correlation all vs unrelated EWAS\ncor=",round(cor(joint_df$beta_all, joint_df$beta_unrelated), 2)))+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

jpeg(filename = paste0(output_path,"all_vs_unrelated_effect_correlation_plot.jpg"),width = 12, height = row_dimensions, units = "in", res = 600)
print(makeplots)
dev.off()


# you may need to change 'coefficient' to the beta/coefficient from the original EWAS
print(paste("correlation between all and unrelated glint standard errors is:",round(cor(joint_df$se_all, joint_df$se_unrelated), 2)))

plot_out <- ggplot() +
  geom_point(data=joint_df, aes(x=se_all,y=se_unrelated), colour="#1F968BFF")+
  labs(title=paste0("Correlation all vs unrelated EWAS\ncor=",round(cor(joint_df$se_all, joint_df$se_unrelated), 2)))+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()

jpeg(filename = paste0(output_path,"all_vs_unrelated_SE_correlation_plot.jpg"),width = 12, height = row_dimensions, units = "in", res = 600)
print(makeplots)
dev.off()
