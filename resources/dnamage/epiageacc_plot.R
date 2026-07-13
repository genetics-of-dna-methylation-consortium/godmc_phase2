###############################discription###################################
# This R script generates a comprehensive set of visualizations for epi age acceleration. 
# The script takes two command-line arguments: an input .RData file and an output filename prefix.
# The outputs include:
# 1.Density plots for overall distributions and sex-stratified comparisons
# 2.Box and violin plots to visualize group-level differences
# 3.Scatterplot matrix showing correlations and standard errors across age metrics

#################################packages####################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))

#################################Sub-functions################################
#### functions for density plot
density.plot.by.sex <- function(CorTable, ClockNames) {
  par(mfrow=c(1,3))
  for (ClockName in ClockNames) {
      PredAgeDansity <- density(CorTable[[ClockName]])
      density_F <- density(CorTable[CorTable$Sex_factor == "F", ClockName])
      density_M <- density(CorTable[CorTable$Sex_factor == "M", ClockName])
      
      plot(PredAgeDansity, xlab = ClockName, col = "white", cex.main=1, cex=0.7,
        xlim = c(min(PredAgeDansity$x) - min(PredAgeDansity$x)/4, max(PredAgeDansity$x) + max(PredAgeDansity$x)/4), 
        ylim = c(0, max(PredAgeDansity$y) + max(PredAgeDansity$y)/4), 
        main = paste("Density plot of", ClockName , "based on sex"))
      polygon(PredAgeDansity, col = alpha("#ffba49", 0.6))
      polygon(density_F, col = alpha("#ef5b5b", 0.6))
      polygon(density_M, col = alpha("#1982c4", 0.6))
      
      legendname = c(paste0(ClockName,"_All"),paste0(ClockName,"_F"), paste0(ClockName,"_M"))
      densitycolor = c("#ffba49", "#ef5b5b", "#1982c4")
      legend("topright", legend = legendname, pch = 19, col = densitycolor, inset = 0.01)
      }
}

#### Function: Combined box and violin plot
box.violin <- function(PlotData, ColorList) {
  p=ggplot(PlotData, aes(x = Group, y = Value, fill = Group)) +
        geom_violin(width=0.6) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + 
        labs(x = "Group", y = "Density") +
        theme_minimal() +
        scale_fill_manual(values=c(ColorList)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
  return(p)
}

#### Function: Box and violin plot grouped by sex
box.violin.by.sex <- function(CorTable, ClockName) {
  clock_sym <- sym(ClockName)
  plot_data <- bind_rows(
        CorTable %>% 
        mutate(Group = paste0(ClockName,"_All"),Value = !!clock_sym) %>% 
        select(Value, Group),
        CorTable %>% 
        filter(Sex_factor == "F") %>% 
        mutate(Group = paste0(ClockName,"_F"),Value = !!clock_sym) %>% 
        select(Value, Group),
        CorTable %>% 
        filter(Sex_factor == "M") %>% 
        mutate(Group = paste0(ClockName,"_M"),Value = !!clock_sym) %>% 
        select(Value, Group)
  )
  p=ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
        geom_violin(width=0.6) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + 
        labs(x = "Group", y = "Density") +
        theme_minimal() +
        scale_fill_manual(values=c("#ffba49", "#ef5b5b", "#1982c4")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
  return(p)
}


#### Panel Functions for Scatter Plot Matrix
# lower panel - correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=3)
  txt <- paste0("cor: ", r)
  text(0.5, 0.5, txt, cex=1.5)
}

# upper panel - scatter plots
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = "grey50")
}
  
# diagnal panel - sd calculation
std <- function(x) round(sd(x)/sqrt(length(x)), digits=4)
panel.se <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  se = round(std(x), digits = 3)
  txt <- paste0("se:", se)
  mtext(txt, side = 1, line = -1.5)
}

#################################Main Function####################################
main <- function()
{
  # Args
  arguments <- commandArgs(T)
  age_RData <- arguments[1]
  age_plot <- arguments[2]

  #load("/lustre/home/sww208/GoDMC/DataSetGoDMC/scz_ab_eur/processed_data/methylation_data/age_prediction.RData")
  load(age_RData)

  message("Generating density plots for age accelerations==========================")
  pdf(file = paste0(age_plot, "_density.pdf"), width=12, height=6)
  par(mfrow=c(1,2))
  if (age_valid){
    message("Plotting density plots on chronological age, DNAmAge and PhenoAge")
    # Chronological age
    legendname = c("Chronological Age")
    densitycolor = c("#ff595e")
    agedensity = density(cortable$Age_numeric)
    plot(agedensity, xlab = "Age", col = "white", cex.main=1, cex=0.7,
        xlim = c(0, max(agedensity$x) + 10), 
        ylim = c(0, max(agedensity$y) + max(agedensity$y)/4), 
        main = "Density plot of chronological age and predicted age")
    polygon(agedensity, col = alpha("#ff595e", 0.6))
    abline(v =mean(cortable$Age_numeric), lty=2, col="#ff595e")
  }else{
    if (dna_valid) {
      dnadensity = density(cortable$DNAmAge)
      plot(dnadensity, xlab = "Predicted Age", col = "white", cex.main=1, cex=0.7,
      xlim = c(0, max(dnadensity$x) + 10), 
      ylim = c(0, max(dnadensity$y) + max(dnadensity$y)/4), 
      main = "Density plot of predicted age")
      polygon(dnadensity, col = alpha("#ffca3a", 0.6))
      abline(v=mean(cortable$DNAmAge), lty=2, col= "#ffca3a")
      legendname = c("DNAmAge")
      densitycolor = c("#ffca3a")
    }else{
      message("Warning: DNAmAge does not exist, please check.")
    }
  }
    
  # PhenoAge
  if (phen_valid){
    phendensity = density(cortable$PhenoAge)
    polygon(phendensity, col = alpha( "#8ac926", 0.6))
    abline(v=mean(cortable$PhenoAge), lty=2, col= "#8ac926")
    legendname = c(legendname, "PhenoAge")
    densitycolor = c(densitycolor, "#8ac926")
  }
  legend("topleft", legend = legendname, pch = 19, col = densitycolor, inset = 0.01)

  # DunedinPACE
  if (pace_valid){
    message("Plotting density plots on DunedinPCAE")
    pacedensity = density(cortable$DunedinPACE)
    plot(pacedensity, xlab = "Pace of Aging", col = "white", cex.main=1, cex=0.7,
        main = "Density plot of DunedinPCAE")
    polygon(pacedensity, col = alpha("#1982c4", 0.6))
    abline(v=mean(cortable$DunedinPACE), lty=2, col="#1982c4")
  }
  
  # Density by sex
  if (sex_valid) {
    par(mfrow=c(1,3))
    if (dna_valid){density.plot.by.sex(cortable, c("DNAmAge", "DNAmAgeSD", "DNAmAgessSD"))}
    if (phen_valid){density.plot.by.sex(cortable, c("PhenoAge", "PhenoAgeSD", "PhenoAgessSD"))}
    if (pace_valid){density.plot.by.sex(cortable, c("DunedinPACE", "DunedinPACESD", "DunedinPACEssSD"))}
  }else{
    message("Warning: Due to the sex variable not being valid, the density plot based sex will not be plotted.")
  }
  
  
  dev.off()
  message(paste0("Density plots saved to ", age_plot, "_density.pdf"))

  #### Box and violin plots
  message("Generating box and violin plots ==========================")
  if (age_valid) {
    if (dna_valid && phen_valid){
      plot_data <- bind_rows(
        cortable %>% mutate(Group = "Age_numeric_All") %>% select(Value = Age_numeric, Group),
        cortable %>% mutate(Group = "DNAmAge_All") %>% select(Value = DNAmAge, Group),
        cortable %>% mutate(Group = "PhenoAge_All") %>% select(Value = PhenoAge, Group))
      p1 = box.violin(plot_data, ColorList=c("#ff595e", "#ffca3a", "#8ac926"))
    }else {
      message("Warning: DNAmAge and PhenoAge do not exist, please check.")
      p1 = NULL
    }
  }else{
    if (dna_valid && phen_valid){
      plot_data <- bind_rows(
        cortable %>% mutate(Group = "DNAmAge_All") %>% select(Value = DNAmAge, Group),
        cortable %>% mutate(Group = "PhenoAge_All") %>% select(Value = PhenoAge, Group))
      p1 = box.violin(plot_data, ColorList=c("#ffca3a", "#8ac926"))
    } else {
      message("Warning: DNAmAge and PhenoAge do not exist, and age variable might has 1 level, please check.")
      p1 = NULL
    }
  }

  plot_data <- bind_rows(cortable %>% mutate(Group = "DunedinPACE_All") %>% select(Value = DunedinPACE, Group))
  p2 = box.violin(plot_data, ColorList=c("#1982c4"))

  arranged_plots1 <- ggarrange(p1, p2, ncol = 2, nrow = 1)

  if (sex_valid){
    fmplot = list()
    if (dna_valid){
      dnafm = box.violin.by.sex(cortable,ClockName="DNAmAge")
      dnasdfm = box.violin.by.sex(cortable,ClockName="DNAmAgeSD")
      dnasdssfm = box.violin.by.sex(cortable,ClockName="DNAmAgessSD")
      fmplot[[length(fmplot)+1]]= list(dnafm, dnasdfm, dnasdssfm)
    }

    if (phen_valid){
      phenfm = box.violin.by.sex(cortable,ClockName="PhenoAge")
      phensdfm = box.violin.by.sex(cortable,ClockName="PhenoAgeSD")
      phensdssfm = box.violin.by.sex(cortable,ClockName="PhenoAgessSD")
      fmplot[[length(fmplot)+1]]= list(phenfm, phensdfm, phensdssfm)
    }

    if (pace_valid){
      pacefm = box.violin.by.sex(cortable,ClockName="DunedinPACE")
      pacesdfm = box.violin.by.sex(cortable,ClockName="DunedinPACESD")
      pacesdssfm = box.violin.by.sex(cortable,ClockName="DunedinPACEssSD")
      fmplot[[length(fmplot)+1]]= list(pacefm, pacesdfm, pacesdssfm)
    }

    arranged_plots=list(arranged_plots1)
    for (i in fmplot){
      arranged_plots[[length(arranged_plots)+1]] <- ggarrange(plotlist = i, ncol = 3, nrow = 1)
    }
    
  }else{message("Warning: Due to the sex variable, the box plot based on the sex will not be plotted.")
    arranged_plots=list(arranged_plots1)
  }
  ggexport(arranged_plots, filename = paste0(age_plot, "_box_violin.pdf"), width = 12, height = 6)

  #### scatter plot of correlation matrix
  message("Generating scatter plot matrix ==========================")
  if (age_valid) {matrixcol = c("Age_numeric")}else{matrixcol = c()}
  if (dna_valid) {matrixcol = c(matrixcol, "DNAmAge", "DNAmAgeSD", "DNAmAgessSD")}
  if (phen_valid) {matrixcol = c(matrixcol, "PhenoAge", "PhenoAgeSD", "PhenoAgessSD")}
  if (pace_valid) {matrixcol = c(matrixcol, "DunedinPACE", "DunedinPACESD", "DunedinPACEssSD")}

  if (length(matrixcol) > 3) {
    png(file = paste0(age_plot, "_correlation.png"), width=1400, height=800)
    matrixtable = subset(cortable, select = matrixcol)
    pairs(matrixtable, 
          lower.panel = panel.cor,
          upper.panel = upper.panel,
          diag.panel = panel.se,
          cex.labels= 1.2, gap = 0.3,
          main = "Scatterplot Matrix of Age Accelerations")
    dev.off()
    message(paste0("Scatterplot plots saved to ", age_plot, "_correlation.png"))
  } else {
    message("Warning: DNAmAge and PhenoAge do not exist, please check, fail on scatter plot.")
  }

}

main()
