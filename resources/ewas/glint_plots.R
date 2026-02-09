# miami and qq plots for GLINT EWAS

# Args
arguments <- commandArgs(T)

glint_ewas <- arguments[1] 
output_path <- arguments[2]
study_name <- arguments[3]
meth_array <- arguments[4]
output_extension <- arguments[5] 


library(ggplot2)
library(ggrepel)
library(data.table)
library(grid)
library(dplyr)
library(meffil)

# function (taken from meffil github) to perform scatter thinning for plot
scatter.thinning <- function(x,y,resolution=100,max.per.cell=100) {
  x.cell <- floor((resolution-1)*(x - min(x,na.rm=T))/diff(range(x,na.rm=T))) + 1
  y.cell <- floor((resolution-1)*(y - min(y,na.rm=T))/diff(range(y,na.rm=T))) + 1
  z.cell <- x.cell * resolution + y.cell
  frequency.table <- table(z.cell)
  frequency <- rep(0,max(z.cell, na.rm=T))
  frequency[as.integer(names(frequency.table))] <- frequency.table
  f.cell <- frequency[z.cell]
  
  big.cells <- length(which(frequency > max.per.cell))
  sort(c(which(f.cell <= max.per.cell),
         sample(which(f.cell > max.per.cell),
                size=big.cells * max.per.cell, replace=F)),
       decreasing=F)
}

# function (adapted from meffil github) to create miami plot
meffil.ewas.miami.plot.1 <- function(results, sig.threshold=0.05/nrow(results_info),
                                     title="miami plot") {
  
  chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")
  chromosomes <- intersect(chromosomes, results$chromosome)
  stats <- results
  stats$chromosome <- factor(as.character(stats$chromosome), levels=chromosomes)
  stats$chr.colour <- "zero"
  stats$chr.colour[stats$chromosome %in% chromosomes[seq(1,length(chromosomes),2)]] <- "one"
  p.values <- stats$p.value
  p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin
  stats$stat <- -log(p.values,10) * sign(stats$beta)
  
  stats <- stats[order(stats$stat, decreasing=T),]
  
  chromosome.lengths <- sapply(chromosomes, function(chromosome)
    max(stats$position[which(stats$chromosome == chromosome)]))
  chromosome.lengths <- as.numeric(chromosome.lengths)
  chromosome.starts <- c(1,cumsum(chromosome.lengths)+1)
  names(chromosome.starts) <- c(chromosomes, "NA")
  stats$global <- stats$position + chromosome.starts[stats$chromosome] - 1
  
  selection.idx <- scatter.thinning(stats$global, stats$stat,
                                    resolution=100, max.per.cell=100)
  
  stats_for_plot <- head(
    stats[order(stats$p.value), ][stats$p.value < sig.threshold, ],
    50
  )
  
  stats_for_plot$gene <- factor(stats_for_plot$gene)
  
  (ggplot(stats[selection.idx,], aes(x=position, y=stat)) +
      geom_point(aes(colour=chr.colour),size=0.8) +
      theme_classic()+
      scale_color_manual(values=c("#481567FF", "#404788FF"))+
      facet_grid(. ~ chromosome, space="free_x", scales="free_x") +
      theme(strip.text.x = element_text(angle = 90))+
      guides(colour=FALSE) +
      labs(x="Position",
           y=bquote(-log[10]("p-value") * sign(beta))) +             
      geom_hline(yintercept=log(sig.threshold,10), colour="red") +
      geom_hline(yintercept=-log(sig.threshold,10), colour="red") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      geom_label_repel(
        data = stats_for_plot,
        aes(label=gene),
        box.padding   = 1, 
        segment.color = 'red',max.overlaps=Inf)+
      coord_cartesian(clip = "off")+
      ggtitle(paste(plot.title))
  ) 
}

# load meta-analysis data and run miami plot
# Load EWAS results
ewas <- read.table(
  glint_ewas,
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

# Inspect
head(ewas)
str(ewas)

meth_array <- tolower(meth_array)

if (meth_array %in% c("450k", "epic", "epic2")) {
  probe.info <- meffil.get.features(meth_array)
} else {
  stop("Unknown methylation array: ", meth_array,
       ". Expected one of: 450k, epic, epic2")
}

plot.title <- paste(study_name)
print(plot.title)
results_info <- merge(ewas,probe.info,by.x="LMM.ID",by.y = "name")
print(head(results_info))
print(dim(results_info))
names(results_info)[names(results_info) == 'p.values'] <- 'p.value'
names(results_info)[names(results_info) == 'chromosome.y'] <- 'chromosome'


results_info <- results_info[order(results_info$p.value),]
print(head(results_info))
results_info$gene <- sapply(strsplit(results_info$gene.symbol,';'), "[", 1)
miami_out <- meffil.ewas.miami.plot.1(results_info,sig.threshold=0.05/nrow(results_info),
                                      title=paste0("Miami plot: ",plot.title))
jpeg(filename = paste0(output_path,"/",study_name,"_glint_ewas_miami_plot_",output_extension,".jpg"),width = 7, height = 4, units = "in", res = 600)
print(miami_out)
dev.off()

# code from meffil to run qq plots
meffil.ewas.qq.plot <- function(ewas.object,
                                sig.threshold=0.05/nrow(results_info),
                                sig.color="red",
                                title="",
                                xlab=bquote(-log[10]("expected p-values")),
                                ylab=bquote(-log[10]("observed p-values")),
                                lambda.method="median") {
  
  p.values <- sort(ewas.object$p.value, decreasing=T)
  p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin
  stats <- data.frame(is.sig=p.values < sig.threshold,
                      expected=-log(sort(ppoints(p.values),decreasing=T),10),
                      observed=-log(p.values, 10))
  lambda <- qq.lambda(p.values[which(p.values > sig.threshold)],
                      method=lambda.method)
  
  label.x <- min(stats$expected) + diff(range(stats$expected))*0.1
  label.y <- min(stats$expected) + diff(range(stats$observed))*0.9
  
  lambda.label <- paste("lambda == ", format(lambda$estimate,digits=3),
                        "%+-%", format(lambda$se, digits=3),
                        "~(", lambda.method, ")", sep="")
  
  selection.idx <- scatter.thinning(stats$observed, stats$expected,
                                    resolution=100, max.per.cell=100)
  
  lim <- range(c(0, stats$expected, stats$observed))
  sig.threshold <- format(sig.threshold, digits=3)
  
  (ggplot(stats[selection.idx,], aes(x=expected, y=observed)) + 
      geom_abline(intercept = 0, slope = 1, colour="black") +              
      geom_point(aes(colour=factor(sign(is.sig)))) +
      scale_colour_manual(values=c("black", "red"),
                          name="Significant",
                          breaks=c("0","1"),
                          labels=c(paste("p-value >", sig.threshold),
                                   paste("p-value <", sig.threshold))) +
      annotate(geom="text", x=label.x, y=label.y, hjust=0,
               label=lambda.label,
               parse=T) +
      xlim(lim) + ylim(lim) + 
      xlab(xlab) + ylab(ylab) +
      coord_fixed() +
      theme(legend.position="none"))# +
}

# function taken from meffil to calculate lambda
qq.lambda <- function(p.values, method="median", B=100) {
  stopifnot(method %in% c("median","regression","robust"))
  p.values <- na.omit(p.values)
  observed <- qchisq(p.values, df=1, lower.tail = FALSE)
  observed <- sort(observed)
  expected <- qchisq(ppoints(length(observed)), df=1, lower.tail=FALSE)
  expected <- sort(expected)
  
  lambda <- se <- NA
  if (method == "median")  {
    lambda <- median(observed)/qchisq(0.5, df=1)
    boot.medians <- sapply(1:B, function(i) median(sample(observed, replace=T)))
    se <- sd(boot.medians/qchisq(0.5,df=1))
  } else if (method %in% c("regression","robust")) {
    if (method == "regression")
      coef.table <- summary(lm(observed ~ 0 + expected))$coeff
    else
      coef.table <- summary(rlm(observed ~ 0 + expected))$coef
    lambda <- coef.table["expected",1]
    se <- coef.table["expected", "Std. Error"]
  }
  list(method=method, estimate=lambda, se=se)
}

# 
qq_out <- meffil.ewas.qq.plot(results_info,sig.threshold=0.05/nrow(results_info),
                              title=paste0(plot.title))

jpeg(filename = paste0(output_path,"/",study_name,"_glint_ewas_qq_plot_",output_extension,".jpg"),width = 4, height = 4, units = "in", res = 600)
print(qq_out)
dev.off()

