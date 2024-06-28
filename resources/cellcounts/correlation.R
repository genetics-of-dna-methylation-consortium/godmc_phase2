#Correlation matrix between measured cell counts and predicted cell counts
library(corrplot)

arguments <- commandArgs(T)
cellcounts_cov <- arguments[1]
measured_cellcounts <- arguments[2]
cor_matrix <- arguments [3]
cor_plot <- arguments [4]

predicted<-read.table(cellcounts_cov,header=T)
predicted$Bcells<-predicted$Bmem+predicted$Bnv
predicted$Tcells<-predicted$CD4Tmem+predicted$CD4Tnv+predicted$Treg+predicted$CD8Tmem+predicted$CD8Tnv
predicted<-subset(predicted,select=c("IID","Tcells","Bcells","Neu","Mono","Eos","Baso"))

measured<-read.table(measured_cellcounts,header=T)

#Check if IID column exist in measured cell counts file 
if (!"IID" %in% colnames(measured)) {
  stop("Set the name 'IID' to the column of individuals identifiers in measured cell count file")
}

data<-merge(measured,predicted,by="IID",all=F)
ids<-data$IID

measured <- measured[match(ids,measured$IID),]
predicted <- predicted[match(ids,predicted$IID),]

partial_names <- c("IID","neu", "mono","T","B","eos","baso")
selected_columns <- measured[, grepl(paste(partial_names, collapse = "|"), 
                                     colnames(measured), ignore.case = TRUE)]


exclude_column <- "IID"
suffix <- "_measured"

colnames(measured) <- ifelse(colnames(measured) != exclude_column, 
                             paste0(colnames(measured), suffix), colnames(measured))

suffix <- "_predicted"
colnames(predicted) <- ifelse(colnames(predicted) != exclude_column, 
                              paste0(colnames(predicted), suffix), colnames(predicted))


correlation_matrix <- cor(predicted[,-which(names(predicted) == "IID")],
                          measured[,-which(names(measured) == "IID")],
                         use = "complete.obs", method="spearman")

write.table(correlation_matrix, file=cor_matrix, row=TRUE, col=TRUE, qu=FALSE, sep="\t")
pdf(height=54,width=87,cor_plot)
#png(height=540, width=870,file=cor_plot)
corrplot(correlation_matrix, method = "circle", type = "full", tl.col = "black",tl.cex=5,cl.cex=5)
dev.off()
