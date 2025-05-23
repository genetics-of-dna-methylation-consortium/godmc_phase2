suppressMessages(library(meffil))


arguments <- commandArgs(T)
beta_file <- arguments[1]
prop_var <- as.numeric(arguments[2])
max_pcs <- as.numeric(arguments[3])
phen_file <- arguments[4]
pc_out <- arguments[5]


message("Loading methylation data")
load(beta_file)


message("Extracting most variable probes and calculate PCs")
#featureset <- meffil:::guess.featureset(rownames(norm.beta))
#autosomal.sites <- meffil.get.autosomal.sites(featureset)
#autosomal.sites <- intersect(autosomal.sites, rownames(norm.beta))
autosomal.sites <- meffil.autosomal.subset(rownames(norm.beta))
norm.beta.aut <- norm.beta[autosomal.sites, ]

message("Calculating variances")
var.sites <- meffil.most.variable.cpgs(norm.beta.aut, n = 20000)
var.idx <- match(var.sites, rownames(norm.beta.aut))

message("Calculating beta PCs")
pc <- prcomp(t(meffil:::impute.matrix(norm.beta.aut[var.idx, ], margin = 1)))

message("Identifying PCs that cumulatively explain ", prop_var, " of variance")
cumvar <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
n_pcs <- which(cumvar > prop_var)[1]
message(n_pcs, " PCs required to explain ", prop_var, " of methylation variance")
n_pcs <- min(n_pcs, max_pcs)

pc_explain_var <- summary(pc)$importance[2,]
n_pc_0.1 <- sum(pc_explain_var > 0.01)
message(n_pc_0.1, " PCs could explain >0.01 variance individually")
n_pcs <- min(n_pcs, n_pc_0.1)

if(n_pcs > 0 ){
	pc <- pc$x[,1:n_pcs]
	pc <- as.matrix(pc)
}else{
	message("no methylation PC will be used")
	q()
}

if (ncol(pc) == 1){colnames(pc) <- "PC1"} 

if(phen_file != "NULL")
{
	message("Removing PCs associated with EWAS phenotypes")
	phen <- read.table(phen_file, he=T)
	rownames(phen) <- phen$IID
	phen <- subset(phen, IID %in% rownames(pc), select=-c(IID))
	pc1 <- pc[rownames(pc) %in% rownames(phen), ]
	phen <- phen[match(rownames(pc1), rownames(phen)), , drop=FALSE]
	stopifnot(all(rownames(phen) == rownames(pc1)))


	l <- lapply(1:ncol(phen), function(i)
	{
		pvals <- coefficients(summary(lm(phen[,i] ~ pc1)))[-1,4]
		which(pvals < 0.05/ncol(phen))
	})
	l <- sort(unique(unlist(l)))

	message("Identified ", length(l), " PC(s) associated with phenotypes")

	if(length(l) > 0)
	{
		pc <- pc[,! 1:ncol(pc) %in% l, drop=FALSE]
	}
}

pc1 <- t(pc)
save(pc, file=paste0(pc_out, ".RData"))
write.table(pc1, file=paste0(pc_out, ".txt"), row=T, col=T, qu=F, sep="\t")
