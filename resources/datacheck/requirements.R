arguments <- commandArgs(T)
related <- as.character(arguments[1])


message("Checking R version")
currentr <- paste0(R.Version()['major'], ".", R.Version()['minor'])
ch <- compareVersion(currentr, "4.0")
if(ch == -1)
{
	stop("You are running R version ", currentr, ". Please upgrade to at least 4.0.")
}



message("Checking that all required packages are present")

pkglist <- c(
	"lattice",
	"ggplot2",
	"data.table",
	"MatrixEQTL",
	"parallel",
	"matrixStats",
	"plyr",
	"Cairo",
	"plotrix",
	"meffil",
	"EasyQC",
	"impute",
	"EpiDISH",
	"ewaff",
	"MASS",
	"dplyr",
	"magrittr",
	"ggrepel",
  "DunedinPACE",
  "RPMM",
  "vioplot",
  "qqman",
  "ggtext",
  "diptest",
  "mixtools",
  "ggpubr",
  "corrplot",
  "janitor",
  "tidyr",
  "readr",
  "scoreInvHap",
  "EpiDISH"
)

index <- pkglist %in% rownames(installed.packages())
if(any(!index))
{
	stop("Before continuing, the following packages need to be installed:\n", paste(pkglist[!index], collapse="\n"))
} else {
	message("All required packages installed")
}

pkglist_related<-c("SNPRelate","GENESIS","GWASTools")

if (related=="yes")
{
message("Checking that all required packages are present for related samples")

index <- pkglist_related %in% rownames(installed.packages())

if(any(!index))
{
	stop("Before continuing, the following packages need to be installed:\n", paste(pkglist_related[!index], collapse="\n"))
} else {
	message("All required packages for related samples are installed")
}
}

l<-list.files("./resources/genetics",pattern="HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz")
if(length(l)==0) 
{
    stop("Before continuing, you need to download HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz from the sftp \n")    
}
l<-list.files("./resources/bin/hase/data",pattern="ref-hrc.ref.gz")
if(length(l)==0) 
{
    stop("Before continuing, you need to download ref-hrc.ref.gz from the sftp \n")    
}

l<-list.files("./resources/bin/hase/data",pattern="ref-hrc.ref_info.h5")
if(length(l)==0)
{
    stop("Before continuing, you need to download ref-hrc.ref_info.h5 from the sftp \n")
}


library(meffil)

y<-packageVersion("meffil")
if(y<"1.3.8")
{
stop ("Meffil warning: please update to latest version")
}

message("All required packages are installed and required files are downloaded \n\n")

library(EpiDISH)
data(cent12CT.m)
if(exists("cent12CT.m")==F)
{
stop ("EpiDISH warning: please install EpiDISH version from github (https://github.com/sjczheng/EpiDISH)")
}

library(DunedinPACE)

x<-packageVersion("DunedinPACE")
if(x<"0.99.0")
{
	stop("DunedinPACE warning: please install DunedinPACE version from github (https://github.com/danbelsky/DunedinPACE) but with remotes package")
}

