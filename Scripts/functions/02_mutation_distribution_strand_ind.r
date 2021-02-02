# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV",
              help="Input directory", metavar="character"),
  make_option(c("-B", "--barplotfile"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Plots/mutation_distribution_strand_independent.png", 
              help="Barplot distribution of mutation [default= %default]", metavar="character"),
  make_option(c("-T", "--tablefile"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Tables/mutation_distribution_strand_independent.txt", 
              help="Table with mutation distribution [default= %default]", metavar="character")
  # make_option(c("-C", "--control_prefix"), type="character", default="control_",
              # help="Prefix identifying control samples", metavar="character"),
  # make_option(c("-R", "--raw_counts"), type="character", default=NULL,
              # help="raw (total) read counts for this starting file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No input directory file specified with '-I' flag.")
} else {  cat ("Input directory is ", opt$indir, "\n")
  indir <- opt$indir  
  #setwd(wd_location)  
  }

if (is.null(opt$barplotfile)) {
  stop("WARNING: No barplot file specified with '-B' flag.")
} else {  cat ("Barplot file is ", opt$barplotfile, "\n")
  barplotfile <- opt$barplotfile  
  #setwd(wd_location)  
  }

  if (is.null(opt$tablefile)) {
  stop("WARNING: No table file specified with '-T' flag.")
} else {  cat ("Table file is ", opt$tablefile, "\n")
  tablefile <- opt$tablefile  
  #setwd(wd_location)  
  }

mutation_distr<-function(indir,barplotfile,tablefile)
{
library("data.table")
library("gplots")
starting_dir <- getwd()
setwd(indir)
myfiles<-dir(pattern=".vcf$")
myfiles<-myfiles[myfiles != "FVG3T.muTect.somatic.snv.vcf"]
png(barplotfile,width=30,height=30,units="cm",res=600,type="cairo")
par(mfrow=c(5,5))
for(aaa in 1:length(myfiles))
{
tm<-scan(myfiles[aaa],what="",sep="\n")
toskip<-max(grep("##",tm))
mysnp<-fread(myfiles[aaa],skip=toskip,sep="\t",data.table=F)
mytum<-unlist(strsplit(myfiles[aaa],"\\."))[1]
mynorm<-gsub("T","A",mytum)
#tum<-unlist(strsplit(myfiles[1],"\\."))[1]
#Columns 10 and 11 have sample names. Sometemies the first is tumor sometimes normal. We use grep to understand the order
myt<-grep("T",names(mysnp)[10:11])
myn<-grep("A",names(mysnp)[10:11])
mysnp$tumor<-unlist(lapply(strsplit(mysnp[,mytum],":"),"[",1))
mysnp$normal<-unlist(lapply(strsplit(mysnp[,mynorm],":"),"[",1))
mysnp$mutation<-NA
mysnp$mutation[mysnp$tumor=="0/1"&mysnp$normal=="0"]<-paste(mysnp$REF,mysnp$ALT,sep=">")
mysnp$mutation[mysnp$mutation=="A>C"]<-"T>G"
mysnp$mutation[mysnp$mutation=="A>G"]<-"T>C"
mysnp$mutation[mysnp$mutation=="A>T"]<-"T>A"
mysnp$mutation[mysnp$mutation=="G>T"]<-"C>A"
mysnp$mutation[mysnp$mutation=="G>C"]<-"C>G"
mysnp$mutation[mysnp$mutation=="G>A"]<-"C>T"
toplot<-round(100*table(mysnp$mutation)/sum(table(mysnp$mutation)),3)
toplot<-toplot[c("T>G","T>C","T>A","C>A","C>G","C>T")]
mycol<-c("deepskyblue","black","firebrick","gainsboro","darkolivegreen3","coral")
#png(mybarplotfile,width=6,height=6,units="cm",res=600,type="cairo")
par(mar=c(2.1,1.5,1.6,0.2),mgp=c(1.7, 0.7, 0))
barplot2(toplot,cex.names=0.6,cex.axis=0.6,las=1,col=mycol,ylim=c(0,40),main=mytum)
#dev.off()
}
dev.off()
setwd(starting_dir)
}
mutation_distr(indir=indir,barplotfile=barplotfile,tablefile=tablefile)
