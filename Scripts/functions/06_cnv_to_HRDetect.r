# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
#Questo si fa sui file di MUTATION
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
make_option(c("-T", "--tumor"), type="character", default="/mnt/vol1/projects/clingen/Control-FREEC/output/FVG9T_chr.bam_CNVs",
              help="Input tumor file [default= %default]", metavar="character"),
make_option(c("-O", "--outfile"), type="character", default="/mnt/vol1/projects/clingen/Control-FREEC/cnv_for_HRDetect/FVG9T_for_HRDetect.txt", 
              help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$tumor)) {
  stop("WARNING: No input tumor file specified with '-T' flag.")
} else {  cat ("Input tumor file is ", opt$tumor, "\n")
  tumor <- opt$tumor  
  #setwd(wd_location)  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No output file specified with '-O' flag.")
} else {  cat ("Output file is ", opt$outfile, "\n")
  outfile <- opt$outfile
  #setwd(wd_location)  
  }

cnv_to_HRDetect<-function(tumor,outfile) #VIENE CARICATA MA NON ESEGUITA, PER ESEGUIRLA SERVE L'ULTIMA RIGA
{
#EVENTUALMENTE AGGIUSTARE PER L'INPUT AGGIORNATO
library("data.table")
mytumor<-fread(tumor,data.table=F,fill=T)
library(stringr)
mytumor$Bcount=str_count(mytumor[,6], pattern="B")
mytumor$Acount=str_count(mytumor[,6], pattern="A")
#mynormal$Bcount=str_count(mynormal[,6], pattern="B")
#mynormal$Acount=str_count(mynormal[,6], pattern="A")
mytumor$minor.copy.number.inTumour=mytumor$Bcount
mytumor$minor.copy.number.inTumour[mytumor$Acount<mytumor$Bcount]=mytumor$Acount[mytumor$Acount<mytumor$Bcount]
mytumor$minor.copy.number.inNormal=1
mytumor$total.copy.number.inNormal=2
setnames(mytumor, "V4", "total.copy.number.inTumour")
#teniamo solo i valori con incertezza <=5%
mytumor<-mytumor[mytumor$V7<=5,]
mytumor$seg_no<-c(1:nrow(mytumor)) #fatta dopo il filtro per avere numeri successivi
names(mytumor)[c(13,1,2,3,12,11,4,10)] <- c("seg_no", "Chromosome", "chromStart", "chromEnd", "total.copy.number.inNormal", "minor.copy.number.inNormal", "total.copy.number.inTumour", "minor.copy.number.inTumour")
mytumor<-mytumor[,c(13,1,2,3,12,11,4,10)]
write.table(mytumor, outfile, sep="\t", quote=F, row.names=F)
}
cnv_to_HRDetect(tumor=tumor,outfile=outfile)
