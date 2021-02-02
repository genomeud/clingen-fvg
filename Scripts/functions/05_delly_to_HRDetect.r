# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SV",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-O", "--outdir"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SV/HRDetect", 
              help="Output directory [default= %default]", metavar="character")
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

if (is.null(opt$outdir)) {
  stop("WARNING: No output directory file specified with '-B' flag.")
} else {  cat ("Output directory is ", opt$outdir, "\n")
  outdir <- opt$outdir
  #setwd(wd_location)  
  }

delly_to_HRDetect<-function(indir,outdir) #VIENE CARICATA MA NON ESEGUITA, PER ESEGUIRLA SERVE L'ULTIMA RIGA
{
library("data.table")
starting_dir <- getwd()
setwd(indir)
myfiles<-dir(pattern=".gff")
for(aaa in 1:length(myfiles)) #l'originale è for(aaa in 1:length(myfiles))
{
mydelly<-fread(myfiles[aaa],data.table=FALSE,header=FALSE,sep="\t")
names(mydelly)<-c("chrom1","V2","svclass","start1","end2","V6","V7","V8","V9")
mydelly<-mydelly[,c("chrom1","svclass","start1","end2","V9")]
mydelly$chrom2<-unlist(lapply(strsplit(mydelly$V9, "TCHR="),"[",2))
mydelly$chrom2<-unlist(lapply(strsplit(mydelly$chrom2, ";"),"[",1))
mydelly$chrom2[mydelly$chrom2=="na"]<-mydelly$chrom1[mydelly$chrom2=="na"]
mydelly$start2<-unlist(lapply(strsplit(mydelly$V9, "TSTART="),"[",2))
mydelly$start2<-unlist(lapply(strsplit(mydelly$start2, ";"),"[",1))
mydelly$end2[mydelly$start2!="na"]<-mydelly$start2[mydelly$start2!="na"]
mydelly$start2[mydelly$start2=="na"]<-mydelly$end2[mydelly$start2=="na"]
mydelly$end1<-mydelly$start1
mydelly$dellySV2<-unlist(lapply(strsplit(mydelly$V9, "SVType="),"[",2))
mydelly$sample=unlist(lapply(strsplit(myfiles[aaa], "\\."), "[",1))
#sum(mydelly$start1[mydelly$dellySV2=="breakpoint"]-mydelly$end2[mydelly$dellySV2=="breakpoint"]) #VERIFICA CHE I BREAKPOINT COMINCINO E FINISCONO NELLA STESSA BASE, DEVE DARE = 0
mydelly<-mydelly[mydelly$dellySV2!="breakpoint",]
#chrom start end sample svclass
mydelly<-mydelly[,c("chrom1","start1","end1","chrom2","start2","end2","sample","svclass")]
mydelly$svclass[mydelly$svclass=="TRA"]<-"translocation"
mydelly$svclass[mydelly$svclass=="INV"]<-"inversion"
mydelly$svclass[mydelly$svclass=="DEL"]<-"deletion"
mydelly$svclass[mydelly$svclass=="DUP"]<-"tandem-duplication"
outfile=gsub("gff", "bed", myfiles[aaa])
write.table(mydelly, paste(outdir, outfile, sep="/"), row.names=FALSE, quote=FALSE, sep="\t")
}
setwd(starting_dir)
}
delly_to_HRDetect(indir=indir,outdir=outdir)
