# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV",
              help="Input directory [default= %default]", metavar="character"),
  make_option(c("-B", "--barplotfile"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Plots/mutation_distribution.png", 
              help="Barplot distribution of mutation [default= %default]", metavar="character"),
  make_option(c("-T", "--tablefile"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Tables/mutation_distribution.txt", 
              help="Table with mutation distribution [default= %default]", metavar="character") #IN REALTA NON VIENE CREATA
			  
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

mutation_distr<-function(indir,barplotfile,tablefile) #VIENE CARICATA MA NON ESEGUITA, PER ESEGUIRLA SERVE L'ULTIMA RIGA
{
library("data.table")
library("gplots")
starting_dir <- getwd()
setwd(indir)
myfiles<-dir(pattern=".vcf")
#myfiles=grep("FVG3T",myfiles,invert=T,value=T) #perchè il 3 è falsato. Senza value = T darebbe come risultato il numero di elementi che soddisfa la condizione (NON è COME LINUX)
png(barplotfile,width=30,height=30,units="cm",res=600,type="cairo")
par(mfrow=c(7,4)) # disposizione dei barplot, quante righe e quante colonne
for(aaa in 1:length(myfiles)) #l'originale è for(aaa in 1:length(myfiles))
{
tm<-scan(myfiles[aaa],what="",sep="\n")
toskip<-max(grep("##",tm)) #si lascia la riga con i nomi delle colonne perchè viene letta come tale dal comando successivo
mysnp<-fread(myfiles[aaa],skip=toskip,sep="\t",data.table=F)#CREA IL FILE SENZA INTESTAZIONE
mytum<-unlist(strsplit(myfiles[aaa],"\\."))[1]#LISTA DEI NOMI CAMPIONI TUMORALI
mynorm<-gsub("T","A",mytum)#LISTA CAMPIONI SANI (CAMBIA LA T CON LA A A QUELLI DI PRIMA)
#tum<-unlist(strsplit(myfiles[1],"\\."))[1]
#Columns 10 and 11 have sample names. Sometemies the first is tumor sometimes normal. We use grep to understand the order
myt<-grep("T",names(mysnp)[10:11]) 
myn<-grep("A",names(mysnp)[10:11])
mysnp$tumor<-unlist(lapply(strsplit(mysnp[,mytum],":"),"[",1)) #PRENDONO LO 0/1 DA TUMOR E LO METTONO IN UNA COLONNA
mysnp$normal<-unlist(lapply(strsplit(mysnp[,mynorm],":"),"[",1))#PRENDONO LO 0 DA NORMAL
mysnp$mutation<-NA
mysnp$mutation[mysnp$tumor=="0/1"&mysnp$normal=="0"]<-paste(mysnp$REF,mysnp$ALT,sep=">") #VEDERE COSA VUOL DIRE CHE IL NORMALE VALE SEMPRE 0
toplot<-round(100*table(mysnp$mutation)/sum(table(mysnp$mutation)),3) #SI VEDE LA PERCENTUALE DI CIASCUNA MUTAZIONE
toplot<-toplot[c("A>C","T>G","A>G","T>C","A>T","T>A","C>A","G>T","C>G","G>C","C>T","G>A")] #CAMBIA L'ORDINE E BASTA?
mycol<-c("deepskyblue","dodgerblue","black","gray45","firebrick","firebrick1","gainsboro","ghostwhite","darkolivegreen3","darkolivegreen4","coral","coral1")
#png(mybarplotfile,width=6,height=6,units="cm",res=600,type="cairo")
#par(mar=c(2,1.4,1.5,0.1),mgp=c(1.6, 0.6, 0))
barplot2(toplot,cex.names=0.3,cex.axis=0.3,las=1,col=mycol,ylim=c(0,20),main=mytum)
#dev.off()
}
dev.off()
setwd(starting_dir)
}
mutation_distr(indir=indir,barplotfile=barplotfile,tablefile=tablefile)
