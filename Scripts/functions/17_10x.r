library(data.table)
vcffile<-"/mnt/vol1/projects/clingen/10X/24A/outs/phased_variants.vcf.gz"
my10x<-fread(cmd=paste("zcat",vcffile,"| grep -v '##'"),data.table=F)
vcffile2<-"/mnt/vol1/projects/clingen/Mutation/Mutation_SNP/FVG24A.GATK.snp.vcf.gz"
myill<-fread(cmd=paste("zcat",vcffile2,"| grep -v '##'"),data.table=F)
my10x<-my10x[my10x$FILTER == "PASS",]
my10x=my10x[nchar(my10x$ALT) == 1,] #eliminano le INDEL da più di una base
my10x=my10x[nchar(my10x$REF) == 1,] #eliminano le INDEL da più di una base
myill<-myill[grep("GL|MT", myill$"#CHROM", invert=T),] #eliminiamo gli scaffold
my10x$"#CHROM"<-gsub("chr","",my10x$"#CHROM")
my10x<-my10x[grep("Y", my10x$"#CHROM", invert=T),] #eliminiamo il cromosoma Y dai 10x
mytot=merge(my10x,myill, by=c("#CHROM", "POS"), sort=F, all=T)
#HANNO ALLINEATO SU DUE REFERENCE DIVERSE, RIP
