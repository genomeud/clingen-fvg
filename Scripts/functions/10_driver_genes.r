#Estraiamo le informazioni sullo strand per CanDrA (comandi ritagliati dallo script 03)
library("data.table")
library("gplots")
library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
starting_dir <- getwd()
setwd("/mnt/vol1/projects/clingen/Somatic/Somatic_SNV")
vcf_files<-dir(pattern=".vcf$") #il $ alla fine fa in modo che prenda solo i file che FINISCONO con vcf
sample_names<-unlist(lapply(strsplit(vcf_files,"\\."),"[",1))#per separare il punto serve il doppio backslash, la quadra serve per lapply
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)#legge tutti i vcf come g-ranges, assegnandogli un nome e un genoma di riferimento
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Ora bisogna estrarre i dati sullo strand dal .gff
gff<-fread("/mnt/vol1/projects/clingen/Reference/hg19.ensGene.gtf.gz", data.table=F, sep="\t")
transgff<-gff[gff$V3=="transcript",] #estrae solo i dati sulle porzioni trascritte

#carichiamo in memoria due funzioni che servono dopo
count.overlap<-function(start.1=1,end.1=4845859,start.2=1,end.2=84945)
{
	max.start<-max(start.1,start.2)
	min.end<-min(end.1,end.2)
	overlap<-max(0,(1+min.end-max.start))
	overlap
}

count.overlap.vector<-function(x)
{
	max.start<-max(x[1],x[3])
	min.end<-min(x[2],x[4])
	overlap<-max(0,(1+min.end-max.start))
	overlap
}

#Get transcriptional strand information for all positions in the first VCF object with mut_strand
for (bbb in 1:length(vcfs))
{ 
bname=names(vcfs)[bbb]
strand = mut_strand(vcfs[[bbb]], genes_hg19)

#Ora bisogna estrarre i dati dai GRanges
mydf<-data.frame(chromosome=seqnames(vcfs[[bbb]]),pos=ranges(vcfs[[bbb]]), strand=strand(vcfs[[bbb]]), ref=elementMetadata(vcfs[[bbb]])$REF, alt=elementMetadata(vcfs[[bbb]])$ALT)
mydf$strand<-strand
mydf$strand<-as.character(mydf$strand)
mydf$strand[mydf$strand=="-"]<-"."


#ora vediamo quali porzioni di mydf sono comprese nei geni e se sono + o -

for (aaa in 1:nrow(mydf))
{
if (mydf$strand[aaa]==".") next
cat("line", aaa, "out of", nrow(mydf),"\n")
stransgff<-transgff[transgff$V1==mydf$chromosome[aaa],]
tdf<-cbind(rep(mydf$"pos.start"[aaa], nrow(stransgff)),rep(mydf$"pos.end"[aaa], nrow(stransgff)),stransgff$V4, stransgff$V5)
overlap<-apply(tdf,1,count.overlap.vector)
stra<-stransgff$V7[overlap>0]
strap<-sum(stra=="+")
stram<-sum(stra=="-")
if (stram>strap) mydf$strand[aaa]<-"-" else mydf$strand[aaa]<-"+"
}

# Per l'input di CanDrA servono:
# -numero cromosoma
# -start/end
# -allele reference
# -allele mutato
# strand

#Questa parte per ora è commentata
# candradf<-data.frame(mydf$chromosome, mydf$pos.start, mydf$ref, mydf$alt.value, mydf$strand)
# candradf$mydf.chromosome<-gsub("chr","",candradf$mydf.chromosome)
# write.table(candradf, paste("/mnt/vol1/projects/clingen/driver_mutations/dfs_for_candra/", bname,".txt", sep =""), sep="\t", col.names=F, row.names=F, quote=F)
# }

#Senza filtrare per il p-value vediamo che VCX3B, che a noi risulta solo in FVG21, si trova mutato anche in 20 e 27 come ci dicono dal CRO. Si tratta di mutazioni diverse (il sito è leggermente diverso) e quella del 21 secondo CanDrA è più candidata ad essere una driver.
#va detto che CanDrA fa semplicemente un confronto con il suo database di mutazioni, non va a vedere cosa succede nei nostri campioni.
#TP53 c'è nel database però non l'ha trovata: perchè? Cerchiamo di scoprirlo

candradf2<-data.frame(mydf$chromosome, mydf$pos.start, mydf$ref, mydf$alt.value, strand=".")
candradf2$mydf.chromosome<-gsub("chr","",candradf2$mydf.chromosome)
write.table(candradf, paste("/mnt/vol1/projects/clingen/driver_mutations/dfs_for_candra/", bname,".txt", sep =""), sep="\t", col.names=F, row.names=F, quote=F)
}

#ora si usa CanDrA su linux

setwd(starting_dir)



