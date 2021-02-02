#Partiamo dai file con gli SNP somatici (vcf), usando count.overlap prendiamo solo le mutazioni in cui non c'è overlap con i trascritti
library(data.table)
#carichiamo in memoria le due funzioni che servono dopo
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

#carichiamo il gff della reference
gff<-fread("/mnt/vol1/projects/clingen/Reference/hg19.ensGene.gtf.gz", data.table=F, sep="\t")
transgff<-gff[gff$V3=="transcript",] #estrae solo i dati sulle porzioni trascritte
#legge i vcf saltando le righe di intestazione E PRENDE SOLO LE INTERGENIC
vcfs<-dir("/mnt/vol1/projects/clingen/Somatic/Somatic_SNV", pattern=".vcf$", full.names=T)
####NON COPIARE QUESTO COME BASE PER ALTRI SCRIPT, USARE QUELLO DELLE INDEL####
for (aaa in 1:length(vcfs))
{
myfile<-fread(cmd=paste("grep -v '##'", vcfs[aaa]),data.table=F)
myfile$intergenic<-0
mychr<-unique(myfile$"#CHROM")
mychror<-mychr[1:24]
mychr<-paste("chr",mychror,sep="")
for (bbb in 1:length(mychr))
{
stransgff<-transgff[transgff$V1==mychr[bbb],]
myfilechr<-myfile[myfile$"#CHROM"==mychror[bbb],]
for (ccc in 1:nrow(myfilechr))
{
tdf<-cbind(rep(myfilechr$POS[ccc], nrow(stransgff)),rep(myfilechr$POS[ccc], nrow(stransgff)),stransgff$V4, stransgff$V5)
overlap<-apply(tdf,1,count.overlap.vector)
if (sum(overlap)== 0) myfilechr$intergenic[ccc]<-1
}
ifelse (bbb==1, isnp<-myfilechr[myfilechr$intergenic==1,], isnp<-rbind(isnp, myfilechr[myfilechr$intergenic==1,]))
}
write.table(isnp, paste("/mnt/vol1/projects/clingen/intergenic_mutations/SNP/", gsub(".muTect.somatic.snv.vcf",".txt",basename(vcfs[aaa])), sep = ""), quote=F, row.names=F, sep="\t")
}

#ora carichiamo tutti i file tranne il 3 per fare un rbind
out<-dir("/mnt/vol1/projects/clingen/intergenic_mutations/SNP", pattern="T.txt$", full.names=T)
out<-out[out!="/mnt/vol1/projects/clingen/intergenic_mutations/SNP/FVG3T.txt"] #togliamo FVG3

for (aaa in 1:length(out))
{
samplename<-gsub(".txt","",basename(out[aaa]))
read<-fread(out[aaa], data.table=F)
read<-read[,1:5]
read$sample<-samplename
ifelse (aaa==1, totalread<-read, totalread<-rbind(totalread,read))
}

#paste di cromosoma e posizione, colonna nuova
totalread$chrpos<-paste(totalread[,1],totalread[,2], sep="_")
compacttotal<-aggregate(totalread$sample, by=list(totalread$chrpos), FUN="paste", sep=";")


#ha funzionato ma dobbiamo sudare un altro po'
lapply(compacttotal$x, length) #da la lunghezza delle liste in compacttotal
compacttotal<-compacttotal[lapply(compacttotal$x, length)>1,] #scartiamo tutti quelli di l=1
compacttotal<-compacttotal[order(unlist(lapply(compacttotal$x, length)),decreasing=T),] #ordiniamo in ordine decrescente
compacttotal$samples<-unlist(lapply(compacttotal$x,paste,sep=";", collapse=";"))
compacttotal$x<-NULL
colnames(compacttotal)[1]<-"chr_pos"

write.table(compacttotal, "/mnt/vol1/projects/clingen/intergenic_mutations/SNP.txt", row.names=F, quote=F, sep="\t")

#In futuro vedremo per tutte le regioni in cui c'è una concentrazione di mutazioni se queste sono conservate o se intorno ci sono geni importanti

#Quanto fatto fino ad ora tiene conto soltanto della posizione (la mutazione potrebbe essere diversa), quindi si potrebbe rifare la stessa cosa ma tenendo conto della specifica mutazione
# totalread$chrposmut<-paste(totalread[,1],totalread[,2],totalread[,5], sep="_")
# compacttotal2<-aggregate(totalread$sample, by=list(totalread$chrposmut), FUN="paste", sep=";")
# lapply(compacttotal2$x, length) #da la lunghezza delle liste in compacttotal
# compacttotal2<-compacttotal2[lapply(compacttotal2$x, length)>1,] #scartiamo tutti quelli di l=1
# compacttotal2<-compacttotal2[order(unlist(lapply(compacttotal2$x, length)),decreasing=T),] #ordiniamo in ordine decrescente

##alla fine i risultati sono praticamente uguali, perchè le transizioni sono molto più comuni delle transversioni##

###Ora partiamo con le piccole INDEL###
#legge i vcf saltando le righe di intestazione
vcfs<-dir("/mnt/vol1/projects/clingen/Somatic/Somatic_INDEL", pattern=".vcf$", full.names=T)
vcfs<-vcfs[vcfs!="/mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/FVG3T.Strelka.somatic.indel.vcf"]

for (aaa in 1:length(vcfs))
{
myfile<-fread(cmd=paste("grep -v '##'", vcfs[aaa]),data.table=F)
myfile$intergenic<-0
mychr<-unique(myfile$"#CHROM")
endchrom<-min(24,length(mychr))
mychror<-mychr[1:endchrom]
mychr<-paste("chr",mychror,sep="")
for (bbb in 1:length(mychr))
{
stransgff<-transgff[transgff$V1==mychr[bbb],]
myfilechr<-myfile[myfile$"#CHROM"==mychror[bbb],]
for (ccc in 1:nrow(myfilechr))
{
tdf<-cbind(rep(myfilechr$POS[ccc], nrow(stransgff)),rep(myfilechr$POS[ccc], nrow(stransgff)),stransgff$V4, stransgff$V5)
overlap<-apply(tdf,1,count.overlap.vector)
if (sum(overlap)== 0) myfilechr$intergenic[ccc]<-1
}
ifelse (bbb==1, isnp<-myfilechr[myfilechr$intergenic==1,], isnp<-rbind(isnp, myfilechr[myfilechr$intergenic==1,]))
}
write.table(isnp, paste("/mnt/vol1/projects/clingen/intergenic_mutations/INDEL/", gsub(".Strelka.somatic.indel.vcf",".txt",basename(vcfs[aaa])), sep = ""), quote=F, row.names=F, sep="\t")
}

out<-dir("/mnt/vol1/projects/clingen/intergenic_mutations/INDEL", pattern="T.txt$", full.names=T)

for (aaa in 1:length(out))
{
samplename<-gsub(".txt","",basename(out[aaa]))
read<-fread(out[aaa], data.table=F)
read<-read[,1:5]
read$sample<-samplename
ifelse (aaa==1, totalread<-read, totalread<-rbind(totalread,read))
}

#paste di cromosoma e posizione, colonna nuova
totalread$chrpos<-paste(totalread[,1],totalread[,2], sep="_")
compacttotal<-aggregate(totalread$sample, by=list(totalread$chrpos), FUN="paste", sep=";")


#ha funzionato ma dobbiamo sudare un altro po'
compacttotal<-compacttotal[lapply(compacttotal$x, length)>1,] #scartiamo tutti quelli di l=1
##Non resta nessuna INDEL in comune tra i campioni. Andiamo a vedere le INDEL vicine tra loro
distance<-diff(totalread$POS)
distance<-which(0<distance & distance<50) #le posizioni interessate sono n e n+1. Abbiamo considerato 50 come il limite per dire che due indel si trovano nella stessa posizione
compacttotal2<-totalread[unique(c(distance, distance +1)),]
compacttotal2<-compacttotal2[order(as.numeric(compacttotal2[,1],compacttotal2[,2])),]

write.table(compacttotal2, "/mnt/vol1/projects/clingen/intergenic_mutations/INDEL.txt", row.names=F, quote=F, sep="\t")
###In realtà sono tutte sostituzioni e quindi le prende come inserzione + delezione, quindi non sono importanti, non ci serve nulla da questa tabella###

#Ora cominciamo con le CNV
#legge i gff saltando le righe di intestazione
gffs<-dir("/mnt/vol1/projects/clingen/Somatic/Somatic_CNV", pattern=".gff$", full.names=T)

for (aaa in 1:length(gffs))
{
myfile<-fread(cmd=paste("grep -v '##'", gffs[aaa]),data.table=F, header=F)
myfile<-myfile[grep("breakpoint", myfile$V9, invert = T),]
myfile$intergenic<-0
mychr<-unique(myfile$V1)
endchrom<-min(24,length(mychr))
mychror<-mychr[1:endchrom]
mychr<-paste("chr",mychror,sep="")
for (bbb in 1:length(mychr))
{
stransgff<-transgff[transgff$V1==mychr[bbb],]
myfilechr<-myfile[myfile$V1==mychror[bbb],]
for (ccc in 1:nrow(myfilechr))
{
tdf<-cbind(rep(myfilechr$V4[ccc], nrow(stransgff)),rep(myfilechr$V5[ccc], nrow(stransgff)),stransgff$V4, stransgff$V5)
overlap<-apply(tdf,1,count.overlap.vector)
if (sum(overlap)== 0) myfilechr$intergenic[ccc]<-1
}
ifelse (bbb==1, isnp<-myfilechr[myfilechr$intergenic==1,], isnp<-rbind(isnp, myfilechr[myfilechr$intergenic==1,]))
}
write.table(isnp, paste("/mnt/vol1/projects/clingen/intergenic_mutations/CNV/", gsub(".final.bam_CNVs.final.gff",".txt",basename(gffs[aaa])), sep = ""), quote=F, row.names=F, sep="\t")
}

#ora carichiamo tutti i file tranne il 3 per fare un rbind
out<-dir("/mnt/vol1/projects/clingen/intergenic_mutations/CNV", pattern="T.txt$", full.names=T)
out<-out[out!="/mnt/vol1/projects/clingen/intergenic_mutations/CNV/FVG3T.txt"] #togliamo FVG3


for (aaa in 1:length(out))
{
cat(out[aaa],"\n")
samplename<-gsub(".txt","",basename(out[aaa]))
read<-fread(out[aaa], data.table=F)
if (nrow(read) == 0) next
read$sample<-samplename
read$cnvid<-unlist(lapply(strsplit(read$V9, ";"), "[", 3))
read$minstart<-0
read$maxend<-0
read$prop_overlap<-0
# paired<-read$cnvid[duplicated(read$cnvid)]
# read<-read[read$cnvid%in%paired,]
# for (bbb in 1:length(unique(read$cnvid)))
# {
# read$V4[read$cnvid==unique(read$cnvid)[bbb]]<-min(read$V4[read$cnvid==unique(read$cnvid)[bbb]])
# read$V5[read$cnvid==unique(read$cnvid)[bbb]]<-max(read$V5[read$cnvid==unique(read$cnvid)[bbb]])
# }
# read<-read[!duplicated(read$cnvid),]
if (aaa==1) 
{
totalread<-read
next
}
for (ccc in 1:nrow(totalread))
{
tempread<-read[read$V1==totalread$V1[ccc],]
tdf<-cbind(rep(totalread$V4[ccc], nrow(tempread)),rep(totalread$V5[ccc], nrow(tempread)), tempread$V4, tempread$V5)
overlap<-apply(tdf,1,count.overlap.vector)
if (sum(overlap)> 0) 
{
minstart<-apply(tdf, 1, "min")
maxend<-apply(tdf, 1, "max")
minstart<-min(tdf[which.max(overlap),])
maxend<-max(tdf[which.max(overlap),])
totalread$minstart[ccc]<-minstart
totalread$maxend[ccc]<-maxend
totalread$prop_overlap[ccc]<-overlap[which.max(overlap)]/(maxend-minstart)
if (totalread$prop_overlap[ccc] > 0.5) ##LAVORI IN CORSO
{
totalread$sample[ccc]<-paste(totalread$sample[ccc],samplename,sep=";")
} else {
totalread <- rbind(totalread,tempread)
totalread <- totalread[!duplicated(totalread),]
}
} else {
totalread <- rbind(totalread,tempread)
totalread <- totalread[!duplicated(totalread),]
}
}
}

write.table(totalread, "/mnt/vol1/projects/clingen/intergenic_mutations/CNV.txt", row.names=F, quote=F, sep="\t")

#Per la ricerca sulle mutazioni dare un'occhiata a CancerAtlas, The Cancer Genome Atlas, Genome Browser (senza troppe speranze quest'ultimo)

#summary(totalread$prop_overlap) per osservare i parametri
#Decidiamo di tenere come soglia di overlap il 20%

###DOBBIAMO RICONTROLLARE CHE NEL LOOP NON SI PERDANO LE CNV CHE SONO SOVRAPPOSTE TRA DUE CAMPIONI QUALSIASI CHE NON SIA IL PRIMO CHE ANALIZZIAMO###

##################SV TIME###################
#prendiamo soltanto le delezioni

svs<-dir("/mnt/vol1/projects/clingen/Somatic/Somatic_SV", pattern=".gff$", full.names=T)

for (aaa in 1:length(svs))
{
myfile<-fread(cmd=paste("grep -v 'SVType=breakpoint'", svs[aaa]),data.table=F, header=F, sep="\t")
myfile<-myfile[grep("breakpoint", myfile$V9, invert = T),]
myfile$tchr <- unlist(lapply(strsplit(myfile$V9, "TCHR="), "[", 2))
myfile$tchr <- unlist(lapply(strsplit(myfile$tchr, ";"), "[", 1))
myfile$tstart <- unlist(lapply(strsplit(myfile$V9, "TSTART="), "[", 2))
myfile$tstart <- unlist(lapply(strsplit(myfile$tstart, ";"), "[", 1))
myfile<-myfile[,c(1:5, ncol(myfile)-1, ncol(myfile))]
myfile$intergenic<-0
mychr<-unique(myfile$V1)
endchrom<-min(24,length(mychr))
mychror<-mychr[1:endchrom]
mychr<-paste("chr",mychror,sep="")
for (bbb in 1:length(mychr))
{
stransgff<-transgff[transgff$V1==mychr[bbb],]
myfilechr<-myfile[myfile$V1==mychror[bbb],]
for (ccc in 1:nrow(myfilechr))
{
tdf<-cbind(rep(myfilechr$V4[ccc], nrow(stransgff)),rep(myfilechr$V5[ccc], nrow(stransgff)),stransgff$V4, stransgff$V5)
overlap<-apply(tdf,1,count.overlap.vector)
if (myfile$tchr != "na")
{
ttransgff<-transgff[transgff$V1==paste("chr",myfilechr$tchr[ccc],sep=""),]
if (nrow(ttransgff) > 0)
{
tdf<-cbind(rep(myfilechr$V4[ccc], nrow(ttransgff)),rep(myfilechr$V5[ccc], nrow(ttransgff)),ttransgff$V4, ttransgff$V5)
overlap<-c(overlap, apply(tdf,1,count.overlap.vector))
}
}
if (sum(overlap)== 0) myfilechr$intergenic[ccc]<-1
}
ifelse (bbb==1, isnp<-myfilechr[myfilechr$intergenic==1,], isnp<-rbind(isnp, myfilechr[myfilechr$intergenic==1,]))
}
idel<-isnp[isnp$V3=="DEL",]
idup<-isnp[isnp$V3=="DUP",]
iinv<-isnp[isnp$V3=="INV",]
itra<-isnp[isnp$V3=="TRA",]
write.table(idel, paste("/mnt/vol1/projects/clingen/intergenic_mutations/SV/DEL/", gsub(".delly.somatic.sv.gff",".txt",basename(svs[aaa])), sep = ""), quote=F, row.names=F, sep="\t")
write.table(idup, paste("/mnt/vol1/projects/clingen/intergenic_mutations/SV/DUP/", gsub(".delly.somatic.sv.gff",".txt",basename(svs[aaa])), sep = ""), quote=F, row.names=F, sep="\t")
write.table(iinv, paste("/mnt/vol1/projects/clingen/intergenic_mutations/SV/INV/", gsub(".delly.somatic.sv.gff",".txt",basename(svs[aaa])), sep = ""), quote=F, row.names=F, sep="\t")
write.table(itra, paste("/mnt/vol1/projects/clingen/intergenic_mutations/SV/TRA/", gsub(".delly.somatic.sv.gff",".txt",basename(svs[aaa])), sep = ""), quote=F, row.names=F, sep="\t")
}
#DA FARE IN SCREEN
for (folder in c("DEL", "DUP", "INV", "TRA"))
{
out<-dir(paste("/mnt/vol1/projects/clingen/intergenic_mutations/SV/",folder, sep=""), pattern="T.txt$", full.names=T)
out<-grep("FVG3T.txt", out, invert=T, value=T)
for (aaa in 1:length(out))
{
cat(out[aaa],"\n")
samplename<-gsub(".txt","",basename(out[aaa]))
read<-fread(out[aaa], data.table=F)
if (nrow(read)==0) next
read$sample<-samplename
read$minstart<-0
read$maxend<-0
read$prop_overlap<-0
if (aaa==1) 
{
totalread<-read
next
}
for (ccc in 1:nrow(totalread))
{
tempread<-read[read$V1==totalread$V1[ccc],]
tdf<-cbind(rep(totalread$V4[ccc], nrow(tempread)),rep(totalread$V5[ccc], nrow(tempread)), tempread$V4, tempread$V5)
overlap<-apply(tdf,1,count.overlap.vector)
if (sum(overlap)> 0) 
{
minstart<-apply(tdf, 1, "min")
maxend<-apply(tdf, 1, "max")
minstart<-min(tdf[which.max(overlap),])
maxend<-max(tdf[which.max(overlap),])
totalread$minstart[ccc]<-minstart
totalread$maxend[ccc]<-maxend
totalread$prop_overlap[ccc]<-overlap[which.max(overlap)]/(maxend-minstart)
if (totalread$prop_overlap[ccc] > 0.5)
{
totalread$sample[ccc]<-paste(totalread$sample[ccc],samplename,sep=";")
} else {
totalread <- rbind(totalread,tempread)
totalread <- totalread[!duplicated(totalread),]
}
} else {
totalread <- rbind(totalread,tempread)
totalread <- totalread[!duplicated(totalread),]
}
}
}
write.table(totalread, paste("/mnt/vol1/projects/clingen/intergenic_mutations/SV_",folder,".txt",sep=""), row.names=F, quote=F, sep="\t")
}

#Abbiamo trovato un breakpoint di traslocazione comune a 5 campioni che è in posizione 12:66451463 e in quella posizione c'è uno SNP di dbSNP che si chiama rs566655613


##############GRAFICI#############
library(stringr)
pdf("/mnt/vol1/projects/clingen/intergenic_mutations/plots/plot.pdf")
###SNP###
#Barplots per ogni tipo di mutazione si plotta quante sono le mutazioni con una certa numerosità:
snpdat<-fread("/mnt/vol1/projects/clingen/intergenic_mutations/SNP.txt",data.table=F)
snpdat$nsamples<-1+str_count(snpdat$samples,";")
snpdat$chrposplot<-paste("chr",gsub("_",":",snpdat$chr_pos),sep="")
#Non serve perché sono già ordinati per numero di campioni
#invdat2<-invdat[order(invdat$nsamples,decreasing=T),]
#Arbitrariamente scegliamo i primi 20

pp<-aggregate(snpdat$nsamples,by=list(snpdat$nsamples),FUN=length) #summary di tutte le entry uniche della colonna di riferimento
#barplot(invdat$nsamples,names=invdat$chrposplot)
#Questo sotto è più bello
barplot(pp$x,names=pp$Group.1, xlab="Numero di ripetizioni delle mutazioni", ylab="Quantità di mutazioni ripetute", main="SNP") #scala logaritmica per rendere i valori paragonabili

###INDEL, CNV, SV###
#LASCIAMO PERDERE SONO TROPPO POCHE E TUTTE IN UN SOLO CAMPIONE, AGGREGANDO QUELLE VICINE SONO POCHE COMUNQUE. VEDERE SOLO LA CNV CON 11 CAMPIONI
dev.off()

##############VEDIAMO QUALI DI QUESTE MUTAZIONI SONO VICINISSIME AI GENI############
gff<-fread("/mnt/vol1/projects/clingen/Reference/hg19.ensGene.gtf.gz", data.table=F, sep="\t")
transgff<-gff[gff$V3=="transcript",] #estrae solo i dati sulle porzioni trascritte
transgff$genename<-unlist(lapply(strsplit(transgff$V9, ";"), "[", 1))
transgff$genename<-gsub("gene_id ","",transgff$genename)
transgff$genename<-gsub("\"","",transgff$genename)
transgff<-transgff[,-9] 

snpdat$chr<-unlist(lapply(strsplit(snpdat$chrposplot,":"), "[", 1))
snpdat$pos<-unlist(lapply(strsplit(snpdat$chrposplot,":"), "[", 2))
snpdat$pos<-as.numeric(snpdat$pos)
snpdat$nearstart<-NA
snpdat$nearend<-NA
snpdat$genestart<-NA
snpdat$geneend<-NA
for(aaa in 1:nrow(snpdat))
{
ppstart<-which.min(abs(transgff$V4[transgff$V1 == snpdat$chr[aaa]]-snpdat$pos[aaa]))
ppend<-which.min(abs(transgff$V5[transgff$V1 == snpdat$chr[aaa]]-snpdat$pos[aaa]))
snpdat$nearstart[aaa]<-abs(transgff$V4[transgff$V1 == snpdat$chr[aaa]][ppstart]-snpdat$pos[aaa])
snpdat$nearend[aaa]<-abs(transgff$V5[transgff$V1 == snpdat$chr[aaa]][ppend]-snpdat$pos[aaa])
snpdat$genestart[aaa]<-transgff$genename[transgff$V1 == snpdat$chr[aaa]][ppstart]
snpdat$geneend[aaa]<-transgff$genename[transgff$V1 == snpdat$chr[aaa]][ppend]
cat(aaa, "\n")
}

#ORA SCRIVIAMO IL FILE CHE RIPORTA I GENI PIU VICINI E LA DISTANZA PER GLI SNP

write.table(snpdat, "/mnt/vol1/projects/clingen/intergenic_mutations/close_genes_SNP.txt", quote=F, sep="\t", row.names=F)

#CNV

cnvdat<-fread("/mnt/vol1/projects/clingen/intergenic_mutations/CNV.txt",data.table=F)
cnvdat$chr<-paste("chr",cnvdat$V1, sep="")
cnvdat$start<-as.numeric(cnvdat$V4)
cnvdat$end<-as.numeric(cnvdat$V5)
cnvdat$nearstart<-NA
cnvdat$nearend<-NA
cnvdat$upstream<-NA
cnvdat$downstream<-NA


#GLI 1 SI RIFERISCONO ALLO START DELLA CNV, I 2 ALL'END DELLA CNV

for(aaa in 1:nrow(cnvdat))
{
ppstart1<-which.min(abs(transgff$V4[transgff$V1 == cnvdat$chr[aaa]]-cnvdat$start[aaa]))
ppend1<-which.min(abs(transgff$V5[transgff$V1 == cnvdat$chr[aaa]]-cnvdat$start[aaa]))
ppstart2<-which.min(abs(transgff$V4[transgff$V1 == cnvdat$chr[aaa]]-cnvdat$end[aaa]))
ppend2<-which.min(abs(transgff$V5[transgff$V1 == cnvdat$chr[aaa]]-cnvdat$end[aaa]))
cnvdat$nearstart[aaa]<-abs(transgff$V4[transgff$V1 == cnvdat$chr[aaa]][ppstart2]-cnvdat$end[aaa]) #DOWNSTREAM
cnvdat$nearend[aaa]<-abs(transgff$V5[transgff$V1 == cnvdat$chr[aaa]][ppend1]-cnvdat$start[aaa])
cnvdat$upstream[aaa]<-transgff$genename[transgff$V1 == cnvdat$chr[aaa]][ppend1]
cnvdat$downstream[aaa]<-transgff$genename[transgff$V1 == cnvdat$chr[aaa]][ppstart2]
cat(aaa, "\n")
}

write.table(cnvdat, "/mnt/vol1/projects/clingen/intergenic_mutations/close_genes_CNV.txt", quote=F, sep="\t", row.names=F)

##INDEL###
invdat<-fread("/mnt/vol1/projects/clingen/intergenic_mutations/INDEL.txt",data.table=F)
invdat$chr<-paste("chr",unlist(lapply(strsplit(invdat$chrpos,"_"), "[", 1)), sep="")
invdat$pos<-as.numeric(unlist(lapply(strsplit(invdat$chrpos,"_"), "[", 2)))
invdat$nearstart<-NA
invdat$nearend<-NA
invdat$genestart<-NA
invdat$geneend<-NA
for(aaa in 1:nrow(invdat))
{
ppstart<-which.min(abs(transgff$V4[transgff$V1 == invdat$chr[aaa]]-invdat$pos[aaa]))
if (length(ppstart) == 0) next
ppend<-which.min(abs(transgff$V5[transgff$V1 == invdat$chr[aaa]]-invdat$pos[aaa]))
invdat$nearstart[aaa]<-abs(transgff$V4[transgff$V1 == invdat$chr[aaa]][ppstart]-invdat$pos[aaa])
invdat$nearend[aaa]<-abs(transgff$V5[transgff$V1 == invdat$chr[aaa]][ppend]-invdat$pos[aaa])
invdat$genestart[aaa]<-transgff$genename[transgff$V1 == invdat$chr[aaa]][ppstart]
invdat$geneend[aaa]<-transgff$genename[transgff$V1 == invdat$chr[aaa]][ppend]
cat(aaa, "\n")
}

write.table(invdat, "/mnt/vol1/projects/clingen/intergenic_mutations/close_genes_INDEL.txt", quote=F, sep="\t", row.names=F)

#SV###
deldat<-fread("/mnt/vol1/projects/clingen/intergenic_mutations/SV_DEL.txt",data.table=F)
dupdat<-fread("/mnt/vol1/projects/clingen/intergenic_mutations/SV_DUP.txt",data.table=F)
invdat<-fread("/mnt/vol1/projects/clingen/intergenic_mutations/SV_INV.txt",data.table=F)
tradat<-fread("/mnt/vol1/projects/clingen/intergenic_mutations/SV_TRA.txt",data.table=F)
deldat$chr<-paste("chr",deldat$V1, sep="")
deldat$start<-as.numeric(deldat$V4)
deldat$end<-as.numeric(deldat$V5)
deldat$nearstart<-NA
deldat$nearend<-NA
deldat$upstream<-NA
deldat$downstream<-NA

dupdat$chr<-paste("chr",dupdat$V1, sep="")
dupdat$start<-as.numeric(dupdat$V4)
dupdat$end<-as.numeric(dupdat$V5)
dupdat$nearstart<-NA
dupdat$nearend<-NA
dupdat$upstream<-NA
dupdat$downstream<-NA

invdat$chr<-paste("chr",invdat$V1, sep="")
invdat$start<-as.numeric(invdat$V4)
invdat$end<-as.numeric(invdat$V5)
invdat$nearstart<-NA
invdat$nearend<-NA
invdat$upstream<-NA
invdat$downstream<-NA

tradat$chr<-paste("chr",tradat$V1, sep="")
tradat$start<-as.numeric(tradat$V4)
tradat$end<-as.numeric(tradat$V5)
tradat$nearstart<-NA
tradat$nearend<-NA
tradat$upstream<-NA
tradat$downstream<-NA
tradat$nearstart_end<-NA
tradat$nearend_end<-NA
tradat$upstream_end<-NA
tradat$downstream_end<-NA
tradat$tchr<-paste("chr",tradat$tchr, sep="")

#GLI 1 SI RIFERISCONO ALLO START DELLA SV, I 2 ALL'END DELLA SV
for(aaa in 1:nrow(deldat))
{
ppend1<-which.min(abs(transgff$V5[transgff$V1 == deldat$chr[aaa]]-deldat$start[aaa]))
if (length(ppend1) == 0) next
ppstart2<-which.min(abs(transgff$V4[transgff$V1 == deldat$chr[aaa]]-deldat$end[aaa]))
deldat$nearstart[aaa]<-abs(transgff$V4[transgff$V1 == deldat$chr[aaa]][ppstart2]-deldat$end[aaa]) #DOWNSTREAM
deldat$nearend[aaa]<-abs(transgff$V5[transgff$V1 == deldat$chr[aaa]][ppend1]-deldat$start[aaa])
deldat$upstream[aaa]<-transgff$genename[transgff$V1 == deldat$chr[aaa]][ppend1]
deldat$downstream[aaa]<-transgff$genename[transgff$V1 == deldat$chr[aaa]][ppstart2]
cat(aaa, "\n")
}

for(aaa in 1:nrow(dupdat))
{
ppend1<-which.min(abs(transgff$V5[transgff$V1 == dupdat$chr[aaa]]-dupdat$start[aaa]))
if (length(ppend1) == 0) next
ppstart2<-which.min(abs(transgff$V4[transgff$V1 == dupdat$chr[aaa]]-dupdat$end[aaa]))
dupdat$nearstart[aaa]<-abs(transgff$V4[transgff$V1 == dupdat$chr[aaa]][ppstart2]-dupdat$end[aaa]) #DOWNSTREAM
dupdat$nearend[aaa]<-abs(transgff$V5[transgff$V1 == dupdat$chr[aaa]][ppend1]-dupdat$start[aaa])
dupdat$upstream[aaa]<-transgff$genename[transgff$V1 == dupdat$chr[aaa]][ppend1]
dupdat$downstream[aaa]<-transgff$genename[transgff$V1 == dupdat$chr[aaa]][ppstart2]
cat(aaa, "\n")
}

for(aaa in 1:nrow(invdat))
{
ppend1<-which.min(abs(transgff$V5[transgff$V1 == invdat$chr[aaa]]-invdat$start[aaa]))
if (length(ppend1) == 0) next
ppstart2<-which.min(abs(transgff$V4[transgff$V1 == invdat$chr[aaa]]-invdat$end[aaa]))
invdat$nearstart[aaa]<-abs(transgff$V4[transgff$V1 == invdat$chr[aaa]][ppstart2]-invdat$end[aaa]) #DOWNSTREAM
invdat$nearend[aaa]<-abs(transgff$V5[transgff$V1 == invdat$chr[aaa]][ppend1]-invdat$start[aaa])
invdat$upstream[aaa]<-transgff$genename[transgff$V1 == invdat$chr[aaa]][ppend1]
invdat$downstream[aaa]<-transgff$genename[transgff$V1 == invdat$chr[aaa]][ppstart2]
cat(aaa, "\n")
}

for(aaa in 1:nrow(tradat))
{
ppend1<-which.min(abs(transgff$V5[transgff$V1 == tradat$chr[aaa]]-tradat$start[aaa]))
ppend_end<-which.min(abs(transgff$V5[transgff$V1 == tradat$tchr[aaa]]-tradat$tstart[aaa]))
if (length(ppend1) == 0 | length(ppend_end) == 0) next
ppstart2<-which.min(abs(transgff$V4[transgff$V1 == tradat$chr[aaa]]-tradat$end[aaa]))
ppstart_end<-which.min(abs(transgff$V4[transgff$V1 == tradat$tchr[aaa]]-tradat$tstart[aaa]))
tradat$nearstart[aaa]<-abs(transgff$V4[transgff$V1 == tradat$chr[aaa]][ppstart2]-tradat$end[aaa]) #DOWNSTREAM
tradat$nearend[aaa]<-abs(transgff$V5[transgff$V1 == tradat$chr[aaa]][ppend1]-tradat$start[aaa])
tradat$upstream[aaa]<-transgff$genename[transgff$V1 == tradat$chr[aaa]][ppend1]
tradat$downstream[aaa]<-transgff$genename[transgff$V1 == tradat$chr[aaa]][ppstart2]
tradat$nearstart_end[aaa]<-abs(transgff$V4[transgff$V1 == tradat$tchr[aaa]][ppstart_end]-tradat$tstart[aaa])
tradat$nearend_end[aaa]<-abs(transgff$V5[transgff$V1 == tradat$tchr[aaa]][ppend_end]-tradat$tstart[aaa])
tradat$upstream_end[aaa]<-transgff$genename[transgff$V1 == tradat$tchr[aaa]][ppend_end]
tradat$downstream_end[aaa]<-transgff$genename[transgff$V1 == tradat$tchr[aaa]][ppstart_end]
cat(aaa, "\n")
}

write.table(deldat, "/mnt/vol1/projects/clingen/intergenic_mutations/close_genes_SVDEL.txt", quote=F, sep="\t", row.names=F)
write.table(dupdat, "/mnt/vol1/projects/clingen/intergenic_mutations/close_genes_SVDUP.txt", quote=F, sep="\t", row.names=F)
write.table(invdat, "/mnt/vol1/projects/clingen/intergenic_mutations/close_genes_SVINV.txt", quote=F, sep="\t", row.names=F)
write.table(tradat, "/mnt/vol1/projects/clingen/intergenic_mutations/close_genes_SVTRA.txt", quote=F, sep="\t", row.names=F)

#ENSG00000256614 



#ANDARE A VEDERE PER LA TESI I GENI CHE SONO DISTANTI MENO DI 1K BASI





