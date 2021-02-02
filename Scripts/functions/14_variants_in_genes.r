library(data.table)
##########
#INDEL
##########
nostriINDEL <- fread("/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsINDEL.txt", data.table=F) #l'elenco delle mutazioni annotate nei nostri campioni
#Prendiamo solo quelle che si trovano nei geni#
genicINDEL <- nostriINDEL[nostriINDEL$GeneName != ".",]
names(genicINDEL)[ncol(genicINDEL)] <- "sample"
#creiamo i summary#
#questo le prende tutte
mut_summary <- aggregate(genicINDEL$sample, by=list(genicINDEL$GeneName), FUN=function(x) (ngenes=length(unique(x))))
mut_summary <- mut_summary[order(mut_summary$x, decreasing=T),]

#questo prende le esoniche
exonINDEL <- genicINDEL[genicINDEL$Func == "exonic",]
mut_summary2 <- aggregate(exonINDEL$sample, by=list(exonINDEL$GeneName), FUN=function(x) (ngenes=length(unique(x))))
mut_summary2 <- mut_summary2[order(mut_summary2$x, decreasing=T),]

write.table(mut_summary, "/mnt/vol1/projects/clingen/driver_mutations/tabs/allINDEL_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary2, "/mnt/vol1/projects/clingen/driver_mutations/tabs/exoINDEL_summary.txt", quote=F, sep="\t", row.names=F)

##########
#CNV
#########

nostriCNV <- fread("/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsCNV.txt", data.table=F)
nostriCNV <- nostriCNV[nostriCNV$CNVType != "breakpoint",] #TENIAMO SOLO GAIN E LOSS
nostriCNV <- nostriCNV[nostriCNV$GeneName != ".",]
names(nostriCNV)[ncol(nostriCNV)] <- "sample"
mylist <- strsplit(nostriCNV$GeneName,",") #Estraiamo i geni interessati dall CNV
#associamo le liste al campione da cui provengono
for (aaa in 1:nrow(nostriCNV))
{
mylist[[aaa]] <- paste(nostriCNV$sample[aaa], mylist[[aaa]], sep=";")
}
myvector <- unique(unlist(mylist)) #togliamo i duplicati dovuti a possibili errori nell'annotazione delle CNV (sottoinsiemi)
allCNV <- data.frame(sample=unlist(lapply(strsplit(myvector, ";"), "[", 1)), gene=unlist(lapply(strsplit(myvector, ";"), "[", 2)))

mut_summary <- aggregate(allCNV$sample, by=list(allCNV$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary <- mut_summary[order(mut_summary$x, decreasing=T),]
##summaries##
#DOBBIAMO SEPARARE I GAIN DA I LOSS#
gainCNV<-nostriCNV[nostriCNV$CNVType == "gain",]
lossCNV<-nostriCNV[nostriCNV$CNVType == "loss",]
##Estraiamo i geni##
mylistg <- strsplit(gainCNV$GeneName,",") #Estraiamo i geni interessati dall CNV
mylistl <- strsplit(lossCNV$GeneName,",")
#associamo le liste al campione da cui provengono
for (aaa in 1:nrow(gainCNV))
{
mylistg[[aaa]] <- paste(gainCNV$sample[aaa], mylistg[[aaa]], sep=";")
}
for (aaa in 1:nrow(lossCNV))
{
mylistl[[aaa]] <- paste(lossCNV$sample[aaa], mylistl[[aaa]], sep=";")
}

myvectorg <- unique(unlist(mylistg)) #togliamo i duplicati dovuti a possibili errori nell'annotazione delle CNV (sottoinsiemi)
myvectorl <- unique(unlist(mylistl)) #togliamo i duplicati dovuti a possibili errori nell'annotazione delle CNV (sottoinsiemi)
gCNV <- data.frame(sample=unlist(lapply(strsplit(myvectorg, ";"), "[", 1)), gene=unlist(lapply(strsplit(myvectorg, ";"), "[", 2)))
lCNV <- data.frame(sample=unlist(lapply(strsplit(myvectorl, ";"), "[", 1)), gene=unlist(lapply(strsplit(myvectorl, ";"), "[", 2)))

mut_summary2 <- aggregate(gCNV$sample, by=list(gCNV$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary2 <- mut_summary2[order(mut_summary2$x, decreasing=T),]

mut_summary3 <- aggregate(lCNV$sample, by=list(lCNV$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary3 <- mut_summary3[order(mut_summary3$x, decreasing=T),]

write.table(mut_summary, "/mnt/vol1/projects/clingen/driver_mutations/tabs/allCNV_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary2, "/mnt/vol1/projects/clingen/driver_mutations/tabs/gainCNV_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary3, "/mnt/vol1/projects/clingen/driver_mutations/tabs/lossCNV_summary.txt", quote=F, sep="\t", row.names=F)

###########
#SV
###########
nostriSV <- fread("/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsSV.txt", data.table=F)
nostriSV <- nostriSV[nostriSV$SVType != "breakpoint",] #TENIAMO SOLO GAIN E LOSS
nostriSV <- nostriSV[nostriSV$GeneName != ".",]
names(nostriSV)[ncol(nostriSV)] <- "sample"
mylist <- strsplit(nostriSV$GeneName,",") #Estraiamo i geni interessati dall CNV
for (aaa in 1:nrow(nostriSV))
{
mylist[[aaa]] <- paste(nostriSV$sample[aaa], mylist[[aaa]], sep=";")
}
myvector <- unique(unlist(mylist)) #togliamo i duplicati dovuti a possibili errori nell'annotazione delle CNV (sottoinsiemi)
allSV <- data.frame(sample=unlist(lapply(strsplit(myvector, ";"), "[", 1)), gene=unlist(lapply(strsplit(myvector, ";"), "[", 2)))

mut_summary <- aggregate(allSV$sample, by=list(allSV$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary <- mut_summary[order(mut_summary$x, decreasing=T),]

#####
#SEPARIAMO LE SV PER TIPO
#####

#read gff#
gff<-fread("/mnt/vol1/projects/clingen/Reference/hg19.ensGene.gtf.gz", data.table=F, sep="\t")
transgff<-gff[gff$V3=="transcript",] #estrae solo i dati sulle porzioni trascritte

#count overlap##
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

#DOBBIAMO SEPARARE I GAIN DA I LOSS#
delSV<-nostriSV[nostriSV$SVType == "DEL",]
dupSV<-nostriSV[nostriSV$SVType == "DUP",]
invSV<-nostriSV[nostriSV$SVType == "INV",]
traSV<-nostriSV[nostriSV$SVType == "TRA",]
traSV$TGENE<-NA

#estraiamo gli endpoint delle traslocazioni#
for (aaa in nrow(traSV))
{
smallgff<-gff[gff$V1 == paste("chr", traSV$TCHR[aaa], sep=""),]
tdf<-cbind(rep(traSV$TSTART[aaa], nrow(smallgff)),rep(traSV$TSTART[aaa], nrow(smallgff)),smallgff$V4, smallgff$V5)
tdf2<-apply(tdf,2,as.numeric)
overlap<-apply(tdf2,1,count.overlap.vector)
if (sum(overlap) > 0) 
{
mygene<-gsub("gene_id ","", unlist(lapply(strsplit(smallgff$V9[overlap > 0], ";"), "[", 1)))
mygene<-gsub("\"", "", mygene)
traSV$TGENE[aaa]<-paste(mygene, collapse = ";")
}
}

##Estraiamo i geni##

mylistdel <- strsplit(delSV$GeneName,",") #Estraiamo i geni interessati dall CNV
mylistdup <- strsplit(dupSV$GeneName,",")
mylistinv <- strsplit(invSV$GeneName,",")
mylisttra <- strsplit(traSV$GeneName,",")
mylisttra2 <- strsplit(traSV$TGENE, ";")
mylisttra <- unique(c(mylisttra,mylisttra2)) 
#associamo le liste al campione da cui provengono
for (aaa in 1:nrow(delSV))
{
mylistdel[[aaa]] <- paste(delSV$sample[aaa], mylistdel[[aaa]], sep=";")
}
for (aaa in 1:nrow(dupSV))
{
mylistdup[[aaa]] <- paste(dupSV$sample[aaa], mylistdup[[aaa]], sep=";")
}
for (aaa in 1:nrow(invSV))
{
mylistinv[[aaa]] <- paste(invSV$sample[aaa], mylistinv[[aaa]], sep=";")
}
for (aaa in 1:nrow(traSV))
{
mylisttra[[aaa]] <- paste(traSV$sample[aaa], mylisttra[[aaa]], sep=";")
}

myvectordel <- unique(unlist(mylistdel)) #togliamo i duplicati dovuti a possibili errori nell'annotazione delle CNV (sottoinsiemi)
myvectordup <- unique(unlist(mylistdup))
myvectorinv <- unique(unlist(mylistinv))
myvectortra <- unique(unlist(mylisttra))
delSV <- data.frame(sample=unlist(lapply(strsplit(myvectordel, ";"), "[", 1)), gene=unlist(lapply(strsplit(myvectordel, ";"), "[", 2)))
dupSV <- data.frame(sample=unlist(lapply(strsplit(myvectordup, ";"), "[", 1)), gene=unlist(lapply(strsplit(myvectordup, ";"), "[", 2)))
invSV <- data.frame(sample=unlist(lapply(strsplit(myvectorinv, ";"), "[", 1)), gene=unlist(lapply(strsplit(myvectorinv, ";"), "[", 2)))
traSV <- data.frame(sample=unlist(lapply(strsplit(myvectortra, ";"), "[", 1)), gene=unlist(lapply(strsplit(myvectortra, ";"), "[", 2)))

mut_summary2 <- aggregate(delSV$sample, by=list(delSV$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary2 <- mut_summary2[order(mut_summary2$x, decreasing=T),]

mut_summary3 <- aggregate(dupSV$sample, by=list(dupSV$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary3 <- mut_summary3[order(mut_summary3$x, decreasing=T),]

mut_summary4 <- aggregate(invSV$sample, by=list(invSV$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary4 <- mut_summary4[order(mut_summary4$x, decreasing=T),]

mut_summary5 <- aggregate(traSV$sample, by=list(traSV$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary5 <- mut_summary5[order(mut_summary5$x, decreasing=T),] #####IN QUESTO CE NE è UNO INTERESSANTE CON 23 CAMPIONI#####

write.table(mut_summary, "/mnt/vol1/projects/clingen/driver_mutations/tabs/allSV_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary2, "/mnt/vol1/projects/clingen/driver_mutations/tabs/delSV_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary3, "/mnt/vol1/projects/clingen/driver_mutations/tabs/dupSV_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary4, "/mnt/vol1/projects/clingen/driver_mutations/tabs/invSV_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary5, "/mnt/vol1/projects/clingen/driver_mutations/tabs/traSV_summary.txt", quote=F, sep="\t", row.names=F)

#Per le traslocazioni facciamo una cosa a parte

svs<-dir("/mnt/vol1/projects/clingen/Somatic/Somatic_SV", pattern=".gff$", full.names=T)
svs<-svs[svs != "/mnt/vol1/projects/clingen/Somatic/Somatic_SV/FVG3T.delly.somatic.sv.gff"]

for (aaa in 1:length(svs))
{
myfile<-fread(cmd=paste("grep -v 'SVType=breakpoint'", svs[aaa]),data.table=F, header=F, sep="\t")
myfile<-myfile[grep("breakpoint", myfile$V9, invert = T),]
myfile<-myfile[myfile$V3 == "TRA",]
myfile$tchr <- unlist(lapply(strsplit(myfile$V9, "TCHR="), "[", 2))
myfile$tchr <- unlist(lapply(strsplit(myfile$tchr, ";"), "[", 1))
myfile$tstart <- unlist(lapply(strsplit(myfile$V9, "TSTART="), "[", 2))
myfile$tstart <- unlist(lapply(strsplit(myfile$tstart, ";"), "[", 1))
myfile<-myfile[,c(1:5, ncol(myfile)-1, ncol(myfile))]
myfile$genicstart<-0
myfile$genicend<-0
mychr<-unique(myfile$V1)
mychrend<-unique(myfile$tchr)
endchrom<-min(24,length(mychr))
mychror<-mychr[1:endchrom]
mychr<-paste("chr",mychror,sep="")
mychrend<-paste("chr",mychror, sep="")
myfile$tstart<-as.numeric(myfile$tstart)
myfile$startname<-NA
myfile$endname<-NA
for (bbb in 1:length(mychr))
{
cat("bbb =",bbb,"\n")
stransgff<-transgff[transgff$V1==mychr[bbb],]
myfilechr<-myfile[myfile$V1==mychror[bbb],]
for (ccc in 1:nrow(myfilechr))
{
cat("ccc =",ccc,"\n")
tdf<-cbind(rep(myfilechr$V4[ccc], nrow(stransgff)),rep(myfilechr$V5[ccc], nrow(stransgff)),stransgff$V4, stransgff$V5)
overlap<-apply(tdf,1,count.overlap.vector)
if (length(overlap)>0 & sum(overlap)>0)
{
tgene <-stransgff[which.max(overlap),]
cat(nrow(tgene),"\n")##
myfilechr$startname[ccc]<-tgene$V9
myfilechr$genicstart[ccc]<-1
}
if (myfilechr$tchr[ccc] != "na")
{
ttransgff<-transgff[transgff$V1==paste("chr",myfilechr$tchr[ccc],sep=""),]
if (nrow(ttransgff) > 0)
{
tdf<-cbind(rep(myfilechr$tstart[ccc], nrow(ttransgff)),rep(myfilechr$tstart[ccc], nrow(ttransgff)),ttransgff$V4, ttransgff$V5)
overlap<-apply(tdf,1,count.overlap.vector)
if (length(overlap)>0 & sum(overlap)>0)
{
tgene <-ttransgff[which.max(overlap),]
cat(nrow(tgene),"\n")##
myfilechr$endname[ccc]<-tgene$V9
myfilechr$genicend[ccc]<-1
}
}
}
}
ifelse (bbb==1, itra<-myfilechr[myfilechr$genicstart==1 | myfilechr$genicend==1,], itra<-rbind(itra, myfilechr[myfilechr$genicstart==1 | myfilechr$genicend==1,]))
}
itra$startname<-unlist(lapply(strsplit(itra$startname, "\""), "[", 2))
itra$endname<-unlist(lapply(strsplit(itra$endname, "\""), "[", 2))
write.table(itra, paste("/mnt/vol1/projects/clingen/driver_mutations/translocations/", gsub(".delly.somatic.sv.gff",".txt",basename(svs[aaa])), sep = ""), quote=F, row.names=F, sep="\t")
}

for (file in dir("/mnt/vol1/projects/clingen/driver_mutations/translocations", full.names=T))
{
itra<-fread(file, data.table=F)
mygenes<-unique(na.omit(c(itra$startname,itra$endname)))
genesample<-data.frame(gene=mygenes,sample=gsub(".txt","",basename(file)), stringsAsFactors=F)
nostrogene<-itra[itra$startname == "ENSG00000230876" | itra$endname == "ENSG00000230876",]
if (file == "/mnt/vol1/projects/clingen/driver_mutations/translocations/FVG10T.txt")
{
fgenesample<-genesample
fnostrogene<-nostrogene
next
}
else 
{
fgenesample<-rbind(genesample,fgenesample)
fnostrogene<-rbind(nostrogene,fnostrogene)
}
}

mut_summary <- aggregate(fgenesample$sample, by=list(fgenesample$gene), FUN=function(x) (ngenes=length(unique(x))))
mut_summary <- mut_summary[order(mut_summary$x, decreasing=T),] #ENSG00000230876 

mut_summary2 <- aggregate(fgenesample$sample, by=list(fgenesample$gene), FUN=function(x) (names=paste(x,sep=";")))

mut_summary3 <- merge(mut_summary, mut_summary2, by="Group.1", sort = F)

mut_summary3$x.y<-unlist(lapply(mut_summary3$x.y,paste,sep=";",collapse=";"))

colnames(mut_summary3)[2]<-"x"

write.table(mut_summary3, "/mnt/vol1/projects/clingen/driver_mutations/tabs/traSV_summary.txt", quote=F, sep="\t", row.names=F)

fgenesample

#DOBBIAMO FARE UN LOOP IN CUI PER CIASCUN CAMPIONE SI PRENDE QUALI GENI SONO TOCCATI DALLA TRASLOCAZIONE E SI CONTA OGNI GENE IN QUANTI CAMPIONI VIENE TROVATO