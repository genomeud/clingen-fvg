library(data.table)
#LEGGIAMO TUTTO MA SOLO LE PRIME 20 RIGHE
allSNP <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/allSNPs_summary.txt", data.table=F)[1:10,]
intSNP <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/intSNPs_summary.txt", data.table=F)[1:10,]
nintSNP <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/nintSNPs_summary.txt", data.table=F)[1:10,]
driSNP <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/driSNPs_summary.txt", data.table=F)[1:10,]
allINDEL <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/allINDEL_summary.txt",data.table=F)[1:10,] ####NON BELLISSIMO MA CI TOCCA###
exoINDEL <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/exoINDEL_summary.txt",data.table=F)[1:10,] ###NON PLOTTABILE TROPPA POCA ROBA###
allCNV <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/allCNV_summary.txt", data.table=F)[1:10,]
gainCNV <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/gainCNV_summary.txt", data.table=F)[1:10,]
lossCNV <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/lossCNV_summary.txt", data.table=F)[1:10,]
allSV <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/allSV_summary.txt", data.table=F)[1:10,]
delSV <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/delSV_summary.txt", data.table=F)[1:10,]
dupSV <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/dupSV_summary.txt", data.table=F)[1:10,]
invSV <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/invSV_summary.txt", data.table=F)[1:10,]
traSV <- fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/traSV_summary.txt", data.table=F)[1:10,]


pdf("/mnt/vol1/projects/clingen/driver_mutations/tabs/plots/global.pdf")
par(mar=c(6.4,3.2,3.8,2.2),las=2,mfrow=c(3,2), cex=0.6) ##MAR CAMBIA I MARGINI, LAS METTE TUTTO PERPENDICOLARE AGLI ASSI (ANCHE I NOMI DEI GENI)
#barplot(allSNP$x, names=allSNP$Group.1, main="All SNPs")
#barplot(intSNP$x, names=intSNP$Group.1, main="Intronic SNPs")
barplot(nintSNP$x, names=nintSNP$Group.1, main="Non-intronic SNPs")
barplot(driSNP$x, names=driSNP$Group.1, main="Driver SNPs")
barplot(allINDEL$x, names=allINDEL$Group.1, main="All INDELs")
#barplot(exoINDEL$x, names=exoINDEL$Group.1, main="Exonic INDELs")
#barplot(allCNV$x, names=allCNV$Group.1, main="All CNVs")
#barplot(gainCNV$x, names=gainCNV$Group.1, main="Gain CNVs")
barplot(lossCNV$x, names=lossCNV$Group.1, main="Loss CNVs")
#barplot(allSV$x, names=allSV$Group.1, main="All SVs")
# barplot(delSV$x, names=delSV$Group.1, main="Deletion SVs")
# barplot(dupSV$x, names=dupSV$Group.1, main="Duplication SVs")
barplot(invSV$x, names=invSV$Group.1, main="Inversion SVs")
barplot(traSV$x, names.arg=c("LINC00486","LAMA5","PTPRN2","RASA3","CHRNA4","SLC12A7","SHANK2","CDH4","CELSR1","RPH3A"), main="Translocation SVs")
#METTERE names.arg al posto di names e mettere un vettore con i nomi dei geni
dev.off()
#GUARDARE I GENI DI RILIEVO E VEDERE SE QUALCUNO è NOTO, LA SCOMMESSA è CHE NON LO SIANO MA SE QUALCUNO LO è DIVENTA INTERESSANTE
#https://revistaclinicapsicologica.com/data-cms/articles/oldissue/20200918054014am.pdf

#Tabella con il numero di mutazioni a carico di geni noti

geneslist<-fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/genes_list.txt", data.table=F, header=F)
names(geneslist)<-"Geni"
geneslist$"Campioni mutati (SNP)"<-0
allsnp<-fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/allSNPs_summary.txt", data.table=F)
for (aaa in 1:nrow(geneslist))
{
if (sum(allsnp$"Group.1" == geneslist[aaa,1]) == 0) next
geneslist[aaa,2]<-allsnp[allsnp$"Group.1" == geneslist[aaa,1],2]
}

write.table(geneslist, "/mnt/vol1/projects/clingen/driver_mutations/tabs/genes_list_SNP.txt", quote=F, sep="\t", row.names=F)

#Ricerca gain/loss noti

locilist<-fread("/mnt/vol1/projects/clingen/driver_mutations/tabs/loci_list.txt", data.table=F, header=F)
annotations<-dir("/mnt/vol1/projects/clingen/Somatic/Somatic_CNV/Annotation/", full.names=T)
annotations<-annotations[annotations != "/mnt/vol1/projects/clingen/Somatic/Somatic_CNV/Annotation/FVG3T.final.bam_CNVs.final.hg19_multianno.xls"]
finalmat<-matrix(NA,ncol=nrow(locilist),nrow=length(annotations),dimnames=list(basename(annotations),locilist[,1]))
rownames(finalmat)<-gsub(".final.bam_CNVs.final.hg19_multianno.xls","",rownames(finalmat))
for (aaa in 1:length(annotations))
{
annot<-fread(annotations[aaa],data.table=F)
annot<-annot[,c("cytoBand", "CNVType")]
for(bbb in 1:nrow(locilist))
{
region<-grep(paste0("^",locilist[bbb,1]),annot[,1])
regiontype<-sum(annot[region,2] == locilist[bbb,2])
finalmat[aaa,bbb]<-regiontype
}
}
finalmat<-finalmat[,1:5]

write.table(finalmat, "/mnt/vol1/projects/clingen/driver_mutations/tabs/gains_losses.txt", quote=F, sep="\t", col.names=NA)





