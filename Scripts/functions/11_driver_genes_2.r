library(data.table)
out<-dir("/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/output", pattern=".annotated", full.names=T)
out<-out[out!="/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/output/FVG3T.annotated"] #togliamo FVG3
#loop per impilare i file in un file unico, si aggiunge una colonna che indica il campione
for (aaa in 1:length(out))
{
samplename<-gsub(".annotated","",basename(out[aaa]))
read<-fread(out[aaa])
read$sample<-samplename
ifelse (aaa==1, totalread<-read, totalread<-rbind(totalread,read))
}

#Ora si filtra per p-value e per tenere solo i driver genes
sig<-totalread[totalread$CanDrA_Significance<0.05,]
sig<-sig[sig$CanDrA_Category == "Driver",]

#però stampiamo anche un file che contiene tutti i risultati, driver e passenger
#aggregated <- totalread[, .(allsample=paste(sample,sep=";")), by=list(Chrom, Coordinate, Ref_Allele, Mut_Allele)]
totalread$aggregate <-paste(totalread$Chrom, totalread$Coordinate, totalread$Ref_Allele, totalread$Mut_Allele, sep = "_")
aggregated <- aggregate(totalread$sample, by=list(totalread$aggregate),FUN="paste",sep=";")
aggregated$samplenumber <- unlist(lapply(aggregated$x,length))
aggregated$x <- unlist(lapply(aggregated$x,paste,sep=";",collapse=";"))
final <- merge(totalread,aggregated,by.x="aggregate",by.y="Group.1")
final <- final[!duplicated(final$aggregate),]
final <- final[,!"aggregate"]
final <- final[order(as.numeric(final$Chrom), final$Coordinate),]
final$sample <- final$x
final <- final[,!"x"]
final_orig<-final

#Si vede che il massimo di condivisione per mutazioni nella stessa posizione è solo 4 campioni quindi niente di che

#salviamo la tabella
#write.table(sig, "/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/output/significant_drivers.txt", quote=F, row.names=F, sep="\t")
#la tabella seguente permette di vedere i geni, spiccano PRAMEF1 (3), PRAMEF2 (2) e NBPF10 (2)
#write.table(table(sig$HGNC), "/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/output/genes.txt", quote=F, row.names=F, sep="\t")

#Questa tabella comprende una colonna che contiene i nomi dei campioni in cui viene riscontrata la mutazione e NON è FILTRATA
#write.table(final,"/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/output/final.txt", quot=F, row.names=F,sep="\t")

#aggiungiamo le stopgain e le splice site

intogen <- fread("/mnt/vol1/projects/clingen/driver_mutations/intogen/mutation_analysis.tsv", data.table=F)
intogen<-intogen[,c(2:6,16,17)]
annotations<-dir("/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation", full.names=T)
annotations<-annotations[annotations != "/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation/FVG3T.muTect.somatic.snv.annovar.hg19_multianno.xls"]
for (aaa in annotations)
{
file<-fread(aaa)
file<-file[,1:10]
samplename<-gsub(".muTect.somatic.snv.annovar.hg19_multianno.xls","",basename(aaa))
file$sample<-samplename
ifelse (aaa=="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation/FVG10T.muTect.somatic.snv.annovar.hg19_multianno.xls", totalfile<-file, totalfile<-rbind(totalfile,file))
}
totalfile$aggregate <-paste(totalfile$CHROM, totalfile$POS, totalfile$REF, totalfile$ALT, sep = "_")
aggregated <- aggregate(totalfile$sample, by=list(totalfile$aggregate),FUN="paste",sep=";")
aggregated$samplenumber <- unlist(lapply(aggregated$x,length))
aggregated$x <- unlist(lapply(aggregated$x,paste,sep=";",collapse=";"))
final <- merge(totalfile,aggregated,by.x="aggregate",by.y="Group.1")
final <- final[!duplicated(final$aggregate),]
final <- final[,!"aggregate"]
final <- final[order(as.numeric(final$CHROM), final$POS),]
final$sample <- final$x
final <- final[,!"x"]
names(intogen)[c(6,7)]<-c("into_protein","into_consequence")
intogen<-intogen[substr(intogen$into_consequence,1,6) == "Splice",]
splicestop<-merge(final,intogen,by.x=c("CHROM","POS","REF","ALT"),by.y=c("chr","pos","ref","alt"))
splicestop<-splicestop[!duplicated(splicestop),]
#splicing<-final[final$Func == "splicing",]
names(splicing)[10]<-"mRNA"
splicing<-splicing[,c(-6,-7)]
names(splicing)[6]<-"HGNC"
names(splicestop)[c(1:4,8,10,13,14,15)]<-c("Chrom","Coordinate","Ref_Allele","Mut_Allele","HGNC","mRNA","Strand","AAS","CONSEQUENCE")
splicestop<-splicestop[,c(-5,-6,-7,-9)]
final2<-merge(final_orig,splicestop,by=names(splicestop),all=T)

write.table(final2,"/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/output/final2.txt", quot=F, row.names=F,sep="\t")