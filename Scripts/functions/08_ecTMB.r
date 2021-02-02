# TMB was defined as the number of somatic, coding,
# base substitution, and indel mutations per megabase of
# genome examined. All base substitutions and indels in
# the coding region of targeted genes, including synonymous
# alterations, are initially counted before filtering as
# described below. Synonymous mutations are counted in
# order to reduce sampling noise. While synonymous mutations
# are not likely to be directly involved in creating
# immunogenicity, their presence is a signal of mutational
# processes that will also have resulted in nonsynonymous
# Chalmers et al. Genome Medicine (2017) 9:34 Page 3 of 14
# mutations and neoantigens elsewhere in the genome.
# Non-coding alterations were not counted. Alterations
# listed as known somatic alterations in COSMIC and
# truncations in tumor suppressor genes were not
# counted, since our assay genes are biased toward genes
# with functional mutations in cancer [63]. Alterations
# predicted to be germline by the somatic-germlinezygosity
# algorithm were not counted [64]. Alterations
# that were recurrently predicted to be germline in our
# cohort of clinical specimens were not counted. Known
# germline alterations in dbSNP were not counted. Germline
# alterations occurring with two or more counts in
# the ExAC database were not counted [65]. To calculate
# the TMB per megabase, the total number of mutations
# counted is divided by the size of the coding region of the
# targeted territory. The nonparametric Mann–Whitney Utest
# was subsequently used to test for significance in
# difference of means between two populations.

# DOI 10.1186/s13073-017-0424-2

#Noi per il TMB abbiamo usato SOLTANTO LE MISSENSO

library(data.table)
library(openxlsx)
#carichiamo il FAI della reference per calcolare la porzione codificante del genoma
fai<-fread("/mnt/vol1/projects/clingen/Reference/GCF_000001405.25_GRCh37.p13_cds_from_genomic.fna.fai", data.table=F)
#cdslength<-sum(fai$V2) #NON USARE, IL PAPER DA UNA STIMA PIU AFFIDABILE
cdslength<-25840698 #secondo il paper https://doi.org/10.1186/s13104-019-4343-8 (UTILIZZIAMO QUESTO)
cdslength<-cdslength/1000000 #risultato in Mbasi
#prepariamo la tabella finale
indelfile=dir("/mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/Annotation", pattern="xls", full.names=T)
snvfile=dir("/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation", pattern="xls", full.names=T)
samplename=unlist(lapply(strsplit(basename(snvfile), "\\."), "[", 1)) #qua vuole i \\., boh
#creiamo un df vuoto
results<-data.frame("sample"=samplename)
#loop
vcfsnp=dir("/mnt/vol1/projects/clingen/Somatic/Somatic_SNV", pattern=glob2rx("*vcf.gz"), full.names=T)
vcfindel=dir("/mnt/vol1/projects/clingen/Somatic/Somatic_INDEL", pattern=glob2rx("*vcf.gz"), full.names=T)
gffcnv=dir("/mnt/vol1/projects/clingen/Somatic/Somatic_CNV", pattern="gff", full.names=T)
gffsv=dir("/mnt/vol1/projects/clingen/Somatic/Somatic_SV", pattern="gff", full.names=T)
for (indelno in 1:length(indelfile))
{
if (indelfile[indelno] == "/mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/Annotation/FVG3T.Strelka.somatic.indel.annovar.hg19_multianno.xls") next
indeltoread<-indelfile[indelno]
snvtoread<-snvfile[indelno]
indel<-fread(indeltoread, data.table=F)
snv<-fread(snvtoread, data.table=F)
vars<-rbind(indel, snv)
exvars<-vars[vars$ExonicFunc!=".",] #perchè non serve il \\.??? boh ok
mut_type<-table(exvars$ExonicFunc)
results[indelno,"synonymous SNV"]<-mut_type["synonymous SNV"]
results[indelno,"missense SNV"]<-mut_type["missense SNV"]
results[indelno,"stoploss"]<-mut_type["stoploss"]
results[indelno,"stopgain"]<-mut_type["stopgain"]
results[indelno,"intergenic"]<-sum(vars$Func=="intergenic")
results[indelno,"intronic"]<-sum(vars$Func=="intronic")
results[indelno,"TMB"]<-results[indelno,"missense SNV"]/cdslength
results[indelno,"TMB_syn"]<-(results[indelno,"missense SNV"] + results[indelno,"synonymous SNV"])/cdslength
scanvcf<-scan(vcfsnp[indelno], what="", sep="\n") #lettura vcf: il risultato indica quanti SNP ci sono nel campione
print(vcfsnp[indelno])
results[indelno,"totSNP"]<-sum(substr(scanvcf,1,1)!="#") #elimina le righe di header che iniziano col cancelletto
scanindel<-scan(vcfindel[indelno], what="", sep="\n") 
results[indelno,"totINDEL"]<-sum(substr(scanindel,1,1)!="#")
scancnv<-fread(gffcnv[indelno],data.table=F)
results[indelno,"totCNVgain"]<-sum(scancnv$V3=="gain")
results[indelno,"totCNVloss"]<-sum(scancnv$V3=="loss")
scansv<-fread(gffsv[indelno],data.table=F,sep="\t")
results[indelno,"totSVdel"]<-sum(scansv[,3]=="DEL")
results[indelno,"totSVdup"]<-sum(scansv[,3]=="DUP")
results[indelno,"totSVtra"]<-sum(scansv[,3]=="TRA")
results[indelno,"totSVinv"]<-sum(scansv[,3]=="INV")
} 
#adesso si incorporano i dati relativi alla percentuale di cellule sane e si valuta quanto questo valore influisca sul TMB (ovviamente con percentuali basse ci si aspetta TMB più bassi)
#l'intuizione è abbastanza ovvia, questa è più che altro una dimostrazione
#lo scopo ultimo è capire cosa è successo con FVG22 in cui non si ritrova il valore di BRCAness atteso
isto<-read.xlsx("/mnt/vol1/projects/clingen/TMB/istologia tumori.xlsx")
isto$ID<-paste(isto$ID,"T", sep="")
merged<-merge(results,isto,by.x="sample", by.y="ID")
cor.test(merged$TMB,merged$"%.cellule.tumorali",method="spearman")
merged$class<-0
merged$class[merged$"%.cellule.tumorali">0.5]<-1

wilcox.test(merged$TMB~merged$class) #la divisione in gruppi risulta significativa quindi l'ipotesi è corretta, quando %cell < 50% si ha un TMB più basso
#ora mettiamo in correlazione il TMB con il grado e con la percentuale di cell. assieme (si cerca di predire TMB da grado e %cell)
#il grado lo dividiamo in classi perchè c'è un solo G1 e creerebbe problemi
merged$grade.class<-"G<3"
merged$grade.class[merged$G=="G3"]<-"G3"
modello<-lm(merged$TMB~merged$class + merged$grade.class)
summary(modello) #serve per guardare il risultato e valutare la significatività sulla base del p-value
#si vede che è significativa solo la percentuale di cellule tumorali, non il grado
#aggiungiamo il valore di T, cioè l'estensione del tumore primitivo, e N, cioè lo stato dei linfonodi
merged$t.class<-"T>1"
merged$t.class[grep("T1",merged$T)]<-"T1"
merged$n.class<-"N>0"
merged$n.class[grep("N0",merged$N)]<-"N0"
modello2<-lm(merged$TMB~merged$class + merged$grade.class + merged$t.class + merged$n.class)
#t value � la statistica usata per calcolare il p value (t di student?)
#Estimate: stima del contributo di quella variabile all'outcome, pi� grosso � estimate/std e pi� alto sar� il tvalue e quindi pi� significativo il valore
#Intercept si pi� togliere e anche il tvalue
#continua ad essere significativo solo class, che è la percentuale di cellule tumorali
modello3<-lm(merged$TMB~merged$class + merged$grade.class + merged$t.class + merged$n.class + merged$ER + merged$PR + merged$HER2)
#a furia di aggiungere roba non significativa si perde anche la significatività di class, quindi diamo per buono solo lui e basta

#write.table(merged, "/mnt/vol1/projects/clingen/TMB/snv_indel/snv_indel.txt", sep="\t", quote=F, row.names=F) #Tabella descrizione campioni per tesi
wilcox.test(merged$TMB~merged$grade.class)
wilcox.test(merged$TMB~merged$t.class)
wilcox.test(merged$TMB~merged$n.class)

pdf("/mnt/vol1/projects/clingen/TMB/tmbplots.pdf")
par(mfrow=c(2,2))

#BOXPLOTS DEL WILCOXON
boxplot(merged$TMB~merged$class, col=rainbow(2), xlab="Percentuale tessuto tumorale", ylab="TMB", names=c(expression(""<=50), ">50%"))
# Add data points
mylevels <- levels(as.factor(merged$class))
levelProportions <- summary(merged$class)/nrow(merged)
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- merged[merged$class==thislevel, "TMB"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}

boxplot(merged$TMB~merged$grade.class, col=heat.colors(2), xlab="Grado del tumore", ylab="TMB")
mylevels <- levels(as.factor(merged$grade.class))
#levelProportions <- summary(merged$grade.$class)/nrow(merged)
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- merged[merged$grade.class==thislevel, "TMB"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=0)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}

boxplot(merged$TMB~merged$t.class, col=terrain.colors(3), xlab="Estensione del tumore", ylab="TMB")
mylevels <- levels(as.factor(merged$t.class))
#levelProportions <- summary(merged$t.class)/nrow(merged)
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- merged[merged$t.class==thislevel, "TMB"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=0)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}

boxplot(merged$TMB~merged$n.class, col=cm.colors(2), xlab="Stato dei linfonodi", ylab="TMB")
mylevels <- levels(as.factor(merged$n.class))
#levelProportions <- summary(merged$n.class)/nrow(merged)
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- merged[merged$n.class==thislevel, "TMB"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=0)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
dev.off()

#PLOT CON CORRELAZIONE TMB E TIPI DI MUTAZIONE
tab<-fread("/mnt/vol1/projects/clingen/TMB/tmb+perc.csv", data.table=F)
tab$TMB<-gsub(",",".",tab$TMB)
tab$TMB<-as.numeric(tab$TMB)
tab$totCNV<- tab$totCNVloss + tab$totCNVgain
tab$totSV<- tab$totSVdel + tab$totSVdup + tab$totSVinv + tab$totSVtra

pdf("/mnt/vol1/projects/clingen/TMB/correlazione.pdf")
par(mfrow=c(2,2))
mod<-lm(tab$TMB ~ tab$totSNP)
pvalue<-format(summary(mod)$coefficients[2,4],digits=3)
rsquared<-round(summary(mod)$adj.r.squared,3)
plot(tab$TMB ~ tab$totSNP, xlab = "Numero di SNP", ylab = "TMB", pch = 19)
abline(mod)
text(13000, 1.2, label=paste("p = ",pvalue, sep=""))
text(13000, 0.8, label=bquote(R^2 == .(rsquared)))

mod2<-lm(tab$TMB ~ tab$totINDEL)
pvalue2<-format(summary(mod2)$coefficients[2,4],digits=3)
rsquared2<-round(summary(mod2)$adj.r.squared,3)
plot(tab$TMB ~ tab$totINDEL, xlab = "Numero di INDEL", ylab = "TMB", pch = 19)
abline(mod2)
text(300, 1.2, label=paste("p = ",pvalue2, sep=""))
text(300, 0.8, label=bquote(R^2 == .(rsquared2)))

mod3<-lm(tab$TMB ~ tab$totCNV)
pvalue3<-format(summary(mod3)$coefficients[2,4],digits=3)
rsquared3<-round(summary(mod3)$adj.r.squared,3)
plot(tab$TMB ~ tab$totCNV, xlab = "Numero di CNV", ylab = "TMB", pch = 19)
abline(mod3)
text(1200, 1.2, label=paste("p = ",pvalue3, sep=""))
text(1200, 0.8, label=bquote(R^2 == .(rsquared3)))

mod4<-lm(tab$TMB ~ tab$totSV)
pvalue4<-format(summary(mod4)$coefficients[2,4],digits=3)
rsquared4<-round(summary(mod4)$adj.r.squared,3)
plot(tab$TMB ~ tab$totSV, xlab = "Numero di SV", ylab = "TMB", pch = 19)
abline(mod4)
text(700, 1.2, label=paste("p = ",pvalue4, sep=""))
text(700, 0.8, label=bquote(R^2 == .(rsquared4)))
dev.off()
