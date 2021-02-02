library(gdsfmt)
library(SNPRelate)
wd<-getwd()
setwd("/mnt/vol1/projects/clingen/SNPRelate")
vcf=dir("/mnt/vol1/projects/clingen/Mutation/Mutation_SNP", pattern="vcf.gz", full.names=T)
#leggere i VCF
for (file in vcf)
{
vcf.fn<-file
basevcf<-basename(file)
basevcf<-gsub("vcf.gz","gds",basevcf)
snpgdsVCF2GDS(vcf.fn, basevcf, method="biallelic.only")
}

#MERGE dei VCF
direct<-dir("/mnt/vol1/projects/clingen/SNPRelate", pattern="GATK.snp.gds",full.names=T)
#snpgdsCombineGeno(c("FVG10T_test.gds", "FVG10A_test.gds"),"mergetest.gds")
snpgdsCombineGeno(direct,"merged.gds") #se il direct non funziona scrivere i nomi a mano

###########FINO A QUA NON SERVE PIU RIPETERLO#############

##PCA##
genofile<-snpgdsOpen("merged.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2) #come filtro va bene
snpset.id <- unlist(unname(snpset))
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
#calcolo percentuale di variazione spiegata dalle componenti principali
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
#Creo il tab senza sample_annot
tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)
#creo una colonna che metta nella stessa classe i campioni dello stesso paziente per usarla come sample_annot
tab$pop<-gsub("T","",tab$sample.id)
tab$pop<-gsub("A","",tab$pop)

#Provo la strada con le info sulla popolazione utilizzando tab$pop come sample_annot
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- tab$pop
head(cbind(sample.id, pop_code))

#creo i simboli per il grafico
mysymbol<-rep(c(21,22,23,24,25), 6)
mysymbol<-mysymbol[1:26]
mycolor<-sort(rep(1:6, 5))
mycolor<-mycolor[1:26]
mypop<-sort(unique(tab$pop))
mydf <- data.frame(pop=mypop, col=mycolor, symbol=mysymbol)

tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

tab <- merge(tab, mydf, sort=F)

#plottiamo cose brutte perchè R non sa contare
plot(sign(tab$EV2)*abs(tab$EV2)^(1/3), sign(tab$EV1)*abs(tab$EV1)^(1/3), col=tab$col, pch=tab$symbol, xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
##non va ancora benissimo

#IBD MoM tutti vs tutti
ibd <- snpgdsIBDMoM(genofile, snp.id=snpset.id,
    maf=0.05, missing.rate=0.05, num.thread=2)
ibd.coeff <- snpgdsIBDSelection(ibd)

#IBD MLE tutti vs tutti
#snp.id <- sample(snpset.id, 1500)  # random 1500 SNPs #noi li usiamo tutti
ibd <- snpgdsIBDMLE(genofile, snp.id=snpset.id,
    maf=0.05, missing.rate=0.05, num.thread=2)
ibd.coeff <- snpgdsIBDSelection(ibd)

###Vediamo a chi è parente FVG3T###
ibd2<-ibd.coeff[,c("ID1","ID2","kinship")]
reibd<-reshape(ibd2, idvar="ID1", timevar="ID2", direction="wide")
colnames(reibd) <- gsub("kinship.","",colnames(reibd)) #1. togliere kinship dai colnames
rownames(reibd) <- reibd[,1] #2. far diventare la prima colonna i nomi delle righe
reibd <- reibd [,-1] #toglie la colonna 1 che ormai è diventata rownames
#Trasformiamo la matrice da triangolare a quadratica
FVG9T<-rep(NA, 51)
reibd<-rbind(reibd, FVG9T) 
FVG10A<-rep(NA, 52)
reibd<-cbind(FVG10A, reibd) 
rownames(reibd)<-colnames(reibd) #,ette il nome alla riga FVG9T, mancava

#Questa funzione converte la matrice tri->squ
tri.to.squ<-function(x)
{
rn<-row.names(x)
cn<-colnames(x)
an<-unique(c(cn,rn))
myval<-x[!is.na(x)]
mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
for(ext in 1:length(cn))
{
for(int in 1:length(rn))
{
if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
}

}
return(mymat)
}
squibd<-tri.to.squ(reibd)
squibd[squibd==1]=0.5 #l'identità deve essere 0.5, non 1
###Si vede che FVG3T è in relazione di "praticamente" identità con 13A-13T, plausibile che la confusione sia avvenuta tra 3 e 13###
#image(1:52,1:52,squibd,xlab=rownames(squibd),ylab=colnames(squibd)) #brutta ma rende l'idea, da migliorare

#trasformiamo la matrice di similarità in una matrice di dissimilarità
squibd<-0.5-squibd

png("/mnt/vol1/projects/clingen/SNPRelate/heatmap.png",width=7,height=7,units="in",res=600,type="cairo")
heatmap(squibd, Rowv=NA, Colv=NA)
dev.off()

#DEVE ESSERE L'ULTIMA COSA
#setwd(wd)

#snpgdsSummary("FVG10A.GATK.snp.gds")

