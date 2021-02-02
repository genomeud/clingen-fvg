library(data.table)
intogen <- fread("/mnt/vol1/projects/clingen/driver_mutations/intogen/mutation_analysis.tsv", data.table=F) #VERSIONE DEL 2016
candra <- fread("/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/sorted_CanDrA_BRCA.gz", data.table=F)
nostri <- fread("/mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsSNP.txt", data.table=F) #l'elenco delle mutazioni annotate nei nostri campioni
intogen <- intogen[,c(2:6,8:ncol(intogen))] #togliamo due colonne inutili
intogen <- intogen[!duplicated(intogen),] #perdiamo circa 1000 righe (da 56k a 55k)
#possiamo buttare via tutte le intergeniche dato che non saranno mai driver
nostri <- nostri[nostri$Func != "intergenic",] #questo dimezza circa le righe
#togliamo FILTER, Gencode, cytoBand, wgRna, targetScanS, e tutto quello che viene dopo
nostri <- nostri[,c(1:6,8:13, ncol(nostri))]
#cambiamo il nome dell'ultima colonna (era buggato)
colnames(nostri)<-c(names(nostri[1:12]), "sample")
tutto1<-merge(nostri,candra,by.x=c("CHROM", "POS", "REF", "ALT"),by.y=c("V1","V2","V3","V4"), all.x=TRUE) #merge tra nostri e candra con solo cromosoma/posizione/ref/alt, intogen lo aggiunge il comando successivo. all.x serve per non fare l'intersezione tra i due file
tutto2<-merge(tutto1,intogen,by.x=c("CHROM", "POS", "REF", "ALT"),by.y=c("chr","pos","ref","alt"), all.x=TRUE) 
#buttiamo via le mutazioni che non sono nè in candra nè in intogen, ma prima creiamo una copia parallela che utilizzeremo dopo per fare le statistiche
tutto3<-tutto2
tutto2<-tutto2[is.na(tutto2$V12)==FALSE | is.na(tutto2$gene)==FALSE | is.na(tutto2$V11)==FALSE | is.na(tutto2$cancer)==FALSE,]
temp<-tutto2[tutto2$gene == "PIK3CA",]
temp[temp$gene %in% "PIK3CA",] #la mutazione E542K di PIK3CA, che le colleghe hanno segnalato come assente, in realtà l'abbiamo trovata ma a detta di Candra il p-value non è significativo
temp2<-tutto2[tutto2$gene == "TP53",]
temp2[temp2$gene %in% "TP53",13:20]
#il database che ci siamo creati sembra abbastanza affidabile, mancano da candra le nonsense e le splicing variants
#Buttare: QUAL,Gene,GeneDetail,V6,AAChange,strand_orig,info,region,strand,protein,known_predisposing exac_af,missense_driver_mut_prediction, disrupting_driver_mut_prediction,inframe_driver_mut_prediction,known_match
tutto2<-tutto2[,!colnames(tutto2) %in% c("QUAL","Gene","GeneDetail","V6","AAChange","strand_orig","info","region","strand","protein","known_predisposing exac_af","missense_driver_mut_prediction", "disrupting_driver_mut_prediction","inframe_driver_mut_prediction","known_match")]
#Con TP53 vediamo che Intogen, in più a Candra, contiene una nonsenso e uno SpliceDonor
sum(is.na(tutto2$V10) | tutto2$V10 %in% "N/A")  #da qui vediamo che delle nostre mutazioni presenti in almeno uno dei due database, solo 90 non vengono annotate da candra (alcune perchè sono stopgain, altre perchè non ci sono proprio. è indicativo solo di quanto sia completo candra rispetto a intogen, non in assoluto
sum(is.na(tutto2$transcript_exons)) #ci sono 4461 mutazioni che non sono presenti in intogen ma in candra si (sempre in riferimento alle nostre). 
sum(!is.na(tutto2$V10) & !is.na(tutto2$transcript_exons)) #soltanto 20 mutazioni del nostro dataset sono in intogen, ma c'è qualcosa di strano, si dovrebbe rifare
browser()
#ora cambio i nomi alle colonne per preparare il file per il cro:
names(tutto2)[names(tutto2) == "V5"] <- "GeneName.1"
names(tutto2)[names(tutto2) == "V7"] <- "MutationType"
names(tutto2)[names(tutto2) == "V8"] <- "AAChange"
names(tutto2)[names(tutto2) == "V9"] <- "AAPosition"
names(tutto2)[names(tutto2) == "V11"] <- "DriveSignificance"
names(tutto2)[names(tutto2) == "V12"] <- "PValue"
tutto2 <- tutto2[,c(-9,-14)] #TOGLIE V10 E SAMPLE
tutto2 <- tutto2[,1:23] #TOGLIE LE COLONNE IN ECCESSO
tutto2 <- tutto2[!duplicated(tutto2),]
write.table(tutto2, "/mnt/vol1/projects/clingen/driver_mutations/intersection_nostre_candra_intogen_4CRO.txt", quote=F, sep="\t", row.names=F) #DA DARE AL CRO

#############################
#ANALISI STATISTICA
#############################

pippo<-!(tutto3$GeneName %in% "." & is.na(tutto3$V5) & is.na(tutto3$gene))
tutto3 <- tutto3[pippo,]
mut_summary <- aggregate(tutto3$sample, by=list(tutto3$GeneName), FUN=function(x) (ngenes=length(unique(x))))
mut_summary <- mut_summary[order(mut_summary$x, decreasing=T),]

#####fare una cosa simile per le introniche, le non-introniche, driver###

tutto3int<-tutto3[tutto3$Func == "intronic",]
mut_summary2 <- aggregate(tutto3int$sample, by=list(tutto3int$GeneName), FUN=function(x) (ngenes=length(unique(x))))
mut_summary2 <- mut_summary2[order(mut_summary2$x, decreasing=T),]

tutto3nint<-tutto3[tutto3$Func != "intronic",]
mut_summary3 <- aggregate(tutto3nint$sample, by=list(tutto3nint$GeneName), FUN=function(x) (ngenes=length(unique(x))))
mut_summary3 <- mut_summary3[order(mut_summary3$x, decreasing=T),]

tutto3dri <- tutto3[tutto3$V11 %in% "Driver",] ##il più rilevante##
mut_summary4 <- aggregate(tutto3dri$sample, by=list(tutto3dri$GeneName), FUN=function(x) (ngenes=length(unique(x))))
mut_summary4 <- mut_summary4[order(mut_summary4$x, decreasing=T),]

##SALVIAMO LE TABELLE CON IL NUMERO DI CAMPIONI IN CUI, PER CIASCUN GENE, VI è UNA MUTAZIONE DEL TIPO SCELTO##
write.table(mut_summary, "/mnt/vol1/projects/clingen/driver_mutations/tabs/allSNPs_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary2, "/mnt/vol1/projects/clingen/driver_mutations/tabs/intSNPs_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary3, "/mnt/vol1/projects/clingen/driver_mutations/tabs/nintSNPs_summary.txt", quote=F, sep="\t", row.names=F)
write.table(mut_summary4, "/mnt/vol1/projects/clingen/driver_mutations/tabs/driSNPs_summary.txt", quote=F, sep="\t", row.names=F)