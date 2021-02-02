library(biomaRt)
#List available marts. Be careful, for plants you have to set a different website, that I have written somewhere!
#If you get an error 500 for the function below, you should exit R and re-enter because you probably have some 
#mess in the cache. You could also clean Marts' cache, but I don't know exaclty how to do that.
listMarts()

regulatory <- useEnsembl(biomart="regulation",dataset="hsapiens_regulatory_feature",GRCh=37)
promoter<- getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'),

tt<-useEnsembl(biomart="ENSEMBL_MART_SNP",GRCh=37)
listDatasets(mart=tt)
tt<-useEnsembl(biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp_som",GRCh=37)
listAttributes(mart=tt)
#Here I get info on exon, but we can get info on several additional features (e.g. predicted effect of variants and several others)
#We can feed this to the enrichment test of mutationalPattenrs and see what happens, just for fun!!!

tutto1<-getBM(attributes=c("chr_name","chrom_start","chrom_end","clinical_significance"),mart=tt)
tutto2<-getBM(attributes=c("chr_name","chrom_start","chrom_end","phenotype_description"),mart=tt)
tutto3<-getBM(attributes=c("chr_name","chrom_start","chrom_end","associated_gene"),mart=tt)
tutto4<-getBM(attributes=c("chr_name","chrom_start","chrom_end","consequence_type_tv"),mart=tt)
tutto5<-getBM(attributes=c("chr_name","chrom_start","chrom_end","polyphen_prediction"), filters = "polyphen_prediction", values = "probably damaging", mart=tt)
tutto6<-getBM(attributes=c("chr_name","chrom_start","chrom_end","polyphen_prediction"), filters = "polyphen_prediction", values = "possibly damaging", mart=tt)
tutto7<-getBM(attributes=c("chr_name","chrom_start","chrom_end","polyphen_prediction"), filters = "polyphen_prediction", values = "benign", mart=tt)
tutto8<-getBM(attributes=c("chr_name","chrom_start","chrom_end","polyphen_prediction"), filters = "polyphen_prediction", values = "unknown", mart=tt)
#listFilterValues(mart=tt,filter ="polyphen_prediction"), mart=tt)
tutto6<-getBM(attributes=c("chr_name","chrom_start","chrom_end","sift_prediction"),mart=tt)

#trasforma il gc content in un dato qualitativo (tutto1)
tutto1$gc_content<-"high"
tutto1$gc_content[tutto1$percentage_gene_gc_content<40]<-"low"
tutto1<-tutto1[,c(1,2,3,5)]

#seleziona solo le righe in cui il gc content è "low"
tutto1<-tutto1[tutto1$gc_content=="low",]
head(tutto1)

#seleziona solo la famiglia con più membri (tutto2)
fam <- table(tutto2$family)
max(fam[2:length(fam)]) #da qui si vede quanti membri ha la famiglia maggiore
names(fam[fam==148]) #qui si vede come si chiama la famiglia con più membri
tutto2<-tutto2[tutto2$family=="PTHR11738",]

#prende solo la description più rappresentata (tutto3)
fam <- table(tutto3$family_description)
max(fam[2:length(fam)]) #da qui si vede quanto è abbondante la description più abbondante
names(fam[fam==700]) #qui si vede qual è la descrizione più abbondante (è "uncharacterized" ma questo è solo un esercizio quindi chissene
tutto3<-tutto3[tutto3$family_description=="UNCHARACTERIZED",]

#identifica genericamente i trascritti (trasforma l'ID in "transcript")
tutto4$transcript<-"transcript"
tutto4<-tutto4[,c(1,2,3,5)]

#crea i GRanges
tutto1_g <- reduce(GRanges(tutto1$chromosome_name, IRanges(tutto1$start_position,tutto1$end_position)))
tutto2_g <- reduce(GRanges(tutto2$chromosome_name, IRanges(tutto2$start_position,tutto2$end_position)))
tutto3_g <- reduce(GRanges(tutto3$chromosome_name, IRanges(tutto3$start_position,tutto3$end_position)))
tutto4_g <- reduce(GRanges(tutto4$chromosome_name, IRanges(tutto4$start_position,tutto4$end_position)))
tutto5_g <- reduce(GRanges(tutto5$chr_name, IRanges(tutto5$chrom_start,tutto5$chrom_end)))
tutto6_g <- reduce(GRanges(tutto6$chr_name, IRanges(tutto6$chrom_start,tutto6$chrom_end)))
tutto7_g <- reduce(GRanges(tutto7$chr_name, IRanges(tutto7$chrom_start,tutto7$chrom_end)))
#li mette insieme
#regions <- GRangesList(tutto1_g, tutto2_g, tutto3_g, tutto4_g)
regions <- GRangesList(tutto5_g, tutto6_g, tutto7_g)
#names(regions) <- c("GC Content", "Family", "Family Description","Transcript")
names(regions) <- c("Probably Damaging", "Possibly Damaging", "Benign")
#prova con FVG10T
library("data.table")
library("gplots")
library(MutationalPatterns)
library(BSgenome)
#Load the needed genome
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") #ESEGUIRE UNA VOLTA SOLA PER INSTALLARE IL GENOMA
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
vcf_files <- "/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/FVG10T.muTect.somatic.snv.vcf"
sample_names<-unlist(lapply(strsplit(vcf_files,"\\."),"[",1))#per separare il punto serve il doppio backslash, la quadra serve per lapply
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
surveyed_file <- "/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/CallableLoci_output/FVG10T.final.bed"
seqlevelsStyle(regions) <- "UCSC"
tissue<-rep("breast",length(sample_names)) #crea un vettore con scritto "breast" tante volte quanti sono i campioni

# Import the file using rtracklayer and use the UCSC naming standard
print("Importing library rtracklayer...")
library(rtracklayer)
surveyed <- import(surveyed_file)
seqlevelsStyle(surveyed) <- "UCSC"

# For this example we use the same surveyed file for each sample.
surveyed_list <- list(surveyed)

# Calculate the number of observed and expected number of mutations in
# each genomic regions for each sample.
distr <- genomic_distribution(vcfs, surveyed_list, regions)
 
# Perform the enrichment/depletion test by tissue type.
distr_test <- enrichment_depletion_test(distr, by = tissue)
#head(distr_test) #la discrepanza tra il numero di mutazioni originali e quello riportato qui indica che probabilmente quelle scartate non erano nei CallableLoci, e quindi non sufficientemente affidabili.
outpng <- gsub(".bed", "_prova.png", surveyed_file)
cat("Now plotting...", outpng, "\n")
png(outpng, width=7,height=7,units="in",res=600,type="cairo")
plot_enrichment_depletion(distr_test)
myplot <- plot_enrichment_depletion (distr_test)
print(myplot) #Quando ggplot è chiamato all'interno di una funzione non stampa a meno che il plot non venga assegnato ad una variabile
dev.off()




