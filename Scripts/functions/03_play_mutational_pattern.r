# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-C", "--coscos"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Plots/cosine_cosmic_heatmap.png", 
              help="Cosine cosmic heatmap plot [default= %default]", metavar="character"),
  make_option(c("-I", "--indir"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV",
              help="Kmers file", metavar="character"),
  make_option(c("-R", "--relcont"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Plots/relative_contribution_heatmap.png", 
              help="Relative contribution of signatures plot [default= %default]", metavar="character"),
  make_option(c("-T", "--tablefile"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Tables/mutation_distribution_strand_independent.txt", 
              help="Table with mutation distribution [default= %default]", metavar="character")
  # make_option(c("-C", "--control_prefix"), type="character", default="control_",
              # help="Prefix identifying control samples", metavar="character"),
  # make_option(c("-R", "--raw_counts"), type="character", default=NULL,
              # help="raw (total) read counts for this starting file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No input directory file specified with '-I' flag.")
} else {  cat ("Input directory is ", opt$indir, "\n")
  indir <- opt$indir  
  #setwd(wd_location)  
  }

if (is.null(opt$coscos)) {
  stop("WARNING: No cosine cosmic heatmap file specified with '-C' flag.")
} else {  cat ("Cosine cosmic heatmap file is ", opt$coscos, "\n")
  coscos <- opt$coscos  
  #setwd(wd_location)  
  }

  if (is.null(opt$relcont)) {
  stop("WARNING: No relative contribution heatmap file specified with '-R' flag.")
} else {  cat ("relative contribution heatmap file is ", opt$relcont, "\n")
  relcont <- opt$relcont  
  #setwd(wd_location)  
  }

  if (is.null(opt$tablefile)) {
  stop("WARNING: No table file specified with '-T' flag.")
} else {  cat ("Table file is ", opt$tablefile, "\n")
  tablefile <- opt$tablefile  
  #setwd(wd_location)  
  }

#I am using MutationalPatterns 
#http://bioconductor.org/packages/release/bioc/vignettes/MutationalPatterns/inst/doc/Introduction_to_MutationalPatterns.pdf
mutation_distr<-function(indir,coscos,relcont,tablefile)
{
library("data.table")
library("gplots")
library(MutationalPatterns)
library(BSgenome)
#Load the needed genome
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") #ESEGUIRE UNA VOLTA SOLA PER INSTALLARE IL GENOMA
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
starting_dir <- getwd()
setwd(indir)
vcf_files<-dir(pattern=".vcf$")
vcf_files<-vcf_files[vcf_files != "FVG3T.muTect.somatic.snv.vcf"]
sample_names<-unlist(lapply(strsplit(vcf_files,"\\."),"[",1))#per separare il punto serve il doppio backslash, la quadra serve per lapply
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)#legge tutti i vcf come g-ranges, assegnandogli un nome e un genoma di riferimento
tissue<-rep("breast",length(sample_names)) #crea un vettore con scritto "breast" tante volte quanti sono i campioni
muts = mutations_from_vcf(vcfs[[1]])#estrae tutte le mutazioni senza informazioni aggiuntive (A>B e basta)
types = mut_type(vcfs[[1]]) # cataloga le mutazioni secondo i 6 tipi convenzionali, giusto?
context = mut_context(vcfs[[1]], ref_genome) #esplicita il contesto (nucleotidi adiacenti) e il cromosoma
type_context = type_context(vcfs[[1]], ref_genome) #contiene i due elementi (types e context) di cui ciascuno contiene ulteriori elementi, type_context[[1]][1] visualizza il primo sottoelemento del primo elemento
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)#tabella che quantifica i tipi di mutazioni nei vari campioni
p1 <- plot_spectrum(type_occurrences) #barplot
p2 <- plot_spectrum(type_occurrences, CT = TRUE)#come sopra ma distingue le C>T in contesto CpG dalle altre, perchè in CpG è abbastanza comune la deaminazione delle citosine (che quindi diventano U)
p3 <- plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)#come sopra ma senza legenda
library("gridExtra") #installata
#VA A SCHERMO grid.arrange(p1, p2, p3, ncol=3, widths=c(3,3,1.75)) #mette insieme i 3 grafici di prima in un'unica immagine
p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE) #esplicita il tessuto nel titolo del grafico (breast)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)#lento, fa il profilo 96 di ciascun campione (calcola il numero di mutazioni in tutti i contesti).
browser()
#VA A SCHERMO plot_96_profile(mut_mat[,c(1,2)], condensed = TRUE,ymax=0.1) #plotta i contesti per ciascuna delle mutazioni nei primi 2 campioni
#Use NMF to estimate mutational pattern (the approach is the one that Alexandrov et al implemented in MatLab)
mut_mat <- mut_mat + 0.0001 #somma 0,0001 a tutti i valori per evitare che il calcolo successivo ritrovi dei valori "0"
library("NMF")#c'è già
 #estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456) #ne parla nel pdf, estrae 2/3/4/5 signatures ma è soltanto una prova
#Extract two mutational signatures from the datasets. We might want to change the nrun parameter

#RIMUOVERE ASSOLUTAMENTE IL CAMPIONE 3#

 nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 10) #estrae le signatures e ne esplicita il contributo. 
 colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
 rownames(nmf_res$contribution) <- c("Signature A", "Signature B")
 #VA A SCHERMO plot_96_profile(nmf_res$signatures, condensed = TRUE)
 
 #Visualize relative contribution of signatures for each sample
pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "relative")
pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute")
 
#Plot signature contribution as a heatmap with sample clustering dendrogram and a specified signature order:
 pch1 <- plot_contribution_heatmap(nmf_res$contribution, sig_order = c("Signature B", "Signature A"))
#Plot signature contribution as a heatmap without sample clustering
 pch2 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=FALSE)
#VA A SCHERMO grid.arrange(pch1, pch2, ncol = 2, widths = c(2,1.6)) #mette insieme i due grafici di prima

#Compare the reconstructed mutational profile with the original mutational profile:
#VA A SCHERMO plot_compare_profiles(mut_mat[,1], nmf_res$reconstructed[,1], profile_names = c("Original", "Reconstructed"), condensed = TRUE) #reconstructed distribuisce le mutazioni nei campioni secondo le signature, praticamente operazione inversa rispetto al calcolo delle contribution

 
 ##Compare with COSMIC data
 
#Download mutational signatures from the COSMIC website:
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96-contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
#Plot mutational profile of the first two COSMIC signatures:
#VA A SCHERMO plot_96_profile(cancer_signatures[,1:2], condensed = TRUE, ymax = 0.3)
 
#Hierarchically cluster the COSMIC signatures based on their similarity with average linkage:
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
#VA A SCHERMO plot(hclust_cosmic)
 #Calculate pairwise cosine similarity between mutational profiles and COSMIC signatures:
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
png(coscos,width=7,height=7,units="in",res=600,type="cairo")
#pdf(coscos)
#par(mar=c(0.8,1.6,0.6,1.4),oma=c(0.2,0.2,0.5,0.1),tck=-0.01)
 plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order, cluster_rows = TRUE) #descrizione delle signature: https://cancer.sanger.ac.uk/cosmic/signatures_v2
dev.off()

#rimuoviamo temporaneamente il campione 3 (cancellare questo comando quando lo si rimuove sul serio)
#mut_mat=mut_mat[,!colnames(mut_mat)%in%"FVG3T"]
#Fit mutation matrix to the COSMIC mutational signatures:
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
# Select signatures with total contribution between samples > 10
select <- which(rowSums(fit_res$contribution) > 10)
# Plot contribution barplot
png(gsub(".png", "_tnp.png", relcont),width=7,height=7,units="in",res=600,type="cairo") #ELIMINARE SE NON SI FA IN SCREEN
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = TRUE, mode = "absolute")
dev.off()#ELIMINARE SE NON SI FA IN SCREEN

png(relcont,width=7,height=7,units="in",res=600,type="cairo")
plot_contribution_heatmap(fit_res$contribution[select,], cluster_samples = TRUE, method = "complete")
dev.off()

#Compare our mutation matrix with the reconstruction obtained by fitting out mutations to COSMIC signatures.
# calculate all pairwise cosine similarities 
cos_sim_ori_rec <- cos_sim_matrix(mut_mat, fit_res$reconstructed)
# extract cosine similarities per sample between original and reconstructed (alla fine si ottiene, per ogni campione, se stesso originale vs ricostruzione)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
# Adjust data frame for plotting with gpplot
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
# Load ggplot2
library(ggplot2)
# Make barplot
png(gsub(".png", "_tnp2.png", relcont),width=7,height=7,units="in",res=600,type="cairo")#ELIMINARE SE NON SI FA IN SCREEN
ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
 geom_bar(stat="identity", fill = "skyblue4") +
 coord_cartesian(ylim=c(0.8, 1)) +
 coord_flip(ylim=c(0.8,1)) +
 ylab("Cosine similarity\n original VS reconstructed") +
 xlab("") +
 # Reverse order of the samples such that first is up
 # xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +
 theme_bw() +
 theme(panel.grid.minor.y=element_blank(),
 panel.grid.major.y=element_blank()) +
 # Add cut.off line
 geom_hline(aes(yintercept=.95)) #SOGLIA DI ACCETTAZIONE
 dev.off()#ELIMINARE SE NON SI FA IN SCREEN

#Analyse transcriptional strand bias
#Get gene definitions for your reference genome:
# For example get known genes table from UCSC for hg19 using
#The line below is commented because it is needed only once.
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#Get transcriptional strand information for all positions in the first VCF object with mut_strand 
 strand = mut_strand(vcfs[[1]], genes_hg19)
#Make mutation count matrix with transcriptional strand information (96 trinucleotides * 2 strands = 192 features). NB: only those mutations that are located within gene bodies are counted
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
#Count the number of mutations on each strand, per tissue, per mutation type:
strand_counts <- strand_occurrences(mut_mat_s, by=tissue)
#Perform Poisson test for strand asymmetry significance testing:
strand_bias <- strand_bias_test(strand_counts) #vengono tutti significativi il che è strano, vedere cosa succede senza il campione 3
#Plot the mutation spectrum with strand distinction:
png(gsub(".png", "_tnp3.png", relcont),width=7,height=7,units="in",res=600,type="cairo") #ELIMINARE SE NON SI FA IN SCREEN
ps1 <- plot_strand(strand_counts, mode = "relative")
#Plot the effect size (log2(untranscribed/transcribed) of the strand bias. Asteriks indicate significant strand bias.
ps2 <- plot_strand_bias(strand_bias)
grid.arrange(ps1, ps2)
dev.off() #ELIMINARE SE NON SI FA IN SCREEN

#Extract 2 signatures from mutation count matrix with strand features:
######################
#Apparently this process takes a long time, I skipped this and went on to studying the genomic distribution of features.
######################
nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2,nrun=10)
# Provide signature names
colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B")
#Plot signatures with 192 features:
a <- plot_192_profile(nmf_res_strand$signatures, condensed = TRUE)
#Plot strand bias per mutation type for each signature with significance test:
b <- plot_signature_strand_bias(nmf_res_strand$signatures)
#Combine the plots into one figure:
grid.arrange(a, b, ncol = 2, widths = c(5, 1.8))


#Study the genomic distribution of features

#Make rainfall plot of sample 1 over all autosomal chromosomes
# Define autosomal chromosomes
chromosomes <- seqnames(get(ref_genome))[1:22]
# Make a rainfall plot (lo fa solo di FVG10)
for (aaa in 1:length(vcfs))
{
print(gsub(".png", paste("_",names(vcfs[aaa]),".png",sep=""), "/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Plots/rainfall/rainplot.png"))
png(gsub(".png", paste("_",names(vcfs[aaa]),".png",sep=""), "/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Plots/rainfall/rainplot.png"),width=7,height=7,units="in",res=600,type="cairo") #ELIMINARE SE NON SI FA IN SCREEN
tmpplot<-plot_rainfall(vcfs[[aaa]], title = names(vcfs[aaa]), chromosomes = chromosomes, cex = 1.5, ylim = 1e+08)
print(tmpplot)
dev.off() #ELIMINARE SE NON SI FA IN SCREEN
}

browser()

#To be sincere, I don't like these rainfall very much. All the samples look quite the same!

# if (!requireNamespace("BiocManager", quietly=TRUE)) 
# install.packages("BiocManager")
# BiocManager::install("biomaRt")

 library(biomaRt)

regulatory <- useEnsembl(biomart="regulation", dataset="hsapiens_regulatory_feature", GRCh = 37)
# Download the regulatory CTCF binding sites and convert them to
# a GRanges object.
#dalla colonna regulatory feature type name selezione solo quelli che rappresentano siti di legame per CTCF (CTCF definisce il confine tra DNA attivo e eterocromatico (inattivo), cercare in letteratura per la tesi
CTCF <- getBM(attributes = c('chromosome_name','chromosome_start','chromosome_end','feature_type_name'),
    filters = "regulatory_feature_type_name",
    values = "CTCF Binding Site",
    mart = regulatory)
CTCF_g <- reduce(GRanges(CTCF$chromosome_name, IRanges(CTCF$chromosome_start,CTCF$chromosome_end)))

# Download the promoter regions and convert them to a GRanges object.
promoter = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'),
    filters = "regulatory_feature_type_name", values = "Promoter", mart = regulatory)
promoter_g = reduce(GRanges(promoter$chromosome_name, IRanges(promoter$chromosome_start, promoter$chromosome_end)))

# Download the promoter flanking regions and convert them to a GRanges object.
 flanking = getBM(attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 'feature_type_name'),
        filters = "regulatory_feature_type_name", values = "Promoter Flanking Region", mart = regulatory)
 flanking_g = reduce(GRanges( flanking$chromosome_name, IRanges(flanking$chromosome_start, flanking$chromosome_end)))

# Combine all genomic regions (GRanges objects) in a named list:
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
# Use the same chromosome naming convention consistently:
seqlevelsStyle(regions) <- "UCSC"


#Write the reference
genome<-BSgenome.Hsapiens.UCSC.hg19
writeXStringSet(genome,"/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Reference/homo_sapiens_hg19-ucsc_nochr.fasta")
}
mutation_distr(indir=indir,coscos=coscos,relcont=relcont,tablefile=tablefile)
