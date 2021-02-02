######################################
# ATTENTION!!!!
# I am currently using a FAKE callable loci file, just for testing purposes.
# I will need to rerun the analysis using a real list of callable loci!
#CallableLoci restituisce tutte le porzioni del genoma in cui si hanno abbastanza info per chiamare uno SNP
######################################

#comandi dallo script03 che servono
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV",
              help="Kmers file", metavar="character"),
	 make_option(c("-V", "--vcf"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/FVG2T.muTect.somatic.snv.vcf",
              help="VCF file (input)", metavar="character"),
	 make_option(c("-B", "--bed"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/CallableLoci_output/FVG2T.final.bed",
              help="BED file (input)", metavar="character"),
	make_option(c("-P", "--outpng"), type="character", default="/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/CallableLoci_output/FVG2T.final.png",
              help="PNG file (output)", metavar="character") #bisogna togliere FVG11 e impostare il loop
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
  }
  
if (is.null(opt$vcf)) {
  stop("WARNING: No input VCF file specified with '-V' flag.")
} else {  cat ("Input VCF is ", opt$vcf, "\n")
  vcf <- opt$vcf
  }
  
if (is.null(opt$bed)) {
  stop("WARNING: No input BED file specified with '-B' flag.")
} else {  cat ("Input BED is ", opt$bed, "\n")
  bed <- opt$bed
  }
  
  if (is.null(opt$outpng)) {
  stop("WARNING: No output PNG file specified with '-P' flag.")
} else {  cat ("Output PNG is ", opt$outpng, "\n")
  outpng <- opt$outpng
  }
  
enrichment_analysis<-function(indir, vcf, bed, outpng)
{ 
#preliminari presi da script precedenti
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
for (myfile in dir(indir, pattern="vcf", full.names=TRUE))
{
vcf_files <- myfile
sample_names<-unlist(lapply(strsplit(vcf_files,"\\."),"[",1))#per separare il punto serve il doppio backslash, la quadra serve per lapply
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)#legge tutti i vcf come g-ranges, assegnandogli un nome e un genoma di riferimento
tissue<-rep("breast",length(sample_names)) #crea un vettore con scritto "breast" tante volte quanti sono i campioni

#Loop per tutti i vcf tranne fvg3t
base_vcf <- basename(vcf_files)
if (base_vcf=="FVG3T.muTect.somatic.snv.vcf") next
base_bed <- gsub("muTect.somatic.snv.vcf", "final.bed", base_vcf)
bed <- paste(indir, "CallableLoci_output", base_bed, sep = "/")
surveyed_file <- bed

print("Importing library biomaRt...")
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
outpng <- gsub("bed", "png", surveyed_file)
cat("Now plotting...", outpng, "\n")
png(outpng, width=7,height=7,units="in",res=600,type="cairo")
plot_enrichment_depletion(distr_test)
myplot <- plot_enrichment_depletion (distr_test)
print(myplot) #Quando ggplot è chiamato all'interno di una funzione non stampa a meno che il plot non venga assegnato ad una variabile
dev.off()
}
setwd(starting_dir)
}
enrichment_analysis(indir=indir, vcf=vcf, bed=bed, outpng=outpng)

#per convertire oggetto GRanges in data frame
#df <- data.frame(chromosome=seqnames(promoter_g),
 # starts=start(promoter_g)-1,
 # ends=end(promoter_g))
