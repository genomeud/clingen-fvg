#07_HRDetect non verrà usato più, cancellarlo quando si riesce a far funzionare il comando)
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
make_option(c("-S", "--samplename"), type="character", default="FVG4T",
              help="Input sample name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$samplename)) {
  stop("WARNING: No input sample name specified with '-S' flag.")
} else {  cat ("Input sample name is ", opt$samplename, "\n")
  samplename <- opt$samplename  
  #setwd(wd_location)  
  }

library(signature.tools.lib)
myoutput<-HRDetect_pipeline(data_matrix=NULL,
                              genome.v="hg19",
                              SNV_vcf_files=c(test=paste("/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/", samplename, ".muTect.somatic.snv.vcf.gz", sep="")),
                              SNV_tab_files=NULL,
                              SNV_catalogues=NULL,
                              Indels_vcf_files=c(test=paste("/mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/",samplename,".Strelka.somatic.indel.vcf.gz", sep="")),
                              CNV_tab_files=c(test=paste("/mnt/vol1/projects/clingen/Control-FREEC/cnv_for_HRDetect/",samplename,"_for_HRDetect.txt",sep="")),
                              SV_bedpe_files=c(test=paste("/mnt/vol1/projects/clingen/Somatic/Somatic_SV/HRDetect/",samplename,".delly.somatic.sv.bed",sep="")),
                              SV_catalogues=NULL,
                              signature_type="COSMIC",
                              cosmic_siglist=NULL,
                              bootstrap_scores=FALSE,
                              nparallel=4)
							  
row.names(myoutput$hrdetect_output)=samplename
if (file.exists("/mnt/vol1/projects/clingen/HRDetect/hrdetect_output.txt"))
{
write.table(myoutput$hrdetect_output,"/mnt/vol1/projects/clingen/HRDetect/hrdetect_output.txt",sep="\t",quote=F,col.names=F,append=T)
} else
{
write.table(myoutput$hrdetect_output,"/mnt/vol1/projects/clingen/HRDetect/hrdetect_output.txt",sep="\t",quote=F,col.names=NA)
}




#Il calcolo che fa per la BRCAness è:
#linear_part <- rep(intercept,nrow(hrdetect_input)) + apply(hrdetect_input*matrix(rep(features_weight,nrow(hrdetect_input)),nrow = nrow(hrdetect_input),byrow = TRUE),1,sum)
#BRCA_prob <- 1/(1+exp(-linear_part))