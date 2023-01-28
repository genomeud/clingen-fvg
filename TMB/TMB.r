#Written by Fabio Marroni and Stefano Pezzella in 2020.


# DOI 10.1186/s13073-017-0424-2

#Noi per il TMB abbiamo usato SOLTANTO LE MISSENSO

library(data.table)
library(openxlsx)

#We use a literature-based proxy for length of the CDS: https://doi.org/10.1186/s13104-019-4343-8 
cdslength<-25840698 
cdslength<-cdslength/1000000 #convert to Mbases
#Files 
snvtoread<-"Sample.muTect.somatic.snv.annovar.hg19_multianno.xls"
vars<-fread(snvtoread, data.table=F)

#Important assumption about the input files
#1) It has one column named ExonicFunc, in which intergenic mutations are annotated as "."
#   if this is not the case you will need to change the symbol in the line below
exvars<-vars[vars$ExonicFunc!=".",] 
mut_type<-table(exvars$ExonicFunc)
TMB<-sum(mut_type)/cdslength

