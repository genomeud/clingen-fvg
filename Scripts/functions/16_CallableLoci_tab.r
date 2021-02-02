library(data.table)
wd<-getwd()
setwd("/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/CallableLoci_output/")
summaries<-dir(pattern=".final.summary.txt")
m<-matrix(0, ncol=length(summaries), nrow=6)
rownames(m)<-c("REF_N","CALLABLE","NO_COVERAGE","LOW_COVERAGE","EXCESSIVE_COVERAGE","POOR_MAPPING_QUALITY")
m<-data.frame(m)
for (aaa in 1:length(summaries))
{
file1<-fread(summaries[aaa], data.table=F)
name1<-gsub(".final.summary.txt","",basename(summaries[aaa]))
names(m)[aaa]<-name1
m[,aaa]<-file1$nBases
m[,aaa]<-as.numeric(m[,aaa])
}
m<-t(m)
m<-data.frame(m)
m$CALLABLE_ON_TOTAL<-m$CALLABLE/(m$REF_N + m$CALLABLE + m$NO_COVERAGE + m$LOW_COVERAGE + m$EXCESSIVE_COVERAGE + m$POOR_MAPPING_QUALITY)
write.table(m, "/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/CallableLoci_output/tab.txt", quote=F, sep="\t")
setwd(wd)