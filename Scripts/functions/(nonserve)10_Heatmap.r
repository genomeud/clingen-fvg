#Perform clustering and heatmpa like my work on mozzarella

do.cluster<-function(infile="/projects/populus/ep/share/marroni/collaborations/Misson/tables/16s_merged_count_by_Family.txt",
                    metafile="/projects/populus/ep/share/marroni/collaborations/Misson/tables/metadata.csv",
					pdffile="/projects/populus/ep/share/marroni/collaborations/Misson/plots/Family_new_heatmap.pdf",
					outtable="/projects/populus/ep/share/marroni/collaborations/Misson/tables/Table_S1.txt",
               use.vst=F,onlytop=T,top=50,remove.extreme.high=F,nhigh=3,variable.name="Salinity",add.col=T,use.new.names=T,
                    new.rotation=T,add.brackets=F)
{
library(DESeq2)
library(data.table)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
mydat<-fread(infile,data.table=F)
#Change names for some classes
mydat$V1[mydat$V1=="Family_XII"]<-"Clostridiales_Family_XII"
mydat$V1[mydat$V1=="Family_XIII"]<-"Clostridiales_Family_XIII"
mydat$V1[mydat$V1=="Family_XI"]<-"Bacillales_Family_XI"
mydat$V1[mydat$V1=="JTB215"]<-"Clostridiales_JTB215"
#Remove all the rows in which the higher taxonomy level is unknown
setnames(mydat,"V1","taxonomy")
row.names(mydat)<-mydat$taxonomy
mydat$taxonomy<-mydat$"16"<-mydat$"17"<-mydat$Tot<-NULL
metadata<-fread(metafile,data.table=F)
#Use DESeq2 to normalize read counts
dds <- DESeqDataSetFromMatrix(countData = mydat,
                              colData = metadata,
                              design = ~ 1)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#We use VST to transform the data for better plotting
if(use.vst)
{
	dds <-varianceStabilizingTransformation(dds, blind=FALSE)                              
	pino<-assay(dds, normalized=TRUE)
}
if(!use.vst)
{
	# dds <- estimateSizeFactors(dds)
	# pino<-counts(dds, normalized=TRUE)                              
	# pino<-log10(0.1+pino)
	pino<-log10(0.1+mydat)
	}

tot<-apply(pino,1,sum)
pino<-pino[order(tot,decreasing=T),]
#I write back the not-log transofrmed counts, because I like it more
write.table((10^pino)-0.1,outtable,quote=F,sep="\t",col.names=NA)

if(onlytop) 
{
if(top>1) pino<-pino[1:min(top,nrow(pino)),]
}


pp<-data.frame(colnames(pino))
names(pp)<-"Sample"
smeta<-metadata[,c("Sample",variable.name)]
ff<-merge(pp,smeta,by="Sample",sort=F)
mycol<-data.frame(var=unique(ff[,variable.name]),col=brewer.pal(n=length(unique(ff[,variable.name])),"Set1"))
setnames(mycol,"var",variable.name)
ff<-merge(ff,mycol,all=T,sort=F)
mycol<-ff$col
names(mycol)<-ff$Sample
pino<-as.matrix(pino)
mycol<-list(Salinity=c("0"="red","4"="blue","9"="darkgreen","18"="purple","35"="orange"))
sal<-HeatmapAnnotation(Salinity=ff[,variable.name],col=mycol,show_legend=F,show_annotation_name=F)
legsal<-Legend(labels=names(mycol$Salinity),title="Salinity",legend_gp=gpar(fill=mycol$Salinity),
		grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
		labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 7, fontface = "bold"))
maincolful<-circlize::colorRamp2(c(min(pino),mean(pino), max(pino)), c("blue", "white", "red"))
at<-round(seq(from=min(pino),to=max(pino),length.out=5),0)
legheat<-Legend(labels=at,title="Normalized\nabundance",legend_gp=gpar(fill=maincolful(at)),
		grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
		labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 7, fontface = "bold"))
pd = packLegend(legsal, legheat)
#png(pdffile,width=13,height=13,units="cm",res=600,type="cairo")
pdf(pdffile)
par(mar=c(1,1,1,1))
cc<-Heatmap(pino, 
	row_names_gp = gpar(fontsize = 7),
	column_names_gp = gpar(fontsize = 8),
	name="Normalized\nabundance",
	width=unit(9,"cm"),
	height=unit(13,"cm"),
	top_annotation=sal,
	show_heatmap_legend = FALSE)
draw(cc,adjust_annotation_extension=TRUE,padding=unit(c(0.1,0.1,0.1,0.1),"mm"))
draw(pd,x=unit(0.92,"npc"),y=unit(0.5,"npc"))
dev.off()

}






