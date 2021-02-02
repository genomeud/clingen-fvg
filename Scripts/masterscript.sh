#Directory con i file: /mnt/vol1/projects/clingen/
#PER AVVIARE R DA CONDA:
export PATH=$PATH:/mnt/vol2/conda/miniconda3/bin
source activate r_3.6.1
R
#PER USCIRE DA CONDA
conda deactivate
#PER CREARE UN AMBIENTE NUOVO
conda create --name NOMEAMBIENTE
#PER ATTIVARE AMBIENTE
source activate nomeambiente

#COMANDI SCREEN:
#Ctrl + D : chiudere screen
#Crtl + A + ...
#C: crea nuovo screen
#N: passa allo screen successivo (o torna al primo)
#P: passa allo screen precedente
#D: detached (si esce da screen)

#Per controllare spazio disco
#df -h

#Write number of somatic SNPs in each sample CONTA GLI SNP (RIGHE) IN OGNI VCF, LA PIPE FUNZIONA SOLO CON EGREP
cd $VCFDIR
[ -f ${VCFDIR}/Tables/SNV_number.txt ] && rm ${VCFDIR}/Tables/SNV_number.txt #rimuove il txt se esiste già
for aaa in *.vcf; do  egrep -v "##|#CHROM" $aaa | echo -e ${aaa/T.muTect.somatic.snv.vcf/} '\t' $(wc -l) >> ${VCFDIR}/Tables/SNV_number.txt; done
#DA QUI SI VEDE CHE IN FVG3 HANNO ACCOPPIATO I CAMPIONI DI 2 PERSONE DIVERSE (TROPPE MUTAZIONI), QUINDI IN MOLTE ANALISI VIENE ESCLUSO

#RICHIAMA LO SCRIPT CHE INDICA LA DISTRIBUZIONE DELLE MUTAZIONI. PARAMETRI I, B, T DOPO IL PATH
Rscript /mnt/vol1/projects/clingen/functions/01_mutation_distribution.r

#Richiama lo stesso script ma utilizza le 6 mutazioni "canoniche" invece che tutte e 12
Rscript /mnt/vol1/projects/clingen/functions/02_mutation_distribution_strand_ind.r

#Analisi su MutationalPatterns
Rscript /mnt/vol1/projects/clingen/functions/03_play_mutational_pattern.r

#Preparare la reference per CallableLoci

gunzip /mnt/vol1/projects/clingen/Reference/human_g1k_v37_decoy.fasta.gz

	#Create fai
module load sw/bio/samtools/0.1.19
samtools faidx human_g1k_v37_decoy.fasta
	#Create dict
module load sw/bio/picard-tools/1.88
java -jar /iga/stratocluster/packages/sw/bio/picard-tools/1.88/CreateSequenceDictionary.jar \
R= /mnt/vol1/projects/clingen/Reference/human_g1k_v37_decoy.fasta \
O= /mnt/vol1/projects/clingen/Reference/human_g1k_v37_decoy.dict


#Test callableloci su FVG11 (coverage: 15-75) conviene non impostare un loop ma farli uno alla voltta perchè ci mette 1-2 gg
module load sw/bio/gatk/3.3-0
java -jar /iga/stratocluster/packages/sw/bio/gatk/3.3-0/GenomeAnalysisTK.jar -T CallableLoci \
-I /mnt/vol1/projects/clingen/Bam/FVG29T.final.bam -R /mnt/vol1/projects/clingen/Reference/human_g1k_v37_decoy.fasta \
-summary /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/CallableLoci_output/FVG29T.final.summary.txt \
-o /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/CallableLoci_output/FVG29T.final.bed \
-minDepth 15 -maxDepth 75 -format BED #-nt 4 #nt = numero di processori


#Valutazione distribuzione mutazioni (necessita dell'output di CallableLoci) in loop per tutti i campioni
Rscript /mnt/vol1/projects/clingen/functions/04_CallableLoci.r

#Produce un output tabulare
Rscript /mnt/vol1/projects/clingen/functions/16_CallableLoci_tab.r

#Creazione BEDPE per HRDetect
Rscript /mnt/vol1/projects/clingen/functions/05_delly_to_HRDetect.r

#Per usare control-FREEC (calcolo cnv). QUI DENTRO FUNZIONA SAMTOOLS SENZA CARICARE MODULI
export PATH=$PATH:/mnt/vol2/conda/miniconda3/bin
source activate controlfreec
#Estrae la lista di SNP noti in umano
gunzip /mnt/vol1/projects/clingen/Reference/hg19_snp142.SingleDiNucl.1based.txt.gz
cd /mnt/vol1/projects/clingen/Control-FREEC/
#conda install -c bioconda samtools openssl=1.0 #downgradare OpenSSL aveva risolto un bug

#è stato necessario rimuovere dal file hg19.len tutti gli scaffold e il cromosoma y (FORSE?)

#Creiamo dei VCF adatti a Control-FREEC
for file in /mnt/vol1/projects/clingen/Mutation/Mutation_SNP/*T.GATK.snp.vcf.gz
do
shortname=$(basename $file)
shortname=${shortname/.gz/} #elimina il .gz alla fine
if [ "$shortname" != "FVG3T.GATK.snp.vcf" ]
then
zcat $file | sed -e 's/##contig=<ID=/##contig=<ID=chr/g' -e 's/^/chr/g' -e 's/chr#/#/g' > /mnt/vol1/projects/clingen/Mutation/Mutation_SNP/SNP_FREEC/$shortname
fi
done

# Try to obtain the full data GIA FATTO PER TUTTI, LA ROBA DOPO NO
#Needed only once to create the reference with names chr...
#sed -e 's/>/>chr/g' /mnt/vol1/projects/clingen/Reference/human_g1k_v37_decoy.fasta > /mnt/vol1/projects/clingen/Reference/human_g1k_v37_decoy_chr.fasta
#samtools faidx /mnt/vol1/projects/clingen/Reference/human_g1k_v37_decoy_chr.fasta

for file in /mnt/vol1/projects/clingen/Bam/*.final.bam
do
shortname=$(basename $file)
shortname=${shortname/.final/_tmp}
if [ "$shortname" != "FVG3A_tmp.bam" ] && [ "$shortname" != "FVG3T_tmp.bam" ]
then
samtools view -b --threads 24 -L /mnt/vol1/projects/clingen/Reference/hg19_onlychr.bed $file > /mnt/vol1/projects/clingen/Control-FREEC/$shortname
fi
done

#Questo ci mette TROPPO tempo, per comodità aprire più screen ed eseguire i comandi parziali riportati sotto

for file in /mnt/vol1/projects/clingen/Control-FREEC/*_tmp.bam
do
shortname=$(basename $file)
shortname=${shortname/_tmp/_chr}
samtools view -h --threads 8 $file | awk 'BEGIN{ FS=OFS="\t"} {if ( $1 !~ /@/ ) {$3="chr"$3}; if ( $1 ~ /@/ ) {gsub("SN:", "SN:chr", $2)}; if (( $1 !~ /@/ ) && ($7!~ /=/ )) {$7="chr"$7};print $0}' | samtools view -b --threads 8  - > /mnt/vol1/projects/clingen/Control-FREEC/$shortname
done

##Comandi parziali##

for number in 9A 9T 10A 10T 11A 11T 12A 12T 13A 13T 14A 14T 15A 15T
do
file=/mnt/vol1/projects/clingen/Control-FREEC/FVG${number}_tmp.bam
shortname=$(basename $file)
shortname=${shortname/_tmp/_chr}
samtools view -h --threads 8 $file | awk 'BEGIN{ FS=OFS="\t"} {if ( $1 !~ /@/ ) {$3="chr"$3}; if ( $1 ~ /@/ ) {gsub("SN:", "SN:chr", $2)}; if (( $1 !~ /@/ ) && ($7!~ /=/ )) {$7="chr"$7};print $0}' | samtools view -b --threads 8  - > /mnt/vol1/projects/clingen/Control-FREEC/$shortname
done

for number in 17A 17T 18A 18T 20A 20T 21A 21T 22A 22T
do
file=/mnt/vol1/projects/clingen/Control-FREEC/FVG${number}_tmp.bam
shortname=$(basename $file)
shortname=${shortname/_tmp/_chr}
samtools view -h --threads 8 $file | awk 'BEGIN{ FS=OFS="\t"} {if ( $1 !~ /@/ ) {$3="chr"$3}; if ( $1 ~ /@/ ) {gsub("SN:", "SN:chr", $2)}; if (( $1 !~ /@/ ) && ($7!~ /=/ )) {$7="chr"$7};print $0}' | samtools view -b --threads 8  - > /mnt/vol1/projects/clingen/Control-FREEC/$shortname
done

for number in 17A 17T 18A 18T 20A 20T 21A 21T 22A 22T
do
file=/mnt/vol1/projects/clingen/Control-FREEC/FVG${number}_tmp.bam
shortname=$(basename $file)
shortname=${shortname/_tmp/_chr}
samtools view -h --threads 8 $file | awk 'BEGIN{ FS=OFS="\t"} {if ( $1 !~ /@/ ) {$3="chr"$3}; if ( $1 ~ /@/ ) {gsub("SN:", "SN:chr", $2)}; if (( $1 !~ /@/ ) && ($7!~ /=/ )) {$7="chr"$7};print $0}' | samtools view -b --threads 8  - > /mnt/vol1/projects/clingen/Control-FREEC/$shortname
done

for number in 23A 23T 24A 24T 25A 25T 26A 26T 27A 27T 28A 28T 29A 29T
do
file=/mnt/vol1/projects/clingen/Control-FREEC/FVG${number}_tmp.bam
shortname=$(basename $file)
shortname=${shortname/_tmp/_chr}
samtools view -h --threads 8 $file | awk 'BEGIN{ FS=OFS="\t"} {if ( $1 !~ /@/ ) {$3="chr"$3}; if ( $1 ~ /@/ ) {gsub("SN:", "SN:chr", $2)}; if (( $1 !~ /@/ ) && ($7!~ /=/ )) {$7="chr"$7};print $0}' | samtools view -b --threads 8  - > /mnt/vol1/projects/clingen/Control-FREEC/$shortname
done

for number in 2A 2T 4A 4T 5A 5T 6A 6T 7A 7T 8A 8T
do
file=/mnt/vol1/projects/clingen/Control-FREEC/FVG${number}_tmp.bam
shortname=$(basename $file)
shortname=${shortname/_tmp/_chr}
samtools view -h --threads 8 $file | awk 'BEGIN{ FS=OFS="\t"} {if ( $1 !~ /@/ ) {$3="chr"$3}; if ( $1 ~ /@/ ) {gsub("SN:", "SN:chr", $2)}; if (( $1 !~ /@/ ) && ($7!~ /=/ )) {$7="chr"$7};print $0}' | samtools view -b --threads 8  - > /mnt/vol1/projects/clingen/Control-FREEC/$shortname
done

#Pulizia file temporanei (NON ESEGUIRE PER ORA)  
#rm /mnt/vol1/projects/clingen/Control-FREEC/*_tmp.bam

#Eseguire Control-FREEC

for file in /mnt/vol1/projects/clingen/Control-FREEC/*T_chr.bam
do
shortname=$(basename $file)
shortname=${shortname/T_chr.bam/}
sed -e "s/FVG9/$shortname/g" /mnt/vol1/projects/clingen/Control-FREEC/config_BAF_clingen.txt > /mnt/vol1/projects/clingen/Control-FREEC/config_BAF_clingen_$shortname.txt
freec --conf /mnt/vol1/projects/clingen/Control-FREEC/config_BAF_clingen_$shortname.txt #Quando si usa sed con una variabile servono le ", non le '
done

########TORNARE IN AMBIENTE CONDA DI R ###########
conda deactivate
source activate r_3.6.1
#conda install -c bioconda samtools
#conda install -c bioconda htslib

#Creazione CNV per HRDetect a partire dai BAF
#/mnt/vol1/projects/clingen/Control-FREEC/output/FVG9T_chr.bam_CNVs
#/mnt/vol1/projects/clingen/Control-FREEC/cnv_for_HRDetect/FVG9T_for_HRDetect.txt

for file in /mnt/vol1/projects/clingen/Control-FREEC/output/*CNVs
do
shortname=$(basename $file)
shortname=${shortname/chr.bam_CNVs/for_HRDetect.txt}
Rscript /mnt/vol1/projects/clingen/functions/06_cnv_to_HRDetect.r -T $file -O /mnt/vol1/projects/clingen/Control-FREEC/cnv_for_HRDetect/$shortname
done

#Usare HRDetect
#Creiamo il tabix per i vcf degli SNP e delle INDEL   e poi lanciamo HRDetect
##NON HA FATTO IL 4T, NON SAPPIAMO PERCHE, LO LASCIAMO INDIETRO

rm /mnt/vol1/projects/clingen/HRDetect/hrdetect_output.txt #è necessario rimuovere l'output vecchio prima di produrre quello nuovo
for file in /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/*T.muTect.somatic.snv.vcf
do
shortname=$(basename $file)
shortname=${shortname/.muTect.somatic.snv.vcf/}
bgzip -c /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/$shortname.muTect.somatic.snv.vcf > /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/$shortname.muTect.somatic.snv.vcf.gz
tabix -p vcf /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/$shortname.muTect.somatic.snv.vcf.gz
bgzip -c /mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/$shortname.Strelka.somatic.indel.vcf > /mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/$shortname.Strelka.somatic.indel.vcf.gz
tabix -p vcf /mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/$shortname.Strelka.somatic.indel.vcf.gz

Rscript /mnt/vol1/projects/clingen/functions/07a_comandi_HRDetect.r -S FVG4T $shortname
done

#Calcolo TMB 
	#creare FAI della reference
	samtools faidx /mnt/vol1/projects/clingen/Reference/GCF_000001405.25_GRCh37.p13_cds_from_genomic.fna
Rscript /mnt/vol1/projects/clingen/functions/08_ecTMB.r #in realtà non usiamo ecTMB

#Dimostrare che il campione 3 è frutto di appaiamento errato tramite SNPRelate

Rscript /mnt/vol1/projects/clingen/functions/09_SNPRelate.r

###############
#Analisi Driver Genes
###############

tar -xvzf /mnt/vol1/projects/clingen/driver_mutations/SPOT-1D-local.tar.gz -C /mnt/vol1/projects/clingen/driver_mutations/ #in realtà per ora lasciamo perdere SPOT1D

#Andiamo avanti con CanDrA
tar -xvjf /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0.tar.bz2 -C /mnt/vol1/projects/clingen/driver_mutations/
tar -xvjf /mnt/vol1/projects/clingen/driver_mutations/BRCA.tar.bz2 -C /mnt/vol1/projects/clingen/driver_mutations/

#compilare tabix
cd /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/tabix-0.2.6/
make
cd

#copiare il dataset di BRCA nella directory giusta
mv /mnt/vol1/projects/clingen/driver_mutations/BRCA /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/

#prova
#perl /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/open_candra.pl BRCA /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/demo_input.txt > /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/demo.annotated

#creazione input per CanDrA
Rscript /mnt/vol1/projects/clingen/functions/10_driver_genes.r

#Ora facciamo CanDrA
for file in /mnt/vol1/projects/clingen/driver_mutations/dfs_for_candra/*.txt
do
shortname=$(basename $file)
shortname=${shortname/.txt/}
perl /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/open_candra.pl BRCA /mnt/vol1/projects/clingen/driver_mutations/dfs_for_candra/$shortname.txt > /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/output/$shortname.annotated
done

#Da FVG10 ne vengono solo due significative (score < 0,05), quello driver significativo è NBPF10 che è noto in letteratura come effettivamente correlato a breast cancer
#CanDrA usa tanti software tutti assieme e fa una media pesata. Non specifica quali features usa per BRCA.
#Dopo faremo merge, prendiamo quelli che hanno pvalue significativo e facciamo table per vedere quali sono e quali sono rappresentati più volte

#Analisi dei risultati
Rscript /mnt/vol1/projects/clingen/functions/11_driver_genes_2.r

#Ricerca mutazioni intergeniche comuni a più campione
Rscript /mnt/vol1/projects/clingen/functions/12_intergenic_mutations.r

##PER IL FUTURO
#1)Cercare mutazioni comuni a più campioni che NON si trovano nei geni (serve a dare un senso all'utilizzo del WGS nel progetto clingen)
#2)Lista di tutti i geni che contengono una mutazione somatica
#3)Cercare mutazioni nei promotori e mutazioni nonsenso (vedere in quali geni sono, se c'è qualche gene interessante etc etc)

#Vediamo se delle mutazioni note non sono presenti nel database di CanDRa
cd /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/
zcat sorted_CanDrA_BRCA.gz | grep -w -f names_intogen.tsv | cut -f5 | sort | uniq > results.txt
#Vediamo che dei 99 geni del nostro elenco ne mancano solo 2, quindi non è quello il problema. Adesso facciamo il confronto con la lista che tiene conto delle mutazioni escluse da CanDRA per vedere se le mutazioni sono presenti nei nostri campioni o no
cat /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/output/final.txt | grep -w -f names_intogen.tsv | cut -f6 | sort | uniq > results2.txt
#ora confrontiamo direttamente con i file di annotazione
for file in /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation/*.xls
do
if [[ $file != "/mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation/FVG3T*.xls" ]]
then
cat $file >> /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/results3.txt
fi
done
cut -f8 results3.txt | grep -w -f names_intogen.tsv | sort | uniq > anno_vs_intogen.txt
#questo ultimo output non va bene, lo sistemiamo
cat <(cut -d"," -f1 anno_vs_intogen.txt) <(cut -d"," -f2 anno_vs_intogen.txt) <(cut -d"," -f3 anno_vs_intogen.txt) | sort | uniq > anno_vs_intogen2.txt
#serve un ultimo ritocco
cat <(cut -d"-" -f1 anno_vs_intogen2.txt) <(cut -d"-" -f2 anno_vs_intogen2.txt) <(cut -d"-" -f3 anno_vs_intogen2.txt) | sort | uniq > anno_vs_intogen3.txt
#ora sono 99 quindi nelle annotazioni i geni sono presenti tutti, quindi non è quello il problema
cut -f8 results3.txt | sort | uniq | wc -l
#non è sufficiente vedere una mutazione somatica in un gene per stabilire se è driver, perchè dall'ultimo comando abbiamo visto che vi sono mutazioni in quasi tutti i geni del trascrittoma e dovrebbero risultare tutti driver se così fosse
#abbiamo delle mutazioni che candra non conosce, nemmeno come passenger, potrebbe essere dovuto all'unicità delle mutazioni tra i vari individui (CandrA richiede la stessa identica mutazione nello stesso identico punto)
#abbiamo il database di intogen di driver mutations del 2016, estraiamo da results3 tutte le mutazioni di PIK3CA, le scriviamo in un file temporaneo e controlliamo dove sono le mutazioni.
cut -f1-8 results3.txt | grep PIK3CA > PIK3CA_muts_intogen.txt
#già guardando le prime 3 sono tutte driver, controlliamo se le abbiamo perse con candra
zcat sorted_CanDrA_BRCA.gz | grep PIK3CA > PIK3CA_muts_candra.txt
cd ~
#per cercarci a mano la corrispondenza tra le nostre mutazioni e quelle di intogen passiamo ad R
Rscript /mnt/vol1/projects/clingen/functions/13_candra_vs_intogen.r
#bisogna creare una copia dei file di annotazione che presenti una colonna con il sample name
for file in /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation/*.xls
do
if [[ $file != /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation/FVG3T*.xls ]]
then
shortname=$(basename $file)
shortname=${shortname/.muTect.somatic.snv.annovar.hg19_multianno.xls/}
if [[ $file != /mnt/vol1/projects/clingen/Somatic/Somatic_SNV/Annotation/FVG10T*.xls ]]
then
sed "1d" $file | sed "s/$/\t$shortname/" >> /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsSNP.txt
else
sed "s/$/\t$shortname/" $file > /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsSNP.txt
fi
fi
done
###Facciamo lo stesso per le INDEL###
for file in /mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/Annotation/*.xls
do
if [[ $file != /mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/Annotation/FVG3T*.xls ]]
then
shortname=$(basename $file)
shortname=${shortname/.Strelka.somatic.indel.annovar.hg19_multianno.xls/}
if [[ $file != /mnt/vol1/projects/clingen/Somatic/Somatic_INDEL/Annotation/FVG10T*.xls ]]
then
sed "1d" $file | sed "s/$/\t$shortname/" >> /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsINDEL.txt
else
sed "s/$/\t$shortname/" $file > /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsINDEL.txt
fi
fi
done
###Ora le CNV###
for file in /mnt/vol1/projects/clingen/Somatic/Somatic_CNV/Annotation/*.xls
do
if [[ $file != /mnt/vol1/projects/clingen/Somatic/Somatic_CNV/Annotation/FVG3T*.xls ]]
then
shortname=$(basename $file)
shortname=${shortname/.final.bam_CNVs.final.hg19_multianno.xls/}
if [[ $file != /mnt/vol1/projects/clingen/Somatic/Somatic_CNV/Annotation/FVG10T*.xls ]]
then
sed "1d" $file | sed "s/$/\t$shortname/" >> /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsCNV.txt
else
sed "s/$/\t$shortname/" $file > /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsCNV.txt
fi
fi
done
###SV###
for file in /mnt/vol1/projects/clingen/Somatic/Somatic_SV/Annotation/*.xls
do
if [[ $file != /mnt/vol1/projects/clingen/Somatic/Somatic_SV/Annotation/FVG3T*.xls ]]
then
shortname=$(basename $file)
shortname=${shortname/.delly.somatic.sv.hg19_multianno.xls/}
if [[ $file != /mnt/vol1/projects/clingen/Somatic/Somatic_SV/Annotation/FVG10T*.xls ]]
then
sed "1d" $file | sed "s/$/\t$shortname/" >> /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsSV.txt
else
sed "s/$/\t$shortname/" $file > /mnt/vol1/projects/clingen/driver_mutations/CanDrA.v1.0/database/BRCA/resultsSV.txt
fi
fi
done

Rscript /mnt/vol1/projects/clingen/functions/14_variants_in_genes.r #RICERCA VARIANTI NEI GENI

Rscript /mnt/vol1/projects/clingen/functions/15_plots_variants_in_genes.r #Grafici varianti nei geni

Rscript /mnt/vol1/projects/clingen/functions/17_10x.r #Verifica traslocazione anomala con 10x
