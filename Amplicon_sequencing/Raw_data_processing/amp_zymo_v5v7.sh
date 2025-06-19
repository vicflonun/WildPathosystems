#!/bin/bash

#Version 3.0
#Script must be run in the main_dir were data is storaged
#main_dir should contain a data folder
#TRIMMOMATIC, FLASH, VSEARCH NEEDED

echo "###### ============================ ######"
echo "###### Setting programs and folders ######"
echo "###### ============================ ######"

#Remove out files so they dont concatenate the same result
rm -r ./stats

#Seting directories
mkdir -p ./trimmed
mkdir -p ./merged
mkdir -p ./process
mkdir -p ./result
mkdir -p ./stats


#Create variables for the tools 
export trimmomatic=~/Documents/ampliconseq/Trimmomatic-0.39/trimmomatic-0.39.jar
export flash=~/Documents/ampliconseq/FLASH-1.2.11-Linux-x86_64/flash
export vsearch=~/Documents/ampliconseq/vsearch-2.22.1-linux-x86_64/bin/vsearch
#export usearch=~/Documents/ampliconseq/usearch11.0.667_i86linux32


echo pipeline for $1

if [ $1 = 16S ];
then
#set variables #avoid $o, $f, $n
#Trimmomatic
hc=20
#cr=280
q=12 #minimum quality 
l=30 #minimun length

#Filtering #more can be added if needed
#sr=19 #strip bases from right 
#sl=20 #strip bases from left
mi=325 #minumum length 240   #REf is 254bp
ma=345 #maximum length 260
ee=2 #expected error

#clustering
i=0.97 #threshold OTU for clustering

#otu classification
b=0.80 #minimum boostrap value
base=~/Documents/data.bases/rdp_16s_v18.fa #database
fi 


if [ $1 = ITS ]
then
#set variables #avoid $o, $f, $n
#Trimmomatic
hc=20
#cr=280
q=12 #minimum quality 
l=30 #minimun length


#Filtering #more can be added if needed
#sr=19 #strip bases from right 
#sl=20 #strip bases from left
mi=165 #minumum length #PREVIOUS 290 - 395. Reference is 310+-67bp aprox 240 to 380 
ma=265 #maximum length
ee=1.5 #expected error

#clustering
i=0.95 #threshold OTU for clustering

#otu table
b=0.80 #minimum boostrap value
base=~/Documents/data.bases/utax_sintax.fasta #database
fi 

if [ $2 = otu ] 
then 

echo "###### ============================ ######"
echo "###### Q trimming and merging pairs ######"
echo "###### ============================ ######"

#TRIM primers and TRIM low quality seqs
for f in $(find ./data/* -type f | grep _R1_001.fastq.gz); do
o=./trimmed/$(awk -F '/|_L001' '{print $3}' <<< "$f")
java -jar $trimmomatic PE -basein $f -baseout $o -summary ./trimmed/out.txt HEADCROP:$hc LEADING:$q TRAILING:$q SLIDINGWINDOW:4:$q MINLEN:$l AVGQUAL:$q
echo $o >> ./stats/trimmed.out.txt
cat ./trimmed/out.txt >> ./stats/trimmed.out.txt
done

#MERGE the pairs
for f in $(find ./trimmed/* -type f | grep _1P.fastq.gz) ; do 
o=./merged/$(awk -F '/|_1P' '{print $3}' <<< "$f").extendedFrags.fastq.gz
$vsearch --fastq_mergepairs $f --reverse ${f/1P/2P} --fastqout $o --log ./merged/out.txt 2> ./merged/out2.txt
cat ./merged/out.txt ./merged/out2.txt >> ./stats/merged.out.txt
done

#--fastq_allowmergestagger 


echo "###### ============================ ######"
echo "###### Q filtering and cropping     ######"
echo "###### ============================ ######"

#RENAME, CROP to size and FILTER merged pairs
for f in $(find ./merged/* -type f | grep extendedFrags.fastq.gz) ; do 
n=$(awk -F '/|.extendedFrags' '{print $3}' <<< "$f")
o=./process/${n}_quality.fastq.gz 
$vsearch --fastx_filter $f --relabel $n. --fastq_maxlen $ma --fastq_minlen $mi --fastq_maxee $ee --fastqout $o --log ./process/out.txt 
cat ./process/out.txt >> ./stats/filter.out.txt
done

echo "###### ============================ ######"
echo "###### Dereplicating sequences      ######"
echo "###### ============================ ######"

#DEREPLICATE by sample 
for f in $(find ./process/* -type f | grep _quality.fastq) ; do 
o=./process/$(awk -F '/|_quality' '{print $3}' <<< "$f")_derep.fasta
$vsearch --fastx_uniques $f --strand plus --fastaout  $o --sizeout --fasta_width 0 --log ./process/out.txt 
cat ./process/out.txt >> ./stats/derep.sample.out.txt
done 

#DEREPLICATE across samples
cat ./process/*derep.fasta > ./result/all_derep.fasta
$vsearch --derep_fulllength ./result/all_derep.fasta --strand plus --output ./result/all_unique.fasta --sizein --sizeout --fasta_width 0 --uc ./result/all-unique.uc --log ./stats/derep.all.out.txt 

echo "###### ============================ ######"
echo "###### Clustering OTUS              ######"
echo "###### ============================ ######"

# OTU CLUSTERING 
$vsearch --cluster_size ./result/all_unique.fasta --id $i --strand plus --sizein  --sizeout  --fasta_width 0 --centroids ./result/centroids.fasta --log ./result/out.txt
echo Clusters generated: $(grep -c "^>" ./result/centroids.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#Remove SINGLETONS
$vsearch --sortbysize ./result/centroids.fasta --sizein --sizeout  --fasta_width 0 --minsize 2 --output ./result/sorted.fasta --log ./result/out.txt
echo Clusters after removal of singletons: $(grep -c "^>" ./result/sorted.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#Remove CHIMERAS de novo 
$vsearch --uchime_denovo ./result/sorted.fasta --sizein --sizeout --fasta_width 0 --qmask none --nonchimeras ./result/denovo.nonchimeras.fasta --log ./result/out.txt
echo Clusters after removal of chimeras de novo: $(grep -c "^>" ./result/denovo.nonchimeras.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#RENAME OTUs
$vsearch --fastx_filter ./result/denovo.nonchimeras.fasta --fasta_width 0 --relabel OTU --fastaout ./result/otus.fasta

echo "###### ============================ ######"
echo "###### Mapping reads and OTU table  ######"
echo "###### ============================ ######"

#Construct concatenated SEMIRAW reads (merged reads before quality filtering)
for f in $(find ./merged/* -type f | grep extendedFrags.fastq.gz) ; do 
n=$(awk -F '/|.extendedFrags' '{print $3}' <<< "$f")
o=./process/${n}_renamed.fasta 
$vsearch --fastx_filter $f --fastaout $o --relabel $n.
done

cat ./process/*renamed.fasta > ./result/all_semiraw.fasta
rm  ./process/*renamed.fasta

#Create OTU TABLE based on semiraw merged pairs
$vsearch --usearch_global ./result/all_semiraw.fasta --db ./result/otus.fasta --id $i --strand plus --sizein --sizeout --fasta_width 0  --qmask none --dbmask none --otutabout ./result/otutab.txt --log ./stats/table.out.txt

echo "###### ============================ ######"
echo "###### Assign taxonomy to OTUs      ######"
echo "###### ============================ ######"

#OTU CLASSIFICATION use the appropriate database for 16S and ITS2 
$vsearch --sintax ./result/otus.fasta --db $base --tabbedout ./result/otus.sintax --strand both --sintax_cutoff $b --log ./stats/taxa.out.txt


#rm ./result/all_semiraw.fasta
rm ./result/all_derep.fasta
fi

if [ $2 = zotu ] 
then 

echo "###### ============================ ######"
echo "###### Q trimming and merging pairs ######"
echo "###### ============================ ######"

#TRIMM primers
for f in $(find ./data/* -type f | grep _R1_001.fastq.gz); do
o=./trimmed/$(awk -F '/|_L001' '{print $3}' <<< "$f")
java -jar $trimmomatic PE -basein $f -baseout $o -summary ./trimmed/out.txt CROP:$cr HEADCROP:$hc
echo $o >> ./stats/trimmed.out.txt
cat ./trimmed/out.txt >> ./stats/trimmed.out.txt
done

for f in $(find ./trimmed/* -type f | grep _1P.fastq.gz) ; do 
o=./merged/$(awk -F '/|_1P' '{print $3}' <<< "$f").extendedFrags.fastq.gz
$vsearch --fastq_mergepairs $f --reverse ${f/1P/2P} --fastqout $o --fastq_allowmergestagger --log ./merged/out.txt
cat ./merged/out.txt >> ./stats/merged.out.txt
done

echo "###### ============================ ######"
echo "###### Q filtering and cropping     ######"
echo "###### ============================ ######"

#RENAME, CROP to size and FILTER merged pairs
for f in $(find ./merged/* -type f | grep extendedFrags.fastq.gz) ; do 
n=$(awk -F '/|.extendedFrags' '{print $3}' <<< "$f")
o=./process/${n}_quality.fastq.gz 
$vsearch --fastx_filter $f --relabel $n. --fastq_maxlen $ma --fastq_minlen $mi --fastq_maxee $ee --fastqout $o --log ./process/out.txt 
cat ./process/out.txt >> ./stats/filter.out.txt
done

echo "###### ============================ ######"
echo "###### Dereplicating sequences      ######"
echo "###### ============================ ######"

#DEREPLICATE by sample 
for f in $(find ./process/* -type f | grep _quality.fastq) ; do 
o=./process/$(awk -F '/|_quality' '{print $3}' <<< "$f")_derep.fasta
$vsearch --fastx_uniques $f --strand plus --fastaout  $o --sizeout --fasta_width 0 --log ./process/out.txt 
cat ./process/out.txt >> ./stats/derep.sample.out.txt
done 

#DEREPLICATE across samples
cat ./process/*derep.fasta > ./result/all_derep.fasta
$vsearch --derep_fulllength ./result/all_derep.fasta --strand plus --output ./result/all_unique.fasta --sizein --sizeout --fasta_width 0 --uc ./result/all-unique.uc --log ./stats/derep.all.out.txt 

echo "###### ============================ ######"
echo "###### Clustering OTUS              ######"
echo "###### ============================ ######"

# OTU CLUSTERING 
$vsearch --cluster_unoise ./result/all_unique.fasta --strand plus --sizein  --sizeout  --fasta_width 0 --centroids ./result/denoised.fasta --log ./result/out.txt
echo Clusters generated: $(grep -c "^>" ./result/denoised.fasta) >> ./stats/zotu.out.txt
cat ./result/out.txt >> ./stats/zotu.out.txt

#Remove SINGLETONS
$vsearch --sortbysize ./result/denoised.fasta --sizein --sizeout  --fasta_width 0 --minsize 2 --output ./result/sorted.denoised.fasta --log ./result/out.txt
echo Clusters after removal of singletons: $(grep -c "^>" ./result/sorted.denoised.fasta) >> ./stats/zotu.out.txt
cat ./result/out.txt >> ./stats/zotu.out.txt

#Remove CHIMERAS de novo 
$vsearch --uchime2_denovo ./result/sorted.denoised.fasta --sizein --sizeout --fasta_width 0 --qmask none --nonchimeras ./result/denovo.denoised.nonchimeras.fasta --log ./result/out.txt
echo Clusters after removal of chimeras de novo: $(grep -c "^>" ./result/denovo.denoised.nonchimeras.fasta) >> ./stats/zotu.out.txt
cat ./result/out.txt >> ./stats/zotu.out.txt


#RENAME OTUs
$vsearch --fastx_filter ./result/denovo.denoised.nonchimeras.fasta --fasta_width 0 --relabel OTU --fastaout ./result/zotus.fasta

echo "###### ============================ ######"
echo "###### Mapping reads and OTU table  ######"
echo "###### ============================ ######"

#Construct concatenated SEMIRAW reads (merged reads before quality filtering)
for f in $(find ./merged/* -type f | grep extendedFrags.fastq.gz) ; do 
n=$(awk -F '/|.extendedFrags' '{print $3}' <<< "$f")
o=./process/${n}_renamed.fasta 
$vsearch --fastx_filter $f --fastaout $o --relabel $n.
done

cat ./process/*renamed.fasta > ./result/all_semiraw.fasta
rm  ./process/*renamed.fasta

#Create OTU TABLE based on semiraw merged pairs
$vsearch --usearch_global ./result/all_semiraw.fasta --db ./result/zotus.fasta --id $i --strand plus --sizein --sizeout --fasta_width 0  --qmask none --dbmask none --otutabout ./result/zotutab.txt --log ./stats/ztable.out.txt

echo "###### ============================ ######"
echo "###### Assign taxonomy to OTUs      ######"
echo "###### ============================ ######"

#OTU CLASSIFICATION use the appropriate database for 16S and ITS2 
$vsearch --sintax ./result/zotus.fasta --db $base --tabbedout ./result/zotus.sintax --strand both --sintax_cutoff $b --log ./stats/taxa.out.txt

rm ./result/all_semiraw.fasta
rm ./result/all_derep.fasta

fi








