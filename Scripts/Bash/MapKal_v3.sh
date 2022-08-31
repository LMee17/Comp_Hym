#!usr/bin/bash

#USAGE: bash MapKal_v2.sh <batch txt file> <location of transcriptome> <directory of reads>

#Batch text file should contain the samples you wish to map, one per line
#In this script I have the sample names with an assumed prefix of fq.gz
#Similarly, the batch script should be named in such a way as it goes RunName_Batch.txt
#ie If I was working with the species Amel or else wanted my project named Proj0 the batch
#files would be named Amel_Batch.txt or Proj0_Batch.txt, respectively.
#The transcriptome will alsed be presumed to be named in the same way, ie Amel_GCF001.fasta etc
#Please ensure all directory paths end in a /

#Author: Lauren Mee, IVES, University of Liverpool 2021

filecheck=$(echo $1 | grep -c "/")
if [ $filecheck -eq 1 ]; then
	run=$(echo $1 | awk -F "/" '{print$NF}' | awk -F "_" '{print$1}')
else
	run=$(echo $1 | awk -F "_" '{print$1}')
fi

echo $run

mkdir $run/
mkdir $run/Counts/
mkdir $run/logs/

kallisto index -i $run/$run.indx "$2""$run"_*fasta

while read s; do
	echo $s
	kallisto quant -i $run/$run.indx $3$s*_1.fq.gz $3$s*_2.fq.gz -o $run/Counts/$s/ -t $4 &> $s.log
	mv $run/Counts/$s/abundance.tsv $run/Counts/$s/$s.abundance.tsv
	mv $s.log $run/logs/
done < $1
