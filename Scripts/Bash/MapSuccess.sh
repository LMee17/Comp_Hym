#!usr/bin/bash

#USAGE: bash MapSuccess.sh <list of targets>

#Assumed that this is ran after MapKal_v3.sh
#Assumed that each species / project being assessed have a directory from the current directory
#Ie ./Amel/ or ./Bter/
#Each of these then contain a folder called logs/ that have the required log files
#List of targets is a text file with the project/species name, one per line, exactly as the directories are named
#Ie to assess Amel/ directory, then one line = Amel

printf "%s\t%s\t%s\t%s\t%s\n" "SampleSet" "Sample" "Reads" "MappedReads" "PercMapped" > PercMapped.tsv

while read t; do
	for l in $t/logs/*log; do
		sample=$(echo $l | awk -F "/" '{print$NF}' | awk -F "." '{print$1}')
		readtot=$(grep "processed" $l | awk -F "reads," '{print$1}' | tr -d -c 0-9)
		maptot=$(grep "processed" $l | awk -F "reads," '{print$2}' | tr -d -c 0-9)
		prop=$(echo "$maptot / $readtot" | bc -l)
		permap=$(echo "$prop * 100" | bc -l)
		echo $permap
		printf "%s\t%s\t%s\t%s\t" $t $sample $readtot $maptot >> PercMapped.tsv
		printf "%.2f\n" $permap >> PercMapped.tsv
	done
done < $1
