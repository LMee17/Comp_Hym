#!usr/bin/bash

#USAGE bash TrimLoop.sh <directory containing raw data directories> <output>

#A script designed to go through the raw_data directories of Novogene transcript reads and trim adapter
#sequences
#Assumed that each set of raw data files are within their own folders and follow the directory pattern
#provided by Novogene
#ie I would save the zip file provided within my own directory (say SetA), the zip file would contain
#two first directories, one of which contains the raw data under the name "raw_data"
#Thus this script would be directed to a master directory (ie Raw) which contains directories that 
#then contain the zip file and its corresponding unzipped directory, which in turn contains a 
#"raw_data" directory which (finally) contains the directories of each fq.gz file per sample.

trim="/pub37/lmee/bin/TrimGalore-0.6.6/trim_galore"

now=$(date '+%c')
echo "Trimming began at $now" >> $2trim.report

for d in $1*/; do
	set=$(echo $d | awk -F "/" '{print$(NF-1)}')
	echo $set
	echo "Preparing to trim read files from $set"
	clock=$(date '+%c')
	echo "$set began at $clock" >> $2trim.report
	for s in $d*/raw_data/*/; do
		tic=$(date '+%c')
		sample=$(echo $s | awk -F "/" '{print$(NF-1)}')
		echo "Trimming $sample ..."
		echo "$sample trim began at $tic" >> $2trim.report
		$trim --fastqc -o $2 --paired "$s$sample""_1.fq.gz" "$s$sample""_2.fq.gz"
		tok=$(date '+%c')
		echo "$sample trim completed at $tok" >> $2trim.report
	done
	time=$(date '+%c')
	echo "$set complete at $time" >> $2trim.report
done

mkdir $2Reports/
mv $2*report.txt

echo "Read trims are complete. Please refer to your output directory to view results."
