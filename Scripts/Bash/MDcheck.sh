#!usr/bin/bash

#USAGE bash MDcheck.sh <directory containing raw data directories>

#A script designed to go through the raw_data directories of Novogene transcriptomic sequence data,
#perform MD5sum checks and compare to those provided by Novogene to assess whether there has been
#any issues in transferring through the files
#Assumed that each set of raw data files are within their own folders and follow the directory pattern
#provided by Novogene
#ie I would save the zip file provided within my own directory (say SetA), the zip file would contain
#two first directories, one of which contains the raw data under the name "raw_data"
#Thus this script would be directed to a master directory (ie Raw) which contains directories that 
#then contain the zip file and its corresponding unzipped directory, which in turn contains a 
#"raw_data" directory which (finally) contains the directories of each fq.gz file per sample.

#Move through the directories, logging what samples are present

printf "%s\t%s\t%s\n" "Set" "Sample" "Decision" >> MD5.report

for d in $1*/; do
	echo "Beginning MD5 check of raw data files"
	set=$(echo $d | awk -F "/" '{print$(NF-1)}')
	for s in $d*/raw_data/*/; do
		sample=$(echo $s | awk -F "/" '{print$(NF-1)}')
		echo "Checking $sample ..."
		md5sum -b $s$sample*1.fq.gz | awk '{print$1}' >> $sample.MD5.check.tmp
		md5sum -b $s$sample*2.fq.gz | awk '{print$1}' >> $sample.MD5.check.tmp
		check=0
		awk '{print$1}' "$s"MD5.txt > $sample.MD5.ori.tmp
		diff -q $sample.MD5.check.tmp $sample.MD5.ori.tmp || check=1
		if [ $check -eq 1 ]; then
			dec=$(echo "Fail")
		else
			dec=$(echo "Pass")
		fi
		printf "%s\t%s\t%s\n" "$set" "$sample" "$dec" >> MD5.report
		rm *tmp
	done
done

echo "MD5sum checks are complete. Please refer to MD5.report to view results."
