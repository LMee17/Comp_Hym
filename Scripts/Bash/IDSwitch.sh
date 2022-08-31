#!usr/bin/bash

#USAGE bash IDSwitch.sh <Folder Transcriptomes> <Folder IDtables>

for f in $1*faa; do
	filecheck=$(echo $1 | grep -c "/")
	if [ $filecheck -eq 1 ]; then
		spec=$(echo $f | awk -F "/" '{print$NF}' | awk -F "." '{print$1}')
	else
		spec=$(echo $f | awk -F "." '{print$1}')
	fi
	echo $spec
	grep ">" $f | awk -F " " '{print$1}' | sed 's/^.//' > $spec.trans.list.tmp
#	count=1
	while read t; do
		gene=$(grep $t $2$spec*tsv | awk '{print$1}' | uniq)
		echo $gene
		let count=count+1
		sed -i "s/$t/$gene/" $f
#		if [ $count -eq 5 ]; then
#			exit
#		fi
	done < $spec.trans.list.tmp
done

rm *tmp
