#usr/bin/bash

# USAGE bash GFFread.sh <gff> <proteome>

# FIRST MUST REMOVE COMMENTED OUT LINES AT THE BEGINNING OF THE FILE SO THAT THIS DOESN'T AFFECT AWK ETC


grep "#" $1 > tmp
check=$(wc -l tmp | awk '{print$1}')

if [ $check -gt 0 ]; then
        sed -n '/^#/!p' $1 > edit.gff3
        rm tmp
else
        cp $1 edit.gff3
        rm tmp
fi


# START WITH GETTING ALL GENE IDS
awk '$3=="gene"' edit.gff3 | sed -n -e 's/^.*gene=//p' | awk -F ";" '{print$1}' | sort -u > genelist.txt

filecheck=$(echo $1 | grep -c "/")
if [ $filecheck -eq 1 ]; then
	species=$(echo $1 | awk -F "/" '{print$NF}' | awk -F "_" '{print$1}')
else
	species=$(echo $1 | awk -F "_" '{print$1}')
fi

echo $species

printf "%s\t%s\t%s\n" "GeneID" "TranscriptID" "ProteinID" >> $species.MasterIDtable.tsv

count=1
while read g; do
	tot=$(wc -l genelist.txt | awk '{print$1}')
	echo "Gene $count / $tot"
        gene=$g
        echo $g
        grep "gene=$g;" edit.gff3 | grep "transcript_id=" | sed -n -e 's/^.*transcript_id=//p' | sort -u | sort > $g.translist
        grep "gene=$g;" edit.gff3 | grep "Parent=" | sed -n -e 's/^.*Parent=rna/rna/p' | awk -F ";" '{print$1}' | sort -u > $g.isolist.txt
#        isono=$(wc -l $g.isolist.txt | awk '{print$1}')
        while read t; do
                i=$(grep "$t" edit.gff3 | grep "Parent=" | sed -n -e 's/^.*Parent=rna/rna/p' | awk -F ";" '{print$1}' | sort -u)
                echo $i
		trans=$t
                grep "$i;" edit.gff3 | grep "protein_id=" | sed -n -e 's/^.*protein_id=//p' | sort -u > $i.pep.check
		check=$(grep -c ";" $i.pep.check)
		if [ $check -gt 0 ]; then
			pep=$(grep "$i;" edit.gff3 | grep "protein_id=" | sed -n -e 's/^.*protein_id=//p' | awk -F ";" '{print $1}' | sort -u)
			echo "$pep"
		else
			pep=$(grep "$i;" edit.gff3 | grep "protein_id=" | sed -n -e 's/^.*protein_id=//p'| sort -u)
			echo "$pep"
		fi
		rm $i.pep.check
                printf "%s\t%s\t%s\n" "$gene" "$trans" "$pep" >> $species.MasterIDtable.tsv
        done < $g.translist
        rm *translist
        rm $g.isolist.txt
	let count=count+1
done < genelist.txt

rm edit.gff3
rm genelist.txt
