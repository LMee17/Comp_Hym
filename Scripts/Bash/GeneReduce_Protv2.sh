#usr/bin/bash

#USAGE bash GeneReduce_Protv2.sh <species.txt> <IDtable directory> <proteome directory>

#A script to reduce protomes to longest isoform per gene
#Requires a table with Gene and corresponding ProteinIDs in two different columns (GENE \t TRANSCRIPT \t PROTEIN)
#Requires the proteome that is to be reduced
#Requires a list of the species of proteomes to be reduced
#Species should be one per line
#Proteomes should begin with species (as appears on line of Species.txt file)
#I.e. Amel as an input should correspond with Amel.faa or Amel_GCF.faa.
#Proteomes must end in .faa
#Produces an output folder of reduced proteomes

while read s; do
	echo "Assessing $s..."
	#take only necessary information from input tables (GENE \t PROTEIN)
        awk '{print$1,$3}' $2/$s*tsv > $s.ID.edit.txt
	#generate a list of just genes
        awk '{print $1}' $2/$s*tsv | sort | uniq > $s.gene.list
	gentot=$(wc -l $s.gene.list | awk '{print$1}')
	echo "$gentot genes to resolve"
        #set the counter
        count=1
        #begin reading through the  entire gene list for this species
	while read g; do
		echo "Gene $count / $gentot : $g..."
                #for each gene, create a temporary text file of the gene and corresponding isoforms(s)
		grep "^$g " $s.ID.edit.txt > $s.$g.tmp
                #count the number of isoforms
                check=$(wc -l $s.$g.tmp | awk '{print$1}')
		echo $check
		#check that there are proteins associated with these genes to be reduced; if not; skip
		protcheck=$(grep -c "P_" $s.$g.tmp)
		echo $protcheck
		if [ $protcheck -lt 1 ]; then
			echo "No Proteins are associated with this gene, skipping..."
			let count=count+1
			rm *tmp
			continue
		else
			#check if there's only one isoform
			if [ $protcheck -eq 1 ]; then
				orth=$(grep "^$g " $s.ID.edit.txt | awk '{print $2}')
                        	echo "Only one isoform, $orth, recording ..."
			else
                        	echo "$g has more than one isoform. Resolving"
                        	#check: NP
                        	#NP are the best annotated isoforms using NCBI transcriptomes. If there are any NP_ isoform(s) these should be assessed first.
                        	npcheck=$(grep -c "NP_" $s.$g.tmp)
                        	#if there are no NP_ isoforms, that we move to taking the longest isoform
                        	if [ $npcheck -eq 0 ]; then
                                	echo "No official isoform, assessing by isoform length..."
                                	#make a temporary list of the transcripts to be assessed
                                	awk '{print$2}' $s.$g.tmp > $s.$g.target.tmp
                                	#run through each transcript, making an individual fasta file for each isoform
                                	while read t; do
                                        	echo $t > $s.$g.$t.tmp
                                        	python ~/Scripts/get_seq2.py $3/$s*fasta $s.$g.$t.tmp $s.$g.$t.fasta.tmp
                                        	#remove the top line of the fasta file
                                        	sed -i '1d' $s.$g.$t.fasta.tmp
                                        	#count what remains - the nucleotides
                                        	charcount=$(wc $s.$g.$t.fasta.tmp | awk '{print$3}')
                                        	#record the size of these isoforms in a comparison file
                                        	printf "%s\t%d\n" "$t" "$charcount" >> $s.$g.comp.tmp
                                	done < $s.$g.target.tmp
                                	#sort the isoforms by size, taking the bottom (longest) value
                                	orth=$(sort -nk 2 $s.$g.comp.tmp | tail -n 1 | awk '{print$1}')
                                	#this transcript is recorded as the preferred isoform for this gene
                                	echo "$orth is the longest isoform, recording..."
                        	else
                                	#if there is NM_ isoforms, we check to see if there are one or multiple
                                	if [ $npcheck -eq 1 ]; then
                                        	#if there is only one, this is recorded as the preferred isoform.
                                        	orth=$(grep "NP_" $s.$g.tmp | awk '{print$2}')
                                        	echo "Official isoform present, $orth, recording ..."
                                		#if there is more than one, we repeat the length comparison above using just the NP_ isoforms
                                	else
                                        	echo "More than one official isoform present, assessing by isoform length..."
                                        	grep "NP_" $s.$g.tmp | awk '{print$2}' > $s.$g.N.tmp
                                        	while read t; do
                                                	echo $t > $s.$g.$t.tmp
                                                	python ~/Scripts/get_seq2.py $3/$s*fasta $s.$g.$t.tmp $s.$g.$t.fasta.tmp
                                                	sed -i '1d' $s.$g.$t.fasta.tmp
                                                	charcount=$(wc $s.$g.$t.fasta.tmp | awk '{print$3}')
                                                	printf "%s\t%d\n" "$t" "$charcount" >> $s.$g.comp.tmp
                                        	done < $s.$g.N.tmp
						orth=$(sort -nk 2 $s.$g.comp.tmp | tail -n 1 | awk '{print$1}')
                                        	#the longest of these isoforms is recorded as the preferred isoform for this gene.
                                        	echo "$orth is the longest official isoform, recording ..."
					fi
				fi
			fi
		fi
		echo $orth
		#the preferred isoform for this gene is recorded
                printf "%s\t%s\n" $g $orth >> $s.IsoResolved.tsv
		#remove the temporary files within the loop to prevent unnecessary file cumilation
              	rm *tmp
                #check for duplicated isoforms within the file. This means that something has gone wrong with the filtering process/input Master table.
                dupcheck=$(awk '{print $2}' $s.IsoResolved.tsv | grep "P_" | sort -k 1 | uniq -d | wc -l)
                #if this goes wrong, the script must be stopped and the input reassessed.
                if [ $dupcheck -gt 0 ]; then
                        echo "Something has gone wrong with the script whilst assessing the $s species, gene $g. Exiting"
                        echo "Script failed to resolve $g, species $s" >> fail.txt
                        exit
                else
                        #if everything is fine, continue
                        echo "..."
                fi
                #as this gene has been completed, we add on to the count
		let count=count+1
	done < $s.gene.list
	#once the loop is finished, we have all the preferred isoforms for all the input genes for this species
        #we use this information to create new reduced fasta files
        echo "Compiling isoform resolved proteome for $s ..."
        #make a list of preferred isoforms
        awk '{print $2}' $s.IsoResolved.tsv > $s.iso.gene.list
        #extract these from this species' original transcriptome
        python ~/Scripts/get_seq2.py $3/$s*fasta $s.iso.gene.list $s.fa
        #remove the isoform list
        rm $s.ID.edit.txt
done < $1

mkdir ReducedProteomes/
mv *.fa ReducedProteomes/

rm *gene.list

echo "Script run complete."
