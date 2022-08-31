#!usr/bin/bash

#USAGE: bash GeneReduce.sh <species text> <ID table folder> <Transcriptome Folder>

#A script to reduce transcriptomes to longest isoform per gene
#Requires a table with Gene and corresponding TranscriptIDs in two different columns (GENE \t TRANSCRIPT)
#Requires the transcriptome that is to be reduced
#Requires a list of the species of transcriptomes to be reduced
#Species should be one per line
#Transcriptomes should begin with species (as appears on line of Species.txt file)
#I.e. Amel as an input should correspond with Amel.fasta or Amel_GCF.fasta.
#Transcriptomes must end in fasta
#Produces an output folder of reduced transcriptomes.

#Begins by reading through the species input file
while read s; do
	#each species is assessed
	echo "Assessing $s...."
	#take only necessary information from input tables (GENE \t TRANSCRIPT \t PROTEIN)
	awk '{print$1,$2}' $2/$s*tsv > $s.ID.edit.txt
	#generate a list of just genes
	awk '{print $1}' $2/$s*tsv | sort | uniq > $s.gene.list
	#set the counter
	count=1
	#begin reading through the  entire gene list for this species
	while read g; do
		echo "Gene $g..."
		#for each gene, create a temporary text file of the gene and corresponding transcript(s)
		grep "^$g " $s.ID.edit.txt > $s.$g.tmp
		#count the number of transcripts
		check=$(wc -l $s.$g.tmp | awk '{print$1}')
		#if there is just one transcript, then this can be recorded as the isoform for this gene and we move to the next gene
		if [ $check -eq 1 ]; then
			orth=$(grep "^$g " $s.ID.edit.txt | awk '{print $2}')
			echo "Only one isoform, $orth, recording ..."
		#if there are more than one isoforms for this gene, we begin to resolve
		else
			echo "$g has more than one isoform. Resolving"
			#check 1: NM
			#NM are the best annotated isoforms using NCBI transcriptomes. If there are any NM_ transcript(s) these should be assessed first.
			nmcheck=$(grep -c "NM_\|NR_" $s.$g.tmp)
			#if there are no NM_ isoforms, that we move to taking the longest isoform
			if [ $nmcheck -eq 0 ]; then
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
				if [ $nmcheck -eq 1 ]; then
					#if there is only one, this is recorded as the preferred isoform.
					orth=$(grep "NM\|NR_" $s.$g.tmp | awk '{print$2}')
					echo "Official isoform present, $orth, recording ..."
				#if there is more than one, we repeat the length comparison above using just the NM_ isoforms
				else
					echo "More than one official isoform present, assessing by isoform length..."
					grep "NM_\|NR_" $s.$g.tmp | awk '{print$2}' > $s.$g.N.tmp
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
		#the preferred isoform for this gene is recorded
		printf "%s\t%s\n" $g $orth >> $s.IsoResolved.tsv
		#remove the temporary files within the loop to prevent unnecessary file cumilation
		rm *tmp
		#check for duplicated isoforms within the file. This means that something has gone wrong with the filtering process/input Master table.
		dupcheck=$(awk '{print $2}' $s.IsoResolved.tsv | sort -k 1 | uniq -d | wc -l)
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
	echo "Compiling isoform resolved transcriptome for $s ..."
	#make a list of preferred isoforms
	awk '{print $2}' $s.IsoResolved.tsv > $s.iso.gene.list
	#extract these from this species' original transcriptome
	python ~/Scripts/get_seq2.py $3/$s*fasta $s.iso.gene.list $s.fa
	#remove the isoform list
	rm $s.ID.edit.txt
done < $1

#once this is done for all species input, we compile the new reduced fastas into one place
mkdir ReducedTranscriptomes/
mv *.fa ReducedTranscriptomes/

rm *gene.list

echo "Script run complete."
