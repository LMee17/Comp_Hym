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
	awk '{print$1,$2}' $2/$s*tsv | sed '1d' > $s.ID.edit.txt
	#generate a list of just genes
	awk '{print $1}' $2/$s*tsv | sed '1d' | sort | uniq > $s.gene.list
	#set the counter
	count=1
	#begin reading through the  entire gene list for this species
	while read g; do
		echo "Gene $g..."
		#for each gene, create a temporary text file of the gene and corresponding transcript(s)
		grep "^$g " $s.ID.edit.txt > $s.$g.tmp
		#CHECK ONE: Check for 1 isoform
		one=$(wc -l $s.$g.tmp | awk '{print$1}')
		if [ $one -eq 1 ]; then
			#one: yes
			#record the single isoform
			echo "Only one transcript. Recording isoform ..."
			orth=$(grep $g $s.$g.tmp | awk '{print$2}')
			echo "Resolved: $g = $orth"
		else
			#one: no
			echo "More than one transcript detected. Assessing..."
			#check for NM/NR annotated refseq references
			nm=$(grep -c "NM_\|NR_" $s.$g.tmp)
			if [ $nm -eq 0 ]; then
				#no annotated transcripts
				echo "No annotated transcripts detected. Assessing..."
				#check for NC
				nc=$(grep -c "XR_" $s.$g.tmp)
				if [ $nc -eq 0 ]; then
					echo "No non-coding transcripts detected. Resolving remaing variants ..."
					#record longest isoform
					awk '{print$2}' $s.$g.tmp > $s.$g.target.tmp
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
					echo "Resolved: $g = $orth"
				else
					echo "Non-coding transcripts detected. Assessing..."
					#all noncoding ?
					if [ $nc -eq $one ]; then
                                        	#all noncoding
                                                echo "All annotated transcripts are non-coding. Resolving..."
                                                #take longest isoform
                                                #resolve remaining isoforms by length
                                                grep "XR_" $s.$g.tmp |awk '{print$2}' > $s.$g.target.tmp
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
                                                echo "Resolved: $g = $orth"
                                      	else
						#non coding and coding both detected
						echo "Both coding and non-coding transcripts detected. Removing non-coding"
						#remove nc
						grep -v "XR_" $s.$g.tmp > $s.$g.nc.tmp
                                                #is there now only one transcript left ?
                                                nc1=$(wc -l $s.$g.nc.tmp | awk '{print$1}')
                                                if [ $nc1 -eq 1 ]; then
                                                	#there is now only one transcript left. Record isoform.
                                                        echo "Only one remaining transcript variant. Recording isoform..."
                                                        orth=$(grep $g $s.$g.nc.tmp | awk '{print$2}')
                           	              	else
                      	         	                #there's still multiple variants left to resolve
                             	                        echo "More than one remaining transcript variant. Resolving..."
                                                        #record longest isform
                                                        #resolve remaining isoforms by length
                                                        awk '{print$2}' $s.$g.nc.tmp > $s.$g.target.tmp
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
                                                        echo "Resolved: $g = $orth"
                                               	fi
					fi
				fi
			else
				#annotated transcripts present
				echo "At least one annotated transcript detected..."
				#is there one or more?
				nm1=$(grep -c "NM_\|NR_" $s.$g.tmp)
				if [ $nm1 -eq 1 ]; then
					echo "Only one annotated transcript detected. Recording isoform..."
					#record single isoform
					orth=$(grep "NM_\|NR_" $s.$g.tmp | awk '{print$2}')
					echo "Resolved: $g = $orth"
				else
					echo "More than one annotated transcript detected. Resolving..."
					#check for NC
					nc2=$(grep -c "NR_" $s.$g.tmp)
					if [ $nc2 -eq 0 ]; then
						#no non coding annotated transcripts present
						echo "No annotated transcripts are non-coding. Resolving remaining coding isoforms..."
						#resolve remaining isoforms by length
						#only consider annotated transcripts
						grep "NM_" $s.$g.tmp | awk '{print$2}' > $s.$g.target.tmp
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
						echo "Resolved: $g = $orth"
					else
						#at least one transcript is non-coding
						echo "At least one annotated transcript is non-coding. Assessing ..."
						#are they all noncoding?
						if [ $nc2 -eq $one ]; then
							#all noncoding
							echo "All annotated transcripts are non-coding. Resolving..."
							#take longest isoform
							#resolve remaining isoforms by length
	                                                grep "NR_" $s.$g.tmp |awk '{print$2}' > $s.$g.target.tmp
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
                                                	echo "Resolved: $g = $orth"
						else
							#not all noncoding
                                                        echo "Both coding and non-coding transcripts detected. Removing non-coding variant..."
                                                        #remove noncoding transcript
                                                        grep -v "NR_" $s.$g.tmp > $s.$g.nc2.tmp
                                                        #is there now only one transcript left ?
                                                        nca1=$(wc -l $s.$g.nc2.tmp | awk '{print$1}')
                                                        if [ $nca1 -eq 1 ]; then
								#there is now only one transcript left. Record isoform.
								echo "Only one remaining transcript variant. Recording isoform..."
								orth=$(grep $g $s.$g.nc2.tmp | awk '{print$2}')
							else
								#there's still multiple variants left to resolve
								echo "More than one remaining transcript variant. Resolving..."
								#record longest isform
								#resolve remaining isoforms by length
		                                                grep "NM_" $s.$g.nc2.tmp | awk '{print$2}' > $s.$g.target.tmp
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
                                                		echo "Resolved: $g = $orth"
							fi
						fi
					fi
				fi
			fi
		fi
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

