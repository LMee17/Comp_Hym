# Comparing the transcriptomic response to immune challenge across 4 Hymenoptera species.

## 1: Assessing Raw Reads and Pseudoalignment

This portion of the analysis was run on remote linux servers at the University of Liverpool. All scripts mentioned in this section are found in `Scripts/Bash/`

### Step 1: Assessing read quality post transfer

The RNA in this experiment was sequenced by NovoGene UK and as such the scripts used here assume that the directory structure reflects that of data that has been received from the above service. 

To check that there had been no corruption during data transfer the script `MDcheck.sh` is used.

```sh
bash MDcheck.sh <directory containing raw reads from NovoGene UK>
```

### Step 2: Trimming reads

Once we are assured that the files have no signs of transfer issues, we can move on to trimming the reads. For this we use `Trim Galore` (a CutAdapt wrapper; Martin, 2011) and the script `TrimLoop.sh`. The script requires to be told where to find the reads and also where to write output.

```sh
bash TrimLoop.sh < directory containing raw reads from NovoGene UK> <output>
```

### Step 3: Mapping using `kallisto`

Build information for the genomes used in this project:
- _Apis mellifera_: AmelHAv3.1, GCF 003254395.2, Wallberg _et al_ (2019)
- _Bombus terrestris_: iyBomTerr1.2, GCF 910591885.1
- _Ceratina australensis_: ASM430768v1, GCA 004307685.1, Rehan _et al_ (2018)
- _Polistes canadensis_: ASM131383v1, GCF 001313835.1, Patalano _et al_ (2015)

As one of our builds (_C. australensis_) lacked adequate available annotations required for other mappers, we used `kallisto` (Bray _et al_, 2016) to pseudoalign k-mers of our reads to the transcriptome builds of the above species. This requires two steps: producing a genome index, and then running the pseudoalignment, both of which is ran using `MapKal_v3.sh`. In order for this script to run, a text file with a list of the samples (one per line) must be produced per species and batch of samples.

```sh
bash MapKal_v3.sh <batch text file> <directory containing transcriptome> \\
<directory of trimmed reads>
```

## 2: Determining Orthology

In order to more directly compare DE across the four species orthogroups must be determined. For this, we used `OrthoFinder` (Emms and Kelly, 2019)., but there were a number of pre-processing steps that were required first. Again, unless otherwise stated, the scripts mentioned in this section are found in `Scripts/Bash/`

### Step 1: Producing tables with Gene / Transcript / Protein IDs per build

We required both the protein and transcript to be reduced to one isoform/variant per gene to optimise our orthogrouping. In order to do this, we needed to use the `.gff` files of the three builds where they were available in order to parse out transcripts and proteins associated with each gene. The _C. australensis_ build had the same ID per gene / transcript / protein (i.e. gene Caust.v210001 would have variants Caust.v210001-RA and Caust.v210001-RB) which made producing ID tables considerably simpler.

Per build, `gff` and `faa` proteome files were downloaded from NCBI and the script `GFFRead_2.sh` was used to produce a `MasterIDTable.tsv` file, which can then be used to reduce isoforms per gene.

```sh
bash GFFRead_2.sh <gff> <proteome>
```

## Step 2: Isoform-Reducing Builds

There was then two slightly different methods of reducing the transcriptomes and proteomes of the four species. For proteomes the script `GeneReduce_Protv2.sh` is used. This requires the creation of another text file which lists the species that are to be reduced in the same way that the build files are prefixed, one per line. I.e. for _A. mellifera_ our proteome would be `Amel_AmelHAv3.1_GCF 003254395.2.faa` or similar, and thus the corresponding line in the species text would simply be Amel. Otherwise, the only file requirement is that the proteome ends in the file extension `.faa`.

ID tables produced in the step above are expected to be kept together in one directory, and the proteomes in the above format also kept in a directory.

```sh
bash GeneReduce_Protv2.sh <species.txt> <location of ID tables> \\
<location of proteomes>
```
This will determine one isoform per gene by:
- Checking to see if there are any annotated isoforms present (NP_). If there are, and there is only oen present, this will be assigned as the gene representative isoform.
- If there are more than one annotated isoform present, then the longest of these is designated the gene isoform.
- If there are no annotated isoforms present then the longest of the isoforms is selected as representative.
- If there is only one isoform associated with that gene, then this automatically is selected. 

For transcriptomes I use the script `GeneReduce_Trans.sh` which runs essentially the same as above, but with the added steps of only taking non-coding variants if there are no coding variants present. It also assumes the fasta extension will be `.fasta`.

Both of these scripts will then produce a reduced fasta file with the selected isoforms per gene.

### Step 3: Running `OrthoFinder`

`OrthoFinder` works best when the phylogeny is filled out, and, as such, I added four other Hymenopterans to the species to be compared, shown below:
- _Atta cephalotes_: AttaCep1.0, GCF 000143395.1, Suen _et al_ (2011)
- _Habropoda laboriosa_: Hlab1.0, GCF 001263275.1, Kapheim _et al_ (2015)
- _Solenopsis invicta_: Sinv3.0, GCA 016802725.1
- _Vespa crabro_: iyVesCrab1.2: GCF 910589235.1

For each of these I ran through the steps outlined above, producing isoform-reduced fasta files ready to feed to `OrthoFinder`.

For protein based sequences:

```sh
orthofinder -f <directory of proteomes>
```

For nucleotide based sequence files the flag -d must be included.

### Step 4: Parsing the Results

`OrthoFinder` produces a `Orthogroups.tsv` file that contains all the orthogrouping information. I downloaded this tables locally and used a R script (`Scripts/R`) to make in a more readable format, removing the superfluous species and combining both the nucleotide and protien results: `OrthoParseNCombine.R`. Note that within this script are a couple of manual additions regarding AMPs.

## 3: Differential Expression Analysis

The following scripts are found in `Scripts/R/`.

There were two levels of differential expression (DE) analyses that were ran: 1) Gene-Level and 2) OrthoGroup Level. `kallisto` produces abundance count files. These files were downloaded locally (`input/Counts_Kal`) and converted to gene-level count files using `tximport` (Soneson _et al_, 2015). From here, DE was ran using DESeq2 (Love _et al_ 2014). All of this are undergone in two scripts. One (`GeneLevelDE_NoOutRem.R`) does not remove outliers whilst the other (`GeneLevelDE.R`) does. Results are saved in respective directories within `output/Ind_DE` once ran by species - those without outliers removed are appended with "_12", for the total number of samples. 

Orthogroup level analyses were ran in much the same way using `OrthoLevelDE_NoOutRem.R` and `OrthoLevelDE.R` respectively. Results are written in `output/Ortho_DE` per grouping of species (ie Aculeata/Anthophila). There are a number of file appendages that mean different things:

- `_sA`: "sans all". All outliers are removed.
- `_sB`: "sans Bombus". Only _B. terrestris_ outliers are removed.
- `_aS`: "all samples". No outliers are removed.

Unfortunately, the Aculeata, Anthophila and Apinae directories without such a suffix as above include all samples, whereas this is the opposite for the AgeControlled grouping, which has all outliers removed. 

## 4: GO Analysis

For GO analysis, use script `GOAnalyses.R`. This script uses GO annotations determined by orthology using `eggNOG`'s sequence mapper (Cantalapiedra _et al_, 2021, Huerta-Cepas _et al_, 2019). This was done for each of the four proteomes with parameters set to default. GO analysis is ran using `topGO` (Alexa and Rahnenfuhrer, 2016). Results are saved in `output/Ind_DE/` per species. No analysis is ran on the Orthogroup level results (as there are not enough genes of interest).

## Program Versions
`DESeq2` v 1.36.0 \
`eggNOG` 2.1.7 \
`kallisto` v.0.43.0 \
`OrthoFinder`  v2.5.4 \
`topGO` 2.48.0 \
`Trim Galore`  v.0.6.6 \
`tximport` v 1.24.0

## References

Alexa, Adrian, and Jörg Rahnenführer. "Gene set enrichment analysis with topGO." _Bioconductor Improv_ 27 (2009): 1-26. \
Bray, Nicolas L., et al. "Near-optimal probabilistic RNA-seq quantification." _Nature biotechnology_ 34.5 (2016): 525-527. \
Cantalapiedra, Carlos P., et al. "eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale." _Molecular biology and evolution_ 38.12 (2021): 5825-5829. \
Emms, David M., and Steven Kelly. "OrthoFinder: phylogenetic orthology inference for comparative genomics." _Genome biology_ 20.1 (2019): 1-14. \
Huerta-Cepas, Jaime, et al. "eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses." _Nucleic acids research_ 47.D1 (2019): D309-D314. \
Kapheim, Karen M., et al. "Genomic signatures of evolutionary transitions from solitary to group living." _Science_ 348.6239 (2015): 1139-1143. \
Love, Michael I., Wolfgang Huber, and Simon Anders. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." _Genome biology_ 15.12 (2014): 1-21 \
Martin, Marcel. "Cutadapt removes adapter sequences from high-throughput sequencing reads." _EMBnet. journal_ 17.1 (2011): 10-12. \
Patalano, Solenn, et al. "Molecular signatures of plastic phenotypes in two eusocial insect species with simple societies." _Proceedings of the National Academy of Sciences_ 112.45 (2015): 13970-13975. \
Rehan, Sandra M., et al. "Conserved genes underlie phenotypic plasticity in an incipiently social bee." _Genome biology and evolution_ 10.10 (2018): 2749-2758. \
Suen, Garret, et al. "The genome sequence of the leaf-cutter ant Atta cephalotes reveals insights into its obligate symbiotic lifestyle." _PLoS genetics_ 7.2 (2011): e1002007. \
Soneson, Charlotte, Michael I. Love, and Mark D. Robinson. "Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences." _F1000Research_ 4 (2015). \
Wallberg, Andreas, et al. "A hybrid de novo genome assembly of the honeybee, Apis mellifera, with chromosome-length scaffolds." _BMC genomics_ 20.1 (2019): 1-19. 











