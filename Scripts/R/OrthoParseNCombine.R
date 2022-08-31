#June 2022
#A script that 1) parses the OrthoFinder result tables (one per protein, one per transcript)
#into a more readable format. Then, these are combined, using the protein results as a base
#but adding the groupings of non-coding RNAs that were determined in the transcript 
#OrthoFinder run. Any manual additions are added at the end (i.e., adding in the apideacin
#grouping)

library("dplyr")
library("stringr")
library("tidyr")

#Parse one: Protein
#End table results in Orthogroup / Gene / Species
#Assumed orthofinder has been run using 8 Hymenopteran species (Acep/Amel/Bter/Caus/Hlab/Pcan/Sinv/Vcra)

df <- read.table("input/Ortho_Files/OrthoFinder_Orthogroups_Prot_Jun22.tsv",
                 sep = "\t", header = T, quote = "")
names(df)

#remove extraneous species
df <- df[,c(1,3:5,7)]

orth.list <- vector(mode = "list", length = 4)

for (i in 1:length(orth.list)){
  orth.list[[i]] <- df[,c(1,(i+1))]
  species <- names(orth.list[[i]][2])
  orth.list[[i]]$Species <- paste(species)
  names(orth.list[[i]]) <- c("Orthogroup", "Description", "Species")
  orth.list[[i]] <- orth.list[[i]][!orth.list[[i]]$Description == "",]
  orth.list[[i]] <- as.data.frame(separate_rows(orth.list[[i]], Description, sep= ","))
}

ortho <- bind_rows(orth.list)
#remove random whitespaces
ortho$Description <- trimws(ortho$Description)
ortho$ProteinID <- sapply(strsplit(ortho$Description, " "), '[', 1)

#Add Gene info
#load the tables for each species
amel.id <- read.table("input/Gene_Info/Amel_HAv3.1_MasterIDtable.tsv", sep = "\t", header = T)
bter.id <- read.table("input/Ortho_Files/Bter.MasterIDtable.tsv", sep = "\t", header = T)
pcan.id <- read.table("input/Ortho_Files/Pcan_MasterIDtable.tsv", sep = "\t", header = T)
#caus needs making
caus.id <- read.table("input/Ortho_Files/Caus.ID.txt", header = F, sep = "\t")
names(caus.id) <- c("GeneID", "TranscriptID")
#protein ids are the same as the transcripts
caus.id$ProteinID <- caus.id$TranscriptID

master.id <- rbind(amel.id, bter.id, pcan.id, caus.id)

#let's get the genes
genes <- vector(length = nrow(ortho))
for (i in 1:length(genes)){
  genes[i] <- master.id$GeneID[master.id$ProteinID == ortho[i,4]]
}

ortho <- as.data.frame(ortho)
ortho$GeneID <- genes

#writeup
write.table(ortho, "input/Ortho_Files/OrthoFinder_Prot_Parsed_Jun22.tsv", 
            col.names = T, row.names = F, sep = "\t", quote =F)
ortho.prot <- ortho

##Repeat for cds

df <- read.table("input/Ortho_Files/OrthoFinder_Orthogroups_CDS_Jun22.tsv",
                 sep = "\t", header = T, quote = "")

#remove extraneous species
df <- df[,c(1,3:5,7)]

orth.list <- vector(mode = "list", length = 4)

for (i in 1:length(orth.list)){
  orth.list[[i]] <- df[,c(1,(i+1))]
  species <- names(orth.list[[i]][2])
  orth.list[[i]]$Species <- paste(species)
  names(orth.list[[i]]) <- c("Orthogroup", "TranscriptID", "Species")
  orth.list[[i]] <- orth.list[[i]][!orth.list[[i]]$TranscriptID == "",]
  orth.list[[i]] <- as.data.frame(separate_rows(orth.list[[i]], TranscriptID, sep= ","))
}


ortho <- bind_rows(orth.list)

#remove random whitespaces
ortho$TranscriptID <- trimws(ortho$TranscriptID)

#Add Gene info
#load the tables for each species
amel.id <- read.table("input/Gene_Info/Amel_HAv3.1_MasterIDtable.tsv", sep = "\t", header = T)
bter.id <- read.table("input/Ortho_Files/Bter.MasterIDtable.tsv", sep = "\t", header = T)
pcan.id <- read.table("input/Ortho_Files/Pcan_MasterIDtable.tsv", sep = "\t", header = T)
#caus needs making
caus.id <- read.table("input/Ortho_Files/Caus.ID.txt", header = F, sep = "\t")
names(caus.id) <- c("GeneID", "TranscriptID")
#protein ids are the same as the transcripts
caus.id$ProteinID <- caus.id$TranscriptID

master.id <- rbind(amel.id, bter.id, pcan.id, caus.id)

#let's get the genes
genes <- vector(length = nrow(ortho))
for (i in 1:length(genes)){
  genes[i] <- master.id$GeneID[master.id$TranscriptID == ortho[i,2]]
}

ortho <- as.data.frame(ortho)
ortho$GeneID <- genes

#writeup
write.table(ortho, "input/Ortho_Files/OrthoFinder_CDS_Parsed_Jun22.tsv", 
            col.names = T, row.names = F, sep = "\t", quote =F)
ortho.cds <- ortho

##Combine the two
#slim down prot ortho to match CDS
ortho.prot2 <- ortho.prot[,c(1,4,3,5)]
head(ortho.prot2)

##Non-coding genes
#The only groupings from the cds run I'm interested in are those that consist of 
#non-coding genes. Each group that contains only non-coding genes will be given an 
#orthogroup ID started with the last protein orthogroup ID number + 1 (ie OG00100 
#becomes OG00101)
#Assumption: Caus transcripts that are in predominantly non-coding orthogroups are 
#also assumed to be non-coding.
#Gonna start with a logic matrix with a 1 for presence of non-coding, a 
#1 for presence of coding in separate columns. 

logmat <- as.data.frame(unique(ortho.cds$Orthogroup))
names(logmat) <- "Orthogroup"

#Functions to assess presence of coding / non-coding transcripts within orthogroups.

ncScan <- function(x){
  trans <- ortho.cds$TranscriptID[ortho.cds$Orthogroup == paste(x)]
  verd <- unique(grepl("R_", trans))
  if ( length(verd) == 1 ){
    ifelse(verd == TRUE, nc <- 1, nc <- 0)
  } else {
    nc <- 1
  }
  return(nc)
}

#apply
logmat$NonCod <- sapply(logmat$Orthogroup, ncScan)

codScan <- function(x){
  trans <- ortho.cds$TranscriptID[ortho.cds$Orthogroup == paste(x)]
  verd <- unique(grepl("M_", trans))
  if ( length(verd) == 1 ){
    ifelse(verd == TRUE, cod <- 1, cod <- 0)
  } else {
    cod <- 1
  }
  return(cod)
}

#apply
logmat$Cod <- sapply(logmat$Orthogroup, codScan)

#Extract groups that are all non-coding.

all.nc <- logmat$Orthogroup[logmat$NonCod == 1 &
                              logmat$Cod == 0]

ortho.nc <- ortho.cds[ortho.cds$Orthogroup  %in% all.nc,]

#Get the last protein ID

last.id <- tail(ortho.prot, n = 1)

nc.res <- as.data.frame(unique(ortho.nc$Orthogroup))
names(nc.res) <- "CDS_Orthogroup"

res.id <- vector(length = nrow(nc.res))

for (i in 1:length(res.id)){
  no <- last.id + i
  res.id[i] <- paste("OG00", no, sep = "")
}

nc.res$ResolvedOrthogroup <- res.id

####Combine
#Reduce the protein dataframe to just gene, species and orthogroup

ortho <- ortho.prot2[,c(1,3:4)]

#massagethe non-coding data so that it can be appended

for (i in 1:nrow(ortho.nc)){
  ortho.nc[i,5] <- nc.res$ResolvedOrthogroup[nc.res$CDS_Orthogroup == ortho.nc[i,1]]
}

toadd <- ortho.nc[,c(5,3,4)]
names(toadd)[1] <- "Orthogroup"

#The Caus Problem: There are cases were Caus sequences have been sorted into groups 
#with NR in the transcript run, and also into a protein orthogroup in the protein run. 
#These have to be assessed.

caus.nc <- toadd$GeneID[grep("Caust", toadd$GeneID)]

ortho.prot2[ortho.prot2$GeneID %in% caus.nc,]

#I went through these occurrences manually
#Caust.v2_009161, Caust.v2_014029, Caust.v2_007962, Caust.v2_017529,
#Caust.v2_007433, Caust.v2_006560, Caust.v2_000693 and Caust.v2_008865 look to be proteins. 
#Any orthogroups with these present in the cds groups will be removed.

protsinnc <- ortho.nc$Orthogroup[ortho.nc$GeneID == "Caust.v2_009161" 
                                 | ortho.nc$GeneID == "Caust.v2_014029" | 
                                   ortho.nc$GeneID == "Caust.v2_007962" | ortho.nc$GeneID == "Caust.v2_017529" |
                                   ortho.nc$GeneID == "Caust.v2_007433" | ortho.nc$GeneID =="Caust.v2_006560" |
                                   ortho.nc$GeneID == "Caust.v2_000693" | ortho.nc$GeneID == "Caust.v2_008865"]

ortho.nc <- ortho.nc[!ortho.nc$Orthogroup %in% protsinnc,]

torem <- c("Caust.v2_009161", "Caust.v2_014029", "Caust.v2_007962", "Caust.v2_017529",
           "Caust.v2_007433", "Caust.v2_006560", "Caust.v2_000693", "Caust.v2_008865")
toadd <- toadd[!toadd$GeneID %in% torem,]

#Caust.v2_012543, Caust.v2_017732, Caust.v2_019456, Caust.v2_000590, 
#Caust.v2_018252, Caust.v2_021656, Caust.v2_006561 all look to be non-coding. 

ncinprot <- ortho$Orthogroup[ortho$GeneID == "Caust.v2_012543" | 
                               ortho$GeneID == "Caust.v2_017732" |
                               ortho$GeneID == "Caust.v2_019456" | ortho$GeneID == "Caust.v2_000590" |
                               ortho$GeneID == "Caust.v2_018252" | ortho$GeneID == "Caust.v2_021656" |
                               ortho$GeneID == "Caust.v2_006561"]

#All of these groupings in the prot run only contain Caus so I can remove them

ortho <- ortho[!ortho$Orthogroup %in% ncinprot,]

#Caust.v2_000414, Caust.v2_008865, Caust.v2_010854, Caust.v2_020433, 
#Caust.v2_006418, Caust.v2_005154, Caust.v2_000151 all have no similarity 
#in either protein or NC blasts to the other bees and wasp.

unknownNC <- ortho.nc$Orthogroup[ortho.nc$GeneID == "Caust.v2_000414" | 
                                   ortho.nc$GeneID == "Caust.v2_010854" | ortho.nc$GeneID == "Caust.v2_020433" |
                                   ortho.nc$GeneID == "Caust.v2_006418" | ortho.nc$GeneID == "Caust.v2_005154" |
                                   ortho.nc$GeneID ==  "Caust.v2_000151"]

ortho.nc[ortho.nc$Orthogroup %in% unknownNC,] %>%
  arrange(Orthogroup)

#There's definitely overlap in the NC runs

unknownProt <- ortho$Orthogroup[ortho$GeneID == "Caust.v2_000414" | 
                                  ortho$GeneID == "Caust.v2_010854" | ortho$GeneID == "Caust.v2_020433" |
                                  ortho$GeneID == "Caust.v2_006418" | ortho$GeneID == "Caust.v2_005154" |
                                  ortho$GeneID ==  "Caust.v2_000151"]

ortho[ortho$Orthogroup %in% unknownProt,] %>%
  arrange(Orthogroup)

#In the protein runs these unknowns are on their own, 
#so I'll leave them in the non-coding groupings

ortho <- ortho[!ortho$Orthogroup %in% unknownProt,]

#Combine the altered orthogroupings
ortho <- rbind(ortho, toadd)
ortho <- ortho[!duplicated(ortho),]

#Write up 
write.table(ortho, "input/Ortho_Files/OrthoFinder_Combined_Parsed_Jun22.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

#manually add Apidaecin grouping for Bter/ Amel
ortho[40557,] <- c("OG0011065", "Bter", "LOC100649867")
ortho[40558,] <- c("OG0011065", "Amel", "Apid1")

#write up (Aug22)
write.table(ortho, "input/Ortho_Files/OrthoFinder_Combined_Parsed_Aug22.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)



