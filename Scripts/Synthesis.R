#5th August 2022
#Comparing the signficantly DE genes between species / orthogroupings

library("dplyr")
library(stringr)

#####All####
#read in and annotate individual results
dir.ind <- "output/Ind_DE"
targ <- list.dirs(dir.ind)
targ <- targ[-(grep("Plots", targ))]
targ <- targ[c(-1, -4, -6, -8)]

ind.list <- vector(mode = "list", length = length(targ))
for ( i in 1:length(ind.list)){
  file <- list.files(path = paste(targ[i], sep ="/"), 
                     pattern = 'DeSeqResults_FDRlt0.1.tsv')
  spec <- strsplit(targ[i], "/")[[1]][3]
  ind.list[[i]] <- read.table(paste(targ[i], file, sep = "/"),
                              header = T, sep = "\t")
  ind.list[[i]]$Source <- paste(spec)
}
ind.res <- bind_rows(ind.list)

#add orthogroup information
ortho <- read.table("input/Ortho_Files/OrthoFinder_Combined_Parsed_Aug22.tsv",
                    sep = "\t", header = T)
#remember to add apidaecin orthogroup info manually
ann <- read.table("input/Gene_Info/Amel_ImmAnnotations_Orthogroup_Jun22.tsv",
                  header = T, sep = "\t")
ann$Orthogroup[ann$GeneID == "Apid1"] <- "OG0011065"
write.table(ann, "input/Gene_Info/Amel_ImmAnnotations_Orthogroup_Aug22.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

ann.slim <- ann[,c(2:3)]

j <- 1
for(i in 1:nrow(ind.res)){
  check <- ind.res$Gene[i] %in% ortho$GeneID
  if (check == FALSE){
    species <- ind.res$Source[i]
    ind.res$OrthoGroup[i] <- paste("Unassigned", species, j, sep = "_")
    j <- j + 1
  } else {
    ind.res$OrthoGroup[i] <- ortho$Orthogroup[ortho$GeneID ==
                                                ind.res$Gene[i]]
  }
}
for(i in 1:nrow(ind.res)){
  if (grepl("Unnasigned", ind.res$OrthoGroup[i]) == TRUE) {
    ind.res$Function[i] <- "Unknown"} else {
      check2 <- ind.res$OrthoGroup[i] %in% ann.slim$Orthogroup
      if (check2 == FALSE){
        ind.res$Function[i] <- "Unknown"
      } else {
        f <- unique(ann.slim$Function[ann.slim$Orthogroup == ind.res$OrthoGroup[i]])
        if (length(f) > 1){
          if(TRUE %in% grepl("Effector|Recognition|Signalling", f)){
            ind.res$Function[i] <- "Immune"
          } else {
            if(TRUE %in% grepl("Putative Immune", f)){
              ind.res$Function[i] <- "Putative Immune"
            } 
          }
        } else {
          f2 <- gsub("Effector", "Immune", 
                     gsub("Recognition", "Immune", 
                          gsub("Signalling", "Immune", 
                               gsub("Background", "None Immune", f))))
          ind.res$Function[i] <- paste(f2)
        }
      }
    }
}


###read in the orthogrouped analyses
dir.ort <- "output/Ortho_DE"
targ <- list.dirs(dir.ort)
targ <- targ[-(grep("Filtered_Counts", targ))]
targ <- targ[c(-1, -2, -4, -6, -8, -9)]

ort.list <- vector(mode = "list", length = length(targ))
for ( i in 1:length(ort.list)){
  file <- list.files(path = paste(targ[i], sep ="/"), 
                     pattern = 'DeSeqResults_FDRlt0.1.tsv')
  spec <- strsplit(targ[i], "/")[[1]][3]
  ort.list[[i]] <- read.table(paste(targ[i], file, sep = "/"),
                              header = T, sep = "\t")
  ort.list[[i]]$Source <- paste(spec)
}
ort.res <- bind_rows(ort.list)

for(i in 1:nrow(ort.res)){
  if (ort.res$OrthoGroup[i] == "Unassigned") {
    ort.res$Function[i] <- "Unknown"} else {
      check2 <- ort.res$OrthoGroup[i] %in% ann.slim$Orthogroup
      if (check2 == FALSE){
        ort.res$Function[i] <- "Unknown"
      } else {
        f <- unique(ann.slim$Function[ann.slim$Orthogroup == ort.res$OrthoGroup[i]])
        if (length(f) > 1){
          if(TRUE %in% grepl("Effector|Recognition|Signalling", f)){
            ort.res$Function[i] <- "Immune"
          } else {
            if(TRUE %in% grepl("Putative Immune", f)){
              ort.res$Function[i] <- "Putative Immune"
            } 
          }
        } else {
          f2 <- gsub("Effector", "Immune", 
                     gsub("Recognition", "Immune", 
                          gsub("Signalling", "Immune", 
                               gsub("Background", "None Immune", f))))
          ort.res$Function[i] <- paste(f2)
        }
      }
    }
}
head(ort.res)

res <- rbind(ind.res, ort.res)
names(ind.res)

#for ind, keep only orthogroup  / function / source / contrast
#for ort, keep orthogroup/function/source /contrast
ind.res <- ind.res[,c(11,12,10,7)]
ort.res <- ort.res[,c(8,10,9,7)]

names(ort.res)

#starting with ort, aggregate by orthogroup, collapsing the Source
#rename runs and treatments
ort.res$Source[ort.res$Source == "Apinae_sB"] <- "Apinae"
ort.res$Source[ort.res$Source == "Anthophila_sA"] <- "Anthophila"
ort.res$Source[ort.res$Source == "Aculeata_sA"] <- "Aculeata"
ort.res$Contrast[grepl("Wound",ort.res$Contrast)] <- "Wound"
ort.res$Contrast[grepl("GramPos",ort.res$Contrast)] <- "GramPositive"
ort.res$Contrast[grepl("GramNeg",ort.res$Contrast)] <- "GramNegative"

ort.res <- aggregate(ort.res[3], ort.res[-3],
          FUN = function(X) paste(unique(X), collapse=", "))
#now aggregate by treatment
ort.res <- aggregate(ort.res[3], ort.res[-3],
                  FUN = function(X) paste(unique(X), collapse=", "))

#now for ind runs
#make it nicer
ind.res$Contrast[grepl("Wound", ind.res$Contrast)] <- "Wound"
ind.res$Contrast[grepl("GramPos", ind.res$Contrast)] <- "GramPositive"
ind.res$Contrast[grepl("GramNeg", ind.res$Contrast)] <- "GramNegative"
ind.res$Source[ind.res$Source == "Amel"] <- "A.melifera"
ind.res$Source[ind.res$Source == "Bter"] <- "B.terrestris"
ind.res$Source[ind.res$Source == "Caus"] <- "C.australensis"
ind.res$Source[ind.res$Source == "Plan"] <- "P.lanio"

#first, aggregate by orthogroup, collapsing source
ind.res <- aggregate(ind.res[3], ind.res[-3],
          FUN = function(X) paste(unique(X), collapse=", "))
#and now by contrast
ind.res <- aggregate(ind.res[3], ind.res[-3],
                  FUN = function(X) paste(unique(X), collapse=", "))

##combine and do again
all.res <- rbind(ind.res, ort.res)
all.res <- aggregate(all.res[3], all.res[-3],
                  FUN = function(X) paste(unique(X), collapse=", "))

all.res$DETot <- str_count(all.res$Source, pattern = ",")+1

write.table(all.res, 
            "output/Synthesis/AllSigOrthogroups_AllDE.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

