##6th August
#Comparing the signficantly DE genes between species / orthogroupings

library("dplyr")
library(stringr)
library("ggplot2")
library("ggrepel")

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

#give the orthogroups a description, if possible
#using amel as an anchor species (and as its in all combinations)

for (i in 1:nrow(ort.res)){
  print(i)
  one <- head(ortho$GeneID[ortho$Orthogroup == ort.res$OrthoGroup[i] &
                        ortho$Species == "Amel"], n= 1)
  ort.res$Gene[i] <- one
  two <- unique(ind.res$Desc[ind.res$Gene == one])
  if (length(two) > 0){
    ort.res$Desc[i] <- two
  } else {
    ort.res$Desc[i] <- "NA"
  }
}

all.res <- rbind(ort.res, ind.res)
tail(all.res)
all.res$Direction <- ifelse(all.res$log2FoldChange > 0, "Up", "Down")
tail(all.res)

ggplot(all.res, aes(x = Direction)) +
       geom_bar() +
  facet_grid(Contrast ~ Source)


###hmm. Ind only

ind.res$Contrast[grepl("Wound", ind.res$Contrast)] <- "Wound"
ind.res$Contrast[grepl("GramPos", ind.res$Contrast)] <- "GramPositive"
ind.res$Contrast[grepl("GramNeg", ind.res$Contrast)] <- "GramNegative"
ind.res$Source[ind.res$Source == "A.melifera"] <- "A.mellifera"
ind.res$Source[ind.res$Source == "Bter"] <- "B.terrestris"
ind.res$Source[ind.res$Source == "Caus"] <- "C.australensis"
ind.res$Source[ind.res$Source == "Plan"] <- "P.lanio"

ind.res$Direction <- ifelse(ind.res$log2FoldChange > 0, "Up", "Down")

ggplot(ind.res, aes(x = Direction)) +
  geom_bar() +
  facet_grid(Contrast ~Source)

ind.res$Contrast <- factor(ind.res$Contrast, levels = c("Wound",
                                                        "GramPositive",
                                                        "GramNegative"))
ind.res$Direction <- factor(ind.res$Direction, levels = c("Up", "Down"))

ggplot(ind.res, aes(x = Direction)) +
  geom_bar(aes(fill = Direction)) +
  scale_fill_manual(values = c("red", "blue")) +
  xlab("Direction of Expression Regulation") +
  ylab("Number of Differentially Expressed Genes (FDR < 0.1)") +
  theme(strip.text.x = element_text(face = "italic")) +
  facet_grid(Contrast ~Source, scales = "free_y")

ggsave("output/Synthesis/AllSig_UpDwnReg_BySpecies.pdf")

#big volcano plot
dir.ind <- "output/Ind_DE"
targ <- list.dirs(dir.ind)
targ <- targ[-(grep("Plots", targ))]
targ <- targ[c(-1, -4, -6, -8)]

ind.list2 <- vector(mode = "list", length = length(targ))
for ( i in 1:length(ind.list2)){
  file <- list.files(path = paste(targ[i], sep ="/"), 
                     pattern = 'DeSeqResults_All.tsv')
  spec <- strsplit(targ[i], "/")[[1]][3]
  ind.list2[[i]] <- read.table(paste(targ[i], file, sep = "/"),
                              header = T, sep = "\t", quote = "")
  ind.list2[[i]]$Source <- paste(spec)
}
ind.res2 <- bind_rows(ind.list2)
j <- 1
for(i in 1:nrow(ind.res2)){
  check <- ind.res2$Gene[i] %in% ortho$GeneID
  if (check == FALSE){
    species <- ind.res2$Source[i]
    ind.res2$OrthoGroup[i] <- paste("Unassigned", species, j, sep = "_")
    j <- j + 1
  } else {
    ind.res2$OrthoGroup[i] <- ortho$Orthogroup[ortho$GeneID ==
                                                ind.res2$Gene[i]]
  }
}
for(i in 1:nrow(ind.res2)){
  if (grepl("Unnasigned", ind.res2$OrthoGroup[i]) == TRUE) {
    ind.res2$Function[i] <- "Unknown"} else {
      check2 <- ind.res2$OrthoGroup[i] %in% ann.slim$Orthogroup
      if (check2 == FALSE){
        ind.res2$Function[i] <- "Unknown"
      } else {
        f <- unique(ann.slim$Function[ann.slim$Orthogroup == ind.res2$OrthoGroup[i]])
        if (length(f) > 1){
          if(TRUE %in% grepl("Effector|Recognition|Signalling", f)){
            ind.res2$Function[i] <- "Immune"
          } else {
            if(TRUE %in% grepl("Putative Immune", f)){
              ind.res2$Function[i] <- "Putative Immune"
            } 
          }
        } else {
          f2 <- gsub("Effector", "Immune", 
                     gsub("Recognition", "Immune", 
                          gsub("Signalling", "Immune", 
                               gsub("Background", "None Immune", f))))
          ind.res2$Function[i] <- paste(f2)
        }
      }
    }
}

ind.res2$Significant <- "No"
ind.res2$Significant[ind.res2$padj < 0.1] <- "Yes"

ind.res2$Direction <- "Not Significant"
ind.res2$Direction[ind.res2$Sig == "Yes" &
              ind.res2$log2FoldChange > 0 ] <- "Up-Regulated"
ind.res2$Direction[ind.res2$Sig == "Yes" &
              ind.res2$log2FoldChange < 0 ] <- "Down-Regulated"
ind.res2$Neglog10adjPvalue <- -log10(ind.res2$padj)
head(ind.res2)


ind.res2$Contrast[grepl("Wound", ind.res2$Contrast)] <- "Wound"
ind.res2$Contrast[grepl("GramPos", ind.res2$Contrast)] <- "Gram Positive"
ind.res2$Contrast[grepl("GramNeg", ind.res2$Contrast)] <- "Gram Negative"

backup <- ind.res2

ind.res2$Source[ind.res2$Source == "A.melifera"] <- "A.mellifera"
ind.res2$Source[ind.res2$Source == "Bter"] <- "B.terrestris"
ind.res2$Source[ind.res2$Source == "Caus"] <- "C.australensis"
ind.res2$Source[ind.res2$Source == "Plan"] <- "P.lanio"

ind.res2$Contrast <- factor(ind.res2$Contrast, levels = c("Wound",
                                                      "Gram Positive",
                                                      "Gram Negative"))

ggplot(ind.res2[!is.na(ind.res2$padj),],
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.1", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(Contrast~Source, scales = "free") 

#label the top 3 most signficant gene per wound/species grid
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                ind.res2$Gene == "LOC724654"] <- "Cytochrome b5"
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                ind.res2$Gene == "LOC102654628"] <- "Uncharacterised LOC102654628"
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                ind.res2$Gene == "LOC113218757"] <- "LOC113218757"

ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC726750"] <- "Fibroin heavy chain"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC411012"] <- "Cactus2"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC724654"] <- "Cytochrome b5"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC726750"] <- "Fibroin heavy chain"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC411012"] <- "Cytochrome b5"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC411807"] <- "Uncharacterised LOC411807"

ind.res2[ind.res2$Source == "B.terrestris" &
           ind.res2$Contrast == "Gram Negative",] %>%
  arrange(padj) %>%
  head(n = 3)

ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "LOC100631061"] <- "Hymenoptaecin"
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "LOC125385347"] <- "Hymenoptaecin-like"
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "LOC100649281"] <- "Maltase A1-like"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC100631061"] <- "Hymenoptaecin"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC100644101"] <- "Uncharacterised LOC100644101"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC100649867"] <- "Apidaecins type 73"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC100631061"] <- "Hymenoptaecin"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC125385347"] <- "Hymenoptaecin-like"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC100649867"] <- "Apidaecins type 73"

ind.res2[ind.res2$Source == "C.australensis" &
           ind.res2$Contrast == "Gram Negative",] %>%
  arrange(padj) %>%
  head(n = 3)

ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "Caust.v2_016347"] <- "ATP-citrate synthase"
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "Caust.v2_020967"] <- "piwi-like protein Siwi"
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "Caust.v2_011640"] <- "TIM50-C-like"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "Caust.v2_018203"] <- "LOC108623428"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "Caust.v2_011034"] <- "Uncharacterised LOC108629120"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "Caust.v2_011210"] <- "LOC108632239"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "Caust.v2_011640"] <- "TIM50-C-like"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "Caust.v2_011034"] <- "Uncharacterised LOC108629120"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "Caust.v2_011550"] <- "Uncharacterised LOC108623077"


ind.res2[ind.res2$Source == "P.lanio" &
           ind.res2$Contrast == "Gram Negative",] %>%
  arrange(padj) %>%
  head(n = 3)


ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "LOC106790552"] <- "Uncharacterised LOC106790552"
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "LOC106787312"] <- "Tyrosine decarboxylase-like"
ind.res2$Label[ind.res2$Contrast == "Wound" & 
                 ind.res2$Gene == "LOC106791417"] <- "Protein croquemort-like"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC106787312"] <- "Tyrosine decarboxylase-like"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC106793670"] <- "LOC106793670"
ind.res2$Label[ind.res2$Contrast == "Gram Positive" & 
                 ind.res2$Gene == "LOC106790552"] <- "Uncharacterised LOC106790552"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC106785761"] <- "Protein L(2)EFL"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC106786250"] <- "Protein bicaudal C"
ind.res2$Label[ind.res2$Contrast == "Gram Negative" & 
                 ind.res2$Gene == "LOC106786125"] <- "Inhibin beta chain"

ggplot(ind.res2[!is.na(ind.res2$padj),],
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction), size = .8) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.1", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(Source~Contrast, scales = "free") +
  theme(strip.text.y = element_text(face = "italic")) 

ggsave("output/Synthesis/Volcano_allIndRuns.pdf")


length(unique(ind.res2$Gene[ind.res2$Source == "P.lanio" & 
                              ind.res2$Contrast == "Gram Negative" &
                              ind.res2$padj < 0.1 &
                              ind.res2$log2FoldChange < 0]))

