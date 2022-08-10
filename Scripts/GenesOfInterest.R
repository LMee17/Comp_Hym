#8th August 2022
#Looking at genes that occur across species

library("dplyr")
library(ggvenn)
library("ggVennDiagram")
library("ggplot2")
library("stringr")

#I'm sick of all these little scripts making an amalgamated version of the individual
#species runs so I'm going to finally save the fucker.

dir.ind <- "output/Ind_DE"
targ <- list.dirs(dir.ind)
targ <- targ[-(grep("Plots", targ))]
#remove unwanted runs = all samples for Bter/Caus/Plan
targ <- targ[c(-1,-4,-6,-8)]

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
                    header = T, sep = "\t")
ann <- read.table("input/Gene_Info/Amel_ImmAnnotations_Orthogroup_Aug22.tsv",
                  header = T, sep = "\t")
ann.slim <- ann[,c(2:3)]

#add immune annotation
for(i in 1:nrow(ind.res)){
  check <- ind.res$Gene[i] %in% ortho$GeneID
  if (check == FALSE){
    ind.res$OrthoGroup[i] <- "Unassigned"
  } else {
    ind.res$OrthoGroup[i] <- ortho$Orthogroup[ortho$GeneID ==
                                                ind.res$Gene[i]]
  }
}
for(i in 1:nrow(ind.res)){
  if (ind.res$OrthoGroup[i] == "Unassigned") {
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
head(ind.res)

write.table(ind.res, "output/Synthesis/AllSpecies_GLDE_Combined_FDRlt1.0.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

#THERE. 

#3 genes that appear as DEG in all 3 species
#1: moduler serine protease: OG0000138
ind.res[ind.res$OrthoGroup == "OG0000138",]

#2: uncharacterised "serine rich" protein
ind.res[ind.res$OrthoGroup == "OG0007622",]

#3: protein lethal(2) essential for life (L(2)EFL)
ind.res[ind.res$OrthoGroup == "OG0000067",]

#15 bee gene overlap
#1: uncharacterised 4238
ortho[ortho$Orthogroup == "OG0004238",]
ind.res[ind.res$OrthoGroup == "OG0004238",]

#2: cactus
ind.res[ind.res$OrthoGroup == "OG0000699",]

#3: serpin-5
ind.res[ind.res$OrthoGroup == "OG0004983",]

#4: malvolio
ind.res[ind.res$OrthoGroup == "OG0000627",]

#5: pug
ind.res[ind.res$OrthoGroup == "OG0002284",]

#6: cytochromeb5
ind.res[ind.res$OrthoGroup == "OG0008945",]

#7: fibroin heavy chain
ind.res[ind.res$OrthoGroup == "OG0008153",]

#8: uncharacterised 2961
ind.res[ind.res$OrthoGroup == "OG0002961",]

#9: MFS-type transporter SLC18B1
ind.res[ind.res$OrthoGroup == "OG0003250",]

#10: four-domain protease inhibitor
ind.res[ind.res$OrthoGroup == "OG0001255",]

#11: B-gluc2
ind.res[ind.res$OrthoGroup == "OG0008514",]

#12: hexamerin
ind.res[ind.res$OrthoGroup == "OG0000151",]

#13: solute carrier family 22 member 21
ind.res[ind.res$OrthoGroup == "OG0001416",]

#14: henna
ind.res[ind.res$OrthoGroup == "OG0002814",]

#15: 3-phosphoinositide-dependent protein kinase 1
ind.res[ind.res$OrthoGroup == "OG0001290",]

##Looking at overlap - venn diagram

j <- 1
for (i in 1:nrow(ind.res)){
  if (ind.res$OrthoGroup[i] == "Unassigned"){
    print(paste(i, "Yup", sep = "="))
    species <- ind.res$Source[i]
    ind.res$OrthoGroup[i] <- paste("Unnassigned", species, j, sep = "_")
    j <- j + 1
  } 
}

#all genes, all treatments
all.venn <- list("A.mellifera" = ind.res$OrthoGroup[ind.res$Source == "Amel"],
                 "B.terrestris" = ind.res$OrthoGroup[ind.res$Source == "Bter"],
                 "C.australensis" = ind.res$OrthoGroup[ind.res$Source == "Caus"],
                 "P.lanio" = ind.res$OrthoGroup[ind.res$Source == "Plan"])
venn <- Venn(all.venn)
data <- process_data(venn)

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = id), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_manual(values = c("#ffa500", "#805380", "#558c8c", "#806980", "#aa3773",
                                 "#80d253", "#aa8c55", "#ff532d", "#005aff", "#0080d2",
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a")) +
  labs(title = "All") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("output/Synthesis/All_Ind_Overlap_Venn.pdf")

##Heatmaps of all overlappy genes 
#ok need all results for this

dir.ind <- "output/Ind_DE"
targ <- list.dirs(dir.ind)
targ <- targ[-(grep("Plots", targ))]
#remove unwanted runs = all samples for Bter/Caus/Plan
targ <- targ[c(-1,-4,-6,-8)]

ind.list <- vector(mode = "list", length = length(targ))
for ( i in 1:length(ind.list)){
  file <- list.files(path = paste(targ[i], sep ="/"), 
                     pattern = 'DeSeqResults_All.tsv')
  spec <- strsplit(targ[i], "/")[[1]][3]
  ind.list[[i]] <- read.table(paste(targ[i], file, sep = "/"),
                              header = T, sep = "\t", quote = "")
  ind.list[[i]]$Source <- paste(spec)
}
ind.all.res <- bind_rows(ind.list)

for(i in 1:nrow(ind.all.res)){
  check <- ind.all.res$Gene[i] %in% ortho$GeneID
  if (check == FALSE){
    ind.all.res$OrthoGroup[i] <- "Unassigned"
  } else {
    ind.all.res$OrthoGroup[i] <- ortho$Orthogroup[ortho$GeneID ==
                                                    ind.all.res$Gene[i]]
  }
}
for(i in 1:nrow(ind.all.res)){
  if (ind.all.res$OrthoGroup[i] == "Unassigned") {
    ind.all.res$Function[i] <- "Unknown"} else {
      check2 <- ind.all.res$OrthoGroup[i] %in% ann.slim$Orthogroup
      if (check2 == FALSE){
        ind.all.res$Function[i] <- "Unknown"
      } else {
        f <- unique(ann.slim$Function[ann.slim$Orthogroup == ind.all.res$OrthoGroup[i]])
        if (length(f) > 1){
          if(TRUE %in% grepl("Effector|Recognition|Signalling", f)){
            ind.all.res$Function[i] <- "Immune"
          } else {
            if(TRUE %in% grepl("Putative Immune", f)){
              ind.all.res$Function[i] <- "Putative Immune"
            } 
          }
        } else {
          f2 <- gsub("Effector", "Immune", 
                     gsub("Recognition", "Immune", 
                          gsub("Signalling", "Immune", 
                               gsub("Background", "None Immune", f))))
          ind.all.res$Function[i] <- paste(f2)
        }
      }
    }
}
head(ind.all.res)

write.table(ind.all.res, "output/Synthesis/AllInd_AllResultsDESeq.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

#pull out only the fifteen orthogroups
fift <- c("OG0004238",
          "OG0004238",
          "OG0000699",
          "OG0004983",
          "OG0000627",
          "OG0002284",
          "OG0008945",
          "OG0008153",
          "OG0002961",
          "OG0003250",
          "OG0001255",
          "OG0008514",
          "OG0000151",
          "OG0001416",
          "OG0002814",
          "OG0001290")

fift.res <- ind.all.res[ind.all.res$OrthoGroup %in% fift,]
fift.res <- fift.res[!fift.res$Source == "Plan",]
fift.res %>%
  arrange(OrthoGroup)
for (i in 1:length(unique(fift.res$OrthoGroup))){
  out <- length(unique(fift.res$Gene[fift.res$OrthoGroup == unique(fift.res$OrthoGroup)[i]]))
  print(out)
}

#there are cases where there are unequal numbers of genes in the orthogroups that need
#to be sorted
#first, the hexamerins
#HexA in Amel = LOC100650745 in Bter = Caust.v2_000045 in Caus
#HexB in Amel = LOC100650380 in Bter = Caust.v2_000046 in Caus
#HexC in Amel = LOC100650863 in Bter = Caust.v2_001442 in Caus
fift.res$Gene[fift.res$Gene == "LOC100650745"] <- "Hex70a"
fift.res$Gene[fift.res$Gene == "Caust.v2_000045"] <- 'Hex70a'
fift.res$Gene[fift.res$Gene == "LOC100650380"] <- "HEX70b"
fift.res$Gene[fift.res$Gene == "Caust.v2_000046"] <- "HEX70b"
fift.res$Gene[fift.res$Gene == "LOC100650863"] <- "Hex70c"
fift.res$Gene[fift.res$Gene == "Caust.v2_001442"] <- "Hex70c"
#cactus1 is the only ortholog present in all 3 species
fift.res$Gene[fift.res$Gene == "cact1"] <- "cactus"
fift.res$Gene[fift.res$Gene == "cact1"] <- "cactus"

#the rest will require a key
key <- unique(fift.res$OrthoGroup)
key <- key[-1]
key <- key[-1]
key2 <- c("uncharacterised OG0004238",
          "SCF22M21",
          "henna",
          "3-PIDPK1",
          "SLC18B1",
          "serpin-5/88Ea",
          "malvolio",
          "pug",
          "4-domain protease inhibitor",
          "cytochrome b5",
          "gnbp-1",
          "fibroin heavy chain",
          "uncharacterised OG0002961")
key <- as.data.frame(cbind(key,key2))
key

for (i in 1:nrow(fift.res)){
  if(fift.res$OrthoGroup[i] %in% key$key){
    fift.res$Gene[i] <- key$key2[key$key == fift.res$OrthoGroup[i]]
  }
}
fift.res %>%
  arrange(Gene)

#remove extra cacti
fift.res2 <- fift.res[!fift.res$Gene == "LOC411012" ,] 
fift.res2 <- fift.res2[!fift.res2$Gene == "LOC725766" ,] 

fift.res2[fift.res2$OrthoGroup == "OG0000699",]

fift.res2$Contrast[grepl("Wound", fift.res2$Contrast)] <- "Wound"
fift.res2$Contrast[grepl("GramPos", fift.res2$Contrast)] <- "Gram Positive"
fift.res2$Contrast[grepl("GramNeg", fift.res2$Contrast)] <- "Gram Negative"

fift.res2$Contrast <- factor(fift.res2$Contrast, levels = c("Wound",
                                                            "Gram Positive",
                                                            "Gram Negative"))

fift.res2$Source[fift.res2$Source == "Amel"] <- "A.mellifera"
fift.res2$Source[fift.res2$Source == "Bter"] <- "B.terrestris"
fift.res2$Source[fift.res2$Source == "Caus"] <- "C.autralensis"

unique(fift.res2$Gene[fift.res2$Function == "None Immune"])

#hex70b isn't actually sig in any of the runs so
fift.res2 <- fift.res2[!fift.res2$Gene == "HEX70b",]

fift.res2$Gene <- factor(fift.res2$Gene, levels = c("uncharacterised OG0004238",
                                                    "SLC18B1",
                                                    "SCF22M21",
                                                    "uncharacterised OG0002961",
                                                    "pug",
                                                    "malvolio",
                                                    "Hex70c",
                                                    "HEX70b",
                                                    "Hex70a",
                                                    "henna",
                                                    "fibroin heavy chain",
                                                    "cytochrome b5",
                                                    "4-domain protease inhibitor",
                                                    "serpin-5/88Ea",
                                                    "gnbp-1",
                                                    "cactus",
                                                    "3-PIDPK1"))


ggplot(data=fift.res2, mapping = aes(x = Contrast, y = Gene,
                                     fill = log2FoldChange, labels = Source)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0080ff", mid = "white", high = "#ff0000") +
  facet_grid(~Source) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  theme(strip.text.x = element_text(face = "italic")) +
  xlab("Condition (vs Naive)")

ggsave("output/Synthesis/3Bees_Heatmap_15genes.pdf")

ind.res[ind.res$OrthoGroup == "OG0000699",]

#notes



#three in aculeata
targ <- c("LOC724703", "LOC100646646", "Caust.v2_008285", "LOC106792475",
          "LOC100577795", "LOC100651554", "Caust.v2_004788", "LOC106790475",
          "LOC724367", "LOC100651433", "LOC106785762", "Caust.v2_010194")
three.res <- ind.all.res[ind.all.res$Gene %in% targ,]
three.res %>%
  arrange(OrthoGroup)

key <- unique(three.res$OrthoGroup)
key2 <- c("S/T-kinase", "L(2)EFL", "modular serine protease")
key <- as.data.frame(cbind(key,key2))

for(i in 1:nrow(three.res)){
  three.res$Gene2[i] <- key$key2[key$key == three.res$OrthoGroup[i]]
}

head(three.res)

ggplot(data=three.res, mapping = aes(x = Contrast, y = Gene,
                                     fill = log2FoldChange, labels = Source)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0080ff", mid = "white", high = "#ff0000") +
  facet_grid(~Source) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  theme(strip.text.x = element_text(face = "italic")) +
  xlab("Condition (vs Naive)")



ind.res[ind.res$OrthoGroup == "OG0002961",]

###Orthogroup Level ####
dir.ort <- "output/Ortho_DE"
targ <- list.dirs(dir.ort)
targ <- targ[-(grep("Filtered_Counts", targ))]
targ <- targ[c(-1, -2, -4, -6, -7, -9, -10)]

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

#leukocyte elastase inhibitor
ort.res[ort.res$OrthoGroup == "OG0003982",]
#ninjurin-1
ort.res[ort.res$OrthoGroup == "OG0004924",]
#cactus
ort.res[ort.res$OrthoGroup == "OG0000699",]
#toll
ort.res[ort.res$OrthoGroup == "OG0000399",]
#art3
ort.res[ort.res$OrthoGroup == "OG0005280",]

#heatmap time
dir.ort <- "output/Ortho_DE"
targ <- list.dirs(dir.ort)
targ <- targ[-(grep("Filtered_Counts", targ))]
targ <- targ[c(-1, -2, -4, -6, -7, -9, -10)]

ort.list <- vector(mode = "list", length = length(targ))
for ( i in 1:length(ort.list)){
  file <- list.files(path = paste(targ[i], sep ="/"), 
                     pattern = 'DeSeqResults_All.tsv')
  spec <- strsplit(targ[i], "/")[[1]][3]
  ort.list[[i]] <- read.table(paste(targ[i], file, sep = "/"),
                              header = T, sep = "\t", quote = "")
  ort.list[[i]]$Source <- paste(spec)
}
ort.all.res <- bind_rows(ort.list)

five.res <- ort.res[ort.all.res$OrthoGroup == "OG0003982" |
                      ort.all.res$OrthoGroup == "OG0004924" |
                      ort.all.res$OrthoGroup == "OG0000699"|
                      ort.all.res[ort.res$OrthoGroup == "OG0000399" |
                      ort.all.res[ort.all.res$OrthoGroup == "OG0005280", ]
                      