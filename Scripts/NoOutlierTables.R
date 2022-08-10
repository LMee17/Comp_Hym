#8th August 2022
#Producing combined results tables for individual and orthogroup DE analyses
#with outliers kept in

library("dplyr")

#read in and annotate individual results
dir.ind <- "output/Ind_DE"
targ <- list.dirs(dir.ind)
targ <- targ[-(grep("Plots", targ))]
targ <- targ[c(-1, -3, -5, -7)]

ind.list <- vector(mode = "list", length = length(targ))
for ( i in 1:length(ind.list)){
  file <- list.files(path = paste(targ[i], sep ="/"), 
                     pattern = 'DeSeqResults_All.tsv')
  spec <- strsplit(targ[i], "/")[[1]][3]
  ind.list[[i]] <- read.table(paste(targ[i], file, sep = "/"),
                              header = T, sep = "\t", quote = "")
  ind.list[[i]]$Source <- paste(spec)
}

#add ccal gene descriptions to caus results
caus <- ind.list[[3]]
ccal <- read.table("~/Downloads/GCF_001652005.1_ASM165200v1_protein.faa",
                   header = F, sep = "\t", quote = "")
names(caus)[9] <- "Ccal_Prot"
j <- 1
for (i in 1:nrow(caus)){
  k <- j/1000
  if ((k == round(k)) == TRUE){
    print(j)
  }
  check <- unique(grepl(caus$Ccal_Prot[i], ccal$V1))
  if (length(check) == 1 ){
    caus$Desc[i] <- "NA"
  } else {
    caus$Desc[i] <- ccal$V2[grepl(caus$Ccal_Prot[i], ccal$V1)]
  }
  j <- j + 1
}
head(caus)
caus$Ccal_Prot <- NULL

ind.res <- bind_rows(ind.list)
ind.res <- ind.res[!ind.res$Source=="Caus_12",]
ind.res <- rbind(ind.res, caus)
ind.res <- ind.res[!ind.res$Source == "Amel",]

summary(as.factor(ind.res$Source))
ind.res$Source[ind.res$Source == "Plan_12"] <- "P.lanio"
ind.res$Source[ind.res$Source == "Bter_12"] <- "B.terrestris"
ind.res$Source[ind.res$Source == "Caus_12"] <- "C.australensis"


write.table(ind.res, "output/Synthesis/AllInd_AllSamples_DESeq2Results.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

#orthogroup level
dir.ort <- "output/Ortho_DE"
targ <- list.dirs(dir.ort)
targ <- targ[-(grep("Filtered_Counts", targ))]
targ <- targ[c(-1, -3, -4, -5, -8, -9, -11)]

ort.list <- vector(mode = "list", length = length(targ))
for ( i in 1:length(ort.list)){
  file <- list.files(path = paste(targ[i], sep ="/"), 
                     pattern = 'DeSeqResults_All.tsv')
  spec <- strsplit(targ[i], "/")[[1]][3]
  ort.list[[i]] <- read.table(paste(targ[i], file, sep = "/"),
                              header = T, sep = "\t", quote = "")
  ort.list[[i]]$Source <- paste(spec)
}
ort.res <- bind_rows(ort.list)

#give the orthogroups a description, if possible
#using amel as an anchor species (and as its in all combinations)

ortho <- read.table("input/Ortho_Files/OrthoFinder_Combined_Parsed_Aug22.tsv",
                    sep = "\t", header =T)

for (i in 1:nrow(ort.res)){
  print(i)
  one <- head(ortho$GeneID[ortho$Orthogroup == ort.res$OrthoGroup[i] &
                             ortho$Species == "Amel"], n= 1)
  ort.res$Gene[i] <- one
  two <- unique(ind.list[[1]]$Desc[ind.list[[1]]$Gene == one])
  if (length(two) > 0){
    ort.res$Desc[i] <- two
  } else {
    ort.res$Desc[i] <- "NA"
  }
}
tail(ort.res)

write.table(ort.res, "output/Synthesis/AllOrt_AllSamples_DESeq2Results.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

summary(as.factor(ort.res$Source))
ort.res %>%
  subset(Source == "Apinae") %>%
  subset(padj < 0.1) %>%
  subset(grepl("GramPos", Contrast)) %>%
  subset(log2FoldChange < 0) %>%
  nrow()
