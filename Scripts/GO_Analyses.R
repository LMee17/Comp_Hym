#28th July 2022
#GO analyses for the individual species runs, significantly up or down regulated genes

#BiocManager::install("topGO")

library("topGO")
library("dplyr")
library("ggplot2")
library("ggwordcloud")


#####Running Analyses####

##Resources
ortho <- read.table("input/Ortho_Files/OrthoFinder_Combined_Parsed_Jul22.tsv",
                    header = T, sep = "\t")
#manually add apidaecin orthogroup
ortho[40559,] <- c("OG0011065", "Amel", "Apid1")

#result tsvs
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
ind.res <- bind_rows(ind.list)

runGO <- function(goiList, species, contrast, runname, deDF){
  goi <- goiList
  if (species == "Bter"){
    geneID2GO <- readMappings(paste("input/Gene_Info/", species, 
                                    "_GOverseTopGO_Jun22.tsv", sep= ""))
  } else {
    geneID2GO <- readMappings(paste("input/Gene_Info/", species, 
                                    "_GOverseTopGO_Apr22.tsv", sep= ""))
  }
  y <- deDF[!is.na(deDF$padj),]
  y <- y[grepl(paste(contrast), y$Contrast),]
  geneUniverse <- unique(y$Gene)
  geneList <- factor(as.integer(geneUniverse %in% goi))
  names(geneList) <- geneUniverse
  myGOdata<-new("topGOdata",
                description=paste(species, "_", runname, sep = ""), 
                ontology="BP",
                allGenes=geneList, 
                annot=annFUN.gene2GO, 
                gene2GO=geneID2GO)
  resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
  bpRes <- GenTable(myGOdata, weightFisher = resultFisher,
                    orderBy = "resultFisher", topNodes = 200)
  bpRes$Ontology <- "BP"
  myGOdata<-new("topGOdata",
                description=paste(species, "_", runname, sep = ""), 
                ontology="CC",
                allGenes=geneList, 
                annot=annFUN.gene2GO, 
                gene2GO=geneID2GO)
  resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
  ccRes <- GenTable(myGOdata, weightFisher = resultFisher,
                    orderBy = "resultFisher", topNodes = 200)
  ccRes$Ontology <- "CC"
  myGOdata<-new("topGOdata",
                description=paste(species, "_", runname, sep = ""), 
                ontology="MF",
                allGenes=geneList, 
                annot=annFUN.gene2GO, 
                gene2GO=geneID2GO)
  resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
  mfRes <- GenTable(myGOdata, weightFisher = resultFisher,
                    orderBy = "resultFisher", topNodes = 200)
  mfRes$Ontology <- "MF"
  goRes <- rbind(bpRes, ccRes, mfRes)
  write.table(goRes, paste("output/Ind_DE/", species,
                          "/", species, runname, "GOResults.tsv", sep = ""),
              col.names = T, row.names = F, quote = F, sep = "\t")
}

getGOI <- function(contrast, regDirection, deDF){
  x <- deDF[!is.na(deDF$padj),]
  x <- x[x$padj < 0.1,]
  if ( regDirection == "Up" ){
    goi <- x$Gene[grepl(paste(contrast), x$Contrast) &
                       x$log2FoldChange > 1]
  }
  if ( regDirection == "Down"){
    goi <- x$Gene[grepl(paste(contrast), x$Contrast) &
                       x$log2FoldChange < 1]
  }
  goi  <- unique(goi)
  return(goi)
}

getGOI2 <- function(contrast, deDF){
  x <- deDF[!is.na(deDF$padj),]
    goi <- x$Gene[grepl(paste(contrast), x$Contrast) &
                          x$padj < 0.1]
    goi  <- unique(goi)
    return(goi)
}

#See if this works. Start with Wound Up, Amel
x <- getGOI("Wound", "Up", ind.list[[1]])
runGO(x, 
      species = "Amel", 
      contrast = "Wound", 
      runname = "Wound_UpReg",  
      deDF = ind.list[[1]])

x <- getGOI("Wound", "Down", ind.list[[1]])
runGO(x, 
      species = "Amel", 
      contrast = "Wound", 
      runname = "Wound_DownReg",  
      deDF = ind.list[[1]])

treats <- c("Wound", "GramPos", "GramNeg")
#Ok, iterate
for (i in 1:length(ind.list)){
  for(j in 1:length(treats)){
    x <- getGOI(treats[j], "Up", ind.list[[i]])
    species <- unique(ind.list[[i]]$Source)
    runGO(x, species = species, contrast = treats[i], 
          runname = paste(treats[j], "UpReg", sep = "_"),
    deDF = ind.list[[i]])
    x <- getGOI(treats[j], "Down", ind.list[[i]])
    species <- unique(ind.list[[i]]$Source)
    runGO(x, species = species, contrast = treats[i], 
          runname = paste(treats[j], "DownReg", sep = "_"),
    deDF = ind.list[[i]])
  }
}

#this isn't working for Polistes and I have no idea why. When I break the function down
#and run it directly it works, but not when its input

#polistes
x <- getGOI(treats[1], "Up", ind.list[[4]])
species <- unique(ind.list[[4]]$Source)
runGO(x, species = species, contrast = treats[1], 
      runname = paste(treats[1], "UpReg", sep = "_"),
      deDF = ind.list[[4]])
x <- getGOI(treats[1], "Down", ind.list[[4]])
runGO(x, species = species, contrast = treats[1], 
      runname = paste(treats[1], "DownReg", sep = "_"),
      deDF = ind.list[[4]])
x <- getGOI(treats[2], "Up", ind.list[[4]])
runGO(x, species = species, contrast = treats[1], 
      runname = paste(treats[2], "UpReg", sep = "_"),
      deDF = ind.list[[4]])
x <- getGOI(treats[2], "Down", ind.list[[4]])
runGO(x, species = species, contrast = treats[1], 
      runname = paste(treats[2], "DownReg", sep = "_"),
      deDF = ind.list[[4]])
x <- getGOI(treats[3], "Up", ind.list[[4]])
runGO(x, species = species, contrast = treats[1], 
      runname = paste(treats[3], "UpReg", sep = "_"),
      deDF = ind.list[[4]])
x <- getGOI(treats[3], "Down", ind.list[[4]])
runGO(x, species = species, contrast = treats[1], 
      runname = paste(treats[3], "DownReg", sep = "_"),
      deDF = ind.list[[4]])

#To make visual plots, just looking at all genes at once

#Again, won't work in a look for some fucking reason.

x <- getGOI2(contrast = "Wound", deDF = ind.list[[1]])
runGO(goiList = x, species = "Amel", 
      contrast = "Wound", 
      runname = "Amel_Wound_All", 
      deDF = ind.list[[1]])

for (i in 1:length(ind.list)){
  for (j in 1:length(treats)){
    x <- getGOI2(contrast = paste(treats[j]), deDF = ind.list[[i]])
    spec <- unique(ind.list[[i]]$Source)
    runGO(goiList = x, 
          species = spec, 
          contrast = treats[j], 
          runname = paste(treats[j], "All", sep = "_"), 
          deDF = ind.list[[i]])
  }
}

####Add Up/ Down Regulated #####

visiDirect <- function(species, experiment, goiList, runname, deDF){
  go.df  <- read.table(paste("output/Ind_DE/", species, "/",
                             species, experiment, "_AllGOResults.tsv", sep = ""),
                       header = T, sep = "\t", quote = "")
  go.df <- go.df[go.df$weightFisher < 0.05,]
  Reg <- vector(length = length(goiList))
  for (i in 1:length(goiList)){
    goi.df <- deDF$log2FoldChange[deDF$Gene == goiList[i] &
                                    grepl(experiment, deDF$Contrast)]
    print(goi.df)
    if(goi.df > 0){
      Reg[i] <- "Up"
    } else {
      Reg[i] <- "Down" 
    }
  }
  Reg
  goiReg <- as.data.frame(cbind(goiList, Reg))
  if (species == "Bter"){
    goverse <- read.table(paste("input/Gene_Info/", species, 
                                "_GOverseTopGO_Jun22.tsv", sep= ""),
                          header = F, sep = "\t")
  } else {
    goverse <- read.table(paste("input/Gene_Info/", species, 
                                "_GOverseTopGO_Apr22.tsv", sep= ""),
                          header = F, sep = "\t")
  }
  goList <- unique(go.df$GO.ID)
  for(i in 1:length(goList)){
    go.genes <- goverse$V1[grepl(goList[i], goverse$V2) &
                             goverse$V1 %in% goiList]
    go.reg <- goiReg$Reg[goiReg$goiList %in% go.genes]
    go.reg.mode <- go.reg[which.max(table(go.reg))]
    go.df$Regulation[go.df$GO.ID == goList[i]] <- go.reg.mode
  }
  write.table(go.df, paste("output/Ind_DE/", species, "/", species, "_",
                           experiment, "_GOAnalysis_Directional.tsv", sep = ""),
              col.names = F, row.names = F, quote = F, sep = "\t")    
}

names(ind.list[[1]])
#hopefully the loop works
for (i in 1:length(ind.list)){
  spec <- unique(ind.list[[i]]$Source)
  print(spec)
  for (j in 1:length(treats)){
    print(treats[j])
    x <- getGOI2(contrast = paste(treats[j]), deDF = ind.list[[i]])
    visiDirect(species = spec, 
               experiment = paste(treats[j]), 
               goiList = x, 
               runname = paste(spec, treats[j], sep = "_"), 
               deDF = ind.list[[i]])
  }
}

sessionInfo()
writeLines(capture.output(sessionInfo()),
            "SessionLogs/GoAnalyses_SesInfo.txt")  
