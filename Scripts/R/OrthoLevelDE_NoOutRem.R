#August2022
#Orthogroup-level DE analyses without removing outliers

set.seed(42)

##Libraries####

library("DESeq2")
library("tidyverse")
library("ggrepel")
library("tximport")
library("pheatmap")
library("reshape")
library("ggvenn")

##Resources####

ortho <- read.table("input/Ortho_Files/OrthoFinder_Combined_Parsed_Jul22.tsv",
                    sep = "\t", header = T)

ann <- read.table("input/Gene_Info/Amel_ImmAnnotations_Orthogroup_Jun22.tsv",
                  header = T, sep = "\t")
ann.slim <- ann[,c(2:3)]

#Fix Kr-h1 issue (distinction not need for this run)

ortho$GeneID[grep("Kr-h1", ortho$GeneID)] <- "Kr-h1"
ortgps <- as.data.frame(unique(ortho$Orthogroup))

for(i in 1:nrow(ortgps)){
  species <- ortho$Species[ortho$Orthogroup == ortgps[i,1]]
  ortgps[i,2] <- ifelse("Amel" %in% species, 1, 0)
  ortgps[i,3] <- ifelse("Bter" %in% species, 1, 0)
  ortgps[i,4] <- ifelse("Caus" %in% species, 1, 0)
  ortgps[i,5] <- ifelse("Pcan" %in% species, 1, 0)
}

names(ortgps) <- c("OrthoGroup", "Amel", "Bter",
                   "Caus", "Plan")

#Make tx object with orthogroup IDs in place of GeneIDs

runs <- c("Amel", "Bter", "Caus", "Plan")

tx.list <- vector(mode = "list", length = length(runs))

for (i in 1:length(runs)){
  tx.list[[i]] <- read.table(paste("input/Gene_Info/", runs[i],
                                   "_TXGene.tsv", sep = ""),
                             header = T, sep = "\t")
  for (j in 1:nrow(tx.list[[i]])){
    check <- tx.list[[i]][j,2] %in% ortho$GeneID
    if (check == FALSE){
      tx.list[[i]]$ORTHOID[j] <- paste("Unassigned", runs[i], sep = "_")
    } else {
      tx.list[[i]]$ORTHOID[j] <- unique(ortho$Orthogroup
                                        [ortho$GeneID == tx.list[[i]][j,2]])
    }
  }
  tx.list[[i]] <- tx.list[[i]][,c(1,3)]
}

ort.tx <- bind_rows(tx.list)

##Functions####

#Function that takes in raw abundance.tsv files and keeps only those transcripts 
#that are members of orthogroups that have representatives of all input species. These are then collected into a folder reader to be read into tximport.

filterCount <- function(specieslist, outname){
  #object read to read in orthogroups that have members in all given species
  ort.grp <- vector(mode = "list", length = length(specieslist))
  #iterate through empty list and fill each with list of orthogroups 
  #that species has been assigend to
  for (i in 1:length(ort.grp)){
    ort.grp[[i]] <- ortgps$OrthoGroup[ortgps[,paste(specieslist[i])] == 1]
  }
  #combine list of orthogroups into one vector
  ort.grp <- unlist(ort.grp)
  #put into table format with orthogroup as dimname and frequency as value
  #ie the number of species assigned to that orthogroup
  x <- table(ort.grp)
  #how many orthogroups are present that have all the species listed as input included
  tot <- length(names(x[x == length(specieslist)]))
  #print to screen
  print(paste("There are", tot, "orthogroups shared between these species", sep = " "))
  #record these target orthogroups
  tar.ort <- names(x[x==length(specieslist)])
  #collect a list of transcripts associated with these orthogroups
  tx.out <- ort.tx$TXNAME[ort.tx$ORTHOID %in% tar.ort]
  #prepare to produce output
  dir.create(paste("output/Ortho_DE/", outname, sep = ""))
  dir.create(paste("output/Ortho_DE/", outname, "/Filtered_Counts/", sep = ""))
  #for each species, read in each samples abundance file and keep only rows
  #associated with transcripts that correspond with target orthogroups
  for (s in 1:length(specieslist)){
    dir <- paste("input/Counts_Kal/", specieslist[s], sep = "")
    files <- list.files(dir)
    for (f in 1:length(files)){
      entrer <- read.table(paste(dir, files[f], sep = "/"),
                           header = T, sep = "\t")
      sortir <- entrer[entrer$target_id %in% tx.out,]
      splitf <- strsplit(files[f], ".", fixed =T)
      fname <- paste(splitf[[1]][1], ".orthoFiltered_", 
                     splitf[[1]][2], ".", splitf[[1]][3],
                     sep = "")
      write.table(sortir,
                  paste("output/Ortho_DE/", outname, "/Filtered_Counts/", 
                        fname, sep =""),
                  col.names = T, row.names = F, quote = F, sep ="\t")
    }
  }
}

#A function that reads all filtered count files in a given orthogroup work directory
#and feed them into tximport in order to produce a list of txobjects (one per input species). Coldata object is produced from file names.

ortCount <- function(ortrun, txgeneobjlist){
  dir <- paste("output/Ortho_DE/", ortrun, "/Filtered_Counts", sep = "")
  files <- list.files(dir)
  samples <- sapply(strsplit(files, ".", fixed = T), "[", 1)
  treat <- vector(length = length(samples))
  species <- vector(length = length(samples))
  for (i in 1:length(samples)){
    if(grepl("_N_", samples[i]) == T){
      treat[i] <- "Naive"
    } else {
      if(grepl("_P_", samples[i]) == T){
        treat[i] <- "Wound"
      } else {
        if(grepl("_SM_", samples[i]) == T){
          treat[i] <- "GramPos"
        } else {
          if(grepl("_SL_", samples[i]) == T){
            treat[i] <- "GramNeg"
          }
        }
      }    
    }
  }
  for (i in 1:length(samples)){
    if(grepl("AM_", samples[i]) == T){
      species[i] <- "Amel"
    } else {
      if(grepl("BT_", samples[i]) == T){
        species[i] <- "Bter"
      } else {
        if(grepl("CA_", samples[i]) == T){
          species[i] <- "Caus"
        } else {
          if(grepl("PL_", samples[i]) == T){
            species[i] <- "Plan"
          }
        }
      }    
    }
  }
  cd <- as.data.frame(cbind(samples, treat, species))
  names(cd) <- c("Samples", "Treatment", "Species")
  #write up coldata
  write.table(cd, 
              paste("output/Ortho_DE/", ortrun, "/", ortrun, "_SampleData.tsv", sep = ""),
              col.names = T, row.names = F, sep = "\t", quote = F)
  #name the files with sample IDS
  names(files) <- samples
  cnts <- file.path(dir, files)
  names(cnts) <- samples
  #each species will have its own txobject that will be combined into one for DE
  sets <- unique(species)
  tx.set <- vector(mode = "list", length = length(sets))
  #run through each species "set", reading in only those abundance files
  for (s in 1:length(sets)){
    print(paste(sets[s], "..." ))
    samp.set <- cd$Samples[cd$Species == paste(sets[s])]
    print(samp.set)
    cd.set <- cd[cd$Species == paste(sets[s]),]
    txgene.set <- txgeneobjlist[[s]]
    cnts.set <- cnts[names(cnts) %in% samp.set]
    tx.set[[s]] <- tximport(cnts.set, type = "kallisto",
                            tx2gene = txgene.set, ignoreAfterBar = T)
  }
  return(tx.set)
}

#Produce DDS object from count table object (includes making coldata object from 
#column names, prefiltering and setting reference condition to Naive)

getDDS <- function(ortrun, txobject){
  cd <- read.table(paste("output/Ortho_DE/", ortrun, "/",
                         ortrun, "_SampleData.tsv", sep = ""), header = T, sep = "\t")
  ddsTx <- DESeqDataSetFromTximport(txobject, colData = cd, design = ~Treatment)
  keep <- rowSums(counts(ddsTx)) >= 10
  ddsTx2 <- ddsTx[keep,]
  ddsTx2$Treatment <- relevel(ddsTx2$Treatment, ref = "Naive")
  ddsTx2 <- DESeq(ddsTx2)
  return(ddsTx2)
}

#Function to plot expression boxplot of a certain orthogroup separating by species. 
#Includes an all species variable.

visiOrt.messy <- function(orthogroup, ortname, dds){
  tmp <- counts(dds, normalized = T)
  samples <- dds$Samples
  counts <- tmp[rownames(tmp) == paste(orthogroup),]
  plot.df <- as.data.frame(cbind(samples, counts))
  names(plot.df) <- c("Samples", "Norm_Counts")
  plot.df$Species[grepl("AM_", plot.df$Samples)] <- "A mellifera"
  plot.df$Species[grepl("BT_", plot.df$Samples)] <- "B terrestris"
  plot.df$Species[grepl("CA_", plot.df$Samples)] <- "C australensis"
  plot.df$Species[grepl("PL_", plot.df$Samples)] <- "P lanio"
  plot.df$Treatment[grepl("_N_", plot.df$Samples)] <- "Naive"
  plot.df$Treatment[grepl("_P_", plot.df$Samples)] <- "Wound"
  plot.df$Treatment[grepl("_SM_", plot.df$Samples)] <- "Gram Negative"
  plot.df$Treatment[grepl("_SL_", plot.df$Samples)] <- "Gram Positive"
  plot.df$Treatment <- factor(plot.df$Treatment, levels = c("Naive",
                                                            "Wound", 
                                                            "Gram Positive", 
                                                            "Gram Negative"))
  plot.df$Norm_Counts <- as.numeric(plot.df$Norm_Counts)
  to.add <-plot.df
  to.add$Species <- "All Species"
  plot.df2 <- rbind(plot.df,to.add)
  plot.df2$Species <- factor(plot.df2$Species, levels = c("All Species",
                                                          "A mellifera",
                                                          "B terrestris",
                                                          "C australensis",
                                                          "P lanio"))
  p <- ggplot(plot.df2, aes(x = Treatment, y = log(Norm_Counts), fill = Species)) +
    geom_boxplot() +
    scale_fill_manual(values = c("grey", "#ffa758",  "#58b0ff", "#58ffa7", "#ff58b0")) +
    xlab("Treatment") +
    ylab("log(DESeq2 Normalised Counts)") +
    labs(title = (paste(ortname, " expression across treatments", sep = ""))) +
    guides(alpha = FALSE, colour = FALSE)
  ggsave(paste(wd, orthogroup, "_ExpressionBoxplot.messy.pdf", sep = ""))
  p
}

#Function to plot expression boxplot of a certain orthogroup separating by treatment. 
#Species are plotted as different shaped geom_points.

visiOrt.clean <- function(orthogroup, ortname, dds){
  tmp <- counts(dds, normalized = T)
  samples <- dds$Samples
  counts <- tmp[rownames(tmp) == paste(orthogroup),]
  plot.df <- as.data.frame(cbind(samples, counts))
  names(plot.df) <- c("Samples", "Norm_Counts")
  plot.df$Species[grepl("AM_", plot.df$Samples)] <- "A mellifera"
  plot.df$Species[grepl("BT_", plot.df$Samples)] <- "B terrestris"
  plot.df$Species[grepl("CA_", plot.df$Samples)] <- "C australensis"
  plot.df$Species[grepl("PL_", plot.df$Samples)] <- "P lanio"
  plot.df$Treatment[grepl("_N_", plot.df$Samples)] <- "Naive"
  plot.df$Treatment[grepl("_P_", plot.df$Samples)] <- "Wound"
  plot.df$Treatment[grepl("_SM_", plot.df$Samples)] <- "Gram Negative"
  plot.df$Treatment[grepl("_SL_", plot.df$Samples)] <- "Gram Positive"
  plot.df$Treatment <- factor(plot.df$Treatment, levels = c("Naive",
                                                            "Wound", 
                                                            "Gram Positive", 
                                                            "Gram Negative"))
  plot.df$Norm_Counts <- as.numeric(plot.df$Norm_Counts)
  plot.df$Species <- factor(plot.df$Species, levels = c(
    "A mellifera",
    "B terrestris",
    "C australensis",
    "P lanio"))
  p <- ggplot(plot.df, aes(x = Treatment, y = log(Norm_Counts), fill = Treatment)) +
    geom_boxplot(alpha = .8, outlier.colour = "red") +
    geom_point(aes(shape = Species), position = position_dodge(width = .2),
               size = 2) +
    scale_fill_manual(values = c("#ffa758",  "#58b0ff", "#58ffa7", "#ff58b0")) +
    scale_shape_manual(values = c(1, 0, 2, 4))+
    xlab("Treatment") +
    ylab("log(DESeq2 Normalised Counts)") +
    labs(title = (paste(ortname, " expression across treatments", sep = ""))) +
    guides(alpha = FALSE, colour = FALSE)
  ggsave(paste(wd, orthogroup, "_ExpressionBoxplot.clean.pdf", sep = ""))
  p
}

#Produce a dataframe with orthogroup and gene class information ready for data visualisation

makeMetaDf <- function(res.all.object){
  x <- res.all.object[!is.na(res.all.object$padj),]
  x$Condition[x$Contrast==paste(treats[1])] <- "Wound"
  x$Condition[x$Contrast==paste(treats[2])] <- "Gram Positive"
  x$Condition[x$Contrast==paste(treats[3])] <- "Gram Negative"
  for(i in 1:nrow(x)){
    check <- x$Gene[i] %in% ortho$GeneID
    if (check == FALSE){
      x$OrthoGroup[i] <- "Unassigned"
    } else {
      x$OrthoGroup[i] <- ortho$Orthogroup[ortho$GeneID ==
                                            x$Gene[i]]
    }
  }
  for(i in 1:nrow(x)){
    if (x$OrthoGroup[i] == "Unassigned") {
      x$Function[i] <- "Unknown"} else {
        check2 <- x$OrthoGroup[i] %in% ann.slim$Orthogroup
        if (check2 == FALSE){
          x$Function[i] <- "Unknown"
        } else {
          f <- unique(ann.slim$Function[ann.slim$Orthogroup == x$OrthoGroup[i]])
          if (length(f) > 1){
            if(TRUE %in% grepl("Effector|Recognition|Signalling", f)){
              x$Function[i] <- "Immune"
            } else {
              if(TRUE %in% grepl("Putative Immune", f)){
                x$Function[i] <- "Putative Immune"
              } 
            }
          } else {
            f2 <- gsub("Effector", "Immune", 
                       gsub("Recognition", "Immune", 
                            gsub("Signalling", "Immune", 
                                 gsub("Background", "None Immune", f))))
            x$Function[i] <- paste(f2)
          }
        }
      }
  }
  x$Sig <- "No"
  x$Sig[x$padj < 0.1] <- "Yes"
  x$Condition <- factor(x$Condition, levels = c("Wound", 
                                                "Gram Negative", 
                                                "Gram Positive"))
  x$Neglog10adjPvalue <- -log10(x$padj)
  x$Direction <- "Not Significant"
  x$Direction[x$Sig == "Yes" &
                x$log2FoldChange > 0 ] <- "Up-Regulated"
  x$Direction[x$Sig == "Yes" &
                x$log2FoldChange < 0 ] <- "Down-Regulated"
  return(x)
}

#Produce a 1/0 binary matrix for use in comparing overlap in genes of interest

makeMatrix <- function(goidf){
  mat <- as.data.frame(unique(goidf$goi))
  names(mat) <- "Gene"
  mat$SigWound <- ifelse(mat$Gene %in% sig.wound, 1, 0)
  mat$SigPos <- ifelse(mat$Gene %in% sig.pos, 1, 0)
  mat$SigNeg <- ifelse(mat$Gene %in% sig.neg, 1, 0)
  mat$TopWound <- ifelse(mat$Gene %in% top.wound, 1, 0)
  mat$TopPos <- ifelse(mat$Gene %in% top.pos, 1, 0)
  mat$TopNeg <- ifelse(mat$Gene %in% top.neg, 1, 0)
  mat$AllSig <- ifelse(mat$SigWound == 1 &
                         mat$SigPos == 1 &
                         mat$SigNeg == 1, 1, 0)
  mat$AllTop <- ifelse(mat$TopWound == 1 &
                         mat$TopPos == 1 &
                         mat$TopNeg == 1, 1, 0)
  mat$BacterialSig <- ifelse(mat$SigWound == 0 &
                               mat$SigPos == 1 &
                               mat$SigNeg == 1, 1, 0)
  mat$BacterialTop <- ifelse(mat$TopWound == 0 &
                               mat$TopPos == 1 &
                               mat$TopNeg == 1, 1, 0)
  mat$GoiWound <- ifelse(mat$TopWound == 1 | mat$SigWound == 1, 1, 0)
  mat$GoiPos <- ifelse(mat$TopPos == 1 | mat$SigPos == 1, 1, 0)
  mat$GoiNeg <- ifelse(mat$TopNeg == 1 | mat$SigNeg == 1, 1, 0)
  mat$AllGoi <- ifelse(mat$GoiWound == 1 &
                         mat$GoiPos == 1 &
                         mat$GoiNeg == 1, 1, 0)
  mat$BacterialGoi <- ifelse(mat$GoiPos == 1 &
                               mat$GoiNeg == 1  &
                               mat$GoiWound == 0, 1, 0)
  return(mat)
}

#Produce and save a boxplot with expression levels across samples with gene input

visiGene <- function(gene, genename, dds){
  p <- plotCounts(dds, gene = paste(gene), intgroup = "Treatment", returnData = T)
  p$Treatment2 <- as.character(p$Treatment)
  p$Treatment2[p$Treatment == "GramPos"] <- "Gram Positive"
  p$Treatment2[p$Treatment == "GramNeg"] <- "Gram Negative"  
  p$Treatment2 <- factor(p$Treatment2, levels = c("Naive",
                                                  "Wound", 
                                                  "Gram Positive", 
                                                  "Gram Negative"))  
  p2 <- ggplot(p, aes(x = Treatment2, y = count)) +
    geom_boxplot(aes(colour = Treatment2, alpha = 0.5)) +
    geom_point(position = position_jitter(w=0.1, h = 0),
               aes(colour = Treatment2), size = 3) +
    scale_colour_manual(values = c("#58ffa7", "#ffa758", "#ff58b0", "#58b0ff")) +
    xlab("Treatment") +
    labs(title = (paste(genename, " expression across treatments", sep = ""))) +
    guides(alpha = FALSE, colour = guide_legend("Treatment"))
  ggsave(paste(wd, "Plots/ExpressionPlots/", gene, "_ExpBoxPlot.pdf", sep = ""))
  p2
}

#Plot all expression boxplots from a group of genes 

visiFacet <- function(genes, genenames, dds, filename){
  p <- vector(mode = "list", length = length(genes))
  for (i in 1:length(p)){
    p[[i]] <- plotCounts(dds, gene = genes[i], 
                         intgroup = "Treatment", returnData = T)
    p[[i]]$Gene <- paste(genes[i])
    p[[i]]$Label <- paste(genenames[i])
  }
  p2 <- bind_rows(p)
  p2$Treatment2 <- as.character(p2$Treatment)
  p2$Treatment2[p2$Treatment == "GramPos"] <- "Gram Positive"
  p2$Treatment2[p2$Treatment == "GramNeg"] <- "Gram Negative"  
  p2$Treatment2 <- factor(p2$Treatment2, levels = c("Naive",
                                                    "Wound", 
                                                    "Gram Positive", 
                                                    "Gram Negative"))  
  p3 <- ggplot(p2, aes(x = Treatment2, y = count)) +
    geom_boxplot(aes(colour = Treatment2, alpha = 0.5)) +
    geom_point(position = position_jitter(w=0.1, h = 0),
               aes(colour = Treatment2), size = 3) +
    scale_colour_manual(values = c("#58ffa7", "#ffa758", "#ff58b0", "#58b0ff")) +
    xlab("Treatment") +
    guides(alpha = FALSE, colour = guide_legend("Treatment")) +
    facet_wrap(~Label, nrow = 5, scale = "free_y")
  ggsave(paste(wd, "Plots/ExpressionPlots/", filename, "_Facet_ExpBoxPlot.pdf", sep = ""))
  p3
}

##OrthoLevel DE: All Species (Aculeata, all samples)####

###Preparing TX Object

#Looking at orthogroups present in all four species.
#First, produce files that have non-orthogroup associated transcripts removed and 
#save them in output/Ortho_DE/Aculeata/ ready to be ready to be fed into tximport.

filterCount(runs, "Aculeata")

wd <- "output/Ortho_DE/Aculeata/"
length(list.files(paste(wd, "/Filtered_Counts", sep = "")))

#And now to make a list of txobjects (one per species) from these filtered count tables.

ort.48.txs <- ortCount("Aculeata", tx.list)

#make combined abundance matrix
abundance <- cbind(ort.48.txs[[1]]$abundance,
                   ort.48.txs[[2]]$abundance,
                   ort.48.txs[[3]]$abundance,
                   ort.48.txs[[4]]$abundance)
nrow(abundance)
#make combined counts matrix
counts <- cbind(ort.48.txs[[1]]$counts,
                   ort.48.txs[[2]]$counts,
                   ort.48.txs[[3]]$counts,
                   ort.48.txs[[4]]$counts)
nrow(counts)
#make combined lengths matrix
length <- cbind(ort.48.txs[[1]]$length,
                   ort.48.txs[[2]]$length,
                   ort.48.txs[[3]]$length,
                   ort.48.txs[[4]]$length)
nrow(length)
ort.48.tx <- list(abundance = abundance,
                  counts = counts, 
                  length = length)

ort.48.tx$countsFromAbundance <- "no"

for (i in 1:length(ort.48.tx)){
  print(head(ort.48.tx[[i]], n = 1))
}

###DE Analysis

ort.48.dds <- getDDS("Aculeata", ort.48.tx)
ort.48.res <- vector(mode = "list", length = 3)
treats <- resultsNames(ort.48.dds)[c(4,3,2)]

for(i in 1:length(ort.48.res)){
  ort.48.res[[i]] <- results(ort.48.dds, name = paste(treats[i]))
  print(treats[i])
  summary(ort.48.res[[i]])
}

###Results

#Save 

for (i in 1:length(ort.48.res)){
  ort.48.res[[i]]$Contrast <- paste(treats[i])
  ort.48.res[[i]]$OrthoGroup <- rownames(ort.48.res[[i]])
  rownames(ort.48.res[[i]]) <- NULL
  ort.48.res[[i]] <- as.data.frame(ort.48.res[[i]])
}

ort.48.res.all <- bind_rows(ort.48.res[1:3])

write.table(ort.48.res.all, paste(wd, "Aculeata_DeSeqResults_All.tsv", sep = ""),
            col.names = T, row.names = F, quote = F, sep = "\t")

ort.48.res.sig <- subset(ort.48.res.all, padj < 0.1)

write.table(ort.48.res.sig, paste(wd, "Aculeata_DeSeqResults_FDRlt0.1.tsv", sep =""), 
            col.names = T, row.names = F, quote = F, sep = "\t")

###ExpBoxPlots

#OG0003913: tyrosine hydrolase

visiOrt.messy("OG0003913", "Tyrosine hydrolase", ort.48.dds)

visiOrt.clean("OG0003913", "Tyrosine hydrolase", ort.48.dds)

#OG0004924: ninjurin-1

visiOrt.messy("OG0004924", "Ninjurin-1", ort.48.dds)

visiOrt.clean("OG0004924", "Ninjurin-1", ort.48.dds)

#OG0003982: leukocyte elastase inhibitor (serpin)

visiOrt.messy("OG0003982", "Leukocyte elastase inhibitor", ort.48.dds)

visiOrt.clean("OG0003982", "Leukocyte elastase inhibitor", ort.48.dds)

##OrthoLevel DE: All Bees (Anthophila, all samples)####

###Preparing TX Object
#Looking at orthogroups present in all three bee species.

#First, produce files that have non-orthogroup associated transcripts removed and 
#save them in output/Ortho_DE/Aculeata/ ready to be ready to be fed into tximport.

runs <- c("Amel", "Bter", "Caus")
filterCount(runs, "Anthophila")

wd <- "output/Ortho_DE/Anthophila/"
length(list.files(paste(wd, "/Filtered_Counts", sep = "")))

#And now to make a list of txobjects (one per species) from these filtered count tables.

ort.36.txs <- ortCount("Anthophila", tx.list)

#make combined abundance matrix
abundance <- cbind(ort.36.txs[[1]]$abundance,
                   ort.36.txs[[2]]$abundance,
                   ort.36.txs[[3]]$abundance)
nrow(abundance)
#make combined counts matrix
counts <- cbind(ort.36.txs[[1]]$counts,
                ort.36.txs[[2]]$counts,
                ort.36.txs[[3]]$counts)
nrow(counts)
#make combined lengths matrix
length <- cbind(ort.36.txs[[1]]$length,
                ort.36.txs[[2]]$length,
                ort.36.txs[[3]]$length)
nrow(length)
ort.36.tx <- list(abundance = abundance,
                  counts = counts, 
                  length = length)

ort.36.tx$countsFromAbundance <- "no"

for (i in 1:length(ort.36.tx)){
  print(head(ort.36.tx[[i]], n = 1))
}

###DE Analysis

ort.36.dds <- getDDS("Anthophila", ort.36.tx)
ort.36.res <- vector(mode = "list", length = 3)
treats <- resultsNames(ort.36.dds)[c(4,3,2)]

for(i in 1:length(ort.36.res)){
  ort.36.res[[i]] <- results(ort.36.dds, name = paste(treats[i]))
  print(treats[i])
  summary(ort.36.res[[i]])
}

###Results

#Save 

for (i in 1:length(ort.36.res)){
  ort.36.res[[i]]$Contrast <- paste(treats[i])
  ort.36.res[[i]]$OrthoGroup <- rownames(ort.36.res[[i]])
  rownames(ort.36.res[[i]]) <- NULL
  ort.36.res[[i]] <- as.data.frame(ort.36.res[[i]])
}

ort.36.res.all <- bind_rows(ort.36.res[1:3])

write.table(ort.36.res.all, paste(wd, "Anthophila_DeSeqResults_All.tsv", sep = ""),
            col.names = T, row.names = F, quote = F, sep = "\t")

ort.36.res.sig <- subset(ort.36.res.all, padj < 0.1)

write.table(ort.36.res.sig, paste(wd, "Anthophila_DeSeqResults_FDRlt0.1.tsv", sep =""), 
            col.names = T, row.names = F, quote = F, sep = "\t")

###ExpBoxPlots

sig.ort <- ort.36.res.sig$OrthoGroup
ann.slim[ann.slim$Orthogroup %in% sig.ort,]

#OG0004924: ninjurin-1

visiOrt.messy("OG0004924", "Ninjurin-1", ort.36.dds)

visiOrt.clean("OG0004924", "Ninjurin-1", ort.36.dds)

#OG0007158: Sorbitol dehydrogenase

visiOrt.messy("OG0007158", "Sorbitol dehydrogenase", ort.36.dds)

visiOrt.clean("OG0007158", "Sorbitol dehydrogenase", ort.36.dds)

#OG0003982: leukocyte elastase inhibitor (serpin)
visiOrt.messy("OG0003982", "Leukocyte elastase inhibitor", ort.36.dds)

visiOrt.clean("OG0003982", "Leukocyte elastase inhibitor", ort.36.dds)

#OG0007195: mt (mitochondrial) lipoyltransferase 1

visiOrt.messy("OG0007195", "Lipoyltransferase 1 (mt)", ort.36.dds)
visiOrt.clean("OG0007195", "Lipoyltransferase 1 (mt)", ort.36.dds)

#OG0008153: fibroin heavy chain

visiOrt.messy("OG0008153", "Fibroin heavy chain", ort.36.dds)
visiOrt.clean("OG0008153", "Fibroin heavy chain", ort.36.dds)

##OrthoLevel DE: Apis and Bombus (Apinae, all samples)####

###Preparing TX Object

#First, produce files that have non-orthogroup associated transcripts removed and 
#save them in output/Ortho_DE/Aculeata/ ready to be ready to be fed into tximport.

runs <- c("Amel", "Bter")
filterCount(runs, "Apinae")

wd <- "output/Ortho_DE/Apinae/"
length(list.files(paste(wd, "/Filtered_Counts", sep = "")))

#And now to make a list of txobjects (one per species) from these filtered count tables.

ort.24.txs <- ortCount("Apinae", tx.list)

#make combined abundance matrix
abundance <- cbind(ort.24.txs[[1]]$abundance,
                   ort.24.txs[[2]]$abundance)
nrow(abundance)
#make combined counts matrix
counts <- cbind(ort.24.txs[[1]]$counts,
                   ort.24.txs[[2]]$counts)
nrow(counts)
#make combined lengths matrix
length <- cbind(ort.24.txs[[1]]$length,
                   ort.24.txs[[2]]$length)
nrow(length)
ort.24.tx <- list(abundance = abundance,
                  counts = counts, 
                  length = length)

ort.24.tx$countsFromAbundance <- "no"

for (i in 1:length(ort.24.tx)){
  print(head(ort.24.tx[[i]], n = 1))
}

###DE Analysis

ort.24.dds <- getDDS("Apinae", ort.24.tx)
ort.24.res <- vector(mode = "list", length = 3)
treats <- resultsNames(ort.24.dds)[c(4,3,2)]

for(i in 1:length(ort.24.res)){
  ort.24.res[[i]] <- results(ort.24.dds, name = paste(treats[i]))
  print(treats[i])
  summary(ort.24.res[[i]])
}

###Results

#Save 

for (i in 1:length(ort.24.res)){
  ort.24.res[[i]]$Contrast <- paste(treats[i])
  ort.24.res[[i]]$OrthoGroup <- rownames(ort.24.res[[i]])
  rownames(ort.24.res[[i]]) <- NULL
  ort.24.res[[i]] <- as.data.frame(ort.24.res[[i]])
}

ort.24.res.all <- bind_rows(ort.24.res[1:3])

write.table(ort.24.res.all, paste(wd, "Apinae_DeSeqResults_All.tsv", sep = ""),
            col.names = T, row.names = F, quote = F, sep = "\t")

ort.24.res.sig <- subset(ort.24.res.all, padj < 0.1)

write.table(ort.24.res.sig, paste(wd, "Apinae_DeSeqResults_FDRlt0.1.tsv", sep =""), 
            col.names = T, row.names = F, quote = F, sep = "\t")

###ExpBoxPlots

ort.24.res.sig

####Upregulated Gram Pos

orts.up.pos <- ort.24.res.sig$OrthoGroup[ort.24.res.sig$log2FoldChange > 0 &
                                           ort.24.res.sig$Contrast == "Treatment_GramPos_vs_Naive"]

ort.24.res.sig[ort.24.res.sig$log2FoldChange > 0 &
                 ort.24.res.sig$Contrast == "Treatment_GramPos_vs_Naive",]

ortho[ortho$Orthogroup == "OG0008078",]

#####OG0000399: Toll
visiOrt.messy("OG0000399", "Toll", ort.24.dds)

visiOrt.clean("OG0000399", "Toll", ort.24.dds)

#####OG0000699: Cactus
ann.slim[ann.slim$Orthogroup == "OG0000699",]
visiOrt.messy("OG0000699", "Cactus", ort.24.dds)

visiOrt.clean("OG0000399", "Cactus", ort.24.dds)
#####OG0001164: Waprin-Thr1 

ann.slim[ann.slim$Orthogroup == "OG0001164",]
visiOrt.messy("OG0001164", "Waprin-Thr1 ", ort.24.dds)

visiOrt.clean("OG0001164", "Waprin-Thr1 ", ort.24.dds)
#####OG0001573: Carcinine transporter

ann.slim[ann.slim$Orthogroup == "OG0001573",]
visiOrt.messy("OG0001573", "Carcinine transporter", ort.24.dds)

visiOrt.clean("OG0001573", "Carcinine transporter", ort.24.dds)
#####OG0001687: Neurotrimin
ann.slim[ann.slim$Orthogroup == "OG0001687",]
visiOrt.messy("OG0001687", "Neurotrimin", ort.24.dds)

visiOrt.clean("OG0001687", "Neurotrimin", ort.24.dds)

#####OG0001733: OG0001733

ann.slim[ann.slim$Orthogroup == "OG0001733",]
visiOrt.messy("OG0001733", "OG0001733", ort.24.dds)

visiOrt.clean("OG0001733", "OG0001733", ort.24.dds)

#####OG0002754: Yellow-e3

ann.slim[ann.slim$Orthogroup == "OG0002754",]

visiOrt.messy("OG0002754", "Yellow-e3", ort.24.dds)

visiOrt.clean("OG0002754", "Yellow-e3", ort.24.dds)

#####OG0002872: CD109 antigen
ann.slim[ann.slim$Orthogroup == "OG0002872",]

visiOrt.messy("OG0002872", "CD109 antigen", ort.24.dds)

visiOrt.clean("OG0002872", "CD109 antigen", ort.24.dds)
#####OG0003143: J domain-containing protein
ann.slim[ann.slim$Orthogroup == "OG0003143",]
visiOrt.messy("OG0003143", "OG0003143", ort.24.dds)

visiOrt.clean("OG0003143", "OG0003143", ort.24.dds)
#####OG0003297: M1Pi

#M1Pi / methylthioribose-1-phosphate isomerase

ann.slim[ann.slim$Orthogroup == "OG0003297",]

visiOrt.messy("OG0003297", "M1Pi", ort.24.dds)

visiOrt.clean("OG0002754", "M1Pi", ort.24.dds)
#####OG0003741: Manganese-transporting ATPase 13A1
ann.slim[ann.slim$Orthogroup == "OG0003741",]
visiOrt.messy("OG0003741", "OG0003741", ort.24.dds)

visiOrt.clean("OG0003741", "OG0003741", ort.24.dds)
#####OG0003982: Leukocyte elastase inhibitor
ann.slim[ann.slim$Orthogroup == "OG0003982",]
visiOrt.messy("OG0003982", "Leukocyte elastase inhibitor", ort.24.dds)

visiOrt.clean("OG0003982", "Leukocyte elastase inhibitor", ort.24.dds)
#####OG0004924: Ninjurin-1
ann.slim[ann.slim$Orthogroup == "OG0004924",]
visiOrt.messy("OG0004924", "Ninjurin-1", ort.24.dds)

visiOrt.clean("OG0004924", "Ninjurin-1", ort.24.dds)

#####OG0005280: Art3
ann.slim[ann.slim$Orthogroup == "OG0005280",]
visiOrt.messy("OG0005280", "Art3", ort.24.dds)

visiOrt.clean("OG0005280", "Art3", ort.24.dds)

#####OG0007008: NHE-RF1
ann.slim[ann.slim$Orthogroup == "OG0007008",]
visiOrt.messy("OG0007008", "NHE-RF1", ort.24.dds)

visiOrt.clean("OG0007008", "NHE-RF1", ort.24.dds)
#####OG0007879: Collagen alpha-2(IV) chain-like
ann.slim[ann.slim$Orthogroup == "OG0007879",]
visiOrt.messy("OG0007879", "OG0007879", ort.24.dds)

visiOrt.clean("OG0007879", "OG0007879", ort.24.dds)
#####OG0008078: threonine--tRNA ligase, cytoplasmic
ann.slim[ann.slim$Orthogroup == "OG0008078",]
visiOrt.messy("OG0008078", "OG0008078", ort.24.dds)

visiOrt.clean("OG0008078", "OG0008078", ort.24.dds)
#####OG0008153: Fibroin heavy chain
ann.slim[ann.slim$Orthogroup == "OG0008153",]
visiOrt.messy("OG0008153", "Fibroin heavy chain", ort.24.dds)

visiOrt.clean("OG0008153", "Fibroin heavy chain", ort.24.dds)

#####OG0008945: Cytochrome b5
ann.slim[ann.slim$Orthogroup == "OG0002754",]
visiOrt.messy("OG0008945", "Cytochrome b5", ort.24.dds)

visiOrt.clean("OG0008945", "Cytochrome b5", ort.24.dds)
#####OG0009885: Uncharacterised OG0009885
ann.slim[ann.slim$Orthogroup == "OG0009885",]
visiOrt.messy("OG0009885", "Uncharacterised OG0009885", ort.24.dds)

visiOrt.clean("OG0009885", "Uncharacterised OG0009885", ort.24.dds)

####Downregulated Gram Pos

#That took forever so I'm only gonna label immune / putative immune genes going forwards

####Upregulated Gram Neg

#That took forever so I'm only gonna label immune / putative immune genes going forwards
orts.up.neg <- ort.24.res.sig$OrthoGroup[ort.24.res.sig$log2FoldChange > 0 &
                                           ort.24.res.sig$Contrast == "Treatment_GramNeg_vs_Naive"]

orts.up.neg[orts.up.neg %in% orts.up.pos]
orts.up.neg <- orts.up.neg[!orts.up.neg %in% orts.up.pos]
orts.up.neg

ann[ann$Orthogroup %in% orts.up.neg,] %>%
  arrange(Function)

#Easy.

labs <- orts.up.neg
labs[labs=="OG0000244"] <- "Uncharacterised LOC100577430"
labs[labs=="OG0001768"] <- "HSDL1"
labs[labs=="OG0009073"] <- "GCWSP"
labs[labs=="OG0000446"] <- "Antitrypsin"
labs

for (i in 1:length(orts.up.neg)){
  visiOrt.messy(orts.up.neg[i], labs[i], ort.24.dds)
  visiOrt.clean(orts.up.neg[i], labs[i], ort.24.dds)
}

##OrthoLevel DE: All Age-Controlled Species (AgeCont) all samples

###Preparing TX Object

runs <- c("Amel", "Bter", "Plan")
filterCount(runs, "AgeControlled_aS")

wd <- "output/Ortho_DE/AgeControlled_aS/"

ort.age.sa.txs <- ortCount("AgeControlled_aS", tx.list)

#make combined abundance matrix
abundance <- cbind(ort.age.sa.txs[[1]]$abundance,
                   ort.age.sa.txs[[2]]$abundance,
                   ort.age.sa.txs[[3]]$abundance)
nrow(abundance)

#make combined counts matrix
counts <- cbind(ort.age.sa.txs[[1]]$counts,
                ort.age.sa.txs[[2]]$counts,
                ort.age.sa.txs[[3]]$counts)
nrow(counts)

#make combined lengths matrix
length <- cbind(ort.age.sa.txs[[1]]$length,
                ort.age.sa.txs[[2]]$length,
                ort.age.sa.txs[[3]]$length)
nrow(length)

ort.age.sa.tx <- list(abundance = abundance,
                      counts = counts, 
                      length = length)

ort.age.sa.tx$countsFromAbundance <- "no"

for (i in 1:length(ort.age.sa.tx)){
  print(head(ort.age.sa.tx[[i]], n = 1))
}

remove(ort.age.sa.txs)

###DE Analysis

ort.age.sa.dds <- getDDS("AgeControlled_aS", ort.age.sa.tx)

ort.age.sa.res <- vector(mode = "list", length = 3)
treats <- resultsNames(ort.age.sa.dds)[c(4,3,2)]

for(i in 1:length(ort.age.sa.res)){
  ort.age.sa.res[[i]] <- results(ort.age.sa.dds, name = paste(treats[i]))
  print(treats[i])
  summary(ort.age.sa.res[[i]])
}

###Results

for (i in 1:length(ort.age.sa.res)){
  ort.age.sa.res[[i]]$Contrast <- paste(treats[i])
  ort.age.sa.res[[i]]$OrthoGroup <- rownames(ort.age.sa.res[[i]])
  rownames(ort.age.sa.res[[i]]) <- NULL
  ort.age.sa.res[[i]] <- as.data.frame(ort.age.sa.res[[i]])
}

ort.age.sa.res.all <- bind_rows(ort.age.sa.res[1:3])

write.table(ort.age.sa.res.all, paste(wd, 
                                      "AgeControlled_aS_DeSeqResults_All.tsv", sep = ""),
            col.names = T, row.names = F, quote = F, sep = "\t")

ort.age.sa.res.sig <- subset(ort.age.sa.res.all, padj < 0.1)

write.table(ort.age.sa.res.sig, paste(wd,
                                      "AgeControlled_aS_DeSeqResults_FDRlt0.1.tsv", sep =""),
            col.names = T, row.names = F, quote = F, sep = "\t")

###ExpBoxPlots

ort.age.sa.res.sig

sig.ort <- ort.age.sa.res.sig$OrthoGroup
ann[ann$Orthogroup %in% sig.ort,] %>%
  arrange(Orthogroup)

ortho[ortho$Orthogroup =="OG0006713",]

####OG0003982: leukocyte elastase inhibitor (serpin)

visiOrt.messy("OG0003982", "Leukocyte elastase inhibitor", ort.age.sa.dds)

visiOrt.clean("OG0003982", "Leukocyte elastase inhibitor", ort.age.sa.dds)

####OG0005698: protein lozenge

visiOrt.messy("OG0005698", "protein lozenge", ort.age.sa.dds)

visiOrt.clean("OG0005698", "protein lozenge", ort.age.sa.dds)

####OG0006713: TSR1 homolog

visiOrt.messy("OG0000699", "TSR1 homolog", ort.age.sa.dds)

visiOrt.clean("OG0000699", "TSR1 homolog", ort.age.sa.dds)


##SessionInfo####

sessionInfo()

