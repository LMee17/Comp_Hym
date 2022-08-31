#Aug 2022
#Individual species' gene-level DE analyses

#These have been converted from .rmd files

set.seed(42)
dir.create("output")
dir.create("output/Ind_DE/")

##Libraries####

```{r}
library("DESeq2")
library("tidyverse")
library("ggrepel")
library("tximport")
library("pheatmap")
library("reshape")
library("ggvenn")
```

##Resources####

#OrthoGroup information / Annotation information

ortho <- read.table("input/Ortho_Files/OrthoFinder_Combined_Parsed_Aug22.tsv",
                    sep = "\t", header = T)

ann <- read.table("input/Gene_Info/Amel_ImmAnnotations_Orthogroup_Aug22.tsv",
                  header = T, sep = "\t")
ann.slim <- ann[,c(2:3)]

##Functions####

#Converts transcript-level counts from input/ directory to a single counts table 
#with gene level estimates. Also generates a coldata object in the output directory

geneCount <- function(species){
  dir <- paste("input/Counts_Kal/", species, sep = "")
  files <- list.files(dir)
  samples <- sapply(strsplit(files, ".", fixed = T), "[", 1)
  treat <- c(rep("Naive", 3), rep("Wound", 3), rep("GramPos", 3),
             rep("GramNeg",3))
  coldata <- as.data.frame(cbind(samples, treat))
  names(coldata) <- c("Sample", "Treatment")
  write.table(coldata, 
              paste("output/Ind_DE/", species, "/", species, "_SampleData.tsv", sep = ""),
              col.names = T, row.names = F, sep = "\t", quote = F)
  names(files) <- samples
  txgene <- read.table(paste("input/Gene_Info/", species, "_TXGene.tsv", sep =""),
                       header = T, sep = "\t")
  cnts <- file.path(dir, files)
  names(cnts) <- samples
  txi.kal.tsv <- tximport(cnts, type = "kallisto", tx2gene = txgene, ignoreAfterBar = T)
  write.table(txi.kal.tsv$counts,
              paste("output/Ind_DE/", species, "/", 
                    species,  "_TxGeneCount.tsv", sep = ""),
              col.names = T, row.names = T, quote = F, sep = "\t")
  return(txi.kal.tsv)
}

#Produce DDS object from tximport object (includes prefiltering and setting 
#reference condition to Naive)

getDDS <- function(species, txobject){
  cd <- read.table(paste("output/Ind_DE/", species, "/",
                         species, "_SampleData.tsv", sep = ""), header = T, sep = "\t")
  ddsTx <- DESeqDataSetFromTximport(txobject, colData = cd, design = ~Treatment)
  keep <- rowSums(counts(ddsTx)) >= 10
  ddsTx2 <- ddsTx[keep,]
  ddsTx2$Treatment <- relevel(ddsTx2$Treatment, ref = "Naive")
  ddsTx2 <- DESeq(ddsTx2)
  return(ddsTx2)
}

#Visually assess the samples if (an) outlier(s) was/were removed

outCheck <- function(outlier_s, species){
  cd <- read.table(paste("output/Ind_DE/", species, "/",
                         species, "_SampleData.tsv", sep = ""), header = T)
  cd <- cd[!cd$Sample %in% outlier_s,]
  dir <- paste("input/Counts_Kal/", species, sep = "")
  files <- list.files(dir)
  samples <- sapply(strsplit(files, ".", fixed = T), "[", 1)
  txgene <- read.table(paste("input/Gene_Info/", species, "_TXGene.tsv", sep =""),
                       header = T, sep = "\t")
  cnts <- file.path(dir, files)
  names(cnts) <- samples
  cnts[outlier_s]
  cnts <- cnts[!names(cnts) %in% outlier_s]
  tx <- tximport(cnts, type = "kallisto", tx2gene = txgene, ignoreAfterBar = T)
  
  ddsTx <- DESeqDataSetFromTximport(tx, colData = cd, design = ~Treatment)
  keep <- rowSums(counts(ddsTx)) >= 10
  ddsTx2 <- ddsTx[keep,]
  ddsTx2$Treatment <- relevel(ddsTx2$Treatment, ref = "Naive")
  ddsTx2 <- DESeq(ddsTx2)
  vsd <- vst(ddsTx2, blind = F)
  select <- order(rowMeans(counts(ddsTx2,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  p1 <- pheatmap(assay(vsd)[select,])
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$type, sep="-")
  colnames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$type, sep="-")
  p2 <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists)
  pca <- plotPCA(vsd, intgroup="Treatment", returnData = TRUE)
  p3 <- ggplot(pca, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
    geom_point() +
    geom_text()
  p1
  p2 
  p3
}

#Return a dds object with specified outlier removed

ddsOmit <- function(outlier_s, species){
  cd <- read.table(paste("output/Ind_DE/", species, "/",
                         species, "_SampleData.tsv", sep = ""), header = T)
  cd <- cd[!cd$Sample %in% outlier_s,]
  dir <- paste("input/Counts_Kal/", species, sep = "")
  files <- list.files(dir)
  samples <- sapply(strsplit(files, ".", fixed = T), "[", 1)
  txgene <- read.table(paste("input/Gene_Info/", species, "_TXGene.tsv", sep =""),
                       header = T, sep = "\t")
  cnts <- file.path(dir, files)
  names(cnts) <- samples
  cnts[outlier_s]
  cnts <- cnts[!names(cnts) %in% outlier_s]
  tx <- tximport(cnts, type = "kallisto", tx2gene = txgene, ignoreAfterBar = T)
  
  ddsTx <- DESeqDataSetFromTximport(tx, colData = cd, design = ~Treatment)
  keep <- rowSums(counts(ddsTx)) >= 10
  ddsTx2 <- ddsTx[keep,]
  ddsTx2$Treatment <- relevel(ddsTx2$Treatment, ref = "Naive")
  ddsTx2 <- DESeq(ddsTx2)
  return(ddsTx2)
}

#Produce a dataframe with orthogroup and gene class information ready 
#for data visualisation

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
  tmp <- counts(dds, normalized = T)
  samples <- dds$Sample
  counts <- tmp[rownames(tmp) == paste(gene),]
  plot.df <- as.data.frame(cbind(samples, counts))
  names(plot.df) <- c("Samples", "Norm_Counts")
  plot.df$Treatment[grepl("_N_", plot.df$Samples)] <- "Naive"
  plot.df$Treatment[grepl("_P_", plot.df$Samples)] <- "Wound"
  plot.df$Treatment[grepl("_SM_", plot.df$Samples)] <- "Gram Negative"
  plot.df$Treatment[grepl("_SL_", plot.df$Samples)] <- "Gram Positive"
  plot.df$Treatment <- factor(plot.df$Treatment, levels = c("Naive",
                                                            "Wound", 
                                                            "Gram Positive", 
                                                            "Gram Negative"))
  plot.df$Norm_Counts <- as.numeric(plot.df$Norm_Counts)
  p <- ggplot(plot.df, aes(x = Treatment, y = log(Norm_Counts), fill = Treatment)) +
    geom_boxplot() +
    geom_point(position = position_jitter(w=0.1, h = 0),
               aes(colour = Treatment), size = 3) +
    scale_colour_manual(values = c("#58ffa7", "#ffa758", "#ff58b0", "#58b0ff")) +
    xlab("Treatment") +
    labs(title = (paste(genename, " expression across treatments", sep = ""))) +
    guides(alpha = FALSE, colour = guide_legend("Treatment"))
  ggsave(paste(wd, "Plots/ExpressionPlots/", gene, "_ExpBoxPlot.pdf", sep = ""))
  p
}

#Plot all expression boxplots from a group of genes 

visiFacet <- function(genes, genenames, dds, filename){
  p <- vector(mode = "list", length = length(genes))
  for (i in 1:length(p)){
    tmp <- counts(dds, normalized = T)
    samples <- dds$Sample
    counts <- tmp[rownames(tmp) == paste(genes[i]),]
    p[[i]] <- as.data.frame(cbind(samples, counts))
    names(p[[i]]) <- c("Samples", "Norm_Counts")
    p[[i]]$Treatment[grepl("_N_", p[[i]]$Samples)] <- "Naive"
    p[[i]]$Treatment[grepl("_P_", p[[i]]$Samples)] <- "Wound"
    p[[i]]$Treatment[grepl("_SM_", p[[i]]$Samples)] <- "Gram Negative"
    p[[i]]$Treatment[grepl("_SL_", p[[i]]$Samples)] <- "Gram Positive"
    p[[i]]$Treatment <- factor(p[[i]]$Treatment, levels = c("Naive",
                                                            "Wound", 
                                                            "Gram Positive", 
                                                            "Gram Negative"))
    p[[i]]$Norm_Counts <- as.numeric(p[[i]]$Norm_Counts)
    p[[i]]$Gene <- paste(genes[i])
    p[[i]]$Label <- paste(genenames[i])
  }
  p2 <- bind_rows(p)
  p3 <- ggplot(p2, aes(x = Treatment, y = log(Norm_Counts))) +
    geom_boxplot(aes(colour = Treatment, alpha = 0.5)) +
    geom_point(position = position_jitter(w=0.1, h = 0),
               aes(colour = Treatment), size = 3) +
    scale_colour_manual(values = c("#58ffa7", "#ffa758", "#ff58b0", "#58b0ff")) +
    xlab("Treatment") +
    guides(alpha = FALSE, colour = guide_legend("Treatment")) +
    facet_wrap(~Label, nrow = 5, scale = "free_y")
  ggsave(paste(wd, "Plots/ExpressionPlots/", filename, "_Facet_ExpBoxPlot.pdf", sep = ""))
  p3
}







##Apis mellifera####

dir.create("output/Ind_DE/Amel/")
wd <- "output/Ind_DE/Amel/"

###DE Analysis

#First, convert transcript-level counts table from Kallisto to gene level 
#counts and produce a coldata object

amel.tx <- geneCount("Amel")
head(amel.tx$counts)

#Build a DESeq DataSet from the above, prefilter (keep only rows that have at least
#10 reads total), and return dds object.

amel.dds <- getDDS("Amel", amel.tx)
amel.res <- vector(mode = "list", length = 3)
treats <- resultsNames(amel.dds)[c(4,3,2)]

for(i in 1:length(amel.res)){
  amel.res[[i]] <- results(amel.dds, name = paste(treats[i]))
  summary(amel.res[[i]])
}

###Check for outliers
#First, check counts

boxplot(log10(assays(amel.dds)[["cooks"]]), range=0, las=2)

#normalise
vsd <- vst(amel.dds, blind = F)
#make heatmap (top 20 genes)
select <- order(rowMeans(counts(amel.dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
pheatmap(assay(vsd)[select,])

pdf("output/Ind_DE/Amel/Amel_Pheatmap.pdf")
  pheatmap(assay(vsd)[select,])
dev.off()

#AM_N_056 is clustering with bacterial and wound challenges... which may be a problem.
#Check with sample to sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)


pdf("output/Ind_DE/Amel/Amel_DistanceMatrix.pdf")
  pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()

#Finally, PCA

amel.pca <- plotPCA(vsd, intgroup="Treatment", returnData = TRUE)
ggplot(amel.pca, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
  geom_point() +
  geom_text()

ggsave("output/Ind_DE/Amel/Amel_PCA.pdf")

#What happens if I remove the potential outlier ? 

torem <- "AM_N_056"
outCheck(torem, "Amel")

#I'm not convinced that I need to remove that sample, but I'll look at the results 
#regardless.

amel.dds2 <- ddsOmit(torem, "Amel")

amel.res2 <- vector(mode = "list", length = 3)
treats2 <- resultsNames(amel.dds2)[c(4,3,2)]

for(i in 1:length(amel.res2)){
  amel.res2[[i]] <- results(amel.dds2, name = paste(treats2[i]))
  print(treats2[i])
  summary(amel.res2[[i]])
}

#Save plots
vsd2 <- vst(amel.dds2, blind = F)
#make heatmap (top 20 genes)
select <- order(rowMeans(counts(amel.dds2,normalized=TRUE)),
                decreasing=TRUE)[1:20]
pheatmap(assay(vsd2)[select,])

pdf("output/Ind_DE/Amel/Amel_OutlierRemoved_Pheatmap.pdf")
  pheatmap(assay(vsd2)[select,])
dev.off()

sampleDists2 <- dist(t(assay(vsd2)))
sampleDistMatrix2 <- as.matrix(sampleDists2)
rownames(sampleDistMatrix2) <- paste(vsd2$Treatment, vsd2$type, sep="-")
colnames(sampleDistMatrix2) <- NULL

pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)

pdf("output/Ind_DE/Amel/Amel_OutlierRemoved_DistanceMatrix.pdf")
  pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)
dev.off()

amel.pca2 <- plotPCA(vsd2, intgroup="Treatment", returnData = TRUE)
ggplot(amel.pca2, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
  geom_point() +
  geom_text()
ggsave("output/Ind_DE/Amel/Amel_OutlierRemoved_PCA.pdf")

#Gonna leave it as it is and not remove any outliers.

###Results
#Save 

for (i in 1:length(amel.res)){
  amel.res[[i]]$Contrast <- paste(treats[i])
  amel.res[[i]]$Gene <- rownames(amel.res[[i]])
  rownames(amel.res[[i]]) <- NULL
  amel.res[[i]] <- as.data.frame(amel.res[[i]])
}
 
amel.res.all <- bind_rows(amel.res[1:3])
 
id <- read.table("input/Gene_Info/Amel_IDwDesc.tsv", 
                 sep = "\t", quote = "", header = T)

for (i in 1:nrow(amel.res.all)){
  d <- unique(id$Description[id$GeneID == amel.res.all$Gene[i]])
  if(length(d) > 1){
    amel.res.all$Desc[i] <- d[1]
  } else {
    amel.res.all$Desc[i] <- d  
  }
}

write.table(amel.res.all, paste(wd, "Amel_DeSeqResults_All.tsv", sep = ""),
            col.names = T, row.names = F, quote = F, sep = "\t")

amel.res.sig <- subset(amel.res.all, padj < 0.1)

write.table(amel.res.sig, paste(wd, "Amel_DeSeqResults_FDRlt0.1.tsv", sep =""), 
            col.names = T, row.names = F, quote = F, sep = "\t")

#Genes of interest?
#I'm going to take the top 10 most significant and top 10 biggest logFC per treatment. 
#Genes that overlap between two (and between the three treatments) will 
#be of particular interest

nona <- amel.res.all[!is.na(amel.res.all$padj),]
nona$LogDif <- sub('-', '', nona$log2FoldChange)
tmp <- nona[nona$Contrast == treats[1],]
sig.wound <- head(tmp$Gene[order(tmp$padj)], n = 10)
top.wound <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)
tmp <- nona[nona$Contrast == treats[2],]
sig.pos <- head(tmp$Gene[order(tmp$padj)], n = 10)
top.pos <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)
tmp <- nona[nona$Contrast == treats[3],]
sig.neg <-  head(tmp$Gene[order(tmp$padj)], n = 10)
top.neg <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)

goi <- c(sig.wound, top.wound, sig.pos, top.pos, sig.neg, top.neg)
stat <- c(rep("Top 10 Sig Wound", 10), rep("Top 10 log2FC Wound", 10),
          rep("Top 10 Sig Gram Positive", 10), rep("Top 10 log2FC Gram Positive", 10),
          rep("Top 10 Sig Gram Negative", 10), rep("Top 10 log2FC Gram Negative", 10))

goi.df <- as.data.frame(cbind(goi, stat))

goi.wound <- unique(c(sig.wound, top.wound))
goi.pos <- unique(c(sig.pos, top.pos))
goi.neg <- unique(c(sig.neg, top.neg))

###Visualise: Venn Diagrams

top.sig <- list("Wound" = sig.wound,
                "Gram Positive" = sig.pos,
                "Gram Negative" = sig.neg)

p <- ggvenn(top.sig,
            set_name_size = 5,
            show_percentage = F,
            fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Top 10 Genes per Treatment by adjPvalue")

dir.create(paste(wd, "Plots/", sep = ""))
ggsave(paste(wd, "Plots/Amel_TopAdjPvalue_Venn.pdf", sep = ""))

top.lfc <- list("Wound" = top.wound,
                "Gram Positive" = top.pos,
                "Gram Negative" = top.neg)

p <- ggvenn(top.lfc,
            set_name_size = 5,
            show_percentage = F,
            fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Top 10 Genes per Treatment by log2FC")

ggsave(paste(wd, "Plots/Amel_Toplog2FC_Venn.pdf", sep = ""))

goi <- list("Wound" = goi.wound,
            "Gram Positive" = goi.pos,
            "Gram Negative" = goi.neg)

p <- ggvenn(goi,
            set_name_size = 5,
            show_percentage = F,
            fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Genes of Interest Across Treatments")

ggsave(paste(wd, "Plots/Amel_GOI_Venn.pdf", sep = ""))

goi.mat <- makeMatrix(goi.df)
goi.mat

###Visualise: Expression BoxPlots

dir.create(paste(wd, "Plots/ExpressionPlots/", sep = ""))

####Significant Wound

goi.mat$Gene[goi.mat$SigWound == 1]
genes <- goi.mat$Gene[goi.mat$SigWound == 1]

visiFacet(genes = genes,
          genenames = genes,
          dds = amel.dds, 
          filename = "Woundlog2FC")

#All Wound: 
visiFacet(genes = genes,
          genenames = genes,
          dds = amel.dds, 
          filename = "WoundSig")

####TopFC Wound

goi.mat$Gene[goi.mat$TopWound == 1]
genes <- goi.mat$Gene[goi.mat$TopWound == 1]

#All Wound: 

visiFacet(genes = genes,
          genenames = genes,
          dds = amel.dds, 
          filename = "Woundlog2FC")

####Significant GramPos

goi.mat$Gene[goi.mat$SigPos == 1]
genes <- goi.mat$Gene[goi.mat$SigPos == 1]

#All Gram Pos Sig: 
visiFacet(genes = genes,
          genenames = genes,
          dds = amel.dds, 
          filename = "PosSig")

####TopFC GramPos
goi.mat$Gene[goi.mat$TopPos == 1]
genes <- goi.mat$Gene[goi.mat$TopPos == 1]

#All log Gram Pos: 

visiFacet(genes = genes,
          genenames = genes,
          dds = amel.dds, 
          filename = "Poslog2FC")

####Significant GramNeg

goi.mat$Gene[goi.mat$SigNeg == 1]
genes <- goi.mat$Gene[goi.mat$SigNeg == 1]

#All Gram Neg Sig: 

visiFacet(genes = genes,
          genenames = genes,
          dds = amel.dds, 
          filename = "NegSig")

####TopFC GramPos

goi.mat$Gene[goi.mat$TopNeg == 1]
genes <- goi.mat$Gene[goi.mat$TopNeg == 1]

###Visualise: PCA
#Decided to not remove potential outliers.

#Resultant PCA:
ggplot(amel.pca, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() 

ggsave(paste(wd, "Plots/Amel_PCA_Simple.pdf", sep = ""))

ggplot(amel.pca, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() +
  geom_label_repel(
    data = amel.pca,
    aes(label = name),
    size = 3
  )

ggsave(paste(wd, "Plots/Amel_PCA_Labelled.pdf", sep = ""))

amel.pca$Treat2 <- c(rep("Naive", 3), rep("Treatment", 9))

ggplot(amel.pca, aes(x = PC1, y = PC2, color = Treat2)) +
  geom_point() +
  stat_ellipse()

ggsave(paste(wd, "Plots/Amel_PCA_TrtVUnTrt.pdf", sep = ""))

###Visualise: Smears

#Amel/Bter consideration : there is a shared gene name that causes issues. To fix:
amel.res.all$Gene[amel.res.all$Gene == "Kr-h1"] <- "Kr-h1_Amel"

ortho$GeneID[ortho$GeneID == "Kr-h1" & ortho$Species == "Amel"] <- "Kr-h1_Amel"

plot.df <- makeMetaDf(amel.res.all)
plot.df

#The function assigns the genes without an orthogroup "unknown" functions. 
#This is fine for the non Amel species but for Amel I can assign function based on 
#gene ID

for (i in 1:nrow(plot.df)){
  if (plot.df$Function[i] == "Unknown"){
    plot.df$Function[i] <- ann$Function[ann$GeneID == plot.df$Gene[i]]
  }
}

plot.df$Function <- gsub("Effector", "Immune", 
                         gsub("Recognition", "Immune", 
                              gsub("Signalling", "Immune", 
                                   gsub("Background", "None Immune", plot.df$Function))))

tail(plot.df)

plot.df$Condition <- factor(plot.df$Condition, levels = c("Wound",
                                                          "Gram Positive",
                                                          "Gram Negative"))

ggplot(plot.df, aes(x = log(baseMean), y = log2FoldChange)) +
  geom_point(aes(colour = Sig)) +
  facet_grid(~Condition) +
  scale_color_manual(values=c("black", "Red"), guide = NULL) +
  theme(legend.position = "bottom")

#Note: remember only those that passed baseMean threshold per comparison are plotted
#(120 for Gram Pos, considerably lower for the others)

ggsave(paste(wd, "Plots/Amel_SmearPlot.pdf", sep = ""))

###Visualise: Rainclouds

sig.genes <- unique(unlist(sig))
plot.sig <- plot.df[plot.df$Gene %in% sig.genes,]

to.add.naive <- plot.sig[plot.sig$Gene %in% sig.genes &
                           plot.sig$Condition == "Wound",]
to.add.naive$Sig <- "No"
to.add.naive$Condition <- "Naive"
to.add.naive$log2FoldChange <- 0

plot.sig <- rbind(plot.sig, to.add.naive)

plot.sig$Condition <- factor(plot.sig$Condition, levels = c("Naive",
                                                            "Wound", 
                                                            "Gram Positive", 
                                                            "Gram Negative"))

plot.sig$Significance[plot.sig$padj > 0.1] <- "Not Significant"
plot.sig$Significance[plot.sig$padj < 0.1 & plot.sig$padj > 0.05] <- 
  "Significant (FDR < 0.1)"
plot.sig$Significance[plot.sig$padj < 0.05 & plot.sig$padj > 0.01] <-
  "Significant (FDR < 0.05)"
plot.sig$Significance[plot.sig$padj < 0.01 & plot.sig$padj > 0.001] <-
  "Significant (FDR < 0.01)"
plot.sig$Significance[plot.sig$padj < 0.001] <-
  "Significant (FDR < 0.001)"

plot.sig$Significance <- factor(plot.sig$Significance, levels =
                                  c("Not Significant",
                                    "Significant (FDR < 0.1)",
                                    "Significant (FDR < 0.05)",
                                    "Significant (FDR < 0.01)",
                                    "Significant (FDR < 0.001)"))

plot.sig$Function <- factor(plot.sig$Function, levels = c("Immune",
                                                          "Putative Immune", 
                                                          "None Immune"))

ggplot(plot.sig[!plot.sig$Condition == "Naive",], 
       (aes(x = Condition, y = log2FoldChange))) +
  geom_jitter(aes(colour = Function, shape = Significance)) +
  scale_colour_manual(values = c("#2c8cd9", "#7a2cd9", "darkgrey")) +
  scale_shape_manual(values = c(1, 5, 18, 19, 17)) 

ggsave(paste(wd, "Plots/Amel_logFCJitter.pdf", sep = ""))

ggplot(plot.sig, 
       (aes(x = Condition, y = log2FoldChange))) +
  geom_point(aes(colour = Function, shape = Significance)) +
  geom_line( linetype = "dashed", alpha = 0.2, aes(group = Gene, colour = Function)) +
  scale_colour_manual(values = c("#2c8cd9", "#7a2cd9", "darkgrey")) +
  scale_shape_manual(values = c(1, 5, 18, 19, 17)) 

ggsave(paste(wd, "Plots/Amel_logFCScatter_GeneLine.pdf", sep = ""))

###Visualise: Volcano Plots

plot.df$Function <- factor(plot.df$Function, levels = c("Immune",
                                                        "Putative Immune",
                                                        "None Immune"))

ggplot(plot.df,
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
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Amel_Simple_FDR0.1_Volcano.pdf", sep = ""))

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1)) +
  guides(colour = guide_legend("FDR < 0.1", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Amel_Class_FDR0.1_Volcano.pdf", sep = ""))

plot.df$Direction2 <- "Not Significant"
plot.df$Direction2[plot.df$padj < 0.05 &
                     plot.df$log2FoldChange > 0 ] <- "Up-Regulated"
plot.df$Direction2[plot.df$padj < 0.05 &
                     plot.df$log2FoldChange < 0 ] <- "Down-Regulated"

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction2), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.05", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Amel_Simple_FDR0.05_Volcano.pdf", sep = ""))

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction2, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1)) +
  guides(colour = guide_legend("FDR < 0.05", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Amel_Class_FDR0.05_Volcano.pdf", sep = ""))

plot.df$Direction3 <- "Not Significant"
plot.df$Direction3[plot.df$padj < 0.001 &
                     plot.df$log2FoldChange > 0 ] <- "Up-Regulated"
plot.df$Direction3[plot.df$padj < 0.001 &
                     plot.df$log2FoldChange < 0 ] <- "Down-Regulated"

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Amel_Simple_FDR0.001_Volcano.pdf", sep = ""))

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1)) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Amel_Class_FDR0.001_Volcano.pdf", sep = ""))

#Label top 3 sig genes per treatment

sig.wound[1:3]
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC724654"] <- "Cytochrome b5"
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC102654628"] <- "Uncharacterised LOC102654628"
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC113218757"] <- "LOC113218757"

sig.pos[1:3]
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC726750"] <- "Fibroin heavy chain"
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC411012"] <- "Cactus2"
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC724654"] <- "Cytochrome b5"

sig.neg[1:3]
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC726750"] <- "Fibroin heavy chain"
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC411012"] <- "Cytochrome b5"
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC411807"] <- "Uncharacterised LOC411807"

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") +
  geom_label_repel(
    data = plot.df[!is.na(plot.df$Label),],
    size = 2,
    aes(label = Label),
    point.padding = unit(0.3, "lines"),
    box.padding = 1.25
  ) 


ggsave(paste(wd, "Plots/Amel_Simple_FDR0.001_Volcano_top3Label.pdf", sep = ""))

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1)) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") +
  geom_label_repel(
    data = plot.df[!is.na(plot.df$Label),],
    size = 2,
    aes(label = Label),
    point.padding = unit(0.3, "lines"),
    box.padding = 1.2
  ) 


ggsave(paste(wd, "Plots/Amel_Class_FDR0.001_Volcano_top3Label.pdf", sep = ""))

##Bombus terrestris####

dir.create("output/Ind_DE/Bter/")
wd <- "output/Ind_DE/Bter/"

###DE Analysis
#First, convert transcript-level counts table from Kallisto to gene level counts and 
#produce a coldata object

bter.tx <- geneCount("Bter")
head(bter.tx$counts)

#Build a DESeq DataSet from the above, prefilter (keep only rows that have at least 10 
#reads total), and return dds object.

bter.dds <- getDDS("Bter", bter.tx)
bter.res <- vector(mode = "list", length = 3)
treats <- resultsNames(bter.dds)[c(4,3,2)]

for(i in 1:length(bter.res)){
  bter.res[[i]] <- results(bter.dds, name = paste(treats[i]))
  summary(bter.res[[i]])
}

###Check for outliers

boxplot(log10(assays(bter.dds)[["cooks"]]), range=0, las=2)

#normalise
vsd <- vst(bter.dds, blind = F)
#make heatmap (top 20 genes)
select <- order(rowMeans(counts(bter.dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
pheatmap(assay(vsd)[select,])

#From this I would suggest P30 and SM32 are outliers - and possibly N001

pdf("output/Ind_DE/Bter/Bter_Pheatmap.pdf")
  pheatmap(assay(vsd)[select,])
dev.off()

#Check with sample to sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

pdf("output/Ind_DE/Bter/Bter_DistanceMatrix.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()

#Yeah that's pretty messy
#Finally, PCA

bter.pca <- plotPCA(vsd, intgroup="Treatment", returnData = TRUE)
ggplot(bter.pca, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
  geom_point() +
  geom_text()

ggsave("output/Ind_DE/Bter/Bter_PCA.pdf")

#What happens if I remove the potential outliers ? 
torem <- c("BT_A_SL_003", "BT_A_P_030", "BT_A_SM_032", "BT_A_N_001")
outCheck(torem, "Bter")

#2 plots are still a mess, but the sample-distance plots are very much improved.
#Now SL 003 is a problem ? Go back and try again....

#Much better.... though I've lost 4 samples!

bter.dds2 <- ddsOmit(torem, "Bter")

bter.res2 <- vector(mode = "list", length = 3)
treats2 <- resultsNames(bter.dds2)[c(4,3,2)]

for(i in 1:length(bter.res2)){
  bter.res2[[i]] <- results(bter.dds2, name = paste(treats2[i]))
  print(treats2[i])
  summary(bter.res2[[i]])
}

#Save plots

vsd2 <- vst(bter.dds2, blind = F)
#make heatmap (top 20 genes)
select <- order(rowMeans(counts(bter.dds2,normalized=TRUE)),
                decreasing=TRUE)[1:20]
pheatmap(assay(vsd2)[select,])

pdf("output/Ind_DE/Bter/Bter_OutlierRemoved_Pheatmap.pdf")
  pheatmap(assay(vsd2)[select,])
dev.off()

sampleDists2 <- dist(t(assay(vsd2)))
sampleDistMatrix2 <- as.matrix(sampleDists2)
rownames(sampleDistMatrix2) <- paste(vsd2$Treatment, vsd2$type, sep="-")
colnames(sampleDistMatrix2) <- NULL

pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)

pdf("output/Ind_DE/Bter/Bter_OutlierRemoved_DistanceMatrix.pdf")
pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)
dev.off()

bter.pca2 <- plotPCA(vsd2, intgroup="Treatment", returnData = TRUE)
ggplot(bter.pca2, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
  geom_point() +
  geom_text()
ggsave("output/Ind_DE/Bter/Bter_OutlierRemoved_PCA.pdf")
#gonna make the call to go forward without the outliers 

bter.dds <- bter.dds2
bter.res <- bter.res2

###Results
#Save 
#(Once this is done once, read from the produced file, annotating takes some time.)

for (i in 1:length(bter.res)){
  bter.res[[i]]$Contrast <- paste(treats[i])
  bter.res[[i]]$Gene <- rownames(bter.res[[i]])
  rownames(bter.res[[i]]) <- NULL
  bter.res[[i]] <- as.data.frame(bter.res[[i]])
}

bter.res.all <- bind_rows(bter.res[1:3])

id <- read.table("input/Gene_Info/Bter_IDwDesc.tsv", 
                 sep = "\t", quote = "", header = T)


for (i in 1:nrow(bter.res.all)){
  d <- unique(id$Description[id$GeneID == bter.res.all$Gene[i]])
  if(length(d) > 1){
    bter.res.all$Desc[i] <- d[1]
  } else {
    bter.res.all$Desc[i] <- d  
  }
}

write.table(bter.res.all, paste(wd, "Bter_DeSeqResults_All.tsv", sep = ""),
            col.names = T, row.names = F, quote = F, sep = "\t")

bter.res.sig <- subset(bter.res.all, padj < 0.1)

write.table(bter.res.sig, paste(wd, "Bter_DeSeqResults_FDRlt0.1.tsv", sep =""), 
            col.names = T, row.names = F, quote = F, sep = "\t")

#Genes of interest?

nona <- bter.res.all[!is.na(bter.res.all$padj),]
nona$LogDif <- sub('-', '', nona$log2FoldChange)
tmp <- nona[nona$Contrast == treats[1],]
sig.wound <- head(tmp$Gene[order(tmp$padj)], n = 10)
top.wound <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)
tmp <- nona[nona$Contrast == treats[2],]
sig.pos <- head(tmp$Gene[order(tmp$padj)], n = 10)
top.pos <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)
tmp <- nona[nona$Contrast == treats[3],]
sig.neg <-  head(tmp$Gene[order(tmp$padj)], n = 10)
top.neg <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)

goi <- c(sig.wound, top.wound, sig.pos, top.pos, sig.neg, top.neg)
stat <- c(rep("Top 10 Sig Wound", 10), rep("Top 10 log2FC Wound", 10),
          rep("Top 10 Sig Gram Positive", 10), rep("Top 10 log2FC Gram Positive", 10),
          rep("Top 10 Sig Gram Negative", 10), rep("Top 10 log2FC Gram Negative", 10))

goi.df <- as.data.frame(cbind(goi, stat))

goi.wound <- unique(c(sig.wound, top.wound))
goi.pos <- unique(c(sig.pos, top.pos))
goi.neg <- unique(c(sig.neg, top.neg))

###Visualise: Venn Diagrams

top.sig <- list("Wound" = sig.wound,
                 "Gram Positive" = sig.pos,
                 "Gram Negative" = sig.neg)

p <- ggvenn(top.sig,
       set_name_size = 5,
       show_percentage = F,
       fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Top 10 Genes per Treatment by adjPvalue")

dir.create(paste(wd, "Plots/", sep = ""))
ggsave(paste(wd, "Plots/Bter_TopAdjPvalue_Venn.pdf", sep = ""))

top.lfc <- list("Wound" = top.wound,
                 "Gram Positive" = top.pos,
                 "Gram Negative" = top.neg)

p <- ggvenn(top.lfc,
       set_name_size = 5,
       show_percentage = F,
       fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Top 10 Genes per Treatment by log2FC")

ggsave(paste(wd, "Plots/Bter_Toplog2FC_Venn.pdf", sep = ""))

goi <- list("Wound" = goi.wound,
                 "Gram Positive" = goi.pos,
                 "Gram Negative" = goi.neg)

p <- ggvenn(goi,
       set_name_size = 5,
       show_percentage = F,
       fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Genes of Interest Across Treatments")

ggsave(paste(wd, "Plots/Bter_GOI_Venn.pdf", sep = ""))

goi.mat <- makeMatrix(goi.df)
goi.mat

###Visualise: Expression BoxPlots
dir.create(paste(wd, "Plots/ExpressionPlots/", sep = ""))

####Significant Wound

goi.mat$Gene[goi.mat$SigWound == 1]
genes <- goi.mat$Gene[goi.mat$SigWound == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = bter.dds)
}

#All Wound: 

visiFacet(genes = genes,
          genenames = genes,
          dds = bter.dds, 
          filename = "WoundSig")

####TopFC Wound

goi.mat$Gene[goi.mat$TopWound == 1]
genes <- goi.mat$Gene[goi.mat$TopWound == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = bter.dds)
}

#All Wound: 

visiFacet(genes = genes,
          genenames = genes,
          dds = bter.dds, 
          filename = "Woundlog2FC")

####Significant GramPos

goi.mat$Gene[goi.mat$SigPos == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = bter.dds)
}

#All Gram Pos Sig: 

visiFacet(genes = genes,
          genenames = genes,
          dds = bter.dds, 
          filename = "PosSig")

####TopFC GramPos

goi.mat$Gene[goi.mat$TopPos == 1]
genes <- goi.mat$Gene[goi.mat$TopPos == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = bter.dds)
}

#All log Gram Pos: 

visiFacet(genes = genes,
          genenames = genes,
          dds = bter.dds, 
          filename = "Poslog2FC")

####Significant GramNeg

goi.mat$Gene[goi.mat$SigNeg == 1]
genes <- goi.mat$Gene[goi.mat$SigNeg == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = bter.dds)
}

visiFacet(genes = genes,
          genenames = genes,
          dds = bter.dds, 
          filename = "NegSig")

####TopFC GramPos

goi.mat$Gene[goi.mat$TopNeg == 1]
genes <- goi.mat$Gene[goi.mat$TopNeg == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = bter.dds)
}

#All log Gram Neg: 

visiFacet(genes = genes,
          genenames = genes,
          dds = bter.dds, 
          filename = "Neglog2FC")

###Visualise: PCA
#Redo after removing outliers.

vsd <- vst(bter.dds, blind = F)

bter.pca <- plotPCA(vsd, intgroup="Treatment", returnData = TRUE)

ggplot(bter.pca, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() 

ggsave(paste(wd, "Plots/Bter_PCA_Simple.pdf", sep = ""))

ggplot(bter.pca, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() +
  geom_label_repel(
    data = bter.pca,
    aes(label = name),
    size = 3
  )

ggsave(paste(wd, "Plots/Bter_PCA_Labelled.pdf", sep = ""))

bter.pca$Treat2 <- c(rep("Naive", 2), rep("Treatment", 6))

ggplot(bter.pca, aes(x = PC1, y = PC2, color = Treat2)) +
  geom_point() +
  stat_ellipse()

ggsave(paste(wd, "Plots/Bter_PCA_TrtVUnTrt.pdf", sep = ""))

###Visualise: Smears

bter.res.all$Gene[bter.res.all$Gene == "Kr-h1"] <- "Kr-h1_Bter"

ortho$GeneID[ortho$GeneID == "Kr-h1" & ortho$Species == "Bter"] <- "Kr-h1_Bter"

plot.df <- makeMetaDf(bter.res.all)
plot.df %>%
  filter(Condition == "Wound") %>%
  arrange(padj)

plot.df$Condition <- factor(plot.df$Condition, levels = c("Wound",
                                                          "Gram Positive",
                                                          "Gram Negative"))

ggplot(plot.df, aes(x = log(baseMean), y = log2FoldChange)) +
  geom_point(aes(colour = Sig)) +
  facet_grid(~Condition) +
  scale_color_manual(values=c("black", "Red"), guide = NULL) +
  theme(legend.position = "bottom")

ggsave(paste(wd, "Plots/Bter_SmearPlot.pdf", sep = ""))

###Visualise: Rainclouds

sig.genes <- unique(plot.df$Gene[plot.df$padj < 0.1])
plot.sig <- plot.df[plot.df$Gene %in% sig.genes,]

to.add.naive <- plot.sig[plot.sig$Gene %in% sig.genes &
                               plot.sig$Condition == "Wound",]
to.add.naive$Sig <- "No"
to.add.naive$Condition <- "Naive"
to.add.naive$log2FoldChange <- 0

plot.sig <- rbind(plot.sig, to.add.naive)

plot.sig$Condition <- factor(plot.sig$Condition, levels = c("Naive",
                                                            "Wound", 
                                                            "Gram Positive", 
                                                            "Gram Negative"))

plot.sig$Significance[plot.sig$padj > 0.1] <- "Not Significant"
plot.sig$Significance[plot.sig$padj < 0.1 & plot.sig$padj > 0.05] <- 
  "Significant (FDR < 0.1)"
plot.sig$Significance[plot.sig$padj < 0.05 & plot.sig$padj > 0.01] <-
  "Significant (FDR < 0.05)"
plot.sig$Significance[plot.sig$padj < 0.01 & plot.sig$padj > 0.001] <-
  "Significant (FDR < 0.01)"
plot.sig$Significance[plot.sig$padj < 0.001] <-
  "Significant (FDR < 0.001)"

plot.sig$Significance <- factor(plot.sig$Significance, levels =
                                 c("Not Significant",
                                   "Significant (FDR < 0.1)",
                                   "Significant (FDR < 0.05)",
                                   "Significant (FDR < 0.01)",
                                   "Significant (FDR < 0.001)"))


plot.sig$Function <- factor(plot.sig$Function, levels = c("Immune",
                                                          "Putative Immune", 
                                                          "None Immune", 
                                                          "Unknown"))

ggplot(plot.sig[!plot.sig$Condition == "Naive",], 
       (aes(x = Condition, y = log2FoldChange))) +
  geom_jitter(aes(colour = Function, shape = Significance)) +
  scale_colour_manual(values = c("#2c8cd9", "#7a2cd9", "darkgrey", "plum")) +
  scale_shape_manual(values = c(1, 5, 18, 19, 17)) 

ggsave(paste(wd, "Plots/Bter_logFCJitter.pdf", sep = ""))

ggplot(plot.sig, 
       (aes(x = Condition, y = log2FoldChange))) +
  geom_point(aes(colour = Function, shape = Significance)) +
  geom_line( linetype = "dashed", alpha = 0.2, aes(group = Gene, colour = Function)) +
  scale_colour_manual(values = c("#2c8cd9", "#7a2cd9", "darkgrey", "plum")) +
  scale_shape_manual(values = c(1, 5, 18, 19, 17)) 

ggsave(paste(wd, "Plots/Bter_logFCScatter_GeneLine.pdf", sep = ""))

###Visualise: Volcano Plots

plot.df$Function <- factor(plot.df$Function, levels = c("Immune",
                                                        "Putative Immune",
                                                        "None Immune",
                                                        "Unknown"))

ggplot(plot.df,
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
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Bter_Simple_FDR0.1_Volcano.pdf", sep = ""))

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.1", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Bter_Class_FDR0.1_Volcano.pdf", sep = ""))

plot.df$Direction2 <- "Not Significant"
plot.df$Direction2[plot.df$padj < 0.05 &
                    plot.df$log2FoldChange > 0 ] <- "Up-Regulated"
plot.df$Direction2[plot.df$padj < 0.05 &
                    plot.df$log2FoldChange < 0 ] <- "Down-Regulated"

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction2), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.05", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Bter_Simple_FDR0.05_Volcano.pdf", sep = ""))

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction2, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.05", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Bter_Class_FDR0.05_Volcano.pdf", sep = ""))

plot.df$Direction3 <- "Not Significant"
plot.df$Direction3[plot.df$padj < 0.001 &
                    plot.df$log2FoldChange > 0 ] <- "Up-Regulated"
plot.df$Direction3[plot.df$padj < 0.001 &
                    plot.df$log2FoldChange < 0 ] <- "Down-Regulated"

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Bter_Simple_FDR0.001_Volcano.pdf", sep = ""))

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Bter_Class_FDR0.001_Volcano.pdf", sep = ""))

#Label top 3 sig genes per treatment

sig.wound[1:3]
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC100631061"] <- "Hymenoptaecin"
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC125385347"] <- "Hymenoptaecin-like"
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC100649281"] <- "Maltase A1-like"

sig.pos[1:3]
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC100631061"] <- "Hymenoptaecin"
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC100644101"] <- "Uncharacterised LOC100644101"
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC100649867"] <- "Apidaecins type 73"

sig.neg[1:3]
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC100631061"] <- "Hymenoptaecin"
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC125385347"] <- "Hymenoptaecin-like"
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC100649867"] <- "Apidaecins type 73"

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") +
  geom_label_repel(
    data = plot.df[!is.na(plot.df$Label),],
    size = 2,
    aes(label = Label),
    point.padding = unit(0.3, "lines"),
    box.padding = 1
  ) 

ggsave(paste(wd, "Plots/Bter_Simple_FDR0.001_Volcano_top3Label.pdf", sep = ""))

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") +
  geom_label_repel(
    data = plot.df[!is.na(plot.df$Label),],
    size = 2,
    aes(label = Label),
    point.padding = unit(0.3, "lines"),
    box.padding = 1
  ) 

ggsave(paste(wd, "Plots/Bter_Class_FDR0.001_Volcano_top3Label.pdf", sep = ""))


##Ceratina australensis####

dir.create("output/Ind_DE/Caus/")
wd <- "output/Ind_DE/Caus/"

###DE Analysis

caus.tx <- geneCount("Caus")
head(caus.tx$counts)

#Build a DESeq DataSet from the above, prefilter (keep only rows that have at least 10
#reads total), and return dds object.

caus.dds <- getDDS("Caus", caus.tx)
caus.res <- vector(mode = "list", length = 3)
treats <- resultsNames(caus.dds)[c(4,3,2)]

for(i in 1:length(caus.res)){
  caus.res[[i]] <- results(caus.dds, name = paste(treats[i]))
  summary(caus.res[[i]])
}

###Check for outliers
#First, check counts

boxplot(log10(assays(caus.dds)[["cooks"]]), range=0, las=2)

#normalise
vsd <- vst(caus.dds, blind = F)
#make heatmap (top 20 genes)
select <- order(rowMeans(counts(caus.dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
pheatmap(assay(vsd)[select,])

#This is messy I'm not sure where to begin.

pdf("output/Ind_DE/Caus/Caus_Pheatmap.pdf")
  pheatmap(assay(vsd)[select,])
dev.off()

#Check with sample to sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

#Yeah that's pretty messy

pdf("output/Ind_DE/Caus/Caus_DistanceMatrix.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()

#Finally, PCA

caus.pca <- plotPCA(vsd, intgroup="Treatment", returnData = TRUE)
ggplot(caus.pca, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
  geom_point() +
  geom_text()

ggsave("output/Ind_DE/Caus/Caus_PCA.pdf")

#P55 and SL56 seem a bit ... away
#What happens if I remove the potential outliers ? 
  
torem <- c("CA_P_055", "CA_SL_056", "CA_N_012")
outCheck(torem, "Caus")

#2 plots are still a mess, but the sample-distance plots are very much improved.
#Now N 12 is a problem ? Go back and try again....

#Much better.... though I've lost 3 samples!

caus.dds2 <- ddsOmit(torem, "Caus")

caus.res2 <- vector(mode = "list", length = 3)
treats2 <- resultsNames(caus.dds2)[c(4,3,2)]

for(i in 1:length(caus.res2)){
  caus.res2[[i]] <- results(caus.dds2, name = paste(treats2[i]))
  print(treats2[i])
  summary(caus.res2[[i]])
}

#Save plots

vsd2 <- vst(caus.dds2, blind = F)
#make heatmap (top 20 genes)
select <- order(rowMeans(counts(caus.dds2,normalized=TRUE)),
                decreasing=TRUE)[1:20]
pheatmap(assay(vsd2)[select,])

pdf("output/Ind_DE/Caus/Caus_OutlierRemoved_Pheatmap.pdf")
  pheatmap(assay(vsd2)[select,])
dev.off()

sampleDists2 <- dist(t(assay(vsd2)))
sampleDistMatrix2 <- as.matrix(sampleDists2)
rownames(sampleDistMatrix2) <- paste(vsd2$Treatment, vsd2$type, sep="-")
colnames(sampleDistMatrix2) <- NULL

pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)

pdf("output/Ind_DE/Caus/Caus_OutlierRemoved_DistanceMatrix.pdf")
  pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)
dev.off()

caus.pca2 <- plotPCA(vsd2, intgroup="Treatment", returnData = TRUE)
ggplot(caus.pca2, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
  geom_point() +
  geom_text()
ggsave("output/Ind_DE/Caus/Caus_OutlierRemoved_PCA.pdf")


#There is quite a significant difference .... gonna make the call to go forward 
#without the outliers 

caus.dds <- caus.dds2
caus.res <- caus.res2

###Results

Save 

for (i in 1:length(caus.res)){
  caus.res[[i]]$Contrast <- paste(treats[i])
  caus.res[[i]]$Gene <- rownames(caus.res[[i]])
  rownames(caus.res[[i]]) <- NULL
  caus.res[[i]] <- as.data.frame(caus.res[[i]])
}

caus.res.all <- bind_rows(caus.res[1:3])

id <- read.table("input/Gene_Info/Caus_IDwCcal.tsv", 
                 sep = "\t", quote = "", header = T)

for (i in 1:nrow(caus.res.all)){
  d <- unique(id$CcalID[id$Gene == caus.res.all$Gene[i]])
  if(length(d) > 1){
    caus.res.all$Desc[i] <- d[1]
  } else {
    caus.res.all$Desc[i] <- d  
  }
}

write.table(caus.res.all, paste(wd, "Caus_DeSeqResults_All.tsv", sep = ""),
            col.names = T, row.names = F, quote = F, sep = "\t")

caus.res.sig <- subset(caus.res.all, padj < 0.1)

write.table(caus.res.sig, paste(wd, "Caus_DeSeqResults_FDRlt0.1.tsv", sep =""), 
            col.names = T, row.names = F, quote = F, sep = "\t")

#Genes of interest?
#I'm going to take the top 10 most significant and top 10 biggest logFC per treatment. 
#Genes that overlap between two (and between the three treatments) will be of particular
#interest

nona <- caus.res.all[!is.na(caus.res.all$padj),]
nona$LogDif <- sub('-', '', nona$log2FoldChange)
tmp <- nona[nona$Contrast == treats[1],]
sig.wound <- head(tmp$Gene[order(tmp$padj)], n = 10)
top.wound <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)
tmp <- nona[nona$Contrast == treats[2],]
sig.pos <- head(tmp$Gene[order(tmp$padj)], n = 10)
top.pos <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)
tmp <- nona[nona$Contrast == treats[3],]
sig.neg <-  head(tmp$Gene[order(tmp$padj)], n = 10)
top.neg <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)

goi <- c(sig.wound, top.wound, sig.pos, top.pos, sig.neg, top.neg)
stat <- c(rep("Top 10 Sig Wound", 10), rep("Top 10 log2FC Wound", 10),
          rep("Top 10 Sig Gram Positive", 10), rep("Top 10 log2FC Gram Positive", 10),
          rep("Top 10 Sig Gram Negative", 10), rep("Top 10 log2FC Gram Negative", 10))

goi.df <- as.data.frame(cbind(goi, stat))


goi.wound <- unique(c(sig.wound, top.wound))
goi.pos <- unique(c(sig.pos, top.pos))
goi.neg <- unique(c(sig.neg, top.neg))

###Visualise: Venn Diagrams

top.sig <- list("Wound" = sig.wound,
                "Gram Positive" = sig.pos,
                "Gram Negative" = sig.neg)

p <- ggvenn(top.sig,
            set_name_size = 5,
            show_percentage = F,
            fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Top 10 Genes per Treatment by adjPvalue")

dir.create(paste(wd, "Plots/", sep = ""))
ggsave(paste(wd, "Plots/Caus_TopAdjPvalue_Venn.pdf", sep = ""))

top.lfc <- list("Wound" = top.wound,
                "Gram Positive" = top.pos,
                "Gram Negative" = top.neg)

p <- ggvenn(top.lfc,
            set_name_size = 5,
            show_percentage = F,
            fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Top 10 Genes per Treatment by log2FC")

ggsave(paste(wd, "Plots/Caus_Toplog2FC_Venn.pdf", sep = ""))

goi <- list("Wound" = goi.wound,
            "Gram Positive" = goi.pos,
            "Gram Negative" = goi.neg)

p <- ggvenn(goi,
            set_name_size = 5,
            show_percentage = F,
            fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Genes of Interest Across Treatments")

ggsave(paste(wd, "Plots/Caus_GOI_Venn.pdf", sep = ""))

goi.mat <- makeMatrix(goi.df)
goi.mat

####Visualise: Expression BoxPlots

dir.create(paste(wd, "Plots/ExpressionPlots/", sep = ""))

####Significant Wound

goi.mat$Gene[goi.mat$SigWound == 1]
genes <- goi.mat$Gene[goi.mat$SigWound == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = caus.dds)
}

#All Wound: 

visiFacet(genes = genes,
          genenames = genes,
          dds = caus.dds, 
          filename = "WoundSig")

####TopFC Wound

goi.mat$Gene[goi.mat$TopWound == 1]
genes <- goi.mat$Gene[goi.mat$TopWound == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = caus.dds)
}

#All Wound: 

visiFacet(genes = genes,
          genenames = genes,
          dds = caus.dds, 
          filename = "Woundlog2FC")

####Significant GramPos

goi.mat$Gene[goi.mat$SigPos == 1]
genes <- goi.mat$Gene[goi.mat$SigPos == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = caus.dds)
}

#All Gram Pos Sig: 

visiFacet(genes = genes,
          genenames = genes,
          dds = caus.dds, 
          filename = "PosSig")

####TopFC GramPos

goi.mat$Gene[goi.mat$TopPos == 1]
genes <- goi.mat$Gene[goi.mat$TopPos == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = caus.dds)
}

#All log Gram Pos: 

visiFacet(genes = genes,
          genenames = genes,
          dds = caus.dds, 
          filename = "Poslog2FC")

####Significant GramNeg

goi.mat$Gene[goi.mat$SigNeg == 1]
genes <- goi.mat$Gene[goi.mat$SigNeg == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = caus.dds)
}

#All Gram Neg Sig: 

visiFacet(genes = genes,
          genenames = genes,
          dds = caus.dds, 
          filename = "NegSig")

####TopFC GramNeg

goi.mat$Gene[goi.mat$TopNeg == 1]
genes <- goi.mat$Gene[goi.mat$TopNeg == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = caus.dds)
}

#All log Gram Neg: 

visiFacet(genes = genes,
          genenames = genes,
          dds = caus.dds, 
          filename = "Neglog2FC")

###Visualise: PCA
#Redo after removing outliers.

vsd <- vst(caus.dds, blind = F)

caus.pca <- plotPCA(vsd, intgroup="Treatment", returnData = TRUE)

ggplot(caus.pca, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() 

ggsave(paste(wd, "Plots/Caus_PCA_Simple.pdf", sep = ""))

ggplot(caus.pca, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() +
  geom_label_repel(
    data = caus.pca,
    aes(label = name),
    size = 3
  )

ggsave(paste(wd, "Plots/Caus_PCA_Labelled.pdf", sep = ""))

caus.pca$Treat2 <- c(rep("Naive", 2), rep("Treatment", 7))

ggplot(caus.pca, aes(x = PC1, y = PC2, color = Treat2)) +
  geom_point() +
  stat_ellipse()

ggsave(paste(wd, "Plots/Caus_PCA_TrtVUnTrt.pdf", sep = ""))

###Visualise: Smears

plot.df <- makeMetaDf(caus.res.all)
plot.df %>%
  filter(Condition == "Wound") %>%
  arrange(padj)

plot.df$Condition <- factor(plot.df$Condition, levels = c("Wound",
                                                          "Gram Positive",
                                                          "Gram Negative"))

ggplot(plot.df, aes(x = log(baseMean), y = log2FoldChange)) +
  geom_point(aes(colour = Sig)) +
  facet_grid(~Condition) +
  scale_color_manual(values=c("black", "Red"), guide = NULL) +
  theme(legend.position = "bottom")

ggsave(paste(wd, "Plots/Caus_SmearPlot.pdf", sep = ""))

###Visualise: Rainclouds

sig.genes <- unique(plot.df$Gene[plot.df$padj < 0.1])
plot.sig <- plot.df[plot.df$Gene %in% sig.genes,]

to.add.naive <- plot.sig[plot.sig$Gene %in% sig.genes &
                           plot.sig$Condition == "Wound",]
to.add.naive$Sig <- "No"
to.add.naive$Condition <- "Naive"
to.add.naive$log2FoldChange <- 0

plot.sig <- rbind(plot.sig, to.add.naive)

plot.sig$Condition <- factor(plot.sig$Condition, levels = c("Naive",
                                                            "Wound", 
                                                            "Gram Positive", 
                                                            "Gram Negative"))

plot.sig$Significance[plot.sig$padj > 0.1] <- "Not Significant"
plot.sig$Significance[plot.sig$padj < 0.1 & plot.sig$padj > 0.05] <- 
  "Significant (FDR < 0.1)"
plot.sig$Significance[plot.sig$padj < 0.05 & plot.sig$padj > 0.01] <-
  "Significant (FDR < 0.05)"
plot.sig$Significance[plot.sig$padj < 0.01 & plot.sig$padj > 0.001] <-
  "Significant (FDR < 0.01)"
plot.sig$Significance[plot.sig$padj < 0.001] <-
  "Significant (FDR < 0.001)"

plot.sig$Significance <- factor(plot.sig$Significance, levels =
                                  c("Not Significant",
                                    "Significant (FDR < 0.1)",
                                    "Significant (FDR < 0.05)",
                                    "Significant (FDR < 0.01)",
                                    "Significant (FDR < 0.001)"))

plot.sig$Function <- factor(plot.sig$Function, levels = c("Immune",
                                                          "Putative Immune", 
                                                          "None Immune", 
                                                          "Unknown"))

ggplot(plot.sig[!plot.sig$Condition == "Naive",], 
       (aes(x = Condition, y = log2FoldChange))) +
  geom_jitter(aes(colour = Function, shape = Significance)) +
  scale_colour_manual(values = c("#2c8cd9", "#7a2cd9", "darkgrey", "plum")) +
  scale_shape_manual(values = c(1, 5, 18, 19, 17)) 

ggsave(paste(wd, "Plots/Caus_logFCJitter.pdf", sep = ""))

ggplot(plot.sig, 
       (aes(x = Condition, y = log2FoldChange))) +
  geom_point(aes(colour = Function, shape = Significance)) +
  geom_line( linetype = "dashed", alpha = 0.2, aes(group = Gene, colour = Function)) +
  scale_colour_manual(values = c("#2c8cd9", "#7a2cd9", "darkgrey", "plum")) +
  scale_shape_manual(values = c(1, 5, 18, 19, 17)) 

ggsave(paste(wd, "Plots/Caus_logFCScatter_GeneLine.pdf", sep = ""))

####Visualise: Volcano Plots

plot.df$Function <- factor(plot.df$Function, levels = c("Immune",
                                                        "Putative Immune",
                                                        "None Immune",
                                                        "Unknown"))
ggplot(plot.df,
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
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Caus_Simple_FDR0.1_Volcano.pdf", sep = ""))

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.1", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Caus_Class_FDR0.1_Volcano.pdf", sep = ""))

plot.df$Direction2 <- "Not Significant"
plot.df$Direction2[plot.df$padj < 0.05 &
                     plot.df$log2FoldChange > 0 ] <- "Up-Regulated"
plot.df$Direction2[plot.df$padj < 0.05 &
                     plot.df$log2FoldChange < 0 ] <- "Down-Regulated"

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction2), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.05", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Caus_Simple_FDR0.05_Volcano.pdf", sep = ""))

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction2, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.05", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Caus_Class_FDR0.05_Volcano.pdf", sep = ""))

plot.df$Direction3 <- "Not Significant"
plot.df$Direction3[plot.df$padj < 0.001 &
                     plot.df$log2FoldChange > 0 ] <- "Up-Regulated"
plot.df$Direction3[plot.df$padj < 0.001 &
                     plot.df$log2FoldChange < 0 ] <- "Down-Regulated"

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Caus_Simple_FDR0.001_Volcano.pdf", sep = ""))

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Caus_Class_FDR0.001_Volcano.pdf", sep = ""))

#Label top 3 sig genes per treatment

sig.wound[1:3]
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "Caust.v2_016347"] <- "ATP-citrate synthase	"
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "Caust.v2_020967"] <- "Piwi-like protein Siwi"
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "Caust.v2_011640"] <- "TIM50-like"

sig.pos[1:3]
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "Caust.v2_018203"] <- "NT5D3-like"
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "Caust.v2_011034"] <- "Uncharacterised LOC108629120"
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "Caust.v2_011210"] <- "Med24"

sig.neg[1:3]
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "Caust.v2_011640"] <- "TIM50-like"
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "Caust.v2_011034"] <- "Uncharacterised LOC108629120"
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "Caust.v2_011550"] <- "Uncharacterised LOC108623077"

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") +
  geom_label_repel(
    data = plot.df[!is.na(plot.df$Label),],
    size = 2,
    aes(label = Label),
    point.padding = unit(0.3, "lines"),
    box.padding = 1
  ) 

ggsave(paste(wd, "Plots/Caus_Simple_FDR0.001_Volcano_top3Label.pdf", sep = ""))

ggplot(plot.df,
       aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
  facet_grid(~Condition, scales = "free_y") +
  geom_label_repel(
    data = plot.df[!is.na(plot.df$Label),],
    size = 2,
    aes(label = Label),
    point.padding = unit(0.3, "lines"),
    box.padding = 2.25
  ) 

ggsave(paste(wd, "Plots/Caus_Class_FDR0.001_Volcano_top3Label.pdf", sep = ""))


##Polistes lanio####

dir.create("output/Ind_DE/Plan/")
wd <- "output/Ind_DE/Plan/"

###DE Analysis
#First, convert transcript-level counts table from Kallisto to gene level counts and 
#produce a coldata object

plan.tx <- geneCount("Plan")
head(plan.tx$counts)

#Build a DESeq DataSet from the above, prefilter (keep only rows that have at least 
#10 reads total), and return dds object.

plan.dds <- getDDS("plan", plan.tx)
plan.res <- vector(mode = "list", length = 3)
treats <- resultsNames(plan.dds)[c(4,3,2)]

for(i in 1:length(plan.res)){
  plan.res[[i]] <- results(plan.dds, name = paste(treats[i]))
  summary(plan.res[[i]])
}

###Check for outliers

boxplot(log10(assays(plan.dds)[["cooks"]]), range=0, las=2)

#normalise
vsd <- vst(plan.dds, blind = F)
#make heatmap (top 20 genes)
select <- order(rowMeans(counts(plan.dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
pheatmap(assay(vsd)[select,])

#N1908 may be an issue.

pdf("output/Ind_DE/Plan/Plan_Pheatmap.pdf")
  pheatmap(assay(vsd)[select,])
dev.off()

#Check with sample to sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

#Yeah that's pretty messy

pdf("output/Ind_DE/Plan/Plan_DistanceMatrix.pdf")
  pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()

#Finally, PCA

plan.pca <- plotPCA(vsd, intgroup="Treatment", returnData = TRUE)
ggplot(plan.pca, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
  geom_point() +
  geom_text()

#And N1909 ...
#I can't remove two from the same treatment, I need an absolute minimum of 2.

ggsave("output/Ind_DE/Plan/Plan_PCA.pdf")

#Going with 1908 and SM1936
#If I could remove 1909 as well, I would
#What happens if I remove the potential outliers ? 

torem <- c("PL_N_1908", "PL_SM_1936", "PL_P_1912")
outCheck(torem, "Plan")

#2 plots are still a mess, but the sample-distance plots are very much improved.
#Hmm. P1912....

plan.dds2 <- ddsOmit(torem, "Plan")

plan.res2 <- vector(mode = "list", length = 3)
treats2 <- resultsNames(plan.dds2)[c(4,3,2)]

for(i in 1:length(plan.res2)){
  plan.res2[[i]] <- results(plan.dds2, name = paste(treats2[i]))
  print(treats2[i])
  summary(plan.res2[[i]])
}

#Wow. That's a lot of genes.

#Save plots

vsd2 <- vst(plan.dds2, blind = F)
#make heatmap (top 20 genes)
select <- order(rowMeans(counts(plan.dds2,normalized=TRUE)),
                decreasing=TRUE)[1:20]
pheatmap(assay(vsd2)[select,])

pdf("output/Ind_DE/Plan/Plan_OutlierRemoved_Pheatmap.pdf")
  pheatmap(assay(vsd2)[select,])
dev.off()

sampleDists2 <- dist(t(assay(vsd2)))
sampleDistMatrix2 <- as.matrix(sampleDists2)
rownames(sampleDistMatrix2) <- paste(vsd2$Treatment, vsd2$type, sep="-")
colnames(sampleDistMatrix2) <- NULL

pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)

pdf("output/Ind_DE/Plan/Plan_OutlierRemoved_DistanceMatrix.pdf")
  pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)
dev.off()

plan.pca2 <- plotPCA(vsd2, intgroup="Treatment", returnData = TRUE)
ggplot(plan.pca2, aes(x = PC1, y = PC2, color = Treatment, label = name)) +
  geom_point() +
  geom_text()
ggsave("output/Ind_DE/Plan/Plan_OutlierRemoved_PCA.pdf")

#Gonna make the call to go forward without the outliers 

plan.dds <- plan.dds2
plan.res <- plan.res2

###Results
#Save 

for (i in 1:length(plan.res)){
  plan.res[[i]]$Contrast <- paste(treats[i])
  plan.res[[i]]$Gene <- rownames(plan.res[[i]])
  rownames(plan.res[[i]]) <- NULL
  plan.res[[i]] <- as.data.frame(plan.res[[i]])
}

plan.res.all <- bind_rows(plan.res[1:3])

id <- read.table("input/Gene_Info/Plan_IDwDesc.tsv", 
                 sep = "\t", quote = "", header = T)

for (i in 1:nrow(plan.res.all)){
  d <- unique(id$Description[id$GeneID == plan.res.all$Gene[i]])
  if(length(d) > 1){
    plan.res.all$Desc[i] <- d[1]
  } else {
    plan.res.all$Desc[i] <- d  
  }
}

write.table(plan.res.all, paste(wd, "Plan_DeSeqResults_All.tsv", sep = ""),
            col.names = T, row.names = F, quote = F, sep = "\t")

plan.res.sig <- subset(plan.res.all, padj < 0.1)

write.table(plan.res.sig, paste(wd, "Plan_DeSeqResults_FDRlt0.1.tsv", sep =""), 
            col.names = T, row.names = F, quote = F, sep = "\t")

#Genes of interest?
  
#I'm going to take the top 10 most significant and top 10 biggest logFC per treatment. 
#Genes that overlap between two (and between the three treatments) will be of particular
#interest

nona <- plan.res.all[!is.na(plan.res.all$padj),]
nona$LogDif <- sub('-', '', nona$log2FoldChange)
tmp <- nona[nona$Contrast == treats[1],]
sig.wound <- head(tmp$Gene[order(tmp$padj)], n = 10)
top.wound <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)
tmp <- nona[nona$Contrast == treats[2],]
sig.pos <- head(tmp$Gene[order(tmp$padj)], n = 10)
top.pos <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)
tmp <- nona[nona$Contrast == treats[3],]
sig.neg <-  head(tmp$Gene[order(tmp$padj)], n = 10)
top.neg <- head(tmp$Gene[order(tmp$LogDif, decreasing = T)], n = 10)

goi <- c(sig.wound, top.wound, sig.pos, top.pos, sig.neg, top.neg)
stat <- c(rep("Top 10 Sig Wound", 10), rep("Top 10 log2FC Wound", 10),
          rep("Top 10 Sig Gram Positive", 10), rep("Top 10 log2FC Gram Positive", 10),
          rep("Top 10 Sig Gram Negative", 10), rep("Top 10 log2FC Gram Negative", 10))

goi.df <- as.data.frame(cbind(goi, stat))


goi.wound <- unique(c(sig.wound, top.wound))
goi.pos <- unique(c(sig.pos, top.pos))
goi.neg <- unique(c(sig.neg, top.neg))

###Visualise: Venn Diagrams

top.sig <- list("Wound" = sig.wound,
                 "Gram Positive" = sig.pos,
                 "Gram Negative" = sig.neg)

p <- ggvenn(top.sig,
       set_name_size = 5,
       show_percentage = F,
       fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Top 10 Genes per Treatment by adjPvalue")

dir.create(paste(wd, "Plots/", sep = ""))
ggsave(paste(wd, "Plots/Plan_TopAdjPvalue_Venn.pdf", sep = ""))

top.lfc <- list("Wound" = top.wound,
                 "Gram Positive" = top.pos,
                 "Gram Negative" = top.neg)

p <- ggvenn(top.lfc,
       set_name_size = 5,
       show_percentage = F,
       fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Top 10 Genes per Treatment by log2FC")

ggsave(paste(wd, "Plots/Plan_Toplog2FC_Venn.pdf", sep = ""))

goi <- list("Wound" = goi.wound,
                 "Gram Positive" = goi.pos,
                 "Gram Negative" = goi.neg)

p <- ggvenn(goi,
       set_name_size = 5,
       show_percentage = F,
       fill_color = c("yellow", "blue", "green"))
p +  ggtitle("Genes of Interest Across Treatments")

ggsave(paste(wd, "Plots/Plan_GOI_Venn.pdf", sep = ""))

goi.mat <- makeMatrix(goi.df)
goi.mat


###Visualise: Expression BoxPlots

dir.create(paste(wd, "Plots/ExpressionPlots/", sep = ""))

####Significant Wound

goi.mat$Gene[goi.mat$SigWound == 1]
genes <- goi.mat$Gene[goi.mat$SigWound == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = plan.dds)
}

#All Wound: 

visiFacet(genes = genes,
          genenames = genes,
          dds = plan.dds, 
          filename = "WoundSig")

####TopFC Wound

goi.mat$Gene[goi.mat$TopWound == 1]
genes <- goi.mat$Gene[goi.mat$TopWound == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = plan.dds)
}

#All Wound: 

visiFacet(genes = genes,
          genenames = genes,
          dds = plan.dds, 
          filename = "Woundlog2FC")

####Significant GramPos

goi.mat$Gene[goi.mat$SigPos == 1]
genes <- goi.mat$Gene[goi.mat$SigPos == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = plan.dds)
}

#All Gram Pos Sig: 

visiFacet(genes = genes,
          genenames = genes,
          dds = plan.dds, 
          filename = "PosSig")

####TopFC GramPos

goi.mat$Gene[goi.mat$TopPos == 1]
genes <- goi.mat$Gene[goi.mat$TopPos == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = plan.dds)
}

#All log Gram Pos: 

visiFacet(genes = genes,
          genenames = genes,
          dds = plan.dds, 
          filename = "Poslog2FC")

####Significant GramNeg

goi.mat$Gene[goi.mat$SigNeg == 1]
genes <- goi.mat$Gene[goi.mat$SigNeg == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = plan.dds)
}

#All Gram Neg Sig: 

visiFacet(genes = genes,
          genenames = genes,
          dds = plan.dds, 
          filename = "NegSig")

####TopFC GramPos

goi.mat$Gene[goi.mat$TopNeg == 1]
genes <- goi.mat$Gene[goi.mat$TopNeg == 1]

for (i in 1:length(genes)){
  visiGene(gene = genes[i], genename = genes[i], dds = plan.dds)
}

#All log Gram Neg: 

visiFacet(genes = genes,
          genenames = genes,
          dds = plan.dds, 
          filename = "Neglog2FC")

###Visualise: PCA
#Redo after removing outliers.

vsd <- vst(plan.dds, blind = F)

plan.pca <- plotPCA(vsd, intgroup="Treatment", returnData = TRUE)

ggplot(plan.pca, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() 

ggsave(paste(wd, "Plots/Plan_PCA_Simple.pdf", sep = ""))

ggplot(plan.pca, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point() +
  geom_label_repel(
    data = plan.pca,
    aes(label = name),
    size = 3
  )

ggsave(paste(wd, "Plots/Plan_PCA_Labelled.pdf", sep = ""))

plan.pca$Treat2 <- c(rep("Naive", 2), rep("Treatment", 7))

ggplot(plan.pca, aes(x = PC1, y = PC2, color = Treat2)) +
  geom_point() +
  stat_ellipse()

ggsave(paste(wd, "Plots/Plan_PCA_TrtVUnTrt.pdf", sep = ""))

###Visualise: Smears

plot.df <- makeMetaDf(plan.res.all)
plot.df %>%
  filter(Condition == "Wound") %>%
  arrange(padj)

plot.df$Condition <- factor(plot.df$Condition, levels = c("Wound",
                                                          "Gram Positive",
                                                          "Gram Negative"))

ggplot(plot.df, aes(x = log(baseMean), y = log2FoldChange)) +
  geom_point(aes(colour = Sig)) +
  facet_grid(~Condition) +
  scale_color_manual(values=c("black", "Red"), guide = NULL) +
  theme(legend.position = "bottom")

ggsave(paste(wd, "Plots/Plan_SmearPlot.pdf", sep = ""))

###Visualise: Rainclouds

sig.genes <- unique(plot.df$Gene[plot.df$padj < 0.1])
plot.sig <- plot.df[plot.df$Gene %in% sig.genes,]

to.add.naive <- plot.sig[plot.sig$Gene %in% sig.genes &
                               plot.sig$Condition == "Wound",]
to.add.naive$Sig <- "No"
to.add.naive$Condition <- "Naive"
to.add.naive$log2FoldChange <- 0

plot.sig <- rbind(plot.sig, to.add.naive)

plot.sig$Condition <- factor(plot.sig$Condition, levels = c("Naive",
                                                            "Wound", 
                                                            "Gram Positive", 
                                                            "Gram Negative"))

plot.sig$Significance[plot.sig$padj > 0.1] <- "Not Significant"
plot.sig$Significance[plot.sig$padj < 0.1 & plot.sig$padj > 0.05] <- 
  "Significant (FDR < 0.1)"
plot.sig$Significance[plot.sig$padj < 0.05 & plot.sig$padj > 0.01] <-
  "Significant (FDR < 0.05)"
plot.sig$Significance[plot.sig$padj < 0.01 & plot.sig$padj > 0.001] <-
  "Significant (FDR < 0.01)"
plot.sig$Significance[plot.sig$padj < 0.001] <-
  "Significant (FDR < 0.001)"

plot.sig$Significance <- factor(plot.sig$Significance, levels =
                                 c("Not Significant",
                                   "Significant (FDR < 0.1)",
                                   "Significant (FDR < 0.05)",
                                   "Significant (FDR < 0.01)",
                                   "Significant (FDR < 0.001)"))


plot.sig$Function <- factor(plot.sig$Function, levels = c("Immune",
                                                          "Putative Immune", 
                                                          "None Immune", 
                                                          "Unknown"))

ggplot(plot.sig[!plot.sig$Condition == "Naive",], 
       (aes(x = Condition, y = log2FoldChange))) +
  geom_jitter(aes(colour = Function, shape = Significance)) +
  scale_colour_manual(values = c("#2c8cd9", "#7a2cd9", "darkgrey", "plum")) +
  scale_shape_manual(values = c(1, 5, 18, 19, 17)) 

ggsave(paste(wd, "Plots/Plan_logFCJitter.pdf", sep = ""))

ggplot(plot.sig, 
       (aes(x = Condition, y = log2FoldChange))) +
  geom_point(aes(colour = Function, shape = Significance)) +
  geom_line( linetype = "dashed", alpha = 0.2, aes(group = Gene, colour = Function)) +
  scale_colour_manual(values = c("#2c8cd9", "#7a2cd9", "darkgrey", "plum")) +
  scale_shape_manual(values = c(1, 5, 18, 19, 17)) 

ggsave(paste(wd, "Plots/Plan_logFCScatter_GeneLine.pdf", sep = ""))

####Visualise: Volcano Plots

plot.df$Function <- factor(plot.df$Function, levels = c("Immune",
                                                        "Putative Immune",
                                                        "None Immune",
                                                        "Unknown"))

ggplot(plot.df,
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
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Plan_Simple_FDR0.1_Volcano.pdf", sep = ""))

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.1", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Plan_Class_FDR0.1_Volcano.pdf", sep = ""))

plot.df$Direction2 <- "Not Significant"
plot.df$Direction2[plot.df$padj < 0.05 &
                    plot.df$log2FoldChange > 0 ] <- "Up-Regulated"
plot.df$Direction2[plot.df$padj < 0.05 &
                    plot.df$log2FoldChange < 0 ] <- "Down-Regulated"

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction2), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.05", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Plan_Simple_FDR0.05_Volcano.pdf", sep = ""))

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction2, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.05", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Plan_Class_FDR0.05_Volcano.pdf", sep = ""))

plot.df$Direction3 <- "Not Significant"
plot.df$Direction3[plot.df$padj < 0.001 &
                    plot.df$log2FoldChange > 0 ] <- "Up-Regulated"
plot.df$Direction3[plot.df$padj < 0.001 &
                    plot.df$log2FoldChange < 0 ] <- "Down-Regulated"

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Plan_Simple_FDR0.001_Volcano.pdf", sep = ""))

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") 

ggsave(paste(wd, "Plots/Plan_Class_FDR0.001_Volcano.pdf", sep = ""))

#Label top 3 sig genes per treatment

sig.wound[1:3]
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC106790552"] <- "Uncharacterised LOC106790552"
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC106787312"] <- "Tyrosine decarboxylase-like"
plot.df$Label[plot.df$Condition == "Wound" & 
                plot.df$Gene == "LOC106791417"] <- "Protein croquemort-like"

sig.pos[1:3]
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC106787312"] <- "Tyrosine decarboxylase-like"
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC106793670"] <- "UDP-GUT 1-3-like"
plot.df$Label[plot.df$Condition == "Gram Positive" & 
                plot.df$Gene == "LOC106790552"] <- "Uncharacterised LOC106790552"

sig.neg[1:3]
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC106785761"] <- "L(2)EFL-like"
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC106786250"] <- "Protein bicaudal C"
plot.df$Label[plot.df$Condition == "Gram Negative" & 
                plot.df$Gene == "LOC106786125"] <- "Inhibin beta chain"

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3), size = 1, 
             position = position_jitter(w=0.1, h = 0.1)) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  guides(colour = guide_legend("FDR < 0.001", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 3)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") +
  geom_label_repel(
    data = plot.df[!is.na(plot.df$Label),],
    size = 2,
    aes(label = Label),
    point.padding = unit(0.3, "lines"),
    box.padding = 1.5
  ) 

ggsave(paste(wd, "Plots/Plan_Simple_FDR0.001_Volcano_top3Label.pdf", sep = ""))

ggplot(plot.df,
     aes(x = log2FoldChange, y = Neglog10adjPvalue)) +
  geom_point(aes(colour = Direction3, shape = Function), size = 1) +
  scale_colour_manual(values=c("Up-Regulated" = "red",
                               "Not Top" = "black",  
                               "Down-Regulated" = "blue" ),
                      breaks = c("Up-Regulated", "Down-Regulated")) +
  scale_shape_manual(values= c(19, 17, 1, 3)) +
  guides(colour = guide_legend("FDR < 0.01", nrow = 2),
         shape = guide_legend("Gene Function", nrow = 4)) +
  ylab("-log10(adjPvalue)") +
  theme(legend.position = "bottom") +
    facet_grid(~Condition, scales = "free_y") +
  geom_label_repel(
    data = plot.df[!is.na(plot.df$Label),],
    size = 2,
    aes(label = Label),
    point.padding = unit(0.3, "lines"),
    box.padding = 1.5
  ) 

ggsave(paste(wd, "Plots/Plan_Class_FDR0.001_Volcano_top3Label.pdf", sep = ""))

##SessionInfo####
sessionInfo()


