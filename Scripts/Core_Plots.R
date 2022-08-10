#28th July 2022
#Plots for Core Project not included in RMD 

library(dplyr)
library(UpSetR)
library(ggvenn)
library("ggVennDiagram")
library("ggplot2")



library("ggVennDiagram")

####Upset Plot of Sig Genes across 4 individual species#####
#Requires: result tsv, ortho information, annotation

#result tsvs
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

summary(as.factor(ind.res$Source))

tail(ind.res)

#add orthogroup/gene class information
ortho <- read.table("input/Ortho_Files/OrthoFinder_Combined_Parsed_Jul22.tsv",
                    sep = "\t", header = T)

ann <- read.table("input/Gene_Info/Amel_ImmAnnotations_Orthogroup_Jun22.tsv",
                  header = T, sep = "\t")
ann.slim <- ann[,c(2:3)]

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

ind.res$OrthoGroup[ind.res$OrthoGroup == "Unassigned" &
                     ind.res$Source == "Amel"] <- "Unassigned_Amel"
ind.res$OrthoGroup[ind.res$OrthoGroup == "Unassigned" &
                     ind.res$Source == "Bter" | ind.res$Source == "Bter_12"] <- "Unassigned_Bter"
ind.res$OrthoGroup[ind.res$OrthoGroup == "Unassigned" &
                     ind.res$Source == "Caus" | ind.res$Source == "Caus_12"] <- "Unassigned_Caus"
ind.res$OrthoGroup[ind.res$OrthoGroup == "Unassigned" &
                     ind.res$Source == "Plan" | ind.res$Source == "Plan_12"] <- "Unassigned_Plan"
summary(as.factor(ind.res$OrthoGroup))




#make matrix for upsetR
res.mat <- as.data.frame(unique(ind.res$OrthoGroup))
names(res.mat) <- "OrthoGroup"

runs <- unique(ind.res$Source)
for ( i in 1:length(runs)){
  x <- unique(ind.res$OrthoGroup[ind.res$Source == paste(runs[i])])
  res.mat[,i+1] <- ifelse(res.mat$OrthoGroup %in% x, 1, 0)
  names(res.mat)[i+1] <- paste(runs[i])
}
head(res.mat)

#immune annotation ? 
for(i in 1:nrow(res.mat)){
  x <- unique(ind.res$Function[ind.res$OrthoGroup == res.mat$OrthoGroup[i]])
  res.mat$Immune[i] <- ifelse(x == "Immune", paste(TRUE), paste(FALSE))
  res.mat$PutativeImmune[i] <- ifelse(x == "Putative Immune", paste(TRUE), paste(FALSE))
}
head(res.mat)

#Plot UpsetR

upset(res.mat,
      query.legend = "top",
      mb.ratio = c(0.55, 0.45),
      nsets = 6)

names(res.mat) <- c("Orthogroup", "Apis mellifera", "Bombus terrestris", 
                    "Ceratina australensis", "Polistes lanio", "Immune",
                    "PutativeImmune")

upset(res.mat,
      query.legend = "top",
      mb.ratio = c(0.7, 0.3),
      sets = c("Apis mellifera", "Bombus terrestris", "Ceratina australensis",
               "Polistes lanio"),
      sets.bar.color = c("#ffa500", "#005aff", "#00ffa5", "#ff005a"),
      order.by = c("degree"), keep.order = T,
      )

####Venn Diagram of all DE Genes####
#Requires all results (not just significant)

ind.list <- vector(mode = "list", length = length(targ))
for ( i in 1:length(ind.list)){
  file <- list.files(path = paste(targ[i], sep ="/"), 
                     pattern = 'DeSeqResults_All.tsv')
  spec <- strsplit(targ[i], "/")[[1]][3]
  ind.list[[i]] <- read.table(paste(targ[i], file, sep = "/"),
                              header = T, sep = "\t", quote = "")
  ind.list[[i]]$Source <- paste(spec)
}
ind.res2 <- bind_rows(ind.list)
#remove omitted genes
ind.res2 <- ind.res2[!is.na(ind.res2$padj),]

#add orthogroup info
#NB: Apideacins between Amel and Bter are in separate orthogroups
#For the purposes of this, I will be combining them into a group 
unique(ind.res2$Gene[grep("apidaecin", ind.res2$Desc)])
ortho[40558,] <- c("OG0011065", "Amel", "Apid1")
tail(ortho)

ann.slim[12123, 1] <- "Effector"
ann.slim[12123, 2] <- "OG0011065"
tail(ann.slim)

j <- 1
for(i in 1:nrow(ind.res2)){
  check <- ind.res2$Gene[i] %in% ortho$GeneID
  if (check == FALSE){
    ind.res2$OrthoGroup[i] <- paste("Unassigned", j, sep ="_")
    j <- j+1
  } else {
    ind.res2$OrthoGroup[i] <- unique(ortho$Orthogroup[ortho$GeneID ==
                                                ind.res2$Gene[i]])
  }
}

##Wound
#upregulated log2FC > 1

wound.up <- ind.res2[ind.res2$Contrast == "Treatment_Wound_vs_Naive" &
                       ind.res2$log2FoldChange > 1,]
#remove unassigned orthogroups
wound.up2 <- wound.up[grep("Unassigned", wound.up$OrthoGroup, invert = T),]
summary(as.factor(wound.up2$Source))

w.u.venn <- list("A.mellifera" = 
                   wound.up$OrthoGroup[wound.up$Source == "Amel"],
     "B.terrestris" = wound.up$OrthoGroup[wound.up$Source == "Bter"],
     "C.australensis" = wound.up$OrthoGroup[wound.up$Source == "Caus"],
     "P.lanio" = wound.up$OrthoGroup[wound.up$Source == "Plan"])

ggvenn(w.u.venn,
       set_name_size = 3,
       show_percentage = F,
       #fill_color = c("yellow", "blue", "green")
       )

#As I can't italicise or otherwise play with the titles here, I'm going to try another
#package: ggVennDiagram

venn <- Venn(w.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
  labs(title = "Wound") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#Plot2: colours in red for upregulated
ggVennDiagram(w.u.venn,
              label = "count") +
  scale_colour_manual(values = c("orange", "orange", "orange", "orange")) +
  scale_fill_gradient(low = "#ffff9d", high = "#ff0000") +
  labs(title = "Wound") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##Wound Down
#upregulated log2FC < 1

wound.dn <- ind.res2[ind.res2$Contrast == "Treatment_Wound_vs_Naive" &
                       ind.res2$log2FoldChange < -1,]


w.d.venn <- list("A.mellifera" = 
                   wound.dn$OrthoGroup[wound.dn$Source == "Amel"],
                 "B.terrestris" = wound.dn$OrthoGroup[wound.dn$Source == "Bter"],
                 "C.australensis" = wound.dn$OrthoGroup[wound.dn$Source == "Caus"],
                 "P.lanio" = wound.dn$OrthoGroup[wound.dn$Source == "Plan"])

venn <- Venn(w.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

#Plot2: colours in red for upregulated
ggVennDiagram(w.d.venn,
              label = "count") +
  scale_colour_manual(values = c("orange", "orange", "orange", "orange")) +
  scale_fill_gradient(low = "#ffff9d", high = "#ff0000") +
  labs(title = "test")

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#d8d8ff", high = "#0000ff")

##Gram Pos
#upregulated log2FC > 1

pos.up <- ind.res2[ind.res2$Contrast == "Treatment_GramPos_vs_Naive" &
                     ind.res2$log2FoldChange > 1,]


p.u.venn <- list("A.mellifera" = 
                   pos.up$OrthoGroup[pos.up$Source == "Amel"],
                 "B.terrestris" = pos.up$OrthoGroup[pos.up$Source == "Bter"],
                 "C.australensis" = pos.up$OrthoGroup[pos.up$Source == "Caus"],
                 "P.lanio" = pos.up$OrthoGroup[pos.up$Source == "Plan"])

venn <- Venn(p.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

#Plot2: colours in red for upregulated

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##pos Down
#upregulated log2FC < 1

pos.dn <- ind.res2[ind.res2$Contrast == "Treatment_GramPos_vs_Naive" &
                     ind.res2$log2FoldChange < -1,]


p.d.venn <- list("A.mellifera" = 
                   pos.dn$OrthoGroup[pos.dn$Source == "Amel"],
                 "B.terrestris" = pos.dn$OrthoGroup[pos.dn$Source == "Bter"],
                 "C.australensis" = pos.dn$OrthoGroup[pos.dn$Source == "Caus"],
                 "P.lanio" = pos.dn$OrthoGroup[pos.dn$Source == "Plan"])

venn <- Venn(p.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

#Plot2: colours in red for upregulated
p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#00bfff", high = "#0000ff")

##Gram Neg
#upregulated log2FC > 1

neg.up <- ind.res2[ind.res2$Contrast == "Treatment_GramNeg_vs_Naive" &
                     ind.res2$log2FoldChange > 1,]


n.u.venn <- list("A.mellifera" = 
                   neg.up$OrthoGroup[neg.up$Source == "Amel"],
                 "B.terrestris" = neg.up$OrthoGroup[neg.up$Source == "Bter"],
                 "C.australensis" = neg.up$OrthoGroup[neg.up$Source == "Caus"],
                 "P.lanio" = neg.up$OrthoGroup[neg.up$Source == "Plan"])

venn <- Venn(n.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##neg.Down
#upregulated log2FC < 1

neg.dn <- ind.res2[ind.res2$Contrast == "Treatment_GramNeg_vs_Naive" &
                     ind.res2$log2FoldChange < -1,]


n.d.venn <- list("A.mellifera" = 
                   neg.dn$OrthoGroup[neg.dn$Source == "Amel"],
                 "B.terrestris" = neg.dn$OrthoGroup[neg.dn$Source == "Bter"],
                 "C.australensis" = neg.dn$OrthoGroup[neg.dn$Source == "Caus"],
                 "P.lanio" = neg.dn$OrthoGroup[neg.dn$Source == "Plan"])

venn <- Venn(n.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#d8d8ff", high = "#0000ff")



####Venn with Apis/Bombus ######

##Wound
#upregulated log2FC > 1

wound.up <- ind.res2[ind.res2$Contrast == "Treatment_Wound_vs_Naive" &
                       ind.res2$log2FoldChange > 1,]

w.u.venn <- list("A.mellifera" = 
                   wound.up$OrthoGroup[wound.up$Source == "Amel"],
                 "B.terrestris" = wound.up$OrthoGroup[wound.up$Source == "Bter"])

venn <- Venn(w.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
p + scale_fill_manual(values = c("#ffa500", "#805380", "#005aff"))

#Plot2: colours in red for upregulated
p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##Wound Down
#upregulated log2FC < 1

wound.dn <- ind.res2[ind.res2$Contrast == "Treatment_Wound_vs_Naive" &
                       ind.res2$log2FoldChange < -1,]


w.d.venn <- list("A.mellifera" = 
                   wound.dn$OrthoGroup[wound.dn$Source == "Amel"],
                 "B.terrestris" = wound.dn$OrthoGroup[wound.dn$Source == "Bter"])

venn <- Venn(w.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
p + scale_fill_manual(values = c("#ffa500", "#805380", "#005aff"))

#Plot2: colours in red for upregulated
p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#d8d8ff", high = "#0000ff")

##Gram Pos
#upregulated log2FC > 1

pos.up <- ind.res2[ind.res2$Contrast == "Treatment_GramPos_vs_Naive" &
                     ind.res2$log2FoldChange > 1,]


p.u.venn <- list("A.mellifera" = 
                   pos.up$OrthoGroup[pos.up$Source == "Amel"],
                 "B.terrestris" = pos.up$OrthoGroup[pos.up$Source == "Bter"])

venn <- Venn(p.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
p + scale_fill_manual(values = c("#ffa500", "#805380", "#005aff"))

#Plot2: colours in red for upregulated

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##pos Down
#upregulated log2FC < 1

pos.dn <- ind.res2[ind.res2$Contrast == "Treatment_GramPos_vs_Naive" &
                     ind.res2$log2FoldChange < -1,]


p.d.venn <- list("A.mellifera" = 
                   pos.dn$OrthoGroup[pos.dn$Source == "Amel"],
                 "B.terrestris" = pos.dn$OrthoGroup[pos.dn$Source == "Bter"])

venn <- Venn(p.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
p + scale_fill_manual(values = c("#ffa500", "#805380", "#005aff"))

#Plot2: colours in red for upregulated
p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#00bfff", high = "#0000ff")

##Gram Neg
#upregulated log2FC > 1

neg.up <- ind.res2[ind.res2$Contrast == "Treatment_GramNeg_vs_Naive" &
                     ind.res2$log2FoldChange > 1,]


n.u.venn <- list("A.mellifera" = 
                   neg.up$OrthoGroup[neg.up$Source == "Amel"],
                 "B.terrestris" = neg.up$OrthoGroup[neg.up$Source == "Bter"])

venn <- Venn(n.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
p + scale_fill_manual(values = c("#ffa500", "#805380", "#005aff"))

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##neg.Down
#upregulated log2FC < 1

neg.dn <- ind.res2[ind.res2$Contrast == "Treatment_GramNeg_vs_Naive" &
                     ind.res2$log2FoldChange < -1,]


n.d.venn <- list("A.mellifera" = 
                   neg.dn$OrthoGroup[neg.dn$Source == "Amel"],
                 "B.terrestris" = neg.dn$OrthoGroup[neg.dn$Source == "Bter"])

venn <- Venn(n.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
p + scale_fill_manual(values = c("#ffa500", "#805380", "#005aff"))

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#d8d8ff", high = "#0000ff")



####Venn Diagram of Significant DE Genes####
#Requires all results (not just significant)

ind.res3 <- ind.res2[ind.res2$padj < 0.1,]

###Wound
#upregulated log2FC > 1

wound.up <- ind.res3[ind.res3$Contrast == "Treatment_Wound_vs_Naive" &
                       ind.res3$log2FoldChange > 1,]
#remove unassigned orthogroups
wound.up2 <- wound.up[grep("Unassigned", wound.up$OrthoGroup, invert = T),]
summary(as.factor(wound.up2$Source))

w.u.venn <- list("A.mellifera" = 
                   wound.up$OrthoGroup[wound.up$Source == "Amel"],
                 "B.terrestris" = wound.up$OrthoGroup[wound.up$Source == "Bter"],
                 "C.australensis" = wound.up$OrthoGroup[wound.up$Source == "Caus"],
                 "P.lanio" = wound.up$OrthoGroup[wound.up$Source == "Plan"])

venn <- Venn(w.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

#Plot2: colours in red for upregulated
p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##Wound Down
#upregulated log2FC < 1

wound.dn <- ind.res3[ind.res3$Contrast == "Treatment_Wound_vs_Naive" &
                       ind.res3$log2FoldChange < -1,]


w.d.venn <- list("A.mellifera" = 
                   wound.dn$OrthoGroup[wound.dn$Source == "Amel"],
                 "B.terrestris" = wound.dn$OrthoGroup[wound.dn$Source == "Bter"],
                 "C.australensis" = wound.dn$OrthoGroup[wound.dn$Source == "Caus"],
                 "P.lanio" = wound.dn$OrthoGroup[wound.dn$Source == "Plan"])

venn <- Venn(w.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

#Plot2: colours in red for upregulated
ggVennDiagram(w.d.venn,
              label = "count") +
  scale_colour_manual(values = c("orange", "orange", "orange", "orange")) +
  scale_fill_gradient(low = "#ffff9d", high = "#ff0000") +
  labs(title = "test")

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#d8d8ff", high = "#0000ff")

##Gram Pos
#upregulated log2FC > 1

pos.up <- ind.res3[ind.res3$Contrast == "Treatment_GramPos_vs_Naive" &
                     ind.res3$log2FoldChange > 1,]


p.u.venn <- list("A.mellifera" = 
                   pos.up$OrthoGroup[pos.up$Source == "Amel"],
                 "B.terrestris" = pos.up$OrthoGroup[pos.up$Source == "Bter"],
                 "C.australensis" = pos.up$OrthoGroup[pos.up$Source == "Caus"],
                 "P.lanio" = pos.up$OrthoGroup[pos.up$Source == "Plan"])

venn <- Venn(p.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

#Plot2: colours in red for upregulated

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##pos Down
#upregulated log2FC < 1

pos.dn <- ind.res3[ind.res3$Contrast == "Treatment_GramPos_vs_Naive" &
                     ind.res3$log2FoldChange < -1,]


p.d.venn <- list("A.mellifera" = 
                   pos.dn$OrthoGroup[pos.dn$Source == "Amel"],
                 "B.terrestris" = pos.dn$OrthoGroup[pos.dn$Source == "Bter"],
                 "C.australensis" = pos.dn$OrthoGroup[pos.dn$Source == "Caus"],
                 "P.lanio" = pos.dn$OrthoGroup[pos.dn$Source == "Plan"])

venn <- Venn(p.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

#Plot2: colours in red for upregulated
p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#00bfff", high = "#0000ff")

##Gram Neg
#upregulated log2FC > 1

neg.up <- ind.res3[ind.res3$Contrast == "Treatment_GramNeg_vs_Naive" &
                     ind.res3$log2FoldChange > 1,]


n.u.venn <- list("A.mellifera" = 
                   neg.up$OrthoGroup[neg.up$Source == "Amel"],
                 "B.terrestris" = neg.up$OrthoGroup[neg.up$Source == "Bter"],
                 "C.australensis" = neg.up$OrthoGroup[neg.up$Source == "Caus"],
                 "P.lanio" = neg.up$OrthoGroup[neg.up$Source == "Plan"])

venn <- Venn(n.u.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#ffff9d", high = "#ff0000")

##neg.Down
#upregulated log2FC < 1

neg.dn <- ind.res3[ind.res3$Contrast == "Treatment_GramNeg_vs_Naive" &
                     ind.res3$log2FoldChange < -1,]


n.d.venn <- list("A.mellifera" = 
                   neg.dn$OrthoGroup[neg.dn$Source == "Amel"],
                 "B.terrestris" = neg.dn$OrthoGroup[neg.dn$Source == "Bter"],
                 "C.australensis" = neg.dn$OrthoGroup[neg.dn$Source == "Caus"],
                 "P.lanio" = neg.dn$OrthoGroup[neg.dn$Source == "Plan"])

venn <- Venn(n.d.venn)
data <- process_data(venn)

#Plot1: with colours per species
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
                                 "#5555aa", "#8000ad", "#00ffa5", "#808080", "#ff005a"))

p <- ggplot() + 
  #region count layer
  geom_sf(aes(fill = count), data = venn_region(data), show.legend = FALSE) +
  #edge layer
  geom_sf(color="black", size = .5, data = venn_setedge(data), show.legend = FALSE) + 
  #label layer
  geom_sf_text(aes(label = name), fontface = "italic", data = venn_setlabel(data)) +  
  #region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()
#manually made from a colour mixing site
p + scale_fill_gradient(low = "#d8d8ff", high = "#0000ff")

######8th August#####



