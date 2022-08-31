ortho <- read.table("input/Ortho_Files/OrthoFinder_Combined_Parsed_Aug22.tsv", 
                    header = T, sep = "\t")

ann <- read.table("input/Gene_Info/Amel_ImmAnnotations_Orthogroup_Aug22.tsv",
                  header = T, sep = "\t")

for(i in 1:nrow(ortho)){
      check <- ortho$Orthogroup[i] %in% ann$Orthogroup
      print(check)
      if (check == FALSE){
        ortho$Function[i] <- "Un-Annotated"
      } else {
        f <- unique(ann$Function[ann$Orthogroup == ortho$Orthogroup[i]])
        print(f)
        if (length(f) > 1){
          if(TRUE %in% grepl("Effector|Recognition|Signalling", f)){
            ortho$Function[i] <- "Immune"
          } else {
            if(TRUE %in% grepl("Putative Immune", f)){
              ortho$Function[i] <- "Putative Immune"
            } 
          }
        } else {
          f2 <- gsub("Effector", "Immune", 
                     gsub("Recognition", "Immune", 
                          gsub("Signalling", "Immune", 
                               gsub("Background", "None Immune", f))))
          print(f2)
          ortho$Function[i] <- paste(f2)
        }
      }
    }

tail(ortho)

write.table(ortho, "output/Synthesis/OrthoGroup_Annotated_Aug22.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)
