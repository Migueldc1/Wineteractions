# Author: M, de Celis Rodriguez
# Date: 29/08/2022
# Project: Wineteractions - Metatranscriptomic Count Reading

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#load("RNA_counts.Rdata")

#
#### LIBRARIES ####
library(reshape2)


#
#### READ COUNT TABLES ####
count.list <- list()
i <- 0

for (count in list.files("Inputs/RNA_data/Count/eggnog/")) {
  if (grepl(".summary", count, fixed = TRUE) == FALSE) {
    
    i <- i+1
    
    table <- read.table(paste("Inputs/RNA_data/Count/eggnog/", count, sep = "/"))
    colnames(table) <- table[1,]
    table <- table[-1,]
    colnames(table)[ncol(table)] <- gsub(".txt", "", count)
    count.list[[i]] <- table[,c(1,6,7)]
    names(count.list)[i] <- gsub(".txt", "", count)
    
  }
  
}

#
#### READ ANNOTATION TABLES ####

attributes <- c("em_target", "em_Preferred_name", "em_KEGG_ko", "em_BRITE", "em_KEGG_Pathway", "em_GOs")

annot.list <- list()
i <- 0

for (annot in list.files("Inputs/RNA_data/Annotation/eggnog/")) {
  
  i <- i+1
  
  table <- read.table(paste("Inputs/RNA_data/Annotation/eggnog/", annot, sep = "/"), 
                      header = FALSE, quote = "", sep = "\t")
  
  for (atr in attributes) {
    
    atr_df <- NULL
    
    for (n in 1:nrow(table)) {
      
      if (sum(grepl(paste("^", atr, sep = ""), unlist(strsplit(table[,9][n], ";")))) == 1) {
        
        atr_df <- c(atr_df, gsub(paste(atr, "=", sep = ""), "", 
                                 grep(paste("^", atr, sep = ""), 
                                      unlist(strsplit(table[,9][n], ";")), value = TRUE)))
        
      }else if (sum(grepl(paste("^", atr, sep = ""), unlist(strsplit(table[,9][n], ";")))) == 0) {
        
        atr_df <- c(atr_df, NA)
        
      }
      
    }
    
    table <- cbind.data.frame(table, atr_df)
    colnames(table)[ncol(table)] <- gsub("em_", "", atr)
    
  }
  
  table <- unique(table[,10:15])
  
  annot.list[[i]] <- table
  names(annot.list)[i] <- gsub(".emapper.decorated.gff", "", annot)
  
  
}

#
#### MERGE TABLES ####
ko_count.list <- list()

for (n in 1:length(count.list)) {
  
  table <- count.list[[n]][,c(1,3)]
  colnames(table)[1] <- "target"
  
  table_anot <- annot.list[[n]][,c(1,3)]
  
  table <- merge(table, table_anot, by = "target", all.x = FALSE)
  
  ko_count.list[[n]] <- table
  names(ko_count.list)[n] <- names(count.list)[n]
  
  
}

comparison.list <- list()

for (n in 1:length(ko_count.list)) {
  
  table <- ko_count.list[[n]]
  table <- table[complete.cases(table),]
  
  table <- aggregate(as.numeric(table[,2]), list(table$KEGG_ko), sum)
  colnames(table)[1] <- "KEGG_ko"
  colnames(table)[2] <- names(ko_count.list)[n]
  
  comparison.list[[n]] <- table
  names(comparison.list)[n] <- names(ko_count.list)[n]
  
  
}

count_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), comparison.list)
count_df[is.na(count_df)] <- 0

#
#### ANNOTATION TABLE ####

annot_df.list <- list()

for (n in 1:length(annot.list)) {
  
  table <- annot.list[[n]][,3:6]
  table <- table[complete.cases(table),]
  
  annot_df.list[[n]] <- table
  names(annot_df.list)[n] <- names(annot.list)[n]
  
  
}

annot_df <- NULL

for (n in 1:length(annot_df.list)) {
  
  table <- annot_df.list[[n]]
  annot_df <- rbind.data.frame(annot_df, table)
  
}

annot_df <- unique(annot_df)

annot.GO_df <- annot_df[annot_df$GOs != "",]

GO_data <- NULL

for (n in 1:nrow(annot.GO_df)) {
  
  GO_data <- rbind.data.frame(GO_data, 
                              cbind(gene_id = annot.GO_df$KEGG_ko[n],
                                    category = unlist(strsplit(as.character(annot.GO_df$GOs[n]),",", fixed = TRUE))))
  
}

GO_data <- unique(GO_data)

#
#### LENGTH TABLE ####
id_df <- NULL

for (n in 1:length(annot.list)) {
  
  table <- annot.list[[n]][,c(1,3)]
  id_df <- rbind.data.frame(id_df, table)
  
}

colnames(id_df)[1] <- "Geneid"
id_df <- id_df[complete.cases(id_df), ]

length_df <- NULL

for (n in 1:length(count.list)) {
  
  table <- count.list[[n]][,-3]
  length_df <- rbind.data.frame(length_df, table)
  
}

length_df <- merge(length_df, id_df, by = "Geneid", all = FALSE)

length_df <- as.data.frame(unique(length_df[,-1]))

length_df <- aggregate(as.numeric(length_df$Length), list(length_df$KEGG_ko), mean)
colnames(length_df) <- c("KEGG_ko", "Length")

#
#### SAVE DATA ####
write.csv(count_df, "Outputs/count.KO_df.csv", row.names = FALSE)
write.csv(GO_data, "Outputs/GO_data.csv", row.names = FALSE)
write.csv(length_df, "Outputs/length_df.csv", row.names = FALSE)

save.image("RNA_counts.Rdata")

#