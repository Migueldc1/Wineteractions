# Author: M, de Celis Rodriguez
# Date: 29/08/2022
# Project: Wineteractions - Metatranscriptomic Count Reading

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#
#### LIBRARIES ####
library(reshape2)


#
#### READ COUNT TABLES ####
count.list <- list()
i <- 0

for (count in list.files("Inputs/RNA_data/Count/")) {
  if (grepl(".summary", count, fixed = TRUE) == FALSE) {
    
    i <- i+1
    
    table <- read.table(paste("Inputs/RNA_data/Count/", count, sep = "/"))
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

for (annot in list.files("Inputs/RNA_data/Annotation")) {
  
  i <- i+1
  
  table <- read.table(paste("Inputs/RNA_data/Annotation/", annot, sep = "/"), 
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

#### SAVE DATA ####
write.csv(count_df, "Outputs/count_df.csv", row.names = FALSE)

#