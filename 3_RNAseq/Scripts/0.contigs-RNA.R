# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - Metatranscriptomics Data Massaging

# Set the project location as working directory
setwd("C:/Users/Migueldc/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#
#### LIBRARIES ####
library(ggplot2)
library(reshape2)

#load("comm.analysis_wnt.Rdata")

#
#### FUNCTIONS ####

# Open Count_dfs from featureCounts
counts2list <- function(path, df.ext) {
  
  open2list <- list()
  i <- 0
  
  for (dfs in list.files(path)) {
    if (grepl(df.ext, dfs, fixed = TRUE) == FALSE) {
      
      i <- i+1
      
      table <- read.table(paste(path, dfs, sep = "/"))
      colnames(table) <- table[1,]
      table <- table[-1,]
      colnames(table)[ncol(table)] <- gsub(".txt", "", dfs)
      open2list[[i]] <- table
      names(open2list)[i] <- gsub(".txt", "", dfs)
      
    }
  }
  
  return(open2list)
  
}

# Open Taxonomy assignation from metaeuk
tax2list <- function(path, df.ext){

  open2list <- list()
  i <- 0
  
  for (annot in list.files(path)) {
    
    for (annot.file in list.files(paste(path, annot, sep = "/"))) {
      
      if (grepl(df.ext, annot.file, fixed = TRUE) == TRUE) {
        
        i <- i+1
        
        table <- read.table(paste(path, annot, annot.file, sep = "/"), 
                            sep = "\t", fill = TRUE)
        
        colnames(table) <- c("IDs", "taxid", "id.level", "Best.taxonomy", "Taxonomy")
        open2list[[i]] <- table
        
        names(open2list)[i] <- annot
        
      }
    }
  }
  
  for (n in 1:length(open2list)) {
    
    table <- open2list[[n]]
    table <- unique(cbind.data.frame(Geneid = apply(as.matrix(table[,1]), 1, 
                                                    function(x) unlist(strsplit(x, "|", fixed = TRUE))[1]),
                                     
                                     Genus = apply(as.matrix(table[,5]), 1, 
                                                   function(x) unlist(strsplit(x, ";", fixed = TRUE))[12])))
    
    table <- table[!duplicated(table$Geneid),]
    
    table[is.na(table)] <- "Unidentified"
    table$Genus <- gsub("g_", "", table$Genus)
    
    open2list[[n]] <- table
    
  }
  
  return(open2list)  

}

# Merge Count and Tax lists

tax.x.count <- function(counts, taxs) {
  
  tax.x.count.list <- list()
  
  for (n in 1:length(counts)) {
    
    table <- merge(counts[[n]], taxs[[n]], by = "Geneid")[,8:7]
    
    table <- aggregate(as.numeric(table[,2]), list(table[,1]), sum)
    colnames(table) <- c("Genus", names(taxs)[n])
    
    tax.x.count.list[[n]] <- table
    names(tax.x.count.list)[n] <- names(taxs)[n]
    
  }
  
  return(tax.x.count.list)
  
}

#
#### DATA LOADING ####

# Counts
path <- "C:/Users/Migueldc/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/RNAseq Experiment/RNAseq/Count/UniprotKB"
count.list <- counts2list(path, ".summary")

# Taxonomy
path <- "C:/Users/Migueldc/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/RNAseq Experiment/RNAseq/Annotation/metaeuk"
tax.list <- tax2list(path, "_tax_per_pred.tsv")

#
#### PREPARE COUNT TABLES FOR TAXONOMIC PROFILING ####

tax_count.list <- tax.x.count(count.list, tax.list)

tax_count.df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Genus", all = TRUE), tax_count.list)

tax_count.df <- tax_count.df[-1,]
tax_count.df[is.na(tax_count.df)] <- 0

row.names(tax_count.df) <- tax_count.df[,1]
tax_count.df <- tax_count.df[,-1]

tax_count.t.df <- apply(tax_count.df, 2, function(x) x/sum(x))

#
#### PREPARE COUNT TABLES FOR DIFFERENTIAL GENE EXPRESSION ANALYSIS ####

path <- "C:/Users/Migueldc/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/RNAseq Experiment/RNAseq/Count/eggnog"
count.list <- counts2list(path, ".summary")

tpm.list <- list()

for (n in 1:length(count.list)) {
  
  table <- count.list[[n]][,c("Geneid", "Length", names(count.list)[n])]
  
  table[,2] <- as.numeric(table[,2])
  table[,3] <- as.numeric(table[,3])
  
  table <- cbind.data.frame(table[,-3], 
                            TPM = (table[,3]/table[,2])/sum(((table[,3])/table[,2]))*10^6)
  
  table <- subset(table, TPM >= 1)
  colnames(table)[3] <- names(count.list)[n]
  
  tpm.list[[n]] <- table
  names(tpm.list)[n] <- names(count.list)[n]
  
  
}

path <- "C:/Users/Migueldc/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/RNAseq Experiment/RNAseq/Annotation/eggnog"
attributes <- c("em_target", "em_Preferred_name", "em_KEGG_ko", "em_BRITE", "em_KEGG_Pathway", "em_GOs")

annot.list <- list()
i <- 0

for (annot in list.files(path)) {
  
  i <- i+1
  
  table <- read.table(paste(path, annot, paste(annot, ".emapper.decorated.gff", sep = ""), 
                            sep = "/"),
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
  names(annot.list)[i] <- annot
  
  
}

ko_tpm.list <- list()

for (n in 1:length(tpm.list)) {
  
  table <- tpm.list[[n]][,c(1,3)]
  colnames(table)[1] <- "target"
  
  table_anot <- annot.list[[n]][,c(1,3)]
  
  table <- merge(table, table_anot, by = "target", all.x = FALSE)
  
  ko_tpm.list[[n]] <- table
  names(ko_tpm.list)[n] <- names(tpm.list)[n]
  
  
}

comparison.list <- list()

for (n in 1:length(ko_tpm.list)) {
  
  table <- ko_tpm.list[[n]]
  table <- table[complete.cases(table),]
  
  table <- aggregate(table[,2], list(table$KEGG_ko), sum)
  colnames(table)[1] <- "KEGG_ko"
  colnames(table)[2] <- names(ko_tpm.list)[n]
  
  comparison.list[[n]] <- table
  names(comparison.list)[n] <- names(ko_tpm.list)[n]
  
  
}


tpm_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), comparison.list)

















#### DATA EXPORTATION ####
write.table(tax_count.t.df, "Inputs/contig_tax.txt", sep = "\t")

write.table(tpm_df, "Inputs/tpm_ko.txt", sep = "\t")





#

