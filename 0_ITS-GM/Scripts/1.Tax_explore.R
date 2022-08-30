# Author: M, de Celis Rodriguez
# Date: 13/07/2022
# Project: Wineteractions - ITS sequence analysis of GM samples

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/0_ITS-GM")

#
#### LIBRARIES ####
library(RColorBrewer)
library(ggplot2)
library(reshape2)

#
#### FUNCTIONS ####
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

set_unid <- function(tax_df) {
  
  for (i in 1:ncol(tax_df)) {
    if (i != 7) {
      tax_df[,i] <- ifelse(is.na(tax_df[,i]),
                           ifelse(grepl("Unidentified", tax_df[,i-1]) == TRUE,
                                  tax_df[,i-1],
                                  paste("Unidentified", tax_df[,i-1], sep = " ")),
                           tax_df[,i])
    } else{
      tax_df[,i] <- ifelse(is.na(tax_df[,i]),
                           ifelse(grepl("Unidentified", tax_df[,i-1]) == TRUE,
                                  tax_df[,i-1],
                                  "sp."),
                           tax_df[,i])
      
    }
    
  }
  
  return(tax_df)
  
}

#
#### DATA LOADING ####
## GM
asv_GM <- readRDS("Outputs/ASV_GM.rds")
asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x))

tax_GM <- readRDS("Outputs/tax_GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- gsub("^[a-z]__", "", as.matrix(tax_GM))

tax_GM.un <- as.data.frame(set_unid(tax_GM))
tax_GM[is.na(tax_GM)] <- "Unidentified"

## SAMPLE DATA
sample_df <- read.table("Inputs/sample_GM.txt", sep = "\t", header = TRUE)
row.names(sample_df) <- sample_df$Seq_ID
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df <- cbind.data.frame(Sample_ID = paste(sample_df$Origin, sample_df$Farming, sample_df$Condition, sep = "-"),
                              sample_df)

#
#### TAXONOMIC EXPLORATION ####
asv.t_GM.p <- melt(asv.t_GM)
colnames(asv.t_GM.p) <- c("Id", "Seq_ID", "value")
asv.t_GM.p <- merge(asv.t_GM.p, tax_GM[,c(6,8)], by = "Id")

asv.t_GM.p <- merge(asv.t_GM.p, sample_df[c(1:6)], by = "Seq_ID")

## Genus

asv.t_plot <- aggregate(asv.t_GM.p$value, list(asv.t_GM.p$Sample_ID, asv.t_GM.p$Genus, asv.t_GM.p$Stage), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "Stage", "value")
asv.t_plot$Genus[asv.t_plot$value < 0.05] <- "Other"

asv.t_plot <- aggregate(asv.t_plot$value, list(asv.t_plot$Sample_ID, asv.t_plot$Genus, asv.t_plot$Stage), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "Stage", "value")

orderG <- levels(factor(asv.t_plot$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

asv.t_plot <- merge(asv.t_plot, sample_df[,1:6], by = c("Sample_ID", "Stage"))
asv.t_plot$Origin <- factor(asv.t_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))

asv.t_plot$Sample_name <- paste(asv.t_plot$Farming, asv.t_plot$Condition, sep = "-")
asv.t_plot$Sample_name <- factor(asv.t_plot$Sample_name, levels = c("CONV-Control", "CONV-18C", "CONV-NH4", "CONV-SO2",
                                                                    "ECO-Control", "ECO-18C", "ECO-NH4", "ECO-SO2"))

set.seed(123)
ggplot(asv.t_plot, 
       aes(x = Sample_name, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = sample(col_vector, 18)) +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black")) +
  xlab("Sample") + ylab("Abundance") +
  guides(fill = guide_legend(nrow = 3)) + facet_wrap(~Stage + Origin, nrow = 3)

#