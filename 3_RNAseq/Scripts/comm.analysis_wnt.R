# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - ITS sequence analysis of T1-SGM samples

# Set the project location as working directory
setwd("C:/Users/Miguel de Celis/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")

#
#### LIBRARIES ####
library(RColorBrewer)
library(ggplot2)
library(reshape2)

setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/Metatrans_experiment/Community/")
#load("comm.analysis_wnt.Rdata")

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
asv_GM <- readRDS("Outputs/ASV_t1-GM.rds")
asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x))

tax_GM <- readRDS("Outputs/tax_t1-GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- gsub("^[a-z]__", "", as.matrix(tax_GM))
tax_GM <- as.data.frame(set_unid(tax_GM))

asv.t_GM.p <- melt(asv.t_GM)
colnames(asv.t_GM.p) <- c("Id", "Seq_ID", "value")
asv.t_GM.p <- merge(asv.t_GM.p, tax_GM[,c(6,8)])
asv.t_GM.p <- merge(asv.t_GM.p, sample_df[c(1,2)], by = "Seq_ID")
asv.t_GM.p <- asv.t_GM.p[,c(2,5,3,4)]

## SGM
asv_SGM <- readRDS("Outputs/ASV_t1-SGM.rds")
asv_SGM <- asv_SGM[-1,]
asv.t_SGM <- apply(asv_SGM, 1, function(x) x/sum(x))

tax_SGM <- readRDS("Outputs/tax_t1-SGM.rds")
tax_SGM <- cbind.data.frame(tax_SGM, Id = row.names(tax_SGM))
tax_SGM <- gsub("^[a-z]__", "", as.matrix(tax_SGM))
tax_SGM <- as.data.frame(set_unid(tax_SGM))

asv.t_SGM.p <- melt(asv.t_SGM)
colnames(asv.t_SGM.p) <- c("Id", "Sample_ID", "value")
asv.t_SGM.p <- merge(asv.t_SGM.p, tax_SGM[,c(6,8)])

## RNA
asv_RNA <- readRDS("Outputs/ASV_t1-RNA.rds")
asv.t_RNA <- apply(asv_RNA, 1, function(x) x/sum(x))

tax_RNA <- readRDS("Outputs/tax_t1-RNA.rds")
tax_RNA <- cbind.data.frame(tax_RNA, Id = row.names(tax_RNA))

asv.t_RNA.p <- melt(asv.t_RNA)
colnames(asv.t_RNA.p) <- c("Id", "Sample_ID", "value")
asv.t_RNA.p <- merge(asv.t_RNA.p, tax_RNA[,c(1,3)])

## SAMPLE DATA
sample_df <- read.table("Inputs/sample_df.txt", sep = "\t", header = TRUE)
row.names(sample_df) <- sample_df$Sample_ID
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))

#
#### TAXONOMIC EXPLORATION ####
## Genus

asv.t_plot <- rbind.data.frame(cbind(asv.t_GM.p, Study = "GM"),
                               cbind(asv.t_SGM.p, Study = "SGM"),
                               cbind(asv.t_RNA.p, Study = "RNA"))

asv.t_plot <- aggregate(asv.t_plot$value, list(asv.t_plot$Sample_ID, asv.t_plot$Genus, asv.t_plot$Study), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "Study", "value")
asv.t_plot$Genus[asv.t_plot$value < 0.05] <- "Other"

asv.t_plot <- aggregate(asv.t_plot$value, list(asv.t_plot$Sample_ID, asv.t_plot$Genus, asv.t_plot$Study), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "Study", "value")

orderG <- levels(factor(asv.t_plot$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

asv.t_plot <- merge(asv.t_plot, sample_df, by = "Sample_ID")
asv.t_plot$Study <- factor(asv.t_plot$Study, levels = c("GM", "SGM", "RNA"))
asv.t_plot$plot_ID <- paste(asv.t_plot$Farming, asv.t_plot$Condition, sep = " ")
asv.t_plot$Origin <- factor(asv.t_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))

ggplot(asv.t_plot, 
       aes(x = plot_ID, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = c("#ffff33", "#caf55d", "#f58d5d", "#5df5cc", "#cc3939",
                                               "#b33cb5", "#9e66d1", "#cab2d6", "#8da0cb", "#d4eb26",
                                               "#b15928", "#1b9e77", "#d4556a", "#e6b42c", "#2c3ce6",
                                               "#f78e4d", "#e5d8bd", "#666666", "#6a3d9a", "#bf5b17")) +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black")) +
  xlab("Sample") + ylab("Abundance") +
  guides(fill = guide_legend(nrow = 3)) + facet_wrap(~Study + Origin, nrow = 3)

#
#### SAVE DATA ####
save.image("comm.analysis_wnt.Rdata")

#