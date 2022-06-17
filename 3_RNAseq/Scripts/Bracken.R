#### LIBRARIES ####

library(reshape2)
library(ggplot2)
library(RColorBrewer)

setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/RNA-seq/")

#
#### FUNCTIONS ####
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#
######################
brack_df <- read.table("Inputs/bracken.report.txt")
colnames(brack_df) <- gsub("\\.", "-", colnames(brack_df))

smr_df <- read.table("Inputs/SMR_log.txt")
colnames(smr_df) <- c("Total reads", "rRNA reads", "rRNA reads (%)", "non-rRNA reads", "non-rRNA reads (%)")

corr_df <- merge(cbind.data.frame(Unidentified = t(brack_df)[,78]), smr_df[,3, drop = FALSE], by = "row.names")
colnames(corr_df)[1] <- "Sample_ID"
row.names(corr_df) <- corr_df[,1]

cor.test(corr_df$Unidentified, corr_df$`rRNA reads (%)`)
plot(corr_df$Unidentified, corr_df$`rRNA reads (%)`)

brack_df$Taxonomy <- row.names(brack_df)
brack_tax <- melt(brack_df)
brack_tax <- cbind.data.frame(brack_tax, Genus = colsplit(brack_tax$Taxonomy, " ", c("Genus", "specie"))[,1])

brack_tax[brack_tax$value < 0.025 & brack_tax$Genus != "Unidentified", "Genus"] <- "Other"
brack_tax <- aggregate(brack_tax$value, list(brack_tax$Genus, brack_tax$variable), sum)
colnames(brack_tax) <- c("Genus", "Sample_ID", "value")

orderT <- levels(factor(brack_tax$Genus))
orderT <- orderT[! orderT %in% c("Other", "Unidentified")]
orderT <- append(orderT, c("Other", "Unidentified"))



set.seed(123)
ggplot(brack_tax, 
       aes(x = Sample_ID, y = value, fill = factor(Genus, levels = orderT))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = sample(col_vector, 19)) +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black")) +
  xlab("Sample") + ylab("Abundance") +
  guides(fill = guide_legend(nrow = 3))
