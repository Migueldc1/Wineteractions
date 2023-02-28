# Author: M, de Celis Rodriguez
# Date: 23/02/2023
# Project: Wineteractions - Yeast community description

library(reshape2)
library(ggplot2)
library(hillR)
library(vegan)
library(cowplot)


rm(list=ls()) #Clear R environment

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/")

#
#### CUSTOM FUNCTIONS ####
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

## SAMPLE DATA
sample_df <- read.table("Data/Metadata/sample_GM.txt", sep = "\t", header = TRUE, encoding = "UTF-8")
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Origin <- factor(sample_df$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
sample_df$Region <- factor(sample_df$Region, levels = c("Ribera del Guadiana", "Valdepeñas", "La Mancha", "Madrid", "Rioja"))

sample_df <- subset(sample_df, Stage == "0_initial")


## Community data
asv_GM <- readRDS("Data/Sequencing/Outputs/ASV_GM.rds")
asv_GM <- asv_GM[row.names(asv_GM) %in% sample_df$Seq_ID, ]
asv_GM <- asv_GM[, colSums(asv_GM != 0) > 0]

asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x))

tax_GM <- readRDS("Data/Sequencing/Outputs/tax_GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- gsub("^[a-z]__", "", as.matrix(tax_GM))
tax_GM <- tax_GM[row.names(tax_GM) %in% colnames(asv_GM), ]

tax_GM.un <- set_unid(tax_GM)
tax_GM[is.na(tax_GM)] <- "Unidentified"


## Colors
col_origin <- c("#2ac219", "#9fffd1", "#dba54d", "#ae0e36", "#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")

#
#### TAXONOMIC EXPLORATION ####
asv.t_GM.p <- melt(asv.t_GM)
colnames(asv.t_GM.p) <- c("Id", "Seq_ID", "value")
asv.t_GM.p <- merge(asv.t_GM.p, tax_GM[,c(6,8)], by = "Id")

asv.t_GM.p <- merge(asv.t_GM.p, sample_df[c(1:4)], by = "Seq_ID")

## Genus
asv.t_plot <- aggregate(asv.t_GM.p$value, list(asv.t_GM.p$Seq_ID, asv.t_GM.p$Genus), sum)
colnames(asv.t_plot) <- c("Seq_ID", "Genus", "value")
asv.t_plot$Genus[asv.t_plot$value < 0.05] <- "Other"

asv.t_plot <- aggregate(asv.t_plot$value, list(asv.t_plot$Seq_ID, asv.t_plot$Genus), sum)
colnames(asv.t_plot) <- c("Seq_ID", "Genus", "value")

orderG <- levels(factor(asv.t_plot$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

asv.t_plot <- merge(asv.t_plot, sample_df[,1:4], by = "Seq_ID")

asv.t_plot$Farming <- ifelse(asv.t_plot$Farming == "ECO", "Ecological", "Conventional")

gg.gen <- ggplot(asv.t_plot, 
       aes(x = Condition, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = c("#ccebc5", "#ffff33", "#e6ab02", "#80b1d3", "#e6f5c9", "#ff7f00", "#fff2ae",
                                               "#8da0cb", "#b15928", "#bebada", "#1b9e77", "#fbb4ae", "#6a3d9a", "#bf5b17")) +
  guides(fill = guide_legend(ncol = 1)) + 
  facet_grid(Farming ~ Origin) +
  xlab("") + ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "right", 
        legend.text.align = 0,
        aspect.ratio = 2,
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"))

gg.gen

#
#### ALPHA DIVERSITY ANALYSIS ####

alpha_GM <- cbind.data.frame(# Hill based taxonomic alpha diversity
                             Richness = hill_taxa(t(asv.t_GM), q = 0),
                             `Hill–Shannon` = hill_taxa(t(asv.t_GM), q = 1),
                             `Hill–Simpson` = hill_taxa(t(asv.t_GM), q = 2))

alpha_GM$Seq_ID <- row.names(alpha_GM)
alpha_GM <- merge(alpha_GM, sample_df, by = "Seq_ID")

alpha_GM.plot <- melt(alpha_GM[,1:7])

wilcox_df <- NULL
for (var in unique(alpha_GM.plot$variable)) {
  
  alpha_sub <- subset(alpha_GM.plot, variable == var)
  wilcox_df <- rbind(wilcox_df,
                     cbind.data.frame(variable = var, p.value = wilcox.test(alpha_sub[alpha_sub$Farming == "CONV", "value"],
                                                                            alpha_sub[alpha_sub$Farming == "ECO", "value"])$p.value))
  
}

alpha_GM.plot <- merge(alpha_GM.plot, wilcox_df, by = "variable")
alpha_GM.plot$p.value <- ifelse(alpha_GM.plot$p.value <= 0.001, "***",
                                ifelse(alpha_GM.plot$p.value <= 0.01, "**",
                                       ifelse(alpha_GM.plot$p.value <= 0.05, "*", "")))

alpha_GM.plot$Farming <- ifelse(alpha_GM.plot$Farming == "ECO", "Ecological", "Conventional")

## Management
gg.alpha_far <- ggplot(data = alpha_GM.plot, aes(x = Farming, y = value)) +
  geom_boxplot(size = 1) +
  geom_text(aes(x = 1.5, y = Inf, label = p.value), hjust = 0.5, vjust = 2, size = 8) +
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  ylab("") + xlab("") + 
  theme_bw() +
  theme(aspect.ratio = 2,
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, angle = -30, color = "black", vjust = 0.5, hjust = 0.20),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.alpha_far

#
#### BETA DIVERSITY - GENUS LEVEL ####

asv.t_gen <- aggregate(asv.t_GM, list(tax_GM.un[,"Genus"]), sum)
row.names(asv.t_gen) <- asv.t_gen[,1]
asv.t_gen <- asv.t_gen[,-1]

## Calculate Bray-Curtis dissimilarity and NMDS
bray_gen <- vegdist(t(asv.t_gen), method = "bray", binary = TRUE)

set.seed(1)
nMDS_gen <- metaMDS(bray_gen)
nMDS_gen$stress

nMDS_gen.plot <- as.data.frame(nMDS_gen[["points"]])
nMDS_gen.plot$Seq_ID <- rownames(nMDS_gen.plot)
nMDS_gen.plot <- merge(nMDS_gen.plot, sample_df, by = "Seq_ID")

gg.nmds.gen_bray <- ggplot(nMDS_gen.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, color = Origin, shape = Farming, size = Farming), stroke = 3) +
  scale_shape_manual(values = c(16, 5)) +
  scale_size_manual(values = c(4, 2), guide = "none") +
  scale_color_manual(values = col_origin) +
  annotate(geom = "text", x = Inf, y = Inf, label = paste("stress =", round(nMDS_gen$stress, 3)), size = 6, vjust = 1.2, hjust = 1.05) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme(aspect.ratio = 0.66,
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.nmds.gen_bray

set.seed(1)
adonis2(bray_gen ~ Origin*Farming, data = sample_df)

#
#### EXPORT FIGURE 1 ####

gg.figure1 <- plot_grid(gg.gen, plot_grid(gg.alpha_far, gg.nmds.gen_bray, rel_widths = c(1, 1.25), labels = c("B", "C"), 
                                          label_size = 18), 
                        ncol = 1, rel_heights = c(1.33, 1), labels = c("A"), label_size = 18)
gg.figure1

ggsave("Figures/Figure_1.png", gg.figure1, bg = "white", width = 14, height = 10)

