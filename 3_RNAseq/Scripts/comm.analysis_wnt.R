#### LIBRARIES ####
library(RColorBrewer)
library(dada2)
library(msa)
library(phangorn)
library(reshape2)
library(ggplot2)
library(hillR)
library(vegan)

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
#### DATA PREPROCESSING ####
asv_wnt <- readRDS("Inputs/ASV_wnt.rds")
asv_wnt <- asv_wnt[-1,]
asv.t_wnt <- apply(asv_wnt, 1, function(x) x/sum(x))

tax_wnt <- readRDS("Inputs/tax_wnt.rds")
tax_wnt <- cbind.data.frame(tax_wnt, Id = row.names(tax_wnt))
tax_wnt <- gsub("^[a-z]__", "", as.matrix(tax_wnt))

tax_wnt <- set_unid(tax_wnt)

sample_df <- read.table("Inputs/sample_df.txt", sep = "\t", header = TRUE)
row.names(sample_df) <- sample_df$Sample_ID
sample_df <- sample_df[sample_df$Sample_ID %in% row.names(asv_wnt),]
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))

#
#### Phylogenetic tree ####
seqs_wnt <- dada2::getSequences(colnames(asv_wnt))
names(seqs_wnt) <- seqs_wnt

mult_wnt <- msa(seqs_wnt, method = "ClustalW", type = "dna", order = "input")

phang.align_wnt <- as.phyDat(mult_wnt, type = "DNA", names = seqs_wnt)
dm_wnt <- dist.ml(phang.align_wnt)
treeNJ_wnt <- NJ(dm_wnt) # Note, tip order != sequence order
fit_wnt <- pml(treeNJ_wnt, data = phang.align_wnt)

fitGTR_wnt <- update(fit_wnt, k = 4, inv = 0.2)
fitGTR_wnt <- optim.pml(fitGTR_wnt, model = "GTR", optInv = TRUE, optGamma = TRUE,
                         rearrangement = "stochastic", control = pml.control(trace = 0))

#
#### TAXONOMIC EXPLORATION ####
##Genus
dataG_wnt <- merge(asv.t_wnt, tax_wnt, by = "row.names")[,c(2:60,66)]
dataG_wnt <- aggregate(dataG_wnt[,1:59], list(dataG_wnt$Genus), sum)

dataG_wnt <- melt(dataG_wnt)
colnames(dataG_wnt) <- c("Genus", "Sample_ID", "value")
dataG_wnt$Genus[dataG_wnt$value < 0.05] <- "Other"

glomG_wnt <- tax_glom(phyt_wnt, taxrank = 'Genus', NArm = FALSE)
dataG_wnt <- psmelt(glomG_wnt)
dataG_wnt$Genus <- as.character(dataG_wnt$Genus)
dataG_wnt$Genus[dataG_wnt$Abundance < 0.05] <- "Other"
dataG_wnt$Genus[is.na(dataG_wnt$Genus)] <- "Unidentified"

orderG <- levels(factor(dataG_wnt$Genus))
orderG <- orderG[! orderG %in% "Other"]
orderG <- append(orderG, "Other")

dataG_wnt <- aggregate(dataG_wnt$value, list(dataG_wnt$Sample_ID, dataG_wnt$Genus), sum)
colnames(dataG_wnt) <- c("Sample_ID", "Genus", "value")

dataG_wnt <- merge(dataG_wnt, sample_df, by = "Sample_ID")

set.seed(123)
ggplot(dataG_wnt, 
       aes(x = Sample_ID, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = sample(col_vector, 18)) +
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


#
#### ALPHA DIVERISTY: Dynamics ####
alpha_wnt <- cbind.data.frame(# Hill based taxonomic alpha diversity
                              t.q0 = hill_taxa(t(asv.t_wnt), q = 0),
                              t.q1 = hill_taxa(t(asv.t_wnt), q = 1),
                              t.q2 = hill_taxa(t(asv.t_wnt), q = 2),
                              
                              # Hill based phylogenetic alpha diversity
                              p.q0 = hill_phylo(t(asv.t_wnt), fitGTR_wnt$tree, q = 0),
                              p.q1 = hill_phylo(t(asv.t_wnt), fitGTR_wnt$tree, q = 1),
                              p.q2 = hill_phylo(t(asv.t_wnt), fitGTR_wnt$tree, q = 2))

alpha_wnt$Sample_ID <- row.names(alpha_wnt)
alpha_wnt <- merge(alpha_wnt, sample_df, by = "Sample_ID")

alpha_wnt.plot <- melt(alpha_wnt)
alpha_wnt.plot$diversity <- ifelse(startsWith(as.character(alpha_wnt.plot$variable), "t."),
                                   "Taxonomic", "Phylogenetic")
alpha_wnt.plot$diversity <- factor(alpha_wnt.plot$diversity, 
                                   levels = c("Taxonomic", "Phylogenetic"))

alpha_wnt.plot$q <- as.factor(substr(alpha_wnt.plot$variable, 
                                     nchar(as.character(alpha_wnt.plot$variable)),
                                     nchar(as.character(alpha_wnt.plot$variable))))

ggplot(subset(alpha_wnt.plot, variable == "t.q0"), 
       aes(x = Condition, y = value, color = Condition)) + 
  geom_boxplot(size = 1) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust = 1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Alpha-diversity")

#
#### BETA DIVERSITY ####
bray_wnt <- vegdist(t(asv.t_wnt), method = "bray")
nMDS_wnt <- metaMDS(bray_wnt)
nMDS_wnt$stress

nMDS_wnt.plot <- as.data.frame(nMDS_wnt$points)
nMDS_wnt.plot$Sample_ID <- row.names(nMDS_wnt.plot)
nMDS_wnt.plot <- merge(nMDS_wnt.plot, sample_df, by = "Sample_ID")

ggplot(nMDS_wnt.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, color = Condition, shape = Farming), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

#
#### SAVE DATA ####
save.image("comm.analysis_wnt.Rdata")

#