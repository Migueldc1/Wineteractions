# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - RNAseq analysis

library(DESeq2)
library(ggplot2)
library(reshape2)
library(limma)
library(ggforce)
library(clusterProfiler)
library(goseq)

rm(list = ls())

# Set the project location as working directory
setwd("C:/Users/Migueldc/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/3_RNAseq/")
load("DESeq_meta.Rdata")

#
#### DATA LOADING ####

## SAMPLE DATA
sample_df <- read.table("Inputs/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df <- cbind.data.frame(Sample_ID = sample_df[,1],
                              colsplit(sample_df[,1], "-",
                                       names = c("Origin", "Farming", "Condition")),
                              sample_df[,-1])
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))

## Bracken Taxonomy
asv_RNA <- readRDS("Outputs/ASV_t1-RNA.rds")
asv_RNA <- apply(asv_RNA, 1, function(x) x/sum(x))

tax_RNA <- readRDS("Outputs/tax_t1-RNA.rds")
tax_RNA <- cbind.data.frame(tax_RNA, Id = row.names(tax_RNA))

asv_RNA <- melt(asv_RNA)
colnames(asv_RNA) <- c("Id", "Sample_ID", "value")
asv_RNA <- merge(asv_RNA, tax_RNA[,c(1,3)])

## Dominant Genus
dom_sp <- aggregate(asv_RNA$value, list(asv_RNA$Sample_ID), max)
colnames(dom_sp) <- c("Sample_ID", "value")
dom_sp <- merge(dom_sp, asv_RNA[,2:4], by = c("Sample_ID", "value"))

sample_df <- merge(sample_df, dom_sp[,c(1,3)], by = "Sample_ID")
sample_df$Genus <- relevel(as.factor(sample_df$Genus), ref = "Saccharomyces")

## COUNT DATA
countData <- read.csv("Outputs/count.KO_df.csv", header = TRUE, sep = ",")
colnames(countData) <- gsub("\\.", "-", colnames(countData))
countData$KEGG_ko <- gsub("ko:", "", countData$KEGG_ko)

## GO DATA
GO_data <- read.csv("Outputs/GO_data.csv")
GO_data$gene_id <- gsub("ko:", "", GO_data$gene_id)

## GENE LENGTH
length_df <- read.csv("Outputs/length_df.csv")
length_df$KEGG_ko <- gsub("ko:", "", length_df$KEGG_ko)

gene_lenth <- length_df$Length
names(gene_lenth) <- length_df$KEGG_ko

## LOAD DESeq2 OBJECT
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = sample_df, 
                              design = ~Genus, tidy = TRUE)

#
#### DESeq2 ####
dds <- DESeq(dds)

#
#### GLOBAL PCA ####
vst <- vst(dds, blind = TRUE)

pcaData <- plotPCA(vst, intgroup = "Genus", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"), 2)
ggplot(pcaData, aes(PC1, PC2, color = Genus)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar[2],"% variance")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = c("#1b9e77", "#caf55d", "#cc3939", "#cab2d6", "#8da0cb")) +
  labs(color = "Dominant Genus")

#
#### GLOBAL CHANGE IN GENE EXPRESSION QUANTIFICATION - vs. SACCHAROMYCES ####

## Lachancea
res_Lt <- results(dds, contrast = c("Genus", "Lachancea", "Saccharomyces"),
                  pAdjustMethod = "fdr", alpha = 0.05, independentFiltering = F)

res_Lt <- res_Lt[complete.cases(res_Lt),]
res_LtOrdered <- res_Lt[order(res_Lt$padj),]

ress_Lt <- as.data.frame(res_LtOrdered)
ress_Lt$Locus.Tag <- row.names(ress_Lt)

DEG_Lt <- as.integer(res_Lt$padj < 0.05 & abs(res_Lt$log2FoldChange) >= 1)
names(DEG_Lt) <- row.names(res_Lt)
table(DEG_Lt)


## Hanseniaspora
res_Hs <- results(dds, contrast = c("Genus", "Hanseniaspora", "Saccharomyces"),
                  pAdjustMethod = "fdr", alpha = 0.05, independentFiltering = F)

res_Hs <- res_Hs[complete.cases(res_Hs),]
res_HsOrdered <- res_Hs[order(res_Hs$padj),]

ress_Hs <- as.data.frame(res_HsOrdered)
ress_Hs$Locus.Tag <- row.names(ress_Hs)

DEG_Hs <- as.integer(res_Hs$padj < 0.05 & abs(res_Hs$log2FoldChange) >= 1)
names(DEG_Hs) <- row.names(res_Hs)
table(DEG_Hs)

#
#### HISTOGRAM ####
DEG_hist <- rbind.data.frame(cbind(ress_Lt[ress_Lt$Locus.Tag %in% names(DEG_Lt[DEG_Lt == 1]),],
                                   Comparison = "Lachancea"),
                             cbind(ress_Hs[ress_Hs$Locus.Tag %in% names(DEG_Hs[DEG_Hs == 1]),],
                                   Comparison = "Hanseniaspora"))


DEG_hist$Comparison <- factor(DEG_hist$Comparison, levels = c("Lachancea", "Hanseniaspora"))

ggplot(DEG_hist, aes(x = abs(log2FoldChange), fill = Comparison)) + 
  geom_histogram(binwidth = 1, position = "dodge", color = "gray30") +
  theme_minimal() + xlab("|LFC|") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  theme(axis.text.y = element_text(size = 18, color  = "black"),
        axis.title.x = element_text(size = 18, color  = "black"),
        axis.title.y = element_text(size = 18, color  = "black"),
        legend.text = element_text(size = 18, color  = "black"),
        legend.title = element_text(size = 18, color  = "black"),
        axis.text.x = element_text(size = 18, color  = "black")) 

aggregate(abs(DEG_hist$log2FoldChange), list(DEG_hist$Comparison), sum)

#
#### VENN DIAGRAM ####
ress_Lt.venn <- ress_Lt
ress_Lt.venn$Lachancea <- 0
ress_Lt.venn$Lachancea[abs(ress_Lt.venn$log2FoldChange) > 1 & ress_Lt.venn$padj <= 0.05] <- 1

ress_Hs.venn <- ress_Hs
ress_Hs.venn$Hanseniaspora <- 0
ress_Hs.venn$Hanseniaspora[abs(ress_Hs.venn$log2FoldChange) > 1 & ress_Hs.venn$padj <= 0.05] <- 1

venn_plot <- merge(ress_Lt.venn, ress_Hs.venn, by = "row.names")
comp <- cbind(Lachancea = (venn_plot$Lachancea == 1), Hanseniaspora = (venn_plot$Hanseniaspora == 1))

vdc <- vennCounts(comp)
class(vdc) <- "matrix"
df.vdc <- as.data.frame(vdc) %>% 
  mutate(x = c(2, 1.4, -1.4, 0), y = c(-1.5, 0, 0, 0))

df.venn <- data.frame(x = c(-0.75, 0.75), y = c(0, 0), labels = c("Lachancea", "Hanseniaspora"))
df.venn$labels <- factor(df.venn$labels, levels = c("Lachancea", "Hanseniaspora"))

ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, size = 1, colour = "gray30") +
  coord_fixed() + theme_void() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 18, color  = "black")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = NULL) + 
  geom_text(data = df.vdc, aes(x = x, y = y, label = Counts), size = 7)


#
#### BIOLOGICAL ENRICHMENT - KEGG ####
## Lachancea
enrich_Lt <- enrichKEGG(names(DEG_Lt[DEG_Lt == 1]), organism = "ko", keyType = "kegg")

LFC_Lt <- res_Lt$log2FoldChange
names(LFC_Lt) <- row.names(res_Lt)


## Hanseniaspora
enrich_Hs <- enrichKEGG(names(DEG_Hs[DEG_Hs == 1]), organism = "ko", keyType = "kegg")

LFC_Hs <- res_Hs$log2FoldChange
names(LFC_Hs) <- row.names(res_Hs)

## PLOT

kk.plot <- rbind(data.frame(enrich_Lt, Comparison = "Lachancea"),
                 data.frame(enrich_Hs, Comparison = "Hanseniaspora"))

to.rm <- c("ko05016", "ko05014", "ko05012", "ko05020", "ko00710", "ko05022", "ko05017", "ko05010",
           "ko04932", "ko04260", "ko01524", "ko05415", "ko04114", "ko04213", "ko05171")

kk.plot <- kk.plot[!kk.plot$ID %in% to.rm,]

kk.plot$Comparison <- factor(kk.plot$Comparison, levels = c("Lachancea", "Hanseniaspora"))
kk.plot <- kk.plot[order(kk.plot$Comparison, kk.plot$p.adjust),]
kk.plot$Description <- factor(kk.plot$Description, levels = unique(kk.plot$Description))

ggplot(kk.plot, aes(x = Description, y = -log(p.adjust), fill = Comparison)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = "Comparison") + xlab("KEGG pathways") + ylab("-log adjusted p value") +
  coord_flip()



# Pathviews
library(pathview)

pathview(gene.data = LFC.Sc_Lt, 
         pathway.id = "ko00620", 
         species = "ko",  
         gene.idtype = "kegg")




#
#### BIOLOGICAL ENRICHMENT - Gene Ontology (BP) ####
## Lachancea
gene_length.Lt <- gene_lenth[names(gene_lenth) %in% names(DEG_Lt)]

pwf_Lt <- nullp(DEG_Lt, bias.data = gene_length.Lt)

GO_Lt <- goseq(pwf_Lt, gene2cat = GO_data)

GO_Lt.BP <- subset(GO_Lt, ontology == "BP")
GO_Lt.BP <- head(GO_Lt.BP, n = 15)


## Hanseniaspora
gene_length.Hs <- gene_lenth[names(gene_lenth) %in% names(DEG_Hs)]

pwf_Hs <- nullp(DEG_Hs, bias.data = gene_length.Hs)

GO_Hs <- goseq(pwf_Hs, gene2cat = GO_data)

GO_Hs.BP <- subset(GO_Hs, ontology == "BP")
GO_Hs.BP <- head(GO_Hs.BP, n = 15)


## PLOT
GO_plot.BP <- rbind.data.frame(cbind(GO_Lt.BP, Comparison = "Lachancea"),
                               cbind(GO_Hs.BP, Comparison = "Hanseniaspora"))

GO_plot.BP$Comparison <- factor(GO_plot.BP$Comparison, levels = c("Lachancea", "Hanseniaspora"))

ggplot(GO_plot.BP) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Comparison), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  coord_flip() +
  theme_bw()+
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("Gene count")  + xlab("")


#
############################################################################ DIVIDE BY DOMINANT SPECIES ####
################################################### SACCHAROMYCES ####
sample_sacc <- subset(sample_df, Genus == "Saccharomyces")

dds_sacc <- DESeqDataSetFromMatrix(countData = cbind(KEGG_ko = countData[,1],
                                                     countData[,colnames(countData) %in% sample_sacc$Sample_ID]), 
                              colData = sample_sacc, 
                              design = ~Condition, tidy = TRUE)

#### DESeq2 ####
dds_sacc <- DESeq(dds_sacc)

#
#### GLOBAL PCA ####
rlog_sacc <- rlog(dds_sacc, blind = TRUE)

pcaData_sacc <- plotPCA(rlog_sacc, intgroup = "Condition", returnData = TRUE)
percentVar_sacc <- round(100 * attr(pcaData_sacc, "percentVar"), 2)
ggplot(pcaData_sacc, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_sacc[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_sacc[2],"% variance")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = c("#1b9e77", "#caf55d", "#cc3939", "#cab2d6", "#8da0cb")) +
  labs(color = "Dominant Genus")

#
################################################### LACHANCEA ####
sample_lach <- subset(sample_df, Genus == "Lachancea")

dds_lach <- DESeqDataSetFromMatrix(countData = cbind(KEGG_ko = countData[,1],
                                                     countData[,colnames(countData) %in% sample_lach$Sample_ID]), 
                                   colData = sample_lach, 
                                   design = ~Condition, tidy = TRUE)

#### DESeq2 ####
dds_lach <- DESeq(dds_lach)

#
#### GLOBAL PCA ####
rlog_lach <- rlog(dds_lach, blind = TRUE)

pcaData_lach <- plotPCA(rlog_lach, intgroup = "Farming", returnData = TRUE)
percentVar_lach <- round(100 * attr(pcaData_lach, "percentVar"), 2)
ggplot(pcaData_lach, aes(PC1, PC2, color = Farming)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_lach[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_lach[2],"% variance")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = c("#1b9e77", "#caf55d", "#cc3939", "#cab2d6", "#8da0cb", "gray")) +
  labs(color = "Dominant Genus")

#
################################################### HANSENIASPORA ####
sample_hans <- subset(sample_df, Genus == "Hanseniaspora")

dds_hans <- DESeqDataSetFromMatrix(countData = cbind(KEGG_ko = countData[,1],
                                                     countData[,colnames(countData) %in% sample_hans$Sample_ID]), 
                                   colData = sample_hans, 
                                   design = ~Condition, tidy = TRUE)

#### DESeq2 ####
dds_hans <- DESeq(dds_hans)

#
#### GLOBAL PCA ####
rlog_hans <- rlog(dds_hans, blind = TRUE)

pcaData_hans <- plotPCA(rlog_hans, intgroup = "Condition", returnData = TRUE)
percentVar_hans <- round(100 * attr(pcaData_hans, "percentVar"), 2)
ggplot(pcaData_hans, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_hans[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_hans[2],"% variance")) + 
  theme_bw() +
  theme(axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = c("#1b9e77", "#caf55d", "#cc3939", "#cab2d6", "#8da0cb", "gray")) +
  labs(color = "Dominant Genus")

#

#
#### SAVE DATA ####
save.image("DESeq2_meta.Rdata")
