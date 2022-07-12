# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - ITS sequence analysis of T0-GM samples

# Set the project location as working directory
setwd("C:/Users/Laboratorio14/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/1_Reg-Farm//")

#
#### LIBRARIES ####
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(vegan)
library(hillR)
library(phyloseq)
library(DESeq2)
library(agricolae)
library(raster)
library(rnaturalearth)
library(sf)
library(ggnewscale)

load("com_analysis.RData")

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
asv_GM <- readRDS("Outputs/ASV_t0-GM.rds")
asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x))

tax_GM <- readRDS("Outputs/tax_t0-GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- gsub("^[a-z]__", "", as.matrix(tax_GM))

tax_GM.un <- set_unid(tax_GM)
tax_GM[is.na(tax_GM)] <- "Unidentified"

## SAMPLE DATA
sample_df <- read.table("Inputs/NGM_initial.txt", sep = "\t", header = TRUE)
row.names(sample_df) <- sample_df$Seq_ID
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Origin <- factor(sample_df$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", 
                                                        "R3B", "R3C"))
sample_df$Region <- factor(sample_df$Region, levels = c("Ribera del Guadiana", "Valdepeñas", "La Mancha",
                                                        "Madrid", "Rioja"))


#
#### SAMPLING MAP (Temperature) ####

# Load WorldClim MAT data and crop the Spanish area
climate1 <- getData("worldclim", var = "bio", res = 0.5, lon = -3.703790, 
                   lat = 40.416775)
climate2 <- getData("worldclim", var = "bio", res = 0.5, lon = 3.703790, 
                   lat = 40.416775)
climate <- mosaic(climate1, climate2)
climate_crop <- crop(climate, extent(-10, 5, 35, 44))

tmean_spain_df <- as.data.frame(climate_crop, xy = TRUE, na.rm = TRUE)

# Load the Spanish community boundaries
spain <- getData("GADM", country = "ESP", level = 1)

sp.lim <- data.frame(ylim = c(35, 44), xlim = c(-10, 5))
sp.gadm <- st_as_sf(spain)
sp.bbox <- st_bbox(c(xmin = sp.lim$xlim[1],
                     xmax = sp.lim$xlim[2],
                     ymin = sp.lim$ylim[1],
                     ymax = sp.lim$ylim[2]))
sp.gadm <- st_crop(sp.gadm, sp.bbox)

# Draw the map
ggplot() +
  geom_raster(data = tmean_spain_df, aes(x = x, y = y, fill = layer.1/10)) +
  scale_fill_gradient2(name = "Temperature (°C)", low = "blue", mid = "white", high = "red", midpoint = 10,
                       limits = c(-4, 19.9)) +
  geom_sf(data = sp.gadm, fill = NA, color = "black") +
  coord_sf(ylim = sp.lim$ylim, xlim = sp.lim$xlim) +
  new_scale("fill") +
  geom_point(data = sample_df, aes(x = Longitude, y = Latitude, fill = Origin),
             colour = "black", shape = 21, size = 2.5, stroke = 1.5) +
  scale_fill_manual(values = c("#2ac219", "#1949c2", "#dba54d", "#e02424", 
                                "#c124e0", "#89209e", "#a6165c", "#750f41", "#5c105e")) +
  theme_void()

#
#### GEOGRAPHIC DISTANCES ####
coord_df <- sample_df[,c(1,9:8)]
row.names(coord_df) <- coord_df[,1]
coord_df <- coord_df[,-1]

coord_dist <- data.frame(matrix(nrow = nrow(coord_df), ncol = nrow(coord_df)), row.names = row.names(coord_df))
colnames(coord_dist) <- row.names(coord_dist)

for (row in 1:nrow(coord_dist)) {
  for (col in 1:ncol(coord_dist)) {
    
    coord_dist[row,col] <- round(distGeo(c(coord_df$Longitude[row], coord_df$Latitude[row]), 
                                         c(coord_df$Longitude[col], coord_df$Latitude[col]))/1000, 2)
  }
  
}

#
#### TAXONOMIC EXPLORATION ####
## Genus

asv.t_GM.p <- melt(asv.t_GM)
colnames(asv.t_GM.p) <- c("Id", "Seq_ID", "value")
asv.t_GM.p <- merge(asv.t_GM.p, tax_GM[,c(6,8)])
asv.t_GM.p <- merge(asv.t_GM.p, sample_df, by = "Seq_ID")

asv.t_plot <- aggregate(asv.t_GM.p$value, list(asv.t_GM.p$Sample_ID, asv.t_GM.p$Genus), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "value")
asv.t_plot$Genus[asv.t_plot$value < 0.05] <- "Other"

asv.t_plot <- aggregate(asv.t_plot$value, list(asv.t_plot$Sample_ID, asv.t_plot$Genus), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "value")

orderG <- levels(factor(asv.t_plot$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

asv.t_plot <- merge(asv.t_plot, sample_df, by = "Sample_ID")
asv.t_plot$plot_ID <- paste(asv.t_plot$Management, asv.t_plot$Condition, sep = " ")
asv.t_plot$Origin <- factor(asv.t_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))

ggplot(asv.t_plot, 
       aes(x = plot_ID, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = c("#f26b9f", "#ffff33", "#37876a", "#f58d5d",
                                               "#b33cb5", "#9e66d1", "#cab2d6", "#8da0cb", "#aec219",
                                               "#b15928", "#1b9e77", "#d4556a", "#6a3d9a", "#bf5b17")) +
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
  guides(fill = guide_legend(nrow = 3)) +
  facet_wrap(~Origin, nrow = 1)

#
#### ALPHA DIVERSITY ANALYSIS ####

alpha_GM <- cbind.data.frame(
  # Hill based taxonomic alpha diversity
  t.q0 = hill_taxa(t(asv.t_GM), q = 0),
  t.q1 = hill_taxa(t(asv.t_GM), q = 1),
  t.q2 = hill_taxa(t(asv.t_GM), q = 2))

alpha_GM$Seq_ID <- row.names(alpha_GM)
alpha_GM <- merge(alpha_GM, sample_df, by = "Seq_ID")

alpha_GM.plot <- melt(alpha_GM)

## Origin

summary(aov(t.q1~Condition, alpha_GM))
LSD.test(aov(t.q0~Origin, alpha_GM), "Origin")$group
LSD.test(aov(t.q1~Origin, alpha_GM), "Origin")$group
LSD.test(aov(t.q2~Origin, alpha_GM), "Origin")$group

ggplot(data = alpha_GM.plot, aes(x = Origin, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y", nrow = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"),
        strip.text = element_text(size = 17, color = "black"))



## Management
ggplot(data = alpha_GM.plot, aes(x = Management, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"),
        strip.text = element_text(size = 17, color = "black"))


t.test(subset(alpha_GM.plot, variable == "t.q0" & Management == "ECO")$value,
       subset(alpha_GM.plot, variable == "t.q0" & Management == "CONV")$value)

t.test(subset(alpha_GM.plot, variable == "t.q1" & Management == "ECO")$value,
       subset(alpha_GM.plot, variable == "t.q1" & Management == "CONV")$value)

t.test(subset(alpha_GM.plot, variable == "t.q2" & Management == "ECO")$value,
       subset(alpha_GM.plot, variable == "t.q2" & Management == "CONV")$value)

#
#### BETA DIVERSITY ANALYSIS ####
bray_GM <- vegdist(t(asv.t_GM), method = "bray")

set.seed(1)
nMDS_GM <- metaMDS(bray_GM)
nMDS_GM$stress

nMDS_GM.plot <- as.data.frame(nMDS_GM[["points"]])
nMDS_GM.plot$Seq_ID <- rownames(nMDS_GM.plot)
nMDS_GM.plot <- merge(nMDS_GM.plot, sample_df, by = "Seq_ID")

ggplot(nMDS_GM.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, color = Origin, shape = Management), size = 4) +
  scale_color_manual(values = c("#2ac219", "#1949c2", "#dba54d", "#e02424", 
                                "#c124e0", "#89209e", "#a6165c", "#750f41", "#5c105e")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

set.seed(1)
adonis2(bray_GM ~ Management, data = sample_df)

sample_df$group <- paste(sample_df$Origin, sample_df$Management)

beta_GM <- betadisper(bray_GM, sample_df$Origin, type = "centroid")
plot(beta_GM)

bdisp.df_GM <- cbind.data.frame(Distance_to_centroid = beta_GM$distances, group = beta_GM$group)
bdisp.df_GM$Seq_ID <- row.names(bdisp.df_GM)
bdisp.df_GM <- merge(bdisp.df_GM, sample_df[,-7], by = "Seq_ID")

ggplot(data = bdisp.df_GM, aes(x = Management, y = Distance_to_centroid, color = Management)) +
  geom_boxplot()


bray.df_GM <- melt(as.matrix(bray_GM))
colnames(bray.df_GM)[1] <- "Seq_ID" 

bray.df_GM <- merge(bray.df_GM, sample_df, by = "Seq_ID")
colnames(bray.df_GM)[1:2] <- c("Seq_ID.x", "Seq_ID") 

bray.df_GM <- merge(bray.df_GM, sample_df[,6:7], by = "Seq_ID")

bray.df_GM <- subset(bray.df_GM, group.x == group.y & value != 0)

ggplot(data = bray.df_GM, aes(x = Management, y = value, color = Management)) +
  geom_boxplot()

t.test(subset(bray.df_GM, Management == "ECO")$value,
       subset(bray.df_GM, Management == "CONV")$value)

set.seed(1)
adonis2(bray_GM ~ group, data = sample_df)


#
#### CORRELATION BETWEEN GEOGRAPHIC DISTANCE AND ENVIRONTMENT AND COMMUNITY ####
env_dist <- as.matrix(dist(apply(sample_df[,12:17], 2, function(x) (x-min(x))/(max(x)-min(x)))))

env_geo <- mantel(env_dist, coord_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
env_geo

env_geo.plot <- merge(melt(as.matrix(env_dist)), melt(as.matrix(coord_dist)), by = c("Var1", "Var2"))
colnames(env_geo.plot) <- c("Var1", "Var2", "env_dist", "coord_dist")

ggplot(env_geo.plot) +
  geom_point(aes(x = coord_dist, y = env_dist)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 24, color = "black"),
        axis.text.x = element_text(size = 22, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        axis.text.y = element_text(size = 22, color = "black")) +
  xlab("Geographic Distance") + ylab("Environmental Distance") +
  annotate("text", label = "p = 0.571, p-value < 0.001", x = 400, y = 2, size = 7)

com_geo <- mantel(bray_GM, coord_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
com_geo

com_geo.plot <- merge(melt(as.matrix(bray_GM)), melt(as.matrix(coord_dist)), by = c("Var1", "Var2"))
colnames(com_geo.plot) <- c("Var1", "Var2", "com_dist", "coord_dist")

ggplot(com_geo.plot) +
  geom_point(aes(x = coord_dist, y = com_dist)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 24, color = "black"),
        axis.text.x = element_text(size = 22, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        axis.text.y = element_text(size = 22, color = "black")) +
  xlab("Geographic Distance") + ylab("Community Distance") +
  annotate("text", label = "p = -0.034, p-value = 0.761", x = 400, y = 0.8, size = 7)

#### DIFERENTIAL ABUNDANCE ANALYSIS ####
asv.t_GM.p <- melt(asv.t_GM)
colnames(asv.t_GM.p) <- c("Id", "Seq_ID", "value")
asv.t_GM.p <- merge(asv.t_GM.p, tax_GM.un[,c(6,8)])
asv.t_GM.p <- merge(asv.t_GM.p, sample_df, by = "Seq_ID")

asv.t_plot.0 <- aggregate(asv.t_GM.p$value, list(asv.t_GM.p$Sample_ID, asv.t_GM.p$Genus), sum)
colnames(asv.t_plot.0) <- c("Sample_ID", "Genus", "value")
asv.t_plot.0 <- merge(asv.t_plot.0, sample_df, by = "Sample_ID")

phy_GM <- phyloseq(otu_table(asv_GM, taxa_are_rows = FALSE), tax_table(tax_GM.un), sample_data(sample_df))
phy_GM.genus <- tax_glom(phy_GM, taxrank = "Genus")

## Management ~ Global

dds_GM_man.tot <- phyloseq_to_deseq2(phy_GM.genus, ~ Management)

keep <- rowSums(counts(dds_GM_man.tot)) >= 1000
dds_GM_man.tot <- dds_GM_man.tot[keep,]

dds_GM_man.tot <- DESeq(dds_GM_man.tot, test = "Wald", fitType = "parametric")

res_GM_man.tot <- results(dds_GM_man.tot, cooksCutoff = FALSE)

sigtab_GM_man.tot <- as.data.frame(res_GM_man.tot[which(res_GM_man.tot$padj < 0.05), ])
sigtab_GM_man.tot <- merge(sigtab_GM_man.tot, tax_GM.un[,6], by = "row.names")
colnames(sigtab_GM_man.tot)[c(1,8)] <- c("Id", "Genus")

x <- tapply(sigtab_GM_man.tot$log2FoldChange, sigtab_GM_man.tot$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_GM_man.tot$Genus = factor(as.character(sigtab_GM_man.tot$Genus), levels=names(x))

ggplot(sigtab_GM_man.tot, aes(x = Genus, y = log2FoldChange)) + 
  geom_point(size = 4) + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"))

ggplot(asv.t_plot.0[asv.t_plot.0$Genus %in% sigtab_GM_man.tot$Genus,], 
       aes(x = Genus, y = value, color = Management)) + 
  geom_boxplot(size = 1) + 
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 15, color = "black")) +
  xlab("Sample") + ylab("Abundance") +
  guides(fill = guide_legend(nrow = 3))



#
#### MUST ANALYSIS ####
must_df <- sample_df[,c(2:6, 12:17)]

must_pca <- prcomp(must_df[,c(6:7,9:11)], scale = TRUE)

must_pca.plot <- must_pca$x
must_pca.plot <- merge(must_df, must_pca.plot, by = "row.names")

ggplot(must_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming), size = 4) +
  scale_color_manual(values = c("#2ac219", "#1949c2", "#dba54d", "#e02424", 
                                "#c124e0", "#89209e", "#a6165c", "#750f41", "#5c105e")) +
  xlab(paste("PC1: ", round(((must_pca$sdev)^2 / sum((must_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((must_pca$sdev)^2 / sum((must_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        legend.text = element_text(size = 24, color = "black"),
        legend.title = element_text(size = 22, color = "black"),
        axis.text.x = element_text(size = 22, color = "black"))

ggplot(must_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Condition, shape = Farming), size = 4) +
  scale_color_manual(values = c("#2ac219", "#1949c2", "#dba54d", "#e02424", 
                                "#c124e0", "#89209e", "#a6165c", "#750f41", "#5c105e")) +
  xlab(paste("PC1: ", round(((must_pca$sdev)^2 / sum((must_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((must_pca$sdev)^2 / sum((must_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        legend.text = element_text(size = 24, color = "black"),
        legend.title = element_text(size = 22, color = "black"),
        axis.text.x = element_text(size = 22, color = "black"))

must_pca

#
#### CONSTRAINED ANALYSIS ####

adonis2(as.matrix(bray_GM) ~ ., sample_df[,c(12:17)], permutations = 1000)
cap_GM <- capscale(as.matrix(bray_GM) ~ ., sample_df[,c(12:17)], na.action = na.omit)

anova(cap_GM)
RsquareAdj(cap_GM)

arrowmat <- rbind.data.frame(scores(cap_GM, display = "bp"),
                             scores(cap_GM, display = "cn"))

arrowmat$Sample <- row.names(arrowmat)

sitemat <- as.data.frame(scores(cap_GM, display = "sites"))
sitemat$Seq_ID <- row.names(sitemat)
sitemat <- merge(sitemat, sample_df, by = "Seq_ID")

summary(cap_GM)$cont

ggplot() + 
  geom_point(data = sitemat, aes(x = CAP1, y = CAP2, color = Origin), size = 3) +
  geom_segment(data = arrowmat, aes(x = 0, y = 0, xend = CAP1*2, yend = CAP2*2), 
               arrow = arrow(length = unit(0.2,"cm")), size = 1, alpha = 0.8) + 
  geom_text(data = arrowmat, aes(x = CAP1*2.5, y = CAP2*2.5, label = Sample), size = 6) +
  scale_color_manual(values = c("#2ac219", "#1949c2", "#dba54d", "#e02424", 
                                "#c124e0", "#89209e", "#a6165c", "#750f41", "#5c105e")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        legend.text = element_text(size = 22, color = "black"),
        legend.title = element_text(size = 22, color = "black"),
        axis.text.x = element_text(size = 24, color = "black")) +
  xlab("CAP1 (7.91%)") + ylab("CAP2 (2.92%)")










#
#### SAVE OUTPUTS ####
save.image("com_analysis.RData")

#