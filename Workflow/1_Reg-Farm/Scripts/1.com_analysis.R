# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - ITS sequence analysis of T0-GM samples

# Set the project location as working directory
setwd("C:/Users/Migueldc/OneDrive - Universidad Complutense de Madrid (UCM)/Wineteractions/GitHub/Wineteractions/1_Reg-Farm//")

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
library(geosphere)

load("com_analysis.RData")

#
#### FUNCTIONS ####
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

r2Log <- function(model) {
  
  summaryLog <- summary(model)
  1 - summaryLog$deviance / summaryLog$null.deviance
  
}

#
#### DATA LOADING ####

## SAMPLE DATA
sample_df <- read.table("../0_ITS-GM/Inputs/sample_GM.txt", sep = "\t", header = TRUE)
row.names(sample_df) <- sample_df$Seq_ID
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Origin <- factor(sample_df$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", 
                                                        "R3B", "R3C"))
sample_df$Region <- factor(sample_df$Region, levels = c("Ribera del Guadiana", "Valdepeñas", "La Mancha",
                                                        "Madrid", "Rioja"))

sample_df <- cbind.data.frame(Sample_ID = paste(sample_df$Origin, sample_df$Farming, sample_df$Condition, sep = "-"),
                              sample_df)

sample_df <- subset(sample_df, Stage == "0_initial")

## Community data
asv_GM <- readRDS("../0_ITS-GM/Outputs/ASV_GM.rds")
asv_GM <- asv_GM[row.names(asv_GM) %in% sample_df$Seq_ID, ]
asv_GM <- asv_GM[, colSums(asv_GM != 0) > 0]

asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x))

tax_GM <- readRDS("../0_ITS-GM/Outputs/tax_GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- gsub("^[a-z]__", "", as.matrix(tax_GM))
tax_GM <- tax_GM[row.names(tax_GM) %in% colnames(asv_GM), ]

tax_GM.un <- set_unid(tax_GM)
tax_GM[is.na(tax_GM)] <- "Unidentified"

#
#### SAMPLING MAP (Temperature) ####

# Load WorldClim MAT data and crop the Spanish area
climate1 <- getData("worldclim", var = "bio", res = 0.5, lon = -3.703790, 
                   lat = 40.416775)
climate2 <- getData("worldclim", var = "bio", res = 0.5, lon = 3.703790, 
                   lat = 40.416775)
climate <- mosaic(climate1, climate2, fun = mean)
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
  scale_fill_distiller(name = "Temperature (°C)", palette = "RdYlBu",
                       limits = c(-4, 19.9),
                       guide = guide_colourbar(frame.colour = "black", ticks.colour = "black")) +
  geom_sf(data = sp.gadm, fill = NA, color = "black", size = 1/2) +
  coord_sf(ylim = sp.lim$ylim, xlim = sp.lim$xlim) +
  new_scale("fill") +
  geom_point(data = sample_df, aes(x = Longitude, y = Latitude, fill = Origin), shape = 21, size = 3.5) +
  scale_fill_manual(values = c("#2ac219", "#19b7c2", "#dba54d", "#ae0e36", 
                                "#433254", "#5900ff", "#7802a3", "#8849d1", "#9d0dd1")) +
  theme_void() +
  theme(legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"))

#
#### SAMPLING MAP (Elevation) ####

# Load WorldClim MAT data and crop the Spanish area
altitude1 <- getData("worldclim", var = "alt", res = 0.5, lon = -3.703790, 
                    lat = 40.416775)
altitude2 <- getData("worldclim", var = "alt", res = 0.5, lon = 3.703790, 
                    lat = 40.416775)
altitude <- mosaic(altitude1, altitude2, fun = mean)
altitude_crop <- crop(altitude, extent(-10, 5, 35, 44))

alt_spain_df <- as.data.frame(altitude_crop, xy = TRUE, na.rm = TRUE)

# Draw the map
ggplot() +
  geom_raster(data = alt_spain_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_distiller(name = "Altitude (m)", palette = "BrBG",
                       guide = "none") +
  geom_sf(data = sp.gadm, fill = NA, color = "black", size = 1/2) +
  coord_sf(ylim = sp.lim$ylim, xlim = sp.lim$xlim) +
  new_scale("fill") +
  geom_point(data = sample_df, aes(x = Longitude, y = Latitude, fill = Origin), shape = 21, size = 3.5) +
  scale_fill_manual(values = c("#2ac219", "#19b7c2", "#dba54d", "#ae0e36", 
                               "#433254", "#5900ff", "#7802a3", "#8849d1", "#9d0dd1")) +
  theme_void() +
  theme(legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"))

#
#### GEOGRAPHIC DISTANCES ####
coord_df <- sample_df[,c(2,9:8)]
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
asv.t_GM.p <- melt(asv.t_GM)
colnames(asv.t_GM.p) <- c("Id", "Seq_ID", "value")
asv.t_GM.p <- merge(asv.t_GM.p, tax_GM[,c(6,8)], by = "Id")

asv.t_GM.p <- merge(asv.t_GM.p, sample_df[c(1:7)], by = "Seq_ID")

## Genus

asv.t_plot <- aggregate(asv.t_GM.p$value, list(asv.t_GM.p$Sample_ID, asv.t_GM.p$Genus), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "value")
asv.t_plot$Genus[asv.t_plot$value < 0.05] <- "Other"

asv.t_plot <- aggregate(asv.t_plot$value, list(asv.t_plot$Sample_ID, asv.t_plot$Genus), sum)
colnames(asv.t_plot) <- c("Sample_ID", "Genus", "value")

orderG <- levels(factor(asv.t_plot$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

asv.t_plot <- merge(asv.t_plot, sample_df[,1:5], by = "Sample_ID")
asv.t_plot$Origin <- factor(asv.t_plot$Origin, 
                            levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))

asv.t_plot$Sample_name <- paste(asv.t_plot$Farming, asv.t_plot$Condition, sep = "-")
asv.t_plot$Sample_name <- factor(asv.t_plot$Sample_name, 
                                 levels = c("CONV-Control", "CONV-18C", "CONV-NH4", "CONV-SO2",
                                            "ECO-Control", "ECO-18C", "ECO-NH4", "ECO-SO2"))

ggplot(asv.t_plot, 
       aes(x = Sample_name, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Genus", values = c("#ccebc5", "#ffff33", "#e6ab02", "#80b1d3",
                                               "#e6f5c9", "#ff7f00", "#fff2ae",
                                               "#8da0cb", "#b15928", "#bebada", "#1b9e77",
                                               "#fbb4ae", "#6a3d9a", "#bf5b17")) +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        strip.text.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black")) +
  xlab("Sample") + ylab("Abundance") +
  guides(fill = guide_legend(nrow = 3)) + facet_wrap(~Origin, nrow = 1)


#
#### ALPHA DIVERSITY ANALYSIS ####

alpha_GM <- cbind.data.frame(
  # Hill based taxonomic alpha diversity
  t.q0 = hill_taxa(t(asv.t_GM), q = 0),
  t.q1 = hill_taxa(t(asv.t_GM), q = 1),
  t.q2 = hill_taxa(t(asv.t_GM), q = 2))

alpha_GM$Seq_ID <- row.names(alpha_GM)
alpha_GM <- merge(alpha_GM, sample_df, by = "Seq_ID")

alpha_GM.plot <- melt(alpha_GM[,1:9])

## Origin

summary(aov(t.q1~Origin, alpha_GM))
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
ggplot(data = alpha_GM.plot, aes(x = Farming, y = value)) +
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


t.test(subset(alpha_GM.plot, variable == "t.q0" & Farming == "ECO")$value,
       subset(alpha_GM.plot, variable == "t.q0" & Farming == "CONV")$value)

t.test(subset(alpha_GM.plot, variable == "t.q1" & Farming == "ECO")$value,
       subset(alpha_GM.plot, variable == "t.q1" & Farming == "CONV")$value)

t.test(subset(alpha_GM.plot, variable == "t.q2" & Farming == "ECO")$value,
       subset(alpha_GM.plot, variable == "t.q2" & Farming == "CONV")$value)

#
#### ALPHA DIVERSITY: ORIGIN vs. FARMING ####

###q = 0
null_t.0 <- glm(log(t.q0) ~ 1,data = alpha_GM)
summary(null_t.0)

fit_t.0_or <- glm(log(t.q0) ~ Origin, data = alpha_GM)
summary(fit_t.0_or)

fit_t.0_fr <- glm(log(t.q0) ~ Farming, data = alpha_GM)
summary(fit_t.0_fr)

r2Log(null_t.0)
r2Log(fit_t.0_or)
r2Log(fit_t.0_fr)

###q = 1
null_t.1 <- glm(log(t.q1) ~ 1,data = alpha_GM)
summary(null_t.1)

fit_t.1_or <- glm(log(t.q1) ~ Origin, data = alpha_GM)
summary(fit_t.1_or)

fit_t.1_fr <- glm(log(t.q1) ~ Farming, data = alpha_GM)
summary(fit_t.1_fr)

r2Log(null_t.1)
r2Log(fit_t.1_or)
r2Log(fit_t.1_fr)

###q = 2
null_t.2 <- glm(log(t.q2) ~ 1,data = alpha_GM)
summary(null_t.2)

fit_t.2_or <- glm(log(t.q2) ~ Origin, data = alpha_GM)
summary(fit_t.2_or)

fit_t.2_fr <- glm(log(t.q2) ~ Farming, data = alpha_GM)
summary(fit_t.2_fr)

r2Log(null_t.2)
r2Log(fit_t.2_or)
r2Log(fit_t.2_fr)

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
  geom_point(aes(x = MDS1, y = MDS2, color = Origin, shape = Farming), size = 4) +
  scale_color_manual(values = c("#2ac219", "#19b7c2", "#dba54d", "#ae0e36", 
                                "#433254", "#5900ff", "#7802a3", "#8849d1", "#9d0dd1")) +
  annotate(geom = "text", x = 0.4, y = 0.35, label = round(nMDS_GM$stress, 3), size = 7) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 22, color = "black"),
        legend.title = element_text(size = 22, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"))

set.seed(1)
adonis2(bray_GM ~ Origin, data = sample_df)

sample_df$group <- paste(sample_df$Origin, sample_df$Farming)

beta_GM <- betadisper(bray_GM, sample_df$Origin, type = "centroid")
plot(beta_GM)

bdisp.df_GM <- cbind.data.frame(Distance_to_centroid = beta_GM$distances, group = beta_GM$group)
bdisp.df_GM$Seq_ID <- row.names(bdisp.df_GM)
bdisp.df_GM <- merge(bdisp.df_GM, sample_df[,-7], by = "Seq_ID")

ggplot(data = bdisp.df_GM, aes(x = Farming, y = Distance_to_centroid, color = Farming)) +
  geom_boxplot()


bray.df_GM <- melt(as.matrix(bray_GM))
colnames(bray.df_GM)[1] <- "Seq_ID" 

bray.df_GM <- merge(bray.df_GM, sample_df, by = "Seq_ID")
colnames(bray.df_GM)[1:2] <- c("Seq_ID.x", "Seq_ID") 

bray.df_GM <- merge(bray.df_GM, sample_df[,c(1:6,33)], by = "Seq_ID")

bray.df_GM <- subset(bray.df_GM, group.x == group.y & value != 0)

ggplot(data = bray.df_GM, aes(x = Farming.x, y = value, color = Farming.x)) +
  geom_boxplot()

t.test(subset(bray.df_GM, Farming.x == "ECO")$value,
       subset(bray.df_GM, Farming.x == "CONV")$value)

set.seed(1)
adonis2(bray_GM ~ group, data = sample_df)


#
#### CORRELATION BETWEEN GEOGRAPHIC DISTANCE AND ENVIRONTMENT AND COMMUNITY ####
env_dist <- as.matrix(dist(apply(sample_df[sample_df$Condition != "NH4", c(12,13,15,16,19,20)], 2, 
                                 function(x) (x-min(x))/(max(x)-min(x)))))

coord_dist.nh4 <- coord_dist[row.names(coord_dist) %in% row.names(env_dist),
                             colnames(coord_dist) %in% colnames(env_dist)]

env_geo <- mantel(env_dist, coord_dist.nh4, method = "spearman", permutations = 9999, na.rm = TRUE)
env_geo

env_geo.plot <- merge(melt(as.matrix(env_dist)), melt(as.matrix(coord_dist.nh4)), by = c("Var1", "Var2"))
colnames(env_geo.plot) <- c("Var1", "Var2", "env_dist", "coord_dist")

ggplot(env_geo.plot) +
  geom_point(aes(x = coord_dist, y = env_dist)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 24, color = "black"),
        axis.text.x = element_text(size = 22, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        axis.text.y = element_text(size = 22, color = "black")) +
  xlab("Geographic Distance") + ylab("Environmental Distance") +
  annotate("text", label = "r = 0.715, p-value < 0.001", x = 400, y = 2, size = 7)

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

#### MUST ANALYSIS ####
must_df <- sample_df[sample_df$Condition != "NH4", c(12,13,15,16,19,20)]
must_df.t <- sample_df[sample_df$Condition != "NH4" & sample_df$Variety != "Garnacha", c(12,13,15,16,19,20)]

must_pca <- prcomp(must_df, scale = TRUE)

must_pca.plot <- must_pca$x
must_pca.plot <- merge(sample_df, must_pca.plot, by = "row.names")

ggplot(must_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming), size = 4) +
  scale_color_manual(values = c("#2ac219", "#19b7c2", "#dba54d", "#ae0e36", 
                                "#433254", "#5900ff", "#7802a3", "#8849d1", "#9d0dd1")) +
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
  scale_color_manual(values = c("#bf2c45", "#1e74eb", "#93bf2c")) +
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

must_df <- merge(must_df, sample_df[,3:5], by = "row.names")

adonis2(as.matrix(dist(must_df[must_df$Condition != "NH4", 2:7])) ~ Condition,
        must_df[must_df$Condition != "NH4", ], permutations = 1000)


#
#### CONSTRAINED ANALYSIS ####

adonis2(as.matrix(bray_GM)[sample_df[sample_df$Condition != "NH4", "Seq_ID"],
                           sample_df[sample_df$Condition != "NH4", "Seq_ID"]] ~ ., 
        sample_df[sample_df$Condition != "NH4", c(12,13,15,16,19,20)], permutations = 1000)

cap_GM <- capscale(as.matrix(bray_GM)[sample_df[sample_df$Condition != "NH4", "Seq_ID"],
                                      sample_df[sample_df$Condition != "NH4", "Seq_ID"]] ~ ., 
                   sample_df[sample_df$Condition != "NH4", c(12,13,15,16,20)], na.action = na.omit)

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
  geom_text(data = arrowmat, aes(x = CAP1*2.25, y = CAP2*2.25, label = Sample), size = 6) +
  scale_color_manual(values = c("#2ac219", "#19b7c2", "#dba54d", "#ae0e36", 
                                "#433254", "#5900ff", "#7802a3", "#8849d1", "#9d0dd1")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        legend.text = element_text(size = 22, color = "black"),
        legend.title = element_text(size = 22, color = "black"),
        axis.text.x = element_text(size = 24, color = "black")) +
  xlab("CAP1 (18.36%)") + ylab("CAP2 (5.78%)")

#
#### SAVE OUTPUTS ####
save.image("com_analysis.RData")

#