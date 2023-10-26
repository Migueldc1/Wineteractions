# Author: M, de Celis Rodriguez
# Date: 17/06/2022
# Project: Wineteractions - ITS sequence analysis of T0-GM samples

library(sf)
library(ggplot2)


rm(list=ls()) #Clear R environment

# Set the project location as working directory
setwd("~/../OneDrive - Universidad Complutense de Madrid (UCM)/Proyecto - Wineteractions/GitHub/Wineteractions/")

#
#### DATA LOADING ####

## SAMPLE DATA
sample_df <- read.table("Data/Metadata/sample_GM.txt", sep = "\t", header = TRUE, encoding = "UTF-8")
sample_df$Origin <- factor(sample_df$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
sample_df$Region <- factor(sample_df$Region, levels = c("Ribera del Guadiana", "Valdepeñas", "La Mancha", "Madrid", "Rioja"))

sample_df <- subset(sample_df, Stage == "0_initial")

#
#### SAMPLING MAP - ACROSS WINE APPELLATIONS ####
# Load the Spanish community boundaries
spain <- geodata::gadm(country = "ESP", level = 1, path = "Data/Map")

sp.gadm <- st_as_sf(spain)
sp.gadm <- sp.gadm[!sp.gadm$NAME_1 %in% c("Ceuta y Melilla", "Islas Canarias", "Islas Baleares"),]


# Draw the map
gg.map_across <- ggplot() +
  geom_sf(data = sp.gadm, aes(fill = NAME_1), color = NA, show.legend = FALSE) +
  scale_fill_manual(values = alpha(c("#7ecf65", "#ccad5e", "#f7de9c", "#f7c67c", "#3c9e43",
                                     "#a16bc7", "#de1e14", "#b00b29", "#ed5a5a", "#88a3bf", 
                                     "#71d1c7", "#f2f291", "#919bf2", "#6c66d1", "#cc74bc"), alpha = 0.75)) +
  ggnewscale::new_scale_fill() +
  geom_point(data = sample_df, aes(x = Longitude, y = Latitude, fill = Origin, shape = Farming), size = 3.5, stroke = 1.5) +
  scale_fill_manual(values = c("#2ac219", "#9fffd1", "#dba54d", "#ae0e36", 
                               "#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")) +
  scale_shape_manual(values = c(21, 23), labels = c("Conventional", "Organic")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "gray20"))) +
  theme_void() +
  theme(legend.position = "top",
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.margin = margin(t = 0, r = 0, b = -0.5, l = 0, unit = "cm"))

gg.map_across

#
#### SAMPLING MAP - WITHIN WINE APPELLATIONS ####

# Draw the map
gg.map_within <- ggplot() +
  geom_sf(data = sp.gadm, aes(fill = NAME_1), color = "gray20", lwd = 0.25, show.legend = FALSE) +
  scale_fill_manual(values = alpha(c("#7ecf65", "#ccad5e", "#f7de9c", "#f7c67c", "#3c9e43",
                                     "#a16bc7", "#de1e14", "#b00b29", "#ed5a5a", "#88a3bf", 
                                     "#71d1c7", "#f2f291", "#919bf2", "#6c66d1", "#cc74bc"), alpha = 0.75)) +
  ggnewscale::new_scale_fill() +
  geom_point(data = subset(sample_df, Origin %in% c("R1", "R2", "R3A", "R3B", "R3C")), show.legend = FALSE, 
             aes(x = Longitude, y = Latitude, fill = Origin, shape = Farming), size = 3.5, stroke = 1.5) +
  scale_fill_manual(values = c("#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")) +
  scale_shape_manual(values = c(21, 23), labels = c("Conventional", "Organic")) +
  theme_void() +
  xlim(-1.896458,-1.859122) + ylim(42.12142,42.14299)

gg.map_within

#
#### EXPORT FIGURES #####
# ACROSS WA MAP
ggsave("Figures/Inkscape/Figure_S1.1.pdf", gg.map_across, 
       width = 10, height = 7, dpi = 300, bg = "white")

# WITHIN WA MAP
ggsave("Figures/Inkscape/Figure_S1.2.pdf", gg.map_within, 
       width = 10, height = 7, dpi = 300, bg = "white")

#