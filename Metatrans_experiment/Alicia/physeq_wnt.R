#### LIBRARIES ####
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(vegan)
library(hillR)
library(igraph)
#library(dada2) #Get nucleotidic sequences
#library(phangorn)
#library(msa) #Perfom multiple alignment

setwd("C:/Users/Laboratorio14/Desktop/Miguel/Wineteractions/Fungi")
load("physeq_fun.RData")

#
#### FUNCTIONS ####
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#
#### DATA PREPROCESSING ####
asv_fun <- readRDS("ASV_cut.rds")
asv_fun <- asv_fun[row.names(asv_fun) %in% meta.ns$Sample_ID, ]
tasv_fun <- asv_fun/rowSums(asv_fun)

tax_fun <- readRDS("tax_cut.rds")
tax_fun <- cbind.data.frame(tax_fun, Id = row.names(tax_fun))
tax_fun$ASV <- paste("ASV", 1:nrow(tax_fun), sep = "-")

meta <- read.table("meta.txt", sep = "\t", header = TRUE)
row.names(meta) <- meta$Sample_ID
meta$name <- paste(meta$Origin, meta$Sample_name)

meta$Time <- factor(meta$Time, levels = c("Initial", "Middle", "Final"))
meta$Origin <- factor(meta$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B",
                                              "R3C"))
meta$Region <- factor(meta$Region, levels = c("RiberaGuadiana", "Madrid", "LaMancha", 
                                              "Valdepe?as", "LaRioja"))

meta.ns <- subset(meta, Condition != "sintetico")

phy_fun <- phyloseq(otu_table(asv_fun, taxa_are_rows = FALSE), tax_table(tax_fun), 
                    sample_data(meta))

phy_fun@tax_table[is.na(phy_fun@tax_table)] <- "Unidentified"

phy_fun <- subset_samples(phy_fun, order > 54)
phyt_fun <- transform_sample_counts(phy_fun, function(x) x / sum(x))

#
#### Phylogenetic tree ####
seqs_fun <- getSequences(colnames(asv_fun))
names(seqs_fun) <- seqs_fun

mult_fun <- msa(seqs_fun, method = "ClustalW", type = "dna", order = "input")

phang.align_fun <- as.phyDat(mult_fun, type = "DNA", names = seqs_fun)
dm_fun <- dist.ml(phang.align_fun)
treeNJ_fun <- NJ(dm_fun) # Note, tip order != sequence order
fit_fun <- pml(treeNJ_fun, data = phang.align_fun)

fitGTR_fun <- update(fit_fun, k = 4, inv = 0.2)
fitGTR_fun <- optim.pml(fitGTR_fun, model = "GTR", optInv = TRUE, optGamma = TRUE,
                         rearrangement = "stochastic", control = pml.control(trace = 0))

#
#### Get most abundant genus at final stage ####
max <- apply(tasv_fun, 1, function(x) which.max(x))
max <- data.frame(max)
max$Sample_ID <- row.names(max)
max$Id <- row.names(xx)[max$max]

max <- merge(max, tax_fun[,c(6,8)])
max <- merge(max, meta, by = "Sample_ID")
max$name <- paste(max$Origin, max$Sample_name)

max <- subset(max, Time == "Final")[,c(14, 4)]

max$most.ab <- ifelse(max$Genus == "g__Lachancea" | max$Genus == "g__Saccharomyces", 
                     as.character(max$Genus), 
                     "Other")

meta <- merge(meta, max[,c(1,3)], by = "name", all.x = TRUE)

#
#### TAXONOMIC ASSIGNATION ####
##Phylum
glomP_fun <- tax_glom(phyt_fun, taxrank = 'Phylum', NArm = FALSE)
dataP_fun <- psmelt(glomP_fun)
dataP_fun$Phylum <- as.character(dataP_fun$Phylum)
dataP_fun$Phylum[dataP_fun$Abundance < 0.01] <- "Other"

orderP <- levels(factor(dataP_fun$Phylum))
orderP <- orderP[! orderP %in% c("Unidentified")]
orderP <- append(orderP, c("Unidentified"))

set.seed(123)
ggplot(dataP_fun, aes(x = Sample, y = Abundance, fill = factor(Phylum, levels = orderP))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Phylum", values = sample(col_vector, 16)) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust=1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Sample") + ylab("Abundance")

##Genus
glomG_fun <- tax_glom(phyt_fun, taxrank = 'Genus', NArm = FALSE)
dataG_fun <- psmelt(glomG_fun)
dataG_fun$Genus <- as.character(dataG_fun$Genus)
dataG_fun$Genus[dataG_fun$Abundance < 0.05] <- "Other"
dataG_fun$Genus[is.na(dataG_fun$Genus)] <- "Unidentified"

orderG <- levels(factor(dataG_fun$Genus))
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

dataG_fun.plot <- aggregate(dataG_fun$Abundance, list(dataG_fun$Sample, dataG_fun$Genus), sum)
colnames(dataG_fun.plot) <- c("Sample_ID", "Genus", "Abundance")

dataG_fun.plot <- merge(dataG_fun.plot, meta, by = "Sample_ID")

set.seed(123)
ggplot(dataG_fun.plot, 
       aes(x = Sample_name, y = Abundance, fill = factor(Genus, levels = orderG))) + 
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
  xlab("Sample") + ylab("Abundance") + facet_wrap(~Time + Origin, nrow = 3) +
  guides(fill = guide_legend(nrow = 3))

#
#### ALPHA DIVERISTY: Dynamics ####
alpha_fun <- cbind.data.frame(# Hill based taxonomic alpha diversity
                              t.q0 = hill_taxa(tasv_fun, q = 0),
                              t.q1 = hill_taxa(tasv_fun, q = 1),
                              t.q2 = hill_taxa(tasv_fun, q = 2),
                              
                              # Hill based phylogenetic alpha diversity
                              p.q0 = hill_phylo(tasv_fun, fitGTR_fun$tree, q = 0),
                              p.q1 = hill_phylo(tasv_fun, fitGTR_fun$tree, q = 1),
                              p.q2 = hill_phylo(tasv_fun, fitGTR_fun$tree, q = 2))

alpha_fun$Sample_ID <- row.names(alpha_fun)

alpha_fun <- merge(alpha_fun, subset(meta, order > 54)[,-c(7,8)], by = "Sample_ID")

alpha_fun.plot <- melt(alpha_fun)
alpha_fun.plot$diversity <- ifelse(startsWith(as.character(alpha_fun.plot$variable), "t."),
                                   "Taxonomic", "Phylogenetic")
alpha_fun.plot$diversity <- factor(alpha_fun.plot$diversity, 
                                   levels = c("Taxonomic", "Phylogenetic"))

alpha_fun.plot$q <- as.factor(substr(alpha_fun.plot$variable, 
                                     nchar(as.character(alpha_fun.plot$variable)),
                                     nchar(as.character(alpha_fun.plot$variable))))

alpha_fun.plot$name <- paste(alpha_fun.plot$Origin, alpha_fun.plot$Sample_name, 
                             alpha_fun.plot$q)

xx <- subset(alpha_fun.plot, diversity == "Taxonomic")
xx <- cbind.data.frame(xx, 
                       value2 = rbind.data.frame(cbind(subset(xx, q == 1)$value),
                                                 cbind(subset(xx, q == 1)$value),
                                                 cbind(subset(xx, q == 1)$value)))
xxx <- cbind.data.frame(xx[,-14], 
                        value2 = rbind.data.frame(cbind(subset(xx, q == 2)$value),
                                                  cbind(subset(xx, q == 2)$value),
                                                  cbind(subset(xx, q == 2)$value)))


ggplot(subset(alpha_fun.plot, diversity == "Taxonomic"), 
       aes(x = Time, y = value, shape = Management, color = Condition, group = name)) + 
  geom_point(position = position_dodge(width = 0.3), size = 2) + 
  geom_line(position = position_dodge(width = 0.3), size = 1) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust = 1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Alpha-diversity") + facet_wrap(~q + Origin, nrow = 3)

#
#### ALPHA DIVERISTY: Normalized Plots ####
alpha_fun.norm <- alpha_fun[,-c(5,6,7)]

alpha_fun.norm$t.q0 <- ave(alpha_fun.norm$t.q0, as.factor(alpha_fun.norm$Origin), 
                           FUN = scale)
alpha_fun.norm$t.q1 <- ave(alpha_fun.norm$t.q1, as.factor(alpha_fun.norm$Origin), 
                           FUN = scale)
alpha_fun.norm$t.q2 <- ave(alpha_fun.norm$t.q2, as.factor(alpha_fun.norm$Origin), 
                           FUN = scale)

alpha_fun.norm_plot <- melt(alpha_fun.norm)

ggplot(alpha_fun.norm_plot, 
       aes(x = Time, y = value, color = Management)) + 
  geom_boxplot(position = "dodge", size = 2) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust = 1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Alpha-diversity") + facet_wrap(~variable + Condition, nrow = 3)

#
#### ALPHA DIVERSITY: Normalized dynamics plot ####
ggplot(alpha_fun.norm_plot, 
       aes(x = Time, y = value, shape = Management, color = Condition, group = Sample_name)) + 
  geom_point(position = position_dodge(width = 0.3), size = 2) + 
  geom_line(position = position_dodge(width = 0.3), size = 1) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust = 1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Alpha-diversity") + facet_wrap(~variable + Origin, nrow = 3)



#
#### ALPHA DIVERISTY: Diversity loss: Raw / Raw.prop ####
alpha_fun.plot <- subset(alpha_fun.plot, diversity == "Taxonomic")

alpha_fun.plot$name2 <- paste(alpha_fun.plot$Origin, alpha_fun.plot$Sample_name,
                              alpha_fun.plot$variable)

alpha_fun.dyn <- subset(alpha_fun.plot, Time == "Initial")
colnames(alpha_fun.dyn)[10] <- "Initial"

alpha_fun.dyn <- merge(alpha_fun.dyn, subset(alpha_fun.plot, Time == "Middle")[,c(10,14)],
                       by = "name2", all.x = TRUE)
colnames(alpha_fun.dyn)[15] <- "Middle"

alpha_fun.dyn <- merge(alpha_fun.dyn, subset(alpha_fun.plot, Time == "Final")[,c(10,14)],
                       by = "name2", all.x = TRUE)
colnames(alpha_fun.dyn)[16] <- "Final"

alpha_fun.dyn$`Initial-Final` <- (alpha_fun.dyn$Initial-alpha_fun.dyn$Final)/alpha_fun.dyn$Initial
alpha_fun.dyn$`Initial-Middle` <- (alpha_fun.dyn$Initial-alpha_fun.dyn$Middle)/alpha_fun.dyn$Initial

alpha_fun.dyn <- alpha_fun.dyn[,c(3,4,6,8,10,12,13,17,18)]
alpha_fun.dyn.plot <- melt(alpha_fun.dyn)
colnames(alpha_fun.dyn.plot)[5] <- "div.index"

alpha_fun.dyn.plot$variable <- factor(alpha_fun.dyn.plot$variable, 
                                      levels = c("Initial-Middle", "Initial-Final"))

ggplot(subset(xx, diversity == "Taxonomic"), 
       aes(x = Condition, y = value, color = Management)) + 
  geom_boxplot(position = "dodge", size = 2) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Diversity Change") + facet_wrap(~div.index + variable, nrow = 3)


xxx <- rbind.data.frame(cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q0"),
                       cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q1"),
                       cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                               div.index == "t.q2")[,-5], div.index = "t.q2"))


#
#### ALPHA DIVERISTY: Diversity loss: Normalized ####
alpha_fun.norm_plot$name <- paste(alpha_fun.norm_plot$Origin, alpha_fun.norm_plot$Sample_name,
                                  alpha_fun.norm_plot$variable)

alpha_fun.dyn_norm <- subset(alpha_fun.norm_plot, Time == "Initial")
colnames(alpha_fun.dyn_norm)[10] <- "Initial"

alpha_fun.dyn_norm <- merge(alpha_fun.dyn_norm, 
                            subset(alpha_fun.norm_plot, Time == "Middle")[,c(10,11)],
                            by = "name", all.x = TRUE)
colnames(alpha_fun.dyn_norm)[12] <- "Middle"

alpha_fun.dyn_norm <- merge(alpha_fun.dyn_norm, 
                            subset(alpha_fun.norm_plot, Time == "Final")[,c(10,11)],
                            by = "name", all.x = TRUE)
colnames(alpha_fun.dyn_norm)[13] <- "Final"

alpha_fun.dyn_norm$`Initial-Final` <- (alpha_fun.dyn_norm$Initial-alpha_fun.dyn_norm$Final)
alpha_fun.dyn_norm$`Initial-Middle` <- (alpha_fun.dyn_norm$Initial-alpha_fun.dyn_norm$Middle)

alpha_fun.dyn_norm <- alpha_fun.dyn_norm[,c(1:10,14,15)]
alpha_fun.dyn_norm.plot <- melt(alpha_fun.dyn_norm)
colnames(alpha_fun.dyn_norm.plot)[10] <- "div.index"

alpha_fun.dyn_norm.plot$variable <- factor(alpha_fun.dyn_norm.plot$variable, 
                                      levels = c("Initial-Middle", "Initial-Final"))

ggplot(alpha_fun.dyn_norm.plot, 
       aes(x = Condition, y = value, color = Management)) + 
  geom_boxplot(position = "dodge", size = 2) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Diversity Change") + facet_wrap(~div.index + variable, nrow = 3)


xxx <- rbind.data.frame(cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q0"),
                        cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q1"),
                        cbind.data.frame(subset(alpha_fun.dyn.plot, 
                                                div.index == "t.q2")[,-5], div.index = "t.q2"))


#
#### BETA DIVERSITY ####
phyt_bray <- phyloseq::distance(phyt_fun, method = "bray")
nMDS_fun <- metaMDS(phyt_bray)
nMDS_fun$stress

nMDS_fun.plot <- as.data.frame(scores(nMDS_fun))
nMDS_fun.plot$Sample_ID <- rownames(nMDS_fun.plot)
nMDS_fun.plot <- merge(nMDS_fun.plot, meta, by = "Sample_ID")
nMDS_fun.plot$group <- paste(nMDS_fun.plot$Origin, nMDS_fun.plot$Sample_name)
nMDS_fun.plot$order <- as.numeric(nMDS_fun.plot$Time)
nMDS_fun.plot <- nMDS_fun.plot[order(nMDS_fun.plot$order),]

ggplot(nMDS_fun.plot) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Time, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))


pcoa_fun <- wcmdscale(phyt_bray)

pcoa_fun.plot <- as.data.frame(scores(pcoa_fun))
pcoa_fun.plot$Sample_ID <- rownames(pcoa_fun.plot)
pcoa_fun.plot <- merge(pcoa_fun.plot, meta, by = "Sample_ID")
pcoa_fun.plot$group <- paste(pcoa_fun.plot$Origin, pcoa_fun.plot$Sample_name)
pcoa_fun.plot$order <- as.numeric(pcoa_fun.plot$Time)
pcoa_fun.plot <- pcoa_fun.plot[order(pcoa_fun.plot$order),]

ggplot(pcoa_fun.plot) + 
  geom_point(aes(x = Dim1, y = Dim2, color = Time, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))


set.seed(1)
adonis(phyt_bray ~ Time + Region + Management:Region, data = meta.ns)

meta.ns <- meta.ns[match(sample_names(phyt_fun), meta.ns$Sample_ID),]

meta.phy <- data.frame(sample_data(phy_fun))


beta <- betadisper(phyt_bray, meta.ns$Time, type = "centroid")
permutest(beta)
plot(beta)

meta.ns$Sample_ID

#
#### BETA DIVERSITY: CONSTRAINED ####
adonis2(as.matrix(phyt_bray) ~ Region + Management + Time + Condition,
        meta.phy, permutations = 1000)

fun_cap <- capscale(as.matrix(phyt_bray) ~ ., meta.phy[,c(9,2,3)], na.action = na.omit)

anova(fun_cap)
RsquareAdj(fun_cap)

arrowmat <- rbind.data.frame(scores(fun_cap, display = "bp")[-c(6:7),],
                             scores(fun_cap, display = "cn"))
arrowmat$Sample_ID <- row.names(arrowmat)

sitemat <- as.data.frame(scores(fun_cap, display = "sites"))
sitemat$Sample_ID <- row.names(sitemat)
sitemat <- merge(sitemat, meta.phy, by = "Sample_ID")

summary(fun_cap)$cont

ggplot() + 
  geom_point(data = sitemat, aes(x = CAP1, y = CAP2, color = Time), size = 3) +
  geom_segment(data = arrowmat[-c(6:8),], aes(x = 0, y = 0, xend = CAP1*2, yend = CAP2*2), 
               arrow = arrow(length = unit(0.2,"cm")), size = 1, alpha = 0.8) + 
  geom_text(data = arrowmat, aes(x = CAP1*2, y = CAP2*2, label = Sample_ID), size = 5) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 17, color = "black")) +
  xlab("CAP1 (38.43%)") + ylab("CAP2 (13.06%)")

#
#### BETA DIVERSITY: Initial ####
phyt_fun.filt <- subset_samples(phyt_fun, Time == "Initial")
meta.filt <- data.frame(sample_data(phyt_fun.filt))

phyt_bray.filt <- phyloseq::distance(phyt_fun.filt, method = "bray")

nMDS_fun <- metaMDS(phyt_bray.filt)
nMDS_fun$stress

nMDS_fun.plot <- as.data.frame(scores(nMDS_fun))
nMDS_fun.plot$Sample_ID <- rownames(nMDS_fun.plot)
nMDS_fun.plot <- merge(nMDS_fun.plot, meta, by = "Sample_ID")
nMDS_fun.plot$group <- paste(nMDS_fun.plot$Origin, nMDS_fun.plot$Sample_name)
nMDS_fun.plot$order <- as.numeric(nMDS_fun.plot$Time)
nMDS_fun.plot <- nMDS_fun.plot[order(nMDS_fun.plot$order),]

ggplot(nMDS_fun.plot) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Management, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

#
## Beta-dispersion boxplots
bd_fun <- betadisper(phyt_bray.filt, meta.filt$Management, type = "centroid")

bd_fun <- data.frame(Distance_to_centroid = bd_fun$distances)
bd_fun$Sample_ID <- row.names(bd_fun)

bd_fun <- merge(bd_fun, meta, by = "Sample_ID")

#BETADISPER BOXPLOT
ggplot(bd_fun, aes(x = Management, y = Distance_to_centroid, color = Management)) + 
  geom_boxplot(size = 1)

summary(aov(Distance_to_centroid ~ Management, data = bd_fun))

set.seed(1)
adonis(phyt_bray ~ Region*Management, data = meta.ns)



#
##PLOTS CON LAS LINEAS
ggplot(nMDS_fun.plot[order(nMDS_fun.plot$order),]) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Condition, shape = Management, size = Time)) +
  geom_line(aes(x = NMDS1, y = NMDS2, group = group, color = Condition)) +
  scale_size_manual(values = c(5, 3, 2)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

ggplot(subset(nMDS_fun.plot, Time == "Initial")) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = Condition, shape = Management), size = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

#
#### NETWORK ANALYSIS: METAWEB - t0 ####
#Transform our dataset into presence/absence table
otu_prob <- tasv_fun[subset(meta, Time == "Initial")$Sample_ID, ]
otu_prob <- otu_prob[, colSums(otu_prob > 0) > 1 & colSums(otu_prob > 0.001) > 0]

otu_prob[otu_prob > 0] <- 1

#Calculate how many times each ASV/OTU pair co-occur (appear in the same sample)
bac_prob <- as.matrix(otu_prob)
occurs_prob <- crossprod((bac_prob)) 

#Calculate probabilitity of co-occurrence for each pair
res <- matrix(NA, ncol(bac_prob), ncol(bac_prob))

for (i in 1:nrow(occurs_prob)){
  for (j in 2:ncol(occurs_prob)){
    if (j <= i) next
    out <- 0
    cooc <- occurs_prob[i,j] - 1 
    for (k in 0:cooc){ 
      out <- out + Veech_calc(nrow(bac_prob), k, occurs_prob[j,j], occurs_prob[i,i]) 
    }
    res[i, j] <- out # out = probability of finding less than k cooccurences. 
  }
}

#Filter only significant pairs (probability >= 0.95)
pos_prob <- subset(melt(res), value >= 0.95)
colnames(pos_prob) <- c("Source", "Target", "value")
pos_prob[, 1] <- colnames(bac_prob)[pos_prob[, 1]]
pos_prob[, 2] <- colnames(bac_prob)[pos_prob[, 2]]

#Construct node tables with taxonomic information of the nodes in the network
node_prob <- unique(rbind(cbind(Id = unique(as.character(pos_prob[,1]))), 
                          cbind(Id = unique(as.character(pos_prob[,2])))))
node_prob <- merge(node_prob, tax_fun, by = "Id")
node_prob[is.na(node_prob)] <- "Unidentified"

##DRAW NETWORK
#Import the egde list (list of significant co-occurring pair) to the igraph package
metanet_prob <- graph_from_data_frame(as.matrix(pos_prob[,c(1:2)]), directed = FALSE)

node_prob <- node_prob[match(names(V(metanet_prob)), node_prob$Id), ]

node_prob$most.ab <- ifelse(node_prob$Genus == "g__Lachancea" | 
                              node_prob$Genus == "g__Saccharomyces", 
                            as.character(node_prob$Genus), "Other")

V(metanet_prob)$ASV <- node_prob$ASV
V(metanet_prob)$Genus <- node_prob$Genus
V(metanet_prob)$color <- as.factor(node_prob$most.ab)

#Set the layout, the position in the plot qhere the nodes will be drawn, attending to their connections
set.seed(1)
l_prob <- layout_with_fr(metanet_prob, grid = "nogrid")
l_prob <- norm_coords(l_prob, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

#Calculate the number of modules (subcommunities) forming the network
cwt_prob <- cluster_walktrap(metanet_prob, steps = 20)
max(cwt_prob$membership)
V(metanet_prob)$community <- cwt_prob$membership
colrs <- adjustcolor( c("yellowgreen", "darkgoldenrod", "deepskyblue", "yellow"), alpha = 0.8)

#Draw the network
plot(metanet_prob, vertex.size = 5, vertex.label = NA, 
     vertex.color = colrs[as.numeric(V(metanet_prob)$community)], 
     rescale = FALSE, layout = l_prob*1.3, edge.color = "gray70")

#
xx <- cbind.data.frame(Id = cwt_prob$names, mod = cwt_prob$membership)
node_prob <- merge(node_prob, xx, by = "Id")

xxx.pro <- subset(node_prob, mod == 3)


xxx.x <- xxx.pro[xxx.pro$ASV %in% xxx$ASV,]
xxx.x <- merge(xxx.x, xxx[,c(9,12)], by = "ASV")

occurs_xxx <- occurs_prob[xxx.x$Id, xxx.x$Id]

node_prob$remove <- ifelse(node_prob$ASV %in% xxx.x$ASV, "yellow", NA)
V(metanet_prob)$remove <- node_prob$remove

otu_xxx <- asv_fun[,colnames(asv_fun) %in% xxx.x$Id]

#
#### NETWORK ANALYSIS: Eco - t0 ####
#Transform our dataset into presence/absence table
otu_eco <- tasv_fun[subset(meta, Time == "Initial" & Management == "eco")$Sample_ID, ]
otu_eco <- otu_eco[, colSums(otu_eco > 0) > 1 & colSums(otu_eco > 0.001) > 0]

otu_eco[otu_eco > 0] <- 1

#Calculate how many times each ASV/OTU pair co-occur (appear in the same sample)
bac_eco <- as.matrix(otu_eco)
occurs_eco <- crossprod((bac_eco)) 

#Calculate probabilitity of co-occurrence for each pair
res <- matrix(NA, ncol(bac_eco), ncol(bac_eco))

for (i in 1:nrow(occurs_eco)){
  for (j in 2:ncol(occurs_eco)){
    if (j <= i) next
    out <- 0
    cooc <- occurs_eco[i,j] - 1 
    for (k in 0:cooc){ 
      out <- out + Veech_calc(nrow(bac_eco), k, occurs_eco[j,j], occurs_eco[i,i]) 
    }
    res[i, j] <- out # out = probability of finding less than k cooccurences. 
  }
}

#Filter only significant pairs (probability >= 0.95)
pos_eco <- subset(melt(res), value >= 0.95)
colnames(pos_eco) <- c("Source", "Target", "value")
pos_eco[, 1] <- colnames(bac_eco)[pos_eco[, 1]]
pos_eco[, 2] <- colnames(bac_eco)[pos_eco[, 2]]

#Construct node tables with taxonomic information of the nodes in the network
node_eco <- unique(rbind(cbind(Id = unique(as.character(pos_eco[,1]))), 
                          cbind(Id = unique(as.character(pos_eco[,2])))))
node_eco <- merge(node_eco, tax_fun, by = "Id")
node_eco[is.na(node_eco)] <- "Unidentified"

##DRAW NETWORK
#Import the egde list (list of significant co-occurring pair) to the igraph package
metanet_eco <- graph_from_data_frame(as.matrix(pos_eco[,c(1:2)]), directed = FALSE)

node_eco <- node_eco[match(names(V(metanet_eco)), node_eco$Id), ]

node_eco$most.ab <- ifelse(node_eco$Genus == "g__Lachancea" | 
                              node_eco$Genus == "g__Saccharomyces", 
                            as.character(node_eco$Genus), "Other")

V(metanet_eco)$ASV <- node_eco$ASV
V(metanet_eco)$Genus <- node_eco$Genus
V(metanet_eco)$color <- as.factor(node_eco$most.ab)

#Set the layout, the position in the plot qhere the nodes will be drawn, attending to their connections
set.seed(1)
l_eco <- layout_with_fr(metanet_eco, grid = "nogrid")
l_eco <- norm_coords(l_eco, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

#Calculate the number of modules (subcommunities) forming the network
cwt_eco <- cluster_walktrap(metanet_eco, steps = 20)
max(cwt_eco$membership)
V(metanet_eco)$community <- cwt_eco$membership
colrs <- adjustcolor( c("yellowgreen", "darkgoldenrod", "deepskyblue"), alpha = 0.8)

#Draw the network
plot(metanet_eco, vertex.size = 5, vertex.label = NA, 
     vertex.color = colrs[as.numeric(V(metanet_eco)$community)], 
     rescale = FALSE, layout = l_eco*1.3, edge.color = "gray70")

#OO
#### NETWORK ANALYSIS: Conv - t0 ####
#Transform our dataset into presence/absence table
otu_conv <- tasv_fun[subset(meta, Time == "Initial" & Management == "conv")$Sample_ID, ]
otu_conv <- otu_conv[, colSums(otu_conv > 0) > 1 & colSums(otu_conv > 0.001) > 0]

otu_conv[otu_conv > 0] <- 1

#Calculate how many times each ASV/OTU pair co-occur (appear in the same sample)
bac_conv <- as.matrix(otu_conv)
occurs_conv <- crossprod((bac_conv)) 

#Calculate probabilitity of co-occurrence for each pair
res <- matrix(NA, ncol(bac_conv), ncol(bac_conv))

for (i in 1:nrow(occurs_conv)){
  for (j in 2:ncol(occurs_conv)){
    if (j <= i) next
    out <- 0
    cooc <- occurs_conv[i,j] - 1 
    for (k in 0:cooc){ 
      out <- out + Veech_calc(nrow(bac_conv), k, occurs_conv[j,j], occurs_conv[i,i]) 
    }
    res[i, j] <- out # out = probability of finding less than k cooccurences. 
  }
}

#Filter only significant pairs (probability >= 0.95)
pos_conv <- subset(melt(res), value >= 0.95)
colnames(pos_conv) <- c("Source", "Target", "value")
pos_conv[, 1] <- colnames(bac_conv)[pos_conv[, 1]]
pos_conv[, 2] <- colnames(bac_conv)[pos_conv[, 2]]

#Construct node tables with taxonomic information of the nodes in the network
node_conv <- unique(rbind(cbind(Id = unique(as.character(pos_conv[,1]))), 
                         cbind(Id = unique(as.character(pos_conv[,2])))))
node_conv <- merge(node_conv, tax_fun, by = "Id")
node_conv[is.na(node_conv)] <- "Unidentified"

##DRAW NETWORK
#Import the egde list (list of significant co-occurring pair) to the igraph package
metanet_conv <- graph_from_data_frame(as.matrix(pos_conv[,c(1:2)]), directed = FALSE)

node_conv <- node_conv[match(names(V(metanet_conv)), node_conv$Id), ]

node_conv$most.ab <- ifelse(node_conv$Genus == "g__Lachancea" | 
                             node_conv$Genus == "g__Saccharomyces", 
                           as.character(node_conv$Genus), "Other")

V(metanet_conv)$ASV <- node_conv$ASV
V(metanet_conv)$Genus <- node_conv$Genus
V(metanet_conv)$color <- as.factor(node_conv$most.ab)

#Set the layout, the position in the plot qhere the nodes will be drawn, attending to their connections
set.seed(1)
l_conv <- layout_with_fr(metanet_conv, grid = "nogrid")
l_conv <- norm_coords(l_conv, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

#Calculate the number of modules (subcommunities) forming the network
cwt_conv <- cluster_walktrap(metanet_conv, steps = 20)
max(cwt_conv$membership)
V(metanet_conv)$community <- cwt_conv$membership
colrs <- adjustcolor(c("yellowgreen", "darkgoldenrod", "deepskyblue", "yellowgreen", 
                        "darkgoldenrod", "deepskyblue", "yellowgreen", "blue", 
                        "deepskyblue", "yellowgreen", "deepskyblue", "yellowgreen",
                       "deepskyblue", "yellowgreen"), alpha = 0.8)

#Draw the network
plot(metanet_conv, vertex.size = 5, vertex.label = V(metanet_conv)$community, 
     vertex.color = colrs[as.numeric(V(metanet_conv)$community)], 
     rescale = FALSE, layout = l_conv*1.3, edge.color = "gray70")

#
#### NETWORK ANALYSIS: Local network - t0 ####
pairs_prob <- rbind.data.frame(pos_prob, 
                               cbind.data.frame(Source = node_prob[,1], 
                                                Target = node_prob[,1], value = 0))

globalmatrix_prob <- dcast(pairs_prob, Source~Target)
row.names(globalmatrix_prob) <- globalmatrix_prob[,1]
globalmatrix_prob <- globalmatrix_prob[,-1]
globalmatrix_prob[is.na(globalmatrix_prob)] <- 0
globalmatrix_prob[globalmatrix_prob != 0] <- 1

bac_prob <- bac_prob[, names(globalmatrix_prob)]

sample.list_prob <- lapply(seq_len(nrow(bac_prob)), function(i) bac_prob[i,])

matrixlist_prob <- vector("list", nrow(bac_prob))

for(i in 1:length(sample.list_prob)) {
  subsetvector <- sample.list_prob[[i]]
  subsetmatrix <- globalmatrix_prob[subsetvector == 1, subsetvector == 1]
  matrixlist_prob[[i]] <- subsetmatrix
}

net.list_prob <- vector("list", nrow(bac_prob))
names(net.list_prob) <- row.names(bac_prob)

for(i in 1:length(matrixlist_prob)) {
  a <- matrixlist_prob[[i]]
  c <- melt(as.matrix(a))
  d <- subset(c, c$value != 0)
  e <- graph_from_data_frame(as.matrix(d[,c(1:2)]), directed = FALSE)
  net.list_prob[[i]] <- e
}

rm(a); rm(c); rm(d); rm(e)

#
#### NETWORK ANALYSIS: Local network properties - t0 ####
## Module completeness
mods_prob <- data.frame(cbind(Id = cwt_prob$names, modularity_class = cwt_prob$membership))

mod.prop_prob <- NULL

for (i in 1:length(net.list_prob)) {
  nodes <- cbind(Id = get.vertex.attribute(net.list_prob[[i]])$name)
  nodes <- merge(nodes, mods_prob, by = "Id")
  mod.prop_prob <- rbind(mod.prop_prob, cbind(`mod 1` = nrow(nodes[nodes$modularity_class == 1, ])/nrow(nodes),
                                              `mod 2` = nrow(nodes[nodes$modularity_class == 2, ])/nrow(nodes),
                                              `mod 3` = nrow(nodes[nodes$modularity_class == 3, ])/nrow(nodes),
                                              `mod 4` = nrow(nodes[nodes$modularity_class == 4, ])/nrow(nodes)))
  row.names(mod.prop_prob)[i] = names(net.list_prob[i])
}

mod.prop_prob <- melt(mod.prop_prob)
colnames(mod.prop_prob) <- c("Sample_ID", "modularity_class", "value")

mod.prop_prob <- merge(mod.prop_prob, meta, by = "Sample_ID")

ggplot(data = mod.prop_prob) +
  geom_boxplot(aes(x = modularity_class, y = value), size = 1) + 
  #Make the plot pretty
  theme_bw() +
  theme(axis.text.y = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black"),
        axis.text.x = element_text(size = 17, color = "black", angle = 60, hjust = 1)) + 
  ylim(0, 1) + xlab("") + ylab("Module completeness") + 
  labs(color = "Modularity\nclass", fill = "Modularity\nclass") 

t.test(subset(mod.prop_prob, Management == "eco" & modularity_class == "mod 4")$value,
       subset(mod.prop_prob, Management == "conv" & modularity_class == "mod 4")$value)

#
## Topological network properties
net.prop_prob <- as.data.frame(t(Reduce(function(...) merge(..., all = TRUE), 
                                        list(
                                          #Modularity
                                          lapply(net.list_prob,function(x){modularity(cluster_walktrap(x, steps = 20))}),
                                          #Clustering
                                          lapply(net.list_prob, function(x){transitivity(x, type = "average")}),
                                          #Edge density
                                          lapply(net.list_prob, function(x){edge_density(x)})))))

colnames(net.prop_prob) <- c("Clustering", "Modularity", "Edge.density")

#Data massaging for pretty plots
row.names(net.prop_prob) <- gsub("\\.", "-", row.names(net.prop_prob))
net.prop_prob$Sample_ID <- row.names(net.prop_prob)

net.prop_prob <- merge(net.prop_prob, meta[,-c(8,9)], by = "Sample_ID")

net.prop_prob.plot <- melt(net.prop_prob)

ggplot(data = subset(net.prop_prob.plot, variable == "Modularity")) +
  geom_boxplot(aes(x = Region, y = value, color = Management), size = 1) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black"),
        axis.text.x = element_text(size = 17, color = "black", angle = 60, hjust = 1)) + 
  xlab("") + ylab("Network Properties")

#
#### NETWORK ANALYSIS: Module completeness - t0 + t0-t20-t1000 ####
asv_fmod <- asv_fun[,colnames(asv_fun) %in% mods_prob$Id]
asv_fmod[asv_fmod > 1] <- 1

asv_fmod <- cbind.data.frame(Id = colnames(asv_fmod), t(asv_fmod))
asv_fmod <- merge(asv_fmod, mods_prob, by = "Id")
mods_fmod <- aggregate(asv_fmod[,2:214], list(asv_fmod$modularity_class), sum)

mods_fmod <- as.data.frame(t(mods_fmod[,-1]))
colnames(mods_fmod) <- c("mod 1", "mod 2", "mod 3", "mod 4")

mods_fmod <- cbind.data.frame(mods_fmod, total = rowSums(mods_fmod))
mods_fmod <- mods_fmod[!(mods_fmod$total == 0),]

mods_fmod.prop <- mods_fmod[,1:4]/mods_fmod$total
mods_fmod.prop$Sample_ID <- row.names(mods_fmod.prop)

mods_fmod.prop <- melt(mods_fmod.prop)
colnames(mods_fmod.prop) <- c("Sample_ID", "modularity_class", "value")

mods_fmod.prop <- merge(mods_fmod.prop, meta, by = "Sample_ID")

ggplot(data = mods_fmod.prop) +
  geom_boxplot(aes(x = modularity_class, y = value, color = Time), size = 1) + 
  #Make the plot pretty
  theme_bw() +
  theme(axis.text.y = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black"),
        axis.text.x = element_text(size = 17, color = "black", angle = 60, hjust = 1)) + 
  ylim(0, 1) + xlab("") + ylab("Module completeness") + 
  labs(color = "Modularity\nclass", fill = "Modularity\nclass") 

t.test(subset(mod.prop_prob, Management == "eco" & modularity_class == "mod 4")$value,
       subset(mod.prop_prob, Management == "conv" & modularity_class == "mod 4")$value)

#
## Topological network properties
net.prop_prob <- as.data.frame(t(Reduce(function(...) merge(..., all = TRUE), 
                                        list(
                                          #Modularity
                                          lapply(net.list_prob,function(x){modularity(cluster_walktrap(x, steps = 20))}),
                                          #Clustering
                                          lapply(net.list_prob, function(x){transitivity(x, type = "average")}),
                                          #Edge density
                                          lapply(net.list_prob, function(x){edge_density(x)})))))

colnames(net.prop_prob) <- c("Clustering", "Modularity", "Edge.density")

#Data massaging for pretty plots
row.names(net.prop_prob) <- gsub("\\.", "-", row.names(net.prop_prob))
net.prop_prob$Sample_ID <- row.names(net.prop_prob)

net.prop_prob <- merge(net.prop_prob, meta[,-c(8,9)], by = "Sample_ID")

net.prop_prob.plot <- melt(net.prop_prob)

ggplot(data = subset(net.prop_prob.plot, variable == "Modularity")) +
  geom_boxplot(aes(x = Region, y = value, color = Management), size = 1) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black"),
        axis.text.x = element_text(size = 17, color = "black", angle = 60, hjust = 1)) + 
  xlab("") + ylab("Network Properties")

#
#### t0 MODULES: Taxonomy ####
mods_gen <- merge(tax_fun, mods_prob, by = "Id")
mods_gen <- cbind.data.frame(mods_gen, ASV.count = 1)
mods_gen[is.na(mods_gen)] <- "Unidentified"

mods_gen <- aggregate(mods_gen$ASV.count, list(mods_gen$Genus, mods_gen$modularity_class),
                       sum)

colnames(mods_gen) <- c("Genus", "Module", "ASV.count")

mods_tot <- aggregate(mods_gen$ASV.count, list(mods_gen$Module), sum)
colnames(mods_tot) <- c("Module", "ASV.tot")

mods_gen <- merge(mods_gen, mods_tot, by = "Module")
mods_gen$ASV.prop <- mods_gen$ASV.count/mods_gen$ASV.tot

set.seed(123)
ggplot(mods_gen, aes(x = Module, y = ASV.prop, fill = Genus)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Phylum", values = sample(col_vector, 66)) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Sample") + ylab("Abundance")


#
#### NETWORK ANALYSIS: METAWEB - total ####
#Transform our dataset into presence/absence table
otu_prob <- tasv_fun[row.names(tasv_fun) %in% meta.ns$Sample_ID, ]
otu_prob <- otu_prob[, colSums(otu_prob > 0) > 1 & colSums(otu_prob > 0.001) > 0]

otu_prob[otu_prob > 0] <- 1

#Calculate how many times each ASV/OTU pair co-occur (appear in the same sample)
bac_prob <- as.matrix(otu_prob)
occurs_prob <- crossprod((bac_prob)) 

#Calculate probabilitity of co-occurrence for each pair
res <- matrix(NA, ncol(bac_prob), ncol(bac_prob))

for (i in 1:nrow(occurs_prob)){
  for (j in 2:ncol(occurs_prob)){
    if (j <= i) next
    out <- 0
    cooc <- occurs_prob[i,j] - 1 
    for (k in 0:cooc){ 
      out <- out + Veech_calc(nrow(bac_prob), k, occurs_prob[j,j], occurs_prob[i,i]) 
    }
    res[i, j] <- out # out = probability of finding less than k cooccurences. 
  }
}

#Filter only significant pairs (probability >= 0.95)
pos_prob <- subset(melt(res), value >= 0.95)
colnames(pos_prob) <- c("Source", "Target", "value")
pos_prob[, 1] <- colnames(bac_prob)[pos_prob[, 1]]
pos_prob[, 2] <- colnames(bac_prob)[pos_prob[, 2]]

#Construct node tables with taxonomic information of the nodes in the network
node_prob <- unique(rbind(cbind(Id = unique(as.character(pos_prob[,1]))), 
                          cbind(Id = unique(as.character(pos_prob[,2])))))
node_prob <- merge(node_prob, tax_fun, by = "Id")
node_prob[is.na(node_prob)] <- "Unidentified"

##DRAW NETWORK
#Import the egde list (list of significant co-occurring pair) to the igraph package
metanet_prob <- graph_from_data_frame(as.matrix(pos_prob[,c(1:2)]), directed = FALSE)

node_prob <- node_prob[match(names(V(metanet_prob)), node_prob$Id), ]

node_prob$most.ab <- ifelse(node_prob$Genus == "g__Lachancea" | 
                              node_prob$Genus == "g__Saccharomyces", 
                            as.character(node_prob$Genus), "Other")

V(metanet_prob)$ASV <- node_prob$ASV
V(metanet_prob)$Genus <- node_prob$Genus
V(metanet_prob)$color <- as.factor(node_prob$most.ab)

#Set the layout, the position in the plot qhere the nodes will be drawn, attending to their connections
set.seed(1)
l_prob <- layout_with_fr(metanet_prob, grid = "nogrid")
l_prob <- norm_coords(l_prob, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

#Calculate the number of modules (subcommunities) forming the network
cwt_prob <- cluster_walktrap(metanet_prob, steps = 20)
max(cwt_prob$membership)
V(metanet_prob)$community <- cwt_prob$membership
colrs <- adjustcolor( c("yellowgreen", "darkgoldenrod"), alpha = 0.8)

#Draw the network
plot(metanet_prob, vertex.size = 5, vertex.label = V(metanet_prob)$community, 
     vertex.color = colrs[as.numeric(V(metanet_prob)$community)], 
     rescale = FALSE, layout = l_prob*1.3, edge.color = "gray70")

#
#### NETWORK ANALYSIS: Local network - total ####
pairs_prob <- rbind.data.frame(pos_prob, 
                               cbind.data.frame(Source = node_prob[,1], 
                                                Target = node_prob[,1], value = 0))

globalmatrix_prob <- dcast(pairs_prob, Source~Target)
row.names(globalmatrix_prob) <- globalmatrix_prob[,1]
globalmatrix_prob <- globalmatrix_prob[,-1]
globalmatrix_prob[is.na(globalmatrix_prob)] <- 0
globalmatrix_prob[globalmatrix_prob != 0] <- 1

bac_prob <- bac_prob[, names(globalmatrix_prob)]

sample.list_prob <- lapply(seq_len(nrow(bac_prob)), function(i) bac_prob[i,])

matrixlist_prob <- vector("list", nrow(bac_prob))
names(matrixlist_prob) <- row.names(bac_prob)

for(i in 1:length(sample.list_prob)) {
  subsetvector <- sample.list_prob[[i]]
  subsetmatrix <- globalmatrix_prob[subsetvector == 1, subsetvector == 1]
  matrixlist_prob[[i]] <- subsetmatrix
}

net.list_prob <- vector("list", nrow(bac_prob))
names(net.list_prob) <- row.names(bac_prob)

for(i in 1:length(matrixlist_prob)) {
  a <- matrixlist_prob[[i]]
  c <- melt(as.matrix(a))
  d <- subset(c, c$value != 0)
  e <- graph_from_data_frame(as.matrix(d[,c(1:2)]), directed = FALSE)
  net.list_prob[[i]] <- e
}

rm(a); rm(c); rm(d); rm(e)

#
#### NETWORK ANALYSIS: Local network properties - total ####
## Module completeness
mods_prob <- data.frame(cbind(Id = cwt_prob$names, modularity_class = cwt_prob$membership))

mod.prop_prob <- NULL

for (i in 1:length(net.list_prob)) {
  if (length(get.vertex.attribute(net.list_prob[[i]])) != 0) {
    
    nodes <- cbind(Id = get.vertex.attribute(net.list_prob[[i]])$name)
    nodes <- merge(nodes, mods_prob, by = "Id")
    mod.prop_prob <- rbind(mod.prop_prob, cbind(`mod 1` = nrow(nodes[nodes$modularity_class == 1, ])/nrow(nodes),
                                                `mod 2` = nrow(nodes[nodes$modularity_class == 2, ])/nrow(nodes)))
  }else {
    mod.prop_prob <- rbind(mod.prop_prob, NA)
  }
  
  row.names(mod.prop_prob)[i] = names(net.list_prob[i])
}

mod.prop_prob <- melt(mod.prop_prob)
colnames(mod.prop_prob) <- c("Sample_ID", "modularity_class", "value")

mod.prop_prob <- merge(mod.prop_prob, meta, by = "Sample_ID")

ggplot(data = mod.prop_prob) +
  geom_boxplot(aes(x = Region, y = value, color = modularity_class), size = 1) + 
  #Make the plot pretty
  scale_color_manual(values = c("yellowgreen", "darkgoldenrod")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black"),
        axis.text.x = element_text(size = 17, color = "black", angle = 60, hjust = 1)) + 
  ylim(0, 1) + xlab("") + ylab("Module completeness") + 
  labs(color = "Modularity\nclass", fill = "Modularity\nclass") 

t.test(subset(mod.prop_prob, Management == "eco" & modularity_class == "mod 4")$value,
       subset(mod.prop_prob, Management == "conv" & modularity_class == "mod 4")$value)

mods_prob <- merge(mods_prob, tax_fun, by = "Id")
View(subset(mods_prob, modularity_class == 2))

nrow(subset(mods_prob, Genus == "g__Lachancea"))

#
## Module Taxonomy
mods_prob$count <- 1
mods_prob <- merge(mods_prob, tax_fun)
mods_prob[is.na(mods_prob)] <- "Unidentified"

mods_prob.m <- aggregate(mods_prob$count, list(mods_prob$modularity_class), sum)
colnames(mods_prob.m) <- c("modularity_class", "total")

#Genus
mods_prob.g <- aggregate(mods_prob$count, list(mods_prob$modularity_class,
                                                   mods_prob$Genus), sum)

colnames(mods_prob.g) <- c("modularity_class", "Genus", "count")
mods_prob.g <- merge(mods_prob.g, mods_prob.m, by = "modularity_class")
mods_prob.g$ab <- mods_prob.g$count/mods_prob.g$total

mods_prob.g$Genus[mods_prob.g$ab < 0.015] <- "Other"

orderG <- levels(factor(mods_prob.g$Genus))
orderG <- orderG[! orderG %in% c("Unidentified")]
orderG <- append(orderG, c("Unidentified"))

set.seed(123)
ggplot(mods_prob.g, aes(x = modularity_class, y = ab, 
                        fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Phylum", values = sample(col_vector, 18)) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(size = 19, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Sample") + ylab("Abundance")

#
#Family
mods_prob.f <- aggregate(mods_prob$count, list(mods_prob$modularity_class,
                                               mods_prob$Family), sum)

colnames(mods_prob.f) <- c("modularity_class", "Family", "count")
mods_prob.f <- merge(mods_prob.f, mods_prob.m, by = "modularity_class")
mods_prob.f$ab <- mods_prob.f$count/mods_prob.f$total

mods_prob.f$Family[mods_prob.f$ab < 0.019] <- "Other"

orderF <- levels(factor(mods_prob.f$Family))
orderF <- orderF[! orderF %in% c("Unidentified")]
orderF <- append(orderF, c("Unidentified"))

set.seed(123)
ggplot(mods_prob.f, aes(x = modularity_class, y = ab, 
                        fill = factor(Family, levels = orderF))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Phylum", values = sample(col_vector, 14)) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(size = 19, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Sample") + ylab("Abundance")


#
#### NETWORK ANALYSIS: METAWEB - Genus ####
#Transform our dataset into presence/absence table
tasv_fun.g <- cbind.data.frame(t(tasv_fun), Id = row.names(t(tasv_fun)))
tasv_fun.g <- merge(tasv_fun.g, tax_fun[,-c(7,9)], by = "Id")

tasv_fun.ug <- tasv_fun.g[!complete.cases(tasv_fun.g),]
tasv_fun.g <- tasv_fun.g[complete.cases(tasv_fun.g),]

tasv_fun.ug$Genus <- apply(tasv_fun.ug, 1, function(r) { r[(which(is.na(r))[1]-1)] } )
tasv_fun.ug$Genus <- paste(tasv_fun.ug$Genus, "Unidentified", sep = ".")

tasv_fun.g <- rbind(tasv_fun.g, tasv_fun.ug)

tasv_fun.g <- aggregate(tasv_fun.g[,2:268], list(tasv_fun.g$Genus), sum)
colnames(tasv_fun.g[,2:268]) <- colnames(tasv_fun.g[,2:268])
row.names(tasv_fun.g) <- tasv_fun.g[,1]

tasv_fun.g <- t(tasv_fun.g[,-1])

otu_prob <- as.matrix(tasv_fun.g[row.names(tasv_fun.g) %in% meta.ns$Sample_ID, ])
otu_prob <- otu_prob[, colSums(otu_prob > 0) > 1 & colSums(otu_prob > 0.001) > 0]

otu_prob[otu_prob > 0] <- 1

#Calculate how many times each ASV/OTU pair co-occur (appear in the same sample)
bac_prob <- as.matrix(otu_prob)
occurs_prob <- crossprod((bac_prob)) 

#Calculate probabilitity of co-occurrence for each pair
res <- matrix(NA, ncol(bac_prob), ncol(bac_prob))

for (i in 1:nrow(occurs_prob)){
  for (j in 2:ncol(occurs_prob)){
    if (j <= i) next
    out <- 0
    cooc <- occurs_prob[i,j] - 1 
    for (k in 0:cooc){ 
      out <- out + Veech_calc(nrow(bac_prob), k, occurs_prob[j,j], occurs_prob[i,i]) 
    }
    res[i, j] <- out # out = probability of finding less than k cooccurences. 
  }
}

#Filter only significant pairs (probability >= 0.95)
pos_prob <- subset(melt(res), value >= 0.95)
colnames(pos_prob) <- c("Source", "Target", "value")
pos_prob[, 1] <- colnames(bac_prob)[pos_prob[, 1]]
pos_prob[, 2] <- colnames(bac_prob)[pos_prob[, 2]]

#Construct node tables with taxonomic information of the nodes in the network
node_prob <- unique(rbind(cbind(Id = unique(as.character(pos_prob[,1]))), 
                          cbind(Id = unique(as.character(pos_prob[,2])))))

##DRAW NETWORK
#Import the egde list (list of significant co-occurring pair) to the igraph package
metanet_prob <- graph_from_data_frame(as.matrix(pos_prob[,c(1:2)]), directed = FALSE)

#Set the layout, the position in the plot qhere the nodes will be drawn, attending to their connections
set.seed(1)
l_prob <- layout_with_fr(metanet_prob, grid = "nogrid")
l_prob <- norm_coords(l_prob, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

#Calculate the number of modules (subcommunities) forming the network
cwt_prob <- cluster_walktrap(metanet_prob, steps = 20)
max(cwt_prob$membership)
V(metanet_prob)$community <- cwt_prob$membership
colrs <- adjustcolor( c("yellowgreen", "darkgoldenrod", "deepskyblue", "yellow"), alpha = 0.8)

#Draw the network
plot(metanet_prob, vertex.size = 5, vertex.label = NA, 
     vertex.color = colrs[as.numeric(V(metanet_prob)$community)], 
     rescale = FALSE, layout = l_prob*1.3, edge.color = "gray70")

#
##################################################################### RANDOM STUFF ####
#### Phylogenetic tree ####
tasv_t20 <- tasv_fun[subset(meta, Time == "Middle")$Sample_ID,]
tasv_t20 <- tasv_t20[, colSums(tasv_t20) > 0]

seqs_t20 <- getSequences(colnames(tasv_t20))
names(seqs_t20) <- seqs_t20

mult_t20 <- msa(seqs_t20, method = "ClustalW", type = "dna", order = "input")

phang.align_t20 <- as.phyDat(mult_t20, type = "DNA", names = seqs_t20)
dm_t20 <- dist.ml(phang.align_t20)
treeNJ_t20 <- NJ(dm_t20) # Note, tip order != sequence order
fit_t20 <- pml(treeNJ_t20, data = phang.align_t20)

fitGTR_t20 <- update(fit_t20, k = 4, inv = 0.2)
fitGTR_t20 <- optim.pml(fitGTR_t20, model = "GTR", optInv = TRUE, optGamma = TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0))

#
#### ALPHA DIVERISTY: t20 ####
alpha_t20 <- cbind.data.frame(# Hill based taxonomic alpha diversity
  t.q0 = hill_taxa(tasv_t20, q = 0),
  t.q1 = hill_taxa(tasv_t20, q = 1),
  t.q2 = hill_taxa(tasv_t20, q = 2),
  
  # Hill based phylogenetic alpha diversity
  p.q0 = hill_phylo(tasv_t20, fitGTR_t20$tree, q = 0),
  p.q1 = hill_phylo(tasv_t20, fitGTR_t20$tree, q = 1),
  p.q2 = hill_phylo(tasv_t20, fitGTR_t20$tree, q = 2))

alpha_t20$Sample_ID <- row.names(alpha_t20)

alpha_t20 <- merge(alpha_t20, subset(meta, order > 54)[,-c(7,8)], by = "Sample_ID")

alpha_t20.plot <- melt(alpha_t20)
alpha_t20.plot$diversity <- ifelse(startsWith(as.character(alpha_t20.plot$variable), "t."),
                                   "Taxonomic", "Phylogenetic")
alpha_t20.plot$diversity <- factor(alpha_t20.plot$diversity, 
                                   levels = c("Taxonomic", "Phylogenetic"))

alpha_t20.plot$q <- as.factor(substr(alpha_t20.plot$variable, 
                                     nchar(as.character(alpha_t20.plot$variable)),
                                     nchar(as.character(alpha_t20.plot$variable))))

alpha_t20.plot$name <- paste(alpha_t20.plot$Origin, alpha_t20.plot$Sample_name, 
                             alpha_t20.plot$q)

xx <- subset(alpha_t20.plot, diversity == "Taxonomic")
xx <- cbind.data.frame(xx, 
                       value2 = rbind.data.frame(cbind(subset(xx, q == 1)$value),
                                                 cbind(subset(xx, q == 1)$value),
                                                 cbind(subset(xx, q == 1)$value)))
xxx <- cbind.data.frame(xx[,-14], 
                        value2 = rbind.data.frame(cbind(subset(xx, q == 2)$value),
                                                  cbind(subset(xx, q == 2)$value),
                                                  cbind(subset(xx, q == 2)$value)))


ggplot(subset(alpha_t20.plot, diversity == "Taxonomic"), 
       aes(x = Time, y = value, shape = Management, color = Condition, group = name)) + 
  geom_point(position = position_dodge(width = 0.3), size = 2) + 
  geom_line(position = position_dodge(width = 0.3), size = 1) + 
  theme_bw() + 
  theme(legend.position = "bottom", legend.text.align = 0,
        axis.text.x = element_text(angle = 60, hjust = 1, size = 17, color = "black"),
        axis.text.y = element_text(size = 19, color = "black"),
        axis.title.x = element_text(size = 19, color = "black"),
        strip.text.x = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black")) +
  xlab("Stage") + ylab("Alpha-diversity") + facet_wrap(~q + Origin, nrow = 3)

#





##################################################################### METABOLITOS TFG RA?L ####
#### METABOLITOS ####
df_meta <- read.table("meta_df.txt", header = TRUE, sep = "\t")
df_meta <- merge(df_meta, meta[,c(6,10,12)], by = "Sample_ID")

df_meta$Time <- factor(df_meta$Time, levels = c("Initial", "Middle", "End"))
df_meta$Region <- factor(df_meta$Region, levels = c("RiberaGuadiana", "LaMancha", "Madrid", "Valdepe?as", "LaRioja"))

df_meta.t0 <- subset(df_meta, Time == "Initial")

df_meta <- merge(df_meta, df_meta.t0[,c(2:4,13:16)], by = c("Origin", "Farming", "Condition"))
colnames(df_meta)[13:16] <- c("Sugars", "Ammonia", "PAN", "Malic_acid")
colnames(df_meta)[27:30] <- c("Sugars.t0", "Ammonia.t0", "PAN.t0", "Malic_acid.t0")

df_meta$Sugars.perc <- (df_meta$Sugars/df_meta$Sugars.t0)*100
df_meta$Sugars.cons <- 100-df_meta$Sugars.perc
df_meta$Sugars.bcons <- (df_meta$Sugars.t0-df_meta$Sugars)

df_meta$Ammonia.perc <- (df_meta$Ammonia/df_meta$Ammonia.t0)*100
df_meta$Ammonia.cons <- 100-df_meta$Ammonia.perc
df_meta$Ammonia.bcons <- (df_meta$Ammonia.t0-df_meta$Ammonia)

df_meta$PAN.perc <- (df_meta$PAN/df_meta$PAN.t0)*100
df_meta$PAN.cons <- 100-df_meta$PAN.perc
df_meta$PAN.bcons <- (df_meta$PAN.t0-df_meta$PAN)

df_meta$Malic_acid.perc <- (df_meta$Malic_acid/df_meta$Malic_acid.t0)*100
df_meta$Malic_acid.cons <- 100-df_meta$Malic_acid.perc
df_meta$Malic_acid.bcons <- (df_meta$Malic_acid.t0-df_meta$Malic_acid)

ggplot(subset(df_meta[complete.cases(df_meta$most.ab),], Time != "Middle")) + 
  geom_boxplot(aes(x = Region, y = Ammonia, color = Time))

ggplot(subset(df_meta[complete.cases(df_meta$most.ab),], Time == "End")) + 
  geom_point(aes(x = Region, y = Ammonia.cons, color = Condition, shape = Farming))



df_meta.plot <- melt(df_meta)


ggplot(subset(df_meta.plot, variable == c("Sugars.cons", "Ammonia.cons", "PAN.cons") & Time == "End")) + 
  geom_boxplot(aes(x = most.ab, y = value, color = variable))

ggplot(subset(df_meta.plot, variable == c("Sugars.cons", "Ammonia.cons", "PAN.cons") & Time == "Initial")) + 
  geom_boxplot(aes(x = most.ab, y = value, color = variable))


df_meta.tf <- subset(df_meta, Time == "End")

ggplot(df_meta.tf, aes(x = log10(Ammonia.bcons/Sugars.bcons), y = log10(PAN.bcons/Sugars.bcons), color = most.ab)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)





ggplot(df_meta.tf, aes(x = (Ammonia.perc), y = (PAN.perc), color = most.ab)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)


cor.test(df_meta.tf$Ammonia, df_meta.tf$PAN)
cor.test(df_meta.tf$Ammonia/df_meta.tf$Sugars, df_meta.tf$PAN/df_meta.tf$Sugars)

cor.test(subset(df_meta.tf, most.ab != "g__Lachancea")$Ammonia.cons/subset(df_meta.tf, most.ab != "g__Lachancea")$Sugars.cons, 
         subset(df_meta.tf, most.ab != "g__Lachancea")$PAN.cons/subset(df_meta.tf, most.ab != "g__Lachancea")$Sugars.cons)
#







#### SAVE DATA ####
save.image("physeq_fun.RData")

#