# Code to plot Protea data in structural trait PC space and color points based on physiological performance


# load libraries
library(ggplot2)
library(cowplot)

# load Protea_data if not loaded
load("Data/Protea_data.RData")

# Define Performance and Structural Traits, as in Protea_data dataframe
Performance_Traits <- c("Ks_scaled", "LSC_scaled", "Photo_scaled", "Total_Assim_scaled",
                        "Cond_scaled", "WUE_Instan_scaled")
Structural_Traits <- c("BW_scaled", "WD_scaled", "Leaf_Area_scaled", "LMA_scaled",
                       "LD_scaled", "LWR_scaled", "Stom_L_scaled", "Stom_D_scaled")

# run PCA
pca_results <- princomp(Protea_data[, Structural_Traits], cor = TRUE)

# Get proprotion variance explained
summary(pca_results)

# PC1: 56.2%
# PC2: 20.7%


##########################################
# Plot with PC1 and PC2 as x- and y-axes #
##########################################

# Each plot will plot the points in the same position, but color based on the performance quntile

# 1. Hydraulic Conductance
color <- cut(Protea_data$Ks_scaled, breaks=quantile(Protea_data$Ks_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC1=pca_results$scores[, 1],
                    PC2=pca_results$scores[, 2],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_Ks <- ggplot(tr_pca_df, aes(x = PC1, y = PC2, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC1, y=PC2, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Sapwood-Specific \n Hydraulic Conductivity") + xlab("PC 1 (56.2%)") + ylab("PC 2 (20.7%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  annotate(geom = "text", x = 3.8, y = 2.4, label = "PRNI", size = 3) +
  annotate(geom = "text", x = 0.4, y = 3.3, label = "PRLA", size = 3) +
  annotate(geom = "text", x = -3.0, y = 0.9, label = "PRRE", size = 3) +
  annotate(geom = "text", x = 0.3, y = 0.3, label = "PRPU", size = 3) +
  annotate(geom = "text", x = 3.7, y = -0.9, label = "PREX", size = 3)
pc_Ks


# 2. Leaf Specific Conductivity
color <- cut(Protea_data$LSC_scaled, breaks=quantile(Protea_data$LSC_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC1=pca_results$scores[, 1],
                    PC2=pca_results$scores[, 2],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_LSC <- ggplot(tr_pca_df, aes(x = PC1, y = PC2, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC1, y=PC2, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Leaf-Specific Conductivity") + xlab("PC 1 (56.2%)") + ylab("PC 2 (20.7%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  annotate(geom = "text", x = 3.8, y = 2.4, label = "PRNI", size = 3) +
  annotate(geom = "text", x = 0.4, y = 3.3, label = "PRLA", size = 3) +
  annotate(geom = "text", x = -3.0, y = 0.9, label = "PRRE", size = 3) +
  annotate(geom = "text", x = 0.3, y = 0.3, label = "PRPU", size = 3) +
  annotate(geom = "text", x = 3.7, y = -0.9, label = "PREX", size = 3)
pc_LSC


# 3. Photosynthetic Rate
color <- cut(Protea_data$Photo_scaled, breaks=quantile(Protea_data$Photo_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC1=pca_results$scores[, 1],
                    PC2=pca_results$scores[, 2],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_Photo <- ggplot(tr_pca_df, aes(x = PC1, y = PC2, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC1, y=PC2, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Photosynthetic Rate") + xlab("PC 1 (56.2%)") + ylab("PC 2 (20.7%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  annotate(geom = "text", x = 3.8, y = 2.4, label = "PRNI", size = 3) +
  annotate(geom = "text", x = 0.4, y = 3.3, label = "PRLA", size = 3) +
  annotate(geom = "text", x = -3.0, y = 0.9, label = "PRRE", size = 3) +
  annotate(geom = "text", x = 0.3, y = 0.3, label = "PRPU", size = 3) +
  annotate(geom = "text", x = 3.7, y = -0.9, label = "PREX", size = 3)
pc_Photo


# 4. Total Assimilation
color <- cut(Protea_data$Total_Assim_scaled, breaks=quantile(Protea_data$Total_Assim_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC1=pca_results$scores[, 1],
                    PC2=pca_results$scores[, 2],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_Total_Assim <- ggplot(tr_pca_df, aes(x = PC1, y = PC2, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC1, y=PC2, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Leaf-Specific \n Photosynthetic Rate") + xlab("PC 1 (56.2%)") + ylab("PC 2 (20.7%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  annotate(geom = "text", x = 3.8, y = 2.4, label = "PRNI", size = 3) +
  annotate(geom = "text", x = 0.4, y = 3.3, label = "PRLA", size = 3) +
  annotate(geom = "text", x = -3.0, y = 0.9, label = "PRRE", size = 3) +
  annotate(geom = "text", x = 0.3, y = 0.3, label = "PRPU", size = 3) +
  annotate(geom = "text", x = 3.7, y = -0.9, label = "PREX", size = 3)
pc_Total_Assim


# 5. Stomatal Conductance
color <- cut(Protea_data$Cond_scaled, breaks=quantile(Protea_data$Cond_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC1=pca_results$scores[, 1],
                    PC2=pca_results$scores[, 2],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_Cond <- ggplot(tr_pca_df, aes(x = PC1, y = PC2, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC1, y=PC2, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Stomatal Conductance") + xlab("PC 1 (56.2%)") + ylab("PC 2 (20.7%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  annotate(geom = "text", x = 3.8, y = 2.4, label = "PRNI", size = 3) +
  annotate(geom = "text", x = 0.4, y = 3.3, label = "PRLA", size = 3) +
  annotate(geom = "text", x = -3.0, y = 0.9, label = "PRRE", size = 3) +
  annotate(geom = "text", x = 0.3, y = 0.3, label = "PRPU", size = 3) +
  annotate(geom = "text", x = 3.7, y = -0.9, label = "PREX", size = 3)
pc_Cond


# 6. Instantaneous Water Use Efficiency
color <- cut(Protea_data$WUE_Instan_scaled, breaks=quantile(Protea_data$WUE_Instan_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC1=pca_results$scores[, 1],
                    PC2=pca_results$scores[, 2],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_WUE <- ggplot(tr_pca_df, aes(x = PC1, y = PC2, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC1, y=PC2, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Instantaneous \n Water-Use Efficiency") + xlab("PC 1 (56.2%)") + ylab("PC 2 (20.7%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  annotate(geom = "text", x = 3.8, y = 2.4, label = "PRNI", size = 3) +
  annotate(geom = "text", x = 0.4, y = 3.3, label = "PRLA", size = 3) +
  annotate(geom = "text", x = -3.0, y = 0.9, label = "PRRE", size = 3) +
  annotate(geom = "text", x = 0.3, y = 0.3, label = "PRPU", size = 3) +
  annotate(geom = "text", x = 3.7, y = -0.9, label = "PREX", size = 3)
pc_WUE


plot_grid(pc_Ks, pc_LSC, pc_Photo, pc_Total_Assim, pc_Cond, pc_WUE,
          labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)


########################################################
# Plot using PC axes most aligned with functional axis #
########################################################

# Using three physiological performance traits for manuscript figure

# 1. Ks: PC4 and PC5 are most aligned
color <- cut(Protea_data$Ks_scaled, breaks=quantile(Protea_data$Ks_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC4=pca_results$scores[, 4],
                    PC5=pca_results$scores[, 5],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_Ks2 <- ggplot(tr_pca_df, aes(x = PC4, y = PC5, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC4, y=PC5, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Sapwood-Specific \n Hydraulic Conductivity") + xlab("PC 4 (6.3%)") + ylab("PC 5 (3.9%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
pc_Ks2


# 3. Photosynthetic Rate: PC5 and PC6 are most aligned
color <- cut(Protea_data$Photo_scaled, breaks=quantile(Protea_data$Photo_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC5=pca_results$scores[, 5],
                    PC6=pca_results$scores[, 6],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_Photo2 <- ggplot(tr_pca_df, aes(x = PC5, y = PC6, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC5, y=PC6, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Photosynthetic Rate") + xlab("PC 5 (3.9%)") + ylab("PC 6 (2.3%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
pc_Photo2


# 5. Stomatal Conductance: PC5 and PC8 are most aligned
color <- cut(Protea_data$Cond_scaled, breaks=quantile(Protea_data$Cond_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC5=pca_results$scores[, 5],
                    PC8=pca_results$scores[, 8],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE),
                    Phys_level_cont = Protea_data$Cond_scaled)

pc_Cond2 <- ggplot(tr_pca_df, aes(x = PC5, y = PC8, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC5, y=PC8, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Stomatal Conductance") + xlab("PC 5 (3.9%)") + ylab("PC 8 (1.0%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
pc_Cond2


# Multi-panel plot for three physiological performance PCs (top is PC1 and PC2, bottom is PC axes aligned w/ function)
plot_grid(pc_Ks, pc_Photo, pc_Cond, pc_Ks2, pc_Photo2, pc_Cond2,
          labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)

ggsave(filename = "PCA_Performance_panel.pdf", 
       path = "Output/Figures",
       height = 8, width = 12)


#########################################################
# Same as above, but for remaining 3 Performance Traits #
#########################################################

# 2. LSC: PC2 and PC8 are most aligned
color <- cut(Protea_data$LSC_scaled, breaks=quantile(Protea_data$LSC_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC2=pca_results$scores[, 2],
                    PC8=pca_results$scores[, 8],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_LSC2 <- ggplot(tr_pca_df, aes(x = PC2, y = PC8, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC2, y=PC8, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Leaf-Specific Conductivity") + xlab("PC 2 (20.7%)") + ylab("PC 8 (1.0%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
pc_LSC2


# 4. Total Assimilation: PC2 and PC7 are most aligned
color <- cut(Protea_data$Total_Assim_scaled, breaks=quantile(Protea_data$Total_Assim_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC2=pca_results$scores[, 2],
                    PC7=pca_results$scores[, 7],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_Total_Assim2 <- ggplot(tr_pca_df, aes(x = PC2, y = PC7, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC2, y=PC7, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Leaf-Specific \n Photosynthetic Rate") + xlab("PC 2 (20.7%)") + ylab("PC 7 (1.5%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
pc_Total_Assim2


# 6. Instantaneous WUE: PC6 and PC8 are most aligned
color <- cut(Protea_data$WUE_Instan_scaled, breaks=quantile(Protea_data$WUE_Instan_scaled, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
             c("Low",
               "Medium low",
               "Medium",
               "Medium high",
               "High"),
             include.lowest=TRUE)
tr_pca_df <- tibble(PC6=pca_results$scores[, 6],
                    PC8=pca_results$scores[, 8],
                    Species=Protea_data$Species,
                    Phys_level = factor(color, ordered = TRUE))

pc_WUE2 <- ggplot(tr_pca_df, aes(x = PC6, y = PC8, color = Phys_level)) +
  geom_point() + stat_ellipse(aes(x=PC6, y=PC8, group=Species), colour = "black") +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  ggtitle("Instantaneous \n Water-Use Efficiency") + xlab("PC 6 (2.3%)") + ylab("PC 8 (1.0%)") +
  theme_classic() +
  guides(color = "none") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))
pc_WUE2


# Multi-panel plot for the other three physiological performance PCs (top is PC1 and PC2, bottom is PC axes aligned w/ function)
plot_grid(pc_LSC, pc_Total_Assim, pc_WUE, pc_LSC2, pc_Total_Assim2, pc_WUE2,
          labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)

ggsave(filename = "PCA_Performance_panel_supplement.pdf", 
       path = "Output/Figures",
       height = 8, width = 12)



