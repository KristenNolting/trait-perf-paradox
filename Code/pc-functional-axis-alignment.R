# Calculating the alignment between each PC axis and the functional axis, 
# for all six physiological performance traits. 


## 1. Run the models to get regression vector
## 2. Pull out regression vector from models
## 3. Calculate consistency between regression vectors (1-6) and PC axes
## 4. Plot

# libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(brms)
library(tidybayes)
library(bayesplot)

# to optimize model run efficiency
options(mc.cores = parallel::detectCores())

# load Protea_data if not loaded
load("Data/Protea_data.RData")

# Define the two functions we'll use to get the cosine between a PC and functional axis
sumsq <- function(x) {
  return(sum(x^2))
}

cosine <- function(x, y) {
  numerator <- x %*% y
  denominator <- sqrt(sumsq(x))*sqrt(sumsq(y))
  return(numerator/denominator)
}

##########################################
# Step 1: Run multiple regression models #
##########################################

# Re-run multiple regressions, but exclude the grouping factors (we know these are importantant)
# but the generalized model and the PCA plot of Protea data do not reflect this information

# 1. Ks
ks_mod_2 <- brm(Ks_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                  LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled,
                data = Protea_data,
                family = gaussian(link = "identity"),
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(normal(0, 1), class = b)),
                control = list(adapt_delta = 0.999, max_treedepth = 12),
                save_pars = save_pars(all = TRUE),
                iter = 4000, warmup = 2000, chains = 4,
                file = "Output/Model_Fits/ks_mod_no_REs")


# 2. LSC
LSC_mod_2 <- brm(LSC_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                   LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled,
                 data = Protea_data,
                 family = gaussian(link = "identity"),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 1), class = b)),
                 control = list(adapt_delta = 0.999, max_treedepth = 12),
                 save_pars = save_pars(all = TRUE),
                 iter = 4000, warmup = 2000, chains = 4,
                 file = "Output/Model_Fits/LSC_mod_no_REs")


# 3. Photo
Photo_mod_2 <- brm(Photo_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                     LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled,
                   data = Protea_data,
                   family = gaussian(link = "identity"),
                   prior = c(prior(normal(0, 1), class = Intercept),
                             prior(normal(0, 1), class = b)),
                   control = list(adapt_delta = 0.999, max_treedepth = 12),
                   save_pars = save_pars(all = TRUE),
                   iter = 4000, warmup = 2000, chains = 4,
                   file = "Output/Model_Fits/Photo_mod_no_REs")


# 4. LSP
LSP_mod_2 <- brm(Total_Assim_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                   LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled,
                 data = Protea_data,
                 family = gaussian(link = "identity"),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 1), class = b)),
                 control = list(adapt_delta = 0.999, max_treedepth = 12),
                 save_pars = save_pars(all = TRUE),
                 iter = 4000, warmup = 2000, chains = 4,
                 file = "Output/Model_Fits/LSP_mod_no_REs")


# 5. Conductance
Cond_mod_2 <- brm(Cond_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                    LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled,
                  data = Protea_data,
                  family = gaussian(link = "identity"),
                  prior = c(prior(normal(0, 1), class = Intercept),
                            prior(normal(0, 1), class = b)),
                  control = list(adapt_delta = 0.999, max_treedepth = 12),
                  save_pars = save_pars(all = TRUE),
                  iter = 4000, warmup = 2000, chains = 4,
                  file = "Output/Model_Fits/Cond_mod_no_REs")


# 6. Instantaneous WUE
WUE_Instan_mod_2 <- brm(WUE_Instan_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                          LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled,
                        data = Protea_data,
                        family = gaussian(link = "identity"),
                        prior = c(prior(normal(0, 1), class = Intercept),
                                  prior(normal(0, 1), class = b)),
                        control = list(adapt_delta = 0.999, max_treedepth = 12),
                        save_pars = save_pars(all = TRUE),
                        iter = 4000, warmup = 2000, chains = 4,
                        file = "Output/Model_Fits/WUE_mod_no_REs")


#################################################
# Run steps 2-4 one performance trait at a time #
#################################################

# Step 2: Pull out model regression vectors
# Step 3: Run Consistency Analysis
# Step 4: Plot

# Make a list of the trait names for the PCA

Structural_Traits <- c("BW_scaled", "WD_scaled", "Leaf_Area_scaled", "LMA_scaled",
                       "LD_scaled", "LWR_scaled", "Stom_L_scaled", "Stom_D_scaled")

######
# Ks #
######

# Pull regression vector
ks2_coefs <- as.data.frame(fixef(ks_mod_2, summary = FALSE))

# Get PCs
dat_for_pca <- Protea_data %>%
  select("Species", all_of(Structural_Traits), 
         "Ks_scaled") %>%
  droplevels()
dat_for_pca <- na.omit(dat_for_pca)
pca_results <- princomp(dat_for_pca[, Structural_Traits], cor = TRUE)

# Distribution of cosines
n_iter <- nrow(ks2_coefs)
cosine_dist <- data.frame(PC1 = numeric(),
                          PC2 = numeric(),
                          PC3 = numeric(),
                          PC4 = numeric(),
                          PC5 = numeric(),
                          PC6 = numeric(),
                          PC7 = numeric(),
                          PC8 = numeric())

cosines <- numeric(ncol(pca_results$loadings))
for(j in 1:n_iter){
  mod_coefs <- as.numeric(ks2_coefs[j,])
  for (i in 1:ncol(pca_results$loadings)) {
    cosines[i] <- cosine(pca_results$loadings[, i], mod_coefs[-1])
  }
  cosines_df <- as.data.frame(cosines)
  cosines_df <- t(cosines_df)
  colnames(cosines_df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  cosine_dist <- rbind(cosine_dist, cosines_df)
}

# Convert to radians using arc cosine
CosDegrees_dist <- (acos(cosine_dist))*(180/pi)

# Convert Degrees to a scale between 0-90 degrees
CosDegrees_Abs_dist <- (acos(abs(cosine_dist)))*(180/pi)

# Plot 
Ks2_hist <-
  CosDegrees_Abs_dist %>%
  pivot_longer(cols = everything(), names_to = "PC_Axis", values_to = "Cosine_Value") %>%
  ggplot(aes(x = Cosine_Value, y = reorder(PC_Axis, desc(PC_Axis)))) +
  xlim(0, 90) +
  stat_histinterval(.width = c(.95, .5)) +
  ggtitle("Sapwood-Specific Hydraulic Conductivity") + 
  xlab("Angle Between Trait PC \n & Functional Axis") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank())
Ks2_hist


#######
# LSC #
#######

# Pull regression vector
LSC2_coefs <- as.data.frame(fixef(LSC_mod_2, summary = FALSE))

# Get PCs
dat_for_pca <- Protea_data %>%
  select("Species", all_of(Structural_Traits), 
         "LSC_scaled") %>%
  droplevels()
dat_for_pca <- na.omit(dat_for_pca)
pca_results <- princomp(dat_for_pca[, Structural_Traits], cor = TRUE)

# Distribution of cosines
n_iter <- nrow(LSC2_coefs)
cosine_dist <- data.frame(PC1 = numeric(),
                          PC2 = numeric(),
                          PC3 = numeric(),
                          PC4 = numeric(),
                          PC5 = numeric(),
                          PC6 = numeric(),
                          PC7 = numeric(),
                          PC8 = numeric())

cosines <- numeric(ncol(pca_results$loadings))
for(j in 1:n_iter){
  mod_coefs <- as.numeric(LSC2_coefs[j,])
  for (i in 1:ncol(pca_results$loadings)) {
    cosines[i] <- cosine(pca_results$loadings[, i], mod_coefs[-1])
  }
  cosines_df <- as.data.frame(cosines)
  cosines_df <- t(cosines_df)
  colnames(cosines_df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  cosine_dist <- rbind(cosine_dist, cosines_df)
}

# Convert to radians using arc cosine
CosDegrees_dist <- (acos(cosine_dist))*(180/pi)

# Convert Degrees to a scale between 0-90 degrees
CosDegrees_Abs_dist <- (acos(abs(cosine_dist)))*(180/pi)

# Plot
LSC2_hist <-
  CosDegrees_Abs_dist %>%
  pivot_longer(cols = everything(), names_to = "PC_Axis", values_to = "Cosine_Value") %>%
  ggplot(aes(x = Cosine_Value, y = reorder(PC_Axis, desc(PC_Axis)))) +
  xlim(0, 90) +
  stat_histinterval(.width = c(.95, .5)) +
  ggtitle("Leaf-Specific Conductivity") +
  xlab("Angle Between Trait PC \n & Functional Axis") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank())
LSC2_hist


#########
# Photo #
#########

# Pull regression vector
Photo2_coefs <- as.data.frame(fixef(Photo_mod_2, summary = FALSE))

# Get PCs
dat_for_pca <- Protea_data %>%
  select("Species", all_of(Structural_Traits), 
         "Photo_scaled") %>%
  droplevels()
dat_for_pca <- na.omit(dat_for_pca)
pca_results <- princomp(dat_for_pca[, Structural_Traits], cor = TRUE)

# Distribution of cosines
n_iter <- nrow(Photo2_coefs)
cosine_dist <- data.frame(PC1 = numeric(),
                          PC2 = numeric(),
                          PC3 = numeric(),
                          PC4 = numeric(),
                          PC5 = numeric(),
                          PC6 = numeric(),
                          PC7 = numeric(),
                          PC8 = numeric())

cosines <- numeric(ncol(pca_results$loadings))
for(j in 1:n_iter){
  mod_coefs <- as.numeric(Photo2_coefs[j,])
  for (i in 1:ncol(pca_results$loadings)) {
    cosines[i] <- cosine(pca_results$loadings[, i], mod_coefs[-1])
  }
  cosines_df <- as.data.frame(cosines)
  cosines_df <- t(cosines_df)
  colnames(cosines_df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  cosine_dist <- rbind(cosine_dist, cosines_df)
}

# Convert to radians using arc cosine
CosDegrees_dist <- (acos(cosine_dist))*(180/pi)

# Convert Degrees to a scale between 0-90 degrees
CosDegrees_Abs_dist <- (acos(abs(cosine_dist)))*(180/pi)

# Plot 
Photo2_hist <-
  CosDegrees_Abs_dist %>%
  pivot_longer(cols = everything(), names_to = "PC_Axis", values_to = "Cosine_Value") %>%
  ggplot(aes(x = Cosine_Value, y = reorder(PC_Axis, desc(PC_Axis)))) +
  xlim(0, 90) +
  stat_histinterval(.width = c(.95, .5)) +
  ggtitle("Photosynthetic Rate") +
  xlab("Angle Between Trait PC \n & Functional Axis") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank())
Photo2_hist


#######
# LSP #
#######

# Pull regression vector
Total_Assim2_coefs <- as.data.frame(fixef(LSP_mod_2, summary = FALSE))

# Get PCs
dat_for_pca <- Protea_data %>%
  select("Species", all_of(Structural_Traits), 
         "Total_Assim_scaled") %>%
  droplevels()
dat_for_pca <- na.omit(dat_for_pca)
pca_results <- princomp(dat_for_pca[, Structural_Traits], cor = TRUE)

# Distribution of cosines
n_iter <- nrow(Total_Assim2_coefs)
cosine_dist <- data.frame(PC1 = numeric(),
                          PC2 = numeric(),
                          PC3 = numeric(),
                          PC4 = numeric(),
                          PC5 = numeric(),
                          PC6 = numeric(),
                          PC7 = numeric(),
                          PC8 = numeric())

cosines <- numeric(ncol(pca_results$loadings))
for(j in 1:n_iter){
  mod_coefs <- as.numeric(Total_Assim2_coefs[j,])
  for (i in 1:ncol(pca_results$loadings)) {
    cosines[i] <- cosine(pca_results$loadings[, i], mod_coefs[-1])
  }
  cosines_df <- as.data.frame(cosines)
  cosines_df <- t(cosines_df)
  colnames(cosines_df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  cosine_dist <- rbind(cosine_dist, cosines_df)
}

# Convert to radians using arc cosine
CosDegrees_dist <- (acos(cosine_dist))*(180/pi)

# Convert Degrees to a scale between 0-90 degrees
CosDegrees_Abs_dist <- (acos(abs(cosine_dist)))*(180/pi)

# Plot 
Total_Assim2_hist <-
  CosDegrees_Abs_dist %>%
  pivot_longer(cols = everything(), names_to = "PC_Axis", values_to = "Cosine_Value") %>%
  ggplot(aes(x = Cosine_Value, y = reorder(PC_Axis, desc(PC_Axis)))) +
  xlim(0, 90) +
  stat_histinterval(.width = c(.95, .5)) +
  ggtitle("Leaf-Specific Photosynthetic Rate") +
  xlab("Angle Between Trait PC \n & Functional Axis") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank())
Total_Assim2_hist


########
# Cond #
########

# Pull regression vector
Cond2_coefs <- as.data.frame(fixef(Cond_mod_2, summary = FALSE))

# Get PCs
dat_for_pca <- Protea_data %>%
  select("Species", all_of(Structural_Traits), 
         "Cond_scaled") %>%
  droplevels()
dat_for_pca <- na.omit(dat_for_pca)
pca_results <- princomp(dat_for_pca[, Structural_Traits], cor = TRUE)

# Distribution of cosines
n_iter <- nrow(Cond2_coefs)
cosine_dist <- data.frame(PC1 = numeric(),
                          PC2 = numeric(),
                          PC3 = numeric(),
                          PC4 = numeric(),
                          PC5 = numeric(),
                          PC6 = numeric(),
                          PC7 = numeric(),
                          PC8 = numeric())

cosines <- numeric(ncol(pca_results$loadings))
for(j in 1:n_iter){
  mod_coefs <- as.numeric(Cond2_coefs[j,])
  for (i in 1:ncol(pca_results$loadings)) {
    cosines[i] <- cosine(pca_results$loadings[, i], mod_coefs[-1])
  }
  cosines_df <- as.data.frame(cosines)
  cosines_df <- t(cosines_df)
  colnames(cosines_df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  cosine_dist <- rbind(cosine_dist, cosines_df)
}

# Convert to radians using arc cosine
CosDegrees_dist <- (acos(cosine_dist))*(180/pi)

# Convert Degrees to a scale between 0-90 degrees
CosDegrees_Abs_dist <- (acos(abs(cosine_dist)))*(180/pi)

# Plot 
Cond2_hist <-
  CosDegrees_Abs_dist %>%
  pivot_longer(cols = everything(), names_to = "PC_Axis", values_to = "Cosine_Value") %>%
  ggplot(aes(x = Cosine_Value, y = reorder(PC_Axis, desc(PC_Axis)))) +
  xlim(0, 90) +
  stat_histinterval(.width = c(.95, .5)) +
  ggtitle("Stomatal Conductance") +
  xlab("Angle Between Trait PC \n & Functional Axis") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank())
Cond2_hist


##############
# WUE_Instan #
##############

# Pull regression vector
WUE_Instan2_coefs <- as.data.frame(fixef(WUE_Instan_mod_2, summary = FALSE))

# Get PCs
dat_for_pca <- Protea_data %>%
  select("Species", all_of(Structural_Traits), 
         "WUE_Instan_scaled") %>%
  droplevels()
dat_for_pca <- na.omit(dat_for_pca)
pca_results <- princomp(dat_for_pca[, Structural_Traits], cor = TRUE)

# Distribution of cosines
n_iter <- nrow(WUE_Instan2_coefs)
cosine_dist <- data.frame(PC1 = numeric(),
                          PC2 = numeric(),
                          PC3 = numeric(),
                          PC4 = numeric(),
                          PC5 = numeric(),
                          PC6 = numeric(),
                          PC7 = numeric(),
                          PC8 = numeric())

cosines <- numeric(ncol(pca_results$loadings))
for(j in 1:n_iter){
  mod_coefs <- as.numeric(WUE_Instan2_coefs[j,])
  for (i in 1:ncol(pca_results$loadings)) {
    cosines[i] <- cosine(pca_results$loadings[, i], mod_coefs[-1])
  }
  cosines_df <- as.data.frame(cosines)
  cosines_df <- t(cosines_df)
  colnames(cosines_df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  cosine_dist <- rbind(cosine_dist, cosines_df)
}

# Convert to radians using arc cosine
CosDegrees_dist <- (acos(cosine_dist))*(180/pi)

# Convert Degrees to a scale between 0-90 degrees
CosDegrees_Abs_dist <- (acos(abs(cosine_dist)))*(180/pi)

# Plot
WUE_Instan2_hist <-
  CosDegrees_Abs_dist %>%
  pivot_longer(cols = everything(), names_to = "PC_Axis", values_to = "Cosine_Value") %>%
  ggplot(aes(x = Cosine_Value, y = reorder(PC_Axis, desc(PC_Axis)))) +
  xlim(0, 90) +
  stat_histinterval(.width = c(.95, .5)) +
  ggtitle("Instantaneous Water-Use Efficiency") +
  xlab("Angle Between Trait PC \n & Functional Axis") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank())
WUE_Instan2_hist 


# Panel Plot
plot_grid(Ks2_hist, LSC2_hist, Photo2_hist, Total_Assim2_hist, Cond2_hist, WUE_Instan2_hist,
          labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)

# save
ggsave(filename = "consistency-analysis_panel-plot_wide.pdf", 
       path = "Output/Figures",
       height = 8, width = 12)
