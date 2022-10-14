# Code to run multiple regressions demonstrating that combos of structural traits 
# predict variation in physiological performance traits, quite well. 

# Code to recreate multi-panel predicted vs. observed plots

# load libraries
library(brms)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(tidybayes)

# optimize computing
options(mc.cores = parallel::detectCores())

# load Protea_data if not loaded
load("Data/Protea_data.RData")

# Run each of the six multiple regressions and pull R2 values

# 1. Stem-Specific Conductivity (Ks)

ks_mod <- brm(Ks_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled +
                (1|Site) + (1|Species) + (1|Site:Species),
              data = Protea_data,
              family = gaussian(link = "identity"),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(normal(0, 1), class = b),
                        prior(cauchy(0, 5), class = sd)),
              control = list(adapt_delta = 0.999, max_treedepth = 12),
              save_pars = save_pars(all = TRUE),
              iter = 4000, warmup = 2000, chains = 4,
              file = "Output/Model_Fits/ks_mod")

bayes_R2(ks_mod)
bayes_R2(ks_mod, re_formula = NA)
write.csv(print(summary(ks_mod), file = "Output/Model_Fits/ks_mod.txt"))


# 2. Leaf-Specific Conductivity (LSC)
LSC_mod <- brm(LSC_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                 LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled +
                 (1|Site) + (1|Species) + (1|Site:Species),
               data = Protea_data,
               family = gaussian(link = "identity"),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(normal(0, 1), class = b),
                         prior(cauchy(0, 5), class = sd)),
               control = list(adapt_delta = 0.999, max_treedepth = 12),
               save_pars = save_pars(all = TRUE),
               iter = 4000, warmup = 2000, chains = 4,
               file = "Output/Model_Fits/LSC_mod")

bayes_R2(LSC_mod)
bayes_R2(LSC_mod, re_formula = NA)


# 3. Area-based photosynthetic rate
Photo_mod <- brm(Photo_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                   LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled +
                   (1|Site) + (1|Species) + (1|Site:Species),
                 data = Protea_data,
                 family = gaussian(link = "identity"),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 1), class = b),
                           prior(cauchy(0, 5), class = sd)),
                 control = list(adapt_delta = 0.999, max_treedepth = 12),
                 save_pars = save_pars(all = TRUE),
                 iter = 4000, warmup = 2000, chains = 4,
                 file = "Output/Model_Fits/Photo_mod")

bayes_R2(photo_mod)
bayes_R2(photo_mod, re_formula = NA)


# 4. Total Assimilation (LSP)
Total_Assim_mod <- brm(Total_Assim_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                         LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled +
                         (1|Site) + (1|Species) + (1|Site:Species),
                       data = Protea_data,
                       family = gaussian(link = "identity"),
                       prior = c(prior(normal(0, 1), class = Intercept),
                                 prior(normal(0, 1), class = b),
                                 prior(cauchy(0, 5), class = sd)),
                       control = list(adapt_delta = 0.999, max_treedepth = 12),
                       save_pars = save_pars(all = TRUE),
                       iter = 4000, warmup = 2000, chains = 4,
                       file = "Output/Model_Fits/LSP_mod")

bayes_R2(Total_Assim_mod)
bayes_R2(Total_Assim_mod, re_formula = NA)


# 5. Stomatal Conductance
Cond_mod <- brm(Cond_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                  LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled +
                  (1|Site) + (1|Species) + (1|Site:Species),
                data = Protea_data,
                family = gaussian(link = "identity"),
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(normal(0, 1), class = b),
                          prior(cauchy(0, 5), class = sd)),
                control = list(adapt_delta = 0.999, max_treedepth = 12),
                save_pars = save_pars(all = TRUE),
                iter = 4000, warmup = 2000, chains = 4,
                file = "Output/Model_Fits/Cond_mod")

bayes_R2(Cond_mod)
bayes_R2(Cond_mod, re_formula = NA)


# 6. Instantaneous WUE
WUE_Instan_mod <- brm(WUE_Instan_scaled ~ BW_scaled + WD_scaled + Leaf_Area_scaled + 
                        LMA_scaled + LD_scaled + LWR_scaled + Stom_L_scaled + Stom_D_scaled +
                        (1|Site) + (1|Species) + (1|Site:Species),
                      data = Protea_data,
                      family = gaussian(link = "identity"),
                      prior = c(prior(normal(0, 1), class = Intercept),
                                prior(normal(0, 1), class = b),
                                prior(cauchy(0, 5), class = sd)),
                      control = list(adapt_delta = 0.999, max_treedepth = 12),
                      save_pars = save_pars(all = TRUE),
                      iter = 4000, warmup = 2000, chains = 4,
                      file = "Output/Model_Fits/WUE_mod")

bayes_R2(WUE_Instan_mod)
bayes_R2(WUE_Instan_mod, re_formula = NA)


# Make a predicted versus observed plot for each model, and combine into multipanel plot

# 1. Ks
p_post <- predict(ks_mod) # get predicted Ks values from original data and fitted model
p_post_Ks <- as.data.frame(p_post)
p_post_Ks$Obs <- Protea_data$Ks_scaled
p_post_Ks$Species <- Protea_data$Species
p_post_Ks$Site <- Protea_data$Site
Ks_pp <- ggplot(data = p_post_Ks, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Sapwood-specific \n Hydraulic Conductivity") + xlab("Predicted Ks (scaled)") + ylab("Observed Ks (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
Ks_pp


# 2. LSC
p_post <- predict(LSC_mod) # get predicted ks values from original data and fitted model
p_post_LSC <- as.data.frame(p_post)
p_post_LSC$Obs <- Protea_data$LSC_scaled
p_post_LSC$Species <- Protea_data$Species
p_post_LSC$Site <- Protea_data$Site
LSC_pp <- ggplot(data = p_post_LSC, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) + 
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Leaf-specific \n Conductivity") + xlab("Predicted LSC (scaled)") + ylab("Observed LSC (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
LSC_pp


# 3. Photo
p_post <- predict(Photo_mod) # get predicted ks values from original data and fitted model
p_post_Photo <- as.data.frame(p_post)
p_post_Photo$Obs <- Protea_data$Photo_scaled
p_post_Photo$Species <- Protea_data$Species
p_post_Photo$Site <- Protea_data$Site
Photo_pp <- ggplot(data = p_post_Photo, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  #scale_color_viridis(discrete = TRUE) +
  xlim(-2, 2) + ylim(-2.8, 2.8) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Photosynthetic \n Rate") + xlab("Predicted Aarea (scaled)") + ylab("Observed Aarea (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
Photo_pp


# 4. Total_Assim (LSP)
p_post <- predict(Total_Assim_mod) # get predicted ks values from original data and fitted model
p_post_Total_Assim <- as.data.frame(p_post)
p_post_Total_Assim$Obs <- Protea_data$Total_Assim_scaled
p_post_Total_Assim$Species <- Protea_data$Species
p_post_Total_Assim$Site <- Protea_data$Site
Total_Assim_pp <- ggplot(data = p_post_Total_Assim, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) + 
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Leaf-specific \n Photosynthetic Rate") + xlab("Predicted LSP (scaled)") + ylab("Observed LSP (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
Total_Assim_pp


# 5. Stomtal Conductance
p_post <- predict(Cond_mod) # get predicted ks values from original data and fitted model
p_post_Cond <- as.data.frame(p_post)
p_post_Cond$Obs <- Protea_data$Cond_scaled
p_post_Cond$Species <- Protea_data$Species
p_post_Cond$Site <- Protea_data$Site
Cond_pp <- ggplot(data = p_post_Cond, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Stomatal \n Conductance") + xlab("Predicted gs (scaled)") + ylab("Observed gs (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
Cond_pp


# 6. WUE_Instan
p_post <- predict(WUE_Instan_mod) # get predicted ks values from original data and fitted model
p_post_WUE_Instan <- as.data.frame(p_post)
p_post_WUE_Instan$Obs <- Protea_data$WUE_Instan_scaled
p_post_WUE_Instan$Species <- Protea_data$Species
p_post_WUE_Instan$Site <- Protea_data$Site
WUE_Instan_pp <- ggplot(data = p_post_WUE_Instan, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) + 
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Instantaneous \n Water-use Efficiency") + xlab("Predicted WUE_Instan (scaled)") + ylab("Observed WUE_Instan (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
WUE_Instan_pp


# Panel Plot
plot_grid(Ks_pp, LSC_pp, Photo_pp, Total_Assim_pp, Cond_pp, WUE_Instan_pp,
          labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)

# save
ggsave(filename = "trait-performance-regression_panel-plot.pdf", 
       path = "Output/Figures",
       height = 6, width = 10)
