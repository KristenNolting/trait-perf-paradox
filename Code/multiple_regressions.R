# Code to run multiple regressions demonstrating that combos of structural traits 
# predict variation in physiological performance traits, quite well. 

# Code to recreate multi-panel predicted vs. observed plots

# load libraries
library(rstan)
library(rstanarm)
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

# 1. Hydraulic Conductance

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
              file = "/Output/Model_Fits/ks_mod")

# Model Checks
summary(ks_mod, digits = 3, probs = c(0.025, 0.975))
brms::bayes_R2(ks_mod)[1] # 0.36, Consistent from path model in Annals paper
hist(brms::bayes_R2(ks_mod, summary = FALSE))
hist(brms::bayes_R2(ks_mod, summary = FALSE, re_formula = NA))
brms::bayes_R2(ks_mod, re_formula = NA)[1] # 0.20
pp_check(ks_mod, type = "dens_overlay")
pp_check(ks_mod, type = "intervals")
pp_check(ks_mod, type = "loo_pit")
ks_mod_loo <- loo(ks_mod, moment_match = TRUE)
plot(ks_mod_loo, label_points = TRUE) # Points 16 and 38 with high Pareto-k values
Protea_data[16,] # PRPU at JK
Protea_data[38,] # PRRE at JK


# Plot Observed versus Predicted
p_post <- predict(ks_mod) # get predicted ks values from original data and fitted model
p_post_Ks <- as.data.frame(p_post)
p_post_Ks$Obs <- Protea_data$Ks_scaled
p_post_Ks$Species <- Protea_data$Species
p_post_Ks$Site <- Protea_data$Site
ks_pp <- ggplot(data = p_post_Ks, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) +
  #scale_color_viridis(discrete = TRUE) +
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Sapwood-specific \n Hydraulic Conductance") + xlab("Predicted Ks (scaled)") + ylab("Observed Ks (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
ks_pp



# 2. Leaf Specific Conductivity
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
               file = "/Users/knolting/Documents/UConn_Projects/MW2PAC/Model_Output/LSC_mod")

# Model Checks
summary(LSC_mod, digits = 3, probs = c(0.025, 0.975))
brms::bayes_R2(LSC_mod)[1] # 0.30, Consistent from path model in Annals paper
hist(brms::bayes_R2(LSC_mod, summary = FALSE))
hist(brms::bayes_R2(LSC_mod, summary = FALSE, re_formula = NA))
brms::bayes_R2(LSC_mod, re_formula = NA)[1] # 0.22
pp_check(LSC_mod, type = "dens_overlay")
pp_check(LSC_mod, type = "intervals")
pp_check(LSC_mod, type = "loo_pit_qq")
LSC_mod_loo <- loo(LSC_mod, moment_match = TRUE) # not working - need to think about model fit
LSC_mod_loo
plot(LSC_mod_loo, label_points = TRUE)


# Plot Observed versus Predicted
p_post <- predict(LSC_mod) # get predicted ks values from original data and fitted model
p_post_LSC <- as.data.frame(p_post)
p_post_LSC$Obs <- Protea_data$LSC_scaled
p_post_LSC$Species <- Protea_data$Species
p_post_LSC$Site <- Protea_data$Site
LSC_pp <- ggplot(data = p_post_LSC, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) + 
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  #scale_color_viridis(discrete = TRUE) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Leaf-specific \n Conductivity") + xlab("Predicted LSC (scaled)") + ylab("Observed LSC (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
LSC_pp


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
                 file = "/Users/knolting/Documents/UConn_Projects/MW2PAC/Model_Output/Photo_mod")

# Model Checks
summary(Photo_mod, digits = 3, probs = c(0.025, 0.975))
brms::bayes_R2(Photo_mod)[1] # 0.44, Consistent from path model in Annals paper
hist(brms::bayes_R2(Photo_mod, summary = FALSE)) 
hist(brms::bayes_R2(Photo_mod, summary = FALSE, re_formula = NA))
brms::bayes_R2(Photo_mod, re_formula = NA)[1] # 0.27
pp_check(Photo_mod, type = "dens_overlay")
pp_check(Photo_mod, type = "intervals")
pp_check(Photo_mod, type = "loo_pit_qq")
Photo_mod_loo <- loo(Photo_mod, moment_match = TRUE)
Photo_mod_loo
plot(Photo_mod_loo, label_points = TRUE) # Points 70 and 126 with high Pareto-k values
Protea_data[c(70, 126),] # PRLA at CB and PRRE at SB


# Plot Observed versus Predicted
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
                       file = "/Users/knolting/Documents/UConn_Projects/MW2PAC/Model_Output/LSP_mod")

# Model Checks
summary(Total_Assim_mod, digits = 3, probs = c(0.025, 0.975))
brms::bayes_R2(Total_Assim_mod)[1] # 0.32, Consistent from path model in Annals paper
hist(brms::bayes_R2(Total_Assim_mod, summary = FALSE))
hist(brms::bayes_R2(Total_Assim_mod, summary = FALSE, re_formula = NA))
brms::bayes_R2(Total_Assim_mod, re_formula = NA)[1] # 0.36
pp_check(Total_Assim_mod, type = "dens_overlay")
pp_check(Total_Assim_mod, type = "intervals")
pp_check(Total_Assim_mod, type = "loo_pit_qq")
Total_Assim_mod_loo <- loo(Total_Assim_mod, moment_match = TRUE)
Total_Assim_mod_loo
plot(Total_Assim_mod_loo, label_points = TRUE)


# Plot Observed versus Predicted
p_post <- predict(Total_Assim_mod) # get predicted ks values from original data and fitted model
p_post_Total_Assim <- as.data.frame(p_post)
p_post_Total_Assim$Obs <- Protea_data$Total_Assim_scaled
p_post_Total_Assim$Species <- Protea_data$Species
p_post_Total_Assim$Site <- Protea_data$Site
Total_Assim_pp <- ggplot(data = p_post_Total_Assim, aes(x = Estimate, y = Obs)) +
  geom_point(size = 2, aes(shape = Site, colour = as.factor(Species))) + 
  scale_fill_distiller(type = "seq", palette = "Spectral", direction = 1) +
  #scale_color_viridis(discrete = TRUE) +
  geom_abline(slope = 1, colour = "grey44", linetype = "dashed") +
  ggtitle("Leaf-specific \n Photosynthetic Rate") + xlab("Predicted LSP (scaled)") + ylab("Observed LSP (scaled)") +
  theme_classic() + labs(colour = "Species") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none")
Total_Assim_pp


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
                file = "/Users/knolting/Documents/UConn_Projects/MW2PAC/Model_Output/Cond_mod")

# Model Checks
summary(Cond_mod, digits = 3, probs = c(0.025, 0.975))
brms::bayes_R2(Cond_mod)[1] # 0.37
hist(brms::bayes_R2(Cond_mod, summary = FALSE))
hist(brms::bayes_R2(Cond_mod, summary = FALSE, re_formula = NA))
brms::bayes_R2(Cond_mod, re_formula = NA)[1] # 0.18
pp_check(Cond_mod, type = "dens_overlay")
pp_check(Cond_mod, type = "intervals")
pp_check(Cond_mod, type = "loo_pit_qq")
Cond_mod_loo <- loo(Cond_mod, moment_match = TRUE)
Cond_mod_loo
plot(Cond_mod_loo, label_points = TRUE) # Points 42 and 51 with high Pareto-k values
Protea_data[c(42, 51),] # BOTH PRLA at JK


