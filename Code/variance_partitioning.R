# Code to perform variance partitioning and plotting. 

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

# For each trait, run an intercept only model, with species as a random intercept
# Among species variance is: (variance of the intercepts)/(variance of the intercepts + square of sigma)

#####################
# Structural Traits #
#####################

# 1. Bark Width (BW)
BW_var <- brm(BW_scaled ~ 1 + (1|Species),
              data = Protea_data,
              family = gaussian(link = "identity"),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(cauchy(0, 5), class = sd)),
              control = list(adapt_delta = 0.99, max_treedepth = 12),
              save_pars = save_pars(all = TRUE),
              iter = 4000, warmup = 2000, chains = 4,
              file = "Output/Model_Fits/BW_var")

BW_var_df <- as.data.frame(BW_var)
BW_among_var <- apply(BW_var_df[,5:9], 1, var)
BW_within_var <- (BW_var_df$sigma)^2
BW_among_var_percent <- (BW_among_var/(BW_among_var + BW_within_var))*100


# 2. Wood Density (WD)
WD_var <- brm(WD_scaled ~ 1 + (1|Species),
              data = Protea_data,
              family = gaussian(link = "identity"),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(cauchy(0, 5), class = sd)),
              control = list(adapt_delta = 0.99, max_treedepth = 12),
              save_pars = save_pars(all = TRUE),
              iter = 4000, warmup = 2000, chains = 4,
              file = "Output/Model_Fits/WD_var")

WD_var_df <- as.data.frame(WD_var)
WD_among_var <- apply(WD_var_df[,5:9], 1, var)
WD_within_var <- (WD_var_df$sigma)^2
WD_among_var_percent <- (WD_among_var/(WD_among_var + WD_within_var))*100


# 3. Leaf Area (Leaf_Area)
LA_var <- brm(Leaf_Area_scaled ~ 1 + (1|Species),
              data = Protea_data,
              family = gaussian(link = "identity"),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(cauchy(0, 5), class = sd)),
              control = list(adapt_delta = 0.99, max_treedepth = 12),
              save_pars = save_pars(all = TRUE),
              iter = 4000, warmup = 2000, chains = 4,
              file = "Output/Model_Fits/LA_var")

LA_var_df <- as.data.frame(LA_var)
LA_among_var <- apply(LA_var_df[,5:9], 1, var)
LA_within_var <- (LA_var_df$sigma)^2
LA_among_var_percent <- (LA_among_var/(LA_among_var + LA_within_var))*100


# 4. Leaf Mass per Area (LMA)
LMA_var <- brm(LMA_scaled ~ 1 + (1|Species),
               data = Protea_data,
               family = gaussian(link = "identity"),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(cauchy(0, 5), class = sd)),
               control = list(adapt_delta = 0.99, max_treedepth = 12),
               save_pars = save_pars(all = TRUE),
               iter = 4000, warmup = 2000, chains = 4,
               file = "Output/Model_Fits/LMA_var")

LMA_var_df <- as.data.frame(LMA_var)
LMA_among_var <- apply(LMA_var_df[,5:9], 1, var)
LMA_within_var <- (LMA_var_df$sigma)^2
LMA_among_var_percent <- (LMA_among_var/(LMA_among_var + LMA_within_var))*100


# 5. Lamina Density (LD)
LD_var <- brm(LD_scaled ~ 1 + (1|Species),
              data = Protea_data,
              family = gaussian(link = "identity"),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(cauchy(0, 5), class = sd)),
              control = list(adapt_delta = 0.99, max_treedepth = 12),
              save_pars = save_pars(all = TRUE),
              iter = 4000, warmup = 2000, chains = 4,
              file = "Output/Model_Fits/LD_var")

LD_var_df <- as.data.frame(LD_var)
LD_among_var <- apply(LD_var_df[,5:9], 1, var)
LD_within_var <- (LD_var_df$sigma)^2
LD_among_var_percent <- (LD_among_var/(LD_among_var + LD_within_var))*100


# 6. Leaf LWR (LWR)
LWR_var <- brm(LWR_scaled ~ 1 + (1|Species),
               data = Protea_data,
               family = gaussian(link = "identity"),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(cauchy(0, 5), class = sd)),
               control = list(adapt_delta = 0.99, max_treedepth = 12),
               save_pars = save_pars(all = TRUE),
               iter = 4000, warmup = 2000, chains = 4,
               file = "Output/Model_Fits/LWR_var")

LWR_var_df <- as.data.frame(LWR_var)
LWR_among_var <- apply(LWR_var_df[,5:9], 1, var)
LWR_within_var <- (LWR_var_df$sigma)^2
LWR_among_var_percent <- (LWR_among_var/(LWR_among_var + LWR_within_var))*100


# 7. Stomatal Length (StomL)
StomL_var <- brm(Stom_L_scaled ~ 1 + (1|Species),
                 data = Protea_data,
                 family = gaussian(link = "identity"),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(cauchy(0, 5), class = sd)),
                 control = list(adapt_delta = 0.99, max_treedepth = 12),
                 save_pars = save_pars(all = TRUE),
                 iter = 4000, warmup = 2000, chains = 4,
                 file = "Output/Model_Fits/StomL_var")

StomL_var_df <- as.data.frame(StomL_var)
StomL_among_var <- apply(StomL_var_df[,5:9], 1, var)
StomL_within_var <- (StomL_var_df$sigma)^2
StomL_among_var_percent <- (StomL_among_var/(StomL_among_var + StomL_within_var))*100


# 8. Stomatal Density (StomD)
StomD_var <- brm(Stom_D_scaled ~ 1 + (1|Species),
                 data = Protea_data,
                 family = gaussian(link = "identity"),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(cauchy(0, 5), class = sd)),
                 control = list(adapt_delta = 0.99, max_treedepth = 12),
                 save_pars = save_pars(all = TRUE),
                 iter = 4000, warmup = 2000, chains = 4,
                 file = "Output/Model_Fits/StomD_var")

StomD_var_df <- as.data.frame(StomD_var)
StomD_among_var <- apply(StomD_var_df[,5:9], 1, var)
StomD_within_var <- (StomD_var_df$sigma)^2
StomD_among_var_percent <- (StomD_among_var/(StomD_among_var + StomD_within_var))*100


# Make a dataframe for plotting
Str_Trait_Var <- cbind(BW_among_var_percent, WD_among_var_percent,
                       LA_among_var_percent, LMA_among_var_percent,
                       LD_among_var_percent, LWR_among_var_percent,
                       StomL_among_var_percent, StomD_among_var_percent)
Str_Trait_Var <- as.data.frame(Str_Trait_Var)
Str_Trait_Var_df <- pivot_longer(data = Str_Trait_Var,
                                 cols = everything(),
                                 names_to = "Structural_Trait",
                                 values_to = "Among_Var")

Str_Trait_plot2 <- ggplot(data = Str_Trait_Var_df, 
                          aes(y = Structural_Trait, x = Among_Var)) +
  xlim(0,100) +
  stat_pointinterval(.width = c(0.95, 0.8)) +
  ggtitle("Structural Traits") + xlab("Among Species Variance (%)") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank()) +
  scale_y_discrete(labels=c("WD_among_var_percent" = "Wood Density", 
                            "StomL_among_var_percent" = "Stomatal Length",
                            "StomD_among_var_percent" = "Stomatal Density",
                            "LWR_among_var_percent" = "LWR",
                            "LMA_among_var_percent" = "LMA",
                            "LD_among_var_percent" = "Lamina Density",
                            "LA_among_var_percent" = "Leaf Area",
                            "BW_among_var_percent" = "Bark Width"))
Str_Trait_plot2


######################
# Performance Traits #
######################

# 1. Hydraulic Conductance (Ks)
Ks_var <- brm(Ks_scaled ~ 1 + (1|Species),
              data = Protea_data,
              family = gaussian(link = "identity"),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(cauchy(0, 5), class = sd)),
              control = list(adapt_delta = 0.99, max_treedepth = 12),
              save_pars = save_pars(all = TRUE),
              iter = 4000, warmup = 2000, chains = 4,
              file = "Output/Model_Fits/Ks_var")

Ks_var_df <- as.data.frame(Ks_var)
Ks_among_var <- apply(Ks_var_df[,5:9], 1, var)
Ks_within_var <- (Ks_var_df$sigma)^2
Ks_among_var_percent <- (Ks_among_var/(Ks_among_var + Ks_within_var))*100


# 2. Leaf Specific Conductivity (LSC)
LSC_var <- brm(LSC_scaled ~ 1 + (1|Species),
               data = Protea_data,
               family = gaussian(link = "identity"),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(cauchy(0, 5), class = sd)),
               control = list(adapt_delta = 0.99, max_treedepth = 12),
               save_pars = save_pars(all = TRUE),
               iter = 4000, warmup = 2000, chains = 4,
               file = "Output/Model_Fits/LSC_var")

LSC_var_df <- as.data.frame(LSC_var)
LSC_among_var <- apply(LSC_var_df[,5:9], 1, var)
LSC_within_var <- (LSC_var_df$sigma)^2
LSC_among_var_percent <- (LSC_among_var/(LSC_among_var + LSC_within_var))*100


# 3. Photosynthetic Rate (Photo)
Photo_var <- brm(Photo_scaled ~ 1 + (1|Species),
                 data = Protea_data,
                 family = gaussian(link = "identity"),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(cauchy(0, 5), class = sd)),
                 control = list(adapt_delta = 0.99, max_treedepth = 12),
                 save_pars = save_pars(all = TRUE),
                 iter = 4000, warmup = 2000, chains = 4,
                 file = "Output/Model_Fits/Photo_var")

Photo_var_df <- as.data.frame(Photo_var)
Photo_among_var <- apply(Photo_var_df[,5:9], 1, var)
Photo_within_var <- (Photo_var_df$sigma)^2
Photo_among_var_percent <- (Photo_among_var/(Photo_among_var + Photo_within_var))*100


# 4. Total Assimilation (Total_Assim)
LSP_var <- brm(Total_Assim_scaled ~ 1 + (1|Species),
               data = Protea_data,
               family = gaussian(link = "identity"),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(cauchy(0, 5), class = sd)),
               control = list(adapt_delta = 0.99, max_treedepth = 12),
               save_pars = save_pars(all = TRUE),
               iter = 4000, warmup = 2000, chains = 4,
               file = "Output/Model_Fits/LSP_var")

LSP_var_df <- as.data.frame(LSP_var)
LSP_among_var <- apply(LSP_var_df[,5:9], 1, var)
LSP_within_var <- (LSP_var_df$sigma)^2
LSP_among_var_percent <- (LSP_among_var/(LSP_among_var + LSP_within_var))*100


# 5. Stomatal Conductance (Cond)
Cond_var <- brm(Cond_scaled ~ 1 + (1|Species),
                data = Protea_data,
                family = gaussian(link = "identity"),
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(cauchy(0, 5), class = sd)),
                control = list(adapt_delta = 0.99, max_treedepth = 12),
                save_pars = save_pars(all = TRUE),
                iter = 4000, warmup = 2000, chains = 4,
                file = "Output/Model_Fits/Cond_var")

Cond_var_df <- as.data.frame(Cond_var)
Cond_among_var <- apply(Cond_var_df[,5:9], 1, var)
Cond_within_var <- (Cond_var_df$sigma)^2
Cond_among_var_percent <- (Cond_among_var/(Cond_among_var + Cond_within_var))*100


# 6. Instantaneous Water Use Efficiency (WUE_Instan)
WUE_Instan_var <- brm(WUE_Instan_scaled ~ 1 + (1|Species),
                      data = Protea_data,
                      family = gaussian(link = "identity"),
                      prior = c(prior(normal(0, 1), class = Intercept),
                                prior(cauchy(0, 5), class = sd)),
                      control = list(adapt_delta = 0.99, max_treedepth = 12),
                      save_pars = save_pars(all = TRUE),
                      iter = 4000, warmup = 2000, chains = 4,
                      file = "Output/Model_Fits/WUE_Instan_var")

WUE_Instan_var_df <- as.data.frame(WUE_Instan_var)
WUE_Instan_among_var <- apply(WUE_Instan_var_df[,5:9], 1, var)
WUE_Instan_within_var <- (WUE_Instan_var_df$sigma)^2
WUE_Instan_among_var_percent <- (WUE_Instan_among_var/(WUE_Instan_among_var + WUE_Instan_within_var))*100


# Make a dataframe for plotting
Perf_Trait_Var <- cbind(Ks_among_var_percent, LSC_among_var_percent,
                        Photo_among_var_percent, LSP_among_var_percent,
                        Cond_among_var_percent, WUE_Instan_among_var_percent)
Perf_Trait_Var <- as.data.frame(Perf_Trait_Var)
Perf_Trait_Var_df <- pivot_longer(data = Perf_Trait_Var,
                                  cols = everything(),
                                  names_to = "Performance_Trait",
                                  values_to = "Among_Var")

Perf_Trait_plot2 <- ggplot(data = Perf_Trait_Var_df, 
                           aes(y = Performance_Trait, x = Among_Var)) +
  xlim(0, 100) +
  stat_pointinterval(.width = c(.95, .5)) +
  ggtitle("Physiological Performance Traits") + xlab("Among Species Variance (%)") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank()) +
  scale_y_discrete(labels=c("WUE_Instan_among_var_percent" = "Instantaneous \n WUE", 
                            "LSP_among_var_percent" = "Leaf-Specific \n Photosynthesis",
                            "Photo_among_var_percent" = "Photosynthetic \n Rate",
                            "LSC_among_var_percent" = "Leaf-Specific \n Conductivity",
                            "Ks_among_var_percent" = "Hydraulic \n Conductance",
                            "Cond_among_var_percent" = "Stomatal \n Conductance"))
Perf_Trait_plot2

# Multi-panel plot
plot_grid(Str_Trait_plot2, Perf_Trait_plot2,
          labels = c("A", "B"), ncol = 1, nrow = 2)

# save
ggsave(filename = "var-partitioning_panel-plot_2.pdf", 
       path = "Output/Figures",
       height = 8, width = 6)
