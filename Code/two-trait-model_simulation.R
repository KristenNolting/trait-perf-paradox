# Running a simulation using our 2-trait model; Plotting output

# load libraries
library(mvtnorm)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

###################################
# Two Traits, with Sampling Error #
###################################

set.seed(1234) # for reproducibility

n_sims <- 1000 # sample 1000 individuals

sigma_2 <- 0.02 # sampling variance of 0.02
y_1 <- numeric(n_sims) # empty vector to store y values under scenario 1
y_2 <- numeric(n_sims) # empty vector to store y values under scenario 2
v_1 <- matrix(c(0.5, 0.5), nrow = 1) # betas for scenario 1
v_2 <- matrix(c(0.5, -0.5), nrow = 1) # betas for scenario 2
rho <- matrix(c(1, 0.8, 0.8, 1), nrow = 2) # correlation matrix w/ rho = 0.8
x <- matrix(nrow = n_sims, ncol = 2) # empty trait matrix
eps_y_1 <- numeric(n_sims) # empty vector to store residuals for scenario 1
eps_y_2 <- numeric(n_sims) # empty vector to store residuals for scenario 2
for (i in 1:n_sims) {
  x[i, ] <- rmvnorm(1, sigma = rho)
  mu_y_1 <- v_1 %*% x[i, ]
  mu_y_2 <- v_2 %*% x[i, ]
  y_1[i] <- rnorm(1, mu_y_1, sqrt(sigma_2))
  y_2[i] <- rnorm(1, mu_y_2, sqrt(sigma_2))
  eps_y_1[i] <- y_1[i] - mu_y_1
  eps_y_2[i] <- y_2[i] - mu_y_2
}

# Take a look at how well the predicted variances match what is observed
rho_sim <- cov(x)
cat("Model 1...\n",
    "   Predicted: ", round(v_1 %*% rho_sim %*% t(v_1) + sigma_2, 3), "\n",
    "   Observed:  ", round(var(y_1), 3), "\n",
    "Model 2...\n",
    "   Predicted: ", round(v_2 %*% rho_sim %*% t(v_2) + sigma_2, 3), "\n",
    "   Observed:  ", round(var(y_2), 3), "\n")

# Scenario 1: predicted = 0.913
# Scenario 1: observed = 0.926

# Scenario 2: predicted = 0.118
# Scenario 2: observed = 0.125

############################################
# Plotting 1 Realization of the Simulation #
############################################

# Make dataframe that has the trait 1 and trait 2 data, and y1 and y2 performance data
dat <- as.data.frame(x)
colnames(dat) <- c("Trait_1", "Trait_2")
dat$Performance_1 <- y_1
dat$Performance_2 <- y_2

# We also want to plot the PC axes associated with the trait data
corre <- as.matrix(cor(x))
eigen <- eigen(corre)  # calculate eigenvectors and values
eigen
eigen$slopes[1] <- eigen$vectors[1,1]/eigen$vectors[2,1]  # calc slopes as ratios
eigen$slopes[2] <- eigen$vectors[1,1]/eigen$vectors[1,2]  # calc slopes as ratios
eigen # check

# And we want to plot the regression vectors (or "functional axes)
lm_1 <- lm(formula = Performance_1 ~ Trait_1 + Trait_2,
           data = dat)
predicted_1 <- data.frame(Perf_pred = predict(lm_1, dat), Trait_2 = dat$Trait_2)

lm_2 <- lm(formula = Performance_2 ~ Trait_1 + Trait_2,
           data = dat) 
predicted_2 <- data.frame(Perf_pred = predict(lm_2, dat), Trait_2 = dat$Trait_2)



# Plot for Scenario 1
p1 <- ggplot(data = dat, aes(x = Trait_1, y = Trait_2)) + 
  stat_ellipse(type = "norm", size = 1.5, linetype = "dashed", colour = "grey24")
p1 <- p1 + geom_point(aes(colour = Performance_1)) +
  scale_colour_distiller(type = "seq", palette = "Spectral", limits = c(-3.1,2.95), direction = 1) +
  ggtitle("Simulation: Scenario 1") + xlab("Trait 1") + ylab("Trait 2") +
  theme_classic() + labs(colour = "Performance") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  geom_smooth(color = 'magenta4', data = predicted_1, aes(x = Perf_pred, y = Trait_2), 
              se = FALSE, linetype = "dotdash", size = 1.5) +
  geom_segment(x = 0.3, y = 0, xend = eigen$values[1] + 0.3, yend = eigen$slopes[1] * eigen$values[1], 
               colour = "grey54", arrow = arrow(length = unit(0.4, "cm")), size = 1) +
  geom_segment(x = 0.3, y = 0, xend = eigen$values[2] + 0.3, yend = eigen$slopes[2] * eigen$values[2], 
               colour = "grey54", arrow = arrow(length = unit(0.4, "cm")), size = 1)

p1

# Plot for Scenario 2
p2 <- ggplot(data = dat, aes(x = Trait_1, y = Trait_2)) + 
  stat_ellipse(type = "norm", size = 1.5, linetype = "dashed", colour = "grey24")
p2 <- p2 + geom_point(aes(colour = Performance_2)) +
  scale_colour_distiller(type = "seq", palette = "Spectral", limits = c(-3.1,2.95), direction = 1) +
  ggtitle("Simulation: Scenario 2") + xlab("Trait 1") + ylab("Trait 2") +
  theme_classic() + labs(colour = "Performance") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  geom_smooth(color = 'magenta4', data = predicted_2, aes(x = Perf_pred, y = Trait_2), 
              se = FALSE, linetype = "dotdash", size = 1.5) +
  geom_segment(x = -0.2, y = -0.2, xend = eigen$values[1] - 0.2, yend = (eigen$slopes[1] * eigen$values[1]) - 0.2, 
               colour = "grey54", arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  geom_segment(x = -0.2, y = -0.2, xend = eigen$values[2] - 0.2, yend = (eigen$slopes[2] * eigen$values[2]) - 0.2, 
               colour = "grey54", arrow = arrow(length = unit(0.3, "cm")), size = 1)

p2


###############################################
# Create a sampling distribution of variances #
###############################################
set.seed(1234)
n_sims <- 1000
n_iters <- 1000
sigma_2 <- 0.02

y_1 <- numeric(n_sims)
y_2 <- numeric(n_sims)
v_1 <- matrix(c(0.5, 0.5), nrow = 1)
v_2 <- matrix(c(0.5, -0.5), nrow = 1)
rho <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
x <- matrix(nrow = n_sims, ncol = 2)
eps_y_1 <- numeric(n_sims)
eps_y_2 <- numeric(n_sims)
var_y_1 <- numeric(n_iters)
var_y_2 <- numeric(n_iters)

for(j in 1:n_iters){
  for (i in 1:n_sims) {
    x[i, ] <- rmvnorm(1, sigma = rho)
    mu_y_1 <- v_1 %*% x[i, ]
    mu_y_2 <- v_2 %*% x[i, ]
    y_1[i] <- rnorm(1, mu_y_1, sqrt(sigma_2))
    y_2[i] <- rnorm(1, mu_y_2, sqrt(sigma_2))
    eps_y_1[i] <- y_1[i] - mu_y_1
    eps_y_2[i] <- y_2[i] - mu_y_2
  } 
  var_y_1[j] <- round(var(y_1), 3)
  var_y_2[j] <- round(var(y_2), 3)
}


################################################
# Plotting Distribution of Simulated Variances #
################################################

p3 <- ggplot(data = as.data.frame(var_y_1), aes(x = var_y_1))
p3 <- p3 + geom_density(colour = "grey4", fill = "grey44") +
  xlim(0.55,1.07) +
  geom_vline(xintercept = mean(var_y_1), size = 1.5) +
  geom_vline(xintercept = quantile(var_y_1, probs = 0.025), linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = quantile(var_y_1, probs = 0.975), linetype = "dashed", size = 1.5) +
  ggtitle("Performance Variance: Scenario 1") + xlab("Variance") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank())
p3

p4 <- ggplot(data = as.data.frame(var_y_2), aes(x = var_y_2))
p4 <- p4 + geom_density(colour = "grey4", fill = "grey44") +
  xlim(0.10, 0.17) +
  geom_vline(xintercept = mean(var_y_2), size = 1.5) +
  geom_vline(xintercept = quantile(var_y_2, probs = 0.025), linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = quantile(var_y_2, probs = 0.975), linetype = "dashed", size = 1.5) +
  ggtitle("Performance Variance: Scenario 2") + xlab("Variance") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank())
p4


#####################
# Save 4-panel plot #
#####################

# Will manually add the broken x-axis in panel D

plot_grid(p1, p3, p2, p4,
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)

# And Saving
ggsave(filename = "two-trait-simulation_panel.pdf",
       path = "Output/Figures",
       height = 8, width = 10)





