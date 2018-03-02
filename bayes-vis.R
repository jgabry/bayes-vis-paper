# *******************************************
# This file contains code for fitting models
# and making figures for the paper 
# "Visualization in Bayesian workflow"
# *******************************************

# Setup -------------------------------------------------------------------
library(sp)
library(spdep)
library(dplyr)
library(rstan)
library(loo)
library(ggplot2)
library(bayesplot)

stopifnot(packageVersion("bayesplot") >= "1.4.0")

options(mc.cores = 4)
options(loo.cores = 4)
theme_set(bayesplot::theme_default(base_size = 14))

# load 'GM' SpatialPointsDataFrame
load("bayes-vis.RData")

GM@data <- GM@data %>% 
  mutate(
    log_pm25 = log(pm25), 
    log_sat = log(sat_2014)
  )

# regions via clustering to compare to WHO super-regions
average <- 
  GM@data %>% 
  group_by(iso3) %>% 
  summarise(pm25 = mean(pm25))
d <- dist(average)
hh <- hclust(d)
clust <- cutree(hh,k = 6)
GM@data$cluster_region <-
  sapply(GM@data$iso3, function(x) clust[which(average$iso3 == x)])



# Exploratory plots -------------------------------------------------------

xylabs <- labs(
  x = expression(log(satellite)), 
  y = expression(log(PM[2.5]))
)

# Plot log(pm2.5) vs log(sat), coloured by super-region
plot1 <- ggplot(GM@data, aes(y = log_pm25, x = log_sat)) +
  geom_point(
    aes(color = super_region_name), 
    alpha = 0.4,
    size = rel(0.8)
  ) + 
  scale_color_manual(
    values = 
      c("E-Eur/C-Eur/C-Asia" = "#00C094",
        "HighIncome" = "#FB61D7",
        "LatAm/Carib" = "#53B400",
        "N-Afr/MidEast" = "#A58AFF",
        "S-Asia" = "#00B6EB",
        "SE-Asia/E-Asia/Oceania" = "#C49A00",
        "Sub-Saharan Afr" = "#F8766D")
  ) +
  geom_smooth(
    method = lm, 
    color = "black", 
    size = 0.5, 
    linetype = 2
  ) + 
  coord_equal() +
  xylabs + 
  guides(color = guide_legend(
    title = NULL, 
    override.aes = list(alpha = 1, size = 2)
  )) + 
  theme(legend.text = element_text(size = rel(0.6)))
  

plot(plot1)
ggsave(filename = "plots/plot1.png", width = 6, height = 3)



# Plot: log(pm2.5) vs log(sat) with super-region trends
plot2 <-
  ggplot(GM@data, aes(
    y = log_pm25,
    x = log_sat
  )) +
  geom_point(aes(colour = super_region_name), alpha = 0.2, size = rel(0.75)) + 
  geom_smooth(
    method = lm, 
    color = "black", 
    size = 0.5, 
    linetype = 2
  ) + 
  scale_color_manual(
    values = 
      c("E-Eur/C-Eur/C-Asia" = "#00C094",
        "HighIncome" = "#FB61D7",
        "LatAm/Carib" = "#53B400",
        "N-Afr/MidEast" = "#A58AFF",
        "S-Asia" = "#00B6EB",
        "SE-Asia/E-Asia/Oceania" = "#C49A00",
        "Sub-Saharan Afr" = "#F8766D")
  ) +
  geom_smooth(aes(colour = super_region_name), method = lm) + 
  coord_equal() +
  xylabs + 
  legend_none() 

plot(plot2)
ggsave(filename = "plots/plot2.png", width = 4.5, height = 3.75)


# Plot: log(pm2.5) vs log(sat) with trends by regions from clustering
plot3 <-
  ggplot(GM@data, aes(
    y = log_pm25,
    x = log_sat
  )) + 
  geom_point(aes(colour = as.factor(cluster_region)), alpha = 0.2, size = rel(0.75)) + 
  scale_color_manual(
    values = 
      c("4" = "#00C094",
        "1" = "#FB61D7",
        "5" = "#53B400",
        "2" = "#A58AFF",
        "3" = "#00B6EB",
        "6" = "#C49A00",
        "7" = "#F8766D")
  ) +
  geom_smooth(
    method = lm, 
    color = "black", 
    size = 0.5, 
    linetype = 2
  ) + 
  geom_smooth(aes(colour = as.factor(cluster_region)), method = lm) + 
  coord_equal() +
  xylabs +
  legend_none()

plot(plot3)
ggsave(filename = "plots/plot3.png", width = 4.5, height = 3.75)



# Prior predictive simulations --------------------------------------------

# Plot: prior predictive with vague priors
set.seed(seed = 1923840483)

tau0 <- 1 / sqrt(rgamma(1, 1, rate = 100))
tau1 <- 1 / sqrt(rgamma(1, 1, rate = 100))
sigma <- 1 / sqrt(rgamma(1, 1, rate = 100))
beta0i <- rnorm(8, 0, tau0)
beta1i <- rnorm(8, 0, tau1)
beta0 <- rnorm(1, 0, 100)
beta1 <- rnorm(1, 0, 100)

Nsim <- length(GM@data$super_region)
xysim_labs <- labs(
  x = expression(paste("Observed ", log(PM[2.5]))),
  y = "Simulated data"
)

data1 <- data.frame(
  log_pm25 = GM$log_pm25,
  sim = beta0 + beta0i[GM$super_region] +
    (beta1 + beta1i[GM$super_region]) * GM$log_sat +
    rnorm(Nsim, mean = 0, sd = sigma)
)

theme_set(bayesplot::theme_default(base_size = 18))
theme_update(axis.text = element_text(size = 20))

ggplot(data1, aes(x = log_pm25, y = sim)) + 
  geom_point(alpha = 0.1, color = "red") + 
  xysim_labs
ggsave(filename = "plots/prior_pred_vague.png", width = 4.5, height = 3.75)


# Plot: prior predictive with weakly informative priors
set.seed(seed = 1923840479)
tau0 <- abs(rnorm(1, 0, 1))
tau1 <- abs(rnorm(1, 0, 1))
sigma <- abs(rnorm(1, 0, 1))
beta0i <- rnorm(8, 0, tau0)
beta1i <- rnorm(8, 0, tau1)
beta0 <- rnorm(1, 0, 1)
beta1 <- rnorm(1, 1, 1)

data2 <- data.frame(
  log_pm25 = GM$log_pm25,
  sim = beta0 + beta0i[GM$super_region] +
    (beta1 + beta1i[GM$super_region]) * GM$log_sat +
    rnorm(Nsim, mean = 0, sd = sigma)
)

ggplot(data2, aes(x = log_pm25, y = sim)) +
  geom_point(alpha = 0.1) + 
  xysim_labs
ggsave(filename = "plots/prior_pred_wip.png", width = 4.5, height = 3.75)


# Plot: prior predictive comparison
data3 <- data.frame(
  log_pm25 = GM$log_pm25, 
  wip = data2$sim, 
  vague = data1$sim
)
ggplot(data3, aes(x=log_pm25, y=wip)) + 
  geom_point(alpha = 0.1) + 
  geom_point(
    aes(y = vague), 
    color = "red", 
    alpha = 0.1
  ) + 
  xysim_labs
ggsave(filename = "plots/prior_pred_compare.png", width = 4.5, height = 3.75)


# Fit Stan models 1, 2, 3 -------------------------------------------------

# Compile Stan programs
# * simple.stan: simple linear regression (Model 1)
# * hier.stan: non-centered parameterization of hierarchical model (Model 2, Model 3)
simple_mod <- stan_model("stan/simple.stan")
hier_mod <- stan_model("stan/hierarchical.stan")

# Data for model 1
standata1 <- with(GM@data, list(
  N = length(log_pm25),
  log_pm = log_pm25,
  log_sat = log_sat
))

# Data for model 2 (using super-regions from WHO)
standata2 <- with(GM@data, list(
  N = length(log_pm25),
  R = length(unique(super_region)),
  log_pm = log_pm25,
  log_sat = log_sat,
  region = super_region
))

# Data for model 3 (using super-regions from clustering)
standata3 <- standata2
standata3$R <- length(unique(GM$cluster_region))
standata3$region <- GM$cluster_region

# Fit the models with Stan
nuts_controls <- list(max_treedepth = 15, adapt_delta = 0.99)
mod1 <- sampling(simple_mod, data = standata1, seed = 2402)
mod2 <- sampling(hier_mod, data = standata2, control = nuts_controls, seed = 2402)
mod3_diverge <- sampling(hier_mod, data = standata3, control = nuts_controls[1], seed = 2402)
mod3 <- sampling(hier_mod, data = standata3, control = nuts_controls, seed = 2402)

save(file = "stan/stanfits.RData", mod1, mod2, mod3, mod3_diverge)

# Extract parameter estimates, pointwise log-lik, and posterior predictive draws
keep_pars <- c("sigma", "beta0", "beta1", "beta0_region", "beta1_region", "tau0", "tau1")
posterior1 <- as.array(mod1, pars = keep_pars[1:3])
posterior2 <- as.array(mod2, pars = keep_pars)
posterior3_diverge <- as.array(mod3_diverge, pars = keep_pars)
posterior3 <- as.array(mod3, pars = keep_pars)



# HMC diagnostics ----------------------------------------------------
theme_update(axis.text = element_text(size = 20))
hmc_diagnostics <- nuts_params(mod3_diverge)

# parallel coordinates with divergences
color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "green", div_size = 0.15, div_alpha = 0.4)

parcoord <- 
  mcmc_parcoord(
    posterior3_diverge,
    regex_pars = c("tau1", "beta1_region"),
    size = 0.15,
    alpha = 0.2,
    np = hmc_diagnostics,
    np_style = div_style,
    transformations = list(tau1 = "log")
  ) +
  scale_x_discrete(
    expand = c(.05, 0),
    labels = parse(text = c("log(tau[1])", paste0("beta[1][", 1:6,"]")))
  )

plot(parcoord)
ggsave(filename = "plots/mcmc_parcoord_divs.png", width = 4.5, height = 3.75)

# scatterplot with divergences
div_style <- scatter_style_np(div_color = "green", div_size = 2.5, div_alpha = 0.75)
scatter <- 
  mcmc_scatter(
    posterior3_diverge,
    size = 1.5,
    alpha = 2/3,
    pars = c("beta1_region[1]", "tau1"), 
    transform = list(tau1 = "log"), 
    np = hmc_diagnostics,
    np_style = div_style
  ) + 
  labs(x = expression(beta[1][1]), y = expression(log(tau[1]))) + 
  xaxis_title(size = rel(1.25)) + 
  yaxis_title(size = rel(1.25))

plot(scatter)
ggsave(filename = "plots/mcmc_scatter_divs.png", width = 4.5, height = 3.75)

# Graphical posterior predictive checks -----------------------------------
theme_set(bayesplot::theme_default(base_size = 14))
y <- standata2$log_pm
yrep1 <- as.matrix(mod1, pars = "log_pm_rep")
yrep2 <- as.matrix(mod2, pars = "log_pm_rep")
yrep3 <- as.matrix(mod3, pars = "log_pm_rep")

samp100 <- sample(nrow(yrep1), 100)

# overlaid densities
color_scheme_set("blue")
ppc_dens_overlay(y, yrep1[samp100, ]) + 
  coord_cartesian(ylim = c(0, 0.7), xlim = c(0, 6)) +
  legend_none()
ggsave(filename = "plots/ppc_dens1.png", width = 4.5, height = 3.75)

color_scheme_set("gray")
ppc_dens_overlay(y, yrep2[samp100, ]) + 
  coord_cartesian(ylim = c(0, 0.7), xlim = c(0, 6)) +
  legend_none()
ggsave(filename = "plots/ppc_dens2.png", width = 4.5, height = 3.75)

color_scheme_set("red")
ppc_dens_overlay(y, yrep3[samp100, ]) + 
  coord_cartesian(ylim = c(0, 0.7), xlim = c(0, 6)) +
  legend_none()
ggsave(filename = "plots/ppc_dens3.png", width = 4.5, height = 3.75)

# stat: skew 
skew <- function(x) {
  xdev <- x - mean(x)
  n <- length(x)
  r <- sum(xdev^3) / sum(xdev^2)^1.5
  return(r * sqrt(n) * (1 - 1/n)^1.5)
}

color_scheme_set("blue")
ppc_stat(y, yrep1, stat = "skew", binwidth = 0.01) + 
  xlim(0, .6) + 
  legend_none()
ggsave(filename = "plots/ppc_skew1.png", width = 4.5, height = 3.75)

color_scheme_set("gray")
ppc_stat(y, yrep2, stat = "skew", binwidth = 0.01) + 
  xlim(0, .6) + 
  legend_none()
ggsave(filename = "plots/ppc_skew2.png", width = 4.5, height = 3.75)

color_scheme_set("red")
ppc_stat(y, yrep3, stat = "skew", binwidth = 0.01) + 
  xlim(0, .6) + 
  legend_none()
ggsave(filename = "plots/ppc_skew3.png", width = 4.5, height = 3.75)


# stat: group medians
superregion <- GM@data$super_region_name
superregion <- factor(superregion, levels = levels(superregion)[c(2,1,3:7)])

color_scheme_set("blue")
ppc_stat_grouped(y, yrep1, 
                 group = superregion, 
                 stat = "median", 
                 facet_args = list(nrow = 2)) + 
  facet_text(size = rel(0.7)) + 
  scale_x_continuous(breaks = function(x) pretty(x, n = 3)) +
  legend_none()
ggsave(filename = "plots/ppc_med_grouped1.png", height = 3, width = 7)

color_scheme_set("gray")
ppc_stat_grouped(y, yrep2, 
                 group = superregion, 
                 stat = "median",
                 facet_args = list(nrow = 2)) +
  facet_text(size = rel(0.7)) + 
  scale_x_continuous(breaks = function(x) pretty(x, n = 3)) +
  legend_none()
ggsave(filename = "plots/ppc_med_grouped2.png", height = 3, width = 7)

color_scheme_set("red")
ppc_stat_grouped(y, yrep3, 
                 group = standata3$region,
                 stat = "median", 
                 facet_args = list(nrow = 2)) + 
  facet_text(size = rel(0.7)) + 
  scale_x_continuous(breaks = function(x) pretty(x, n = 3)) +
  legend_none()
ggsave(filename = "plots/ppc_med_grouped3.png", height = 3, width = 7)


# LOO-PIT plots ----------------------------------------------------------
loglik1 <- as.matrix(mod1, pars = "log_lik")
loglik2 <- as.matrix(mod2, pars = "log_lik")
loglik3 <- as.matrix(mod3, pars = "log_lik")
loo1 <- loo(loglik1)
loo2 <- loo(loglik2)
loo3 <- loo(loglik3)
psis1 <- psislw(-loglik1)
psis2 <- psislw(-loglik2)
psis3 <- psislw(-loglik3)
pit1 <- rstantools::loo_pit(object = yrep1, y = y, lw = psis1$lw_smooth)
pit2 <- rstantools::loo_pit(object = yrep2, y = y, lw = psis2$lw_smooth)
pit3 <- rstantools::loo_pit(object = yrep3, y = y, lw = psis3$lw_smooth)

unifs <- matrix(runif(length(pit1) * 100), nrow = 100)
coords <- coord_cartesian(ylim = c(0, 2), xlim = c(0.1, 0.9))

color_scheme_set("blue")
ppc_dens_overlay(pit1, unifs) + legend_none() + coords
ggsave(filename = "plots/ppc_loo_pit_overlay1.png", width = 4.5, height = 3.75)

color_scheme_set("darkgray")
ppc_dens_overlay(pit2, unifs) + legend_none() + coords
ggsave(filename = "plots/ppc_loo_pit_overlay2.png", width = 4.5, height = 3.75)

color_scheme_set("red")
ppc_dens_overlay(pit3, unifs) + legend_none() + coords
ggsave(filename = "plots/ppc_loo_pit_overlay3.png", width = 4.5, height = 3.75)



# Pareto k diagnostics ----------------------------------------------------
theme_update(axis.text = element_text(size = 16))
k2 <- psis2$pareto_k
kdata <- data.frame(Index = seq_along(k2), khat = k2)
ggplot(kdata, aes(x = Index, y = khat)) + 
  hline_at(c(0, 0.5, 1), size = 0.25, linetype = 2, color = "gray") +
  geom_point(size = 0.5, alpha = 0.9) + 
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(y = "k-hat") + 
  annotate(
    "text", 
    x = pareto_k_ids(psis2, 0.7), 
    y = kdata$khat[pareto_k_ids(psis2, 0.7)] - 0.04, 
    label = pareto_k_ids(psis2, 0.7), 
    color = "red", size = 4
  ) 
ggsave(filename = "plots/loo_khat.png", width = 4.5, height = 3.75)



# Pointwise ELPD difference -----------------------------------------------
elpdi1 <- loo1$pointwise[, "elpd_loo"]
elpdi2 <- loo2$pointwise[, "elpd_loo"]
elpdi3 <- loo3$pointwise[, "elpd_loo"]

elpd_diffs <-
  data.frame(
    country = GM@data$iso3,
    who_cluster = GM@data$super_region_name,
    our_cluster = GM@data$cluster_region,
    diff12 = elpdi2 - elpdi1,
    diff23 = elpdi3 - elpdi2
  ) %>%
  arrange(who_cluster, country) %>%
  mutate(
    big_diff12 = abs(diff12) > 6,
    big_diff23 = abs(diff23) > 6, 
    Index = 1:n()
  )

outlier_ids <- with(elpd_diffs, Index[big_diff23])
outlier_labs <- with(elpd_diffs, as.character(country[outlier_ids]))

elpd_diff_plot <- ggplot(elpd_diffs, aes(x = Index,y = diff23)) + 
  geom_point(
    aes(colour = factor(who_cluster)), 
    size = 1, 
    alpha = 0.8
  ) + 
  scale_color_manual(
    values = 
      c("E-Eur/C-Eur/C-Asia" = "#00C094",
        "HighIncome" = "#FB61D7",
        "LatAm/Carib" = "#53B400",
        "N-Afr/MidEast" = "#A58AFF",
        "S-Asia" = "#00B6EB",
        "SE-Asia/E-Asia/Oceania" = "#C49A00",
        "Sub-Saharan Afr" = "#F8766D")
  ) +
  hline_0(size = 0.25) +
  ylab(expression(ELPD[i][3] - ELPD[i][2])) +
  legend_none() +
  annotate(
    "text", 
    x = outlier_ids,
    y = elpd_diffs$diff23[outlier_ids],
    label = outlier_labs, 
    size = 4
  )
plot(elpd_diff_plot)
ggsave(filename = "plots/loo_elpd_diff_23.png", width = 4.5, height = 3.75)





# Code for supplementary material ---------------------------------------

# fit 8-schools model
schools_data <- list(
  y = c(28,  8,-3,  7,-1,  1, 18, 12),
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18),
  J = 8
)
schools_mod <- stan("stan/eight_schools.stan", data = schools_data, seed = 1000)
schools_draws <- as.array(schools_mod)
schools_diagnostics <- nuts_params(schools_mod)


# parallel coordinates
color_scheme_set("darkgray")
theme_update(axis.text = element_text(size = 20))
div_style <- parcoord_style_np(div_color = "green", div_size = 0.15, div_alpha = 0.4)
parcoord_schools <-
  mcmc_parcoord(
    schools_draws,
    size = 0.15,
    alpha = 0.2,
    np = schools_diagnostics,
    np_style = div_style
  ) +
  scale_x_discrete(
    expand = c(0, 0),
    labels = parse(text = dimnames(schools_draws)[[3]])
  )

plot(parcoord_schools)
ggsave(filename = "plots/mcmc_parcoord_divs_8school.png", width = 4.5, height = 3.75)

# scatterplot with divergences
div_style <- scatter_style_np(div_color = "green", div_size = 2.5, div_alpha = 0.75)
scatter_schools <-
  mcmc_scatter(
    schools_draws,
    size = 1.5,
    alpha = 2/3,
    pars = c("theta[1]", "tau"),
    transform = list(tau = "log"),
    np = schools_diagnostics,
    np_style = div_style
  ) +
  labs(x = expression(theta[1]), y = expression(log(tau))) +
  xaxis_title(size = rel(1.25)) +
  yaxis_title(size = rel(1.25))

plot(scatter_schools)
ggsave(filename = "plots/mcmc_scatter_divs_8school.png", width = 4.5, height = 3.75)

