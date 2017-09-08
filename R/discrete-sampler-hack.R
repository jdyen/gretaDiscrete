# hacking discrete variables into greta
## set working directory
setwd("~/Dropbox/Post-doc stuff/gretaDiscrete/")

## load packages
library(greta)

## source functions
source("./R/mcmcDiscrete.R")

## simulate data
# general settings
set.seed(1252105)
n_obs <- 100
num_vars <- 5

# predictors
pred_vals <- matrix(rnorm(n_obs * num_vars), ncol = num_vars)
pred_vals <- cbind(rep(1, n_obs), pred_vals)

# coefficients
beta_vals <- rnorm(num_vars + 1)

# set some betas to zero
num_zero <- ceiling(num_vars / 2)
beta_vals[sample(seq_len(num_vars), size = num_zero, replace = FALSE)] <- 0.0

# response
y <- pred_vals %*% beta_vals + rnorm(n_obs, mean = 0.0, sd = 0.5)

## setup greta model
I <- variable(dim = (num_vars + 1)) # dummy variable for a discrete variable (with no prior)
beta_est <- normal(mean = 0.0, sd = 10.0, dim = (num_vars + 1))
mu <- pred_vals %*% (I * beta_est)
distribution(y) <- normal(mu, sd = (sd(y) / 2))

# get the model
m <- model(mu, beta_est, I)

# sample from model
samples <- mcmcDiscrete(m$dag,
                        init = NULL,
                        discrete_vars = "I",
                        n_samples = 100,
                        thin = 1,
                        verbose = TRUE,
                        pb = greta:::create_progress_bar("sampling", c(100, 100), 10),
                        control = list(Lmin = 10,
                                       Lmax = 20,
                                       lower = -5.0,
                                       upper = 5.0,
                                       w_size = 0.5,
                                       max_iter = 100,
                                       epsilon = 0.005))
