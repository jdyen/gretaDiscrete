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
n_obs <- 1000
num_vars <- 20

# predictors
pred_vals <- matrix(rnorm(n_obs * num_vars), ncol = num_vars)
pred_vals <- cbind(rep(1, n_obs), pred_vals)

# coefficients
beta_vals <- rnorm(num_vars + 1)

# set some betas to zero
num_zero <- num_vars - ceiling(num_vars / 3)
beta_vals[sample(seq_len(num_vars), size = num_zero, replace = FALSE)] <- 0.0

# response
y <- pred_vals %*% beta_vals + rnorm(n_obs, mean = 0.0, sd = 0.5)

## setup greta model
I <- variable(dim = (num_vars + 1)) # dummy variable for a discrete variable (with no prior)
beta_sd <- I * 100.0 + (1 - I) * 0.01
beta_est <- normal(mean = zeros(num_vars + 1), sd = beta_sd, dim = (num_vars + 1))
mu <- pred_vals %*% beta_est
distribution(y) <- normal(mu, sd = (sd(y) / 2))

# get the model
m <- model(mu, beta_est, I)

# sample from model
samples <- mcmcDiscrete(model = m,
                        discrete_vars = "I",
                        n_samples = 500,
                        warmup = 500,
                        control = list(lower = -5.0,
                                       upper = 5.0,
                                       w_size = 0.5,
                                       max_iter = 100))
