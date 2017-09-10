# hacking discrete variables into greta
## set working directory
setwd("~/Dropbox/Post-doc stuff/gretaDiscrete/")

## load packages
library(greta)

## source functions
source("./R/mcmcDiscrete.R")

## simulate data
# general settings
set.seed(1252100)
n_obs <- 100
num_vars <- 10

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
                        n_samples = 100,
                        warmup = 100,
                        control = list(lower = -5.0,
                                       upper = 5.0,
                                       w_size = 0.5,
                                       max_iter = 100))

# summarise fitted model
fitted_mean <- apply(samples[[1]], 2, mean)
beta_est <- cbind(round(fitted_mean[grep("beta", names(fitted_mean))], 2),
                  round(beta_vals, 2),
                  fitted_mean[grep("I", names(fitted_mean))])
colnames(beta_est) <- c("Estimated", "Real", "Pr(inclusion)")
print(beta_est)
cat(paste0("The fitted r2 is ",
           round(cor(y, fitted_mean[grep("mu", names(fitted_mean))]) ** 2, 2),
           "\n"))
