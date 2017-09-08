# test SSVS in greta using a binary threshold of a standard Gaussian

## load packages
library(greta)

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

## setup model
# priors
ind_vars <- normal(mean = 0.0, sd = 1.0, dim = (num_vars + 1))
beta_ind <- round(1 / (1 + exp(-20000 * ind_vars)))
beta_tmp <- normal(mean = zeros(num_vars + 1), sd = ones(num_vars + 1), dim = (num_vars + 1))
beta_est <- beta_ind * beta_tmp

# linear predictor
mu <- pred_vals %*% beta_est

# define likelihood
distribution(y) = normal(mu, sd = 1.0)

# compile model
mod <- model(mu, beta_est, ind_vars, beta_ind)

# draw samples from compiled model
samples <- mcmc(mod, n_samples = 200, warmup = 200)

# summarise samples
post_mean <- apply(samples[[1]], 2, mean)

# print some random stuff
compare_beta <- cbind(beta_vals,
                      post_mean[grep("beta_est", names(post_mean))],
                      post_mean[grep("beta_ind", names(post_mean))])
colnames(compare_beta) <- c("RealBeta", "EstimatedBeta", "Pr(inclusion)")
cat(paste0("The model has r2 = ", round(cor(post_mean[grep("mu", names(post_mean))], y) ** 2, 2), "\n"))
print(round(compare_beta, 2))
