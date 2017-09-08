# hacking discrete variables into greta
y <- rnorm(10)

library (greta)
x <- normal(0, 1)
I <- variable() # dummy variable for a discrete variable (with no prior)
mu <- x + 3 * I
distribution(y) <- normal(mu, 1)

# get the model
m <- model(mu)

# hack to find out which element corresponds to I
# will need finessing for vectors
discrete_name <- m$dag$tf_name(I$node)
discrete <- names(m$dag$example_parameters()) == discrete_name

# update discrete and continuous parameters separately
params <- m$dag$example_parameters()
params[!discrete] <- rnorm(sum(!discrete))
params[discrete] <- rbinom(sum(discrete), 1, 0.5)

# send parameters and get log density back
m$dag$send_parameters(params[1:2])
m$dag$log_density()

# you can still use the gradients of the continuous variables for HMC, if you want
m$dag$gradients()[!discrete]
