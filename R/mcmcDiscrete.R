mcmcDiscrete <- function(m, discrete_vars, n_samples, n_burnin, w_size, lower, upper, max_iter,
                         verbose = TRUE, pb_update = 10) {
  # setup progress bar
  if (verbose) {
    pb <- greta:::create_progress_bar("warmup", c(n_burnin, n_samples), pb_update)
    greta:::iterate_progress_bar(pb = pb, it = 0, rejects = 0)
  }

  # find out which element corresponds to discrete variables
  discrete_name <- sapply(discrete_vars, function(x) m$dag$tf_name(get(x)$node))
  params <- m$dag$example_parameters()
  discrete <- grepl(paste(discrete_name, collapse = "|"), names(params))
  
  # initialise parameters
  params[!discrete] <- rnorm(sum(!discrete))
  params[discrete] <- rbinom(sum(discrete), 1, 0.5)
  
  # create output matrix
  samples <- matrix(NA, nrow = (n_samples + n_burnin), ncol  = length(params))
  for (i in 1:(n_burnin + n_samples)) {
    # update continuous parameters
    for (j in seq_len(sum(!discrete))) {
      # send parameters and get log density back
      x1 <- params[!discrete][j]
      m$dag$send_parameters(params)
      logy <- m$dag$log_density()
      logz <- logy - rexp(1)
      u <- runif(1) * w_size
      L <- x1 - u
      R <- x1 + (w_size - u);
      params[!discrete][j] <- L
      m$dag$send_parameters(params)
      log_dens_tmp <- m$dag$log_density()
      while ((L > lower) & (log_dens_tmp > logz)) {
        L <- L - w_size
        params[!discrete][j] <- L
        m$dag$send_parameters(params)
        log_dens_tmp <- m$dag$log_density()
      }
      
      params[!discrete][j] <- R
      m$dag$send_parameters(params)
      log_dens_tmp <- m$dag$log_density()
      while ((R < upper) & (log_dens_tmp > logz)) {
        R <- R + w_size
        params[!discrete][j] <- R
        m$dag$send_parameters(params)
        log_dens_tmp <- m$dag$log_density()
      }

      # sample until draw is in the correct range
      r0 <- max(L, lower)
      r1 <- min(R, upper)
      
      for (k in seq_len(max_iter)) {
        xs <- runif(1, r0, r1)
        params[!discrete][j] <- xs
        m$dag$send_parameters(params)
        logys <- m$dag$log_density()
        if (logys > logz)
          break
        if (xs < x1) {
          r0 <- xs
        } else {
          r1 <- xs
        }
        if ((r1 - r0) < 0.0001) {
          xs <- r0
          break
        }
        if (k == max_iter) {
          xs <- x1
          break
        }
      }
      x1 <- xs
      params[!discrete][j] <- x1    
    }
    
    # update discrete parameters
    for (j in seq_len(sum(discrete))) {
      # send parameters and get log density back
      x1 <- params[discrete][j]
      m$dag$send_parameters(params)
      logy <- m$dag$log_density()
      logz <- logy - rexp(1)
      u <- runif(1) * w_size
      L <- x1 - u
      R <- x1 + (w_size - u);
      params[discrete][j] <- ifelse(L > 0, 1, 0)
      m$dag$send_parameters(params)
      log_dens_tmp <- m$dag$log_density()
      while ((L > lower) & (log_dens_tmp > logz)) {
        L <- L - w_size
        params[discrete][j] <- ifelse(L > 0, 1, 0)
        m$dag$send_parameters(params)
        log_dens_tmp <- m$dag$log_density()
      }
      
      params[discrete][j] <- ifelse(R > 0, 1, 0)
      m$dag$send_parameters(params)
      log_dens_tmp <- m$dag$log_density()
      while ((R < upper) & (log_dens_tmp > logz)) {
        R <- R + w_size
        params[discrete][j] <- ifelse(R > 0, 1, 0)
        m$dag$send_parameters(params)
        log_dens_tmp <- m$dag$log_density()
      }
      
      # sample until draw is in the correct range
      r0 <- max(L, lower)
      r1 <- min(R, upper)
      
      for (k in seq_len(max_iter)) {
        xs <- runif(1, r0, r1)
        params[discrete][j] <- ifelse(xs > 0, 1, 0)
        m$dag$send_parameters(params)
        logys <- m$dag$log_density()
        if (logys > logz)
          break
        if (xs < x1) {
          r0 <- xs
        } else {
          r1 <- xs
        }
        if ((r1 - r0) < 0.0001) {
          xs <- r0
          break
        }
        if (k == max_iter) {
          xs <- x1
          break
        }
      }
      x1 <- xs
      params[discrete][j] <- ifelse(x1 > 0, 1, 0)
    }    
    
    numerical_rejections <- 0
    if (verbose) {
      it_state <- ifelse(i > n_burnin, i - n_burnin, i)
      greta:::iterate_progress_bar(pb = pb, it = it_state, rejects = numerical_rejections)
      if (i == n_burnin) {
        pb <- greta:::create_progress_bar("sampling", c(n_burnin, n_samples), pb_update)
        greta:::iterate_progress_bar(pb = pb, it = 0, rejects = numerical_rejections)
      }
    }

    samples[i, ] <- params
  }
  
  # you can still use the gradients of the continuous variables for HMC, if you want
  #m$dag$gradients()[!discrete]
  
  # return samples
  samples[(n_burnin + 1):nrow(samples), ]
}
