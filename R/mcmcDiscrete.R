mcmcDiscrete <- function (model,
                          discrete_vars,
                          n_samples = 1000,
                          thin = 1,
                          warmup = 100,
                          verbose = TRUE,
                          pb_update = 10,
                          control = list(),
                          initial_values = NULL) {
  
  # find variable names to label samples
  target_greta_arrays <- model$target_greta_arrays
  names <- names(target_greta_arrays)
  
  # check they're not data nodes, provide a useful error message if they are
  are_data <- vapply(target_greta_arrays,
                     function (x) inherits(x$node, 'data_node'),
                     FUN.VALUE = FALSE)
  
  if (any(are_data)) {
    is_are <- ifelse(sum(are_data) == 1, 'is a data greta array', 'are data greta arrays')
    bad_greta_arrays <- paste(names[are_data], collapse = ', ')
    msg <- sprintf('%s %s, data greta arrays cannot be sampled',
                   bad_greta_arrays,
                   is_are)
    stop (msg, call. = FALSE)
  }
  
  # get the dag containing the target nodes
  dag <- model$dag
  
  # random starting locations
  if (is.null(initial_values)) {
    
    # try several times
    valid <- FALSE
    attempts <- 1
    while (!valid & attempts < 10) {
      
      # find out which element corresponds to discrete variables
      initial_values <- dag$example_parameters()
      discrete_name <- sapply(discrete_vars, function(x) dag$tf_name(get(x)$node))
      discrete <- grepl(paste(discrete_name, collapse = "|"), names(initial_values))
      initial_values[!discrete] <- rnorm(sum(!discrete))
      initial_values[discrete] <- rbinom(sum(discrete), 1, 0.5)

      # increase the jitter each time
      initial_values[!discrete] <- rnorm(sum(!discrete), 0, 1 + attempts / 5)
      
      # test validity of values
      valid <- greta:::valid_parameters(dag, initial_values)
      attempts <- attempts + 1
      
    }
    
    if (!valid) {
      stop ('Could not find reasonable starting values after ', attempts,
            ' attempts. Please specify initial values manually via the ',
            'initial_values argument to mcmc',
            call. = FALSE)
    }
  } else {
    if (!greta:::valid_parameters(dag, initial_values)) {
      stop ('The log density and gradients could not be evaluated at these ',
            'initial values.',
            call. = FALSE)
    }
  }
  
  # get default control options
  con <- list(Lmin = 10,
              Lmax = 20,
              lower = -10.0,
              upper = 10.0,
              w_size = 1.0,
              max_iter = 10000,
              epsilon = 0.005)
  
  
  # update them with user overrides
  con[names(control)] <- control
  
  # if warmup is required, do that now and update init
  if (warmup > 0) {
    
    if (verbose)
      pb_warmup <- greta:::create_progress_bar('warmup', c(warmup, n_samples), pb_update)
    else
      pb_warmup <- NULL
    
    # run it
    warmup_draws <- samplerDiscrete(dag = dag,
                                    init = initial_values,
                                    discrete_vars = discrete_vars,
                                    n_samples = warmup,
                                    thin = thin,
                                    verbose = verbose,
                                    pb = pb_warmup,
                                    tune = TRUE,
                                    stash = FALSE,
                                    control = con)
    
    # use the last draw of the full parameter vector as the init
    initial_values <- attr(warmup_draws, 'last_x')
    con <- attr(warmup_draws, 'control')
    
  }
  
  if (verbose)
    pb_sampling <- greta:::create_progress_bar('sampling', c(warmup, n_samples), pb_update)
  else
    pb_sampling <- NULL
  
  # run the sampler
  draws <- samplerDiscrete(dag = dag,
                           init = initial_values,
                           discrete_vars = discrete_vars,
                           n_samples = n_samples,
                           thin = thin,
                           verbose = verbose,
                           pb = pb_sampling,
                           tune = FALSE,
                           stash = TRUE,
                           control = con)
  
  # if this was successful, trash the stash, prepare and return the draws
  rm('trace_stash', envir = greta_stash)
  greta:::prepare_draws(draws)
}

samplerDiscrete <- function(dag,
                            init,
                            discrete_vars,
                            n_samples,
                            thin,
                            verbose,
                            pb,
                            tune = FALSE,
                            stash = FALSE,
                            control = list(Lmin = 10,
                                           Lmax = 20,
                                           lower = -10.0,
                                           upper = 10.0,
                                           w_size = 1.0,
                                           max_iter = 10000,
                                           epsilon = 0.005)) {
  # setup progress bar
  if (verbose) {
    greta:::iterate_progress_bar(pb = pb, it = 0, rejects = 0)
  }
  
  # initialise parameters
  params <- init
  discrete_name <- sapply(discrete_vars, function(x) dag$tf_name(get(x)$node))
  discrete <- grepl(paste(discrete_name, collapse = "|"), names(params))
  
  ## setup continuous sampler
  # unpack options
  Lmin <- control$Lmin
  Lmax <- control$Lmax
  epsilon <- control$epsilon
  lower <- control$lower
  upper <- control$upper
  w_size <- control$w_size
  max_iter <- control$max_iter
  # tuning parameters
  accept_group = 50
  target_acceptance = 0.651
  kappa = 0.75
  gamma = 0.1
  numerical_rejections <- 0
  
  # get initial gradients
  dag$send_parameters(params)
  grad <- dag$gradients()[!discrete]
  logprob <- dag$log_density()
  
  if (tune)
    epsilon_trace <- rep(NA, n_samples)
  
  # set up trace store (grab values of target variables from graph to get
  # dimension and names)
  init_trace <- dag$trace_values()
  n_target <- length(init_trace)
  trace <- matrix(NA,
                  nrow = n_samples %/% thin,
                  ncol = n_target)
  colnames(trace) <- names(init_trace)
  
  # track acceptance
  accept_trace <- rep(0, n_samples)
  
  # if anything goes awry, stash the trace so far
  if (stash)
    on.exit(stash_trace(trace))
  
  accept_count <- 0
  for (i in seq_len(n_samples)) {
    # update continuous parameters
    # send parameters and get log density back
    x_old <- params[!discrete]
    logprob_old <- logprob
    grad_old <- grad
    p <- p_old <- rnorm(sum(!discrete))

    # start leapfrog steps
    reject <- FALSE
    n_steps <- base::sample(Lmin:Lmax, 1)
    for (l in seq_len(n_steps)) {
      # step
      p <- p + 0.5 * epsilon * grad
      params[!discrete] <- params[!discrete] + epsilon * p
      
      # send parameters
      dag$send_parameters(params)
      grad <- dag$gradients()[!discrete]
      
      # check gradients are finite
      if (any(!is.finite(grad))) {
        reject <- TRUE
        break
      }
      p <- p + 0.5 * epsilon * grad
    }
    
    # if the step was bad, reject it out of hand
    if (reject) {
      numerical_rejections <- numerical_rejections + 1
      params[!discrete] <- x_old
      logprob <- logprob_old
      grad <- grad_old
    } else {
      # otherwise do the Metropolis accept/reject step
      # inner products
      p_prod <- 0.5 * sum(p ^ 2)
      p_prod_old <- 0.5 * sum(p_old ^ 2)
      
      # acceptance ratio
      logprob <- dag$log_density()
      log_accept_ratio = logprob - p_prod - logprob_old + p_prod_old
      log_u = log(runif(1))
      
      if (log_u < log_accept_ratio) {
        # on acceptance, iterate the counter and leave the parameters in the dag
        # to be put in the trace
        accept_count <- accept_count + 1
        accept_trace[i] <- 1
      } else {
        # on rejection, reset all the parameters and push old parameters to the
        # graph for the trace
        params[!discrete] <- x_old
        logprob <- logprob_old
        grad <- grad_old
      }
    }

    # update discrete parameters
    for (j in seq_len(sum(discrete))) {
      # send parameters and get log density back
      x1 <- params[discrete][j]
      dag$send_parameters(params)
      logy <- dag$log_density()
      logz <- logy - rexp(1)
      u <- runif(1) * w_size
      L <- x1 - u
      R <- x1 + (w_size - u);
      params[discrete][j] <- ifelse(L > 0, 1, 0)
      dag$send_parameters(params)
      log_dens_tmp <- dag$log_density()
      while ((L > lower) & (log_dens_tmp > logz)) {
        L <- L - w_size
        params[discrete][j] <- ifelse(L > 0, 1, 0)
        dag$send_parameters(params)
        log_dens_tmp <- dag$log_density()
      }
      
      params[discrete][j] <- ifelse(R > 0, 1, 0)
      dag$send_parameters(params)
      log_dens_tmp <- dag$log_density()
      while ((R < upper) & (log_dens_tmp > logz)) {
        R <- R + w_size
        params[discrete][j] <- ifelse(R > 0, 1, 0)
        dag$send_parameters(params)
        log_dens_tmp <- dag$log_density()
      }
      
      # sample until draw is in the correct range
      r0 <- max(L, lower)
      r1 <- min(R, upper)
      
      for (k in seq_len(max_iter)) {
        xs <- runif(1, r0, r1)
        params[discrete][j] <- ifelse(xs > 0, 1, 0)
        dag$send_parameters(params)
        logys <- dag$log_density()
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
    
    # either way, store density and location of target parameters straight from the graph
    # reset dag parameters for extracting the trace
    if (i %% thin == 0) {
      dag$send_parameters(params)
      trace[i / thin, ] <- dag$trace_values()
    }
    
    # optionally tune epsilon
    if (tune) {
      # acceptance rate over the last accept_group runs
      start <- max(1, i - accept_group)
      end <- i
      accept_rate <- mean(accept_trace[start:end], na.rm = TRUE)
      
      # decrease the adaptation rate as we go
      adapt_rate <- min(1, gamma * i ^ (-kappa))
      
      # shift epsilon in the right direction, making sure it never goes negative
      epsilon <- epsilon + pmax(-(epsilon + sqrt(.Machine$double.eps)),
                                adapt_rate * (accept_rate - target_acceptance))
      
      # keep track of epsilon
      epsilon_trace[i] <- epsilon
    }
    
    if (verbose)
      greta:::iterate_progress_bar(pb = pb, it = i, rejects = numerical_rejections)
    
  }
  
  # store the tuned epsilon as the mean of the last half
  if (tune) {
    start <- floor(n_samples / 2)
    end <- n_samples
    control$epsilon <- mean(epsilon_trace[start:end], na.rm = TRUE)
  }
  
  # return samples
  attr(trace, 'last_x') <- params
  attr(trace, 'control') <- control
  trace
}
