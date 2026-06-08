#Retrieving model evaluation stats
calc_mod_stats <- function(mod, data = NULL, response = NULL) {
  
  if (inherits(mod, "brmsfit")) {
    
    # extract data and response
    if (is.null(data)) {
      data <- mod$data
      response <- 1
    } else if (is.null(response)) {
      stop("You must enter the response column")
    }
    
    resp <- data[, response]
    
    # handle binomial vs multinomial
    fam <- mod$family$family
    
    if (fam == "bernoulli") {
      y <- as.numeric(resp) - 1
      fits_tab <- fitted(mod)
      fits <- fits_tab[, "Estimate"]
      preds <- ifelse(fits > .5, 1, 0)
      
    } else if (fam == "categorical") {
      fits_tab <- fitted(mod)
      
      # predicted class (1, 2, 3, ...)
      preds <- apply(fits_tab[, , "Estimate"], 1, which.max)
      y <- as.numeric(resp)
      
      # For now: use probability of true class for scoring
      fits <- fits_tab[cbind(seq_len(nrow(fits_tab)), y), "Estimate"]
    }
    
    # predictive metrics 
    brier_score <- mean((fits - y)^2, na.rm = TRUE) ## Brier score measures how close your predicted probabilities are to the true outcomes (0-perfect, >0.25 worse than guessing).
    log_score <- mean(abs(y*log(fits) + (1 - y) * log(1 - fits)), na.rm = TRUE) ##log predictive density (0-perfect, higher - worse)-- it detects overconfidence
    
    mean.rank <- mean(rank(fits)[y == 1])
    n <- length(y)
    n1 <- sum(y == 1)
    c.index <- (mean.rank - (n1 + 1)/2)/(n - n1) ## concordance index: how well the model ranks cases (1=perfect discrimination)
    
    # information criteria
    waic <- mod$criteria$waic$estimates[3, 1] ## predictive fit
    loo_estimates <- mod$criteria$loo$estimates ##leave-one-out test
    
    # convergence diagnostics
    summ <- posterior_summary(mod)
    max_rhat <- max(summ[, "Rhat"], na.rm = TRUE) ## should be lower than 1.01
    min_bulk_ess <- min(summ[, "Bulk_ESS"], na.rm = TRUE) ## should be larger than 400 or 1000 for better convergence
    min_tail_ess <- min(summ[, "Tail_ESS"], na.rm = TRUE)
    
    # sampler diagnostics
    sampler_params <- rstan::get_sampler_params(mod$fit, inc_warmup = FALSE)
    divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
    ebfmi <- min(sapply(sampler_params, function(x) {
      sum(x[, "energy__"]) / sum(x[, "energy__"]^2)
    }))
    
    output <- c(
      N = as.integer(n),
      baseline = max(table(y)/length(y)),
      predicted.corr = mean(preds == y),
      Brier = brier_score,
      C = c.index,
      LogScore = log_score,
      WAIC = waic,
      elpd_loo = loo_estimates[1, 1],
      p_loo = loo_estimates[2, 1],
      looic = loo_estimates[3, 1],
      max_Rhat = max_rhat,
      min_bulk_ESS = min_bulk_ess,
      min_tail_ESS = min_tail_ess,
      divergences = divergences,
      EBFMI = ebfmi
    )
    
    return(output)
  }
}

  