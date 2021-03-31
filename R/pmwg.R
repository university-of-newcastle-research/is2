is2 <- function(x, importance_samples, n_particles, ...) {
  UseMethod("is2", x)
}

is2.pmwgs <- function(x, ...) {
  if (!requireNamespace("pmwg", quietly = TRUE)) {
    stop(
      "Package \"pmwg\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  # grab the sampled stage of PMwG
  # store the random effects
  alpha <- sampled$samples$alpha[, , sampled$samples$stage == "sample"]
  # store the mu
  theta <- sampled$samples$theta_mu[, sampled$samples$stage == "sample"]
  # store the cholesky transformed sigma
  sig <- sampled$samples$theta_sig[, , sampled$samples$stage == "sample"]
  # the a-half is used in calculating the Huang and Wand (2013) prior.
  # The a is a random sample from inv gamma which weights the inv wishart.
  # The mix of inverse wisharts is the prior on the correlation matrix
  a_half <- log(sampled$samples$a_half[, sampled$samples$stage == "sample"])

  pts2_unwound <- apply(sig, 3, unwind)
  n_params <- nrow(pts2_unwound) + x$n_pars + x$n_pars
  mu_tilde <- array(dim = c(x$n_subjects, n_params))
  sigma_tilde <- array(dim = c(x$n_subjects, n_params, n_params))
  v_alpha <- attr(x, "v_half")
  
  for (j in 1:x$n_subjects) {
    subject_samples <- rbind(alpha[, j, ], theta[, ], pts2_unwound[, ])
    # calculate the mean for re, mu and sigma
    mu_tilde[j, ] <- apply(subject_samples, 1, mean)
    # calculate the covariance matrix for random effects, mu and sigma
    sigma_tilde[j, , ] <- cov(t(subject_samples))
  }

  parvector <- cbind(t(theta), t(pts2_unwound), t(a_half))
  parvector2 <- cbind(
    as_mcmc(sampled, filter="sample"),
    apply(as_mcmc(sampled, selection = "theta_sig", filter="sample"), 3, unwind),
    t(log(sampled$samples$a_half[, sampled$samples$stage == "sample"]))
  )
}


