#' Generic function for calculating importance samples
#'
#' Should be used by creating a S3 method (`is2.object`) for your specific
#' object.
#'
#' @param x The vector of parameter estimates you are working with
#' @param n_isamples The number of importance samples to generate
#' @param n_particles The number of particles to generate for each iteration of
#'   the is2 process
#' @param n_cpus The number of cpus available for parallel processing
#' @param ... Other arguments passed to methods
#'
#' @return An isample object with the importance samples.
#'
#' @export
is2 <- function(x, n_isamples, n_particles, n_cpus = 1, ...) {
  UseMethod("is2", x)
}


#' S3 method for calculating importance samples for pmwgs objects
#'
#' Calculates importance samples for a pmwgs object (from the `pmwg` package)
#'
#' @param x The vector of parameter estimates you are working with
#' @param n_isamples The number of importance samples to generate
#' @param n_particles The number of particles to generate for each iteration of
#'   the is2 process
#' @param n_cpus The number of cpus available for parallel processing
#' @param ... Other arguments passed to methods
#'
#' @return An isample object with the importance samples.
#'
#' @export
is2.pmwgs <- function(x, n_isamples, n_particles, n_cpus = 1, ...) {
  if (!requireNamespace("pmwg", quietly = TRUE)) {
    stop(
      "Package \"pmwg\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  samples <- extract_pmwgs(x)
  mix <- mix_gaussian(samples$parvector)
  prop_theta <- get_proposals(mix, n_isamples)
  # Calculate importance samples for pmwgs object
  samples <- unlist(
    parallel::mclapply(
      X = 1:n_isamples,
      mc.cores = n_cpus,
      FUN = compute_lw,
      prop_theta = prop_theta,
      data = x$data,
      n_subjects = x$n_subjects,
      n_particles = n_particles,
      n_randeffect = x$n_pars,
      mu_tilde = samples$mu_tilde,
      sigma_tilde = samples$sigma_tilde,
      prior_dist = prior_dist_pmwg,
      prior = x$prior,
      group_dist = group_dist_pmwg,
      mix = mix,
      n_params = samples$n_params,
      par_names = x$par_names,
      ll_func = x$ll_func,
      show = TRUE
    )
  )
  is2 <- list(
    samples = samples
  )
  class(is2) <- "isample"
  is2
}


#' Extract samples from a pmwgs object (from pmwg package)
#'
#' Currently in rough form, must be from sample stage, too much hard coded.
#'
#' @param x The pmwgs object passed from is2 generic function
#'
#' @return A list with various samples, possibly reshaped
extract_pmwgs <- function(x) {
  # grab the sampled stage of PMwG
  # store the random effects
  alpha <- x$samples$alpha[, , x$samples$stage == "sample"]
  # store the mu
  theta <- x$samples$theta_mu[, x$samples$stage == "sample"]
  # store the cholesky transformed sigma
  sig <- x$samples$theta_sig[, , x$samples$stage == "sample"]
  # the a-half is used in calculating the Huang and Wand (2013) prior.
  # The a is a random sample from inv gamma which weights the inv wishart.
  # The mix of inverse wisharts is the prior on the correlation matrix
  a_half <- log(x$samples$a_half[, x$samples$stage == "sample"])

  pts2_unwound <- apply(sig, 3, pmwg:::unwind)
  n_params <- nrow(pts2_unwound) + x$n_pars + x$n_pars
  mu_tilde <- array(dim = c(x$n_subjects, n_params))
  sigma_tilde <- array(dim = c(x$n_subjects, n_params, n_params))
  v_alpha <- attr(x, "v_half")

  for (j in 1:x$n_subjects) {
    subject_samples <- rbind(alpha[, j, ], theta[, ], pts2_unwound[, ])
    # calculate the mean for re, mu and sigma
    mu_tilde[j, ] <- apply(subject_samples, 1, mean)
    # calculate the covariance matrix for random effects, mu and sigma
    sigma_tilde[j, , ] <- stats::cov(t(subject_samples))
  }

  parvector <- cbind(t(theta), t(pts2_unwound), t(a_half))
  list(
    mu_tilde = mu_tilde,
    sigma_tilde = sigma_tilde,
    parvector = parvector,
    n_params = n_params
  )
}
