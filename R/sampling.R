#' Calculate importance samples
#'
#' @param prop_theta - one proposals for theta_mu, the group level model
#'   parameter estimates
#' @param n_particles - the number of particles to draw for each importance
#'   sample.
#' @param subj_est - A named list containing the mean `mu` vector and covariance
#'   matrices `sigma` of each subjects random effects as well as the number of
#'   parameters used for the random effect estimates.
#' @param dist_funcs - A named list with the specific calculation for the log
#'   density of the `prior` and `group` distribution
#' @param samples - The object containing original data, model design (number of
#'   subjects, number of parameters and parameter names), the specification of
#'   the prior and thelog likelihood function used to generate the model
#'   estimates.
#' @param mix - The calculated values for the importance mixing vars
#' @param show - Set to TRUE to show feedback after each call
#'
#' @return The logweight of the importance samples
#'
#' @export
compute_lw <- function(prop_theta,
                       n_particles,
                       subj_est,
                       dist_funcs,
                       mix,
                       samples,
                       show = FALSE) {
  logp_out <- get_logp(
    prop_theta,
    samples,
    n_particles,
    subj_est,
    dist_funcs$group
  )
  ## do equation 10
  logw_num <- logp_out[1] + dist_funcs$prior(
    parameters = prop_theta,
    samples$prior,
    samples$n_pars,
    samples$par_names
  )
  logw_den <- log(
    mix$lambda[1] * mvtnorm::dmvnorm(
      prop_theta,
      mix$mu[[1]],
      mix$sigma[[1]]
    ) +
    mix$lambda[2] * mvtnorm::dmvnorm(
      prop_theta,
      mix$mu[[2]],
      mix$sigma[[2]]
    )
  ) # density of proposed params under the means
  logw <- logw_num - logw_den # this is the equation 10
  if (show) {
    cat(".")
  }
  return(c(logw))
  # NOTE: we should leave a note if variance is shit -
  # variance is given by the logp function (currently commented out)
}


#' Calculate the log probability of the data given the xxx
#'
#' Given th list of X parameters do some calucltions and return stuff
#'
#' @param prop_theta - the proposals for theta_mu, the group level model
#'   parameter estimates
#' @param sampler - the object which holds the information about the sampling
#'   process, including log likelihood function, data, n_subjects etc
#' @param n_particles - the number of particles to draw for each importance
#'   sample.
#' @param group_dist - The specific calculation for the log density for the
#'   group distribution
#'
#' @return Something or ther - not sure yet.
#'
#' @export
get_logp <- function(prop_theta,
                     sampler,
                     n_particles,
                     subj_est,
                     group_dist) {
  # make an array for the density
  logp <- array(dim = c(n_particles, sampler$n_subjects))
  # ∀ subjects, get 1000 IS samples (particles) and find log weight of each
  for (j in 1:sampler$n_subjects) {
    # generate the particles from the conditional MVnorm AND
    # a mix of group level proposals
    wmix <- 0.95
    n1 <- stats::rbinom(n = 1, size = n_particles, prob = wmix)
    if (n1 < 2) n1 <- 2
    ## These just avoid degenerate arrays.
    if (n1 > (n_particles - 2)) n1 <- n_particles - 2
    n2 <- n_particles - n1
    # do conditional MVnorm based on the proposal distribution
    conditional <- condMVNorm::condMVN(
      mean = subj_est$mu_tilde[j, ],
      sigma = subj_est$sigma_tilde[j, , ],
      dependent.ind = 1:sampler$n_pars,
      given.ind = (sampler$n_pars + 1):subj_est$n_params,
      X.given = prop_theta[1:(subj_est$n_params - sampler$n_pars)]
    )
    particles1 <- mvtnorm::rmvnorm(
      n1,
      conditional$condMean,
      conditional$condVar
    )
    # mix of proposal params and conditional
    particles2 <- group_dist(
      n_samples = n2,
      parameters = prop_theta,
      sample = TRUE,
      n_randeffect = sampler$n_pars
    )
    particles <- rbind(particles1, particles2)

    for (k in 1:n_particles) {
      x <- particles[k, ]
      # names for ll function to work
      # mod notes: this is the bit the prior effects
      names(x) <- sampler$par_names
      # do lba log likelihood with given parameters for each subject,
      # gets density of particle from ll func
      logw_first <- sampler$ll_func(
        x,
        data = sampler$data[as.numeric(factor(sampler$data$subject)) == j, ]
      ) # mod notes: do we pass this in or the whole sampled object????
      # below gets second part of eq'n 5 numerator ie density under prop_theta
      # particle k and big vector of things
      logw_second <- group_dist(
        random_effect = particles[k, ],
        parameters = prop_theta,
        sample = FALSE,
        n_randeffect = sampler$n_pars
      ) # mod notes: group dist
      # below is the denominator - ie mix of density under conditional and
      # density under pro_theta
      logw_third <- log(
        wmix * mvtnorm::dmvnorm(
          particles[k, ],
          conditional$condMean,
          conditional$condVar
        ) + (1 - wmix) * exp(logw_second)
      ) # mod notes: fine?
      # does equation 5
      logw <- (logw_first + logw_second) - logw_third
      # assign to correct row/column
      logp[k, j] <- logw
    }
  }
  # we use this part to centre the logw before addign back on at the end.
  # This avoids inf and -inf values
  sub_max <- apply(logp, 2, max)
  logw <- t(t(logp) - sub_max)
  w <- exp(logw)
  subj_logp <- log(apply(w, 2, mean)) + sub_max # means

  # sum the logp and return
  return(sum(subj_logp))
}

group_dist_de <- function(random_effect = NULL,
                          parameters,
                          sample = FALSE,
                          n_samples = NULL,
                          n_randeffect) {
  theta_mu <- parameters[(1:n_randeffect) * 2 - 1]
  theta_sig <- parameters[(1:n_randeffect) * 2]

  if (sample) {
    out <- msm::rtnorm(n_samples * n_randeffect,
      mean = theta_mu,
      sd = theta_sig,
      lower = 0
    )
    out <- t(array(out, dim = c(n_randeffect, n_samples)))
    return(out)
  } else {
    logw_second <- msm::dtnorm(random_effect,
      theta_mu,
      theta_sig,
      lower = rep(0, n_randeffect),
      log = TRUE
    )
    return(sum(logw_second))
  }
}

prior_dist_de <- function(parameters, prior, n_randeffect, par_names) {
  ### mod notes: the sampled$prior needs to be fixed/passed in some other time
  theta_mu <- parameters[(1:n_randeffect) * 2 - 1]
  names(theta_mu) <- par_names
  # needs par names
  theta_sig <- parameters[(1:n_randeffect) * 2]
  names(theta_sig) <- par_names
  output <- 0
  for (pars in par_names) {
    output <- output + msm::dtnorm(
      theta_mu[pars],
      mean = prior[[pars]]$mu[1],
      sd = prior[[pars]]$mu[2],
      lower = 0,
      log = TRUE
    )
    output <- output + stats::dgamma(
      theta_sig[pars],
      shape = prior[[pars]]$sig[1],
      rate = prior[[pars]]$sig[2],
      log = TRUE
    )
  }
  return(output)
}

group_dist_pmwg <- function(random_effect = NULL,
                            parameters,
                            sample = FALSE,
                            n_samples = NULL,
                            n_randeffect) {
  theta_mu <- parameters[1:n_randeffect]
  theta_sig_unwound <- parameters[
    (n_randeffect + 1):(length(parameters) - n_randeffect)
  ]
  theta_sig <- pmwg:::wind(theta_sig_unwound)
  if (sample) {
    return(mvtnorm::rmvnorm(n_samples, theta_mu, theta_sig))
  } else {
    logw_second <- mvtnorm::dmvnorm(
      random_effect,
      theta_mu,
      theta_sig,
      log = TRUE
    )
    return(logw_second)
  }
}

prior_dist_pmwg <- function(parameters,
                            prior_parameters = NULL,
                            n_randeffect,
                            par_names = NULL) {
  theta_mu <- parameters[1:n_randeffect]
  theta_sig_unwound <- parameters[
    (n_randeffect + 1):(length(parameters) - n_randeffect)
  ] ## scott would like it to ask for n(unwind)
  theta_sig <- pmwg:::wind(theta_sig_unwound)
  param_a <- exp(
    parameters[((length(parameters) - n_randeffect) + 1):(length(parameters))]
  )
  v_alpha <- 2

  log_prior_mu <- mvtnorm::dmvnorm(
    theta_mu,
    mean = prior_parameters$theta_mu_mean,
    sigma = prior_parameters$theta_mu_var,
    log = TRUE
  )
  log_prior_sigma <- log(MCMCpack::diwish(
    theta_sig,
    v = v_alpha + n_randeffect - 1,
    S = 2 * v_alpha * diag(1 / param_a)
  )) # exp of a-half -> positive only
  log_prior_a <- sum(invgamma::dinvgamma(
    param_a,
    scale = 0.5,
    shape = 1,
    log = TRUE
  ))

  # jacobian determinant of transformation of log of the a-half
  logw_den2 <- sum(log(1 / param_a))
  # jacobian determinant of cholesky factors of cov matrix
  logw_den3 <- log(2^n_randeffect) +
    sum((n_randeffect:1 + 1) * log(diag(theta_sig)))

  return(log_prior_mu + log_prior_sigma + log_prior_a + logw_den3 - logw_den2)
}
