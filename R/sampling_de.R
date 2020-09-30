group_dist <- function(random_effect = NULL,
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

prior_dist <- function(parameters, prior_parameters = prior, n_randeffect) { ### mod notes: the sampled$prior needs to be fixed/passed in some other time
  param.theta.mu <- parameters[(1:n_randeffect) * 2 - 1]
  names(param.theta.mu) <- par.names
  # needs par names
  param.theta.sig <- parameters[(1:n_randeffect) * 2]
  names(param.theta.sig) <- par.names
  output <- 0
  for (pars in names(param.theta.mu)) {
    output <- output + msm::dtnorm(param.theta.mu[pars], mean = prior[[pars]]$mu[1], sd = prior[[pars]]$mu[2], lower = 0, log = TRUE)
    output <- output + dgamma(param.theta.sig[pars], shape = prior[[pars]]$sig[1], rate = prior[[pars]]$sig[2], log = TRUE)
  }
  return(output)
}

get_logp <- function(prop_theta, data, n_subjects, n_particles, n_randeffect, mu_tilde, sigma_tilde, i, group_dist = group_dist) {
  # make an array for the density
  logp <- array(dim = c(n_particles, n_subjects))
  # for each subject, get 1000 IS samples (particles) and find log weight of each
  for (j in 1:n_subjects) {
    # generate the particles from the conditional MVnorm AND mix of group level proposals
    wmix <- 0.95
    n1 <- rbinom(n = 1, size = n_particles, prob = wmix)
    if (n1 < 2) n1 <- 2
    if (n1 > (n_particles - 2)) n1 <- n_particles - 2 ## These just avoid degenerate arrays.
    n2 <- n_particles - n1
    # do conditional MVnorm based on the proposal distribution
    conditional <- condMVN(
      mean = mu_tilde[j, ], sigma = sigma_tilde[j, , ], dependent.ind = 1:n_randeffect,
      given.ind = (n_randeffect + 1):n.params, X.given = prop_theta[i, 1:(n.params - n_randeffect)]
    )
    particles1 <- mvtnorm::rmvnorm(n1, conditional$condMean, conditional$condVar)
    # mix of proposal params and conditional
    particles2 <- group_dist(n_samples = n2, parameters = prop_theta[i, ], sample = TRUE, n_randeffect = n_randeffect)
    particles <- rbind(particles1, particles2)

    for (k in 1:n_particles) {
      x <- particles[k, ]
      # names for ll function to work
      # mod notes: this is the bit the prior effects
      names(x) <- par.names
      #   do lba log likelihood with given parameters for each subject, gets density of particle from ll func
      logw_first <- log.dens.like(x, data = data[[j]], par.names = par.names) # mod notes: do we pass this in or the whole sampled object????
      # below gets second part of equation 5 numerator ie density under prop_theta
      # particle k and big vector of things
      logw_second <- group_dist(random_effect = particles[k, ], parameters = prop_theta[i, ], sample = FALSE, n_randeffect = n_randeffect) # mod notes: group dist
      # below is the denominator - ie mix of density under conditional and density under pro_theta
      logw_third <- log(wmix * dmvnorm(particles[k, ], conditional$condMean, conditional$condVar) + (1 - wmix) * exp(logw_second)) # mod notes: fine?
      # does equation 5
      logw <- (logw_first + logw_second) - logw_third
      # assign to correct row/column
      logp[k, j] <- logw
    }
  }
  # we use this part to centre the logw before addign back on at the end. This avoids inf and -inf values
  sub_max <- apply(logp, 2, max)
  logw <- t(t(logp) - sub_max)
  w <- exp(logw)
  subj_logp <- log(apply(w, 2, mean)) + sub_max # means

  # this part gets the variance the same as David
  # has been commented out here, but could be an option to use this
  # var_numerator = apply(w^2, 2, sum)
  # var_denominator = apply(w, 2, sum)^2
  # variance
  # logp_variance = (var_numerator/var_denominator) - (1/n_particles) #for each subject

  # sum the logp and return
  return(sum(subj_logp))
}

compute_lw <- function(prop_theta, data, n_subjects, n_particles, n_randeffect, mu_tilde, sigma_tilde, i, prior_dist = prior_dist) {
  logp.out <- get_logp(prop_theta, data, n_subjects, n_particles, n_randeffect, mu_tilde, sigma_tilde, i, group_dist = group_dist)
  ## do equation 10
  logw_num <- logp.out[1] + prior_dist(parameters = prop_theta[i, ], prior_parameters = prior, n_randeffect)
  logw_den <- log(mix_weight[1] * mvtnorm::dmvnorm(prop_theta[i, ], mix_mu[[1]], mix_sigma[[1]]) + mix_weight[2] * mvtnorm::dmvnorm(prop_theta[i, ], mix_mu[[2]], mix_sigma[[2]])) # density of proposed params under the means
  logw <- logw_num - logw_den # this is the equation 10
  return(c(logw))
  # NOTE: we should leave a note if variance is shit - variance is given by the logp function (currently commented out)
}
