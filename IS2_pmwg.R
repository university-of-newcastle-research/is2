#######       IS2 Psuedo Code       #########

## set up environment and packages
rm(list = ls())
library(mixtools)
library(mvtnorm)
library(MCMCpack)
library(rtdists)
library(invgamma)
library(mixtools)
library(condMVNorm)
library(parallel)
library(msm)
devtools::load_all()

load("forstmann_long.Rdata")
message("Setup")
cpus <- 1
###### set up variables #####
# number of particles, samples, subjects, random effects etc
n_randeffect <- sampled$n_pars
n_subjects <- sampled$n_subjects
n_iter <- length(sampled$samples$stage[sampled$samples$stage == "sample"])
# length of the full transformed random effect vector and/or parameter vector
length_draws <- sampled$samples$idx
importance_samples <- 100 # number of importance samples
n_particles <- 10 # number of particles
v_alpha <- 2 # ?
pars <- sampled$pars

message("Extract necessary samples")
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

# unwound sigma
pts2_unwound <- apply(sig, 3, unwind)
n_params <- nrow(pts2_unwound) + n_randeffect + n_randeffect
all_samples <- array(dim = c(n_subjects, n_params, n_iter))
mu_tilde <- array(dim = c(n_subjects, n_params))
sigma_tilde <- array(dim = c(n_subjects, n_params, n_params))

for (j in 1:n_subjects) {
  all_samples[j, , ] <- rbind(alpha[, j, ], theta[, ], pts2_unwound[, ])
  # calculate the mean for re, mu and sigma
  mu_tilde[j, ] <- apply(all_samples[j, , ], 1, mean)
  # calculate the covariance matrix for random effects, mu and sigma
  sigma_tilde[j, , ] <- cov(t(all_samples[j, , ]))
}

message("Create parameter vector")
parvector <- cbind(t(theta), t(pts2_unwound), t(a_half))

# do k=2, for a mixture of 2 gaussians
# (Davids suggestions, could be between 1-5 really)
k <- 2 # number of dists
message("Estimate mix of gaussians")
# mvnormalmixEM is a weak point - function can fail. needs a note or output to
# show if it doesn't work. Should restart if it fails
mix <- NULL
save.image(file = "pmwg_line71.RData")
while (is.null(mix)) {
  tryCatch(
    mix <- mixtools::mvnormalmixEM(parvector, k = k, maxit = 5000),
    error = function(e) {
    },
    finally = {
    }
  )
}

save.image(file = "pmwg_line76.RData")
#### generate the proposal parameters from the mix of importance samples  ####
message("Get samples by weight")

# use the weight to get samples for n1. n2 = samples-n1 (i.e 9000 and 1000)
n1 <- rbinom(n = 1, size = importance_samples, prob = max(mix$lambda))
n1 <- pmax(n1, 2)
n1 <- pmin(n1, importance_samples - 2)
n2 <- importance_samples - n1

# generates the 10,000 IS proposals given the mix
proposals1 <- rmvnorm(n1, mix$mu[[1]], mix$sigma[[1]])
proposals2 <- rmvnorm(n2, mix$mu[[2]], mix$sigma[[2]])
prop_theta <- rbind(proposals1, proposals2)

# makes an array to store the IS samples
tmp <- array(dim = c(importance_samples))
message("Do importance Sampling")

# do the sampling
if (cpus > 1) {
  tmp <- mclapply(
    X = 1:importance_samples,
    mc.cores = cpus,
    FUN = compute_lw,
    prop_theta = prop_theta,
    data = data,
    n_subjects = n_subjects,
    n_particles = n_particles,
    n_randeffect = n_randeffect,
    mu_tilde = mu_tilde,
    sigma_tilde = sigma_tilde,
    prior_dist = prior_dist_pmwg,
    prior = sampled$prior,
    group_dist = group_dist_pmwg,
    mix = mix,
    n_params = n_params,
    par_names = sampled$par_names,
    ll_func = sampled$ll_func
  )
} else {
  for (i in 1:importance_samples) {
    cat(i)
    tmp[i] <- compute_lw(
      prop_theta,
      data,
      n_subjects,
      n_particles,
      n_randeffect,
      mu_tilde,
      sigma_tilde,
      i,
      prior_dist = prior_dist_pmwg,
      prior = sampled$prior,
      group_dist = group_dist_pmwg,
      mix = mix,
      n_params = n_params,
      par_names = sampled$par_names,
      ll_func = sampled$ll_func
    )
  }
}


save.image(file = "pmwg_line138.RData")
# get the ML value
finished <- tmp
tmp <- unlist(tmp)
max_lw <- max(tmp)
# takes off the max and gets mean (avoids infs)
mean_centred_lw <- mean(exp(tmp - max_lw))
lw <- log(mean_centred_lw) + max_lw # puts max back on to get the lw


message("Bootstrapping for Standard error")

##### bootstrapping for SE ######
bootstrap <- 10000
log_marglik_boot <- array(dim = bootstrap)
for (i in 1:bootstrap) {
  # resample with replacement from the lw
  log_weight_boot <- sample(tmp, importance_samples, replace = TRUE)
  max_boot <- max(log_weight_boot)
  # takes off the max and gets mean (avoids infs)
  centred_boot <- mean(exp(log_weight_boot - max_boot))
  log_marglik_boot[i] <- log(centred_boot) + max_boot # puts max back on
}
var(log_marglik_boot)


save.image("IS2_v2.Rdata")
