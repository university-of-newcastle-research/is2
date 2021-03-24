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

load("b_fit.RData")
message("Setup")
cpus <- 1
###### set up variables #####
# number of particles, samples, subjects, random effects etc
burnin <- round(nmc / 2)
n_randeffect <- length(theta[1, , 1, 1])
n_subjects <- length(theta[1, 1, , 1])
n_iter <- length(theta[1, 1, 1, ]) - burnin
#length of the full transformed random effect vector and/or parameter vector
importance_samples <- 100 # number of importance samples
n_particles <- 100 # number of particles
par_names <- names(theta[1, , 1, 1])

message("Extract necessary samples")
# grab the sampled stage of PMwG
# store the random effects - theta
theta_samples <- theta[, , , burnin:5000]
theta_samples <- aperm(theta_samples, c(1, 4, 2, 3))
dim(theta_samples) <- c((2501 * 20), 5, 40)
theta_samples <- theta_samples[seq(1, length(theta_samples[, 1, 1]), 10), , ]

# store the mu - phi
phi_samples <- phi[, , burnin:5000]
phi_samples <- aperm(phi_samples, c(1, 3, 2))
dim(phi_samples) <- c(((n_iter + 1) * n.chains), n_randeffect * 2)
phi_samples <- t(phi_samples[seq(1, length(phi_samples[, 1]), 10), ]) # thinning


n_iter <- length(phi_samples[1, ])
n_params <- n_randeffect * 2 + n_randeffect
all_samples <- array(dim = c(n_subjects, n_params, n_iter))
mu_tilde <- array(dim = c(n_subjects, n_params))
sigma_tilde <- array(dim = c(n_subjects, n_params, n_params))

# This is pretty much the same. Get theta and phi together for each subject and
# then get the mu_tilde as the means and the sigma as the covariance of those
# 21 means
for (j in 1:n_subjects) {
  theta_new <- theta_samples[, , j]
  all_samples[j, , ] <- rbind(t(theta_new), phi_samples[, ])
  # calculate the mean for re, mu and sigma
  mu_tilde[j, ] <- apply(all_samples[j, , ], 1, mean)
  # calculate the covariance matrix for random effects, mu and sigma
  sigma_tilde[j, , ] <- cov(t(all_samples[j, , ]))
}

message("Create parameter vector")
parvector <- t(phi_samples)

# do k=2, for a mixture of 2 gaussians
# (Davids suggestions, could be between 1-5 really)
k <- 2 # number of dists
message("Estimate mix of gaussians")
# mvnormalmixEM is a weak point - function can fail. needs a note or output to
# show if it doesn't work. Should restart if it fails
mix <- NULL
save.image(file = "de_line74.RData")
while (is.null(mix)) {
  tryCatch(
    mix <- mixtools::mvnormalmixEM(parvector, k = k, maxit = 5000),
    error = function(e) {
    },
    finally = {
    }
  )
}

save.image(file = "de_line85.RData")
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

# Autofill the par.names for log.dens.like
# log.dens.like comes from load("b_fit.RData")
adj_log_dens_like <- function(x, data) {
  log.dens.like(x, data, par.names = par_names)
}

data <- do.call(rbind, data)

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
    prior_dist = prior_dist_de,
    prior = prior,
    group_dist = group_dist_de,
    mix = mix,
    n_params = n_params,
    par_names = par_names,
    ll_func = adj_log_dens_like
  )
} else {
  for (i in 1:importance_samples) {
    cat(i, sep = "_")
    tmp[i] <- compute_lw(
      prop_theta,
      data,
      n_subjects,
      n_particles,
      n_randeffect,
      mu_tilde,
      sigma_tilde,
      i,
      prior_dist = prior_dist_de,
      prior = prior,
      group_dist = group_dist_de,
      mix = mix,
      n_params = n_params,
      par_names = par_names,
      ll_func = adj_log_dens_like
    )
  }
}


save.image(file = "de_line157.RData")
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


save.image("IS2_de_bs.RData")
