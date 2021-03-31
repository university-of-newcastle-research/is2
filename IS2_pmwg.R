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

load("forstmann_long.RData")
message("Setup")
cpus <- 1
importance_samples <- 100 # number of importance samples
n_particles <- 10 # number of particles

message("Extract necessary samples")
samples <- extract_pmwgs(x)

save.image(file = "pmwg2_line25.RData")
message("Get mix of gaussians")
mix <- mix_gaussian(samples$parvector)

#### generate the proposal parameters from the mix of importance samples  ####
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
    n_subjects = x$n_subjects,
    n_particles = n_particles,
    n_randeffect = x$n_pars,
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
    cat(i, " ")
    tmp[i] <- compute_lw(
      prop_theta,
      data,
      x$n_subjects,
      n_particles,
      x$n_pars,
      mu_tilde,
      sigma_tilde,
      i,
      prior_dist = prior_dist_pmwg,
      prior = sampled$prior,
      group_dist = group_dist_pmwg,
      mix = mix,
      n_params = n_params,
      par_names = x$par_names,
      ll_func = x$ll_func
    )
  }
}
cat("\n")

message("Get Maximum Likelihood and Bootstrap for Standard error")

save.image(file = "pmwg2_line92.RData")

summary_like <- summarise(tmp)
print(summary_like)

save.image("IS2_pmwg2.RData")
