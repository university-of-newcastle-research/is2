#######       IS2 Psuedo Code       #########

## set up environment and packages 
rm(list=ls())
library(mvtnorm)
library(MCMCpack)
library(rtdists)
library(invgamma)
library(mixtools)
library(condMVNorm)
library(parallel)
setwd("~/Documents/Research/Modelling Project/Work/recovery/de-mcmc")

load("3bs_fit.RData")

cpus = 1
###### set up variables #####
# number of particles, samples, subjects, random effects etc
n_randeffect=length(theta[1,,1,1])
n_subjects = length(theta[1,1,,1])
n_iter = length(theta[1,1,1,])- burnin
#length_draws = sampled$samples$idx #length of the full transformed random effect vector and/or parameter vector
IS_samples = 54 #number of importance samples
n_particles = 60 #number of particles
par.names = names(theta[1,,1,1])
hpar.names = names(phi[1,,1])


# grab the sampled stage of PMwG
# store the random effects - theta
theta_samples <- theta[,,,burnin:(burnin+n_iter)]
theta_samples <- aperm(theta_samples, c(1,4,2,3))
dim(theta_samples)<- c(((n_iter+1)*n.chains),n_randeffect,n_subjects)
theta_samples <- theta_samples[seq(1, length(theta_samples[,1,1]), 10),,] #thinning

#sampled$samples$alpha[,,sampled$samples$stage=="sample"]
# store the mu - phi
phi_samples <- phi[,,burnin:(burnin+n_iter)]
phi_samples <- aperm(phi_samples, c(1,3,2))
dim(phi_samples)<- c(((n_iter+1)*n.chains),n_randeffect*2)
phi_samples <- t(phi_samples[seq(1, length(phi_samples[,1]), 10),]) #thinning

#sampled$samples$theta_mu[,sampled$samples$stage=="sample"]
# store the cholesky transformed sigma
# sig <- sampled$samples$theta_sig[,,sampled$samples$stage=="sample"]
# the a-hlaf is used in  calculating the Huang and Wand (2013) prior. 
# The a is a random sample from inv gamma which weights the inv wishart. The mix of inverse wisharts is the prior on the correlation matrix
# a_half <- sampled$samples$a_half[,sampled$samples$stage=="sample"]


n_iter <- length(phi_samples[1,])
n.params<- n_randeffect*2+n_randeffect # theta+phi
all_samples=array(dim=c(n_subjects,n.params,n_iter))
mu_tilde=array(dim = c(n_subjects,n.params))
sigma_tilde=array(dim = c(n_subjects,n.params,n.params))

#This is pretty much the same. Get theta and phi together for each subject and 
#then get the mu_tilde as the means and the sigma as the covariance of those 21 means
for (j in 1:n_subjects){
  theta_new <- theta_samples[,,j]
  all_samples[j,,] = rbind(t(theta_new),phi_samples[,])
  # calculate the mean for re, mu and sigma
  mu_tilde[j,] =apply(all_samples[j,,],1,mean)
  # calculate the covariance matrix for random effects, mu and sigma
  sigma_tilde[j,,] = cov(t(all_samples[j,,]))
}


#### this is the same, just bind theta and phi together in the right order?
X <- t(phi_samples)
# do k=2, for a mixture of 2 gaussians (Davids suggestions, could be between 1-5 really)
k = 2 #number of dists

#mvnormalmixEM is a weak point - function can fail. needs a note or output to show if it doesn't work. Should restart if it fails
mix = NULL
while(is.null(mix)) {
  tryCatch(mix<-mvnormalmixEM(X,k=k, maxit = 5000),error=function(e){
  },finally={})
}


mix_weight <- mix$lambda
mix_mu <- mix$mu
mix_sigma <- mix$sigma

###### generate the proposal parameters from the mix of importance samples  #####

# use the weight to get samples for n1. n2 = samples-n1 (i.e 9000 and 1000)
n1=rbinom(n=1, size=IS_samples,prob=max(mix_weight))
n1=pmax(n1,2)
n1=pmin(n1,IS_samples-2)
n2=IS_samples-n1

# generates the 10,000 IS proposals given the mix
proposals1=rmvnorm(n1,mix_mu[[1]],mix_sigma[[1]])
proposals2=rmvnorm(n2,mix_mu[[2]],mix_sigma[[2]])
prop_theta=rbind(proposals1,proposals2)




### main functions

# get_params = function(parameters, n_randeffect){
#   param.theta.mu <- parameters[1:n_randeffect]
#   param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)] ##scott would like it to ask for n(unwind)
#   param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE)
#   param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
# }



##### make it work

#makes an array to store the IS samples
tmp<-array(dim=c(IS_samples))

#do the sampling
if (cpus>1){
  tmp <- mclapply(X=1:IS_samples,mc.cores = cpus, FUN = compute_lw, prop_theta = prop_theta,data = data,n_subjects= n_subjects,n_particles = n_particles,
                  n_randeffect = n_randeffect,mu_tilde=mu_tilde,sigma_tilde = sigma_tilde, prior_dist=prior_dist)
} else{
  for (i in 1:IS_samples){
    cat(i, sep="_")
    tmp[i]<-compute_lw(prop_theta,data,n_subjects,n_particles, n_randeffect,mu_tilde,sigma_tilde,i,prior_dist=prior_dist)
  }
}


# get the ML value
finished <- tmp
tmp<-unlist(tmp)
max.lw <- max(tmp)
mean.centred.lw <- mean(exp(tmp-max.lw)) #takes off the max and gets mean (avoids infs)
lw <-log(mean.centred.lw)+max.lw #puts max back on to get the lw



##### bootstrapping for SE ######
bootstrap = 10000
log_marglik_boot= array(dim = bootstrap)
for (i in 1:bootstrap){
  log_weight_boot = sample(tmp, IS_samples, replace = TRUE) #resample with replacement from the lw
  max.boot <- max(log_weight_boot)
  centred.boot <- mean(exp(log_weight_boot-max.boot)) #takes off the max and gets mean (avoids infs)
  log_marglik_boot[i] <-log(centred.boot)+max.boot #puts max back on 
}
var(log_marglik_boot)


save.image("IS2_de.Rdata")

