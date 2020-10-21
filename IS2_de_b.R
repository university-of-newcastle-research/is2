#######       IS2 Psuedo Code       #########

## set up environment and packages 
rm(list=ls())
library(mixtools)
library(mvtnorm)
library(MCMCpack)
library(rtdists)
library(invgamma)
library(condMVNorm)
library(parallel)
library(msm)
#setwd("~/Documents/Research/Modelling Project/Work/IS2/de-mcmc")

load("b_fit.RData")
burnin=round(nmc/2)
cpus = 1
###### set up variables #####
# number of particles, samples, subjects, random effects etc
n_randeffect=length(theta[1,,1,1])
n_subjects = length(theta[1,1,,1])
n_iter = length(theta[1,1,1,])- burnin
#length_draws = sampled$samples$idx #length of the full transformed random effect vector and/or parameter vector
IS_samples = 100 #number of importance samples
n_particles = 100 #number of particles
par.names = names(theta[1,,1,1])
hpar.names = names(phi[1,,1])


# grab the sampled stage of PMwG
# store the random effects - theta
theta_samples <- theta[,,,burnin:5000]
theta_samples <- aperm(theta_samples, c(1,4,2,3))
dim(theta_samples)<- c((2501*20),5,40)
theta_samples <- theta_samples[seq(1, length(theta_samples[,1,1]), 10),,]

#sampled$samples$alpha[,,sampled$samples$stage=="sample"]
# store the mu - phi
phi_samples <- phi[,,burnin:5000]
phi_samples <- aperm(phi_samples, c(1,3,2))
dim(phi_samples)<- c((2501*20),10) #iterations*chains - and then 10 = hparams
phi_samples <- t(phi_samples[seq(1, length(phi_samples[,1]), 10),])

  #sampled$samples$theta_mu[,sampled$samples$stage=="sample"]
# store the cholesky transformed sigma
# sig <- sampled$samples$theta_sig[,,sampled$samples$stage=="sample"]
# the a-hlaf is used in  calculating the Huang and Wand (2013) prior. 
# The a is a random sample from inv gamma which weights the inv wishart. The mix of inverse wisharts is the prior on the correlation matrix
# a_half <- sampled$samples$a_half[,sampled$samples$stage=="sample"]


n_iter <- length(phi_samples[1,])
n.params<- n_randeffect*2+n_randeffect
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
  tryCatch(mix<-mixtools::mvnormalmixEM(X,k=k, maxit = 5000),error=function(e){
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


group_dist = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, n_randeffect){
  names(parameters)<- hpar.names
  param.theta.mu <- parameters[(1:n_randeffect)*2-1]
  param.theta.sig <- parameters[(1:n_randeffect)*2]

  
  if (sample){
    out<- msm::rtnorm(n_samples*n_randeffect, mean=param.theta.mu, sd=param.theta.sig, lower=0)
    out <- t(array(out, dim=c(n_randeffect,n_samples)))
    return(out)
  }else{
    logw_second<-msm::dtnorm(random_effect, param.theta.mu,param.theta.sig,lower=rep(0,n_randeffect),log=TRUE)
    return(sum(logw_second))
  }
}



prior_dist = function(parameters, prior = prior, n_randeffect){ ###mod notes: the sampled$prior needs to be fixed/passed in some other time
  param.theta.mu <- parameters[(1:n_randeffect)*2-1]
  names(param.theta.mu)<-par.names
  #needs par names
  param.theta.sig <- parameters[(1:n_randeffect)*2]
  names(param.theta.sig)<-par.names
  output=0
  for (pars in names(param.theta.mu)){
    output = output+ msm::dtnorm(param.theta.mu[pars], mean=prior[[pars]]$mu[1], sd= prior[[pars]]$mu[2], lower=0, log=TRUE)
    output = output+stats::dgamma(param.theta.sig[pars], shape=prior[[pars]]$sig[1], rate=prior[[pars]]$sig[2], log=TRUE)
  }
  return(output)
}


get_logp=function(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, group_dist=group_dist){
  #make an array for the density
  logp=array(dim=c(n_particles,n_subjects))
  # for each subject, get 1000 IS samples (particles) and find log weight of each
  for (j in 1:n_subjects){
    #generate the particles from the conditional MVnorm AND mix of group level proposals
    wmix <- 0.95
    n1=rbinom(n=1,size=n_particles,prob=wmix)
    if (n1<2) n1=2
    if (n1>(n_particles-2)) n1=n_particles-2 ## These just avoid degenerate arrays.
    n2=n_particles-n1
    # do conditional MVnorm based on the proposal distribution
    conditional = condMVNorm::condMVN(mean=mu_tilde[j,],sigma=sigma_tilde[j,,],dependent.ind=1:n_randeffect,
                          given.ind=(n_randeffect+1):n.params,X.given=prop_theta[i,1:(n.params-n_randeffect)])
    particles1 <- mvtnorm::rmvnorm(n1, conditional$condMean,conditional$condVar)
    # mix of proposal params and conditional
    particles2 <- group_dist(n_samples=n2, parameters = prop_theta[i,],sample=TRUE, n_randeffect=n_randeffect)
    particles <- rbind(particles1,particles2)
    
    for (k in 1:n_particles){
      x <-particles[k,]
      #names for ll function to work
      #mod notes: this is the bit the prior effects
      names(x)<-par.names
      #   do lba log likelihood with given parameters for each subject, gets density of particle from ll func
      logw_first=log.dens.like(x,data = data[[j]],par.names = par.names) #mod notes: do we pass this in or the whole sampled object????
      # below gets second part of equation 5 numerator ie density under prop_theta
      # particle k and big vector of things
      logw_second<-group_dist(random_effect = particles[k,], parameters = prop_theta[i,], sample= FALSE, n_randeffect = n_randeffect) #mod notes: group dist
      # below is the denominator - ie mix of density under conditional and density under pro_theta
      logw_third <- log(wmix*mvtnorm::dmvnorm(particles[k,], conditional$condMean,conditional$condVar)+(1-wmix)*exp(logw_second)) #mod notes: fine?
      #does equation 5
      logw=(logw_first+logw_second)-logw_third
      #assign to correct row/column
      logp[k,j]=logw 
    }
  }
  #we use this part to centre the logw before addign back on at the end. This avoids inf and -inf values
  sub_max = apply(logp,2,max)
  logw = t(t(logp) - sub_max)
  w = exp(logw)
  subj_logp = log(apply(w,2,mean))+sub_max #means
  
  # this part gets the variance the same as David
  # has been commented out here, but could be an option to use this
  #var_numerator = apply(w^2, 2, sum)
  #var_denominator = apply(w, 2, sum)^2
  #variance
  #logp_variance = (var_numerator/var_denominator) - (1/n_particles) #for each subject
  
  # sum the logp and return 
  return(sum(subj_logp))
}

compute_lw=function(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, prior_dist=prior_dist){
  
  logp.out <- get_logp(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, group_dist=group_dist)
  ##do equation 10
  logw_num <- logp.out[1]+prior_dist(parameters = prop_theta[i,], prior = prior, n_randeffect)
  logw_den <- log(mix_weight[1]* mvtnorm::dmvnorm(prop_theta[i,], mix_mu[[1]], mix_sigma[[1]])+ mix_weight[2]* mvtnorm::dmvnorm(prop_theta[i,], mix_mu[[2]], mix_sigma[[2]])) #density of proposed params under the means
  logw <- logw_num-logw_den #this is the equation 10
  return(c(logw))
  #NOTE: we should leave a note if variance is shit - variance is given by the logp function (currently commented out)
}



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

save.image("IS2_de_bs.Rdata")
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


save.image("IS2_de_bs.Rdata")

