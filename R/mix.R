#' Create a mixture of gaussians for importance sampling
#'
#' A slow process, that fails regularly - uses `mixtools::mvnormalmixEM` to
#' create a mixture of `k` gaussians used in the importance sampling to generate
#' proposals (particles) from the importance sampling distribution.
#'
#' The parvector can be a difficult array to create, as its exact structure
#' depends on the prior_dist
#'
#' @param parvector A 2D array with rows = sampler iterations and columns being
#'   the model parameter estimates and other samples necessary.
#' @param k The number of multivariate gaussians to estimate in the mixture
#' @param maxit The maximum number of iterations to pass to `mvnormalmixEM`
#'
#' @return A mix object containing mu and sigma for each gaussian in mix
#'
#' @export
mix_gaussian <- function(parvector, k = 2, maxit = 5000) {
  # do k=2, for a mixture of 2 gaussians
  # (Davids suggestions, could be between 1-5 really)
  # mvnormalmixEM is a weak point - function can fail. needs a note or output to
  # show if it doesn't work. Should restart if it fails
  mix <- NULL
  while (is.null(mix)) {
    tryCatch(
      mix <- mixtools::mvnormalmixEM(
        parvector,
        k = k,
        maxit = maxit,
        verb = TRUE
      ),
      error = function(e) {
      },
      finally = {
      }
    )
  }
  mix
}


#' Gets proposal mvnormal samples based on mix of gaussians
#'
#' Currently assumes that the mix of gaussians was generated with k = 2
#'
#' @param parvector A 2D array with rows = sampler iterations and columns being
#'   the model parameter estimates and other samples necessary.
#' @param n_isamples The number of importance samples to generate
#' @param prop_args The extra arguments passed from the is2 function
#'
#' @return An array with n_isamples proposal importance sample distributions
#'   drawn from the mixture of gaussians using rmvnorm.
get_proposals <- function(parvector, n_isamples, prop_args) {
  if (is.null(prop_args$use_mix)) {
    mu <- apply(parvector, 2, mean)
    sigma <- var(parvector)

    # generate the IS proposals using the multivariate t dist
    proposals <- mvtnorm::rmvt(n_isamples, sigma = sigma, df = 1, delta = mu)
    attr(proposals, "mu") <- mu
    attr(proposals, "sigma") <- sigma
  } else {
    if (is.null(prop_args$precomputed_mix)) {
      mix <- mix_gaussian(parvector)
    } else {
      mix <- prop_args$precomputed_mix
    }
    # generate the proposal parameters from the mix of importance samples
    # use the weight to get samples for n1. n2 = samples-n1 (i.e 9000 and 1000)
    if (!inherits(mix, "mixEM")) {
      stop("mix object in `get_poposals` should be mixEM from mixtools")
    }
    if (length(mix$lambda) != 2) {
      stop("`get_proposals` only implemented for a mixture of two gaussians")
    }
    n1 <- stats::rbinom(n = 1, size = n_isamples, prob = max(mix$lambda))
    n1 <- pmax(n1, 2)
    n1 <- pmin(n1, n_isamples - 2)
    n2 <- n_isamples - n1

    proposals1 <- mvtnorm::rmvnorm(n1, mix$mu[[1]], mix$sigma[[1]])
    proposals2 <- mvtnorm::rmvnorm(n2, mix$mu[[2]], mix$sigma[[2]])
    proposals <- rbind(proposals1, proposals2)
    attr(proposals, "mix") <- mix
    proposals
  }
}
