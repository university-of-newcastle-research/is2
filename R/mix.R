#' Create a mixture of gaussians for importance sampling
#'
#' A slow process, that fails regularly - uses `mixtools::mvnormalmixEM` to
#' create a mixture of `k` gaussians used in the importance sampling to generate
#' proposals (particles) from the importance sampling distribution.
#'
#' @param parvector A vector containing
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
      mix <- mixtools::mvnormalmixEM(parvector, k = k, maxit = maxit),
      error = function(e) {
      },
      finally = {
      }
    )
  }
  mix
}

get_proposals <- function(mix, importance_samples) {
  #### generate the proposal parameters from the mix of importance samples  ####
  # use the weight to get samples for n1. n2 = samples-n1 (i.e 9000 and 1000)
  n1 <- rbinom(n = 1, size = importance_samples, prob = max(mix$lambda))
  n1 <- pmax(n1, 2)
  n1 <- pmin(n1, importance_samples - 2)
  n2 <- importance_samples - n1

  proposals1 <- mvtnorm::rmvnorm(n1, mix$mu[[1]], mix$sigma[[1]])
  proposals2 <- mvtnorm::rmvnorm(n2, mix$mu[[2]], mix$sigma[[2]])
  rbind(proposals1, proposals2)
}
