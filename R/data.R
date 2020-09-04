#' Small pmwg set of samples for SDT model of LDT
#'
#' A pmwgs object from the pmwg package that contains samples from a short run
#' of the sampling stages using a low number of particles. Primarily used for
#' testing the functionality of the package.
#'
#'
#' @format A pmwgs object with 427 samples, 19 subjects random effects and 
#' 5 model parameters. A pmwgs object is a thin wrapper on a named list which
#' contains the following (among other elements)
#' \describe{
#'   \item{data}{The original data which is the basis of the model estimates}
#'   \item{samples}{Samples for model estimation, including group means,
#'     covariance matrix and random effect samples}
#'   \item{ll_func}{The log likelihood function used by pmwg for estimation}
#'   \item{prior}{Priors for the model parameter means}
#' }
"sdt_samples"
