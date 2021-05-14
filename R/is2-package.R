#' is2: Importance Sampling Squared
#'
#' The is2 package provides a general purpose implementation of the
#' marginal likelihood calculation methods outlined in Tran et al. (2020)
#' \doi{10.3758/s13428-020-01348-w}. This package will work with pmwg sampler
#' objects, and can be also be used with the output of other MCMC samplers.
#'
#' @section Documentation:
#' The documentation found at \url{https://newcastlecl.github.io/is2/} and
#' contains background information and motivation for the approach used in
#' this package and several detailed examples of the package in action.
#'
#' @section User input:
#' The user is expected to provide a structured object with MCMC samples, a
#' log likelihood function. The initial version works well with a pmwg sampler
#' object.
#'
#' @references
#' Tran, M. N., Scharth, M., Gunawan, D., Kohn, R., Brown, S. D., &
#' Hawkins, G. E. (2020). Robustly estimating the marginal likelihood for
#' cognitive models via importance sampling.
#' \emph{Behavior Research Methods, 1-18.}
#' @keywords internal
#' @aliases is2-package
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
