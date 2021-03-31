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
