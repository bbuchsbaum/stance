# Performance profiling for ContinuousBayesianDecoder (S1-T12)
#
# This script profiles a full VB fit on a medium sized dataset
# using profvis. Results are summarised in cbd_profile_results.md

library(stance)
if (!requireNamespace("profvis", quietly = TRUE)) {
  stop("Package 'profvis' required for profiling")
}

#' Profile ContinuousBayesianDecoder fit
#'
#' @param V Number of voxels
#' @param T Number of time points
#' @param K Number of hidden states
#' @param r Rank of spatial factorisation
#' @param hrf_basis HRF basis specification
#' @return A profvis object
#' @export
profile_cbd_fit <- function(V = 5000, T = 500, K = 5,
                           r = 10, hrf_basis = "canonical") {
  sim <- simulate_fmri_data(V = V, T = T, K = K,
                            algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(
    Y = sim$Y,
    K = K,
    r = r,
    hrf_basis = hrf_basis
  )
  profvis::profvis({
    cbd$fit(max_iter = 50, verbose = FALSE)
  })
}

if (interactive()) {
  p <- profile_cbd_fit()
  print(p)
}
