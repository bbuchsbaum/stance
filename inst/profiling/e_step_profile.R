# E-step profiling (S2-T04)
#
# Profiles compute_log_likelihoods and forward/backward passes
# using profvis, lineprof, and Rprofmem.

library(stance)

if (!requireNamespace("profvis", quietly = TRUE)) {
  stop("Package 'profvis' required for profiling")
}
if (!requireNamespace("lineprof", quietly = TRUE)) {
  stop("Package 'lineprof' required for profiling")
}

#' Profile E-step components of ContinuousBayesianDecoder
#'
#' @param V Number of voxels
#' @param T Number of time points
#' @param K Number of states
#' @param r Rank of spatial factorisation
#' @return list with profvis, lineprof, and Rprofmem outputs
#' @export
profile_e_step <- function(V = 2000, T = 300, K = 4, r = 5) {
  sim <- simulate_fmri_data(V = V, T = T, K = K, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = K, r = r)

  prof <- profvis::profvis({
    loglik <- cbd$.__enclos_env__$private$.compute_log_likelihoods()
    cbd$.__enclos_env__$private$.forward_backward(loglik)
  })

  lp <- lineprof::lineprof({
    loglik <- cbd$.__enclos_env__$private$.compute_log_likelihoods()
    cbd$.__enclos_env__$private$.forward_backward(loglik)
  })

  tmp <- tempfile()
  utils::Rprofmem(tmp)
  loglik <- cbd$.__enclos_env__$private$.compute_log_likelihoods()
  cbd$.__enclos_env__$private$.forward_backward(loglik)
  utils::Rprofmem(NULL)
  mem <- readLines(tmp)
  unlink(tmp)

  list(profvis = prof, lineprof = lp, mem = mem)
}

if (interactive()) {
  res <- profile_e_step()
  print(res$profvis)
}
