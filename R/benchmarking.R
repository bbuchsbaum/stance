#' Benchmark performance against proposal targets
#'
#' Runs small simulations for each scenario in proposal.md and fits the
#' Continuous Bayesian Decoder using all optimisations. Timing and memory
#' usage are compared against the targets.
#'
#' @param max_iter Number of VB iterations to run for each scenario
#'
#' @return Named list with time (minutes), memory usage (GB), logical flag
#'   indicating whether the target was met and processing efficiency in
#'   voxel-TRs per second.
#' @export
benchmark_against_proposal <- function(max_iter = 50) {
  scenarios <- list(
    roi = list(V = 1000, T = 500, K = 3, r = 10,
               target = "5-15 minutes", ram = 2),
    parcellated = list(V = 5000, T = 500, K = 4, r = 15,
                       target = "30-90 minutes", ram = 8),
    full_wb = list(V = 50000, T = 1000, K = 5, r = 20,
                   target = "4-12 hours", ram = 32)
  )

  results <- list()
  for (name in names(scenarios)) {
    sc <- scenarios[[name]]

    sim <- simulate_fmri_data(V = sc$V, T = sc$T, K = sc$K,
                              rank = sc$r, algorithm = "CBD",
                              verbose = FALSE)
    mask <- if (name != "roi") array(TRUE, dim = determine_spatial_dims(sc$V)) else NULL

    cbd <- ContinuousBayesianDecoder$new(
      Y = sim$Y,
      K = sc$K,
      r = sc$r,
      engine = "cpp",
      use_gmrf = !is.null(mask),
      mask = mask
    )

    elapsed <- system.time({
      cbd$fit(max_iter = max_iter, verbose = FALSE)
    })[["elapsed"]]

    mem_used <- as.numeric(lobstr::mem_used())

    results[[name]] <- list(
      time_minutes = elapsed / 60,
      memory_gb = mem_used / 1e9,
      meets_target = elapsed / 60 <= parse_time(sc$target),
      efficiency = sc$V * sc$T / elapsed
    )
  }

  results
}
