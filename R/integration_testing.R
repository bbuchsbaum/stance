#' Integration testing helpers for stance
#'
#' These utilities provide lightweight integration tests combining
#' features developed across Sprints 1-3. They are primarily used in
#' unit tests and benchmarking scripts.
#'
#' @name integration_testing
NULL

#' Run integration tests across configurations
#'
#' @param configs Named list of configuration lists. Each configuration
#'   can contain elements `use_gmrf`, `use_rcpp`, `lambda_h` and
#'   optional `batch_size`.
#' @param sim_data Simulation output from [simulate_fmri_data()].
#' @param max_iter Number of VB iterations for each configuration.
#'
#' @return Named list with timing information for each configuration.
#' @export
run_integration_tests <- function(configs, sim_data, max_iter = 5) {
  results <- list()
  dims <- determine_spatial_dims(sim_data$params$V)
  mask <- array(TRUE, dim = dims)

  for (name in names(configs)) {
    cfg <- configs[[name]]
    engine <- if (isTRUE(cfg$use_rcpp)) "cpp" else "R"
    start <- proc.time()[3]
    cbd <- ContinuousBayesianDecoder$new(
      Y = sim_data$Y,
      K = sim_data$params$K,
      r = sim_data$params$rank,
      hrf_basis = sim_data$params$hrf_spec,
      engine = engine,
      use_gmrf = isTRUE(cfg$use_gmrf),
      lambda_h = cfg$lambda_h %||% 1,
      mask = if (isTRUE(cfg$use_gmrf)) mask else NULL
    )
    cbd$fit(max_iter = max_iter, verbose = FALSE,
            batch_size = cfg$batch_size %||% NULL)
    elapsed <- proc.time()[3] - start
    results[[name]] <- list(
      time = elapsed,
      converged = cbd$get_convergence()$converged
    )
  }

  results
}

#' Sprint 3 integration test
#'
#' Executes a small end-to-end test incorporating spatial priors,
#' Rcpp optimisations and optional stochastic updates.
#'
#' @return Named list of results as produced by
#'   [run_integration_tests()].
#' @export
#' @examples
#' res <- test_sprint3_integration()
#' str(res)
#' 
#' @seealso [run_integration_tests()]
#' @keywords internal
test_sprint3_integration <- function() {
  set.seed(42)
  sim_data <- simulate_fmri_data(
    V = 200, T = 80, K = 3,
    algorithm = "CBD",
    snr = 1.5,
    spatial_smooth = TRUE
  )

  configs <- list(
    sprint1_baseline = list(
      use_gmrf = FALSE,
      use_rcpp = FALSE
    ),
    sprint2_rcpp = list(
      use_gmrf = FALSE,
      use_rcpp = TRUE
    ),
    sprint3_gmrf = list(
      use_gmrf = TRUE,
      use_rcpp = TRUE,
      lambda_h = 10
    ),
    sprint3_full = list(
      use_gmrf = TRUE,
      use_rcpp = TRUE,
      lambda_h = 10,
      batch_size = 20
    )
  )

  run_integration_tests(configs, sim_data, max_iter = 2)
}

