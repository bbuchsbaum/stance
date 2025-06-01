#' Model diagnostic utilities
#'
#' Helper functions to evaluate model fit and performance. These are
#' primarily used by `diagnose_model_fit()` within the
#' `ContinuousBayesianDecoder` class.
#'
#' @name diagnostics
NULL

#' Check convergence diagnostics
#'
#' Computes simple convergence information based on the ELBO history of a
#' fitted decoder.
#'
#' @param cbd A `ContinuousBayesianDecoder` instance
#'
#' @return List with fields `converged`, `iterations` and `delta`.
#' @keywords internal
check_convergence_diagnostics <- function(cbd) {
  conv <- cbd$get_convergence()
  elbo <- conv$elbo_history
  if (length(elbo) < 2) {
    delta <- NA_real_
  } else {
    delta <- abs(elbo[length(elbo)] - elbo[length(elbo) - 1]) /
      abs(elbo[length(elbo) - 1])
  }
  list(
    converged = conv$converged,
    iterations = conv$iterations,
    delta = delta
  )
}

#' Assess parameter recovery
#'
#' Compares estimated spatial maps to a set of ground truth parameters if
#' provided. This is mainly intended for simulations.
#'
#' @param cbd A `ContinuousBayesianDecoder` instance
#' @param true_params List with true `U` and `V` matrices
#'
#' @return Numeric correlation value or `NA` if not available
#' @keywords internal
assess_parameter_recovery <- function(cbd, true_params) {
  if (is.null(true_params) ||
      !all(c("U", "V") %in% names(true_params))) {
    return(NA_real_)
  }
  est <- cbd$get_parameters()
  W_est <- est$U %*% t(est$V)
  W_true <- true_params$U %*% t(true_params$V)
  stats::cor(c(W_est), c(W_true))
}

#' Assess HRF smoothness
#'
#' Computes the average GMRF roughness of the voxel-wise HRFs.
#'
#' @param cbd A `ContinuousBayesianDecoder` instance
#'
#' @return Numeric roughness value or `NA`
#' @keywords internal
assess_hrf_smoothness <- function(cbd) {
  L <- cbd$.__enclos_env__$private$.L_gmrf
  if (is.null(L)) return(NA_real_)
  H <- cbd$get_hrf_estimates()
  compute_roughness(H, L)
}

#' Compute effective degrees of freedom for the GMRF prior
#'
#' Approximates model complexity by the rank of the Laplacian matrix.
#'
#' @param cbd A `ContinuousBayesianDecoder` instance
#'
#' @return Numeric degrees of freedom or `NA`
#' @keywords internal
compute_effective_df_gmrf <- function(cbd) {
  L <- cbd$.__enclos_env__$private$.L_gmrf
  if (is.null(L)) return(NA_real_)
  as.numeric(Matrix::rankMatrix(L))
}

#' Compute spatial autocorrelation of HRF estimates
#'
#' Uses simple pairwise correlations between neighbouring voxels.
#'
#' @param cbd A `ContinuousBayesianDecoder` instance
#'
#' @return Mean correlation across neighbour pairs or `NA`
#' @keywords internal
compute_spatial_autocorrelation <- function(cbd) {
  info <- cbd$.__enclos_env__$private$.spatial_info
  if (is.null(info)) return(NA_real_)
  H <- cbd$get_hrf_estimates()
  vals <- numeric(length(info$neighbors))
  for (i in seq_along(info$neighbors)) {
    nb <- info$neighbors[[i]]
    if (length(nb) > 0) {
      corrs <- stats::cor(t(cbind(H[i, ], H[nb, , drop = FALSE])))
      vals[i] <- mean(corrs[1, -1], na.rm = TRUE)
    } else {
      vals[i] <- NA_real_
    }
  }
  mean(vals, na.rm = TRUE)
}

#' Check performance targets
#'
#' Uses `profile_vb_iteration()` to verify runtime against proposal
#' benchmarks.
#'
#' @param cbd A `ContinuousBayesianDecoder` instance
#'
#' @return Logical, `TRUE` if targets met
#' @keywords internal
check_performance_targets <- function(cbd) {
  res <- profile_vb_iteration(cbd, target = "parcellated")
  isTRUE(res$meets_target)
}

#' Generate diagnostic HTML report
#'
#' Writes a very simple HTML report summarising the diagnostic list.
#'
#' @param diagnostics List returned by `diagnose_model_fit()`
#' @param output_dir Directory where the report should be written
#'
#' @return Path to generated HTML file (invisibly)
#' @keywords internal
generate_diagnostic_report <- function(diagnostics, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  out_file <- file.path(output_dir, "diagnostics.html")
  lines <- c(
    "<html><head><title>Model Diagnostics</title></head><body>",
    "<h1>Model Diagnostics</h1>",
    "<pre>",
    paste(capture.output(str(diagnostics)), collapse = "\n"),
    "</pre>",
    "</body></html>"
  )
  writeLines(lines, out_file)
  invisible(out_file)
}

