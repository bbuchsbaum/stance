#' Cross-validation for ContinuousBayesianDecoder
#'
#' Performs a simple time-blocked cross-validation over a grid of
#' hyperparameters. The function fits the model on training blocks and
#' records the final ELBO as the scoring metric.
#'
#' @param Y Data matrix (voxels x time)
#' @param mask Optional mask for spatial smoothing
#' @param rank_values Candidate ranks for the spatial decomposition
#' @param lambda_h_values Candidate GMRF precision values
#' @param n_folds Number of folds
#' @param metric Scoring metric (currently ignored)
#' @param use_parallel Use parallel computation via the future package
#'
#' @return List with the CV results, best parameter combination and a
#'   textual recommendation.
#' @export
cbd_cross_validate <- function(Y, mask = NULL,
                               rank_values = c(5, 10, 15, 20),
                               lambda_h_values = c(0.1, 1, 10, 100),
                               n_folds = 5,
                               metric = "predictive_likelihood",
                               use_parallel = TRUE) {
  T <- ncol(Y)
  fold_size <- floor(T / n_folds)
  fold_starts <- seq(1, T - fold_size + 1, by = fold_size)[1:n_folds]

  if (use_parallel &&
      requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)) {
    future::plan(future::multisession, workers = min(n_folds, 4))
    lapply_fun <- future.apply::future_lapply
  } else {
    if (use_parallel) {
      warning("future packages not available; running sequentially")
    }
    lapply_fun <- lapply
  }

  grid <- expand.grid(
    rank = rank_values,
    lambda_h = lambda_h_values,
    fold = seq_len(n_folds),
    score = NA_real_
  )

  cv_fun <- function(i) {
    cfg <- grid[i, ]
    start <- fold_starts[cfg$fold]
    idx <- start:(start + fold_size - 1)
    Y_train <- Y[, -idx, drop = FALSE]

    cbd <- ContinuousBayesianDecoder$new(
      Y = Y_train,
      K = 3,
      r = cfg$rank,
      use_gmrf = !is.null(mask),
      lambda_h = cfg$lambda_h,
      mask = mask,
      hrf_basis = "canonical"
    )
    cbd$fit(max_iter = 5, verbose = FALSE)
    tail(cbd$get_convergence()$elbo_history, 1)
  }

  grid$score <- unlist(lapply_fun(seq_len(nrow(grid)), cv_fun))

  best <- aggregate(score ~ rank + lambda_h, data = grid, mean)
  best <- best[which.max(best$score), c("rank", "lambda_h")]

  recommendation <- sprintf(
    "Use rank %s and lambda_h %s for data with %d timepoints",
    best$rank, best$lambda_h, T)

  list(
    results = grid,
    best_params = best,
    recommendation = recommendation
  )
}

#' Select best parameters from CV grid
#'
#' @param cv_results Data frame returned by `cbd_cross_validate()`
#'
#' @return Two-column data frame with the best rank and lambda_h
#' @keywords internal
select_best_parameters <- function(cv_results) {
  best <- aggregate(score ~ rank + lambda_h, data = cv_results, mean)
  best[which.max(best$score), c("rank", "lambda_h")]
}

#' Produce a simple textual recommendation
#'
#' @param Y Data matrix
#' @param best_params Data frame with `rank` and `lambda_h`
#'
#' @return Character recommendation
#' @keywords internal
recommend_parameters <- function(Y, best_params) {
  sprintf("Recommended rank %s and lambda_h %s for dataset of size %d x %d",
          best_params$rank, best_params$lambda_h, nrow(Y), ncol(Y))
}

