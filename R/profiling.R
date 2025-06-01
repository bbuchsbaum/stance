#' Profile a single VB iteration
#'
#' Runs timing benchmarks on the major components of one VB iteration for a
#' fitted `ContinuousBayesianDecoder` object.
#'
#' @param cbd_object A fitted `ContinuousBayesianDecoder` instance
#'
#' @return Named list of timings in seconds
#' @keywords internal
profile_cbd_components <- function(cbd_object) {
  priv <- cbd_object$.__enclos_env__$private
  timings <- list()
  timings$log_lik <- system.time({
    ll <- priv$.compute_log_likelihoods()
  })[3]
  timings$forward_backward <- system.time({
    fb <- priv$.forward_backward(ll)
  })[3]
  timings$update_Pi <- system.time({
    priv$.update_Pi()
  })[3]
  timings$update_UV <- system.time({
    priv$.update_U_V()
  })[3]
  timings$update_sigma2 <- system.time({
    priv$.update_sigma2()
  })[3]
  timings
}

#' Parse human readable time strings
#'
#' Converts strings like "5-15 minutes" or "4-12 hours" to a numeric value in
#' minutes. The midpoint of the range is returned.
#'
#' @param x Character string describing a time range
#'
#' @return Numeric time in minutes
#' @keywords internal
parse_time <- function(x) {
  nums <- as.numeric(unlist(regmatches(x, gregexpr("[0-9]+", x))))
  if (length(nums) == 0) return(Inf)
  val <- mean(nums)
  mult <- if (grepl("hour", x)) 60 else 1
  val * mult
}

#' Suggest optimization targets
#'
#' Prints a simple message highlighting the slowest components.
#'
#' @param timings Named list of timings
#' @param target Target scenario name
#'
#' @return Invisible `NULL`
#' @keywords internal
suggest_optimizations <- function(timings, target) {
  ord <- order(unlist(timings), decreasing = TRUE)
  top <- names(timings)[ord[1:min(3, length(ord))]]
  message("Optimization suggestions for ", target, ": ", paste(top, collapse = ", "))
  invisible(NULL)
}

#' Profile full VB performance
#'
#' Provides timing and memory usage for a VB run and checks against proposal
#' performance targets.
#'
#' @param cbd_object A `ContinuousBayesianDecoder` instance
#' @param target Character string of performance scenario (`"roi"`,
#'   `"parcellated"`, or `"full_wb"`)
#'
#' @return List with timing, memory usage and logical `meets_target`
#' @export
profile_vb_iteration <- function(cbd_object, target = "parcellated") {
  targets <- list(
    roi = list(voxels = 1000, time = "5-15 minutes", ram = "<2 GB"),
    parcellated = list(voxels = 5000, time = "30-90 minutes", ram = "<8GB"),
    full_wb = list(voxels = 50000, time = "4-12 hours", ram = "<16-32GB")
  )
  current_target <- targets[[target]]

  timings <- profile_cbd_components(cbd_object)
  iter <- cbd_object$.__enclos_env__$private$.config$max_iter
  total_time <- sum(unlist(timings)) * iter
  meets <- total_time <= parse_time(current_target$time)
  if (!meets) {
    suggest_optimizations(timings, target)
  }
  list(
    timings = timings,
    memory = lobstr::mem_used(),
    meets_target = meets
  )
}
