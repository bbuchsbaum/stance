#' Optimize convolution strategy thresholds
#'
#' Determine when FFT-based convolution should be used based on the
#' time series length and HRF duration. Optionally performs a small
#' benchmark to adapt the threshold for the current machine.
#'
#' @param T Integer length of the time series
#' @param L_h Length of the HRF kernel
#' @param profile Logical; run micro-benchmarks to tune the threshold
#'
#' @return List with elements `fft` and `batch_fft` representing the
#'   minimum time points for each strategy.
#' @export
optimize_convolution_strategy <- function(T, L_h, profile = FALSE) {
  thresholds <- list(
    fft = 128L,
    batch_fft = 64L
  )

  if (profile) {
    times_direct <- numeric(5)
    times_fft <- numeric(5)

    for (i in seq_len(5)) {
      test_signal <- rnorm(T)
      test_kernel <- exp(-seq_len(min(L_h, 32)) / 6)
      times_direct[i] <- system.time({
        convolve_with_hrf(test_signal, test_kernel, method = "direct")
      })[3]
      times_fft[i] <- system.time({
        convolve_with_hrf(test_signal, test_kernel, method = "fft")
      })[3]
    }

    if (median(times_fft) < median(times_direct)) {
      thresholds$fft <- T
    }
  }

  thresholds
}

