library(stance)
library(microbenchmark)

#' Benchmark spatial_smooth_3d implementations
#'
#' Compares the legacy loop-based smoother with the new FFT-based version
#' for a typical 10x10x10 volume.
benchmark_spatial_smooth <- function(size = c(10, 10, 10)) {
  set.seed(123)
  arr <- array(rnorm(prod(size)), dim = size)

  time_loop <- microbenchmark(
    loop = stance:::spatial_smooth_3d_loop(arr),
    times = 20
  )

  time_fft <- microbenchmark(
    fft = stance:::spatial_smooth_3d(arr),
    times = 20
  )

  data.frame(
    method = c("loop", "fft"),
    median_ms = c(median(time_loop$time) / 1e6,
                  median(time_fft$time) / 1e6)
  )
}

# Example run:
# print(benchmark_spatial_smooth())
