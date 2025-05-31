#' HRF Utilities
#'
#' Functions for handling HRF kernels using fmrireg.
#'
#' @name hrf_utils
NULL

#' Setup HRF Kernel
#'
#' @param spec HRF specification compatible with fmrireg
#' @return Numeric vector kernel
#' @export
setup_hrf_kernel <- function(spec = "canonical") {
  if (is.numeric(spec)) {
    return(spec / sum(spec))
  }
  basis <- fmrireg::create_hrf(spec)
  basis[,1] / sum(basis[,1])
}

#' Convolve with HRF
#'
#' @param X Matrix with signals in rows
#' @param hrf Kernel vector
#' @return Matrix after convolution
#' @export
convolve_with_hrf <- function(X, hrf) {
  k <- length(hrf)
  apply(X, 1, function(row) {
    stats::convolve(row, rev(hrf), type = "open")[seq_len(length(row))]
  }) %>% t()
}
