#' Utility functions for stance
#'
#' Miscellaneous helper utilities used internally by the package.
#'
#' @name utils
NULL

#' Null coalescing operator
#'
#' Returns the first argument if it is not `NULL`, otherwise returns the
#' second.
#'
#' @param a First value
#' @param b Fallback value used when `a` is `NULL`
#'
#' @return `a` if not `NULL`, otherwise `b`
#' @export
#'
#' @examples
#' NULL %||% 1
"%||%" <- function(a, b) {
  if (is.null(a)) b else a
}
#' Normalize rows of a matrix to sum to one
#'
#' Convenience function used in stochastic VB updates.
#'
#' @param X Numeric matrix
#'
#' @return Matrix with rows scaled to sum to one
#' @keywords internal
normalize_rows <- function(X) {
  rs <- rowSums(X)
  rs[rs == 0] <- 1
  sweep(X, 1, rs, '/')
}
