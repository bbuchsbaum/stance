#' Visualization Utilities
#'
#' Basic plotting for decoder outputs.
#'
#' @name visualization
NULL

#' Plot Spatial Maps
#'
#' @param W Spatial map matrix
#' @export
plot_spatial_maps <- function(W) {
  image(t(W), main = "Spatial Maps", xlab = "Voxel", ylab = "State")
}

#' Plot State Sequence
#'
#' @param X State activation matrix
#' @export
plot_state_sequence <- function(X) {
  matplot(t(X), type = "l", lty = 1, ylab = "Activation", xlab = "Time")
}
