#' Data Structure Utilities
#'
#' Helper functions for converting between neuroim2 objects and matrices.
#'
#' @name data_structures
NULL

#' Validate fMRI Input
#'
#' @param Y Input data matrix or NeuroVec object.
#' @return Matrix representation of Y.
#' @export
validate_fmri_input <- function(Y) {
  if (inherits(Y, "NeuroVec")) {
    return(as.matrix(Y))
  }
  if (!is.matrix(Y)) {
    stop("Y must be a matrix or NeuroVec")
  }
  Y
}

#' Extract Data Matrix
#'
#' @param obj NeuroVol/NeuroVec object.
#' @return Numeric matrix
#' @export
extract_data_matrix <- function(obj) {
  if (inherits(obj, "NeuroVec")) {
    return(as.matrix(obj))
  }
  stop("Unsupported object type")
}

#' Restore Spatial Structure
#'
#' @param mat Data matrix
#' @param reference NeuroVol/NeuroVec providing spatial metadata
#' @return NeuroVec with data from mat
#' @export
restore_spatial_structure <- function(mat, reference) {
  if (!inherits(reference, "NeuroVec")) {
    stop("reference must be a NeuroVec")
  }
  neuroim2::NeuroVec(mat, space(reference))
}
