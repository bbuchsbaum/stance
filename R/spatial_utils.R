#' Spatial Neighborhood Utilities
#'
#' Helper functions to work with voxel-wise spatial neighbourhoods and
#' construct sparse graph Laplacians. These utilities handle both
#' `neuroim2::NeuroVol` objects and plain logical arrays or coordinate
#' matrices so that the spatial smoothing infrastructure can operate on
#' generic data.
#'
#' @name spatial_utils
NULL

#' Extract voxel neighbours from a mask
#'
#' Returns the adjacency information for a 3D mask using the specified
#' connectivity pattern. When a `neuroim2::NeuroVol` is supplied the
#' voxel coordinates are extracted via `neuroim2::coords()` so that
#' spatial metadata are preserved.
#'
#' @param mask Logical array or `neuroim2::NeuroVol` defining the spatial
#'   mask. Non-zero voxels are considered part of the graph.
#' @param connectivity Number of neighbouring voxels to use. One of 6, 18
#'   or 26.
#'
#' @return A list with elements `neighbors` (adjacency list) and
#'   `coords` (voxel indices).
#' @export
get_spatial_neighbors <- function(mask, connectivity = 6) {
  if (inherits(mask, "NeuroVol")) {
    mask_arr <- as.logical(mask > 0)
  } else {
    mask_arr <- as.logical(mask)
  }

  dims <- dim(mask_arr)
  res <- spatial_neighbors_cpp(as.vector(mask_arr), dims, as.integer(connectivity))

  coords <- if (inherits(mask, "NeuroVol")) {
    neuroim2::coords(mask)
  } else {
    res$coords
  }

  list(neighbors = res$neighbors, coords = coords)
}

#' Create a sparse GMRF Laplacian
#'
#' Builds a sparse matrix representation of the graph Laplacian from an
#' adjacency list as returned by `get_spatial_neighbors()`.
#'
#' @param neighbors List of integer vectors giving neighbour indices.
#' @param n_voxels Number of voxels in the graph.
#'
#' @return A `dgCMatrix` sparse Laplacian.
#' @export
create_gmrf_laplacian <- function(neighbors, n_voxels) {
  lap <- laplacian_from_neighbors_cpp(neighbors)
  Matrix::sparseMatrix(i = lap$i, j = lap$j, x = lap$x,
                       dims = c(n_voxels, n_voxels))
}
