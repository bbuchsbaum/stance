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
    coords <- neuroim2::coords(mask)
    mask_arr <- as.logical(mask > 0)
  } else {
    mask_arr <- as.logical(mask)
    coords <- which(mask_arr, arr.ind = TRUE)
  }

  offsets <- switch(as.character(connectivity),
    "6"  = rbind(c(-1,0,0), c(1,0,0), c(0,-1,0), c(0,1,0), c(0,0,-1), c(0,0,1)),
    "18" = as.matrix(expand.grid(-1:1,-1:1,-1:1))[-14,],
    "26" = as.matrix(expand.grid(-1:1,-1:1,-1:1))[-14,],
    stop("Unsupported connectivity: ", connectivity)
  )

  n_vox <- nrow(coords)
  neighbors <- vector("list", n_vox)
  dims <- dim(mask_arr)

  # precompute linear indices for all voxel coordinates
  coord_idx <- coords[,1] + (coords[,2]-1) * dims[1] +
    (coords[,3]-1) * dims[1] * dims[2]

  for (i in seq_len(n_vox)) {
    nb_coords <- sweep(offsets, 2, coords[i, ], "+")
    valid <- nb_coords[,1] > 0 & nb_coords[,1] <= dims[1] &
             nb_coords[,2] > 0 & nb_coords[,2] <= dims[2] &
             nb_coords[,3] > 0 & nb_coords[,3] <= dims[3] &
             mask_arr[ nb_coords ]
    if (any(valid)) {
      val_coords <- nb_coords[valid, , drop = FALSE]
      nb_idx <- val_coords[,1] + (val_coords[,2]-1) * dims[1] +
        (val_coords[,3]-1) * dims[1] * dims[2]
      neighbors[[i]] <- match(nb_idx, coord_idx)
      neighbors[[i]] <- neighbors[[i]][!is.na(neighbors[[i]])]
    } else {
      neighbors[[i]] <- integer(0)
    }
  }

  list(neighbors = neighbors, coords = coords)
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
  i <- integer(0)
  j <- integer(0)
  x <- numeric(0)

  for (v in seq_len(n_voxels)) {
    nb <- neighbors[[v]]
    n_nb <- length(nb)
    if (n_nb > 0) {
      i <- c(i, v)
      j <- c(j, v)
      x <- c(x, n_nb)

      i <- c(i, rep(v, n_nb))
      j <- c(j, nb)
      x <- c(x, rep(-1, n_nb))
    }
  }

  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n_voxels, n_voxels))
}
