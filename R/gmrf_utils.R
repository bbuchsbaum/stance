#' GMRF Utilities for stance
#'
#' Functions to construct graph Laplacians used for spatial smoothing
#' of HRF coefficients. These helpers work with neuroim2 spatial
#' structures or generic matrices.
#'
#' @name gmrf_utils
NULL

#' Create GMRF Laplacian from neuroim2 Space
#'
#' Builds a sparse graph Laplacian using the voxel grid defined by a
#' `neuroim2::NeuroSpace` object. A 6-neighbour connectivity pattern is
#' used.
#'
#' @param space A `neuroim2::NeuroSpace` describing the 3D grid.
#' @param n_voxels Optional number of voxels. Defaults to the product of
#'   the first three dimensions of `space`.
#'
#' @return A sparse matrix representing the graph Laplacian.
#' @export
create_gmrf_laplacian_neuroim2 <- function(space, n_voxels = prod(dim(space)[1:3])) {
  dims <- dim(space)[1:3]
  coords <- arrayInd(seq_len(n_voxels), .dim = dims)
  offsets <- matrix(c(-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1), ncol = 3, byrow = TRUE)

  neighbors <- vector("list", n_voxels)
  for (v in seq_len(n_voxels)) {
    nb_coords <- sweep(offsets, 2, coords[v, ], "+")
    inside <- nb_coords[,1] >= 1 & nb_coords[,1] <= dims[1] &
              nb_coords[,2] >= 1 & nb_coords[,2] <= dims[2] &
              nb_coords[,3] >= 1 & nb_coords[,3] <= dims[3]
    nb_coords <- nb_coords[inside, , drop = FALSE]
    if (nrow(nb_coords) > 0) {
      nb_idx <- (nb_coords[,1] - 1) + (nb_coords[,2] - 1) * dims[1] +
                (nb_coords[,3] - 1) * dims[1] * dims[2] + 1
      neighbors[[v]] <- nb_idx
    } else {
      neighbors[[v]] <- integer(0)
    }
  }

  i <- integer(0)
  j <- integer(0)
  x <- numeric(0)

  for (v in seq_len(n_voxels)) {
    nb <- neighbors[[v]]
    n_nb <- length(nb)
    # Diagonal entry: degree
    i <- c(i, v)
    j <- c(j, v)
    x <- c(x, n_nb)
    if (n_nb > 0) {
      i <- c(i, rep(v, n_nb))
      j <- c(j, nb)
      x <- c(x, rep(-1, n_nb))
    }
  }

  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n_voxels, n_voxels))
}

#' Create Chain Graph Laplacian
#'
#' Generates a simple tridiagonal Laplacian matrix for a chain graph
#' with an arbitrary number of voxels.
#'
#' @param n_voxels Number of voxels in the chain.
#'
#' @return A sparse matrix Laplacian.
#' @export
create_chain_laplacian <- function(n_voxels) {
  if (n_voxels <= 0) {
    stop("n_voxels must be positive")
  }
  diag_vals <- rep(2, n_voxels)
  diag_vals[c(1, n_voxels)] <- 1

  i <- c(seq_len(n_voxels), seq_len(n_voxels - 1), seq_len(n_voxels - 1) + 1)
  j <- c(seq_len(n_voxels), seq_len(n_voxels - 1) + 1, seq_len(n_voxels - 1))
  x <- c(diag_vals, rep(-1, 2 * (n_voxels - 1)))

  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n_voxels, n_voxels))
}

