library(stance)
library(microbenchmark)

# legacy implementations preserved for benchmarking
get_spatial_neighbors_loop <- function(mask, connectivity = 6) {
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

create_gmrf_laplacian_loop <- function(neighbors, n_voxels) {
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

benchmark_spatial_neighbors <- function(size = c(30,30,30), connectivity = 6) {
  mask <- array(TRUE, dim = size)

  time_loop <- microbenchmark(
    {
      info <- get_spatial_neighbors_loop(mask, connectivity)
      create_gmrf_laplacian_loop(info$neighbors, length(info$neighbors))
    }, times = 5
  )

  time_cpp <- microbenchmark(
    {
      info <- get_spatial_neighbors(mask, connectivity)
      create_gmrf_laplacian(info$neighbors, length(info$neighbors))
    }, times = 5
  )

  data.frame(
    method = c("loop", "cpp"),
    median_ms = c(median(time_loop$time) / 1e6,
                  median(time_cpp$time) / 1e6)
  )
}

# Example:
# print(benchmark_spatial_neighbors())
