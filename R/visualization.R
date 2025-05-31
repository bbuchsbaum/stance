#' Shared Visualization Framework for stance
#'
#' Plotting functions for both CLD and CBD algorithms, providing
#' consistent visualization of spatial maps, state activations,
#' and convergence diagnostics.
#'
#' @name visualization
NULL

#' Plot Spatial Maps
#'
#' Visualizes spatial patterns from decoder output. Supports both
#' matrix format and neuroimaging objects.
#'
#' @param W Spatial map matrix (V x K) or list of NeuroVol objects
#' @param layout Grid layout c(rows, cols) for multiple maps
#' @param zlim Range for color scale (default: data range)
#' @param col Color palette (default: heat colors)
#' @param titles Optional vector of titles for each map
#' @param mask Optional logical mask for voxels
#' 
#' @return Invisible NULL, plots are displayed
#' @export
#' @importFrom graphics par image layout
plot_spatial_maps <- function(W, layout = NULL, zlim = NULL, 
                              col = heat.colors(100), titles = NULL, mask = NULL) {
  # Handle NeuroVol input
  if (is.list(W) && all(sapply(W, inherits, "NeuroVol"))) {
    # Extract matrices from NeuroVol objects
    W_mat <- sapply(W, as.vector)
  } else if (is.matrix(W)) {
    W_mat <- W
  } else {
    stop("W must be a matrix or list of NeuroVol objects")
  }
  
  V <- nrow(W_mat)
  K <- ncol(W_mat)
  
  # Apply mask if provided
  if (!is.null(mask)) {
    W_mat[!mask, ] <- NA
  }
  
  # Determine layout
  if (is.null(layout)) {
    layout <- c(ceiling(sqrt(K)), ceiling(K / ceiling(sqrt(K))))
  }
  
  # Set up plot layout
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = layout, mar = c(2, 2, 3, 1))
  
  # Determine color scale
  if (is.null(zlim)) {
    zlim <- range(W_mat, na.rm = TRUE)
  }
  
  # Compute matrix dimensions once
  dims <- try_reshape_vector(V)

  # Plot each spatial map
  for (k in 1:K) {
    map_matrix <- matrix(W_mat[, k], dims[1], dims[2])
    
    image(map_matrix, zlim = zlim, col = col, axes = FALSE,
          main = if (!is.null(titles)) titles[k] else paste("State", k))
  }
  
  invisible(NULL)
}

#' Plot State Time Course
#'
#' Visualizes temporal evolution of state activations.
#'
#' @param X State activation matrix (K x T)
#' @param type Plot type ("line", "heatmap", "stacked")
#' @param col Colors for each state
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param main Plot title
#' @param TR Repetition time for x-axis scaling
#' @param ... Additional arguments passed to plotting functions
#' 
#' @return ggplot object (if ggplot2 available) or NULL
#' @export
#' @importFrom graphics matplot image
plot_state_timecourse <- function(X, type = c("line", "heatmap", "stacked"),
                                  col = NULL, xlab = "Time", ylab = "Activation",
                                  main = "State Activations", TR = NULL, ...) {
  type <- match.arg(type)
  
  K <- nrow(X)
  T <- ncol(X)
  
  # Time axis
  if (!is.null(TR)) {
    time_axis <- (1:T) * TR
    xlab <- paste(xlab, "(s)")
  } else {
    time_axis <- 1:T
  }
  
  # Default colors
  if (is.null(col)) {
    col <- rainbow(K)
  }
  
  if (type == "line") {
    matplot(time_axis, t(X), type = "l", lty = 1, col = col,
            xlab = xlab, ylab = ylab, main = main, ...)
    legend("topright", legend = paste("State", 1:K), 
           col = col, lty = 1, bty = "n")
    
  } else if (type == "heatmap") {
    image(time_axis, 1:K, t(X), col = heat.colors(100),
          xlab = xlab, ylab = "State", main = main, ...)
    
  } else if (type == "stacked") {
    # Normalize to proportions for stacked plot
    X_prop <- X / colSums(X)
    X_prop[is.nan(X_prop)] <- 0
    
    # Cumulative sums for stacking
    X_cumsum <- apply(X_prop, 2, cumsum)
    
    plot(time_axis, rep(1, T), type = "n", ylim = c(0, 1),
         xlab = xlab, ylab = "Proportion", main = main, ...)
    
    # Plot stacked areas
    for (k in K:1) {
      y_top <- X_cumsum[k, ]
      y_bottom <- if (k > 1) X_cumsum[k-1, ] else rep(0, T)
      
      polygon(c(time_axis, rev(time_axis)), 
              c(y_top, rev(y_bottom)),
              col = col[k], border = NA)
    }
    
    legend("topright", legend = paste("State", 1:K), 
           fill = col, bty = "n")
  }
  
  invisible(NULL)
}

#' Plot Convergence
#'
#' Visualizes optimization convergence for both FISTA (CLD) and ELBO (CBD).
#'
#' @param values Vector of objective/ELBO values
#' @param type Type of values ("objective" or "elbo")
#' @param log_scale Logical, use log scale for y-axis
#' @param highlight_converged Iteration where convergence detected (optional)
#' @param ... Additional plotting arguments
#' 
#' @return NULL (plot displayed)
#' @export
#' @importFrom graphics plot lines abline
plot_convergence <- function(values, type = c("objective", "elbo"),
                             log_scale = FALSE, highlight_converged = NULL, ...) {
  type <- match.arg(type)
  
  if (length(values) == 0) {
    warning("No convergence values to plot")
    return(invisible(NULL))
  }
  
  n_iter <- length(values)
  
  # Set up plot
  ylab <- if (type == "objective") "Objective Value" else "ELBO"
  
  if (log_scale && all(values > 0)) {
    plot(1:n_iter, log(values), type = "l", 
         xlab = "Iteration", ylab = paste("log", ylab), ...)
  } else {
    plot(1:n_iter, values, type = "l",
         xlab = "Iteration", ylab = ylab, ...)
  }
  
  # Highlight convergence point
  if (!is.null(highlight_converged) && highlight_converged <= n_iter) {
    abline(v = highlight_converged, col = "red", lty = 2)
    text(highlight_converged, mean(range(values)), "Converged", 
         col = "red", pos = 4)
  }
  
  invisible(NULL)
}

#' Plot HRF
#'
#' Visualizes hemodynamic response functions.
#'
#' @param hrf HRF vector or matrix (multiple HRFs in columns)
#' @param TR Repetition time in seconds
#' @param col Colors for HRFs
#' @param main Plot title
#' @param ... Additional plotting arguments
#' 
#' @return NULL (plot displayed)
#' @export
plot_hrf <- function(hrf, TR = 2, col = NULL, main = "HRF", ...) {
  if (is.vector(hrf)) {
    hrf <- matrix(hrf, ncol = 1)
  }
  
  n_time <- nrow(hrf)
  n_hrf <- ncol(hrf)
  
  time_axis <- (0:(n_time-1)) * TR
  
  if (is.null(col)) {
    col <- if (n_hrf == 1) "black" else rainbow(n_hrf)
  }
  
  matplot(time_axis, hrf, type = "l", lty = 1, col = col,
          xlab = "Time (s)", ylab = "Response", main = main, ...)
  
  if (n_hrf > 1) {
    legend("topright", legend = paste("HRF", 1:n_hrf),
           col = col, lty = 1, bty = "n")
  }
  
  invisible(NULL)
}

#' Helper Function to Reshape Vector to Approximate Square
#'
#' @param n Length of vector
#' @return c(rows, cols) dimensions
#' @keywords internal
try_reshape_vector <- function(n) {
  sqrt_n <- sqrt(n)
  if (sqrt_n == floor(sqrt_n)) {
    return(c(sqrt_n, sqrt_n))
  }
  
  # Find factors closest to square root
  factors <- which(n %% 1:n == 0)
  diffs <- abs(factors - sqrt_n)
  best_factor <- factors[which.min(diffs)]
  
  c(best_factor, n / best_factor)
}

#' Create Diagnostic Plot Panel
#'
#' Creates a multi-panel diagnostic plot for decoder output.
#'
#' @param decoder_output Output from CLD or CBD
#' @param which Vector of plot numbers to display (1-4)
#' @param ... Additional arguments passed to individual plots
#' 
#' @return NULL (plots displayed)
#' @export
plot_diagnostics <- function(decoder_output, which = 1:4, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  n_plots <- length(which)
  layout_matrix <- matrix(c(1:n_plots, rep(0, 4-n_plots)), 2, 2, byrow = TRUE)
  layout(layout_matrix)
  
  for (i in which) {
    if (i == 1 && !is.null(decoder_output$W)) {
      # Spatial maps
      plot_spatial_maps(decoder_output$W, ...)
    } else if (i == 2 && !is.null(decoder_output$X_hat)) {
      # State timecourse
      plot_state_timecourse(decoder_output$X_hat, ...)
    } else if (i == 3 && !is.null(decoder_output$diagnostics$objective_values)) {
      # Convergence
      plot_convergence(decoder_output$diagnostics$objective_values, ...)
    } else if (i == 4 && !is.null(decoder_output$hrf)) {
      # HRF
      plot_hrf(decoder_output$hrf, ...)
    }
  }
  
  invisible(NULL)
}
