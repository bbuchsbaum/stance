#' Data Structure Utilities for stance
#'
#' Helper functions for seamless conversion between neuroim2 objects and matrices.
#' These utilities provide a bridge between neuroimaging data structures and
#' computational algorithms, ensuring spatial metadata preservation throughout
#' the analysis pipeline.
#'
#' @name data_structures
NULL

#' Validate fMRI Input Data
#'
#' Validates and converts fMRI data to matrix format, checking dimensions
#' and data integrity. Handles both standard matrices and neuroim2 objects.
#'
#' @param Y Input data (matrix, NeuroVec, or NeuroVol object)
#' @param expected_dims Optional named list with 'V' (voxels) and/or 'T' (timepoints)
#' @param check_finite Logical, whether to check for finite values (default TRUE)
#' @param verbose Logical, whether to print diagnostic information (default FALSE)
#' 
#' @return A list with components:
#'   \item{data}{V x T matrix of fMRI data}
#'   \item{type}{Character string indicating input type}
#'   \item{space}{NeuroSpace object if available, NULL otherwise}
#'   \item{mask}{Logical mask if available, NULL otherwise}
#'   \item{dims}{Named vector with V and T dimensions}
#' 
#' @export
#' @examples
#' \dontrun{
#' # With matrix input
#' Y_mat <- matrix(rnorm(1000 * 200), 1000, 200)
#' validated <- validate_fmri_input(Y_mat)
#' 
#' # With NeuroVec input
#' Y_nvec <- neuroim2::NeuroVec(Y_mat, space_info)
#' validated <- validate_fmri_input(Y_nvec)
#' }
validate_fmri_input <- function(Y, expected_dims = NULL, check_finite = TRUE, verbose = FALSE) {
  # Initialize return structure
  result <- list(
    data = NULL,
    type = NULL,
    space = NULL,
    mask = NULL,
    dims = NULL
  )
  
  # Handle NeuroVec input
  if (inherits(Y, "NeuroVec")) {
    result$type <- "NeuroVec"
    result$data <- as.matrix(Y)
    result$space <- neuroim2::space(Y)
    
    # Extract mask if it's a SparseNeuroVec
    if (inherits(Y, "SparseNeuroVec")) {
      result$mask <- Y@mask
      if (verbose) message("Input is SparseNeuroVec with ", sum(result$mask), " active voxels")
    }
    
  # Handle NeuroVol input (3D, needs to be reshaped)
  } else if (inherits(Y, "NeuroVol")) {
    stop("NeuroVol input not supported. Please provide 4D data as NeuroVec or matrix.")
    
  # Handle matrix input
  } else if (is.matrix(Y)) {
    result$type <- "matrix"
    result$data <- Y
    if (verbose) message("Input is standard matrix")
    
  # Handle data.frame that can be converted
  } else if (is.data.frame(Y)) {
    result$type <- "data.frame"
    # ensure all columns are numeric before conversion
    non_numeric_cols <- !vapply(Y, is.numeric, logical(1))
    if (any(non_numeric_cols)) {
      stop("data.frame contains non-numeric columns: ",
           paste(names(Y)[non_numeric_cols], collapse = ", "))
    }
    result$data <- as.matrix(Y, mode = "numeric")
    if (verbose) message("Converted data.frame to numeric matrix")
    
  } else {
    stop("Input must be a matrix, NeuroVec, or convertible data.frame")
  }
  
  # Get dimensions
  result$dims <- c(V = nrow(result$data), T = ncol(result$data))
  
  # Validate dimensions if expected
  if (!is.null(expected_dims)) {
    if (!is.null(expected_dims$V) && result$dims["V"] != expected_dims$V) {
      stop(sprintf("Expected %d voxels but got %d", expected_dims$V, result$dims["V"]))
    }
    if (!is.null(expected_dims$T) && result$dims["T"] != expected_dims$T) {
      stop(sprintf("Expected %d timepoints but got %d", expected_dims$T, result$dims["T"]))
    }
  }
  
  # Check for finite values
  if (check_finite) {
    non_finite <- !is.finite(result$data)
    if (any(non_finite)) {
      n_na <- sum(is.na(result$data[non_finite]))
      n_inf <- sum(is.infinite(result$data[non_finite]))
      stop(sprintf("Data contains non-finite values: %d NA, %d Inf/-Inf", n_na, n_inf))
    }
  }
  
  # Informative message
  if (verbose) {
    message(sprintf("Validated %s: %d voxels × %d timepoints", 
                    result$type, result$dims["V"], result$dims["T"]))
  }
  
  result
}

#' Extract Data Matrix from Neuroimaging Objects
#'
#' Converts various neuroimaging objects to matrix format while preserving
#' metadata for later reconstruction.
#'
#' @param obj NeuroVec, NeuroVol, or list of neuroimaging objects
#' @param flatten_space Logical, whether to flatten spatial dimensions (default TRUE)
#' @param preserve_attributes Logical, whether to preserve object attributes (default TRUE)
#' @param force_matrix Logical, if TRUE always return a matrix by flattening
#'   3-D data (default FALSE)
#' 
#' @return A list with components:
#'   \item{data}{Numeric matrix. If \code{flatten_space = FALSE} and
#'     \code{force_matrix = FALSE} for 3-D input, a 3-D array is returned}
#'   \item{metadata}{List of preserved metadata for reconstruction}
#' 
#' @export
extract_data_matrix <- function(obj, flatten_space = TRUE,
                                preserve_attributes = TRUE,
                                force_matrix = FALSE) {
  metadata <- list()
  
  # NeuroVec (4D: space x time)
  if (inherits(obj, "NeuroVec")) {
    data <- as.matrix(obj)
    
    if (preserve_attributes) {
      metadata$class <- class(obj)
      metadata$space <- neuroim2::space(obj)
      metadata$mask <- if (inherits(obj, "SparseNeuroVec")) obj@mask else NULL
    }
    
  # NeuroVol (3D: space only)
  } else if (inherits(obj, "NeuroVol")) {
    if (flatten_space || force_matrix) {
      data <- matrix(as.vector(obj), ncol = 1)
      metadata$original_dim <- dim(obj)
    } else {
      data <- as.array(obj)
    }
    
    if (preserve_attributes) {
      metadata$class <- class(obj)
      metadata$space <- neuroim2::space(obj)
    }
    
  # List of NeuroVols (e.g., multiple spatial maps)
  } else if (is.list(obj) && all(sapply(obj, inherits, "NeuroVol"))) {
    data <- sapply(obj, as.vector)
    
    if (preserve_attributes) {
      metadata$class <- "list_NeuroVol"
      metadata$space <- neuroim2::space(obj[[1]])
      metadata$n_volumes <- length(obj)
    }
    
  # Already a matrix
  } else if (is.matrix(obj)) {
    data <- obj
    metadata$class <- "matrix"
    
  } else {
    stop("Unsupported object type: ", paste(class(obj), collapse = ", "))
  }

  if (force_matrix && !is.matrix(data)) {
    metadata$original_dim <- dim(data)
    data <- matrix(as.vector(data), ncol = 1)
  }

  list(data = data, metadata = metadata)
}

#' Restore Spatial Structure to Data
#'
#' Converts matrix data back to neuroimaging objects with preserved spatial metadata.
#' Supports both deterministic (CLD) and probabilistic (CBD) outputs.
#'
#' @param mat Data matrix (V x T for time series, V x K for spatial maps)
#' @param reference Reference neuroimaging object or metadata list
#' @param output_type Character, type of output ("temporal" or "spatial")
#' @param probabilistic Logical, whether data contains uncertainty estimates
#' 
#' @return NeuroVec, NeuroVol, or list of neuroimaging objects
#' 
#' @export
restore_spatial_structure <- function(mat, reference, output_type = c("temporal", "spatial"), 
                                      probabilistic = FALSE) {
  output_type <- match.arg(output_type)
  
  # Extract space information
  if (inherits(reference, c("NeuroVec", "NeuroVol"))) {
    space_obj <- neuroim2::space(reference)
    mask <- if (inherits(reference, "SparseNeuroVec")) reference@mask else NULL
  } else if (is.list(reference) && !is.null(reference$space)) {
    space_obj <- reference$space
    mask <- reference$mask
  } else if (is.list(reference) && !is.null(reference$metadata$space)) {
    space_obj <- reference$metadata$space
    mask <- reference$metadata$mask
  } else {
    stop("Cannot extract spatial information from reference")
  }
  
  # Determine expected voxel count from mask or space
  expected_voxels <- NULL
  if (!is.null(mask)) {
    expected_voxels <- sum(mask)
  } else {
    space_dim <- try(neuroim2::dim(space_obj), silent = TRUE)
    if (!inherits(space_dim, "try-error")) {
      if (length(space_dim) > 3) {
        expected_voxels <- prod(space_dim[1:3])
      } else {
        expected_voxels <- prod(space_dim)
      }
    }
  }

  if (!is.null(expected_voxels) && nrow(mat) != expected_voxels) {
    stop(sprintf(
      "Matrix has %d rows but spatial reference implies %d voxels",
      nrow(mat), expected_voxels
    ))
  }

  # Restore based on output type
  if (output_type == "temporal") {
    # Time series data (V x T)
    if (!is.null(mask)) {
      # Create SparseNeuroVec
      result <- neuroim2::SparseNeuroVec(mat, space_obj, mask)
    } else {
      # Create DenseNeuroVec
      result <- neuroim2::NeuroVec(mat, space_obj)
    }
    
  } else {  # spatial
    # Spatial maps (V x K)
    if (ncol(mat) == 1) {
      # Single map
      result <- neuroim2::NeuroVol(mat[, 1], space_obj)
    } else {
      # Multiple maps
      result <- lapply(seq_len(ncol(mat)), function(k) {
        neuroim2::NeuroVol(mat[, k], space_obj)
      })
    }
  }
  
  # Add uncertainty information if probabilistic
  if (probabilistic) {
    if (is.list(result)) {
      result <- lapply(result, function(elem) {
        attr(elem, "probabilistic") <- TRUE
        elem
      })
      attr(result, "probabilistic") <- TRUE
    } else {
      attr(result, "probabilistic") <- TRUE
    }
  }
  
  result
}

#' Check Temporal Alignment
#'
#' Validates temporal alignment between fMRI data and experimental design.
#' Ensures TR consistency and appropriate temporal coverage.
#'
#' @param Y_info Validated fMRI data info from validate_fmri_input()
#' @param design Design matrix or event timing information
#' @param TR Repetition time in seconds
#' @param tolerance Numeric, tolerance for timing mismatches (default 0.1)
#' 
#' @return Logical TRUE if aligned, otherwise stops with informative error
#' 
#' @export
check_temporal_alignment <- function(Y_info, design, TR, tolerance = 0.1) {
  n_timepoints <- Y_info$dims["T"]
  
  # Check design matrix dimensions
  if (is.matrix(design)) {
    if (ncol(design) != n_timepoints) {
      stop(sprintf("Design matrix has %d columns but data has %d timepoints", 
                   ncol(design), n_timepoints))
    }
  }
  
  # Check event timing if provided as data.frame
  if (is.data.frame(design) && all(c("onset", "duration") %in% names(design))) {
    max_time <- max(design$onset + design$duration)
    data_duration <- n_timepoints * TR
    
    if (max_time > data_duration + tolerance) {
      stop(sprintf("Events extend to %.1fs but data only covers %.1fs", 
                   max_time, data_duration))
    }
  }
  
  # Check TR consistency if space object available
  if (!is.null(Y_info$space)) {
    # Try to extract TR from NeuroSpace if it has temporal information
    space_tr <- try({
      axes <- neuroim2::axes(Y_info$space)
      if (length(axes) >= 4) {
        neuroim2::spacing(Y_info$space)[4]
      } else {
        NULL
      }
    }, silent = TRUE)
    
    if (!inherits(space_tr, "try-error") && !is.null(space_tr)) {
      if (abs(space_tr - TR) > tolerance) {
        warning(sprintf("TR mismatch: data suggests %.3fs but analysis uses %.3fs", 
                        space_tr, TR))
      }
    }
  }
  
  TRUE
}

#' Create Neuroimaging-Compatible Output List
#'
#' Helper function to create standardized output format that can work
#' with both CLD and CBD algorithms.
#'
#' @param spatial_maps Matrix of spatial maps (V x K)
#' @param temporal_activations Matrix of temporal activations (K x T)
#' @param metadata List of metadata from original data
#' @param algorithm Character, "CLD" or "CBD"
#' @param additional_info List of algorithm-specific information
#' 
#' @return List with standardized structure
#' 
#' @export
create_output_structure <- function(spatial_maps, temporal_activations,
                                    metadata, algorithm = c("CLD", "CBD"),
                                    additional_info = list()) {
  algorithm <- match.arg(algorithm)

  if (!is.matrix(spatial_maps) || !is.numeric(spatial_maps)) {
    stop("spatial_maps must be a numeric matrix")
  }

  if (!is.matrix(temporal_activations) || !is.numeric(temporal_activations)) {
    stop("temporal_activations must be a numeric matrix")
  }

  if (nrow(temporal_activations) != ncol(spatial_maps)) {
    stop(sprintf(
      "Number of rows in temporal_activations (%d) must equal number of columns in spatial_maps (%d)",
      nrow(temporal_activations), ncol(spatial_maps)
    ))
  }
  
  output <- list(
    algorithm = algorithm,
    spatial_maps = spatial_maps,
    temporal_activations = temporal_activations,
    metadata = metadata,
    dims = list(
      V = nrow(spatial_maps),
      K = ncol(spatial_maps),
      T = ncol(temporal_activations)
    )
  )
  
  # Add algorithm-specific fields
  if (algorithm == "CBD") {
    output$uncertainty <- additional_info$uncertainty
    output$posterior_params <- additional_info$posterior_params
  } else {  # CLD
    output$objective_values <- additional_info$objective_values
    output$convergence_info <- additional_info$convergence_info
  }
  
  class(output) <- c(paste0(algorithm, "_output"), "stance_output")
  output
}
