#' @title stance: Continuous Bayesian Decoder for fMRI Analysis
#'
#' @description
#' The stance package provides a Continuous Bayesian Decoder (CBD) for analyzing
#' fMRI data using variational Bayes inference. Built on established R neuroimaging
#' infrastructure (neuroim2, fmrireg), it offers efficient spatial pattern learning
#' and temporal state inference with voxel-specific HRF estimation and spatial
#' smoothing priors.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{ContinuousBayesianDecoder}}}{Main R6 class for CBD analysis}
#'   \item{\code{\link{run_cbd_analysis}}}{High-level wrapper function}
#'   \item{\code{\link{qc_report}}}{Generate diagnostic reports}
#'   \item{\code{\link{simulate_cbd_data}}}{Generate synthetic data for testing}
#' }
#'
#' @section Integration with Neuroimaging Packages:
#' The package seamlessly integrates with:
#' \describe{
#'   \item{neuroim2}{For spatial data structures (NeuroVol, NeuroVec, NeuroSpace)}
#'   \item{fmrireg}{For HRF basis functions and regression components}
#' }
#'
#' @section Key Features:
#' \itemize{
#'   \item Variational Bayes inference for efficient computation
#'   \item Voxel-specific HRF estimation using established basis functions
#'   \item Spatial smoothing priors (GMRF) for robust parameter estimation
#'   \item Low-rank spatial decompositions for computational efficiency
#'   \item Native support for BIDS-formatted neuroimaging data
#'   \item Comprehensive diagnostic and visualization tools
#' }
#'
#' @docType package
#' @name stance-package
#' @aliases stance
#' @useDynLib stance, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Suppress R CMD check notes for global variables used in data.table/dplyr operations
utils::globalVariables(c(".", "time", "voxel", "state", "value", "iteration", "elbo"))
