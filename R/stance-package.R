#' @title stance: Continuous Decoders for fMRI Analysis
#'
#' @description
#' The stance package provides two related algorithms for continuous fMRI decoding:
#' the fast \emph{Continuous Linear Decoder} (CLD) and the probabilistic
#' \emph{Continuous Bayesian Decoder} (CBD). Both build on established R
#' neuroimaging infrastructure (neuroim2, fmrireg) for efficient spatial pattern
#' learning and temporal state inference with optional voxel-specific HRF
#' estimation and spatial smoothing priors.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{ContinuousLinearDecoder}}}{R6 class implementing the deterministic CLD algorithm}
#'   \item{\code{\link{ContinuousBayesianDecoder}}}{Main R6 class for CBD analysis}
#'   \item{\code{\link{run_cbd_analysis}}}{High-level wrapper for CBD}
#'   \item{\code{\link{qc_report}}}{Generate diagnostic reports}
#'   \item{\code{\link{simulate_fmri_data}}}{Generate synthetic data for testing}
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
#'   \item Fast deterministic CLD algorithm with GLM+SVD and FISTA optimisation
#'   \item Variational Bayes inference for CBD with uncertainty estimates
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
