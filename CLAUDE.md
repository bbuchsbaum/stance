# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`stance` is an R package implementing a Continuous Bayesian Decoder (CBD) for fMRI analysis using variational Bayes inference. It provides spatial pattern learning and temporal state inference with voxel-specific HRF estimation.

## Common Development Commands

### Build and Check
```bash
# Build the package
R CMD build .

# Check the package (replace with actual version)
R CMD check stance_*.tar.gz

# Install the package
R CMD INSTALL .
```

### Testing
```r
# Run all tests
testthat::test_package("stance")

# Run a specific test file
testthat::test_file("tests/testthat/test-basic.R")
```

### Documentation
```r
# Generate documentation from roxygen2 comments
roxygen2::roxygenise()

# Build vignettes
devtools::build_vignettes()
```

### Development Workflow
```r
# Load package for interactive development
devtools::load_all()

# Check package
devtools::check()

# Run tests interactively
devtools::test()
```

## Architecture

The package uses an R6 class-based architecture centered around the `ContinuousBayesianDecoder` class, which implements:

1. **Spatial Decomposition**: Low-rank factorization of voxel patterns (U, V matrices)
2. **Temporal Modeling**: Hidden Markov Model for state transitions with learnable transition matrix
3. **HRF Modeling**: Voxel-specific hemodynamic response functions using spline basis
4. **Inference**: Variational Bayes updates for all model parameters

Key dependencies:
- `neuroim2`: Provides neuroimaging data structures (`NeuroVec`, `NeuroVol`)
- `fmrireg`: fMRI regression modeling and HRF utilities
- `Rcpp`/`RcppArmadillo`: C++ integration for performance-critical computations

### neuroim2 Usage

The package relies heavily on neuroim2 for neuroimaging data handling:

- **Input Data**: Expects `NeuroVec` objects (4D fMRI data) with proper `NeuroSpace` metadata
- **Masks**: Uses `LogicalNeuroVol` for brain masks to define analysis regions
- **Sparse Data**: Leverages `SparseNeuroVec` for memory-efficient masked data storage
- **ROI Operations**: May use searchlight functions (`spherical_roi`, `searchlight`) for local pattern analysis
- **Key Classes**:
  - `NeuroVec`: 4D fMRI time series data
  - `NeuroVol`: 3D volumetric data (for spatial patterns)
  - `LogicalNeuroVol`: Binary masks
  - `NeuroSpace`: Spatial reference information (dimensions, voxel size, orientation)

### fmrireg Usage

The package leverages fmrireg for HRF modeling and convolution operations:

- **HRF Basis Functions**: Uses pre-defined HRF basis sets:
  - `HRF_SPMG1`: SPM canonical HRF (single basis) - default
  - `HRF_SPMG2`: SPM canonical + temporal derivative (2 basis)
  - `HRF_BSPLINE`: B-spline basis (5 basis functions)
  - `HRF_FIR`: Finite impulse response basis (12 basis)
- **Convolution**: Leverages fmrireg's FFT-based convolution utilities for efficient HRF convolution
- **Event Models**: May use `event_model()` for design matrix construction during simulation
- **Key Functions**:
  - `gen_hrf()`: Generate HRF with specific parameters
  - `evaluate()`: Evaluate HRF at specific time points
  - `sampling_frame()`: Define temporal sampling grid for HRF evaluation

## Important Notes

- The package uses `renv` for dependency management. Run `renv::restore()` to set up the environment.
- C++ implementations now exist in `src/`, including `fista_solver.cpp` for optimization and `log_likelihood.cpp` for computing likelihoods.
- The main implementation is in `R/continuous_bayesian_decoder.R`.
- Sprint planning documents in `data-raw/` outline the development roadmap.
- When implementing new features, follow the existing R6 method structure and update the variational Bayes inference accordingly.