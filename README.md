# stance

`stance` provides continuous decoding tools for fMRI analysis. Two algorithms are available:

* **Continuous Linear Decoder (CLD)** – deterministic GLM+SVD approach with FISTA optimization.
* **Continuous Bayesian Decoder (CBD)** – probabilistic model using a hidden Markov prior and variational Bayes.

Use **CLD** when you have a reliable design matrix and want fast estimation. Use **CBD** when state boundaries are unknown or you need uncertainty estimates for latent states. See `vignette("getting-started", package = "stance")` for worked examples.

## Parallel Processing

Several low-level routines such as row-wise convolution and total variation denoising are implemented in C++ with optional OpenMP support. When compiled with OpenMP, these functions can run in parallel across rows and the number of threads can be controlled via the `n_threads` argument (use `0` to auto-detect).
