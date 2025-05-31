# stance

TODO: Describe what your package does.

## Parallel Processing

Several low-level routines such as row-wise convolution and total variation
denoising are implemented in C++ with optional OpenMP support. When compiled
with OpenMP, these functions can run in parallel across rows and the number of
threads can be controlled via the `n_threads` argument (use `0` to auto-detect).
