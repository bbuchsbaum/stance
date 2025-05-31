# stance

Internal helper functions `fast_matmul`, `convolve_rows`, and
`process_blocks` were removed as they were unused. This keeps the code base
leaner while retaining the main API.


## Parallel Processing

Several low-level routines such as row-wise convolution and total variation
denoising are implemented in C++ with optional OpenMP support. When compiled
with OpenMP, these functions can run in parallel across rows and the number of
threads can be controlled via the `n_threads` argument (use `0` to auto-detect).
