#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Helper used by the CBD implementation to convolve state regressors
// with voxel-specific HRFs generated via \pkg{fmrireg}.  Accepts
// matrices or data extracted from `neuroim2` objects and supports
// optional OpenMP parallelism.

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


// Forward declaration of the FFT helper implemented in
// convolution_fft.cpp.  This allows us to reuse the optimised FFT path
// when the time series is long enough.
arma::cube convolve_voxel_hrf_fft_rcpp(const arma::mat& design,
                                       const arma::mat& hrfs,
                                       int n_threads = 0);

//' Convolve design matrix with voxel-specific HRFs
//'
//' Each voxel can have a different HRF kernel. The design matrix has
//' state regressors in rows. This function performs convolution of each
//' state's regressor with the corresponding voxel HRF for all voxels.
//'

//' Uses direct time-domain convolution by default but switches to an
//' FFT-based implementation when the number of time points exceeds
//' `fft_threshold`.  This mirrors the behaviour of the higher level R
//' helpers and improves performance for long time series.  The
//' threshold can be increased in unit tests to force the direct path.
//' When compiled with OpenMP the convolution loops run in parallel.
//'
//' @param design Matrix of regressors (K x T)
//' @param hrfs   Matrix of HRF kernels (V x L)
//' @param fft_threshold Integer threshold for using FFT;
//'   FFT is used when the number of time points exceeds this value
//' @param n_threads Number of OpenMP threads (0 = all available)
//'
//' @return A cube with dimensions (V x K x T)
//'
//' @keywords internal
// [[Rcpp::export]]
arma::cube convolve_voxel_hrf_rcpp(const arma::mat& design,
                                   const arma::mat& hrfs,
                                   int fft_threshold = 256,
                                   int n_threads = 0) {
  if (design.is_empty() || hrfs.is_empty()) {
    stop("design and hrfs cannot be empty");
  }

  int K = design.n_rows;
  int T = design.n_cols;
  int V = hrfs.n_rows;
  int L = hrfs.n_cols;

#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n_threads <= 0 || n_threads > max_threads) {
    n_threads = max_threads;
  }
  omp_set_num_threads(n_threads);
#else
  (void)n_threads; // suppress warning
#endif

  bool use_fft = (T > fft_threshold);
  if (use_fft) {
    // Delegate to the specialised FFT implementation for speed
    return convolve_voxel_hrf_fft_rcpp(design, hrfs, n_threads);
  }

  // Output cube for direct convolution

  arma::cube result(V, K, T, fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int v = 0; v < V; ++v) {
    const arma::rowvec h = hrfs.row(v);
    for (int k = 0; k < K; ++k) {
      const arma::rowvec s = design.row(k);
      arma::vec conv(T, fill::zeros);
      for (int t = 0; t < T; ++t) {
        double sum = 0.0;
        int max_tau = std::min(L, t + 1);
        for (int tau = 0; tau < max_tau; ++tau) {
          int idx = t - tau;
          if (idx >= 0) {
            sum += s[idx] * h[tau];
          }
        }
        conv[t] = sum;
      }
      for (int t = 0; t < T; ++t) {
        result(v, k, t) = conv[t];
      }
    }
  }

  return result;
}

