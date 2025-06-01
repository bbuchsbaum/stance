#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//' Batched FFT convolution for voxel-specific HRFs
//'
//' Performs convolution of state regressors with voxel-specific HRFs
//' using FFTs. The design matrix has states in rows and time in
//' columns. HRFs are supplied per voxel. The result is returned as a
//' cube with dimensions (V x K x T).
//'
//' @param design Matrix of regressors (K x T)
//' @param hrfs   Matrix of HRF kernels (V x L)
//' @param n_threads Number of threads (0 = all available)
//'
//' @return Cube of convolved regressors (V x K x T)
//' @keywords internal
// [[Rcpp::export]]
arma::cube convolve_voxel_hrf_fft_rcpp(const arma::mat& design,
                                       const arma::mat& hrfs,
                                       int n_threads = 0) {
  if (design.is_empty() || hrfs.is_empty()) {
    stop("design and hrfs cannot be empty");
  }

  int K = design.n_rows;
  int T = design.n_cols;
  int V = hrfs.n_rows;
  int L = hrfs.n_cols;

  int n_fft = 1;
  while (n_fft < T + L - 1) n_fft <<= 1;

#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n_threads <= 0 || n_threads > max_threads) n_threads = max_threads;
  omp_set_num_threads(n_threads);
#endif

  // Precompute FFT of regressors
  arma::cx_mat design_fft(K, n_fft);
  for (int k = 0; k < K; ++k) {
    arma::vec x_pad(n_fft, fill::zeros);
    x_pad.head(T) = design.row(k).t();
    design_fft.row(k) = arma::fft(x_pad).t();
  }

  arma::cube result(V, K, T, fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int v = 0; v < V; ++v) {
    arma::vec h_pad(n_fft, fill::zeros);
    h_pad.head(L) = hrfs.row(v).t();
    arma::cx_vec h_fft = arma::fft(h_pad);

    for (int k = 0; k < K; ++k) {
      arma::cx_vec Y_fft = design_fft.row(k).t() % h_fft;
      arma::vec y = arma::real(arma::ifft(Y_fft));
      result.tube(v, k) = y.head(T);
    }
  }

  return result;
}

