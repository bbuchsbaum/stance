#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Batched HRF Coefficient Update
//'
//' Estimates HRF coefficients for each voxel using batched least squares.
//' The design matrix cross-product is solved with `solve_sympd` when
//' possible and falls back to a general solver with ridge regularisation
//' if the system is not positive definite.
//'
//' @param Y Data matrix (V x T)
//' @param hrf_basis HRF basis matrix (T x L)
//' @param ridge Small ridge penalty added to the diagonal
//' @param block_size Number of voxels per batch
//'
//' @return Matrix of HRF coefficients (V x L)
//' @keywords internal
// [[Rcpp::export]]
arma::mat update_hrf_coefficients_batched_cpp(const arma::mat& Y,
                                              const arma::mat& hrf_basis,
                                              double ridge = 1e-6,
                                              int block_size = 64) {
  if (Y.is_empty() || hrf_basis.is_empty()) {
    stop("Y and hrf_basis must not be empty");
  }
  if (hrf_basis.n_rows != Y.n_cols) {
    stop("hrf_basis rows must match number of time points");
  }

  int V = Y.n_rows;
  int B = hrf_basis.n_cols;
  arma::mat H(V, B, fill::zeros);

  arma::mat XtX = hrf_basis.t() * hrf_basis;
  XtX.diag() += ridge;

  arma::mat XtX_inv;
  bool ok = arma::solve_sympd(XtX_inv, XtX, arma::eye(B,B));
  if (!ok) {
    arma::mat XtX_reg = XtX + ridge * arma::eye(B,B);
    ok = arma::solve(XtX_inv, XtX_reg, arma::eye(B,B));
    if (!ok) {
      XtX_inv = arma::inv(XtX_reg);
    }
  }

  for (int start = 0; start < V; start += block_size) {
    int end = std::min(start + block_size - 1, V - 1);
    int n = end - start + 1;
    arma::mat Yblock = Y.rows(start, end).t(); // T x n
    arma::mat XTy = hrf_basis.t() * Yblock;    // B x n
    arma::mat Hblock = XtX_inv * XTy;          // B x n
    H.rows(start, end) = Hblock.t();
  }

  return H;
}

