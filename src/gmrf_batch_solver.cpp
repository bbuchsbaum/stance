#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#ifdef _OPENMP
#include <omp.h>
#endif

using Eigen::MatrixXd;
using Eigen::SparseMatrix;

// [[Rcpp::export]]
MatrixXd solve_gmrf_batched(const SparseMatrix<double>& XtX,
                             const SparseMatrix<double>& L_gmrf,
                             const MatrixXd& XtY,
                             double lambda_h,
                             int block_size = 64,
                             double tol = 1e-6) {
  int V = XtY.rows();
  int L_basis = XtY.cols();
  int n_blocks = (V + block_size - 1) / block_size;

  MatrixXd H_v(V, L_basis);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int b = 0; b < n_blocks; ++b) {
    int start = b * block_size;
    int end = std::min(start + block_size, V);
    int block_v = end - start;

    SparseMatrix<double> XtX_block = XtX.block(start, start, block_v, block_v);
    SparseMatrix<double> L_block = L_gmrf.block(start, start, block_v, block_v);

    SparseMatrix<double> Q_block = XtX_block + lambda_h * L_block;

    Eigen::SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(Q_block);

    if (solver.info() != Eigen::Success) {
      // Retry with diagonal regularization for numerical stability
      SparseMatrix<double> I(block_v, block_v);
      I.setIdentity();
      solver.compute(Q_block + tol * I);
      if (solver.info() != Eigen::Success) {
        Rcpp::stop("Cholesky failed in block %d", b);
      }
    }

    MatrixXd RHS = XtY.block(start, 0, block_v, L_basis);
    MatrixXd H_block = solver.solve(RHS);

    H_v.block(start, 0, block_v, L_basis) = H_block;
  }

  return H_v;
}

