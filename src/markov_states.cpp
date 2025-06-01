#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Generate discrete Markov states (C++)
//'
//' Samples a sequence of states from a first-order Markov chain
//' defined by `Pi` and `pi0` and returns a one-hot encoded matrix.
//'
//' @param Pi   Transition matrix (K x K)
//' @param pi0  Initial state probabilities (length K)
//' @param T_len Number of time points to generate
//'
//' @return K x T matrix of one-hot encoded states
//' @keywords internal
// [[Rcpp::export]]
arma::mat generate_markov_states_cpp(const arma::mat& Pi,
                                     const arma::vec& pi0,
                                     int T_len) {
  if (Pi.is_empty() || pi0.is_empty()) {
    stop("Pi and pi0 must not be empty");
  }
  int K = Pi.n_rows;
  if (Pi.n_cols != K) {
    stop("Pi must be square");
  }
  if (pi0.n_elem != (unsigned)K) {
    stop("pi0 must have length K");
  }
  if (T_len <= 0) {
    stop("T must be positive");
  }

  RNGScope scope; // ensure RNG initialised

  arma::uvec states(T_len);
  // sample first state
  double u = R::runif(0.0, 1.0);
  double cum = 0.0;
  for (int k = 0; k < K; ++k) {
    cum += pi0[k];
    if (u <= cum) {
      states[0] = k;
      break;
    }
  }

  // subsequent states
  for (int t = 1; t < T_len; ++t) {
    u = R::runif(0.0, 1.0);
    cum = 0.0;
    for (int k = 0; k < K; ++k) {
      cum += Pi(states[t-1], k);
      if (u <= cum) {
        states[t] = k;
        break;
      }
    }
  }

  arma::mat S(K, T_len, fill::zeros);
  for (int t = 0; t < T_len; ++t) {
    S(states[t], t) = 1.0;
  }

  return S;
}

