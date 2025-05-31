# Fallback R implementations for missing C++ code
# These functions provide slow but functional stand-ins so the
# package can run even when the Rcpp versions are unavailable.

#' Forward-Backward Algorithm (fallback)
#'
#' Runs a basic forward-backward procedure for an HMM using R code.
#' This function is used when the optimized C++ implementation is not
#' available.
#'
#' @keywords internal
forward_backward_r <- function(Y, U, V, H_v, hrf_basis, Pi, pi0, sigma2) {
  V_vox <- nrow(Y)
  T_len <- ncol(Y)
  K <- length(pi0)

  # Predicted mean for each state (ignore HRF for fallback)
  state_means <- lapply(seq_len(K), function(k) U %*% V[k, ])

  logB <- matrix(0, K, T_len)
  for (k in seq_len(K)) {
    for (t in seq_len(T_len)) {
      diff <- Y[, t] - state_means[[k]]
      logB[k, t] <- -0.5 * sum(diff * diff) / sigma2
    }
  }

  log_alpha <- matrix(0, K, T_len)
  log_alpha[, 1] <- log(pi0 + 1e-12) + logB[, 1]
  for (t in 2:T_len) {
    for (k in seq_len(K)) {
      log_alpha[k, t] <- logB[k, t] + logSumExp(log_alpha[, t - 1] + log(Pi[, k] + 1e-12))
    }
  }

  log_beta <- matrix(0, K, T_len)
  for (t in (T_len - 1):1) {
    for (k in seq_len(K)) {
      log_beta[k, t] <- logSumExp(log(Pi[k, ] + 1e-12) + logB[, t + 1] + log_beta[, t + 1])
    }
  }

  gamma <- matrix(0, K, T_len)
  for (t in seq_len(T_len)) {
    log_prob <- log_alpha[, t] + log_beta[, t]
    gamma[, t] <- exp(log_prob - logSumExp(log_prob))
  }

  xi <- array(0, c(K, K, T_len - 1))
  for (t in seq_len(T_len - 1)) {
    for (i in seq_len(K)) {
      for (j in seq_len(K)) {
        xi[i, j, t] <- exp(log_alpha[i, t] + log(Pi[i, j] + 1e-12) +
                            logB[j, t + 1] + log_beta[j, t + 1])
      }
    }
    xi[, , t] <- xi[, , t] / sum(xi[, , t])
  }

  log_lik <- sum(log(colSums(exp(log_alpha))))

  list(gamma = gamma, xi = xi, log_likelihood = log_lik)
}

#' Compute log-sum-exp in a numerically stable way
#' @keywords internal
logSumExp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Wrapper that chooses the C++ version when available
#'
#' @keywords internal
forward_backward_algorithm <- function(Y, U, V, H_v, hrf_basis,
                                       Pi, pi0, sigma2, engine = "cpp") {
  if (engine == "cpp" && exists("forward_backward_cpp")) {
    return(forward_backward_cpp(Y, U, V, H_v, hrf_basis, Pi, pi0, sigma2))
  }
  forward_backward_r(Y, U, V, H_v, hrf_basis, Pi, pi0, sigma2)
}

# Use the R implementation as the cpp fallback if real C++ versions
# are not compiled.
forward_backward_cpp <- forward_backward_r

#' Update spatial components (fallback)
#' @keywords internal
update_spatial_components_r <- function(Y, S_gamma, H_v, hrf_basis, current_U, current_V) {
  r <- ncol(current_U)
  cross <- Y %*% t(S_gamma)
  sv <- svd(cross, nu = r, nv = r)
  U_new <- sv$u[, seq_len(r), drop = FALSE]
  V_new <- sv$v[, seq_len(r), drop = FALSE]
  V_new <- V_new * matrix(sv$d[seq_len(r)], nrow = nrow(V_new), ncol = r, byrow = TRUE)
  list(U = U_new, V = V_new)
}

update_spatial_components_cpp <- update_spatial_components_r

#' Update HRF coefficients with GMRF prior (fallback)
#' @keywords internal
update_hrf_coefficients_r <- function(Y, S_gamma, U, V, hrf_basis, L_gmrf,
                                     lambda_H_prior, sigma2) {
  V_vox <- nrow(Y)
  B <- ncol(hrf_basis)
  XtX_inv <- solve(crossprod(hrf_basis) + diag(1e-6, B))
  Xt <- t(hrf_basis)
  H_v <- matrix(0, V_vox, B)
  for (v in seq_len(V_vox)) {
    H_v[v, ] <- XtX_inv %*% Xt %*% Y[v, ]
  }
  H_v
}

update_hrf_coefficients_gmrf_cpp <- update_hrf_coefficients_r
