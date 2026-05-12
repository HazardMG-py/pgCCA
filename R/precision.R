# =============================================================================
# pgCCA — Stage 1: Precision Matrix Estimation
# File:   R/precision.R
# =============================================================================

#' Ledoit-Wolf Precision Matrix Estimator
#'
#' Analytic shrinkage toward scaled identity. No tuning required.
#'
#' @param X  n x p numeric matrix (centred internally if centre = TRUE).
#' @param centre Logical; default TRUE.
#' @return   p x p positive-definite precision matrix.
#' @references Ledoit & Wolf (2004) J. Multivariate Anal. 88, 365-411.
#' @export
lw_precision <- function(X, centre = TRUE) {
  if (centre) X <- sweep(X, 2, colMeans(X))
  n <- nrow(X); p <- ncol(X)
  S  <- crossprod(X) / n
  mu <- sum(diag(S)) / p

  # Frobenius distance S vs mu*I
  delta2 <- sum((S - mu * diag(p))^2) / p

  # Oracle shrinkage intensity
  beta2 <- 0
  for (i in seq_len(n)) {
    Si    <- outer(X[i, ], X[i, ])
    beta2 <- beta2 + sum((Si - S)^2)
  }
  beta2 <- beta2 / (n^2 * p)
  alpha <- min(beta2 / delta2, 1)

  solve((1 - alpha) * S + alpha * mu * diag(p))
}


#' Graphical Lasso Precision Matrix Estimator
#'
#' CV-selected sparsity penalty; falls back to Ledoit-Wolf if glasso
#' is not installed.
#'
#' @param X        n x p numeric matrix.
#' @param rho_seq  Penalty grid. Default c(0.005,0.01,0.02,0.05,0.10,0.20).
#' @param val_frac Validation fraction for CV. Default 0.2.
#' @param centre   Logical; default TRUE.
#' @return p x p positive-definite precision matrix.
#' @references Friedman, Hastie & Tibshirani (2008) Biostatistics 9, 432-441.
#' @export
gl_precision <- function(X,
                         rho_seq  = c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20),
                         val_frac = 0.2,
                         centre   = TRUE) {
  if (!requireNamespace("glasso", quietly = TRUE)) {
    message("glasso not installed — falling back to Ledoit-Wolf.")
    return(lw_precision(X, centre = centre))
  }
  if (centre) X <- sweep(X, 2, colMeans(X))
  n <- nrow(X); p <- ncol(X)

  idx    <- sample(n)
  n_val  <- max(floor(n * val_frac), p + 2)
  n_tr   <- n - n_val
  S_tr   <- crossprod(X[idx[(n_val + 1):n], ]) / n_tr
  S_val  <- crossprod(X[idx[1:n_val],        ]) / n_val

  best_ll    <- -Inf
  best_Theta <- lw_precision(X, centre = FALSE)

  for (rho in rho_seq) {
    tryCatch({
      gl <- glasso::glasso(S_tr, rho = rho, penalize.diagonal = FALSE)
      Th <- gl$wi
      ld <- as.numeric(determinant(Th, logarithm = TRUE)$modulus)
      if (!is.finite(ld)) return(invisible(NULL))
      ll <- ld - sum(diag(S_val %*% Th))
      if (ll > best_ll) { best_ll <- ll; best_Theta <- Th }
    }, error = function(e) invisible(NULL))
  }
  best_Theta
}
