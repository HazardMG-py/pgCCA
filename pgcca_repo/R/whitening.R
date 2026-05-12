# =============================================================================
# pgCCA — Stage 2: Whitening
# File:   R/whitening.R
# =============================================================================

#' Symmetric Matrix Square Root
#'
#' Computes M^{1/2} via eigendecomposition. Eigenvalues are floored at
#' a small positive value for numerical stability.
#'
#' @param M  Symmetric positive-definite matrix.
#' @param floor_val  Minimum eigenvalue floor. Default 1e-6.
#' @return   Matrix of same dimensions as M.
#' @export
mat_sqrt_sym <- function(M, floor_val = 1e-6) {
  e <- eigen(M, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, floor_val))) %*% t(e$vectors)
}


#' Whiten a Dataset Using an Estimated Precision Matrix
#'
#' Computes X_tilde = X * Theta^{1/2} and centres each column.
#' Does NOT rescale to unit variance — rescaling would undo the covariance
#' change induced by whitening.
#'
#' @param X      n x p numeric matrix.
#' @param Theta  p x p precision matrix (e.g. from lw_precision).
#' @return       n x p whitened and centred matrix.
#'
#' @details
#' The whitening transformation maps the covariance of X from Sigma to
#' approximately the identity, making previously correlated variables
#' independent. After whitening, the true cross-modal signal is amplified
#' (0.30 -> 0.64) while block-mate false positives are suppressed
#' (0.18 -> 0.05) for typical methylation parameters (rho_w = 0.6,
#' rho_c = 0.3, block size b = 10, threshold eps = 0.08).
#'
#' @examples
#' set.seed(1)
#' X     <- matrix(rnorm(200 * 40), 200, 40)
#' Theta <- lw_precision(X)
#' Xw    <- whiten(X, Theta)
#' # Within-set correlations are now close to zero
#' mean(abs(cor(Xw)[lower.tri(cor(Xw))]))
#'
#' @export
whiten <- function(X, Theta) {
  Th_sqrt <- mat_sqrt_sym(Theta)
  Xw      <- X %*% Th_sqrt
  sweep(Xw, 2, colMeans(Xw))   # centre; do NOT scale
}


#' Whitening Diagnostic
#'
#' Reports the mean absolute within-set pairwise correlation before and
#' after whitening, and the Ledoit-Wolf shrinkage intensity.
#' Use to confirm that whitening is having the expected effect.
#'
#' @param X        n x p numeric matrix (original).
#' @param Xw       n x p whitened matrix (from whiten()).
#' @param Theta    p x p precision matrix used for whitening.
#' @param n_sample Number of variables to sample for the correlation check.
#'   Default: min(100, p).
#' @param label    Character label printed in the output. Default "X".
#' @return Invisibly returns a list with raw_mean_abscor and white_mean_abscor.
#' @export
whitening_diagnostic <- function(X, Xw, Theta,
                                 n_sample = NULL, label = "X") {
  p <- ncol(X)
  n_sample <- min(n_sample %||% 100L, p)
  idx <- sample(p, n_sample)

  R_raw   <- cor(X[,  idx])
  R_white <- cor(Xw[, idx])
  raw_m   <- mean(abs(R_raw[lower.tri(R_raw)]))
  white_m <- mean(abs(R_white[lower.tri(R_white)]))
  pct_red <- 100 * (1 - white_m / raw_m)

  # Ledoit-Wolf shrinkage intensity
  n <- nrow(X); S <- crossprod(sweep(X, 2, colMeans(X))) / n
  mu <- sum(diag(S)) / p
  delta2 <- sum((S - mu * diag(p))^2) / p
  Xc <- sweep(X, 2, colMeans(X)); beta2 <- 0
  for (i in seq_len(n)) { Si <- outer(Xc[i,], Xc[i,]); beta2 <- beta2 + sum((Si-S)^2) }
  alpha <- min((beta2 / (n^2 * p)) / delta2, 1)

  cat(sprintf("--- Whitening diagnostic: %s ---\n", label))
  cat(sprintf("  Mean |corr| raw:      %.4f\n", raw_m))
  cat(sprintf("  Mean |corr| whitened: %.4f  (%.1f%% reduction)\n",
              white_m, pct_red))
  cat(sprintf("  LW shrinkage alpha:   %.4f  (0=no shrink, 1=full shrink)\n",
              alpha))
  if (pct_red > 20)
    cat("  -> Meaningful block structure removed by whitening\n")
  else
    cat("  -> Whitening effect small (weak within-set structure)\n")

  invisible(list(raw_mean_abscor   = raw_m,
                 white_mean_abscor = white_m,
                 reduction_pct     = pct_red,
                 lw_alpha          = alpha))
}

# Internal null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b
