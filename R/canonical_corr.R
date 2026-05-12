# =============================================================================
# pgCCA — Canonical Correlation
# File:   R/canonical_corr.R
# =============================================================================

#' Canonical Correlation for an Identified Biclique
#'
#' Computes the first sample canonical correlation between the subsets
#' X[, ix] and Y[, iy] using the formula:
#'   rho = sigma_max( C_xx^{-1/2} C_xy C_yy^{-1/2} )
#' This always returns a value in [0, 1], unlike raw SVD of the
#' cross-covariance which can exceed 1.
#'
#' @param X   n x p matrix (original, unwhitened).
#' @param Y   n x q matrix (original, unwhitened).
#' @param ix  Integer vector of X-variable indices.
#' @param iy  Integer vector of Y-variable indices.
#' @param eps Small regularisation added to diagonal for numerical stability.
#'   Default 1e-6.
#'
#' @return Scalar in [0, 1].
#'
#' @details
#' Always called on the ORIGINAL (unwhitened) X and Y after the variable
#' sets have been identified by pgCCA or gCCA. Whitening is used only for
#' variable selection; the canonical correlation is reported in the original
#' data space.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200 * 10), 200, 10)
#' Y <- matrix(rnorm(200 *  8), 200,  8)
#' rho <- canonical_corr(X, Y, ix = 1:3, iy = 1:3)
#' rho  # should be close to 0 (independent random data)
#'
#' @export
canonical_corr <- function(X, Y, ix, iy, eps = 1e-6) {
  if (length(ix) == 0 || length(iy) == 0) return(0)

  Xs <- X[, ix, drop = FALSE]
  Ys <- Y[, iy, drop = FALSE]
  n  <- nrow(X)
  p0 <- length(ix); q0 <- length(iy)

  Cxx <- crossprod(Xs) / n + eps * diag(p0)
  Cyy <- crossprod(Ys) / n + eps * diag(q0)
  Cxy <- crossprod(Xs, Ys) / n

  isq <- function(M) {
    e <- eigen(M, symmetric = TRUE)
    d <- 1 / sqrt(pmax(e$values, eps))
    if (length(d) == 1)
      matrix(d * e$vectors^2, 1, 1)
    else
      e$vectors %*% diag(d) %*% t(e$vectors)
  }

  sv <- tryCatch(
    svd(isq(Cxx) %*% Cxy %*% isq(Cyy))$d[1],
    error = function(e) 0
  )
  min(max(sv, 0), 1)
}
