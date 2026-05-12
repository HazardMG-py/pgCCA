# =============================================================================
# pgCCA — Main Pipeline
# File:   R/pgcca.R
# =============================================================================

#' Precision-Guided Graph Canonical Correlation Analysis
#'
#' Three-stage pipeline for recovering cross-modal biclique associations
#' between two high-dimensional datasets:
#' \enumerate{
#'   \item Estimate precision matrices (Ledoit-Wolf or graphical Lasso).
#'   \item Whiten each dataset using the precision square root.
#'   \item Run gCCA greedy biclique detection on the whitened data.
#' }
#' The canonical correlation is reported on the ORIGINAL data using the
#' identified variable sets.
#'
#' @param X       n x p numeric matrix (n observations, p variables).
#'   Should be mean-centred and scaled to unit variance before calling.
#' @param Y       n x q numeric matrix. Same subjects as X.
#' @param method  Precision estimator: "ledoitwolf" (default), "glasso",
#'   or "none" (runs standard gCCA without whitening).
#' @param eps     Detection threshold for the cross-correlation biclique
#'   graph. Default 0.08. Recommended range: 0.05-0.15.
#' @param lambda_seq Lambda grid for the KL-divergence selector.
#'   Default seq(0.5, 0.85, by = 0.05).
#' @param rho_seq    Penalty grid for graphical Lasso CV (method = "glasso").
#'   Default c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20).
#' @param verbose Logical. If TRUE, prints progress. Default FALSE.
#'
#' @return A list with class "pgcca" containing:
#' \describe{
#'   \item{ix}{Integer vector of identified X-variable indices.}
#'   \item{iy}{Integer vector of identified Y-variable indices.}
#'   \item{rho}{Canonical correlation on the original data (in [0,1]).}
#'   \item{lambda_opt}{Optimal lambda chosen by KL divergence.}
#'   \item{method}{Precision method used.}
#'   \item{eps}{Detection threshold used.}
#'   \item{n_x}{|IX| = number of identified X-variables.}
#'   \item{n_y}{|IY| = number of identified Y-variables.}
#'   \item{valid}{Logical: TRUE if n_x + n_y < n (not overfitting).}
#' }
#'
#' @details
#' \strong{Why whitening helps.}
#' Within-set block correlation (common in DNA methylation CpG islands,
#' gene co-expression modules, financial sector clusters) creates
#' false-positive associations in the raw cross-correlation matrix.
#' If site A truly regulates gene X with cross-correlation rho_c = 0.3,
#' and site B sits in the same CpG island (within-block correlation
#' rho_w = 0.6), then site B has spurious cross-correlation
#' rho_w * rho_c = 0.18 > eps = 0.08. gCCA includes site B as a false
#' positive. After whitening with the precision square root, the true
#' signal becomes alpha_d * rho_c * alpha_d ~= 0.642 >> eps, and the
#' false positive becomes alpha_o * rho_c * alpha_d ~= 0.052 < eps.
#' Only the true variable survives the threshold.
#'
#' \strong{Validity check.}
#' When n_x + n_y >= n, the canonical correlation is 1.0 by rank
#' deficiency and the result is not meaningful. The $valid field flags
#' this. If method = "none" (gCCA) and valid = FALSE, try method =
#' "ledoitwolf".
#'
#' @references
#' Park, H., Bai, S., Ye, Z., Lee, H., Ma, T., and Chen, S. (2025).
#' Graph canonical correlation analysis. arXiv:2502.01780.
#'
#' Ledoit, O. and Wolf, M. (2004). A well-conditioned estimator for
#' large-dimensional covariance matrices. J. Multivariate Anal. 88, 365-411.
#'
#' @examples
#' # Simulate a small biclique
#' set.seed(42)
#' n <- 200; p <- 40; q <- 60
#' IX_true <- c(1, 11); IY_true <- c(1, 16, 31)
#'
#' # Build covariance with block structure and biclique signal
#' Sigma_x <- diag(p); for (b in 0:3) Sigma_x[(b*10+1):(b*10+10), (b*10+1):(b*10+10)] <- 0.6
#' diag(Sigma_x) <- 1
#' Sigma_y <- diag(q); for (b in 0:5) Sigma_y[(b*10+1):(b*10+10), (b*10+1):(b*10+10)] <- 0.6
#' diag(Sigma_y) <- 1
#' Sigma_xy <- matrix(0, p, q); Sigma_xy[IX_true, IY_true] <- 0.3
#' Sigma <- rbind(cbind(Sigma_x, Sigma_xy), cbind(t(Sigma_xy), Sigma_y))
#'
#' dat <- MASS::mvrnorm(n, mu = rep(0, p + q), Sigma = Sigma)
#' X   <- scale(dat[, 1:p])
#' Y   <- scale(dat[, (p+1):(p+q)])
#'
#' # Run pgCCA (with Ledoit-Wolf whitening)
#' res_pg <- pgcca(X, Y, method = "ledoitwolf", eps = 0.08)
#' cat("pgCCA: |IX|=", res_pg$n_x, " |IY|=", res_pg$n_y,
#'     " rho=", round(res_pg$rho, 4), "\n")
#'
#' # Run standard gCCA (no whitening)
#' res_gc <- pgcca(X, Y, method = "none", eps = 0.08)
#' cat("gCCA:  |IX|=", res_gc$n_x, " |IY|=", res_gc$n_y,
#'     " rho=", round(res_gc$rho, 4), "\n")
#'
#' @export
pgcca <- function(X, Y,
                  method     = c("ledoitwolf", "glasso", "none"),
                  eps        = 0.08,
                  lambda_seq = seq(0.5, 0.85, by = 0.05),
                  rho_seq    = c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20),
                  verbose    = FALSE) {

  method <- match.arg(method)
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  stopifnot(nrow(Y) == n, is.numeric(X), is.numeric(Y))

  # ── Stage 1: Precision estimation ─────────────────────────────────────────
  if (verbose) cat("Stage 1: precision estimation (method =", method, ")...\n")

  if (method == "ledoitwolf") {
    Tx <- lw_precision(X); Ty <- lw_precision(Y)
  } else if (method == "glasso") {
    Tx <- gl_precision(X, rho_seq = rho_seq)
    Ty <- gl_precision(Y, rho_seq = rho_seq)
  } else {
    Tx <- diag(p); Ty <- diag(q)
  }

  # ── Stage 2: Whitening ────────────────────────────────────────────────────
  if (verbose) cat("Stage 2: whitening...\n")
  Xw <- whiten(X, Tx)
  Yw <- whiten(Y, Ty)

  # ── Stage 3: Biclique detection ───────────────────────────────────────────
  if (verbose) cat("Stage 3: greedy biclique detection...\n")
  res <- run_gcca(Xw, Yw, eps = eps, lambda_seq = lambda_seq)

  # Canonical correlation on ORIGINAL data
  rho <- canonical_corr(X, Y, res$ix, res$iy)

  out <- list(
    ix         = res$ix,
    iy         = res$iy,
    rho        = rho,
    lambda_opt = res$lambda_opt,
    method     = method,
    eps        = eps,
    n_x        = length(res$ix),
    n_y        = length(res$iy),
    valid      = (length(res$ix) + length(res$iy)) < n
  )
  class(out) <- "pgcca"
  out
}


#' Print Method for pgcca Objects
#' @param x   A pgcca object.
#' @param ... Further arguments (ignored).
#' @export
print.pgcca <- function(x, ...) {
  cat("pgCCA result\n")
  cat("  Method:      ", x$method, "\n")
  cat("  Threshold:   ", x$eps,    "\n")
  cat("  |IX|:        ", x$n_x,    "\n")
  cat("  |IY|:        ", x$n_y,    "\n")
  cat("  rho:         ", round(x$rho, 4), "\n")
  cat("  lambda_opt:  ", round(x$lambda_opt, 3), "\n")
  cat("  Valid (n_x+n_y < n):", x$valid, "\n")
  if (!x$valid)
    cat("  WARNING: total vars >= n. rho=1.0 may be trivial.\n",
        "  Consider using method='ledoitwolf' or increasing eps.\n")
  invisible(x)
}


#' Summary Method for pgcca Objects
#' @param object A pgcca object.
#' @param X      Original X matrix (for named variable output). Optional.
#' @param Y      Original Y matrix. Optional.
#' @param ...    Further arguments (ignored).
#' @export
summary.pgcca <- function(object, X = NULL, Y = NULL, ...) {
  print(object)
  if (!is.null(X) && !is.null(colnames(X)) && length(object$ix) > 0) {
    cat("\n  X variables identified:\n  ")
    cat(paste(head(colnames(X)[object$ix], 10), collapse = ", "))
    if (length(object$ix) > 10)
      cat(sprintf(" ... (%d total)", length(object$ix)))
    cat("\n")
  }
  if (!is.null(Y) && !is.null(colnames(Y)) && length(object$iy) > 0) {
    cat("  Y variables identified:\n  ")
    cat(paste(head(colnames(Y)[object$iy], 10), collapse = ", "))
    if (length(object$iy) > 10)
      cat(sprintf(" ... (%d total)", length(object$iy)))
    cat("\n")
  }
  invisible(object)
}
