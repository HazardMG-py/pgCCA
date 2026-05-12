# =============================================================================
# pgCCA — Epsilon Sweep Diagnostic
# File:   R/eps_sweep.R
#
# Runs gCCA and pgCCA across a range of detection thresholds.
# Use to diagnose whether gCCA overfits (total vars >= n) on real data
# and to confirm pgCCA produces valid results across all thresholds.
# =============================================================================

#' Epsilon Sweep: Run pgCCA/gCCA Across Detection Thresholds
#'
#' For each value of eps in eps_seq, runs the specified method and records
#' the number of selected X-variables, Y-variables, total variables, the
#' canonical correlation, and whether the result is valid (total < n).
#'
#' @param X       n x p numeric matrix.
#' @param Y       n x q numeric matrix.
#' @param eps_seq Numeric vector of threshold values to sweep.
#'   Default seq(0.08, 0.30, by = 0.02).
#' @param method  Precision method: "ledoitwolf", "glasso", or "none".
#'   Default "ledoitwolf".
#' @param verbose Logical. If TRUE, prints a progress line per eps value.
#'   Default TRUE.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{eps}{Detection threshold.}
#'   \item{n_x}{Number of selected X-variables.}
#'   \item{n_y}{Number of selected Y-variables.}
#'   \item{total}{n_x + n_y.}
#'   \item{rho}{Estimated canonical correlation.}
#'   \item{valid}{Logical: total < n.}
#'   \item{method}{Method used.}
#' }
#'
#' @details
#' The primary use case is to demonstrate that standard gCCA (method =
#' "none") overfits at every threshold on real methylation data, while
#' pgCCA (method = "ledoitwolf") produces valid results. See Table 4
#' (eps sweep) in the paper.
#'
#' @examples
#' \dontrun{
#' # Compare gCCA vs pgCCA across eps on TCGA-GBM data
#' sweep_gc <- eps_sweep(X_sub, Y_sub, method = "none")
#' sweep_pg <- eps_sweep(X_sub, Y_sub, method = "ledoitwolf")
#' print(sweep_gc)
#' print(sweep_pg)
#' }
#'
#' @export
eps_sweep <- function(X, Y,
                      eps_seq = seq(0.08, 0.30, by = 0.02),
                      method  = "ledoitwolf",
                      verbose = TRUE) {
  n   <- nrow(X)
  out <- vector("list", length(eps_seq))

  for (k in seq_along(eps_seq)) {
    eps_k <- eps_seq[k]
    res   <- tryCatch(
      pgcca(X, Y, method = method, eps = eps_k, verbose = FALSE),
      error = function(e)
        list(n_x = NA, n_y = NA, rho = NA, valid = NA)
    )
    total <- res$n_x + res$n_y
    out[[k]] <- data.frame(
      eps    = eps_k,
      n_x    = res$n_x,
      n_y    = res$n_y,
      total  = total,
      rho    = round(res$rho, 4),
      valid  = isTRUE(total < n),
      method = method,
      stringsAsFactors = FALSE
    )
    if (verbose) {
      status <- if (isTRUE(total < n)) "OK " else "OVERFIT"
      cat(sprintf("  eps=%.2f  |IX|=%5d  |IY|=%5d  total=%6d  rho=%.4f  %s\n",
                  eps_k, res$n_x, res$n_y, total, res$rho, status))
    }
  }
  do.call(rbind, out)
}


#' Print Eps Sweep Comparison Table
#'
#' Prints a side-by-side comparison of gCCA and pgCCA eps sweep results.
#'
#' @param sweep_gcca  Data.frame from eps_sweep(method = "none").
#' @param sweep_pgcca Data.frame from eps_sweep(method = "ledoitwolf").
#' @param n           Sample size (for validity check label).
#'
#' @export
print_eps_comparison <- function(sweep_gcca, sweep_pgcca, n) {
  cat(sprintf("\n%s\n", strrep("=", 72)))
  cat(sprintf("  eps sweep comparison  (n = %d)\n", n))
  cat(sprintf("%s\n", strrep("=", 72)))
  cat(sprintf("  %-6s  %-22s  %-22s\n",
              "eps",
              "gCCA (no whitening)",
              "pgCCA-LW"))
  cat(sprintf("  %-6s  %-10s %-10s  %-10s %-10s\n",
              "", "total", "rho", "total", "rho"))
  cat(sprintf("%s\n", strrep("-", 72)))

  for (i in seq_len(nrow(sweep_gcca))) {
    gc_tot  <- sweep_gcca$total[i]
    pg_tot  <- sweep_pgcca$total[i]
    gc_flag <- if (gc_tot >= n) " OVERFIT" else " ok"
    pg_flag <- if (pg_tot >= n) " OVERFIT" else " ok"
    cat(sprintf("  %.2f    %6d%s     %.4f    %6d%s     %.4f\n",
                sweep_gcca$eps[i],
                gc_tot, gc_flag, sweep_gcca$rho[i],
                pg_tot, pg_flag, sweep_pgcca$rho[i]))
  }
  cat(sprintf("%s\n", strrep("=", 72)))
}
