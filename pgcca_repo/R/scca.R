# =============================================================================
# pgCCA — sCCA Wrapper (PMA package)
# File:   R/scca.R
#
# Provides a clean interface to sparse CCA via the PMA package
# (Witten, Tibshirani & Hastie 2009) for use in simulations and
# method comparisons.
# =============================================================================

#' Sparse CCA via PMA Package
#'
#' Runs sparse CCA using permutation-based cross-validation to select
#' the L1 penalties c1 and c2. Returns identified variable sets and the
#' canonical correlation on the original data.
#'
#' @param X         n x p numeric matrix (original, unwhitened).
#' @param Y         n x q numeric matrix.
#' @param nperms    Number of permutations for CV. Default 5.
#' @param K         Number of canonical pairs. Default 1.
#' @param thr       Coefficient threshold for variable selection. Default 1e-6.
#'
#' @return List with:
#' \describe{
#'   \item{ix}{Integer vector of selected X-variable indices.}
#'   \item{iy}{Integer vector of selected Y-variable indices.}
#'   \item{rho}{Canonical correlation on the original data.}
#'   \item{n_x}{Number of selected X-variables.}
#'   \item{n_y}{Number of selected Y-variables.}
#' }
#'
#' @details
#' Falls back gracefully to an empty result if PMA is not installed
#' or if cross-validation fails. This wrapper is used for the simulation
#' comparison in the paper (Tables 2a/2b).
#'
#' @references
#' Witten, D.M., Tibshirani, R., and Hastie, T. (2009). A penalized matrix
#' decomposition, with applications to sparse principal components and
#' canonical correlation analysis. Biostatistics, 10(3), 515-534.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat <- generate_data(200)
#' res <- run_scca(dat$X, dat$Y)
#' cat("|IX|=", res$n_x, " |IY|=", res$n_y, " rho=", round(res$rho, 3))
#' }
#'
#' @export
run_scca <- function(X, Y, nperms = 5, K = 1, thr = 1e-6) {
  empty <- list(ix = integer(0), iy = integer(0), rho = 0,
                n_x = 0L, n_y = 0L)

  if (!requireNamespace("PMA", quietly = TRUE)) {
    message("PMA package not installed. Install with: install.packages('PMA')")
    return(empty)
  }

  tryCatch({
    cv <- PMA::CCA.permute(X, Y,
                           typex  = "standard",
                           typez  = "standard",
                           nperms = nperms,
                           trace  = FALSE)
    cr <- PMA::CCA(X, Y,
                   typex     = "standard",
                   typez     = "standard",
                   penaltyx  = cv$bestpenaltyx,
                   penaltyz  = cv$bestpenaltyz,
                   K         = K,
                   trace     = FALSE)

    ix <- which(abs(cr$u[, 1]) > thr)
    iy <- which(abs(cr$v[, 1]) > thr)

    if (length(ix) == 0 || length(iy) == 0) return(empty)

    rho <- canonical_corr(X, Y, ix, iy)
    list(ix = ix, iy = iy, rho = rho,
         n_x = length(ix), n_y = length(iy))

  }, error = function(e) {
    message("sCCA failed: ", conditionMessage(e))
    empty
  })
}
