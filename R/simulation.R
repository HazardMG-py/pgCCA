# =============================================================================
# pgCCA — Simulation Utilities
# File:   R/simulation.R
# =============================================================================

#' Generate Simulation Data with Block Covariance Structure
#'
#' Generates paired (X, Y) datasets with block-diagonal within-set covariance
#' and a biclique cross-covariance structure, matching the design used in
#' the paper (Design II, Scenarios A/B/C).
#'
#' @param n          Sample size.
#' @param rho_within Within-block correlation. Default 0.6.
#' @param rho_cross  Cross-modal biclique correlation. Default 0.3.
#' @param ix_true    Integer vector of true X-biclique indices.
#'   Default c(1, 11).
#' @param iy_true    Integer vector of true Y-biclique indices.
#'   Default c(1, 16, 31).
#' @param p          Total number of X-variables. Default 40.
#' @param q          Total number of Y-variables. Default 60.
#' @param n_blocks   Number of equal-size blocks. Default 4.
#' @param df         Degrees of freedom for t_df distribution. If NULL
#'   (default), Gaussian data are generated.
#' @param seed       Random seed. Default NULL.
#'
#' @return A list with:
#' \describe{
#'   \item{X}{n x p scaled matrix.}
#'   \item{Y}{n x q scaled matrix.}
#'   \item{ix_true}{True X-biclique indices.}
#'   \item{iy_true}{True Y-biclique indices.}
#'   \item{rho_canonical}{True first canonical correlation.}
#' }
#'
#' @examples
#' dat <- generate_data(200)
#' dim(dat$X)  # 200 x 40
#' dat$ix_true # c(1, 11)
#'
#' @export
generate_data <- function(n,
                          rho_within = 0.6,
                          rho_cross  = 0.3,
                          ix_true    = c(1L, 11L),
                          iy_true    = c(1L, 16L, 31L),
                          p          = 40L,
                          q          = 60L,
                          n_blocks   = 4L,
                          df         = NULL,
                          seed       = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Build block-diagonal within-set covariance
  make_block_cov <- function(d, nb, rw) {
    bsz <- d %/% nb
    S   <- matrix(0, d, d)
    for (b in seq_len(nb)) {
      s <- (b - 1) * bsz + 1; e <- b * bsz
      S[s:e, s:e] <- rw
    }
    diag(S) <- 1
    # Ensure positive definiteness
    lmin <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
    if (lmin <= 1e-6) S <- S + (1e-4 - lmin) * diag(d)
    S
  }

  Sx  <- make_block_cov(p, n_blocks, rho_within)
  Sy  <- make_block_cov(q, n_blocks, rho_within)
  Sxy <- matrix(0, p, q)
  for (i in ix_true) for (j in iy_true) Sxy[i, j] <- rho_cross

  Sigma <- rbind(cbind(Sx, Sxy), cbind(t(Sxy), Sy))
  lmin  <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  if (lmin <= 1e-6) Sigma <- Sigma + (1e-3 - lmin) * diag(p + q)

  # Generate data
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * (p + q)), n, p + q)
  if (!is.null(df)) {
    chi2 <- rchisq(n, df)
    Z    <- Z * sqrt(df / chi2)
  }
  D <- Z %*% L
  X <- scale(D[, seq_len(p)])
  Y <- scale(D[, (p + 1):(p + q)])

  # True canonical correlation
  rho_can <- rho_cross * sqrt(length(ix_true) * length(iy_true))

  list(X = X, Y = Y,
       ix_true = ix_true, iy_true = iy_true,
       rho_canonical = rho_can)
}


#' Evaluate Biclique Recovery Against Ground Truth
#'
#' Computes sensitivity, specificity, exact recovery, near-exact recovery,
#' and Jaccard index for an estimated biclique against the true one.
#'
#' @param ix_hat   Estimated X-variable indices.
#' @param iy_hat   Estimated Y-variable indices.
#' @param ix_true  True X-variable indices.
#' @param iy_true  True Y-variable indices.
#' @param p        Total number of X-variables.
#' @param q        Total number of Y-variables.
#'
#' @return Named numeric vector with elements:
#' \describe{
#'   \item{sens}{Sensitivity: fraction of true biclique variables recovered.}
#'   \item{spec}{Specificity: fraction of noise variables correctly excluded.}
#'   \item{exact}{1 if exact recovery, 0 otherwise.}
#'   \item{near_exact}{1 if sens >= 0.8 AND spec >= 0.8, 0 otherwise.}
#'   \item{jaccard}{Mean Jaccard index (X and Y averaged).}
#' }
#'
#' @examples
#' eval_recovery(ix_hat = 1:2, iy_hat = c(1, 16, 31),
#'               ix_true = c(1, 11), iy_true = c(1, 16, 31),
#'               p = 40, q = 60)
#'
#' @export
eval_recovery <- function(ix_hat, iy_hat,
                          ix_true, iy_true,
                          p, q) {
  # Sensitivity
  tp_x <- length(intersect(ix_hat, ix_true))
  tp_y <- length(intersect(iy_hat, iy_true))
  sens  <- (tp_x + tp_y) / (length(ix_true) + length(iy_true))

  # Specificity
  noise_x <- setdiff(seq_len(p), ix_true)
  noise_y <- setdiff(seq_len(q), iy_true)
  tn_x    <- length(intersect(setdiff(seq_len(p), ix_hat), noise_x))
  tn_y    <- length(intersect(setdiff(seq_len(q), iy_hat), noise_y))
  spec    <- (tn_x + tn_y) / (length(noise_x) + length(noise_y))

  # Exact recovery
  exact <- as.integer(setequal(ix_hat, ix_true) && setequal(iy_hat, iy_true))

  # Near-exact (sens >= 0.8 AND spec >= 0.8)
  near_exact <- as.integer(sens >= 0.8 && spec >= 0.8)

  # Jaccard
  jac_one <- function(h, t) {
    u <- length(union(h, t))
    if (u == 0) 1.0 else length(intersect(h, t)) / u
  }
  jaccard <- (jac_one(ix_hat, ix_true) + jac_one(iy_hat, iy_true)) / 2

  c(sens = sens, spec = spec, exact = exact,
    near_exact = near_exact, jaccard = jaccard)
}


#' Run a Full Simulation Scenario
#'
#' Runs a specified number of replications for one simulation scenario,
#' comparing multiple methods.
#'
#' @param methods    Named list where each element has fields:
#'   \describe{
#'     \item{method}{Character: "ledoitwolf", "glasso", or "none".}
#'     \item{eps}{Numeric threshold.}
#'   }
#' @param n_list     Integer vector of sample sizes. Default c(200, 500, 1000).
#' @param n_reps     Number of replications. Default 50.
#' @param rho_within Within-block correlation. Default 0.6.
#' @param rho_cross  Biclique correlation. Default 0.3.
#' @param df         Degrees of freedom (NULL for Gaussian). Default NULL.
#' @param seed       Base random seed. Default 42.
#' @param verbose    Logical. Default TRUE.
#' @param ...        Additional arguments passed to generate_data().
#'
#' @return Named list (one element per method), each containing a named
#'   list (one per n), each containing a matrix of n_reps x 6 results
#'   (sens, spec, exact, near_exact, jaccard, rho).
#'
#' @examples
#' \dontrun{
#' methods <- list(
#'   pgcca_lw = list(method = "ledoitwolf", eps = 0.08),
#'   gcca     = list(method = "none",       eps = 0.08)
#' )
#' store <- run_scenario(methods, n_list = c(200, 500), n_reps = 10)
#' }
#'
#' @export
run_scenario <- function(methods,
                         n_list     = c(200L, 500L, 1000L),
                         n_reps     = 50L,
                         rho_within = 0.6,
                         rho_cross  = 0.3,
                         df         = NULL,
                         seed       = 42L,
                         verbose    = TRUE,
                         ...) {
  store <- lapply(methods, function(m)
    setNames(lapply(n_list, function(n) matrix(NA_real_, n_reps, 6)),
             as.character(n_list))
  )

  for (n in n_list) {
    if (verbose) cat(sprintf("  n = %d ", n))
    for (rep in seq_len(n_reps)) {
      dat <- generate_data(n, rho_within = rho_within,
                           rho_cross = rho_cross, df = df,
                           seed = seed + rep * 1000L + n, ...)
      for (mn in names(methods)) {
        res <- tryCatch(
          pgcca(dat$X, dat$Y,
                method = methods[[mn]]$method,
                eps    = methods[[mn]]$eps,
                verbose = FALSE),
          error = function(e) list(ix = integer(0), iy = integer(0), rho = 0)
        )
        ev  <- eval_recovery(res$ix, res$iy,
                              dat$ix_true, dat$iy_true,
                              p = ncol(dat$X), q = ncol(dat$Y))
        store[[mn]][[as.character(n)]][rep, ] <-
          c(ev, rho = if (inherits(res, "pgcca")) res$rho else res$rho)
      }
    }
    if (verbose) cat("ok\n")
  }
  store
}


#' Summarise Simulation Results
#'
#' Prints a formatted table of mean performance metrics across sample sizes.
#'
#' @param store    Output from run_scenario().
#' @param n_list   Sample sizes used.
#' @param rho_true True canonical correlation (for MSE computation).
#' @param title    Title string printed above the table.
#'
#' @export
summarise_scenario <- function(store, n_list, rho_true, title = "") {
  metrics <- list(
    "Sensitivity"   = function(m) m[, 1],
    "Specificity"   = function(m) m[, 2],
    "Near-exact %"  = function(m) m[, 4] * 100,
    "Jaccard"       = function(m) m[, 5],
    "rho MSE x1e3"  = function(m) (m[, 6] - rho_true)^2 * 1e3
  )
  cat(sprintf("\n%s\n  %s\n%s\n", strrep("=", 60), title, strrep("=", 60)))
  for (lbl in names(metrics)) {
    cat(sprintf("\n  %s\n", lbl))
    hdr <- sprintf("  %-20s", "Method")
    for (n in n_list) hdr <- paste0(hdr, sprintf("  n=%5d", n))
    cat(hdr, "\n")
    for (mn in names(store)) {
      row <- sprintf("  %-20s", mn)
      for (n in n_list) {
        v   <- metrics[[lbl]](store[[mn]][[as.character(n)]])
        row <- paste0(row, sprintf("  %7.2f", mean(v, na.rm = TRUE)))
      }
      cat(row, "\n")
    }
  }
}
