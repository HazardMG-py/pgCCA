# =============================================================================
# pgCCA — Unit Tests
# File:   tests/test_pgcca.R
#
# Run with: source("tests/test_pgcca.R")
# All tests print PASS or FAIL with a description.
# =============================================================================

source("R/precision.R")
source("R/whitening.R")
source("R/canonical_corr.R")
source("R/biclique.R")
source("R/pgcca.R")
source("R/simulation.R")
source("R/scca.R")
source("R/eps_sweep.R")

# ── Test helper ───────────────────────────────────────────────────────────────
n_pass <- 0L; n_fail <- 0L

expect <- function(desc, cond) {
  if (isTRUE(cond)) {
    cat(sprintf("  PASS  %s\n", desc)); n_pass <<- n_pass + 1L
  } else {
    cat(sprintf("  FAIL  %s\n", desc)); n_fail <<- n_fail + 1L
  }
}

cat("\n=== pgCCA Unit Tests ===\n\n")

# ─────────────────────────────────────────────────────────────────────────────
cat("-- precision.R --\n")

set.seed(1)
X10 <- matrix(rnorm(100 * 10), 100, 10)
Th  <- lw_precision(X10)
expect("lw_precision returns 10x10 matrix",   all(dim(Th) == c(10, 10)))
expect("lw_precision output is symmetric",     max(abs(Th - t(Th))) < 1e-10)
expect("lw_precision is positive definite",
       min(eigen(Th, only.values = TRUE)$values) > 0)

# Verify LW shrinks toward identity for iid data
X_iid <- matrix(rnorm(500 * 50), 500, 50)
Th_iid <- lw_precision(X_iid)
expect("lw_precision near identity for iid data",
       max(abs(Th_iid / mean(diag(Th_iid)) - diag(50))) < 0.5)

# ─────────────────────────────────────────────────────────────────────────────
cat("\n-- whitening.R --\n")

set.seed(2)
X20  <- matrix(rnorm(300 * 20), 300, 20)
Th20 <- lw_precision(X20)
Xw   <- whiten(X20, Th20)

expect("whiten returns same dimensions",  all(dim(Xw) == dim(X20)))
expect("whitened columns are centred",
       max(abs(colMeans(Xw))) < 1e-10)
expect("whitening reduces mean |corr|",
       mean(abs(cor(Xw)[lower.tri(cor(Xw))])) <
       mean(abs(cor(X20)[lower.tri(cor(X20))])))

# mat_sqrt_sym is symmetric
M    <- crossprod(matrix(rnorm(25), 5, 5)) + diag(5)
Msq  <- mat_sqrt_sym(M)
expect("mat_sqrt_sym is symmetric",   max(abs(Msq - t(Msq))) < 1e-10)
expect("mat_sqrt_sym squares to M",   max(abs(Msq %*% Msq - M)) < 1e-8)

# ─────────────────────────────────────────────────────────────────────────────
cat("\n-- canonical_corr.R --\n")

set.seed(3)
X5 <- matrix(rnorm(200 * 5), 200, 5)
Y5 <- matrix(rnorm(200 * 5), 200, 5)

rho_noise <- canonical_corr(X5, Y5, 1:3, 1:3)
expect("canonical_corr in [0,1] for noise",
       rho_noise >= 0 && rho_noise <= 1)
expect("canonical_corr near 0 for independent data", rho_noise < 0.4)

# Perfectly correlated subset
Yp <- cbind(X5[, 1:2], matrix(rnorm(200 * 3), 200, 3))
rho_high <- canonical_corr(X5, Yp, 1:2, 1:2)
expect("canonical_corr near 1 for identical variables", rho_high > 0.95)

# Empty index sets return 0
expect("canonical_corr returns 0 for empty ix",
       canonical_corr(X5, Y5, integer(0), 1:3) == 0)
expect("canonical_corr returns 0 for empty iy",
       canonical_corr(X5, Y5, 1:3, integer(0)) == 0)

# ─────────────────────────────────────────────────────────────────────────────
cat("\n-- simulation.R --\n")

set.seed(4)
dat <- generate_data(100, rho_within = 0.6, rho_cross = 0.3, seed = 4)

expect("generate_data X dimensions correct",   all(dim(dat$X) == c(100, 40)))
expect("generate_data Y dimensions correct",   all(dim(dat$Y) == c(100, 60)))
expect("generate_data ix_true correct",        identical(dat$ix_true, c(1L, 11L)))
expect("generate_data iy_true correct",        identical(dat$iy_true, c(1L, 16L, 31L)))
expect("generate_data X is scaled",
       max(abs(colMeans(dat$X))) < 0.1)

# eval_recovery: exact recovery
ev_exact <- eval_recovery(c(1L, 11L), c(1L, 16L, 31L),
                           c(1L, 11L), c(1L, 16L, 31L),
                           p = 40, q = 60)
expect("eval_recovery exact = 1 for perfect recovery",  ev_exact["exact"] == 1)
expect("eval_recovery sens = 1 for perfect recovery",   ev_exact["sens"] == 1)
expect("eval_recovery spec = 1 for perfect recovery",   ev_exact["spec"] == 1)
expect("eval_recovery near_exact = 1 for perfect",      ev_exact["near_exact"] == 1)
expect("eval_recovery jaccard = 1 for perfect",         ev_exact["jaccard"] == 1)

# eval_recovery: miss everything
ev_miss <- eval_recovery(integer(0), integer(0),
                          c(1L, 11L), c(1L, 16L, 31L),
                          p = 40, q = 60)
expect("eval_recovery sens = 0 for empty result",   ev_miss["sens"] == 0)
expect("eval_recovery spec = 1 for empty result",   ev_miss["spec"] == 1)
expect("eval_recovery exact = 0 for empty result",  ev_miss["exact"] == 0)

# ─────────────────────────────────────────────────────────────────────────────
cat("\n-- pgcca.R (main pipeline) --\n")

set.seed(5)
dat2 <- generate_data(300, rho_within = 0.6, rho_cross = 0.3, seed = 5)

res_lw <- pgcca(dat2$X, dat2$Y, method = "ledoitwolf", eps = 0.08)
expect("pgcca returns pgcca class",          inherits(res_lw, "pgcca"))
expect("pgcca ix is integer vector",         is.integer(res_lw$ix))
expect("pgcca iy is integer vector",         is.integer(res_lw$iy))
expect("pgcca rho in [0,1]",                 res_lw$rho >= 0 && res_lw$rho <= 1)
expect("pgcca $method stored correctly",     res_lw$method == "ledoitwolf")
expect("pgcca $valid field exists",          !is.null(res_lw$valid))

# gCCA (no whitening)
res_gc <- pgcca(dat2$X, dat2$Y, method = "none", eps = 0.08)
expect("gCCA (method=none) runs without error", inherits(res_gc, "pgcca"))
expect("gCCA rho in [0,1]",                     res_gc$rho >= 0 && res_gc$rho <= 1)

# pgCCA should recover the biclique at n=300 with rho_within=0.6
ev_lw <- eval_recovery(res_lw$ix, res_lw$iy,
                        dat2$ix_true, dat2$iy_true,
                        p = ncol(dat2$X), q = ncol(dat2$Y))
expect("pgCCA-LW achieves near-exact at n=300",
       ev_lw["near_exact"] == 1)

# Scenario C: weak within-set — all methods should work
dat_c <- generate_data(200, rho_within = 0.05, rho_cross = 0.3, seed = 77)
res_c  <- pgcca(dat_c$X, dat_c$Y, method = "ledoitwolf", eps = 0.08)
ev_c   <- eval_recovery(res_c$ix, res_c$iy,
                         dat_c$ix_true, dat_c$iy_true,
                         p = ncol(dat_c$X), q = ncol(dat_c$Y))
expect("pgCCA works when rho_within is small",
       ev_c["sens"] > 0.5)

# ─────────────────────────────────────────────────────────────────────────────
cat("\n-- Mechanism check (Lemma 1 values) --\n")

# Verify alpha_d and alpha_o analytically
rho_w <- 0.6; b <- 10
lam1  <- 1 + (b - 1) * rho_w   # 6.4
lam2  <- 1 - rho_w              # 0.4
alpha_d <- 1/sqrt(lam2) + (1/sqrt(lam1) - 1/sqrt(lam2)) / b
alpha_o <- (1/sqrt(lam1) - 1/sqrt(lam2)) / b
R_signal <- alpha_d * 0.3 * alpha_d
R_fp     <- abs(alpha_o) * 0.3 * alpha_d

expect("R_signal > eps=0.08",       R_signal > 0.08)
expect("R_fp < eps=0.08",           R_fp < 0.08)
expect("R_signal approx 0.642",     abs(R_signal - 0.642) < 0.01)
expect("R_fp approx 0.052",         abs(R_fp - 0.052) < 0.01)
expect("eta = min gap > 0",
       min((R_signal - 0.08)/2, (0.08 - R_fp)/2) > 0)

# ─────────────────────────────────────────────────────────────────────────────
cat(sprintf("\n=== Results: %d passed, %d failed ===\n\n", n_pass, n_fail))
if (n_fail > 0) {
  cat("Some tests failed. Check output above.\n")
} else {
  cat("All tests passed.\n")
}
