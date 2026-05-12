# =============================================================================
# pgCCA — Reproduce Simulation Results (Tables 2a and 2b)
# File:   vignettes/reproduce_simulation.R
#
# Reproduces the Design II simulation from the paper.
# Runtime: ~60 min for n_reps = 50. Set n_reps = 10 for a quick test.
# =============================================================================

# ── Load all pgCCA modules ────────────────────────────────────────────────────
source("R/precision.R")
source("R/whitening.R")
source("R/canonical_corr.R")
source("R/biclique.R")
source("R/pgcca.R")
source("R/simulation.R")

# ── Settings ──────────────────────────────────────────────────────────────────
N_REPS  <- 50    # set to 10 for a quick test
N_LIST  <- c(200, 500, 1000)
EPS     <- 0.08

METHODS <- list(
  pgcca_lw = list(method = "ledoitwolf", eps = EPS),
  pgcca_gl = list(method = "glasso",     eps = EPS),
  gcca     = list(method = "none",       eps = EPS)
)

IX_TRUE <- c(1L, 11L)
IY_TRUE <- c(1L, 16L, 31L)
RHO_CAN <- 0.3 * sqrt(length(IX_TRUE) * length(IY_TRUE))
cat(sprintf("True canonical correlation: %.4f\n", RHO_CAN))

# ── Scenario A: Gaussian, rho_within = 0.6 ───────────────────────────────────
cat("\nRunning Scenario A (Gaussian, rho_within=0.6)...\n")
store_A <- run_scenario(
  methods    = METHODS,
  n_list     = N_LIST,
  n_reps     = N_REPS,
  rho_within = 0.6,
  rho_cross  = 0.3,
  df         = NULL,
  seed       = 42L,
  verbose    = TRUE,
  ix_true    = IX_TRUE,
  iy_true    = IY_TRUE
)
summarise_scenario(store_A, N_LIST, RHO_CAN,
                   "Scenario A — Gaussian rho_within=0.6")

# ── Scenario B: t5, rho_within = 0.6 ─────────────────────────────────────────
cat("\nRunning Scenario B (t5, rho_within=0.6)...\n")
store_B <- run_scenario(
  methods    = METHODS,
  n_list     = N_LIST,
  n_reps     = N_REPS,
  rho_within = 0.6,
  rho_cross  = 0.3,
  df         = 5,
  seed       = 99L,
  verbose    = TRUE,
  ix_true    = IX_TRUE,
  iy_true    = IY_TRUE
)
summarise_scenario(store_B, N_LIST, RHO_CAN,
                   "Scenario B — t5  rho_within=0.6")

# ── Scenario C: Gaussian, rho_within = 0.05 ──────────────────────────────────
cat("\nRunning Scenario C (Gaussian, rho_within=0.05)...\n")
store_C <- run_scenario(
  methods    = METHODS,
  n_list     = N_LIST,
  n_reps     = N_REPS,
  rho_within = 0.05,
  rho_cross  = 0.3,
  df         = NULL,
  seed       = 77L,
  verbose    = TRUE,
  ix_true    = IX_TRUE,
  iy_true    = IY_TRUE
)
summarise_scenario(store_C, N_LIST, RHO_CAN,
                   "Scenario C — Gaussian rho_within=0.05")

cat("\nDone. store_A, store_B, store_C contain full results.\n")
cat("Use summarise_scenario() to print any scenario.\n")
