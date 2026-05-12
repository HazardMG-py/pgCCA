# =============================================================================
# pgCCA — Real Data Vignette: TCGA-GBM
# File:   vignettes/tcga_gbm.R
#
# Reproduces Tables 3 and 4 (eps sweep) from the paper.
# Data from LinkedOmics: http://linkedomics.org  (dataset: TCGA-GBM)
#
# Download these two files before running:
#   Human__TCGA_GBM__JHU_USC__Methylation__Meth27K__...__cct.gz
#   Human__TCGA_GBM__UNC__RNAseq__GA_RNA__...__cct.gz
# Then set DATA_DIR below.
# =============================================================================

source("R/precision.R")
source("R/whitening.R")
source("R/canonical_corr.R")
source("R/biclique.R")
source("R/pgcca.R")
source("R/eps_sweep.R")

# ── Settings ──────────────────────────────────────────────────────────────────
DATA_DIR  <- "data/"          # path to downloaded LinkedOmics files
TOP_P     <- 6000L            # top methylation sites by variance
TOP_Q     <- 8000L            # top genes by variance
NA_THRESH <- 0.80             # max proportion of NAs per feature
EPS       <- 0.08

# ── Data loading helper ───────────────────────────────────────────────────────
load_linkedomics <- function(path, sep = "\t") {
  df  <- read.table(path, header = TRUE, sep = sep,
                    row.names = 1, check.names = FALSE,
                    comment.char = "", stringsAsFactors = FALSE)
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  t(mat)    # samples x features
}

# ── Load methylation and expression ──────────────────────────────────────────
cat("Loading methylation data...\n")
meth_files <- list.files(DATA_DIR, pattern = "Methylation.*cct", full.names = TRUE)
expr_files <- list.files(DATA_DIR, pattern = "RNAseq.*cct",      full.names = TRUE)

if (length(meth_files) == 0 || length(expr_files) == 0) {
  stop("Data files not found in ", DATA_DIR,
       "\nDownload from http://linkedomics.org (TCGA-GBM dataset)\n",
       "and set DATA_DIR to the correct path.")
}

meth_raw <- load_linkedomics(meth_files[1])
expr_raw <- load_linkedomics(expr_files[1])
cat(sprintf("  Methylation: %d samples x %d sites\n",  nrow(meth_raw), ncol(meth_raw)))
cat(sprintf("  Expression:  %d samples x %d genes\n",  nrow(expr_raw), ncol(expr_raw)))

# ── Match subjects ────────────────────────────────────────────────────────────
common <- intersect(rownames(meth_raw), rownames(expr_raw))
cat(sprintf("  Matched samples: n = %d\n", length(common)))
meth <- meth_raw[common, ]; expr <- expr_raw[common, ]

# ── Preprocessing ─────────────────────────────────────────────────────────────
# Remove features with too many NAs
na_frac_m <- colMeans(is.na(meth)); meth <- meth[, na_frac_m <= NA_THRESH]
na_frac_e <- colMeans(is.na(expr)); expr <- expr[, na_frac_e <= NA_THRESH]
cat(sprintf("  After NA filter (>%.0f%%): meth=%d sites, expr=%d genes\n",
            NA_THRESH * 100, ncol(meth), ncol(expr)))

# Impute remaining NAs with column medians
impute_median <- function(M) {
  for (j in seq_len(ncol(M))) {
    nas <- is.na(M[, j])
    if (any(nas)) M[nas, j] <- median(M[, j], na.rm = TRUE)
  }
  M
}
meth <- impute_median(meth); expr <- impute_median(expr)

# Select top-variance features
var_m <- apply(meth, 2, var); top_m <- order(var_m, decreasing = TRUE)[seq_len(TOP_P)]
var_e <- apply(expr, 2, var); top_e <- order(var_e, decreasing = TRUE)[seq_len(TOP_Q)]
X_sub <- scale(meth[, top_m])   # n x TOP_P, zero mean unit variance
Y_sub <- scale(expr[, top_e])   # n x TOP_Q
n     <- nrow(X_sub)
cat(sprintf("  Final: n=%d, p=%d, q=%d\n", n, ncol(X_sub), ncol(Y_sub)))

# ── Table 3: pgCCA vs gCCA at eps=0.08 ───────────────────────────────────────
cat("\n--- Table 3: Main comparison at eps =", EPS, "---\n")

cat("Running pgCCA-LW...\n")
res_pg <- pgcca(X_sub, Y_sub, method = "ledoitwolf", eps = EPS, verbose = TRUE)

cat("\nRunning gCCA (no whitening)...\n")
res_gc <- pgcca(X_sub, Y_sub, method = "none", eps = EPS, verbose = FALSE)

cat(sprintf("\n%-25s  %6s  %6s  %8s  %7s  %s\n",
            "Method", "|IX|", "|IY|", "Total", "rho", "Valid?"))
cat(strrep("-", 65), "\n")
cat(sprintf("%-25s  %6d  %6d  %8d  %7.4f  %s\n",
            "pgCCA-LW",
            res_pg$n_x, res_pg$n_y, res_pg$n_x + res_pg$n_y,
            res_pg$rho,
            if (res_pg$valid) "YES" else "NO (OVERFIT)"))
cat(sprintf("%-25s  %6d  %6d  %8d  %7.4f  %s\n",
            "gCCA (no whitening)",
            res_gc$n_x, res_gc$n_y, res_gc$n_x + res_gc$n_y,
            res_gc$rho,
            if (res_gc$valid) "YES" else "NO (OVERFIT)"))

# ── Table 4: Epsilon sweep ────────────────────────────────────────────────────
cat("\n--- Table 4: Epsilon sweep ---\n")

eps_seq <- c(0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25, 0.28, 0.30)

cat("\ngCCA sweep:\n")
sweep_gc <- eps_sweep(X_sub, Y_sub, eps_seq = eps_seq, method = "none")

cat("\npgCCA-LW sweep:\n")
sweep_pg <- eps_sweep(X_sub, Y_sub, eps_seq = eps_seq, method = "ledoitwolf")

print_eps_comparison(sweep_gc, sweep_pg, n)

# ── Save results ──────────────────────────────────────────────────────────────
results <- list(
  n          = n,
  eps        = EPS,
  pgcca_lw   = res_pg,
  gcca       = res_gc,
  sweep_gcca = sweep_gc,
  sweep_pgcca= sweep_pg
)
saveRDS(results, file.path(DATA_DIR, "tcga_gbm_pgcca_results.rds"))
cat("\nResults saved to", file.path(DATA_DIR, "tcga_gbm_pgcca_results.rds"), "\n")

cat("\n=== pgCCA identified genes ===\n")
if (!is.null(colnames(Y_sub)) && length(res_pg$iy) > 0) {
  genes <- colnames(Y_sub)[res_pg$iy]
  cat(sprintf("  %d genes identified:\n  ", length(genes)))
  cat(paste(head(genes, 20), collapse = ", "))
  if (length(genes) > 20) cat(sprintf(" ... (%d total)", length(genes)))
  cat("\n")
  write.csv(data.frame(gene = genes), "data/pgcca_identified_genes.csv",
            row.names = FALSE)
  cat("  Saved to data/pgcca_identified_genes.csv\n")
  cat("  Submit to MSigDB Hallmark enrichment at:\n")
  cat("  https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp\n")
}
