# pgCCA: Precision-Guided Graph Canonical Correlation Analysis

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)

R implementation of pgCCA, a three-stage method for recovering
cross-modal biclique associations between two high-dimensional datasets.
Accompanies the paper:

> **Precision-Guided Graph Canonical Correlation Analysis**
> *Statistics and Innovation, 2025 (under review)*

---

## What pgCCA does

When you have two datasets measured on the same subjects (e.g. DNA
methylation and gene expression from cancer patients), pgCCA finds which
variables in dataset X are genuinely associated with which variables in
dataset Y — as a biclique subgraph.

The key improvement over standard gCCA (Park et al. 2025): pgCCA first
*whitens* each dataset using its estimated precision matrix. This removes
within-set block correlation (e.g. co-methylation within CpG islands) that
otherwise floods the cross-correlation matrix with false positives.

| | Raw corr | Whitened corr |
|--|--|--|
| True signal | 0.30 | **0.64** |
| Block-mate (false positive) | 0.18 | **0.05** |
| Detection threshold ε | 0.08 | 0.08 |

---

## Installation

```r
# Install dependencies
install.packages(c("glasso"))

# Install pgCCA from GitHub
# install.packages("remotes")
remotes::install_github("YOUR_USERNAME/pgCCA")
```

Or clone and source locally:

```r
source("R/precision.R")
source("R/whitening.R")
source("R/canonical_corr.R")
source("R/biclique.R")
source("R/pgcca.R")
source("R/simulation.R")
```

---

## Quick start

```r
library(pgCCA)

# Simulate data with block structure and biclique signal
set.seed(42)
dat <- generate_data(n = 200, rho_within = 0.6, rho_cross = 0.3)
X   <- dat$X   # 200 x 40
Y   <- dat$Y   # 200 x 60

# Run pgCCA (Ledoit-Wolf whitening)
res_pg <- pgcca(X, Y, method = "ledoitwolf", eps = 0.08)
print(res_pg)
# pgCCA result
#   Method:      ledoitwolf
#   |IX|:        2
#   |IY|:        3
#   rho:         0.7312
#   Valid:       TRUE

# Run standard gCCA (no whitening)
res_gc <- pgcca(X, Y, method = "none", eps = 0.08)
print(res_gc)
```

---

## Real data (TCGA-GBM)

Download the data from LinkedOmics (http://linkedomics.org,
dataset: TCGA-GBM), then:

```r
# Load and align (see vignettes/tcga_gbm.R for full preprocessing)
X_sub <- scale(meth_data)   # n x 6000 methylation
Y_sub <- scale(expr_data)   # n x 8000 expression

# Run pgCCA
res <- pgcca(X_sub, Y_sub, method = "ledoitwolf", eps = 0.08, verbose = TRUE)
print(res)
# |IX| = 128, |IY| = 137, rho = 0.9961, valid = TRUE

# Compare: standard gCCA always overfits on this data
res_gc <- pgcca(X_sub, Y_sub, method = "none", eps = 0.08)
# |IX| = 4003, |IY| = 4664, rho = 1.0000, valid = FALSE (4667 > n=278)
```

---

## Reproduce paper results

```r
# Design II simulation (Scenarios A, B, C)
source("vignettes/reproduce_simulation.R")

# Real data analysis
source("vignettes/tcga_gbm.R")
```

---

## File structure

```
pgCCA/
├── R/
│   ├── precision.R       # Stage 1: Ledoit-Wolf and GLasso estimators
│   ├── whitening.R       # Stage 2: Whitening transformation + diagnostics
│   ├── biclique.R        # Stage 3: Greedy biclique + KL-divergence selector
│   ├── canonical_corr.R  # Canonical correlation on identified variable sets
│   ├── pgcca.R           # Main pgcca() function + print/summary methods
│   └── simulation.R      # Data generation, evaluation metrics, scenario runner
├── tests/
│   └── test_pgcca.R      # Unit tests
├── vignettes/
│   ├── reproduce_simulation.R   # Reproduce Tables 2a/2b from paper
│   └── tcga_gbm.R               # Reproduce Table 3 (real data)
└── README.md
```

---

## Key functions

| Function | Purpose |
|----------|---------|
| `pgcca(X, Y, method, eps)` | Main pipeline: estimate, whiten, detect, report |
| `lw_precision(X)` | Ledoit-Wolf precision matrix |
| `gl_precision(X)` | Graphical Lasso precision matrix (CV) |
| `whiten(X, Theta)` | Whiten X using precision matrix |
| `canonical_corr(X, Y, ix, iy)` | Canonical correlation on identified subsets |
| `generate_data(n, ...)` | Simulate block-structured biclique data |
| `eval_recovery(ix_hat, iy_hat, ...)` | Sensitivity, specificity, near-exact, Jaccard |
| `run_scenario(methods, n_list, ...)` | Run full simulation scenario |

---

## Dependencies

- **R** >= 4.0
- **glasso** (optional; required for method = "glasso")
- **MASS** (for mvrnorm in simulation, optional)

---

## Citation

If you use pgCCA, please cite:

```
@article{pgcca2025,
  title   = {Precision-Guided Graph Canonical Correlation Analysis},
  author  = {[Authors]},
  journal = {Statistics and Innovation},
  year    = {2025},
  note    = {Under review. arXiv:[XXXX.XXXXX]}
}
```

---

## License

MIT License. See LICENSE file.
