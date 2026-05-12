# =============================================================================
# pgCCA — Stage 3a: Greedy Biclique Detection (Park et al. 2025)
# File:   R/biclique.R
# =============================================================================

#' Greedy Bipartite Dense Subgraph Extraction
#'
#' Internal function implementing the Park et al. (2025) greedy biclique
#' elimination algorithm. At each step removes the lowest-degree vertex
#' (row or column of W) and records the density of the remaining subgraph.
#' Returns the subgraph with maximum density.
#'
#' @param W      Non-negative adjacency matrix (p x q).
#' @param lambda Density-vs-size trade-off parameter in (0, 1).
#' @param index  List with $x (row indices) and $y (col indices) active
#'               in the current submatrix.
#' @return List with $s_out (row indices kept) and $t_out (col indices kept).
#' @keywords internal
greedy_bipar3 <- function(W, lambda, index) {
  W1  <- W
  deleteS <- list(); deleteT <- list()
  nr  <- nrow(W); nc <- ncol(W)
  dens <- matrix(0, nrow = nr + nc - 1, ncol = 8)
  R   <- rowSums(W); C <- colSums(W)
  rs  <- integer(0); cs <- integer(0)

  for (i in seq_len(nr + nc - 2)) {
    if (length(rs) > nr - 2 || length(cs) > nc - 2) break

    row <- nr - length(rs); col <- nc - length(cs)
    IS  <- order(R)[(length(rs) + 1):(length(rs) + 1)]
    IT  <- order(C)[(length(cs) + 1):(length(cs) + 1)]
    dS  <- mean(R[IS]); dT <- mean(C[IT])

    if ((col / row) * dS <= dT) {
      C   <- C - if (length(IS) > 1) colSums(W[IS, ]) else W[IS, ]
      rs  <- c(rs, IS); W[IS, ] <- 0; R[IS] <- 0
      dens[i, 1] <- sum(R != 0); dens[i, 2] <- sum(C != 0); dens[i, 3] <- 1
    } else {
      R   <- R - if (length(IT) > 1) rowSums(W[, IT]) else W[, IT]
      cs  <- c(cs, IT); W[, IT] <- 0; C[IT] <- 0
      dens[i, 1] <- row; dens[i, 2] <- col; dens[i, 3] <- 2
    }
    deleteS[[i]] <- IS; deleteT[[i]] <- IT
    d12 <- dens[i, 1] * dens[i, 2]
    dens[i, 8] <- if (is.na(d12) || d12 == 0) -Inf else sum(R) / d12^lambda
  }

  dens[is.nan(dens[, 8]), 8] <- -Inf
  best <- which.max(dens[, 8])
  rm_s <- integer(0); rm_t <- integer(0)
  if (length(best) > 0 && best > 0)
    for (j in seq_len(best)) {
      if (dens[j, 3] == 1) rm_s <- union(rm_s, deleteS[[j]])
      else                  rm_t <- union(rm_t, deleteT[[j]])
    }

  base_dens <- sum(W1) / (nrow(W1) * ncol(W1))^lambda
  if (base_dens > max(dens[, 8], na.rm = TRUE))
    list(s_out = c(), t_out = c())
  else
    list(s_out = index$x[rm_s], t_out = index$y[rm_t])
}


#' Extract Dense Biclique Subgraph (Park et al. 2025)
#'
#' Applies the greedy biclique elimination for a given lambda.
#'
#' @param W       Non-negative p x q correlation matrix (thresholded).
#' @param lambda  Density parameter in (0.5, 0.85).
#' @param num_subg Number of nested subgraphs to extract. Default 1.
#' @return List of subgraphs; each element has $x and $y index sets.
#' @references Park et al. (2025) arXiv:2502.01780.
#' @keywords internal
subg <- function(W, lambda = 0.5, num_subg = 1) {
  W0 <- W; result <- list(); i <- 1; subg0 <- list(); index <- list()
  index[[i]] <- list(x = seq_len(nrow(W0)), y = seq_len(ncol(W0)))
  result[[i]] <- greedy_bipar3(W0, lambda = lambda, index = index[[i]])

  while (i <= num_subg &&
         length(result[[i]]$s_out) > 1 &&
         length(result[[i]]$t_out) > 1 &&
         length(result[[i]]$s_out) < length(index[[i]]$x) &&
         length(result[[i]]$t_out) < length(index[[i]]$y)) {
    subg0[[i]] <- list(
      x = index[[i]]$x[!index[[i]]$x %in% result[[i]]$s_out],
      y = index[[i]]$y[!index[[i]]$y %in% result[[i]]$t_out]
    )
    index[[i + 1]] <- list(x = result[[i]]$s_out, y = result[[i]]$t_out)
    W11 <- W0[index[[i + 1]]$x, index[[i + 1]]$y]
    if (sum(W11) == 0) break
    result[[i + 1]] <- greedy_bipar3(W = W11, lambda = lambda,
                                      index = index[[i + 1]])
    if (length(index[[i+1]]$x[!index[[i+1]]$x %in% result[[i+1]]$s_out]) <= 2 &&
        length(index[[i]]$y[!index[[i]]$y %in% result[[i]]$t_out]) <= 2) break
    i <- i + 1
  }
  subg0
}


#' KL-Divergence Criterion for Lambda Selection
#'
#' Computes the KL divergence between within-biclique and background
#' correlation densities. Used to select the optimal lambda.
#'
#' @param W     Binary adjacency matrix (0/1), p x q.
#' @param index List of subgraph index sets (from subg()).
#' @return Scalar KL divergence value (-Inf if invalid).
#' @references Park et al. (2025) arXiv:2502.01780.
#' @keywords internal
kld_cross <- function(W, index) {
  if (length(index) == 0) return(-Inf)
  p_  <- sum(W > 0) / (nrow(W) * ncol(W))
  W0  <- lapply(index, function(z) W[z$x, z$y])
  p1  <- mean(unlist(W0)); n_in <- length(unlist(W0))
  p0  <- (sum(W > 0) - sum(unlist(W0) > 0)) / (nrow(W) * ncol(W) - n_in)

  if (length(index[[1]]$x) == 0 || length(index[[1]]$y) == 0)   return(-Inf)
  if (length(index[[1]]$x) + length(index[[1]]$y) >=
      nrow(W) + ncol(W)) return(-Inf)

  p1 <- pmax(pmin(p1, 1 - 1e-9), 1e-9)
  p0 <- pmax(pmin(p0, 1 - 1e-9), 1e-9)
  p_ <- pmax(pmin(p_,  1 - 1e-9), 1e-9)

  sum(unlist(W0) > 0)   * p1 * log(p1 / p_) +
  (n_in - sum(unlist(W0))) * (1 - p1) * log((1 - p1) / p_) +
  (sum(W) - sum(unlist(W0))) * p0 * log(p0 / p_) +
  (length(which(W == 0)) - length(which(unlist(W0) == 0))) *
    (1 - p0) * log((1 - p0) / p_)
}


#' Run gCCA Greedy Biclique Detection with Automatic Lambda Selection
#'
#' Applies the Park et al. (2025) algorithm to a (possibly whitened)
#' cross-correlation matrix. The lambda grid is swept and the value
#' maximising the KL divergence is chosen automatically.
#'
#' @param X         n x p matrix (possibly whitened).
#' @param Y         n x q matrix (possibly whitened).
#' @param eps       Detection threshold. Default 0.08.
#' @param lambda_seq Lambda grid to search. Default seq(0.5, 0.85, by = 0.05).
#' @return List with:
#'   \describe{
#'     \item{ix}{Integer vector of selected X-variable indices.}
#'     \item{iy}{Integer vector of selected Y-variable indices.}
#'     \item{rho}{Canonical correlation (on the matrices passed in).}
#'     \item{lambda_opt}{Optimal lambda chosen by KL divergence.}
#'   }
#' @references Park et al. (2025) arXiv:2502.01780.
#' @export
run_gcca <- function(X, Y,
                     eps        = 0.08,
                     lambda_seq = seq(0.5, 0.85, by = 0.05)) {
  W00 <- abs(cor(X, Y)); W00[W00 < eps] <- 0
  Wb  <- (W00 > 0) * 1L

  res_all <- lapply(lambda_seq, function(lam) subg(W00, lambda = lam, 1))
  klds    <- sapply(res_all, function(z) kld_cross(Wb, z))
  best    <- res_all[[which.max(klds)]]
  lam_opt <- lambda_seq[which.max(klds)]

  if (length(best) == 0 || length(best[[1]]$x) == 0 || length(best[[1]]$y) == 0)
    return(list(ix = integer(0), iy = integer(0), rho = 0, lambda_opt = NA))

  ix  <- sort(unique(unlist(lapply(best, `[[`, "x"))))
  iy  <- sort(unique(unlist(lapply(best, `[[`, "y"))))
  rho <- canonical_corr(X, Y, ix, iy)

  list(ix = ix, iy = iy, rho = rho, lambda_opt = lam_opt)
}
