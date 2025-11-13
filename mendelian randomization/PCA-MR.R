pca_mr <- function(bx, by, sey, ld, neff, overlap_frac = 1) {
  #Sanity check
  if (length(neff) != length(sey))
    stop(sprintf("Length mismatch in neff vs sey: %d vs %d", length(neff), length(sey)))
  if (!is.numeric(neff))
    stop("neff must be numeric")

#Spectral decomposition
  eig <- eigen(ld, symmetric = TRUE)
  Q <- eig$vectors
  Lambda <- pmax(eig$values, 1e-8)

  #Pairwise SNP sample overlap scaling
  N_mat  <- outer(as.numeric(neff), as.numeric(neff),
                  function(Ni, Nj) overlap_frac * sqrt(Ni * Nj))
  Nij_mat <- outer(as.numeric(neff), as.numeric(neff),
                   function(Ni, Nj) overlap_frac * pmin(Ni, Nj))
  overlap_scale <- Nij_mat / N_mat
  stopifnot(is.matrix(overlap_scale), all(dim(overlap_scale) == dim(ld)))

  #Approximate sampling covariance matrix of GWAS outcome effect sizes
  cov_gwas <- (sey %o% sey) * ld * overlap_scale
  cov_gwas_tilde <- t(Q) %*% cov_gwas %*% Q
  se_gwas_tilde <- sqrt(diag(cov_gwas_tilde))

  #Rotate effect sizes
  bx_tilde <- as.numeric(crossprod(Q, bx))
  by_tilde <- as.numeric(crossprod(Q, by))

  #Truncate PCs (We recommend staying between 95-99%)
  var_frac <- Lambda / sum(Lambda)
  keep <- which(cumsum(var_frac) <= 0.99)
  if (length(keep) == 0) keep <- which.max(var_frac)

  bx_tilde <- bx_tilde[keep]
  by_tilde <- by_tilde[keep]
  se_gwas_tilde <- se_gwas_tilde[keep]
  Lambda_keep <- Lambda[keep]

  Multiplicative random-effects IVW regression
  w <- Lambda_keep / (se_gwas_tilde^2)
  slope <- sum(w * bx_tilde * by_tilde) / sum(w * bx_tilde^2)

  resid <- by_tilde - slope * bx_tilde
  Qstat <- sum(w * resid^2)
  df <- length(bx_tilde) - 1
  se_fe <- sqrt(1 / sum(bx_tilde^2 / se_gwas_tilde^2))
  phi <- max(1, Qstat / df)
  se <- se_fe * sqrt(phi)
  p <- 2 * pnorm(-abs(slope / se))
  Q_pval <- pchisq(Qstat, df = df, lower.tail = FALSE)
              

  list(
    slope = slope,
    se = se,
    p = p,
    Q_pval = Q_pval, 
    n_modes = length(keep),
    var_expl = 100 * sum(var_frac[keep])
  )
}

results <- pca_mr_gene(bx, by, sey, ld, neff)
