#' PCA-MR: PCA-based MR using correlated instruments
#'
#' Performs PCA-based Mendelian randomization using spectral decomposition of LD
#' and multiplicative random-effects IVW in rotated space
#'
#' @param bx Numeric vector of SNP exposure effects.
#' @param by Numeric vector of SNP outcome effects.
#' @param sey Numeric vector of standard errors of outcome effects
#' @param ld SNP correlation matrix
#' @param neff Vector of effective sample sizes for each SNP. If unavailable, use the GWAS sample size for every SNP
#' @param overlap_frac Proportion of SNP-level sample overlap (default=1).
#'
#' @return A named list with slope, se, p, Q_pval, n_modes, var_expl.
#'
#' @export
pca_mr <- function(bx, by, sey, ld, neff, overlap_frac = 1) {

  if (length(neff) != length(sey))
    stop(sprintf("Length mismatch in neff vs sey: %d vs %d", length(neff), length(sey)))

  if (!is.numeric(neff))
    stop("neff must be numeric")
  
  ld <- (ld + t(ld)) / 2
  diag(ld) <- 1
  eig <- eigen(ld, symmetric = TRUE)
  Q <- eig$vectors
  Lambda <- pmax(eig$values, 1e-8)

  N_mat  <- outer(as.numeric(neff), as.numeric(neff),
                  function(Ni, Nj) sqrt(Ni * Nj))
  Nij_mat <- outer(as.numeric(neff), as.numeric(neff),
                   function(Ni, Nj) overlap_frac * pmin(Ni, Nj))
  overlap_scale <- Nij_mat / N_mat

  cov_gwas <- (sey %o% sey) * ld * overlap_scale
  cov_gwas_tilde <- t(Q) %*% cov_gwas %*% Q
  se_gwas_tilde <- sqrt(diag(cov_gwas_tilde))

  bx_tilde <- as.numeric(crossprod(Q, bx))
  by_tilde <- as.numeric(crossprod(Q, by))

  var_frac <- Lambda / sum(Lambda)
  keep <- which(cumsum(var_frac) <= 0.99)
  if (length(keep) == 0) keep <- which.max(var_frac)

  bx_tilde <- bx_tilde[keep]
  by_tilde <- by_tilde[keep]
  se_gwas_tilde <- se_gwas_tilde[keep]
  Lambda_keep <- Lambda[keep]

  if (length(Lambda_keep) <= 1) {
    return(list(
      slope = NA,
      se = NA,
      p = NA,
      Q_pval = NA,
      n_modes = NA,
      var_expl = NA
    ))
  }

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

  return(list(
    slope = slope,
    se = se,
    p = p,
    Q_pval = Q_pval,
    n_modes = length(keep),
    var_expl = 100 * sum(var_frac[keep])
  ))
}
