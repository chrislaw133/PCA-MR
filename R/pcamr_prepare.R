#' Prepare projected PC-space quantities for PCA-MR
#'
#' @param bx SNP-exposure effects.
#' @param by SNP-outcome effects.
#' @param bxse Standard errors of bx.
#' @param byse Standard errors of by.
#' @param ld LD correlation matrix.
#' @param keep_var Cumulative variance explained threshold for PC retention.
#' @param min_modes Minimum number of retained PCs.
#' @param kappa_max Maximum allowed condition number.
#'
#' @return A list of projected quantities used by PCA-MR estimators.
#' @keywords internal
.pcamr_prepare <- function(bx, by, bxse, byse, ld,
                           keep_var = 0.999,
                           min_modes = 1,
                           kappa_max = 100) {
  args <- .pcamr_validate_inputs(bx, by, bxse, byse, ld)
  bx   <- args$bx
  by   <- args$by
  bxse <- args$bxse
  byse <- args$byse
  ld   <- args$ld

  eig <- svd(ld)
  Q <- eig$u
  lambda <- eig$d

  bx_pc <- as.numeric(crossprod(Q, bx))
  by_pc <- as.numeric(crossprod(Q, by))

  SigX <- (bxse %o% bxse) * ld
  SigY <- (byse %o% byse) * ld
  SigX_pc <- t(Q) %*% SigX %*% Q
  SigY_pc <- t(Q) %*% SigY %*% Q

  cumvar <- cumsum(lambda) / sum(lambda)
  K0 <- which(cumvar >= keep_var)[1]
  if (is.na(K0)) K0 <- length(lambda)

  keep <- seq_len(K0)

  if (length(keep) < min_modes) {
    keep <- seq_len(min(min_modes, length(lambda)))
  }
  if (length(keep) < min_modes) return(NULL)

  Omega <- SigY_pc[keep, keep, drop = FALSE]
  kappa <- .pcamr_condnum(Omega)

  while (length(keep) > min_modes && is.finite(kappa) && kappa > kappa_max) {
    keep <- keep[-length(keep)]
    Omega <- SigY_pc[keep, keep, drop = FALSE]
    kappa <- .pcamr_condnum(Omega)
  }

  if (length(keep) < min_modes) return(NULL)
  if (!is.finite(kappa) || kappa > kappa_max) return(NULL)

  cumvar_explained <- sum(lambda[keep]) / sum(lambda)

  SigmaX <- SigX_pc[keep, keep, drop = FALSE]
  SigmaY <- SigY_pc[keep, keep, drop = FALSE]

  list(
    bx = bx_pc[keep],
    by = by_pc[keep],
    bxse = sqrt(diag(SigmaX)),
    byse = sqrt(diag(SigmaY)),
    Omega = Omega,
    SigmaX = SigmaX,
    SigmaY = SigmaY,
    keep = keep,
    n_modes = length(keep),
    var_explained = cumvar_explained,
    kappa = kappa,
    lambda = lambda[keep],
    erank = .pcamr_effective_rank(lambda[keep])
  )
}
