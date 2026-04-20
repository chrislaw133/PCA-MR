pcamr_fstat <- function(bx, bxse, ld, N_X, keep_var = 0.999, min_modes = 1) {
  bx   <- as.numeric(bx)
  bxse <- as.numeric(bxse)
  ld   <- as.matrix(ld)

  eig <- svd(ld)
  Q <- eig$u
  lambda <- eig$d

  bx_pc <- as.numeric(crossprod(Q, bx))
  SigX <- (bxse %o% bxse) * ld
  SigX_pc <- t(Q) %*% SigX %*% Q

  cumvar <- cumsum(lambda) / sum(lambda)
  K0 <- which(cumvar >= keep_var)[1]
  if (is.na(K0)) K0 <- length(lambda)

  keep <- seq_len(K0)
  if (length(keep) < min_modes) {
    keep <- seq_len(min(min_modes, length(lambda)))
  }

  se_bx <- sqrt(diag(SigX_pc[keep, keep, drop = FALSE]))

  R2 <- bx_pc[keep]^2 / (bx_pc[keep]^2 + N_X * se_bx^2)
  Fpc <- (N_X - 2) * R2 / pmax(1 - R2, 1e-12)

  list(
    keep = keep,
    bx_pc = bx_pc[keep],
    se_bx_pc = se_bx,
    R2 = R2,
    Fpc = Fpc,
    mean_F = mean(Fpc, na.rm = TRUE),
    min_F = min(Fpc, na.rm = TRUE),
    median_F = median(Fpc, na.rm = TRUE)
  )
}
