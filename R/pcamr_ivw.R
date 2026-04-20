#' PCA-MR-IVW
#'
#' @param bx SNP-exposure effects.
#' @param by SNP-outcome effects.
#' @param bxse Standard errors of bx.
#' @param byse Standard errors of by.
#' @param ld LD correlation matrix.
#' @param keep_var Cumulative variance explained threshold for PC retention.
#' @param min_modes Minimum number of retained PCs.
#' @param kappa_max Maximum allowed condition number.
#' @param model Either "random" or "fixed".
#' @param phi_floor Lower bound for random-effects dispersion.
#'
#' @return A list containing the PCA-MR-IVW estimate and diagnostics.
#' @export
pcamr_ivw <- function(bx, by, bxse, byse, ld,
                      keep_var = 0.999,
                      min_modes = 2,
                      kappa_max = 100,
                      model = c("random", "fixed"),
                      phi_floor = 1) {
  model <- match.arg(model)

  pc <- .pcamr_prepare(
    bx = bx, by = by, bxse = bxse, byse = byse, ld = ld,
    keep_var = keep_var,
    min_modes = min_modes,
    kappa_max = kappa_max
  )
  if (is.null(pc)) return(NULL)

  Omega_inv <- tryCatch(solve(pc$Omega), error = function(e) NULL)
  if (is.null(Omega_inv)) return(NULL)

  Lambda <- diag(pc$lambda, nrow = length(pc$lambda), ncol = length(pc$lambda))
  W <- Lambda %*% Omega_inv

  denom <- as.numeric(t(pc$bx) %*% W %*% pc$bx)
  if (!is.finite(denom) || denom <= 0) return(NULL)

  slope <- as.numeric((t(pc$bx) %*% W %*% pc$by) / denom)

  resid <- as.numeric(pc$by - slope * pc$bx)
  Q <- as.numeric(t(resid) %*% W %*% resid)
  df <- length(pc$bx) - 1
  if (!is.finite(df) || df <= 0) return(NULL)

  phi <- if (model == "fixed") 1 else max(phi_floor, Q / df)
  se <- sqrt(phi / denom)
  p <- 2 * pnorm(-abs(slope / se))
  Q_pvalue <- pchisq(Q, df = df, lower.tail = FALSE)

  list(
    method = "PCA-MR-IVW",
    theta = slope,
    se = se,
    p = p,
    Q = Q,
    df = df,
    phi = phi,
    rse = sqrt(Q / df),
    model = model,
    Q_pvalue = Q_pvalue,
    n_modes = pc$n_modes,
    var_explained = pc$var_explained,
    kappa = pc$kappa,
    erank = pc$erank,
    lambda = pc$lambda,
    keep = pc$keep
  )
}
