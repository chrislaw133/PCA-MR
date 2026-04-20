#' PCA-MR-Egger
#'
#' @param bx SNP-exposure effects.
#' @param by SNP-outcome effects.
#' @param bxse Standard errors of bx.
#' @param byse Standard errors of by.
#' @param ld LD correlation matrix.
#' @param keep_var Cumulative variance explained threshold for PC retention.
#' @param min_modes Minimum number of retained PCs. Defaults to 3.
#' @param kappa_max Maximum allowed condition number.
#'
#' @return A list containing the PCA-MR-Egger estimate and diagnostics.
#' @export
pcamr_egger <- function(bx, by, bxse, byse, ld,
                        keep_var = 0.999,
                        min_modes = 3,
                        kappa_max = 100) {
  pc <- .pcamr_prepare(
    bx = bx, by = by, bxse = bxse, byse = byse, ld = ld,
    keep_var = keep_var,
    min_modes = min_modes,
    kappa_max = kappa_max
  )
  if (is.null(pc)) return(NULL)
  if (pc$n_modes < 3) return(NULL)

  fit <- .pcamr_gls_egger(
    bx_tilde = pc$bx,
    by_tilde = pc$by,
    Omega = pc$Omega,
    Lambda = pc$lambda
  )

  list(
    method = "PCA-MR-Egger",
    intercept = fit$intercept,
    slope = fit$slope,
    se_intercept = fit$se_intercept,
    se_slope = fit$se_slope,
    p_intercept = fit$p_intercept,
    p_slope = fit$p_slope,
    Q = fit$Q,
    df = fit$df,
    rse = fit$rse,
    Q_pvalue = fit$Q_pval,
    n_modes = pc$n_modes,
    var_explained = pc$var_explained,
    kappa = pc$kappa,
    erank = pc$erank,
    lambda = pc$lambda,
    keep = pc$keep
  )
}
