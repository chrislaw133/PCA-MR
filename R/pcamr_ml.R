#' PCA-MR-ML
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
#' @param alpha Confidence interval tail probability.
#'
#' @return A list containing the PCA-MR-ML estimate and diagnostics.
#' @export
pcamr_ml <- function(bx, by, bxse, byse, ld,
                     keep_var = 0.999,
                     min_modes = 1,
                     kappa_max = 100,
                     model = c("random", "fixed"),
                     alpha = 0.05) {
  model <- match.arg(model)

  pc <- .pcamr_prepare(
    bx = bx, by = by, bxse = bxse, byse = byse, ld = ld,
    keep_var = keep_var,
    min_modes = min_modes,
    kappa_max = kappa_max
  )
  if (is.null(pc)) return(NULL)

  bx_pc <- as.numeric(pc$bx)
  by_pc <- as.numeric(pc$by)
  SigmaX <- pc$SigmaX
  SigmaY <- pc$SigmaY

  k <- length(bx_pc)
  if (k < 1) return(NULL)
  if (model == "random" && k < 2) return(NULL)

  Lambda <- diag(pc$lambda, nrow = length(pc$lambda), ncol = length(pc$lambda))

  Taux <- tryCatch(solve(SigmaX), error = function(e) NULL)
  Tauy <- tryCatch(solve(SigmaY), error = function(e) NULL)
  if (is.null(Taux) || is.null(Tauy)) return(NULL)

  Taux <- Lambda %*% Taux
  Tauy <- Lambda %*% Tauy

  init_ivw <- pcamr_ivw(
    bx = bx, by = by, bxse = bxse, byse = byse, ld = ld,
    keep_var = keep_var,
    min_modes = min_modes,
    kappa_max = kappa_max,
    model = model
  )
  theta0 <- if (!is.null(init_ivw) && is.finite(init_ivw$slope)) init_ivw$slope else 0

  loglikelihoodcorrel <- function(param, x, Taux, y, Tauy) {
    k <- length(x)
    b <- param[1:k]
    theta <- param[k + 1]

    as.numeric(
      0.5 * t(x - b) %*% Taux %*% (x - b) +
      0.5 * t(y - theta * b) %*% Tauy %*% (y - theta * b)
    )
  }

  par0 <- c(bx_pc, theta0)

  opt <- tryCatch(
    optim(
      par = par0,
      fn = loglikelihoodcorrel,
      x = bx_pc,
      Taux = Taux,
      y = by_pc,
      Tauy = Tauy,
      hessian = TRUE,
      method = "BFGS",
      control = list(maxit = 25000, reltol = 1e-10)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || is.null(opt$par) || is.null(opt$hessian)) return(NULL)
  if (!is.finite(opt$value)) return(NULL)

  V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  if (is.null(V)) return(NULL)

  theta <- as.numeric(opt$par[k + 1])
  theta_var_fixed <- as.numeric(V[k + 1, k + 1])
  if (!is.finite(theta_var_fixed) || theta_var_fixed < 0) return(NULL)

  if (model == "fixed") {
    theta_se <- sqrt(theta_var_fixed)
  } else {
    phi <- max(2 * opt$value / (k - 1), 1)
    theta_se <- sqrt(theta_var_fixed * phi)
  }

  if (!is.finite(theta_se) || theta_se <= 0) return(NULL)

  zval <- theta / theta_se
  pval <- 2 * pnorm(-abs(zval))
  ci_lo <- theta - qnorm(1 - alpha / 2) * theta_se
  ci_hi <- theta + qnorm(1 - alpha / 2) * theta_se

  heter_stat <- 2 * opt$value
  heter_df <- k - 1
  heter_p <- if (heter_df > 0) pchisq(heter_stat, df = heter_df, lower.tail = FALSE) else NA_real_
  rse <- if (heter_df > 0) sqrt(heter_stat / heter_df) else NA_real_

  list(
    method = "PCA-MR-ML",
    theta = theta,
    se = theta_se,
    p = pval,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    K = k,
    model = model,
    rse = rse,
    heter_stat = heter_stat,
    heter_df = heter_df,
    heter_p = heter_p,
    conv = opt$convergence,
    value = opt$value,
    n_modes = pc$n_modes,
    kappa = pc$kappa,
    erank = pc$erank,
    lambda = pc$lambda,
    keep = pc$keep
  )
}
