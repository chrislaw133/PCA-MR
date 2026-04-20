#' GLS Egger in projected PC space
#'
#' @param bx_tilde Projected exposure effects.
#' @param by_tilde Projected outcome effects.
#' @param Omega Projected outcome covariance matrix.
#'
#' @return Egger fit object.
#' @keywords internal
.pcamr_gls_egger <- function(bx_tilde, by_tilde, Omega, Lambda) {
  n <- length(bx_tilde)

  if (length(by_tilde) != n) stop("Length mismatch between bx_tilde and by_tilde.")
  if (!is.matrix(Omega) || any(dim(Omega) != c(n, n))) {
    stop("Omega has incorrect dimensions.")
  }

  if (n < 3) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      se_intercept = NA_real_,
      se_slope = NA_real_,
      p_intercept = NA_real_,
      p_slope = NA_real_,
      Q = NA_real_,
      df = NA_real_,
      rse = NA_real_,
      Q_pval = NA_real_
    ))
  }

  sgn <- ifelse(bx_tilde >= 0, 1, -1)
  bx_tilde <- abs(bx_tilde)
  by_tilde <- sgn * by_tilde
  S <- diag(sgn, nrow = n, ncol = n)
  Omega <- S %*% Omega %*% S

  X <- cbind(1, bx_tilde)
  y <- matrix(by_tilde, ncol = 1)

  Omega_inv <- tryCatch(solve(Omega), error = function(e) NULL)
  if (is.null(Omega_inv)) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      se_intercept = NA_real_,
      se_slope = NA_real_,
      p_intercept = NA_real_,
      p_slope = NA_real_,
      Q = NA_real_,
      df = NA_real_,
      rse = NA_real_,
      Q_pval = NA_real_
    ))
  }

  Omega_inv <- Lambda %*% Omega_inv

  XtOiX <- t(X) %*% Omega_inv %*% X
  bread <- tryCatch(solve(XtOiX), error = function(e) NULL)
  if (is.null(bread)) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      se_intercept = NA_real_,
      se_slope = NA_real_,
      p_intercept = NA_real_,
      p_slope = NA_real_,
      Q = NA_real_,
      df = NA_real_,
      rse = NA_real_,
      Q_pval = NA_real_
    ))
  }

  beta_hat <- bread %*% t(X) %*% Omega_inv %*% y
  resid <- y - X %*% beta_hat

  Q <- as.numeric(t(resid) %*% Omega_inv %*% resid)
  df <- n - 2

  if (!is.finite(df) || df <= 0) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      se_intercept = NA_real_,
      se_slope = NA_real_,
      p_intercept = NA_real_,
      p_slope = NA_real_,
      Q = Q,
      df = df,
      rse = NA_real_,
      Q_pval = NA_real_
    ))
  }

  rse <- sqrt(Q / df)
  infl <- max(1, rse)

  se_intercept <- sqrt(bread[1, 1]) * infl
  se_slope <- sqrt(bread[2, 2]) * infl

  intercept_hat <- as.numeric(beta_hat[1, 1])
  slope_hat <- as.numeric(beta_hat[2, 1])

  t_int <- intercept_hat / se_intercept
  t_slope <- slope_hat / se_slope

  p_intercept <- 2 * pt(abs(t_int), df = df, lower.tail = FALSE)
  p_slope <- 2 * pt(abs(t_slope), df = df, lower.tail = FALSE)
  Q_pval <- pchisq(Q, df = df, lower.tail = FALSE)

  list(
    intercept = intercept_hat,
    slope = slope_hat,
    se_intercept = se_intercept,
    se_slope = se_slope,
    p_intercept = p_intercept,
    p_slope = p_slope,
    Q = Q,
    df = df,
    rse = rse,
    Q_pval = Q_pval
  )
}
