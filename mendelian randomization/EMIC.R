library(MASS)
library(data.table)
library(dplyr)

set.seed(1234)

calculate_eta_em <- function(eta, eta_cor, w, n) {
  order <- order(eta)
  eta <- eta[order]
  eta_cor <- eta_cor[order,][, order]
  w <- w[order]

  ne <- c()
  for (i in seq(n)) {
    eigenvals <- eigen(eta_cor[seq(i), seq(i)])$values
    eigenvals[eigenvals > 1] <- 1
    ne <- c(ne, sum(eigenvals))
  }
  ne1 <- c(0, ne[seq(n - 1)])
  delta <- ne - ne1

  tau <- delta * w
  tau <- tau / sum(tau)

  s <- c()
  for (i in seq(n)) s <- c(s, sum(tau[seq(i)]))
  p <- c()
  for (i in seq(n)) p <- c(p, (s[i] - tau[i] / 2))

  t <- length(p[p < 0.5])
  k <- t + 1
  extrapolate <- (0.5 - p[t]) / (p[k] - p[t])
  eta_em <- eta[t] + (eta[k] - eta[t]) * extrapolate
  return(list(eta_em = eta_em, delta = delta))
}

emic <- function(beta_y, sigma_y, beta_x, sigma_x, genotype_cor) {
  boot_n <- 1000
  c <- 1.1

  eta <- beta_y / beta_x
  eta_cor <- 0.4428 * (genotype_cor^2) + 0.5665 * abs(genotype_cor)
  diag(eta_cor) <- 1
  w <- (beta_x / sigma_y)^4
  n <- length(beta_y)
  eta_em_result <- calculate_eta_em(eta, eta_cor, w, n)
  eta_em <- eta_em_result$eta_em
  delta <- eta_em_result$delta

  beta_y_cov <- sigma_y %o% sigma_y * genotype_cor
  beta_x_cov <- sigma_x %o% sigma_x * genotype_cor
  beta_y_boot <- mvrnorm(boot_n, beta_y, beta_y_cov)
  beta_x_boot <- mvrnorm(boot_n, beta_x, beta_x_cov)
  eta_boot <- beta_y_boot / beta_x_boot

  eta_em_boot <- c()
  for (i in seq(boot_n)) {
    eta_em_boot <- c(eta_em_boot, calculate_eta_em(eta_boot[i,], eta_cor, w, n)$eta_em)
  }
  sigma_em <- sd(eta_em_boot)
pval <- 2 * pnorm(-abs((eta_em / sigma_em * c)))

# --- Add Cochran's Q ---
# ratio estimates
eta <- beta_y / beta_x

# weights (inverse variance approx)
w_Q <- beta_x^2 / sigma_y^2

# weighted median function
weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  w <- w / sum(w)
  cum_w <- cumsum(w)
  return(x[which(cum_w >= 0.5)[1]])
}

# weighted median causal estimate
eta_wm <- weighted_median(eta, w_Q)

# Cochran’s Q with delta adjustment (EMIC definition)
Q <- sum((delta^2) * w_Q * (eta - eta_wm)^2)
Q_pval <- pchisq(Q, df = n - 1, lower.tail = FALSE)


return(c(n, eta_em, sigma_em, pval, Q, Q_pval))

}
