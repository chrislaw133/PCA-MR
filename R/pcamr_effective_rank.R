#' Effective rank of a positive spectrum
#'
#' @param vals Numeric vector of positive eigenvalues.
#'
#' @return Effective rank.
#' @keywords internal
.pcamr_effective_rank <- function(vals) {
  vals <- as.numeric(vals)
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) return(NA_real_)
  p_i <- vals / sum(vals)
  H <- -sum(p_i * log(p_i))
  exp(H)
}
