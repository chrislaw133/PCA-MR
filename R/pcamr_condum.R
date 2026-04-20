#' Condition number
#'
#' @param M Numeric matrix.
#'
#' @return Condition number.
#' @keywords internal
.pcamr_condnum <- function(M) {
  ev <- tryCatch(
    svd(M, nu = 0, nv = 0)$d,
    error = function(e) NA_real_
  )
  if (any(!is.finite(ev))) return(Inf)
  ev <- pmax(ev, 1e-12)
  max(ev) / min(ev)
}
