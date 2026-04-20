#' Validate PCA-MR inputs
#'
#' @param bx SNP-exposure effects.
#' @param by SNP-outcome effects.
#' @param bxse Standard errors of bx.
#' @param byse Standard errors of by.
#' @param ld LD correlation matrix.
#'
#' @return A validated list of inputs.
#' @keywords internal
.pcamr_validate_inputs <- function(bx, by, bxse, byse, ld) {
  bx   <- as.numeric(bx)
  by   <- as.numeric(by)
  bxse <- as.numeric(bxse)
  byse <- as.numeric(byse)
  ld   <- as.matrix(ld)

  p <- length(bx)

  if (length(by) != p || length(bxse) != p || length(byse) != p) {
    stop("bx, by, bxse, and byse must all have the same length.")
  }

  if (!all(dim(ld) == c(p, p))) {
    stop("ld must be a square matrix with dimensions length(bx) x length(bx).")
  }

  if (any(!is.finite(bx)) || any(!is.finite(by)) ||
      any(!is.finite(bxse)) || any(!is.finite(byse)) ||
      any(!is.finite(ld))) {
    stop("Inputs contain non-finite values.")
  }

  if (any(bxse <= 0) || any(byse <= 0)) {
    stop("bxse and byse must be strictly positive.")
  }

  list(
    bx = bx,
    by = by,
    bxse = bxse,
    byse = byse,
    ld = ld
  )
}
