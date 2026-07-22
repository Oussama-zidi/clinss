# core-quantiles.R
# Shared building blocks used by every procedure.
# Almost every closed-form sample size formula reduces to a function of
# z_{1-alpha} (or z_{1-alpha/2}) and z_{1-beta}.  Writing them once means
# a fix propagates everywhere automatically.

#' @keywords internal
z_alpha <- function(alpha, sides = 2L) {
  stats::qnorm(1 - alpha / sides)
}

#' @keywords internal
z_beta <- function(power) {
  stats::qnorm(power)
}

#' Inflate a sample size for anticipated dropout
#'
#' Returns \code{ceiling(n / (1 - dropout))}.
#'
#' @param n Evaluable sample size (single positive integer).
#' @param dropout Dropout proportion in [0, 1).
#'
#' @references
#' Chow, S.C., Shao, J., Wang, H. (2008). \emph{Sample Size Calculations
#' in Clinical Research}, 2nd ed. Chapman & Hall/CRC. Section 1.2.5.
#'
#' @keywords internal
adjust_dropout <- function(n, dropout) {
  if (dropout <= 0) return(n)
  ceiling(n / (1 - dropout))
}

# The variance of the difference in means under unequal allocation is
# sigma^2 * (1/n1 + 1/n2) = sigma^2 / n1 * (1 + 1/ratio).
# This factor (1 + 1/ratio) appears in every two-sample closed-form formula.
#' @keywords internal
allocation_factor <- function(ratio) {
  1 + 1 / ratio
}
