# binary-two-proportions.R
#
# Sample size and power for the comparison of two independent proportions
# (parallel-group design, binary endpoint), on the difference scale.
#
# Hypotheses supported:
#   superiority         H0: p1 - p2  = 0        vs H1: p1 - p2 != 0 (or >)
#   noninferiority      H0: p1 - p2 <= -margin  vs H1: p1 - p2  > -margin
#   superiority_margin  H0: p1 - p2 <=  margin  vs H1: p1 - p2  >  margin
#   equivalence (TOST)  H0: |p1 - p2| >= margin vs H1: |p1 - p2| < margin
#
# Tests implemented:
#   superiority: two-sample z-test with variance pooled under H0
#     (equivalent to the uncorrected chi-square test).
#   margin-based hypotheses: Farrington-Manning score test, whose null
#     variance uses the restricted maximum likelihood estimates of
#     (p1, p2) under the constraint p1 - p2 = d0.
#
# Formula source (primary literature):
#   Farrington, C.P., Manning, G. (1990). Test statistics and sample size
#     formulae for comparative binomial trials with null hypothesis of
#     non-zero risk difference or non-unity relative risk.
#     Statistics in Medicine, 9(12), 1447-1454.
#   Fleiss, J.L., Levin, B., Paik, M.C. (2003). Statistical Methods for
#     Rates and Proportions, 3rd ed. Wiley. Chapter 4.
#   Chow, S.C., Shao, J., Wang, H. (2008). Sample Size Calculations in
#     Clinical Research, 2nd ed. Chapman & Hall/CRC. Chapter 4.
#
# Conventions:
#   Group 1 is the treatment group, group 2 the control/reference.
#   Higher proportions are assumed favourable. If your endpoint is one
#   where LOWER is better (e.g., failure, relapse), work with the
#   complementary event: p <- 1 - p for both groups.
#   ratio = n2 / n1; margin is always a positive number.

#' @keywords internal
REF_TWO_PROPS <- paste(
  "Farrington, C.P., Manning, G. (1990). Statistics in Medicine, 9(12),",
  "1447-1454;",
  "Fleiss, J.L., Levin, B., Paik, M.C. (2003). Statistical Methods for",
  "Rates and Proportions, 3rd ed. Wiley, Chapter 4;",
  "Chow, S.C., Shao, J., Wang, H. (2008). Sample Size Calculations in",
  "Clinical Research, 2nd ed. Chapman & Hall/CRC, Chapter 4.")

# ---------------------------------------------------------------------------
# Farrington-Manning restricted MLE (closed-form cubic solution)
# ---------------------------------------------------------------------------

# Restricted MLE of (p1, p2) under H0: p1 - p2 = d0, given working values
# (p1, p2) with weights (n1, n2). Closed-form solution of the cubic score
# equation from Farrington & Manning (1990), Section 2.
# Validated in tests/tests.R against direct numerical maximization of the
# restricted binomial likelihood.
#' @keywords internal
fm_restricted_mle <- function(p1, p2, n1, n2, d0) {
  if (d0 == 0) {
    pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
    return(c(p1t = pbar, p2t = pbar))
  }
  theta <- n2 / n1
  a  <- 1 + theta
  b  <- -(1 + theta + p1 + theta * p2 + d0 * (theta + 2))
  cc <- d0^2 + d0 * (2 * p1 + theta + 1) + p1 + theta * p2
  d  <- -p1 * d0 * (1 + d0)
  v  <- b^3 / (3 * a)^3 - b * cc / (6 * a^2) + d / (2 * a)
  u  <- sqrt(max(0, b^2 / (3 * a)^2 - cc / (3 * a)))
  if (v < 0) u <- -u
  w  <- (pi + acos(min(1, max(-1, v / u^3)))) / 3
  p1t <- 2 * u * cos(w) - b / (3 * a)
  c(p1t = p1t, p2t = p1t - d0)
}

# Null-hypothesis standard error of (p1hat - p2hat) at boundary d0, using
# the FM restricted MLEs evaluated at the true proportions.
#' @keywords internal
fm_null_se <- function(p1, p2, n1, n2, d0) {
  pt <- fm_restricted_mle(p1, p2, n1, n2, d0)
  sqrt(pt["p1t"] * (1 - pt["p1t"]) / n1 +
       pt["p2t"] * (1 - pt["p2t"]) / n2)
}

# ---------------------------------------------------------------------------
# Power engine
# ---------------------------------------------------------------------------

#' Power for the difference of two independent proportions
#'
#' Computes the power of a two-sample test comparing two independent
#' proportions at fixed sample sizes, on the risk difference scale.
#' Superiority uses the z-test with variance pooled under the null
#' (uncorrected chi-square test). Non-inferiority, superiority by a margin,
#' and equivalence use the Farrington-Manning score test, whose null
#' variance is evaluated at the restricted maximum likelihood estimates.
#'
#' Group 1 is the treatment group and higher proportions are assumed
#' favourable. If lower is better for your endpoint (e.g., relapse), use the
#' complementary event probabilities \code{1 - p} for both groups.
#'
#' @param n1 Sample size in group 1 (treatment).
#' @param n2 Sample size in group 2 (control). Defaults to \code{n1}.
#' @param p1 True proportion in group 1.
#' @param p2 True proportion in group 2.
#' @param alpha Significance level (one-sided for margin-based hypotheses).
#' @param hypothesis One of \code{"superiority"}, \code{"noninferiority"},
#'   \code{"superiority_margin"}, or \code{"equivalence"}.
#' @param margin Positive margin on the risk difference scale, required for
#'   all hypotheses except \code{"superiority"}.
#' @param sides 1 or 2. Used only for \code{hypothesis = "superiority"}.
#'
#' @return A single numeric value: the power at the specified sample sizes.
#'
#' @references
#' Farrington, C.P., Manning, G. (1990). Test statistics and sample size
#' formulae for comparative binomial trials with null hypothesis of
#' non-zero risk difference or non-unity relative risk.
#' \emph{Statistics in Medicine}, 9(12), 1447--1454.
#'
#' Fleiss, J.L., Levin, B., Paik, M.C. (2003). \emph{Statistical Methods
#' for Rates and Proportions}, 3rd ed. Wiley. Chapter 4.
#'
#' @examples
#' # Power of a non-inferiority comparison at n = 200 per group
#' power_two_proportions(n1 = 200, p1 = 0.85, p2 = 0.85, margin = 0.10,
#'                       alpha = 0.025, hypothesis = "noninferiority")
#'
#' @export
power_two_proportions <- function(n1, n2 = n1, p1, p2,
                                  alpha = 0.025,
                                  hypothesis = c("superiority",
                                                 "noninferiority",
                                                 "superiority_margin",
                                                 "equivalence"),
                                  margin = NULL,
                                  sides = if (hypothesis == "superiority") 2L else 1L) {
  hypothesis <- match.arg(hypothesis)
  check_positive(n1, "n1"); check_positive(n2, "n2")
  check_prob(p1, "p1"); check_prob(p2, "p2"); check_prob(alpha, "alpha")
  if (hypothesis != "superiority") {
    check_positive(margin, "margin")
    if (margin >= 1) stop("`margin` must be below 1 on the risk difference scale.",
                          call. = FALSE)
    sides <- 1L
  }

  delta <- p1 - p2
  se1   <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)  # SE under H1
  za    <- z_alpha(alpha, sides)

  switch(
    hypothesis,
    superiority = {
      se0 <- fm_null_se(p1, p2, n1, n2, 0)  # pooled SE under H0
      p <- stats::pnorm((delta - za * se0) / se1)
      if (sides == 2L) p <- p + stats::pnorm((-delta - za * se0) / se1)
      p
    },
    noninferiority = {
      se0 <- fm_null_se(p1, p2, n1, n2, -margin)
      stats::pnorm(((delta + margin) - za * se0) / se1)
    },
    superiority_margin = {
      se0 <- fm_null_se(p1, p2, n1, n2, margin)
      stats::pnorm(((delta - margin) - za * se0) / se1)
    },
    equivalence = {
      se0_lo <- fm_null_se(p1, p2, n1, n2, -margin)
      se0_hi <- fm_null_se(p1, p2, n1, n2,  margin)
      max(0, stats::pnorm(((margin - delta) - za * se0_hi) / se1) +
             stats::pnorm(((margin + delta) - za * se0_lo) / se1) - 1)
    }
  )
}

# ---------------------------------------------------------------------------
# Sample size
# ---------------------------------------------------------------------------

#' Sample size for the difference of two independent proportions
#'
#' Returns the smallest per-group sample size \eqn{n_1} (with
#' \eqn{n_2 = \lceil r\, n_1 \rceil}) whose power reaches \code{power},
#' on the risk difference scale. Superiority uses the pooled z-test
#' (uncorrected chi-square); non-inferiority, superiority by a margin, and
#' equivalence use the Farrington-Manning score test.
#'
#' Group 1 is the treatment group and higher proportions are assumed
#' favourable. If lower is better for your endpoint, use the complementary
#' event probabilities \code{1 - p} for both groups.
#'
#' @param p1 True (anticipated) proportion in group 1 (treatment).
#' @param p2 True (anticipated) proportion in group 2 (control).
#' @param power Target power (e.g., 0.80 or 0.90).
#' @param alpha Significance level. One-sided for margin-based hypotheses
#'   (e.g., 0.025); two-sided for superiority (e.g., 0.05).
#' @param ratio Allocation ratio \eqn{r = n_2 / n_1}.
#' @param hypothesis One of \code{"superiority"}, \code{"noninferiority"},
#'   \code{"superiority_margin"}, or \code{"equivalence"}.
#' @param margin Positive margin on the risk difference scale, required for
#'   all hypotheses except \code{"superiority"}.
#' @param sides 1 or 2. Relevant for \code{"superiority"} only.
#' @param dropout Anticipated dropout proportion in \eqn{[0, 1)}.
#'
#' @return A \code{\link{clinss_result}} object.
#'
#' @references
#' Farrington, C.P., Manning, G. (1990). Test statistics and sample size
#' formulae for comparative binomial trials with null hypothesis of
#' non-zero risk difference or non-unity relative risk.
#' \emph{Statistics in Medicine}, 9(12), 1447--1454.
#'
#' Fleiss, J.L., Levin, B., Paik, M.C. (2003). \emph{Statistical Methods
#' for Rates and Proportions}, 3rd ed. Wiley. Chapter 4.
#'
#' Chow, S.C., Shao, J., Wang, H. (2008). \emph{Sample Size Calculations
#' in Clinical Research}, 2nd ed. Chapman & Hall/CRC. Chapter 4.
#'
#' @seealso \code{\link{power_two_proportions}},
#'   \code{\link{ss_one_proportion}}
#'
#' @examples
#' # Superiority: 60% vs 40% response, 90% power, two-sided 0.05
#' ss_two_proportions(p1 = 0.60, p2 = 0.40, power = 0.90, alpha = 0.05,
#'                    sides = 2)
#'
#' # Non-inferiority: both arms 85%, margin 10 percentage points
#' res <- ss_two_proportions(p1 = 0.85, p2 = 0.85, margin = 0.10,
#'                           alpha = 0.025, power = 0.90,
#'                           hypothesis = "noninferiority")
#' report(res)
#'
#' @export
ss_two_proportions <- function(p1, p2,
                               power = 0.90,
                               alpha = 0.025,
                               ratio = 1,
                               hypothesis = c("superiority", "noninferiority",
                                              "superiority_margin",
                                              "equivalence"),
                               margin = NULL,
                               sides = if (hypothesis == "superiority") 2L else 1L,
                               dropout = 0) {
  hypothesis <- match.arg(hypothesis)
  check_prob(p1, "p1"); check_prob(p2, "p2")
  check_prob(power, "power"); check_prob(alpha, "alpha")
  check_positive(ratio, "ratio")
  if (dropout < 0 || dropout >= 1) {
    stop("`dropout` must be in [0, 1).", call. = FALSE)
  }

  delta <- p1 - p2
  dist <- switch(
    hypothesis,
    superiority = {
      if (delta == 0) stop("`p1` and `p2` must differ for a superiority test.",
                           call. = FALSE)
      abs(delta)
    },
    noninferiority = {
      check_positive(margin, "margin"); sides <- 1L
      d <- delta + margin
      if (d <= 0) stop(
        "`p1` - `p2` + `margin` must be positive: the true difference must ",
        "lie on the non-inferior side of the margin.", call. = FALSE)
      d
    },
    superiority_margin = {
      check_positive(margin, "margin"); sides <- 1L
      d <- delta - margin
      if (d <= 0) stop(
        "`p1` - `p2` must exceed `margin` for a superiority-by-a-margin test.",
        call. = FALSE)
      d
    },
    equivalence = {
      check_positive(margin, "margin"); sides <- 1L
      d <- margin - abs(delta)
      if (d <= 0) stop(
        "|`p1` - `p2`| must be smaller than `margin` for an equivalence test.",
        call. = FALSE)
      d
    }
  )

  # Normal-approximation seed using unpooled variance (adequate as a
  # starting point; the walk below uses the exact engine).
  zb  <- if (hypothesis == "equivalence" && delta == 0) {
    stats::qnorm(1 - (1 - power) / 2)
  } else {
    z_beta(power)
  }
  vbar <- p1 * (1 - p1) + p2 * (1 - p2) / ratio
  n1 <- vbar * (z_alpha(alpha, sides) + zb)^2 / dist^2
  n1 <- max(2L, ceiling(n1))

  pw_at <- function(n1) {
    power_two_proportions(n1 = n1, n2 = ceiling(ratio * n1), p1 = p1, p2 = p2,
                          alpha = alpha, hypothesis = hypothesis,
                          margin = margin, sides = sides)
  }

  while (pw_at(n1) < power) n1 <- n1 + 1L
  while (n1 > 2L && pw_at(n1 - 1L) >= power) n1 <- n1 - 1L

  n2    <- ceiling(ratio * n1)
  n_vec <- c(`group 1` = as.integer(n1), `group 2` = as.integer(n2))

  method <- if (hypothesis == "superiority") {
    "Two-Sample Z-Test for the Difference of Two Proportions (Pooled)"
  } else {
    "Farrington-Manning Score Test for the Difference of Two Proportions"
  }

  new_clinss_result(
    n = n_vec, n_total = n1 + n2,
    power_target = power, power_achieved = pw_at(n1),
    alpha = alpha, sides = sides, hypothesis = hypothesis,
    method = method, reference = REF_TWO_PROPS,
    parameters = list(p1 = p1, p2 = p2, margin = margin, ratio = ratio),
    dropout = dropout,
    n_enrolled = if (dropout > 0) adjust_dropout(n1 + n2, dropout) else NULL,
    endpoint = "two_proportions"
  )
}
