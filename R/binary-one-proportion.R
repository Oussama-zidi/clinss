# binary-one-proportion.R
#
# Sample size and power for a single proportion against a reference value,
# using either the exact binomial test or the z-test (score statistic, with
# the variance evaluated at the null boundary value).
#
# Hypotheses supported (p = true proportion, p0 = reference):
#   superiority         H0: p  = p0        vs H1: p != p0 (or >)
#   noninferiority      H0: p <= p0 - m    vs H1: p  > p0 - m
#   superiority_margin  H0: p <= p0 + m    vs H1: p  > p0 + m
#   equivalence (TOST)  H0: |p - p0| >= m  vs H1: |p - p0| < m
#
# Formula source (primary literature):
#   Chow, S.C., Shao, J., Wang, H. (2008). Sample Size Calculations in
#     Clinical Research, 2nd ed. Chapman & Hall/CRC. Chapter 4 (one-sample
#     z-test formulations).
#   Fleiss, J.L., Levin, B., Paik, M.C. (2003). Statistical Methods for
#     Rates and Proportions, 3rd ed. Wiley. Chapter 2 (exact binomial test).
#
# Conventions:
#   Higher proportions are assumed favourable for margin-based hypotheses.
#   If lower is better, work with the complementary event (p <- 1 - p,
#   p0 <- 1 - p0). Margins are always positive.
#
# Note on the exact test:
#   Exact binomial power is a step function that is NOT monotone in n
#   (the "sawtooth" phenomenon). ss_one_proportion() therefore scans n
#   upward and returns the SMALLEST n whose exact power reaches the target;
#   power may dip below the target again for some larger n. This is
#   standard behaviour for exact tests and is documented in the help page.

#' @keywords internal
REF_ONE_PROP <- paste(
  "Chow, S.C., Shao, J., Wang, H. (2008). Sample Size Calculations in",
  "Clinical Research, 2nd ed. Chapman & Hall/CRC, Chapter 4;",
  "Fleiss, J.L., Levin, B., Paik, M.C. (2003). Statistical Methods for",
  "Rates and Proportions, 3rd ed. Wiley, Chapter 2.")

# ---------------------------------------------------------------------------
# Exact one-sided binomial rejection regions
# ---------------------------------------------------------------------------

# Smallest k such that P(X >= k | n, pb) <= alpha  (upper-tailed test).
#' @keywords internal
exact_crit_upper <- function(n, pb, alpha) {
  k <- stats::qbinom(1 - alpha, n, pb) + 1
  # qbinom guarantees P(X <= k-1) >= 1 - alpha, so P(X >= k) <= alpha.
  k
}

# Largest k such that P(X <= k | n, pb) <= alpha  (lower-tailed test).
#' @keywords internal
exact_crit_lower <- function(n, pb, alpha) {
  k <- stats::qbinom(alpha, n, pb) - 1
  # Guard: qbinom(alpha,...) may itself have cdf > alpha; step down if needed.
  while (k >= 0 && stats::pbinom(k, n, pb) > alpha) k <- k - 1
  k
}

# ---------------------------------------------------------------------------
# Power engine
# ---------------------------------------------------------------------------

#' Power for one proportion against a reference value
#'
#' Computes the power of a one-sample test of a proportion at a fixed
#' sample size, using either the exact binomial test (default) or the
#' z-test with the variance evaluated at the null boundary value.
#'
#' Higher proportions are assumed favourable for margin-based hypotheses.
#' If lower is better for your endpoint, use the complementary event:
#' \code{p <- 1 - p}, \code{p0 <- 1 - p0}.
#'
#' @param n Sample size.
#' @param p True (anticipated) proportion.
#' @param p0 Reference (null) proportion.
#' @param alpha Significance level (one-sided for margin-based hypotheses).
#' @param hypothesis One of \code{"superiority"}, \code{"noninferiority"},
#'   \code{"superiority_margin"}, or \code{"equivalence"}.
#' @param margin Positive margin, required for margin-based hypotheses.
#' @param sides 1 or 2. Used only for \code{hypothesis = "superiority"};
#'   the exact two-sided test is performed as two one-sided tests at
#'   \code{alpha/2}.
#' @param test \code{"exact"} for the exact binomial test (default) or
#'   \code{"z"} for the normal approximation.
#'
#' @return A single numeric value: the power at sample size \code{n}.
#'
#' @references
#' Chow, S.C., Shao, J., Wang, H. (2008). \emph{Sample Size Calculations
#' in Clinical Research}, 2nd ed. Chapman & Hall/CRC. Chapter 4.
#'
#' Fleiss, J.L., Levin, B., Paik, M.C. (2003). \emph{Statistical Methods
#' for Rates and Proportions}, 3rd ed. Wiley. Chapter 2.
#'
#' @examples
#' # Exact power at n = 50 for 35% response vs a 20% reference
#' power_one_proportion(n = 50, p = 0.35, p0 = 0.20, alpha = 0.05, sides = 1)
#'
#' @export
power_one_proportion <- function(n, p, p0,
                                 alpha = 0.025,
                                 hypothesis = c("superiority",
                                                "noninferiority",
                                                "superiority_margin",
                                                "equivalence"),
                                 margin = NULL,
                                 sides = if (hypothesis == "superiority") 2L else 1L,
                                 test = c("exact", "z")) {
  hypothesis <- match.arg(hypothesis)
  test <- match.arg(test)
  check_positive(n, "n"); check_prob(p, "p"); check_prob(p0, "p0")
  check_prob(alpha, "alpha")
  if (hypothesis != "superiority") {
    check_positive(margin, "margin")
    sides <- 1L
  }

  # Null boundary value(s) at which the test is carried out.
  pb <- switch(hypothesis,
               superiority        = p0,
               noninferiority     = p0 - margin,
               superiority_margin = p0 + margin,
               equivalence        = c(lo = p0 - margin, hi = p0 + margin))
  if (any(pb <= 0 | pb >= 1)) {
    stop("The null boundary `p0` +/- `margin` must lie strictly in (0, 1).",
         call. = FALSE)
  }

  if (test == "z") {
    za  <- z_alpha(alpha, sides)
    se1 <- sqrt(p * (1 - p) / n)
    switch(
      hypothesis,
      superiority = {
        se0 <- sqrt(p0 * (1 - p0) / n)
        pw  <- stats::pnorm(((p - p0) - za * se0) / se1)
        if (sides == 2L) pw <- pw + stats::pnorm((-(p - p0) - za * se0) / se1)
        pw
      },
      noninferiority = {
        se0 <- sqrt(pb * (1 - pb) / n)
        stats::pnorm(((p - pb) - za * se0) / se1)
      },
      superiority_margin = {
        se0 <- sqrt(pb * (1 - pb) / n)
        stats::pnorm(((p - pb) - za * se0) / se1)
      },
      equivalence = {
        se_lo <- sqrt(pb["lo"] * (1 - pb["lo"]) / n)
        se_hi <- sqrt(pb["hi"] * (1 - pb["hi"]) / n)
        max(0, stats::pnorm(((pb["hi"] - p) - za * se_hi) / se1) +
               stats::pnorm(((p - pb["lo"]) - za * se_lo) / se1) - 1)
      }
    )
  } else {
    switch(
      hypothesis,
      superiority = {
        if (sides == 2L) {
          ku <- exact_crit_upper(n, p0, alpha / 2)
          kl <- exact_crit_lower(n, p0, alpha / 2)
          stats::pbinom(ku - 1, n, p, lower.tail = FALSE) +
            (if (kl >= 0) stats::pbinom(kl, n, p) else 0)
        } else if (p >= p0) {
          ku <- exact_crit_upper(n, p0, alpha)
          stats::pbinom(ku - 1, n, p, lower.tail = FALSE)
        } else {
          kl <- exact_crit_lower(n, p0, alpha)
          if (kl >= 0) stats::pbinom(kl, n, p) else 0
        }
      },
      noninferiority = {
        ku <- exact_crit_upper(n, pb, alpha)
        stats::pbinom(ku - 1, n, p, lower.tail = FALSE)
      },
      superiority_margin = {
        ku <- exact_crit_upper(n, pb, alpha)
        stats::pbinom(ku - 1, n, p, lower.tail = FALSE)
      },
      equivalence = {
        ku <- exact_crit_upper(n, pb["lo"], alpha)  # must exceed lower bound
        kl <- exact_crit_lower(n, pb["hi"], alpha)  # must fall below upper
        if (kl < ku) 0 else
          stats::pbinom(kl, n, p) - stats::pbinom(ku - 1, n, p)
      }
    )
  }
}

# ---------------------------------------------------------------------------
# Sample size
# ---------------------------------------------------------------------------

#' Sample size for one proportion against a reference value
#'
#' Returns the smallest sample size whose power reaches \code{power} for a
#' one-sample test of a proportion, using either the exact binomial test
#' (default) or the z-test approximation.
#'
#' Because exact binomial power is a step function that is not monotone in
#' \code{n} (the "sawtooth" phenomenon), the returned \code{n} is the
#' smallest sample size at which the target power is reached; power may dip
#' slightly below the target for some larger sample sizes. This is inherent
#' to exact tests. The achieved power and the exact rejection threshold are
#' both reported.
#'
#' Higher proportions are assumed favourable for margin-based hypotheses.
#' If lower is better for your endpoint, use the complementary event:
#' \code{p <- 1 - p}, \code{p0 <- 1 - p0}.
#'
#' @param p True (anticipated) proportion.
#' @param p0 Reference (null) proportion.
#' @param power Target power.
#' @param alpha Significance level (one-sided for margin-based hypotheses).
#' @param hypothesis One of \code{"superiority"}, \code{"noninferiority"},
#'   \code{"superiority_margin"}, or \code{"equivalence"}.
#' @param margin Positive margin, required for margin-based hypotheses.
#' @param sides 1 or 2. Relevant for \code{"superiority"} only.
#' @param test \code{"exact"} (default) or \code{"z"}.
#' @param dropout Anticipated dropout proportion in \eqn{[0, 1)}.
#'
#' @return A \code{\link{clinss_result}} object.
#'
#' @references
#' Chow, S.C., Shao, J., Wang, H. (2008). \emph{Sample Size Calculations
#' in Clinical Research}, 2nd ed. Chapman & Hall/CRC. Chapter 4.
#'
#' Fleiss, J.L., Levin, B., Paik, M.C. (2003). \emph{Statistical Methods
#' for Rates and Proportions}, 3rd ed. Wiley. Chapter 2.
#'
#' @seealso \code{\link{power_one_proportion}},
#'   \code{\link{ss_two_proportions}}
#'
#' @examples
#' # Exact test: 35% anticipated response vs 20% reference, one-sided 0.05
#' ss_one_proportion(p = 0.35, p0 = 0.20, power = 0.80, alpha = 0.05,
#'                   sides = 1)
#'
#' # z-approximation for comparison
#' ss_one_proportion(p = 0.35, p0 = 0.20, power = 0.80, alpha = 0.05,
#'                   sides = 1, test = "z")
#'
#' @export
ss_one_proportion <- function(p, p0,
                              power = 0.90, alpha = 0.025,
                              hypothesis = c("superiority", "noninferiority",
                                             "superiority_margin",
                                             "equivalence"),
                              margin = NULL,
                              sides = if (hypothesis == "superiority") 2L else 1L,
                              test = c("exact", "z"), dropout = 0) {
  hypothesis <- match.arg(hypothesis)
  test <- match.arg(test)
  check_prob(p, "p"); check_prob(p0, "p0")
  check_prob(power, "power"); check_prob(alpha, "alpha")
  if (dropout < 0 || dropout >= 1) {
    stop("`dropout` must be in [0, 1).", call. = FALSE)
  }

  dist <- switch(
    hypothesis,
    superiority = {
      if (p == p0) stop("`p` and `p0` must differ for a superiority test.",
                        call. = FALSE)
      abs(p - p0)
    },
    noninferiority = {
      check_positive(margin, "margin"); sides <- 1L
      d <- p - (p0 - margin)
      if (d <= 0) stop("`p` must exceed `p0` - `margin`.", call. = FALSE)
      d
    },
    superiority_margin = {
      check_positive(margin, "margin"); sides <- 1L
      d <- p - (p0 + margin)
      if (d <= 0) stop("`p` must exceed `p0` + `margin`.", call. = FALSE)
      d
    },
    equivalence = {
      check_positive(margin, "margin"); sides <- 1L
      d <- margin - abs(p - p0)
      if (d <= 0) stop("|`p` - `p0`| must be smaller than `margin`.",
                       call. = FALSE)
      d
    }
  )

  # Normal-approximation seed.
  zb <- if (hypothesis == "equivalence" && p == p0) {
    stats::qnorm(1 - (1 - power) / 2)
  } else {
    z_beta(power)
  }
  n_seed <- p * (1 - p) * (z_alpha(alpha, sides) + zb)^2 / dist^2
  n_seed <- max(5L, ceiling(n_seed))

  pw_at <- function(n) power_one_proportion(n, p, p0, alpha, hypothesis,
                                            margin, sides, test)

  if (test == "z") {
    n <- n_seed
    while (pw_at(n) < power) n <- n + 1L
    while (n > 2L && pw_at(n - 1L) >= power) n <- n - 1L
  } else {
    # Exact power is a non-monotone step function of n: scan upward from a
    # low start and take the FIRST n reaching the target.
    n_cap <- 4L * n_seed + 1000L
    n <- 2L
    found <- FALSE
    while (n <= n_cap) {
      if (pw_at(n) >= power) { found <- TRUE; break }
      n <- n + 1L
    }
    if (!found) stop("No sample size up to ", n_cap,
                     " reaches the target power; check the inputs.",
                     call. = FALSE)
  }

  nn <- c(subjects = as.integer(n))
  method <- if (test == "exact") "Exact Binomial Test for One Proportion"
            else "One-Sample Z-Test for One Proportion"

  new_clinss_result(
    n = nn, n_total = as.integer(n),
    power_target = power, power_achieved = pw_at(n),
    alpha = alpha, sides = sides, hypothesis = hypothesis,
    method = method, reference = REF_ONE_PROP,
    parameters = list(p = p, p0 = p0, margin = margin, test = test),
    dropout = dropout,
    n_enrolled = if (dropout > 0) adjust_dropout(n, dropout) else NULL,
    endpoint = "one_proportion"
  )
}
