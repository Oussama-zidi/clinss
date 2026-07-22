# continuous-one-mean.R
#
# Sample size and power for:
#   - a single mean against a reference value (one-sample t / z test)
#   - the mean of paired differences (paired t / z test)
#
# Both reduce to the same one-sample t-test calculation.  They are exposed
# as two separate functions because the user-facing framing differs: in
# ss_one_mean() `sd` is the population SD; in ss_paired_means() `sd_diff`
# is the SD of within-pair differences, which the user may derive from the
# endpoint SD and the within-pair correlation via
# sd_diff = sd * sqrt(2 * (1 - rho)).
#
# Formula source (primary literature):
#   Julious, S.A. (2004). Statistics in Medicine, 23(12), 1921-1986.
#   Chow, S.C., Shao, J., Wang, H. (2008). Chapter 3, Sections 3.1, 3.3.

#' @keywords internal
REF_ONE_MEAN <- paste(
  "Julious, S.A. (2004). Sample sizes for clinical trials with Normal data.",
  "Statistics in Medicine, 23(12), 1921-1986;",
  "Chow, S.C., Shao, J., Wang, H. (2008). Sample Size Calculations in",
  "Clinical Research, 2nd ed. Chapman & Hall/CRC, Chapter 3.")

# ---------------------------------------------------------------------------
# Shared power engine (not exported; users go through ss_one_mean etc.)
# ---------------------------------------------------------------------------

#' @keywords internal
power_one_mean <- function(n, effect, sd, alpha = 0.025,
                           hypothesis = "superiority", margin = NULL,
                           sides = if (hypothesis == "superiority") 2L else 1L,
                           test = "t") {
  se <- sd / sqrt(n)
  if (test == "z") {
    za <- z_alpha(alpha, sides)
    switch(hypothesis,
      superiority = {
        p <- stats::pnorm(effect / se - za)
        if (sides == 2L) p <- p + stats::pnorm(-effect / se - za)
        p
      },
      noninferiority     = stats::pnorm((effect + margin) / se - za),
      superiority_margin = stats::pnorm((effect - margin) / se - za),
      equivalence = max(0, stats::pnorm((margin - effect) / se - za) +
                           stats::pnorm((margin + effect) / se - za) - 1))
  } else {
    df    <- n - 1
    tcrit <- stats::qt(1 - alpha / sides, df)
    switch(hypothesis,
      superiority = {
        ncp <- effect / se
        p <- stats::pt(tcrit, df, ncp = ncp, lower.tail = FALSE)
        if (sides == 2L) p <- p + stats::pt(-tcrit, df, ncp = ncp)
        p
      },
      noninferiority = stats::pt(tcrit, df, ncp = (effect + margin) / se,
                                 lower.tail = FALSE),
      superiority_margin = stats::pt(tcrit, df, ncp = (effect - margin) / se,
                                     lower.tail = FALSE),
      equivalence = max(0,
        stats::pt(-tcrit, df, ncp = (effect - margin) / se) -
        stats::pt( tcrit, df, ncp = (effect + margin) / se)))
  }
}

# ---------------------------------------------------------------------------
# Shared sample size engine (used by both ss_one_mean and ss_paired_means)
# ---------------------------------------------------------------------------

#' @keywords internal
ss_one_mean_engine <- function(effect, sd, power, alpha, hypothesis, margin,
                               sides, test, dropout, method, reference,
                               group_label) {
  check_number(effect, "effect"); check_positive(sd, "sd")
  check_prob(power, "power"); check_prob(alpha, "alpha")
  if (dropout < 0 || dropout >= 1) {
    stop("`dropout` must be in [0, 1).", call. = FALSE)
  }

  delta <- switch(
    hypothesis,
    superiority = {
      if (effect == 0) stop("`effect` must be non-zero for a superiority test.",
                            call. = FALSE)
      abs(effect)
    },
    noninferiority = {
      check_positive(margin, "margin"); sides <- 1L
      d <- effect + margin
      if (d <= 0) stop("`effect` + `margin` must be positive.", call. = FALSE)
      d
    },
    superiority_margin = {
      check_positive(margin, "margin"); sides <- 1L
      d <- effect - margin
      if (d <= 0) stop("`effect` must exceed `margin`.", call. = FALSE)
      d
    },
    equivalence = {
      check_positive(margin, "margin"); sides <- 1L
      d <- margin - abs(effect)
      if (d <= 0) stop("|`effect`| must be smaller than `margin`.",
                       call. = FALSE)
      d
    }
  )

  zb <- if (hypothesis == "equivalence" && effect == 0) {
    stats::qnorm(1 - (1 - power) / 2)
  } else {
    z_beta(power)
  }
  n <- (sd / delta)^2 * (z_alpha(alpha, sides) + zb)^2
  n <- max(2L, ceiling(n))

  pw_at <- function(n) power_one_mean(n, effect, sd, alpha, hypothesis,
                                      margin, sides, test)
  while (pw_at(n) < power) n <- n + 1L
  while (n > 2L && pw_at(n - 1L) >= power) n <- n - 1L

  nn <- c(as.integer(n)); names(nn) <- group_label
  new_clinss_result(
    n = nn, n_total = as.integer(n),
    power_target = power, power_achieved = pw_at(n),
    alpha = alpha, sides = sides, hypothesis = hypothesis,
    method = method, reference = reference,
    parameters = list(effect = effect, sd = sd, margin = margin, test = test),
    dropout = dropout,
    n_enrolled = if (dropout > 0) adjust_dropout(n, dropout) else NULL
  )
}

# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

#' Sample size for one mean against a reference value
#'
#' Returns the smallest sample size whose exact power reaches \code{power}
#' for a one-sample test of a mean against a fixed reference value, under
#' superiority, non-inferiority, superiority by a margin, or equivalence
#' (TOST) hypotheses.
#'
#' @param effect Difference between the true mean and the reference value.
#' @param sd Standard deviation of the endpoint.
#' @param power Target power.
#' @param alpha Significance level (one-sided for margin-based hypotheses).
#' @param hypothesis One of \code{"superiority"}, \code{"noninferiority"},
#'   \code{"superiority_margin"}, or \code{"equivalence"}.
#' @param margin Positive margin; required for margin-based hypotheses.
#' @param sides 1 or 2. Relevant for \code{"superiority"} only.
#' @param test \code{"t"} (default, exact noncentral t) or \code{"z"}
#'   (normal approximation).
#' @param dropout Anticipated dropout proportion in \eqn{[0, 1)}.
#'
#' @return A \code{\link{clinss_result}} object.
#'
#' @references
#' Julious, S.A. (2004). Sample sizes for clinical trials with Normal
#' data. \emph{Statistics in Medicine}, 23(12), 1921--1986.
#'
#' Chow, S.C., Shao, J., Wang, H. (2008). \emph{Sample Size Calculations
#' in Clinical Research}, 2nd ed. Chapman & Hall/CRC. Chapter 3.
#'
#' @seealso \code{\link{ss_two_means}}, \code{\link{ss_paired_means}}
#'
#' @examples
#' ss_one_mean(effect = 5, sd = 10, power = 0.90, alpha = 0.05, sides = 2)
#'
#' @export
ss_one_mean <- function(effect, sd,
                        power = 0.90, alpha = 0.025,
                        hypothesis = c("superiority", "noninferiority",
                                       "superiority_margin", "equivalence"),
                        margin = NULL,
                        sides = if (hypothesis == "superiority") 2L else 1L,
                        test = c("t", "z"), dropout = 0) {
  hypothesis <- match.arg(hypothesis)
  test <- match.arg(test)
  ss_one_mean_engine(effect, sd, power, alpha, hypothesis, margin, sides,
                     test, dropout,
                     method    = sprintf("One-Sample %s-Test for One Mean",
                                         toupper(test)),
                     reference = REF_ONE_MEAN,
                     group_label = "subjects")
}

#' Sample size for paired means
#'
#' Returns the smallest number of pairs whose exact power reaches
#' \code{power} for a paired test (one-sample t or z test on within-pair
#' differences), under superiority, non-inferiority, superiority by a
#' margin, or equivalence hypotheses.
#'
#' @param effect True mean of the within-pair differences.
#' @param sd_diff Standard deviation of the within-pair differences.
#'   If you know the endpoint SD \eqn{\sigma} and the within-pair
#'   correlation \eqn{\rho}, use
#'   \eqn{\sigma_d = \sigma \sqrt{2(1 - \rho)}}.
#' @param power Target power.
#' @param alpha Significance level (one-sided for margin-based hypotheses).
#' @param hypothesis One of \code{"superiority"}, \code{"noninferiority"},
#'   \code{"superiority_margin"}, or \code{"equivalence"}.
#' @param margin Positive margin; required for margin-based hypotheses.
#' @param sides 1 or 2. Relevant for \code{"superiority"} only.
#' @param test \code{"t"} (default) or \code{"z"}.
#' @param dropout Anticipated dropout proportion in \eqn{[0, 1)}.
#'
#' @return A \code{\link{clinss_result}} object. The \code{n} component
#'   gives the number of evaluable pairs; \code{n_total} equals \code{n}
#'   because there is a single sampling unit (the pair).
#'
#' @references
#' Julious, S.A. (2004). Sample sizes for clinical trials with Normal
#' data. \emph{Statistics in Medicine}, 23(12), 1921--1986.
#'
#' Chow, S.C., Shao, J., Wang, H. (2008). \emph{Sample Size Calculations
#' in Clinical Research}, 2nd ed. Chapman & Hall/CRC. Chapter 3.
#'
#' @seealso \code{\link{ss_two_means}}, \code{\link{ss_one_mean}}
#'
#' @examples
#' ss_paired_means(effect = 2, sd_diff = 6, power = 0.80, alpha = 0.05,
#'                 sides = 2)
#'
#' # Derive sd_diff from endpoint SD and within-pair correlation
#' sd_endpoint <- 10; rho <- 0.60
#' sd_diff <- sd_endpoint * sqrt(2 * (1 - rho))
#' ss_paired_means(effect = 3, sd_diff = sd_diff, power = 0.90, alpha = 0.05,
#'                 sides = 2)
#'
#' @export
ss_paired_means <- function(effect, sd_diff,
                            power = 0.90, alpha = 0.025,
                            hypothesis = c("superiority", "noninferiority",
                                           "superiority_margin", "equivalence"),
                            margin = NULL,
                            sides = if (hypothesis == "superiority") 2L else 1L,
                            test = c("t", "z"), dropout = 0) {
  hypothesis <- match.arg(hypothesis)
  test <- match.arg(test)
  ss_one_mean_engine(effect, sd_diff, power, alpha, hypothesis, margin, sides,
                     test, dropout,
                     method    = sprintf("Paired %s-Test for the Mean Difference",
                                         toupper(test)),
                     reference = REF_ONE_MEAN,
                     group_label = "pairs")
}
