# continuous-two-means.R
#
# Sample size and power for the comparison of two independent means
# (parallel-group design, continuous endpoint).
#
# Hypotheses supported:
#   superiority         H0: mu1 - mu2  = 0        vs H1: mu1 - mu2 != 0 (or >)
#   noninferiority      H0: mu1 - mu2 <= -margin  vs H1: mu1 - mu2  > -margin
#   superiority_margin  H0: mu1 - mu2 <=  margin  vs H1: mu1 - mu2  >  margin
#   equivalence (TOST)  H0: |mu1 - mu2| >= margin vs H1: |mu1 - mu2| < margin
#
# Formula source (primary literature, not software documentation):
#   Julious, S.A. (2004). Sample sizes for clinical trials with Normal data.
#     Statistics in Medicine, 23(12), 1921-1986.
#     Section 2: superiority; Section 4: non-inferiority; Section 6: TOST.
#   Chow, S.C., Shao, J., Wang, H. (2008). Sample Size Calculations in
#     Clinical Research, 2nd ed. Chapman & Hall/CRC. Chapter 3.
#
# Conventions:
#   effect = mu1 - mu2  (positive means group 1 is higher)
#   ratio  = n2 / n1
#   margin is always a positive number

#' @keywords internal
REF_TWO_MEANS <- paste(
  "Julious, S.A. (2004). Sample sizes for clinical trials with Normal data.",
  "Statistics in Medicine, 23(12), 1921-1986;",
  "Chow, S.C., Shao, J., Wang, H. (2008). Sample Size Calculations in",
  "Clinical Research, 2nd ed. Chapman & Hall/CRC, Chapter 3.")

# ---------------------------------------------------------------------------
# Power engine
# ---------------------------------------------------------------------------

#' Power for the difference of two independent means
#'
#' Computes the power of a two-sample test comparing two independent means
#' at a fixed sample size. Use \code{\link{ss_two_means}} to find the sample
#' size that achieves a target power.
#'
#' @param n1 Sample size in group 1.
#' @param n2 Sample size in group 2. Defaults to \code{n1} (equal allocation).
#' @param effect True difference in means, \eqn{\mu_1 - \mu_2}.
#' @param sd Common standard deviation of the endpoint.
#' @param alpha Significance level (one-sided for margin-based hypotheses).
#' @param hypothesis One of \code{"superiority"}, \code{"noninferiority"},
#'   \code{"superiority_margin"}, or \code{"equivalence"}.
#' @param margin Positive margin, required for all hypotheses except
#'   \code{"superiority"}.
#' @param sides 1 or 2. Used only for \code{hypothesis = "superiority"}.
#' @param test \code{"t"} for the noncentral t distribution (default) or
#'   \code{"z"} for the normal approximation.
#'
#' @return A single numeric value: the power at the specified sample sizes.
#'
#' @references
#' Julious, S.A. (2004). Sample sizes for clinical trials with Normal
#' data. \emph{Statistics in Medicine}, 23(12), 1921--1986.
#'
#' Chow, S.C., Shao, J., Wang, H. (2008). \emph{Sample Size Calculations
#' in Clinical Research}, 2nd ed. Chapman & Hall/CRC. Chapter 3.
#'
#' @examples
#' # Power at n = 86 per group (should be just over 0.90)
#' power_two_means(n1 = 86, effect = 5, sd = 10, alpha = 0.05, sides = 2)
#'
#' # Power curve: how does power grow with n?
#' n_seq <- seq(20, 150, by = 10)
#' pwr   <- sapply(n_seq, function(n)
#'   power_two_means(n1 = n, effect = 5, sd = 10, alpha = 0.05, sides = 2))
#' plot(n_seq, pwr, type = "l", xlab = "n per group", ylab = "Power")
#' abline(h = 0.90, lty = 2)
#'
#' @export
power_two_means <- function(n1, n2 = n1, effect, sd,
                            alpha = 0.025,
                            hypothesis = c("superiority", "noninferiority",
                                           "superiority_margin", "equivalence"),
                            margin = NULL,
                            sides = if (hypothesis == "superiority") 2L else 1L,
                            test = c("t", "z")) {
  hypothesis <- match.arg(hypothesis)
  test <- match.arg(test)
  check_positive(n1, "n1"); check_positive(n2, "n2")
  check_number(effect, "effect"); check_positive(sd, "sd")
  check_prob(alpha, "alpha")
  if (hypothesis != "superiority") {
    check_positive(margin, "margin")
    sides <- 1L
  }

  se <- sd * sqrt(1 / n1 + 1 / n2)

  if (test == "z") {
    za <- z_alpha(alpha, sides)
    pw <- switch(
      hypothesis,
      superiority = {
        p <- stats::pnorm(effect / se - za)
        if (sides == 2L) p <- p + stats::pnorm(-effect / se - za)
        p
      },
      noninferiority     = stats::pnorm((effect + margin) / se - za),
      superiority_margin = stats::pnorm((effect - margin) / se - za),
      equivalence = max(0,
        stats::pnorm((margin - effect) / se - za) +
        stats::pnorm((margin + effect) / se - za) - 1)
    )
  } else {
    df    <- n1 + n2 - 2
    tcrit <- stats::qt(1 - alpha / sides, df)
    pw <- switch(
      hypothesis,
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
        stats::pt( tcrit, df, ncp = (effect + margin) / se))
    )
  }
  pw
}

# ---------------------------------------------------------------------------
# Sample size
# ---------------------------------------------------------------------------

#' Sample size for the difference of two independent means
#'
#' Returns the smallest per-group sample size \eqn{n_1} (with
#' \eqn{n_2 = \lceil r\, n_1 \rceil}) whose exact power reaches
#' \code{power}, for superiority, non-inferiority, superiority by a margin,
#' or equivalence (TOST) hypotheses.
#'
#' The search is seeded by the closed-form normal approximation
#' \deqn{n_1 = (1 + 1/r)\,\sigma^2 (z_{1-\alpha} + z_{1-\beta})^2 / \delta^2}
#' where \eqn{\delta} is the distance between the true effect and the null
#' boundary: \eqn{|\text{effect}|} for superiority,
#' \eqn{\text{effect} + m} for non-inferiority,
#' \eqn{\text{effect} - m} for superiority by a margin, and
#' \eqn{m - |\text{effect}|} for equivalence.
#' The function then walks to the exact integer minimum using the noncentral
#' t distribution (when \code{test = "t"}).
#'
#' @param effect True difference in means, \eqn{\mu_1 - \mu_2}. Positive
#'   values favour group 1.
#' @param sd Common standard deviation of the endpoint.
#' @param power Target power (e.g., 0.80 or 0.90).
#' @param alpha Significance level. For one-sided hypotheses (all
#'   margin-based types) this is the one-sided alpha; for two-sided
#'   superiority it is the two-sided alpha.
#' @param ratio Allocation ratio \eqn{r = n_2 / n_1}. Default 1 gives
#'   equal allocation.
#' @param hypothesis One of \code{"superiority"} (test of no difference),
#'   \code{"noninferiority"}, \code{"superiority_margin"}, or
#'   \code{"equivalence"} (TOST).
#' @param margin Positive non-inferiority or equivalence margin \eqn{m}.
#'   Required for all hypotheses except \code{"superiority"}.
#' @param sides 1 or 2. Relevant for \code{"superiority"} only; all
#'   margin-based hypotheses are inherently one-sided.
#' @param test \code{"t"} for the exact noncentral t calculation (default)
#'   or \code{"z"} for the normal approximation.
#' @param dropout Anticipated dropout proportion in \eqn{[0, 1)}. The
#'   evaluable sample size is inflated to
#'   \code{ceiling(n / (1 - dropout))} and reported as \code{n_enrolled}.
#'
#' @return An object of class \code{\link{clinss_result}} with components:
#' \describe{
#'   \item{n}{Named integer vector of per-group evaluable sample sizes.}
#'   \item{n_total}{Total evaluable sample size.}
#'   \item{n_enrolled}{Total enrolled size after dropout inflation, or
#'     \code{NULL} if \code{dropout = 0}.}
#'   \item{power_achieved}{Exact power at the returned sample sizes.}
#'   \item{power_target}{The requested \code{power}.}
#'   \item{alpha, sides, hypothesis, method, reference}{Design metadata.}
#'   \item{parameters}{Named list of design assumptions.}
#' }
#'
#' @references
#' Julious, S.A. (2004). Sample sizes for clinical trials with Normal
#' data. \emph{Statistics in Medicine}, 23(12), 1921--1986.
#'
#' Chow, S.C., Shao, J., Wang, H. (2008). \emph{Sample Size Calculations
#' in Clinical Research}, 2nd ed. Chapman & Hall/CRC. Chapter 3.
#'
#' @seealso \code{\link{power_two_means}}, \code{\link{ss_one_mean}},
#'   \code{\link{ss_paired_means}}, \code{\link{report}}
#'
#' @examples
#' # Superiority: detect a 5-unit difference (SD 10), 90% power, alpha 0.05
#' ss_two_means(effect = 5, sd = 10, power = 0.90, alpha = 0.05, sides = 2)
#'
#' # Non-inferiority: margin 10, assuming truly equal means
#' ss_two_means(effect = 0, sd = 20, margin = 10, alpha = 0.025,
#'              power = 0.90, hypothesis = "noninferiority")
#'
#' # Equivalence (TOST): 2:1 allocation, 10% dropout
#' res <- ss_two_means(effect = 0, sd = 10, margin = 5, alpha = 0.05,
#'                     power = 0.80, hypothesis = "equivalence",
#'                     ratio = 2, dropout = 0.10)
#' summary(res)
#' report(res)
#'
#' # Sensitivity table over a range of standard deviations
#' do.call(rbind, lapply(c(8, 10, 12), function(s)
#'   as.data.frame(ss_two_means(effect = 5, sd = s, power = 0.90,
#'                              alpha = 0.05, sides = 2))))
#'
#' @export
ss_two_means <- function(effect, sd,
                         power = 0.90,
                         alpha = 0.025,
                         ratio = 1,
                         hypothesis = c("superiority", "noninferiority",
                                        "superiority_margin", "equivalence"),
                         margin = NULL,
                         sides = if (hypothesis == "superiority") 2L else 1L,
                         test = c("t", "z"),
                         dropout = 0) {
  hypothesis <- match.arg(hypothesis)
  test <- match.arg(test)
  check_number(effect, "effect"); check_positive(sd, "sd")
  check_prob(power, "power"); check_prob(alpha, "alpha")
  check_positive(ratio, "ratio")
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
      if (d <= 0) stop(
        "`effect` + `margin` must be positive: the true effect must lie on ",
        "the non-inferior side of the margin.", call. = FALSE)
      d
    },
    superiority_margin = {
      check_positive(margin, "margin"); sides <- 1L
      d <- effect - margin
      if (d <= 0) stop(
        "`effect` must exceed `margin` for a superiority-by-a-margin test.",
        call. = FALSE)
      d
    },
    equivalence = {
      check_positive(margin, "margin"); sides <- 1L
      d <- margin - abs(effect)
      if (d <= 0) stop(
        "|`effect`| must be smaller than `margin` for an equivalence test.",
        call. = FALSE)
      d
    }
  )

  # Closed-form normal-approximation seed (Julious 2004, eq. 5 / Chow 3.2.3).
  # For TOST with effect = 0 the z_{1-beta} term uses z_{1-beta/2}.
  zb <- if (hypothesis == "equivalence" && effect == 0) {
    stats::qnorm(1 - (1 - power) / 2)
  } else {
    z_beta(power)
  }
  n1 <- allocation_factor(ratio) * (sd / delta)^2 *
        (z_alpha(alpha, sides) + zb)^2
  n1 <- max(2L, ceiling(n1))

  pw_at <- function(n1) {
    power_two_means(n1 = n1, n2 = ceiling(ratio * n1), effect = effect,
                    sd = sd, alpha = alpha, hypothesis = hypothesis,
                    margin = margin, sides = sides, test = test)
  }

  while (pw_at(n1) < power) n1 <- n1 + 1L
  while (n1 > 2L && pw_at(n1 - 1L) >= power) n1 <- n1 - 1L

  n2    <- ceiling(ratio * n1)
  n_vec <- c(`group 1` = as.integer(n1), `group 2` = as.integer(n2))

  method <- sprintf("Two-Sample %s-Test for the Difference of Two Means",
                    toupper(test))

  new_clinss_result(
    n = n_vec, n_total = n1 + n2,
    power_target = power, power_achieved = pw_at(n1),
    alpha = alpha, sides = sides, hypothesis = hypothesis,
    method = method, reference = REF_TWO_MEANS,
    parameters = list(effect = effect, sd = sd, margin = margin,
                      ratio = ratio, test = test),
    dropout = dropout,
    n_enrolled = if (dropout > 0) adjust_dropout(n1 + n2, dropout) else NULL
  )
}
