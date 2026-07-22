# core-report.R
# report() turns a clinss_result into protocol-ready text suitable for a
# statistical analysis plan or study protocol. The wording is dispatched on
# the endpoint type stored in the result object, so each outcome module gets
# appropriate phrasing without duplicating the surrounding machinery.

#' Generate protocol-ready text from a sample size result
#'
#' Produces a sentence suitable for a study protocol or statistical analysis
#' plan, describing the sample size, its justification, all key assumptions,
#' and, when applicable, the enrolled size after dropout inflation.
#'
#' @param x A \code{clinss_result} object returned by any \code{ss_*()}
#'   function.
#' @param ... Ignored.
#'
#' @return A single character string.
#'
#' @examples
#' res <- ss_two_means(effect = 5, sd = 10, power = 0.90, alpha = 0.05,
#'                     sides = 2, dropout = 0.10)
#' cat(report(res))
#'
#' @export
report <- function(x, ...) UseMethod("report")

#' @rdname report
#' @export
report.clinss_result <- function(x, ...) {
  endpoint <- if (is.null(x$endpoint)) "continuous" else x$endpoint

  txt <- switch(endpoint,
                continuous      = report_continuous(x),
                two_proportions = report_two_proportions(x),
                one_proportion  = report_one_proportion(x),
                report_generic(x))

  if (!is.null(x$n_enrolled)) {
    txt <- paste0(
      txt,
      sprintf(paste0(" Allowing for %.0f%% dropout, %d participants will be ",
                     "enrolled."), 100 * x$dropout, x$n_enrolled))
  }
  txt
}

#' @keywords internal
report_generic <- function(x) {
  sprintf("A total sample size of %d participants provides %.0f%% power (%s).",
          x$n_total, 100 * x$power_achieved, x$method)
}

#' @keywords internal
report_continuous <- function(x) {
  p <- x$parameters
  sided <- if (x$sides == 2L) "two-sided" else "one-sided"
  switch(
    x$hypothesis,
    superiority = sprintf(
      paste0("A total sample size of %d participants (%s) provides %.0f%% ",
             "power to detect a difference in means of %s units, assuming a ",
             "common standard deviation of %s units, using a %s significance ",
             "level of %s (%s)."),
      x$n_total, describe_groups(x$n), 100 * x$power_achieved,
      format(p$effect), format(p$sd), sided, format(x$alpha), x$method),
    noninferiority = sprintf(
      paste0("A total sample size of %d participants (%s) provides %.0f%% ",
             "power to establish non-inferiority with a margin of %s units, ",
             "assuming a true difference in means of %s units and a common ",
             "standard deviation of %s units, using a one-sided significance ",
             "level of %s (%s)."),
      x$n_total, describe_groups(x$n), 100 * x$power_achieved,
      format(p$margin), format(p$effect), format(p$sd), format(x$alpha),
      x$method),
    superiority_margin = sprintf(
      paste0("A total sample size of %d participants (%s) provides %.0f%% ",
             "power to establish superiority by a margin of %s units, ",
             "assuming a true difference in means of %s units and a common ",
             "standard deviation of %s units, using a one-sided significance ",
             "level of %s (%s)."),
      x$n_total, describe_groups(x$n), 100 * x$power_achieved,
      format(p$margin), format(p$effect), format(p$sd), format(x$alpha),
      x$method),
    equivalence = sprintf(
      paste0("A total sample size of %d participants (%s) provides %.0f%% ",
             "power to establish equivalence within margins of \u00b1%s ",
             "units, assuming a true difference in means of %s units and a ",
             "common standard deviation of %s units, using two one-sided ",
             "tests each at significance level %s (%s)."),
      x$n_total, describe_groups(x$n), 100 * x$power_achieved,
      format(p$margin), format(p$effect), format(p$sd), format(x$alpha),
      x$method),
    report_generic(x)
  )
}

#' @keywords internal
report_two_proportions <- function(x) {
  p <- x$parameters
  sided <- if (x$sides == 2L) "two-sided" else "one-sided"
  pct <- function(v) sprintf("%.1f%%", 100 * v)
  switch(
    x$hypothesis,
    superiority = sprintf(
      paste0("A total sample size of %d participants (%s) provides %.0f%% ",
             "power to detect a difference between proportions of %s in the ",
             "treatment group and %s in the control group, using a %s ",
             "significance level of %s (%s)."),
      x$n_total, describe_groups(x$n), 100 * x$power_achieved,
      pct(p$p1), pct(p$p2), sided, format(x$alpha), x$method),
    noninferiority = sprintf(
      paste0("A total sample size of %d participants (%s) provides %.0f%% ",
             "power to establish non-inferiority with a margin of %s ",
             "percentage points, assuming true proportions of %s in the ",
             "treatment group and %s in the control group, using a ",
             "one-sided significance level of %s (%s)."),
      x$n_total, describe_groups(x$n), 100 * x$power_achieved,
      format(100 * p$margin), pct(p$p1), pct(p$p2), format(x$alpha),
      x$method),
    superiority_margin = sprintf(
      paste0("A total sample size of %d participants (%s) provides %.0f%% ",
             "power to establish superiority by a margin of %s percentage ",
             "points, assuming true proportions of %s in the treatment group ",
             "and %s in the control group, using a one-sided significance ",
             "level of %s (%s)."),
      x$n_total, describe_groups(x$n), 100 * x$power_achieved,
      format(100 * p$margin), pct(p$p1), pct(p$p2), format(x$alpha),
      x$method),
    equivalence = sprintf(
      paste0("A total sample size of %d participants (%s) provides %.0f%% ",
             "power to establish equivalence within margins of \u00b1%s ",
             "percentage points, assuming true proportions of %s in the ",
             "treatment group and %s in the control group, using two ",
             "one-sided tests each at significance level %s (%s)."),
      x$n_total, describe_groups(x$n), 100 * x$power_achieved,
      format(100 * p$margin), pct(p$p1), pct(p$p2), format(x$alpha),
      x$method),
    report_generic(x)
  )
}

#' @keywords internal
report_one_proportion <- function(x) {
  p <- x$parameters
  sided <- if (x$sides == 2L) "two-sided" else "one-sided"
  pct <- function(v) sprintf("%.1f%%", 100 * v)
  switch(
    x$hypothesis,
    superiority = sprintf(
      paste0("A sample size of %d participants provides %.0f%% power to ",
             "detect a true proportion of %s against a reference value of ",
             "%s, using a %s significance level of %s (%s)."),
      x$n_total, 100 * x$power_achieved, pct(p$p), pct(p$p0), sided,
      format(x$alpha), x$method),
    noninferiority = sprintf(
      paste0("A sample size of %d participants provides %.0f%% power to ",
             "establish non-inferiority to a reference value of %s with a ",
             "margin of %s percentage points, assuming a true proportion of ",
             "%s, using a one-sided significance level of %s (%s)."),
      x$n_total, 100 * x$power_achieved, pct(p$p0), format(100 * p$margin),
      pct(p$p), format(x$alpha), x$method),
    superiority_margin = sprintf(
      paste0("A sample size of %d participants provides %.0f%% power to ",
             "establish superiority over a reference value of %s by a ",
             "margin of %s percentage points, assuming a true proportion of ",
             "%s, using a one-sided significance level of %s (%s)."),
      x$n_total, 100 * x$power_achieved, pct(p$p0), format(100 * p$margin),
      pct(p$p), format(x$alpha), x$method),
    equivalence = sprintf(
      paste0("A sample size of %d participants provides %.0f%% power to ",
             "establish equivalence to a reference value of %s within ",
             "margins of \u00b1%s percentage points, assuming a true ",
             "proportion of %s, using two one-sided tests each at ",
             "significance level %s (%s)."),
      x$n_total, 100 * x$power_achieved, pct(p$p0), format(100 * p$margin),
      pct(p$p), format(x$alpha), x$method),
    report_generic(x)
  )
}

#' @keywords internal
describe_groups <- function(n) {
  if (length(n) == 1L) {
    sprintf("%d in a single group", n[1])
  } else if (length(unique(n)) == 1L) {
    sprintf("%d per group", n[1])
  } else {
    paste(sprintf("%d in the %s group", n, names(n)), collapse = ", ")
  }
}
