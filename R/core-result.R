# core-result.R
# Every procedure returns a `clinss_result` object rather than a bare number.
# The object carries sample sizes, achieved power, all assumptions, the method
# name, and the literature reference, making results self-documenting and
# reproducible.

#' @keywords internal
new_clinss_result <- function(n, n_total, power_target, power_achieved,
                              alpha, sides, hypothesis, method, reference,
                              parameters, dropout = 0, n_enrolled = NULL,
                              endpoint = "continuous") {
  structure(
    list(
      n              = n,
      n_total        = n_total,
      n_enrolled     = n_enrolled,
      power_target   = power_target,
      power_achieved = power_achieved,
      alpha          = alpha,
      sides          = sides,
      hypothesis     = hypothesis,
      method         = method,
      reference      = reference,
      parameters     = parameters,
      dropout        = dropout,
      endpoint       = endpoint,
      call           = sys.call(-1)
    ),
    class = "clinss_result"
  )
}

#' Print a clinss_result object
#'
#' @param x A \code{clinss_result} object.
#' @param digits Number of significant digits for the achieved power.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.clinss_result <- function(x, digits = 4, ...) {
  cat(x$method, "\n")
  cat(strrep("-", nchar(x$method)), "\n")
  cat("Hypothesis:      ", x$hypothesis, "\n")
  for (i in seq_along(x$n)) {
    cat(sprintf("n (%s):%s%d\n", names(x$n)[i],
                strrep(" ", max(1, 11 - nchar(names(x$n)[i]))), x$n[i]))
  }
  cat("n (total):       ", x$n_total, "\n")
  if (!is.null(x$n_enrolled)) {
    cat(sprintf("n (enrolled):     %d  (allowing for %.0f%% dropout)\n",
                x$n_enrolled, 100 * x$dropout))
  }
  cat("Achieved power:  ", format(round(x$power_achieved, digits)), "\n")
  cat("Target power:    ", format(x$power_target), "\n")
  cat(sprintf("Alpha:            %s (%s-sided)\n", format(x$alpha),
              if (x$sides == 2L) "two" else "one"))
  invisible(x)
}

#' Summarise a clinss_result object
#'
#' Prints the standard \code{print()} output followed by the full list of
#' design assumptions and the literature reference.
#'
#' @param object A \code{clinss_result} object.
#' @param ... Ignored.
#'
#' @return \code{object}, invisibly.
#' @export
summary.clinss_result <- function(object, ...) {
  print(object)
  cat("\nAssumptions\n")
  for (nm in names(object$parameters)) {
    cat(sprintf("  %-22s %s\n", nm, format(object$parameters[[nm]])))
  }
  cat("\nReference\n  ", object$reference, "\n", sep = "")
  invisible(object)
}

#' Coerce a clinss_result to a data frame
#'
#' Returns a one-row data frame containing all key outputs and inputs.
#' Useful for building sensitivity-analysis tables over a grid of assumptions
#' with \code{lapply()} and \code{do.call(rbind, ...)}.
#'
#' @param x A \code{clinss_result} object.
#' @param ... Ignored.
#'
#' @return A one-row \code{data.frame}.
#' @export
as.data.frame.clinss_result <- function(x, ...) {
  base <- data.frame(
    method         = x$method,
    hypothesis     = x$hypothesis,
    n_total        = x$n_total,
    power_target   = x$power_target,
    power_achieved = x$power_achieved,
    alpha          = x$alpha,
    sides          = x$sides,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(x$n)) {
    base[[paste0("n_", names(x$n)[i])]] <- x$n[i]
  }
  base
}
