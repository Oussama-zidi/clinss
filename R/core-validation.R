# core-validation.R
# Internal input validation helpers. Every public procedure calls these so
# that error messages are consistent across the whole package.
# None of these functions are exported.

#' @keywords internal
check_prob <- function(x, name, open = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single numeric value.", name), call. = FALSE)
  }
  if (open && (x <= 0 || x >= 1)) {
    stop(sprintf("`%s` must be strictly between 0 and 1.", name), call. = FALSE)
  }
  invisible(x)
}

#' @keywords internal
check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive number.", name), call. = FALSE)
  }
  invisible(x)
}

#' @keywords internal
check_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single numeric value.", name), call. = FALSE)
  }
  invisible(x)
}

#' @keywords internal
check_choice <- function(x, name, choices) {
  if (!is.character(x) || length(x) != 1L || !(x %in% choices)) {
    stop(sprintf("`%s` must be one of: %s.", name,
                 paste(shQuote(choices), collapse = ", ")), call. = FALSE)
  }
  invisible(x)
}
