# Internal summary kernels shared by the descriptive functions. Shape estimators
# are written out rather than taken from a package: both are closed-form
# expressions over the sample moments, and `40_dependency_policy.md` rules out
# adding an Import for a wrapper around arithmetic this small.

#' Column layout of one descriptive summary row
#'
#' Named in one place so that the row builder and the all-NA fallback can never
#' disagree about which columns exist or in what order.
#'
#' @keywords internal
#' @noRd
sa_describe_columns <- function() {
  c("n", "n_missing", "mean", "sd", "var", "se", "cv",
    "min", "q1", "median", "q3", "max", "iqr",
    "out_lower_bound", "out_upper_bound", "mad",
    "skewness", "excess_kurtosis")
}


#' Sample skewness, G1
#'
#' The bias-corrected estimator SAS and SPSS report, matching
#' `e1071::skewness(type = 2)`. It needs at least three observations and a
#' non-zero spread; outside that the correction divides by zero, so the result
#' is `NA` rather than a `NaN` that would travel on into a summary table.
#'
#' @keywords internal
#' @noRd
sa_skewness <- function(v) {
  n <- length(v)
  m <- mean(v)
  m2 <- sum((v - m)^2) / n
  if (n < 3L || m2 <= 0) {
    return(NA_real_)
  }
  g1 <- (sum((v - m)^3) / n) / m2^1.5
  g1 * sqrt(n * (n - 1)) / (n - 2)
}


#' Sample excess kurtosis, G2
#'
#' The counterpart of `sa_skewness()`, matching `e1071::kurtosis(type = 2)`.
#' Excess, so a normal sample sits near zero. Needs at least four observations.
#'
#' @keywords internal
#' @noRd
sa_kurtosis <- function(v) {
  n <- length(v)
  m <- mean(v)
  m2 <- sum((v - m)^2) / n
  if (n < 4L || m2 <= 0) {
    return(NA_real_)
  }
  g2 <- (sum((v - m)^4) / n) / m2^2 - 3
  ((n + 1) * g2 + 6) * (n - 1) / ((n - 2) * (n - 3))
}


#' Describe one numeric vector
#'
#' Values and names are written as a single `c(name = value)` call, so a
#' quantity cannot end up under the wrong label.
#'
#' @param x Numeric vector, missing and non-finite values included. They are
#'   counted into `n_missing` and left out of every other statistic, which
#'   keeps a single `Inf` from turning the whole row into `Inf` or `NaN`.
#'
#' @return Named numeric vector with the columns of `sa_describe_columns()`.
#'
#' @keywords internal
#' @noRd
sa_describe_vector <- function(x) {
  v <- x[is.finite(x)]
  n <- length(v)

  if (n == 0L) {
    out <- sa_na_row(sa_describe_columns())
    out[["n"]] <- 0
    out[["n_missing"]] <- length(x)
    return(out)
  }

  m <- mean(v)
  s <- stats::sd(v)
  q <- unname(stats::quantile(v, c(0.25, 0.5, 0.75)))
  iqr <- q[3] - q[1]

  c(n               = n,
    n_missing       = length(x) - n,
    mean            = m,
    sd              = s,
    var             = stats::var(v),
    se              = s / sqrt(n),
    cv              = s / m,
    min             = min(v),
    q1              = q[1],
    median          = q[2],
    q3              = q[3],
    max             = max(v),
    iqr             = iqr,
    out_lower_bound = q[1] - 1.5 * iqr,
    out_upper_bound = q[3] + 1.5 * iqr,
    mad             = stats::mad(v),
    skewness        = sa_skewness(v),
    excess_kurtosis = sa_kurtosis(v))
}
