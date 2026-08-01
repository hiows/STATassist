# Assumption-check kernels. Two callers use them: `diagnose_distribution()`,
# where they are the analysis, and the comparison scenarios, where they fill the
# `diagnostics` slot so that an assumption the test rests on is never silently
# ignored.
#
# A failed check never changes what gets run. It changes which member of the
# reported test family deserves the most weight, and that judgement stays with
# the user.


#' Shapiro-Wilk normality test
#'
#' @param v Numeric vector without missing values.
#'
#' @return Named numeric vector: `shapiro_stat`, `shapiro_pval`.
#'
#' @references
#' Shapiro, S. S. and Wilk, M. B. (1965). An analysis of variance test for
#' normality (complete samples). *Biometrika*, 52(3-4), 591-611.
#'
#' @keywords internal
#' @noRd
sa_shapiro <- function(v) {
  n <- length(v)
  # The engine's own limits, checked here so the message names the sample size
  # rather than reporting a failure from inside stats.
  if (n < 3L || n > 5000L) {
    stop("Shapiro-Wilk needs between 3 and 5000 observations, got ", n, ".",
         call. = FALSE)
  }
  res <- stats::shapiro.test(v)
  c(shapiro_stat = unname(res$statistic),
    shapiro_pval = res$p.value)
}


#' Kolmogorov-Smirnov goodness-of-fit test against a fitted normal
#'
#' The reference distribution is the normal with the sample's own mean and
#' standard deviation. Estimating the parameters from the same data the test
#' judges makes the p-value anti-conservative: it is too large, so the test
#' rejects normality less often than its nominal level says. That is a real
#' limitation of the test and the reason both normality checks are reported
#' together rather than one of them alone.
#'
#' @param v Numeric vector without missing values.
#'
#' @return Named numeric vector: `ks_stat`, `ks_pval`.
#'
#' @references
#' Massey, F. J. (1951). The Kolmogorov-Smirnov test for goodness of fit.
#' *Journal of the American Statistical Association*, 46(253), 68-78.
#'
#' @keywords internal
#' @noRd
sa_ks_normal <- function(v) {
  n <- length(v)
  if (n < 2L) {
    stop("the Kolmogorov-Smirnov test needs at least 2 observations, got ", n,
         ".", call. = FALSE)
  }
  spread <- stats::sd(v)
  if (!is.finite(spread) || spread <= 0) {
    stop("the sample is constant, so no normal reference distribution can be ",
         "fitted to it.", call. = FALSE)
  }
  # Ties make the exact p-value unavailable and stats warns once per call. The
  # warning says nothing the caller can act on here, and the grouped reporting
  # in sa_feature_table() would repeat it for every feature.
  res <- suppressWarnings(
    stats::ks.test(v, "pnorm", mean(v), spread)
  )
  c(ks_stat = unname(res$statistic),
    ks_pval = res$p.value)
}


#' Levene test for homogeneity of variance
#'
#' A one-way ANOVA on how far each observation sits from its own group centre.
#' The default centre is the median, which is the Brown-Forsythe variant: it
#' keeps the test honest when the groups are skewed, whereas centring on the
#' mean makes the test itself sensitive to the non-normality it is meant to
#' tolerate.
#'
#' @param samples List of numeric vectors, one per group level, no missing
#'   values.
#' @param center `"median"`, `"mean"` or `"trimmed"`.
#' @param trim Trimming proportion used when `center = "trimmed"`.
#'
#' @return Named numeric vector: `levene_stat`, `levene_df1`, `levene_df2`,
#'   `levene_pval`.
#'
#' @references
#' Brown, M. B. and Forsythe, A. B. (1974). Robust tests for the equality of
#' variances. *Journal of the American Statistical Association*, 69(346),
#' 364-367.
#'
#' @keywords internal
#' @noRd
sa_levene <- function(samples, center = "median", trim = 0.1) {
  centre_of <- switch(center,
    median  = function(v) stats::median(v),
    mean    = function(v) mean(v),
    trimmed = function(v) mean(v, trim = trim),
    stop("`center` must be one of: median, mean, trimmed.", call. = FALSE)
  )

  deviations <- lapply(samples, function(v) abs(v - centre_of(v)))
  res <- sa_oneway_anova(deviations)

  c(levene_stat = res[["f_stat"]],
    levene_df1  = res[["df1"]],
    levene_df2  = res[["df2"]],
    levene_pval = res[["pval"]])
}


#' Bartlett test for homogeneity of variance
#'
#' More powerful than the Levene test when the groups really are normal, and
#' misleading when they are not: it cannot tell unequal variances apart from
#' heavy tails. Read it next to the Levene result rather than instead of it.
#'
#' @inheritParams sa_levene
#'
#' @return Named numeric vector: `bartlett_stat`, `bartlett_df`,
#'   `bartlett_pval`.
#'
#' @references
#' Bartlett, M. S. (1937). Properties of sufficiency and statistical tests.
#' *Proceedings of the Royal Society A*, 160(901), 268-282.
#'
#' @keywords internal
#' @noRd
sa_bartlett <- function(samples) {
  res <- stats::bartlett.test(samples)
  c(bartlett_stat = unname(res$statistic),
    bartlett_df   = unname(res$parameter),
    bartlett_pval = res$p.value)
}


#' Grubbs test for a single outlier
#'
#' Tests whether the observation furthest from the mean is further out than a
#' normal sample of that size would produce. It assumes the rest of the sample
#' is normal and it looks at one observation, so it is the weakest of the three
#' screening rules and the only one that produces a p-value.
#'
#' @param v Numeric vector without missing values.
#'
#' @return Named numeric vector: `grubbs_stat`, `grubbs_pval`, `grubbs_index`,
#'   the position in `v` of the most extreme observation.
#'
#' @references
#' Grubbs, F. E. (1969). Procedures for detecting outlying observations in
#' samples. *Technometrics*, 11(1), 1-21.
#'
#' @keywords internal
#' @noRd
sa_grubbs <- function(v) {
  n <- length(v)
  if (n < 3L) {
    stop("the Grubbs test needs at least 3 observations, got ", n, ".",
         call. = FALSE)
  }
  spread <- stats::sd(v)
  if (!is.finite(spread) || spread <= 0) {
    stop("the sample is constant, so no observation can be called extreme.",
         call. = FALSE)
  }

  distance <- abs(v - mean(v))
  index <- which.max(distance)
  g <- distance[index] / spread

  # Inverting the Grubbs statistic gives a t on n - 2 degrees of freedom. The
  # factor n is the Bonferroni correction for having looked at whichever of the
  # n observations turned out to be furthest out, so it can exceed 1.
  denominator <- (n - 1)^2 - n * g^2
  pval <- if (denominator <= 0) {
    0
  } else {
    t_stat <- sqrt(n * (n - 2) * g^2 / denominator)
    min(1, n * 2 * stats::pt(-t_stat, n - 2L))
  }

  c(grubbs_stat  = g,
    grubbs_pval  = pval,
    grubbs_index = index)
}


#' Flag outlying observations in one numeric vector
#'
#' Returns a logical flag per observation rather than a cleaned vector. Which
#' observations to keep is a decision about the experiment, not about the
#' arithmetic, so the package never makes it.
#'
#' @param v Numeric vector, missing values allowed and never flagged.
#' @param criterion `"iqr"`, `"robust_z"` or `"grubbs"`.
#' @param iqr_multiplier Fence width for `criterion = "iqr"`.
#' @param z_threshold Cut-off for `criterion = "robust_z"`.
#' @param alpha Significance level for `criterion = "grubbs"`.
#'
#' @return List with `flag`, a logical vector the length of `v`, and `score`,
#'   the numeric quantity the rule thresholded.
#'
#' @references
#' Iglewicz, B. and Hoaglin, D. C. (1993). *How to Detect and Handle Outliers*.
#'
#' @keywords internal
#' @noRd
sa_flag_outliers <- function(v, criterion = "iqr", iqr_multiplier = 1.5,
                             z_threshold = 3.5, alpha = 0.05) {
  usable <- is.finite(v)
  flag <- rep(FALSE, length(v))
  score <- rep(NA_real_, length(v))
  clean <- v[usable]
  if (length(clean) < 3L) {
    return(list(flag = flag, score = score))
  }

  if (criterion == "iqr") {
    q <- unname(stats::quantile(clean, c(0.25, 0.75)))
    iqr <- q[2] - q[1]
    if (iqr > 0) {
      # How far past the nearer quartile the value sits, measured in IQR units,
      # so the score itself does not depend on `iqr_multiplier` and the two can
      # be compared across calls that used different fences. Values inside the
      # box score negative.
      score[usable] <- pmax(q[1] - clean, clean - q[2]) / iqr
      flag[usable] <- score[usable] > iqr_multiplier
    }
  } else if (criterion == "robust_z") {
    # The median and MAD are used instead of the mean and SD because a single
    # extreme value inflates the SD enough to hide itself.
    spread <- stats::mad(clean)
    if (spread > 0) {
      score[usable] <- abs(clean - stats::median(clean)) / spread
      flag[usable] <- score[usable] > z_threshold
    }
  } else if (criterion == "grubbs") {
    res <- tryCatch(sa_grubbs(clean), error = function(e) NULL)
    if (!is.null(res)) {
      score[usable][res[["grubbs_index"]]] <- res[["grubbs_stat"]]
      if (res[["grubbs_pval"]] <= alpha) {
        flag[usable][res[["grubbs_index"]]] <- TRUE
      }
    }
  } else {
    stop("`criterion` must be one of: iqr, robust_z, grubbs.", call. = FALSE)
  }

  list(flag = flag, score = score)
}
