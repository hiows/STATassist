# Robust test kernels, written out rather than taken from a package.
#
# These used to be `lawstat::brunner.munzel.test()` and `WRS2::yuend()`. Both
# were dropped for the same reason: the Python port has to reproduce these
# numbers, and `scipy.stats.brunnermunzel()` reports no interval while no scipy
# or statsmodels function covers the dependent-samples Yuen test at all. Wrapping
# a different third party in each language would make the two implementations
# disagree over defaults rather than over arithmetic, which is the hardest kind
# of difference to track down. Written here once, the formula is the spec.
#
# Every function takes plain numeric vectors and returns a named numeric vector,
# so the port is a transcription rather than a redesign.

#' Two-sided or one-sided p-value from a t-distributed statistic
#'
#' `"greater"` always means the first sample exceeds the second, so the caller
#' must hand over a statistic that is positive in that case.
#'
#' @keywords internal
#' @noRd
sa_t_pval <- function(stat, df, alternative) {
  switch(alternative,
    two.sided = 2 * stats::pt(-abs(stat), df = df),
    greater   = stats::pt(stat, df = df, lower.tail = FALSE),
    less      = stats::pt(stat, df = df, lower.tail = TRUE)
  )
}


#' Confidence interval for an estimate with a t-distributed pivot
#'
#' A one-sided alternative leaves the side it does not test open, which is what
#' [stats::t.test()] and [stats::wilcox.test()] do, so the three tables of a
#' comparison can be read the same way. `bounds` sets what "open" means: an
#' unbounded quantity runs to infinity, a probability stops at 0 and 1.
#'
#' @keywords internal
#' @noRd
sa_t_ci <- function(est, se, df, alternative, conf_level,
                    bounds = c(-Inf, Inf)) {
  alpha <- 1 - conf_level
  switch(alternative,
    two.sided = est + c(-1, 1) * stats::qt(1 - alpha / 2, df = df) * se,
    greater   = c(est - stats::qt(1 - alpha, df = df) * se, bounds[2]),
    less      = c(bounds[1], est + stats::qt(1 - alpha, df = df) * se)
  )
}


#' Winsorise the tails of a sample
#'
#' The `floor(tr * n)` smallest values are pulled up to the next one and the same
#' number of largest values pulled down, leaving the length unchanged. Order is
#' preserved, which is what lets the result be used for a covariance as well as
#' for a variance.
#'
#' @param v Numeric vector without missing values.
#' @param tr Proportion winsorised at each tail, in `[0, 0.5)`.
#'
#' @keywords internal
#' @noRd
sa_winsorize <- function(v, tr) {
  n <- length(v)
  g <- floor(tr * n)
  sorted <- sort(v)
  lower <- sorted[g + 1L]
  upper <- sorted[n - g]
  v[v < lower] <- lower
  v[v > upper] <- upper
  v
}


#' Variance of a winsorised standard normal sample
#'
#' Winsorising shortens the tails, so a winsorised variance underestimates the
#' variance of the underlying normal. Dividing by this factor rescales it back,
#' which is what makes the robust effect size of `sa_yuen_paired()` comparable to
#' an ordinary standardised difference.
#'
#' Closed form of `integrate(z^2 * dnorm(z), qnorm(tr), qnorm(1 - tr))` plus the
#' two winsorised tails, using `integral(z^2 * dnorm(z)) = pnorm(z) - z * dnorm(z)`.
#' At `tr = 0.2` it gives 0.4120867, whose square root is the 0.642 constant
#' quoted for the 20 percent trimmed case.
#'
#' @keywords internal
#' @noRd
sa_winsorized_normal_var <- function(tr) {
  if (tr <= 0) {
    return(1)
  }
  z <- stats::qnorm(tr)
  (1 - 2 * tr) + 2 * z * stats::dnorm(z) + 2 * z^2 * tr
}


#' Brunner-Munzel test for two independent samples
#'
#' The nonparametric Behrens-Fisher problem: unlike the Wilcoxon rank-sum test it
#' does not assume the two distributions share a shape, so it stays valid when
#' the groups differ in spread.
#'
#' Orientation follows the rest of the package: `relative_effect` is
#' `P(X > Y) + 0.5 * P(X = Y)`, above 0.5 when `x` tends to be the larger sample,
#' and `bm_stat` is positive in that same case. Published presentations of the
#' test, and `lawstat::brunner.munzel.test()`, state the estimate the other way
#' round as `P(X < Y) + 0.5 * P(X = Y)`; the two are complements and the p-value
#' is identical either way. This direction is used so that `relative_effect`,
#' `mean_diff`, `hl_shift` and `trim_diff` all point the same way in the result
#' tables.
#'
#' @param x,y Numeric vectors without missing values, at least 2 long.
#' @param alternative `"two.sided"`, `"less"` or `"greater"`, where `"greater"`
#'   tests whether `x` exceeds `y`.
#' @param conf_level Confidence level of the reported interval.
#'
#' @return Named numeric vector: `relative_effect`, `bm_stat`, `df`, `pval`,
#'   `lower_conf`, `upper_conf`.
#'
#' @references
#' Brunner, E. and Munzel, U. (2000). The nonparametric Behrens-Fisher problem:
#' asymptotic theory and a small-sample approximation. *Biometrical Journal*,
#' 42(1), 17-25.
#'
#' @keywords internal
#' @noRd
sa_brunner_munzel <- function(x, y, alternative = "two.sided",
                              conf_level = 0.95) {
  n_x <- length(x)
  n_y <- length(y)

  # Ranks within each sample and within the pooled sample. The placement
  # r_pooled - r_within is what carries the information about the other group.
  r_x <- rank(x)
  r_y <- rank(y)
  r_pooled <- rank(c(x, y))
  r_pooled_x <- r_pooled[seq_len(n_x)]
  r_pooled_y <- r_pooled[n_x + seq_len(n_y)]

  m_x <- mean(r_pooled_x)
  m_y <- mean(r_pooled_y)

  v_x <- sum((r_pooled_x - r_x - m_x + (n_x + 1) / 2)^2) / (n_x - 1)
  v_y <- sum((r_pooled_y - r_y - m_y + (n_y + 1) / 2)^2) / (n_y - 1)

  pooled_var <- n_x * v_x + n_y * v_y
  if (pooled_var <= 0) {
    stop("the groups do not overlap, leaving the Brunner-Munzel variance ",
         "estimate at zero and the statistic undefined.", call. = FALSE)
  }

  # (m_x - m_y) rather than (m_y - m_x), so a positive statistic means x is the
  # larger sample.
  bm_stat <- n_x * n_y * (m_x - m_y) / (n_x + n_y) / sqrt(pooled_var)
  df <- pooled_var^2 /
    ((n_x * v_x)^2 / (n_x - 1) + (n_y * v_y)^2 / (n_y - 1))

  relative_effect <- 1 - (m_y - (n_y + 1) / 2) / n_x
  se <- sqrt(v_x / (n_x * n_y^2) + v_y / (n_y * n_x^2))
  ci <- sa_t_ci(relative_effect, se, df, alternative, conf_level,
                bounds = c(0, 1))

  c(relative_effect = relative_effect,
    bm_stat         = bm_stat,
    df              = df,
    pval            = sa_t_pval(bm_stat, df, alternative),
    lower_conf      = ci[1],
    upper_conf      = ci[2])
}


#' Yuen's trimmed mean test for two dependent samples
#'
#' Compares trimmed means using a standard error built from the winsorised
#' variances and their covariance, so the pairing is kept while outliers in
#' either sample lose their leverage.
#'
#' `robust_dz` is the trimmed mean difference over a robust estimate of the
#' standard deviation of the paired differences, rescaled by
#' `sa_winsorized_normal_var()` so that it reads on the same scale as Cohen's
#' `dz` when the differences are normal. `WRS2::yuend()` reports an explanatory
#' power measure here instead, which is computed as though the two samples were
#' independent and is returned unsigned; a signed, pairing-aware quantity is more
#' use in a table where every other estimate is signed.
#'
#' @param x,y Equal-length numeric vectors of complete pairs, no missing values.
#' @param tr Proportion trimmed at each tail, in `[0, 0.5)`.
#' @param alternative `"two.sided"`, `"less"` or `"greater"`, where `"greater"`
#'   tests whether `x` exceeds `y`.
#' @param conf_level Confidence level of the reported interval.
#'
#' @return Named numeric vector: `x_trim_mean`, `y_trim_mean`, `trim_diff`,
#'   `stderr`, `yuen_stat`, `df`, `pval`, `lower_conf`, `upper_conf`,
#'   `robust_dz`.
#'
#' @references
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population
#' variances. *Biometrika*, 61(1), 165-170.
#'
#' Algina, J., Keselman, H. J. and Penfield, R. D. (2005). An alternative to
#' Cohen's standardized mean difference effect size: a robust parameter and
#' confidence interval in the two independent groups case. *Psychological
#' Methods*, 10(3), 317-328.
#'
#' @keywords internal
#' @noRd
sa_yuen_paired <- function(x, y, tr = 0.2, alternative = "two.sided",
                           conf_level = 0.95) {
  n_pairs <- length(x)
  h <- n_pairs - 2L * floor(tr * n_pairs)

  win_x <- sa_winsorize(x, tr)
  win_y <- sa_winsorize(y, tr)

  # Sums of squared deviations, not variances: the h in the denominator below
  # replaces n, which is the whole point of the trimmed test.
  ss_x <- (n_pairs - 1) * stats::var(win_x)
  ss_y <- (n_pairs - 1) * stats::var(win_y)
  ss_xy <- (n_pairs - 1) * stats::var(win_x, win_y)

  # ss_x + ss_y - 2 * ss_xy is the sum of squares of (win_x - win_y), so it can
  # only reach zero when the winsorised differences are constant.
  stderr <- sqrt((ss_x + ss_y - 2 * ss_xy) / (h * (h - 1)))
  if (!is.finite(stderr) || stderr <= 0) {
    stop("the winsorised paired differences have zero variance, leaving the ",
         "standard error at zero and the statistic undefined.", call. = FALSE)
  }

  df <- h - 1
  x_trim_mean <- mean(x, trim = tr)
  y_trim_mean <- mean(y, trim = tr)
  trim_diff <- x_trim_mean - y_trim_mean
  yuen_stat <- trim_diff / stderr
  ci <- sa_t_ci(trim_diff, stderr, df, alternative, conf_level)

  win_diff_var <- stats::var(sa_winsorize(x - y, tr))
  robust_dz <- if (win_diff_var > 0) {
    trim_diff / sqrt(win_diff_var / sa_winsorized_normal_var(tr))
  } else {
    NA_real_
  }

  c(x_trim_mean = x_trim_mean,
    y_trim_mean = y_trim_mean,
    trim_diff   = trim_diff,
    stderr      = stderr,
    yuen_stat   = yuen_stat,
    df          = df,
    pval        = sa_t_pval(yuen_stat, df, alternative),
    lower_conf  = ci[1],
    upper_conf  = ci[2],
    robust_dz   = robust_dz)
}
