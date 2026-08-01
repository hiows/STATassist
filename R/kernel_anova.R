# Omnibus kernels for three or more groups, written to the same rule as
# `kernel_robust.R`: plain numeric vectors in, a named numeric vector out, no
# fitted model kept anywhere. `stats::oneway.test()` and friends would each
# return a different htest shape, and the Python port would then have to
# reproduce four different unpacking rules instead of four formulas.
#
# None of the omnibus tests reports a confidence interval. An omnibus test says
# that the conditions are not all alike; it does not say by how much, and there
# is no single quantity for an interval to be about. The intervals of a
# multi-group comparison live in the post-hoc table, where each row is one
# contrast and does have a scale of its own. The `lower_conf` and `upper_conf`
# columns are therefore present and `NA` in every table here, which the result
# contract allows: it requires that the columns exist, not that they are finite.


#' Split a numeric vector into one sample per group level
#'
#' Missing values are dropped inside each level, which is the independent-sample
#' rule the rest of the package uses. Levels left with too few observations are
#' an error rather than a silently smaller design.
#'
#' @param values Numeric vector, missing values included.
#' @param group Factor with the levels to split on, same length as `values`.
#' @param n_min Smallest acceptable number of usable observations per level.
#'
#' @keywords internal
#' @noRd
sa_split_groups <- function(values, group, n_min = 2L) {
  samples <- lapply(levels(group), function(lv) {
    v <- values[group == lv]
    v[!is.na(v)]
  })
  names(samples) <- levels(group)

  sizes <- lengths(samples)
  short <- names(samples)[sizes < n_min]
  if (length(short) > 0L) {
    stop("needs at least ", n_min, " usable observation(s) per group; ",
         paste0(short, " = ", sizes[short], collapse = ", "), ".",
         call. = FALSE)
  }
  samples
}


#' One-way analysis of variance
#'
#' The equal-variance omnibus F test. Written out rather than taken from
#' [stats::oneway.test()] because the sums of squares are needed anyway: Tukey's
#' post-hoc test runs on the same mean square error, and computing it twice from
#' two different code paths is how the two end up disagreeing.
#'
#' @param samples List of numeric vectors, one per group level, no missing
#'   values.
#'
#' @return Named numeric vector: `n_used`, `n_groups`, `f_stat`, `df1`, `df2`,
#'   `eta_sq`, `omega_sq`, `pval`, `lower_conf`, `upper_conf`.
#'
#' @references
#' Fisher, R. A. (1925). *Statistical Methods for Research Workers*.
#'
#' Okada, K. (2013). Is omega squared less biased? A comparison of three major
#' effect size indices in one-way ANOVA. *Behaviormetrika*, 40(2), 129-147.
#'
#' @keywords internal
#' @noRd
sa_oneway_anova <- function(samples) {
  k <- length(samples)
  n <- lengths(samples)
  total <- sum(n)
  df1 <- k - 1L
  df2 <- total - k
  if (df2 < 1L) {
    stop("needs more observations than groups to leave any residual degrees ",
         "of freedom.", call. = FALSE)
  }

  means <- vapply(samples, mean, numeric(1))
  grand <- sum(n * means) / total

  ss_between <- sum(n * (means - grand)^2)
  ss_within <- sum(vapply(seq_len(k), function(i) {
    sum((samples[[i]] - means[i])^2)
  }, numeric(1)))
  ss_total <- ss_between + ss_within

  if (ss_within <= 0) {
    stop("every group has zero variance, leaving the mean square error at ",
         "zero and the F statistic undefined.", call. = FALSE)
  }

  ms_within <- ss_within / df2
  f_stat <- (ss_between / df1) / ms_within

  c(n_used     = total,
    n_groups   = k,
    f_stat     = f_stat,
    df1        = df1,
    df2        = df2,
    eta_sq     = ss_between / ss_total,
    # Omega squared subtracts the variance the grouping would explain by chance
    # alone, so it can go negative when the groups are indistinguishable. That
    # is not clipped here: a negative value is the estimate saying so.
    omega_sq   = (ss_between - df1 * ms_within) / (ss_total + ms_within),
    pval       = stats::pf(f_stat, df1, df2, lower.tail = FALSE),
    lower_conf = NA_real_,
    upper_conf = NA_real_)
}


#' Welch's heteroscedastic one-way analysis of variance
#'
#' Weights each group by its own precision instead of pooling the variances, so
#' unequal spreads or unequal group sizes no longer distort the test. The effect
#' size columns are the ordinary sums-of-squares ones: they describe how far
#' apart the group means are, which does not change with the choice of test.
#'
#' @inheritParams sa_oneway_anova
#'
#' @return The same columns as `sa_oneway_anova()`.
#'
#' @references
#' Welch, B. L. (1951). On the comparison of several mean values: an alternative
#' approach. *Biometrika*, 38(3-4), 330-336.
#'
#' @keywords internal
#' @noRd
sa_welch_anova <- function(samples) {
  k <- length(samples)
  n <- lengths(samples)
  if (any(n < 2L)) {
    stop("Welch's ANOVA needs at least 2 observations per group to estimate ",
         "a within-group variance.", call. = FALSE)
  }

  means <- vapply(samples, mean, numeric(1))
  vars <- vapply(samples, stats::var, numeric(1))
  if (any(vars <= 0)) {
    zero <- names(samples)[vars <= 0]
    stop("group(s) with zero variance leave the Welch weight infinite: ",
         paste(zero, collapse = ", "), ".", call. = FALSE)
  }

  w <- n / vars
  sum_w <- sum(w)
  weighted_mean <- sum(w * means) / sum_w

  lambda <- sum((1 - w / sum_w)^2 / (n - 1))
  numerator <- sum(w * (means - weighted_mean)^2) / (k - 1)
  denominator <- 1 + 2 * (k - 2) / (k^2 - 1) * lambda
  f_stat <- numerator / denominator

  df1 <- k - 1L
  df2 <- 1 / (3 / (k^2 - 1) * lambda)

  ss_reference <- sa_oneway_anova(samples)

  c(n_used     = sum(n),
    n_groups   = k,
    f_stat     = f_stat,
    df1        = df1,
    df2        = df2,
    eta_sq     = ss_reference[["eta_sq"]],
    omega_sq   = ss_reference[["omega_sq"]],
    pval       = stats::pf(f_stat, df1, df2, lower.tail = FALSE),
    lower_conf = NA_real_,
    upper_conf = NA_real_)
}


#' Yuen's trimmed mean one-way analysis of variance
#'
#' The robust member of the independent omnibus family: Welch's construction
#' applied to trimmed means and winsorised variances, so heavy tails and stray
#' observations lose the leverage they have over an ordinary F test.
#'
#' `robust_eta_sq` is defined here rather than borrowed. It is the spread of the
#' trimmed means around their unweighted centre, divided by that spread plus the
#' mean rescaled winsorised variance, where the rescaling by `(1 - 2 * tr)^2`
#' puts a winsorised variance back on the scale of the trimmed mean it belongs
#' to. It is zero when the trimmed means coincide, approaches one as they
#' separate, and does not change if every observation is multiplied by a
#' constant. It is not Wilcox's explanatory measure and does not reproduce it.
#'
#' @inheritParams sa_oneway_anova
#' @param tr Proportion trimmed at each tail, in `[0, 0.5)`.
#'
#' @return Named numeric vector: `n_used`, `n_groups`, `f_stat`, `df1`, `df2`,
#'   `robust_eta_sq`, `pval`, `lower_conf`, `upper_conf`.
#'
#' @references
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population
#' variances. *Biometrika*, 61(1), 165-170.
#'
#' @keywords internal
#' @noRd
sa_yuen_anova <- function(samples, tr = 0.2) {
  k <- length(samples)
  n <- lengths(samples)

  h <- n - 2L * floor(tr * n)
  if (any(h < 2L)) {
    short <- names(samples)[h < 2L]
    stop("fewer than 2 observations survive trimming ", tr,
         " from each tail in group(s): ", paste(short, collapse = ", "), ".",
         call. = FALSE)
  }

  trim_means <- vapply(samples, mean, numeric(1), trim = tr)
  win_vars <- vapply(samples, function(v) stats::var(sa_winsorize(v, tr)),
                     numeric(1))
  if (any(win_vars <= 0)) {
    flat <- names(samples)[win_vars <= 0]
    stop("group(s) whose winsorised values are constant leave the trimmed ",
         "weight infinite: ", paste(flat, collapse = ", "), ".", call. = FALSE)
  }

  # d is the squared standard error of one trimmed mean; its reciprocal is the
  # Welch weight, exactly as `sa_welch_anova()` uses n / var.
  d <- (n - 1) * win_vars / (h * (h - 1))
  w <- 1 / d
  sum_w <- sum(w)
  weighted_mean <- sum(w * trim_means) / sum_w

  lambda <- sum((1 - w / sum_w)^2 / (h - 1))
  numerator <- sum(w * (trim_means - weighted_mean)^2) / (k - 1)
  f_stat <- numerator / (1 + 2 * (k - 2) / (k^2 - 1) * lambda)

  df1 <- k - 1L
  df2 <- 1 / (3 / (k^2 - 1) * lambda)

  between <- sum((trim_means - mean(trim_means))^2) / k
  within <- mean(win_vars / (1 - 2 * tr)^2)

  c(n_used        = sum(n),
    n_groups      = k,
    f_stat        = f_stat,
    df1           = df1,
    df2           = df2,
    robust_eta_sq = between / (between + within),
    pval          = stats::pf(f_stat, df1, df2, lower.tail = FALSE),
    lower_conf    = NA_real_,
    upper_conf    = NA_real_)
}


#' Kruskal-Wallis rank sum test
#'
#' The rank-based omnibus test. It asks whether one group tends to produce
#' larger values than another, which is only a statement about medians when the
#' distributions share a shape.
#'
#' @inheritParams sa_oneway_anova
#'
#' @return Named numeric vector: `n_used`, `n_groups`, `h_stat`, `df`,
#'   `epsilon_sq`, `eta_sq_rank`, `pval`, `lower_conf`, `upper_conf`.
#'
#' @references
#' Kruskal, W. H. and Wallis, W. A. (1952). Use of ranks in one-criterion
#' variance analysis. *Journal of the American Statistical Association*,
#' 47(260), 583-621.
#'
#' Tomczak, M. and Tomczak, E. (2014). The need to report effect size estimates
#' revisited. *Trends in Sport Sciences*, 1(21), 19-25.
#'
#' @keywords internal
#' @noRd
sa_kruskal <- function(samples) {
  k <- length(samples)
  total <- sum(lengths(samples))
  pooled <- unlist(samples, use.names = FALSE)
  # Every value tied means the tie correction divides by zero, and kruskal.test
  # returns NaN rather than complaining. Refusing here is what turns it into an
  # NA row with a named reason, like every other test that cannot run.
  if (length(unique(pooled)) < 2L) {
    stop("every observation takes the same value, so the ranks carry no ",
         "information and the tie correction is undefined.", call. = FALSE)
  }
  res <- stats::kruskal.test(samples)
  h <- unname(res$statistic)

  c(n_used      = total,
    n_groups    = k,
    h_stat      = h,
    df          = unname(res$parameter),
    # Both rescale H onto [0, 1]; epsilon squared divides by the largest value H
    # could take, eta squared by the residual degrees of freedom.
    epsilon_sq  = h * (total + 1) / (total^2 - 1),
    eta_sq_rank = (h - k + 1) / (total - k),
    pval        = res$p.value,
    lower_conf  = NA_real_,
    upper_conf  = NA_real_)
}


#' One-way repeated measures analysis of variance
#'
#' Removes the between-subject variation before testing the conditions, which is
#' what makes a within-subject design more sensitive than the same number of
#' independent observations.
#'
#' The uncorrected F test assumes sphericity: that every pair of conditions has
#' the same variance of differences. Mauchly's test reports whether that holds
#' and the Greenhouse-Geisser and Huynh-Feldt epsilons say how badly it fails.
#' Both corrected p-values are returned alongside the uncorrected one instead of
#' one being chosen, since which to trust is a judgement about the design.
#'
#' @param mat Subjects-by-conditions numeric matrix, complete, at least 2 rows
#'   and 3 columns.
#'
#' @return Named numeric vector: `n_used`, `n_groups`, `f_stat`, `df1`, `df2`,
#'   `partial_eta_sq`, `gen_eta_sq`, `mauchly_w`, `mauchly_pval`, `gg_eps`,
#'   `pval_gg`, `hf_eps`, `pval_hf`, `pval`, `lower_conf`, `upper_conf`.
#'
#' @references
#' Mauchly, J. W. (1940). Significance test for sphericity of a normal
#' n-variate distribution. *Annals of Mathematical Statistics*, 11(2), 204-209.
#'
#' Greenhouse, S. W. and Geisser, S. (1959). On methods in the analysis of
#' profile data. *Psychometrika*, 24(2), 95-112.
#'
#' Huynh, H. and Feldt, L. S. (1976). Estimation of the Box correction for
#' degrees of freedom from sample data in randomized block and split-plot
#' designs. *Journal of Educational Statistics*, 1(1), 69-82.
#'
#' Bakeman, R. (2005). Recommended effect size statistics for repeated measures
#' designs. *Behavior Research Methods*, 37(3), 379-384.
#'
#' @keywords internal
#' @noRd
sa_rm_anova <- function(mat) {
  n <- nrow(mat)
  k <- ncol(mat)
  if (n < 2L) {
    stop("needs at least 2 complete subjects, got ", n, ".", call. = FALSE)
  }

  grand <- mean(mat)
  condition_means <- colMeans(mat)
  subject_means <- rowMeans(mat)

  ss_condition <- n * sum((condition_means - grand)^2)
  ss_subject <- k * sum((subject_means - grand)^2)
  ss_total <- sum((mat - grand)^2)
  ss_error <- ss_total - ss_condition - ss_subject

  df1 <- k - 1L
  df2 <- (n - 1L) * (k - 1L)
  if (ss_error <= 0) {
    stop("the subject-by-condition residuals are all zero, leaving the F ",
         "statistic undefined.", call. = FALSE)
  }
  ms_error <- ss_error / df2
  f_stat <- (ss_condition / df1) / ms_error

  sphericity <- sa_sphericity(mat)
  gg <- sphericity[["gg_eps"]]
  hf <- sphericity[["hf_eps"]]

  c(n_used         = n,
    n_groups       = k,
    f_stat         = f_stat,
    df1            = df1,
    df2            = df2,
    partial_eta_sq = ss_condition / (ss_condition + ss_error),
    # Generalised eta squared keeps the between-subject variance in the
    # denominator, which is what makes it comparable with the eta squared of an
    # independent design measuring the same thing.
    gen_eta_sq     = ss_condition / (ss_condition + ss_subject + ss_error),
    mauchly_w      = sphericity[["mauchly_w"]],
    mauchly_pval   = sphericity[["mauchly_pval"]],
    gg_eps         = gg,
    pval_gg        = stats::pf(f_stat, df1 * gg, df2 * gg, lower.tail = FALSE),
    hf_eps         = hf,
    pval_hf        = stats::pf(f_stat, df1 * hf, df2 * hf, lower.tail = FALSE),
    pval           = stats::pf(f_stat, df1, df2, lower.tail = FALSE),
    lower_conf     = NA_real_,
    upper_conf     = NA_real_)
}


#' Mauchly's sphericity test and the two epsilon corrections
#'
#' All three quantities come from the same orthonormal contrast of the
#' condition covariance matrix, so they are computed together.
#'
#' The Huynh-Feldt epsilon can exceed one, which is meaningless as a correction
#' factor, so it is capped. Both epsilons are also floored at the lower bound
#' `1 / (k - 1)`, the value they take when sphericity fails as badly as it can.
#'
#' @param mat Subjects-by-conditions numeric matrix.
#'
#' @return Named numeric vector: `mauchly_w`, `mauchly_pval`, `gg_eps`,
#'   `hf_eps`.
#'
#' @keywords internal
#' @noRd
sa_sphericity <- function(mat) {
  n <- nrow(mat)
  k <- ncol(mat)
  p <- k - 1L
  f <- n - 1L

  lower_bound <- 1 / p
  na_out <- c(mauchly_w = NA_real_, mauchly_pval = NA_real_,
              gg_eps = lower_bound, hf_eps = lower_bound)
  if (n <= k) {
    # The condition covariance is singular below this point, so both the
    # determinant and its eigenvalues stop meaning anything. Falling back to the
    # lower bound applies the most conservative correction available rather than
    # reporting no correction at all.
    return(na_out)
  }

  covariance <- stats::cov(mat)

  # Orthonormal contrasts spanning the k - 1 differences between conditions.
  # qr.Q() of the deviation basis gives a basis orthogonal to the unit vector,
  # which is what sphericity is defined relative to.
  contrasts <- qr.Q(qr(diag(k) - 1 / k), complete = FALSE)[, seq_len(p),
                                                           drop = FALSE]
  transformed <- t(contrasts) %*% covariance %*% contrasts
  eigenvalues <- eigen(transformed, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= 0)) {
    return(na_out)
  }

  w <- prod(eigenvalues) / (sum(eigenvalues) / p)^p

  # Mauchly's statistic is only asymptotically chi-square, so the p-value uses
  # the two-term expansion rather than the leading term alone.
  #
  # `weight` is written with p, the rank of the contrast, throughout.
  # stats::mauchly.test() evaluates the same expression with the number of
  # conditions in the `3 * p` term and the contrast rank everywhere else, so the
  # two disagree in the fourth decimal place. The single-symbol form is the
  # published one and is what is used here.
  rho <- 1 - (2 * p^2 + p + 2) / (6 * p * f)
  chi_sq <- -f * rho * log(w)
  chi_df <- p * (p + 1) / 2 - 1
  weight <- (p + 2) * (p - 1) * (p - 2) * (2 * p^3 + 6 * p^2 + 3 * p + 2) /
    (288 * (f * p * rho)^2)
  lead <- stats::pchisq(chi_sq, chi_df, lower.tail = FALSE)
  correction <- stats::pchisq(chi_sq, chi_df + 4, lower.tail = FALSE)
  mauchly_pval <- lead + weight * (correction - lead)

  gg <- sum(eigenvalues)^2 / (p * sum(eigenvalues^2))
  gg <- min(max(gg, lower_bound), 1)
  hf <- (n * p * gg - 2) / (p * (f - p * gg))
  hf <- min(max(hf, lower_bound), 1)

  c(mauchly_w    = w,
    mauchly_pval = mauchly_pval,
    gg_eps       = gg,
    hf_eps       = hf)
}


#' Friedman rank sum test
#'
#' The rank-based counterpart of repeated measures ANOVA: observations are
#' ranked within each subject, so only the ordering of the conditions inside a
#' subject matters and no distributional assumption is made across subjects.
#'
#' @param mat Subjects-by-conditions numeric matrix, complete.
#'
#' @return Named numeric vector: `n_used`, `n_groups`, `chi_sq`, `df`,
#'   `kendalls_w`, `pval`, `lower_conf`, `upper_conf`.
#'
#' @references
#' Friedman, M. (1937). The use of ranks to avoid the assumption of normality
#' implicit in the analysis of variance. *Journal of the American Statistical
#' Association*, 32(200), 675-701.
#'
#' Kendall, M. G. and Babington Smith, B. (1939). The problem of m rankings.
#' *Annals of Mathematical Statistics*, 10(3), 275-287.
#'
#' @keywords internal
#' @noRd
sa_friedman <- function(mat) {
  n <- nrow(mat)
  k <- ncol(mat)
  # A subject who gives every condition the same value contributes only ties.
  # With no subject ranking anything, friedman.test returns NaN in silence; the
  # refusal here turns that into an NA row with a reason attached.
  if (all(apply(mat, 1L, function(row) length(unique(row)) < 2L))) {
    stop("no subject distinguishes the conditions, so the within-subject ",
         "ranks carry no information.", call. = FALSE)
  }
  res <- stats::friedman.test(mat)
  chi_sq <- unname(res$statistic)

  c(n_used     = n,
    n_groups   = k,
    chi_sq     = chi_sq,
    df         = unname(res$parameter),
    # Kendall's W is the Friedman statistic expressed as agreement between
    # subjects: 0 when the subjects rank the conditions independently, 1 when
    # they all produce the same ranking.
    kendalls_w = chi_sq / (n * (k - 1)),
    pval       = res$p.value,
    lower_conf = NA_real_,
    upper_conf = NA_real_)
}
