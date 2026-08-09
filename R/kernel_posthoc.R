# Pairwise kernels, each paired with the omnibus test that shares its
# assumptions. A rank-based omnibus test is never followed by a parametric
# comparison, which is the whole reason these are not interchangeable.
#
# Every function here returns a matrix with one row per pair, in the row order
# of `sa_level_pairs()`, and the columns `sa_posthoc_table()` expects. The
# estimate reads as `group_lv[j] - group_lv[i]` with `i < j`, the same direction
# rule `group_lv` fixes everywhere else: the reference is the first level, so it
# is the one being subtracted.


#' Column layout of one post-hoc pair
#'
#' @keywords internal
#' @noRd
sa_posthoc_columns <- function() {
  c("n1", "n2", "estimate", "stderr", "statistic", "df", "pval", "lower_conf",
    "upper_conf")
}


#' Assemble a pairwise result matrix from a per-pair function
#'
#' @param group_lv Group levels, fixing the pair order.
#' @param fun Function of the two level indices returning a named numeric vector
#'   with the columns of `sa_posthoc_columns()`.
#'
#' @keywords internal
#' @noRd
sa_pair_matrix <- function(group_lv, fun) {
  pairs <- sa_level_pairs(group_lv)
  rows <- lapply(seq_len(nrow(pairs)), function(p) {
    fun(pairs$i[p], pairs$j[p])[sa_posthoc_columns()]
  })
  matrix(unlist(rows, use.names = FALSE),
         nrow = nrow(pairs), byrow = TRUE,
         dimnames = list(NULL, sa_posthoc_columns()))
}


#' Tukey's honestly significant difference
#'
#' All pairwise mean differences judged against the studentised range, which
#' controls the error rate over the whole set of comparisons at once. The
#' p-values are therefore already family-wise and must not be adjusted again.
#'
#' Assumes equal variances, since every pair is judged against the same pooled
#' mean square error. That is the assumption it shares with the one-way ANOVA it
#' follows.
#'
#' @param samples List of numeric vectors, one per group level, no missing
#'   values, named by and ordered as `group_lv`.
#' @param conf_level Confidence level for the reported intervals.
#'
#' @return Matrix of `sa_posthoc_columns()`, one row per pair.
#'
#' @references
#' Tukey, J. W. (1949). Comparing individual means in the analysis of variance.
#' *Biometrics*, 5(2), 99-114.
#'
#' @keywords internal
#' @noRd
sa_tukey <- function(samples, conf_level = 0.95) {
  k <- length(samples)
  n <- lengths(samples)
  means <- vapply(samples, mean, numeric(1))
  df <- sum(n) - k

  ss_within <- sum(vapply(seq_len(k), function(i) {
    sum((samples[[i]] - means[i])^2)
  }, numeric(1)))
  ms_within <- ss_within / df
  if (ms_within <= 0) {
    stop("the pooled mean square error is zero, so no pairwise comparison can ",
         "be scaled.", call. = FALSE)
  }

  q_crit <- stats::qtukey(conf_level, k, df)

  sa_pair_matrix(names(samples), function(i, j) {
    estimate <- means[i] - means[j]
    # The studentised range is the range of k means over the standard error of
    # one mean, so the divisor carries a 1/2 that a two-sample t does not.
    stderr <- sqrt(ms_within / 2 * (1 / n[i] + 1 / n[j]))
    q_stat <- estimate / stderr
    sa_row(n1         = n[i],
           n2         = n[j],
           estimate   = estimate,
           stderr     = stderr,
           statistic  = q_stat,
           df         = df,
           pval       = stats::ptukey(abs(q_stat), k, df, lower.tail = FALSE),
           lower_conf = estimate - q_crit * stderr,
           upper_conf = estimate + q_crit * stderr)
  })
}


#' Games-Howell pairwise comparisons
#'
#' Tukey's procedure with the pooled variance replaced by a per-pair Welch
#' standard error and Welch degrees of freedom, so it stays valid when the
#' groups differ in spread or in size. The post-hoc partner of Welch's ANOVA.
#'
#' Like Tukey's test the p-values come from the studentised range and are
#' already family-wise, so they must not be adjusted again.
#'
#' @inheritParams sa_tukey
#'
#' @return Matrix of `sa_posthoc_columns()`, one row per pair.
#'
#' @references
#' Games, P. A. and Howell, J. F. (1976). Pairwise multiple comparison
#' procedures with unequal n's and/or variances: a Monte Carlo study.
#' *Journal of Educational Statistics*, 1(2), 113-125.
#'
#' @keywords internal
#' @noRd
sa_games_howell <- function(samples, conf_level = 0.95) {
  k <- length(samples)
  n <- lengths(samples)
  means <- vapply(samples, mean, numeric(1))
  vars <- vapply(samples, stats::var, numeric(1))

  sa_pair_matrix(names(samples), function(i, j) {
    v_i <- vars[i] / n[i]
    v_j <- vars[j] / n[j]
    if (v_i + v_j <= 0) {
      stop("both groups of the pair ", names(samples)[i], " - ",
           names(samples)[j], " have zero variance.", call. = FALSE)
    }
    estimate <- means[i] - means[j]
    stderr <- sqrt((v_i + v_j) / 2)
    df <- (v_i + v_j)^2 / (v_i^2 / (n[i] - 1) + v_j^2 / (n[j] - 1))
    q_stat <- estimate / stderr
    q_crit <- stats::qtukey(conf_level, k, df)
    sa_row(n1         = n[i],
           n2         = n[j],
           estimate   = estimate,
           stderr     = stderr,
           statistic  = q_stat,
           df         = df,
           pval       = stats::ptukey(abs(q_stat), k, df, lower.tail = FALSE),
           lower_conf = estimate - q_crit * stderr,
           upper_conf = estimate + q_crit * stderr)
  })
}


#' Dunn's pairwise rank comparisons
#'
#' Compares mean ranks taken from the pooled ranking the Kruskal-Wallis test
#' already computed, rather than re-ranking each pair on its own. Using the
#' pooled ranks is what keeps the post-hoc conclusions consistent with the
#' omnibus one; a set of separate rank-sum tests can contradict it.
#'
#' The variance carries a tie correction, so midranks do not inflate the
#' statistic. Unlike Tukey's test these p-values are not family-wise on their
#' own and are meant to be adjusted by the caller.
#'
#' @inheritParams sa_tukey
#'
#' @return Matrix of `sa_posthoc_columns()`, one row per pair. `estimate` is the
#'   mean rank difference and `df` is `NA`, the statistic being standard normal.
#'
#' @references
#' Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics*,
#' 6(3), 241-252.
#'
#' @keywords internal
#' @noRd
sa_dunn <- function(samples, conf_level = 0.95) {
  n <- lengths(samples)
  total <- sum(n)
  pooled <- unlist(samples, use.names = FALSE)
  ranks <- rank(pooled)
  split_at <- rep(seq_along(samples), n)
  mean_ranks <- vapply(seq_along(samples), function(i) {
    mean(ranks[split_at == i])
  }, numeric(1))

  tie_sizes <- table(pooled)
  tie_term <- sum(tie_sizes^3 - tie_sizes) / (12 * (total - 1))
  base_var <- total * (total + 1) / 12 - tie_term
  if (base_var <= 0) {
    stop("every observation is tied, leaving the rank variance at zero.",
         call. = FALSE)
  }

  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  sa_pair_matrix(names(samples), function(i, j) {
    estimate <- mean_ranks[i] - mean_ranks[j]
    stderr <- sqrt(base_var * (1 / n[i] + 1 / n[j]))
    z_stat <- estimate / stderr
    sa_row(n1         = n[i],
           n2         = n[j],
           estimate   = estimate,
           stderr     = stderr,
           statistic  = z_stat,
           df         = NA_real_,
           pval       = 2 * stats::pnorm(-abs(z_stat)),
           lower_conf = estimate - z_crit * stderr,
           upper_conf = estimate + z_crit * stderr)
  })
}


#' Yuen's trimmed mean test for two independent samples
#'
#' The independent counterpart of `sa_yuen_paired()`. Trimming both tails before
#' comparing means, and building the standard error from winsorised variances,
#' keeps a few extreme observations from deciding the result.
#'
#' @param x,y Numeric vectors without missing values.
#' @param tr Proportion trimmed at each tail, in `[0, 0.5)`.
#' @param alternative `"two.sided"`, `"less"` or `"greater"`, where `"greater"`
#'   tests whether `x` exceeds `y`.
#' @param conf_level Confidence level of the reported interval.
#'
#' @return Named numeric vector: `x_trim_mean`, `y_trim_mean`, `trim_diff`,
#'   `stderr`, `yuen_stat`, `df`, `pval`, `lower_conf`, `upper_conf`.
#'
#' @references
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population
#' variances. *Biometrika*, 61(1), 165-170.
#'
#' @keywords internal
#' @noRd
sa_yuen_independent <- function(x, y, tr = 0.2, alternative = "two.sided",
                                conf_level = 0.95) {
  n_x <- length(x)
  n_y <- length(y)
  h_x <- n_x - 2L * floor(tr * n_x)
  h_y <- n_y - 2L * floor(tr * n_y)
  if (h_x < 2L || h_y < 2L) {
    stop("fewer than 2 observations survive trimming ", tr,
         " from each tail (", h_x, " and ", h_y, ").", call. = FALSE)
  }

  d_x <- (n_x - 1) * stats::var(sa_winsorize(x, tr)) / (h_x * (h_x - 1))
  d_y <- (n_y - 1) * stats::var(sa_winsorize(y, tr)) / (h_y * (h_y - 1))
  stderr <- sqrt(d_x + d_y)
  if (!is.finite(stderr) || stderr <= 0) {
    stop("both winsorised samples are constant, leaving the standard error at ",
         "zero and the statistic undefined.", call. = FALSE)
  }

  df <- (d_x + d_y)^2 / (d_x^2 / (h_x - 1) + d_y^2 / (h_y - 1))
  x_trim_mean <- mean(x, trim = tr)
  y_trim_mean <- mean(y, trim = tr)
  trim_diff <- x_trim_mean - y_trim_mean
  yuen_stat <- trim_diff / stderr
  ci <- sa_t_ci(trim_diff, stderr, df, alternative, conf_level)

  c(x_trim_mean = x_trim_mean,
    y_trim_mean = y_trim_mean,
    trim_diff   = trim_diff,
    stderr      = stderr,
    yuen_stat   = yuen_stat,
    df          = df,
    pval        = sa_t_pval(yuen_stat, df, alternative),
    lower_conf  = ci[1],
    upper_conf  = ci[2])
}


#' Pairwise Yuen comparisons between independent groups
#'
#' The post-hoc partner of the trimmed mean ANOVA, run pair by pair on the same
#' trimming proportion the omnibus test used. Its p-values are not family-wise
#' and are meant to be adjusted by the caller.
#'
#' @inheritParams sa_tukey
#' @param tr Proportion trimmed at each tail, in `[0, 0.5)`.
#'
#' @return Matrix of `sa_posthoc_columns()`, one row per pair.
#'
#' @keywords internal
#' @noRd
sa_pairwise_yuen <- function(samples, tr = 0.2, conf_level = 0.95) {
  n <- lengths(samples)
  sa_pair_matrix(names(samples), function(i, j) {
    res <- sa_yuen_independent(samples[[i]], samples[[j]], tr = tr,
                               conf_level = conf_level)
    sa_row(n1         = n[i],
           n2         = n[j],
           estimate   = res[["trim_diff"]],
           stderr     = res[["stderr"]],
           statistic  = res[["yuen_stat"]],
           df         = res[["df"]],
           pval       = res[["pval"]],
           lower_conf = res[["lower_conf"]],
           upper_conf = res[["upper_conf"]])
  })
}


#' Pairwise paired t-tests between repeated conditions
#'
#' The post-hoc partner of repeated measures ANOVA. Each pair is tested on its
#' own differences rather than against the pooled residual error, so a pair
#' whose difference happens to be far more variable than the rest is not judged
#' as though it were not. Its p-values are meant to be adjusted by the caller.
#'
#' @param mat Subjects-by-conditions numeric matrix, complete, columns named by
#'   and ordered as `group_lv`.
#' @param conf_level Confidence level for the reported intervals.
#'
#' @return Matrix of `sa_posthoc_columns()`, one row per pair.
#'
#' @keywords internal
#' @noRd
sa_pairwise_paired_t <- function(mat, conf_level = 0.95) {
  n <- nrow(mat)
  sa_pair_matrix(colnames(mat), function(i, j) {
    res <- stats::t.test(mat[, i], mat[, j], paired = TRUE,
                         conf.level = conf_level)
    sa_row(n1         = n,
           n2         = n,
           estimate   = mean(mat[, i]) - mean(mat[, j]),
           stderr     = res$stderr,
           statistic  = res$statistic,
           df         = res$parameter,
           pval       = res$p.value,
           lower_conf = res$conf.int[1],
           upper_conf = res$conf.int[2])
  })
}


#' Conover's pairwise comparisons after a Friedman test
#'
#' Compares the within-subject rank sums the Friedman test already formed,
#' scaled by the residual variability of those same ranks and judged against a
#' t distribution. Working from the within-block ranking is what keeps the
#' post-hoc conclusions consistent with the omnibus one.
#'
#' Its p-values are not family-wise and are meant to be adjusted by the caller.
#'
#' @inheritParams sa_pairwise_paired_t
#'
#' @return Matrix of `sa_posthoc_columns()`, one row per pair. `estimate` is the
#'   rank sum difference.
#'
#' @references
#' Conover, W. J. (1999). *Practical Nonparametric Statistics*, 3rd edition.
#'
#' @keywords internal
#' @noRd
sa_conover <- function(mat, conf_level = 0.95) {
  n <- nrow(mat)
  k <- ncol(mat)
  # Ranked within each subject, so only the ordering a subject produces counts.
  ranks <- t(apply(mat, 1, rank))
  rank_sums <- colSums(ranks)

  # a - b is the sum of squares of the ranks left after the condition rank sums
  # are accounted for, which is the residual the pairwise scale is built from.
  a <- sum(ranks^2)
  b <- sum(rank_sums^2) / n

  df <- (n - 1L) * (k - 1L)
  variance <- 2 * n * (a - b) / df
  if (!is.finite(variance) || variance <= 0) {
    stop("every subject ranks the conditions identically, leaving the ",
         "residual rank variance at zero.", call. = FALSE)
  }
  stderr <- sqrt(variance)
  t_crit <- stats::qt(1 - (1 - conf_level) / 2, df)

  sa_pair_matrix(colnames(mat), function(i, j) {
    estimate <- rank_sums[i] - rank_sums[j]
    t_stat <- estimate / stderr
    sa_row(n1         = n,
           n2         = n,
           estimate   = estimate,
           stderr     = stderr,
           statistic  = t_stat,
           df         = df,
           pval       = 2 * stats::pt(-abs(t_stat), df),
           lower_conf = estimate - t_crit * stderr,
           upper_conf = estimate + t_crit * stderr)
  })
}
