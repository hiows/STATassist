#' Check the assumptions a comparison rests on
#'
#' Runs the normality tests, the homogeneity of variance tests and the outlier
#' screen together, and reports them as one object. The same checks are attached
#' automatically to [compare_two_groups()] and [compare_multiple_groups()], so
#' an assumption is never silently ignored; this function is for looking at them
#' on their own, before a test has been chosen.
#'
#' A failed assumption never blocks an analysis and never causes a test to be
#' swapped for another one. It changes which member of the reported test family
#' deserves the most weight, and that judgement stays with the user. Skewed
#' groups make the rank-based and robust members more trustworthy than the
#' parametric one; unequal variances make Welch's and Brunner-Munzel's
#' treatments of the same data more trustworthy than the pooled ones.
#'
#' Each assumption is checked twice on purpose, by tests that fail differently:
#'
#' \describe{
#'   \item{Normality}{Shapiro-Wilk is the more powerful of the two and is the
#'     one to read first. The Kolmogorov-Smirnov test is fitted against a normal
#'     with the sample's own mean and standard deviation, which makes its
#'     p-value anti-conservative, so it disagreeing with Shapiro-Wilk usually
#'     means the departure is in the tails.}
#'   \item{Homogeneity of variance}{The Levene test is centred on the median,
#'     the Brown-Forsythe variant, and tolerates skew. Bartlett's test is more
#'     powerful when the groups really are normal and cannot tell unequal
#'     variances apart from heavy tails when they are not. The two disagreeing
#'     is itself evidence about normality.}
#' }
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to check.
#' @param group Optional grouping vector with one entry per row of `data`. When
#'   supplied, normality is checked within each level and the variance tests are
#'   run across them. Without it there is only one sample per feature, so the
#'   variance table is empty.
#' @param group_lv Group levels to keep, in display order. Defaults to the
#'   sorted unique values of `group`.
#' @param alpha Threshold applied to the p-values when setting the `normal_ok`
#'   and `variance_ok` flags of the summary table.
#' @param criterion,iqr_multiplier,z_threshold Passed to [screen_outliers()].
#' @param center,trim Centre used by the Levene test and its trimming
#'   proportion, passed through to the Brown-Forsythe variant.
#'
#' @return A `sa_diagnosis` object: a plain list carrying `schema_version`,
#'   `analysis`, `features`, `design`, `parameters`, `metadata` and four tables.
#'
#'   \describe{
#'     \item{`normality`}{One row per feature and group level: `features`,
#'       `group`, `n_used`, `shapiro_stat`, `shapiro_pval`, `ks_stat`,
#'       `ks_pval`, `skewness` and `excess_kurtosis`. The level column is named
#'       `group` to match [summarize_descriptive_stats()], and is `NA` when no
#'       grouping was supplied.}
#'     \item{`variance`}{One row per feature: `features`, `n_used`, `n_groups`,
#'       `levene_stat`, `levene_df1`, `levene_df2`, `levene_pval`,
#'       `bartlett_stat`, `bartlett_df` and `bartlett_pval`. Zero rows when no
#'       `group` was supplied.}
#'     \item{`outliers`}{The [screen_outliers()] table, one row per flagged
#'       observation.}
#'     \item{`summary`}{One row per feature: `features`, `n_levels`,
#'       `n_outliers`, `min_shapiro_pval`, `normal_ok` and `variance_ok`. The
#'       two flags are `NA` when the corresponding test could not be run.}
#'   }
#'
#'   The result is deliberately not an `sa_comparison`. Normality is a property
#'   of one sample and homogeneity a property of a set of them, so the two
#'   tables have different numbers of rows and would not fit a contract built
#'   around one row per feature.
#'
#' @seealso [screen_outliers()] for the outlier stage on its own, and
#'   [compare_multiple_groups()], whose `$diagnostics` slot holds the same
#'   tables for the data it tested.
#'
#' @references
#' Shapiro, S. S. and Wilk, M. B. (1965). An analysis of variance test for
#' normality (complete samples). *Biometrika*, 52(3-4), 591-611.
#'
#' Massey, F. J. (1951). The Kolmogorov-Smirnov test for goodness of fit.
#' *Journal of the American Statistical Association*, 46(253), 68-78.
#'
#' Brown, M. B. and Forsythe, A. B. (1974). Robust tests for the equality of
#' variances. *Journal of the American Statistical Association*, 69(346),
#' 364-367.
#'
#' Bartlett, M. S. (1937). Properties of sufficiency and statistical tests.
#' *Proceedings of the Royal Society A*, 160(901), 268-282.
#'
#' @examples
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#'
#' ## Within species, so a genuine species difference is not read as skew
#' d <- diagnose_distribution(iris, feats, iris$Species)
#' d
#' d$normality
#' d$variance
#' d$summary
#'
#' ## Petal length pooled over species is strongly bimodal, which every
#' ## normality test picks up once the grouping is taken away.
#' diagnose_distribution(iris, "Petal.Length")$normality
#'
#' @export
diagnose_distribution <- function(data,
                                  feats,
                                  group = NULL,
                                  group_lv = NULL,
                                  alpha = 0.05,
                                  criterion = c("iqr", "robust_z", "grubbs"),
                                  iqr_multiplier = 1.5,
                                  z_threshold = 3.5,
                                  center = c("median", "mean", "trimmed"),
                                  trim = 0.1) {

  criterion <- match.arg(criterion)
  center <- match.arg(center)
  sa_check_scalar_num(alpha, "alpha", 0, 1, lower_open = TRUE)
  sa_check_scalar_num(trim, "trim", 0, 0.5, upper_open = TRUE)

  split <- sa_split_for_screening(data, feats, group, group_lv)
  groups_used <- if (split$grouped) names(split$rows) else NA_character_

  per_feature <- lapply(feats, function(f) {
    samples <- lapply(split$rows, function(rows) {
      v <- split$data[[f]][rows]
      v[is.finite(v)]
    })
    stats::setNames(samples, groups_used)
  })
  names(per_feature) <- feats

  outliers <- screen_outliers(data, feats, group, group_lv, criterion,
                              iqr_multiplier, z_threshold, alpha)

  sa_new_diagnosis(
    features   = feats,
    design     = list(group_lv = if (split$grouped) groups_used else NULL,
                      grouped  = split$grouped),
    parameters = list(alpha = alpha, criterion = criterion,
                      iqr_multiplier = iqr_multiplier,
                      z_threshold = z_threshold, center = center,
                      trim = trim),
    normality  = sa_normality_table(per_feature, feats),
    variance   = sa_variance_table(per_feature, feats, split$grouped, center,
                                   trim),
    outliers   = outliers,
    alpha      = alpha
  )
}


#' Normality and shape, one row per feature and group level
#'
#' A level that cannot be tested, because it is too small or constant, yields a
#' row of `NA` rather than being left out. Its absence would otherwise be
#' indistinguishable from the level not existing.
#'
#' The level column is named `group` rather than `level` so that it lines up
#' with the output of `summarize_descriptive_stats()`.
#'
#' @param per_feature Named list of named lists of numeric vectors.
#' @param feats Feature names, fixing the block order.
#'
#' @keywords internal
#' @noRd
sa_normality_table <- function(per_feature, feats) {
  blocks <- lapply(feats, function(f) {
    samples <- per_feature[[f]]
    # Indexed by position, not by name: an ungrouped diagnosis labels its single
    # sample `NA`, and `samples[[NA_character_]]` is not a lookup failure but a
    # silent empty vector.
    rows <- lapply(seq_along(samples), function(j) {
      group <- names(samples)[j]
      v <- samples[[j]]
      normality <- tryCatch(sa_shapiro(v), error = function(e) {
        c(shapiro_stat = NA_real_, shapiro_pval = NA_real_)
      })
      ks <- tryCatch(sa_ks_normal(v), error = function(e) {
        c(ks_stat = NA_real_, ks_pval = NA_real_)
      })
      data.frame(
        features        = f,
        group           = group,
        n_used          = length(v),
        shapiro_stat    = normality[["shapiro_stat"]],
        shapiro_pval    = normality[["shapiro_pval"]],
        ks_stat         = ks[["ks_stat"]],
        ks_pval         = ks[["ks_pval"]],
        skewness        = if (length(v) > 0L) sa_skewness(v) else NA_real_,
        excess_kurtosis = if (length(v) > 0L) sa_kurtosis(v) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}


#' Homogeneity of variance, one row per feature
#'
#' @param grouped Whether more than one sample per feature exists. Without a
#'   grouping there is nothing to compare variances across, so the table is
#'   empty rather than full of `NA`.
#'
#' @keywords internal
#' @noRd
sa_variance_table <- function(per_feature, feats, grouped, center, trim) {
  columns <- c("features", "n_used", "n_groups", "levene_stat", "levene_df1",
               "levene_df2", "levene_pval", "bartlett_stat", "bartlett_df",
               "bartlett_pval")
  if (!grouped) {
    empty <- data.frame(features = character(0), stringsAsFactors = FALSE)
    for (nm in setdiff(columns, "features")) {
      empty[[nm]] <- numeric(0)
    }
    return(empty)
  }

  rows <- lapply(feats, function(f) {
    samples <- per_feature[[f]]
    levene <- tryCatch(sa_levene(samples, center, trim), error = function(e) {
      c(levene_stat = NA_real_, levene_df1 = NA_real_, levene_df2 = NA_real_,
        levene_pval = NA_real_)
    })
    bartlett <- tryCatch(sa_bartlett(samples), error = function(e) {
      c(bartlett_stat = NA_real_, bartlett_df = NA_real_,
        bartlett_pval = NA_real_)
    })
    data.frame(features = f, n_used = sum(lengths(samples)),
               n_groups = length(samples),
               as.list(levene), as.list(bartlett),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[columns]
}


#' Assemble a diagnosis object
#'
#' @keywords internal
#' @noRd
sa_new_diagnosis <- function(features, design, parameters, normality, variance,
                             outliers, alpha) {
  summary_table <- data.frame(
    features = features,
    n_levels = vapply(features, function(f) {
      sum(normality$features == f)
    }, integer(1)),
    n_outliers = vapply(features, function(f) {
      sum(outliers$features == f)
    }, integer(1)),
    stringsAsFactors = FALSE
  )
  # The worst level decides: a comparison is only as normal as its least normal
  # group, so taking the minimum is what makes the flag mean something.
  summary_table$min_shapiro_pval <- vapply(features, function(f) {
    p <- normality$shapiro_pval[normality$features == f]
    if (all(is.na(p))) NA_real_ else min(p, na.rm = TRUE)
  }, numeric(1))
  summary_table$normal_ok <- summary_table$min_shapiro_pval > alpha
  summary_table$variance_ok <- if (nrow(variance) == 0L) {
    rep(NA, length(features))
  } else {
    variance$levene_pval[match(features, variance$features)] > alpha
  }
  rownames(summary_table) <- NULL

  structure(
    list(
      schema_version = sa_schema_version(),
      analysis       = "distribution_diagnosis",
      features       = features,
      design         = design,
      parameters     = parameters,
      normality      = normality,
      variance       = variance,
      outliers       = outliers,
      summary        = summary_table,
      metadata       = sa_metadata()
    ),
    class = c("sa_diagnosis", "sa_result")
  )
}


#' Attach the assumption checks a comparison rests on
#'
#' Called by the comparison scenarios with the samples they actually tested,
#' rather than with the original data. That matters: a paired design keeps
#' complete cases only, and a diagnosis run on the full column would describe a
#' different set of observations than the p-value it sits next to.
#'
#' @param per_feature Per feature, a named list of samples or, when `paired`, a
#'   subjects-by-conditions matrix.
#' @param feats Feature names.
#' @param group_lv Group levels.
#' @param paired Whether `per_feature` holds matrices.
#'
#' @keywords internal
#' @noRd
sa_diagnose_samples <- function(per_feature, feats, group_lv, paired,
                                alpha = 0.05) {
  samples_by_feature <- lapply(feats, function(f) {
    if (paired) {
      mat <- per_feature[[f]]
      stats::setNames(lapply(seq_len(ncol(mat)), function(j) mat[, j]),
                      group_lv)
    } else {
      per_feature[[f]]
    }
  })
  names(samples_by_feature) <- feats

  normality <- sa_normality_table(samples_by_feature, feats)
  # Homogeneity of variance across independent groups is not the assumption a
  # within-subject test makes; sphericity is, and the repeated measures ANOVA
  # row already carries Mauchly's test and both epsilon corrections.
  variance <- sa_variance_table(samples_by_feature, feats, !paired, "median",
                                0.1)

  list(normality = normality, variance = variance,
       summary = sa_new_diagnosis(feats,
                                  design = list(group_lv = group_lv,
                                                grouped = TRUE),
                                  parameters = list(alpha = alpha),
                                  normality = normality, variance = variance,
                                  outliers = data.frame(features = character(0),
                                                        stringsAsFactors = FALSE),
                                  alpha = alpha)$summary)
}


#' Print a distribution diagnosis
#'
#' Reports how many features failed each check rather than printing the tables.
#' Reach into `$normality`, `$variance` and `$outliers` for those.
#'
#' @param x A diagnosis, as returned by [diagnose_distribution()].
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' diagnose_distribution(iris, c("Sepal.Length", "Petal.Length"), iris$Species)
#'
#' @export
print.sa_diagnosis <- function(x, ...) {
  alpha <- x$parameters$alpha
  s <- x$summary

  cat("<sa_diagnosis> ", x$analysis, "\n", sep = "")
  cat("  features : ", length(x$features), "\n", sep = "")
  cat("  groups   : ",
      if (isTRUE(x$design$grouped)) {
        paste(x$design$group_lv, collapse = ", ")
      } else {
        "none, so no variance test"
      }, "\n", sep = "")
  cat("  settings : alpha = ", alpha, ", outlier criterion = ",
      x$parameters$criterion, "\n", sep = "")

  count <- function(flag) sum(!is.na(flag) & !flag)
  cat("\n  checks\n")
  cat("    normality  ", count(s$normal_ok), " of ", nrow(s),
      " feature(s) have a group failing Shapiro-Wilk at ", alpha, "\n",
      sep = "")
  if (nrow(x$variance) > 0L) {
    cat("    variance   ", count(s$variance_ok), " of ", nrow(s),
        " feature(s) fail Levene at ", alpha, "\n", sep = "")
  }
  cat("    outliers   ", nrow(x$outliers), " observation(s) flagged across ",
      sum(s$n_outliers > 0L), " feature(s)\n", sep = "")

  cat("\n  A failed check never changes which tests run. It changes which of\n")
  cat("  them deserves the most weight, and that judgement stays with you.\n")

  invisible(x)
}
