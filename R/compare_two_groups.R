#' Run every applicable two-group test at once
#'
#' Compares exactly two group levels across any number of numeric features and
#' returns a parametric, a rank-based and a robust test side by side, together
#' with the fold change between the two groups. Nothing is chosen on the user's
#' behalf: reporting all of them makes disagreement between them visible, which
#' is the situation where the choice of test actually matters.
#'
#' Which member of each family is used depends on `paired`:
#'
#' | family  | `paired = FALSE`          | `paired = TRUE`                 |
#' |---------|---------------------------|---------------------------------|
#' | t       | Welch's t-test            | Paired t-test                   |
#' | Wilcoxon| Rank-sum (Mann-Whitney U) | Signed-rank                     |
#' | robust  | Brunner-Munzel            | Yuen's trimmed mean (dependent) |
#'
#' Direction is set once, by the order of `group_lv`, and every quantity in the
#' result follows it. `alternative = "greater"` tests whether `group_lv[1]`
#' exceeds `group_lv[2]` in all three families, and `mean_diff`, `hl_shift`,
#' `trim_diff`, `fold_change` and `relative_effect` are all above their null
#' value when `group_lv[1]` is the larger group.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to test.
#'   One output row per entry.
#' @param group Grouping vector with one entry per row of `data`.
#' @param group_lv Character vector of exactly two group levels. The first is
#'   treated as `x`, the second as `y`, so all differences read as
#'   `group_lv[1] - group_lv[2]` and all ratios as `group_lv[1] / group_lv[2]`.
#'   Rows belonging to any other level are dropped.
#' @param id Optional pairing key with one entry per row of `data`, used only
#'   when `paired = TRUE`. Supplying it matches observations by key instead of
#'   by row order, which is the safer choice whenever the rows may have been
#'   reordered or a subject may be missing from one group. Ids present in only
#'   one group are dropped.
#' @param alternative Direction of the alternative hypothesis, one of
#'   `"two.sided"`, `"less"` or `"greater"`.
#' @param paired Logical. If `TRUE`, observations are treated as matched. See
#'   `id` for how the pairs are formed.
#' @param conf_level Confidence level for all reported intervals.
#' @param tr Trimming proportion for Yuen's test, in `[0, 0.5)`. Ignored when
#'   `paired = FALSE`.
#' @param fc_mean Which centre the fold change divides, `"arith"` for the
#'   arithmetic mean or `"geom"` for the geometric mean. The geometric mean
#'   requires strictly positive values and is the usual choice for
#'   concentration-like data.
#' @param p_adjust Multiplicity adjustment applied across `feats` within each
#'   test table, passed to [stats::p.adjust()]. Use `"none"` to disable.
#'
#' @return A `sa_comparison` object: a plain list, so it survives being written
#'   out as JSON and read back in another language, with an S3 class on top that
#'   only supplies [print()]. Its elements are
#'
#'   \describe{
#'     \item{`schema_version`}{Version of this layout, `"0.1.0"`.}
#'     \item{`analysis`}{`"two_group_comparison"`.}
#'     \item{`features`}{The feature names, in the row order every table uses.}
#'     \item{`design`}{`group_lv`, `paired`, `pairing` (`"order"`, `"id"` or
#'       `NA` when not paired), `n_dropped` (rows removed for belonging to a
#'       level outside `group_lv`) and `unmatched_ids`.}
#'     \item{`parameters`}{`alternative`, `conf_level`, `tr`, `fc_mean` and
#'       `p_adjust`, as used.}
#'     \item{`effect`}{One row per feature: `x_center`, `y_center` (the two
#'       centres `fc_mean` selected), `fold_change` and `log2fc`.}
#'     \item{`tests`}{`t_test`, `wilcox_test` and `robust_test`, described
#'       below.}
#'     \item{`test_info`}{Per test, the method `id`, a readable `label` and
#'       whether it was the paired variant.}
#'     \item{`metadata`}{`package_version`, `r_version`, `platform` and an
#'       ISO-8601 `timestamp`.}
#'   }
#'
#'   Every table in `tests` has one row per feature and starts with `features`,
#'   the per-group sample sizes `n_x` / `n_y` and `n_used` (total observations
#'   for independent samples, complete pairs for paired samples), and carries
#'   `pval`, `pval_adj`, `lower_conf` and `upper_conf`. The remaining columns
#'   are:
#'
#'   \describe{
#'     \item{`t_test`}{`x_mean`, `y_mean`, `mean_diff`, `stderr`, `t_stat`,
#'       `df`. Group means are computed directly, so the columns are identical
#'       for paired and independent designs.}
#'     \item{`wilcox_test`}{`hl_shift` (Hodges-Lehmann location shift, the
#'       pseudo-median of differences when paired) and `w_stat`.}
#'     \item{`robust_test`, independent}{`relative_effect`
#'       (`P(X > Y) + 0.5 * P(X = Y)`, above 0.5 when `group_lv[1]` is the
#'       larger group), `bm_stat` and `df`. The interval is on the probability
#'       scale, so a one-sided alternative leaves it open at 0 or at 1 rather
#'       than at infinity.}
#'     \item{`robust_test`, paired}{`x_trim_mean`, `y_trim_mean`, `trim_diff`,
#'       `stderr`, `yuen_stat`, `df` and `robust_dz`, a robust counterpart of
#'       Cohen's `dz`.}
#'   }
#'
#' @section Pairing:
#' With `paired = TRUE` and no `id`, the only available information is row
#' order: the first row of `group_lv[1]` is matched with the first row of
#' `group_lv[2]`, and so on. Both groups must then have the same number of
#' rows. Note that this cannot detect rows that have been reordered, so a data
#' set sorted differently in each group would produce wrong pairs and no
#' complaint. Passing `id` removes that failure mode: pairs are matched on the
#' key, ids appearing in only one group are dropped with a message, and an id
#' repeated within a group is an error.
#'
#' @details
#' Features that cannot be tested do not abort the run. Their row is filled
#' with `NA` and all such features are reported together in a single warning.
#' Informational engine warnings, such as an exact p-value being unavailable
#' because of ties, are grouped into one `message()`.
#'
#' Missing values are handled per feature: independent samples drop `NA` within
#' each group, paired samples keep only complete pairs. `n_x`, `n_y` and
#' `n_used` therefore vary between features when the data are incomplete. The
#' fold change is computed from those same observations, so it never rests on a
#' different subset of the data than the p-value beside it.
#'
#' @seealso [estimate_significance()] to reduce the result to one significance
#'   verdict per feature, and [draw_grouped_boxplot()] to visualise the same
#'   input.
#'
#' @references
#' Welch, B. L. (1947). The generalization of Student's problem when several
#' different population variances are involved. *Biometrika*, 34(1-2), 28-35.
#'
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods.
#' *Biometrics Bulletin*, 1(6), 80-83.
#'
#' Mann, H. B. and Whitney, D. R. (1947). On a test of whether one of two
#' random variables is stochastically larger than the other. *Annals of
#' Mathematical Statistics*, 18(1), 50-60.
#'
#' Brunner, E. and Munzel, U. (2000). The nonparametric Behrens-Fisher problem:
#' asymptotic theory and a small-sample approximation. *Biometrical Journal*,
#' 42(1), 17-25.
#'
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population
#' variances. *Biometrika*, 61(1), 165-170.
#'
#' @examples
#' ## Independent samples: two of the three iris species
#' iris2 <- iris[iris$Species != "setosa", ]
#' res <- compare_two_groups(
#'   data     = iris2,
#'   feats    = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   group    = iris2$Species,
#'   group_lv = c("versicolor", "virginica")
#' )
#' res
#' res$tests$t_test
#' res$tests$robust_test
#' res$effect
#'
#' ## setosa is perfectly separated from the others on petal size, which leaves
#' ## the Brunner-Munzel variance at zero. Those rows come back NA with a
#' ## warning rather than aborting the other features.
#' suppressWarnings(
#'   compare_two_groups(
#'     data     = iris[iris$Species != "virginica", ],
#'     feats    = c("Sepal.Length", "Petal.Length"),
#'     group    = iris$Species[iris$Species != "virginica"],
#'     group_lv = c("setosa", "versicolor")
#'   )$tests$robust_test
#' )
#'
#' ## Paired samples: sleep holds 10 subjects under 2 drugs, listed in the same
#' ## subject order within each group.
#' paired_res <- compare_two_groups(
#'   data     = sleep["extra"],
#'   feats    = "extra",
#'   group    = sleep$group,
#'   group_lv = c("1", "2"),
#'   paired   = TRUE,
#'   alternative = "less"
#' )
#' paired_res$tests$robust_test
#'
#' ## Same data with the second group shuffled. Row order pairing silently uses
#' ## the wrong partners, while `id` recovers the correct result.
#' shuffled <- rbind(sleep[1:10, ], sleep[10 + c(4, 9, 1, 7, 2, 10, 3, 6, 8, 5), ])
#' by_order <- compare_two_groups(
#'   data = shuffled["extra"], feats = "extra",
#'   group = shuffled$group, group_lv = c("1", "2"), paired = TRUE
#' )
#' by_id <- compare_two_groups(
#'   data = shuffled["extra"], feats = "extra",
#'   group = shuffled$group, group_lv = c("1", "2"),
#'   id = shuffled$ID, paired = TRUE
#' )
#' ## The means agree but the paired standard error, and so the p-value, do not.
#' rbind(order = by_order$tests$t_test,
#'       id    = by_id$tests$t_test)
#'
#' @export
compare_two_groups <- function(data,
                               feats,
                               group,
                               group_lv,
                               id = NULL,
                               alternative = c("two.sided", "less", "greater"),
                               paired = FALSE,
                               conf_level = 0.95,
                               tr = 0.2,
                               fc_mean = c("arith", "geom"),
                               p_adjust = "BH") {

  alternative <- match.arg(alternative)
  fc_mean <- match.arg(fc_mean)
  sa_check_flag(paired, "paired")
  sa_check_scalar_num(conf_level, "conf_level", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)
  sa_check_scalar_num(tr, "tr", 0, 0.5, upper_open = TRUE)
  sa_check_p_adjust(p_adjust, "p_adjust")

  if (!is.null(id) && !paired) {
    warning("`id` is only used to form pairs and is ignored when ",
            "`paired = FALSE`. Set `paired = TRUE` if the observations are ",
            "matched.", call. = FALSE)
  }

  input <- sa_validate_wide_input(data, feats, group, group_lv, id = id,
                                  n_levels = 2L)
  data <- input$data
  feats <- input$feats
  group <- input$group
  id <- input$id
  group_lv <- levels(group)

  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }

  if (!paired) {
    idx_x <- which(group == group_lv[1])
    idx_y <- which(group == group_lv[2])
    unmatched <- character(0)
  } else {
    pairing <- if (is.null(id)) {
      sa_pair_by_order(group, group_lv)
    } else {
      sa_pair_by_id(id, group, group_lv)
    }
    idx_x <- pairing$idx_x
    idx_y <- pairing$idx_y
    unmatched <- pairing$unmatched
    if (length(unmatched) > 0L) {
      message("Dropped ", length(unmatched),
              " id(s) present in only one group: ",
              paste(unmatched, collapse = ", "), ".")
    }
  }

  # Samples are taken straight out of the wide columns. Reshaping to long format
  # first and then subsetting was what let the group column be picked up from
  # the calling frame instead of the data.
  samples <- lapply(feats, function(f) {
    x <- data[[f]][idx_x]
    y <- data[[f]][idx_y]
    if (paired) {
      keep <- !is.na(x) & !is.na(y)
      list(x = x[keep], y = y[keep])
    } else {
      list(x = x[!is.na(x)], y = y[!is.na(y)])
    }
  })
  names(samples) <- feats

  effect <- sa_fold_change(samples, feats, group_lv, fc_mean)


  # T-test
  t_label <- if (paired) "Paired t-test" else "Welch's t-test"
  t_columns <- c("n_x", "n_y", "n_used", "x_mean", "y_mean", "mean_diff",
                 "stderr", "t_stat", "df", "pval", "lower_conf", "upper_conf")

  t_result <- sa_feature_table(feats, t_columns, t_label, p_adjust = p_adjust,
    fun = function(i) {
      s <- samples[[i]]
      n_x <- length(s$x)
      n_y <- length(s$y)
      if (n_x < 2L || n_y < 2L) {
        stop("needs at least 2 usable observations per group, got ",
             n_x, " and ", n_y, ".", call. = FALSE)
      }
      res <- stats::t.test(s$x, s$y, alternative = alternative,
                          paired = paired, conf.level = conf_level)
      # t.test()$estimate is two group means when independent but a single mean
      # difference when paired, so the means are computed here instead to keep
      # one column layout for both designs.
      sa_row(n_x         = n_x,
             n_y         = n_y,
             n_used      = if (paired) n_x else n_x + n_y,
             x_mean      = mean(s$x),
             y_mean      = mean(s$y),
             mean_diff   = mean(s$x) - mean(s$y),
             stderr      = res$stderr,
             t_stat      = res$statistic,
             df          = res$parameter,
             pval        = res$p.value,
             lower_conf  = res$conf.int[1],
             upper_conf  = res$conf.int[2])
    })


  # Wilcoxon Test
  w_label <- if (paired) {
    "Wilcoxon signed-rank test"
  } else {
    "Wilcoxon rank sum test (Mann-Whitney U test)"
  }
  w_columns <- c("n_x", "n_y", "n_used", "hl_shift", "w_stat", "pval",
                 "lower_conf", "upper_conf")

  w_result <- sa_feature_table(feats, w_columns, w_label, p_adjust = p_adjust,
    fun = function(i) {
      s <- samples[[i]]
      n_x <- length(s$x)
      n_y <- length(s$y)
      if (n_x < 1L || n_y < 1L) {
        stop("needs at least 1 usable observation per group, got ",
             n_x, " and ", n_y, ".", call. = FALSE)
      }
      res <- stats::wilcox.test(s$x, s$y, alternative = alternative,
                                paired = paired, conf.int = TRUE,
                                conf.level = conf_level)
      sa_row(n_x        = n_x,
             n_y        = n_y,
             n_used     = if (paired) n_x else n_x + n_y,
             hl_shift   = res$estimate,
             w_stat     = res$statistic,
             pval       = res$p.value,
             lower_conf = res$conf.int[1],
             upper_conf = res$conf.int[2])
    })


  # Robust Test
  if (!paired) {
    ## Brunner-Munzel Test
    robust_label <- "Brunner-Munzel test"
    robust_columns <- c("n_x", "n_y", "n_used", "relative_effect", "bm_stat",
                        "df", "pval", "lower_conf", "upper_conf")

    robust_result <- sa_feature_table(feats, robust_columns, robust_label,
      p_adjust = p_adjust,
      fun = function(i) {
        s <- samples[[i]]
        n_x <- length(s$x)
        n_y <- length(s$y)
        if (n_x < 2L || n_y < 2L) {
          stop("needs at least 2 usable observations per group, got ",
               n_x, " and ", n_y, ".", call. = FALSE)
        }
        c(sa_row(n_x = n_x, n_y = n_y, n_used = n_x + n_y),
          sa_brunner_munzel(s$x, s$y, alternative = alternative,
                            conf_level = conf_level))
      })
  } else {
    ## Yuen's trimmed mean test
    robust_label <- "Yuen's trimmed mean test for dependent samples"
    robust_columns <- c("n_x", "n_y", "n_used", "x_trim_mean", "y_trim_mean",
                        "trim_diff", "stderr", "yuen_stat", "df", "pval",
                        "lower_conf", "upper_conf", "robust_dz")

    robust_result <- sa_feature_table(feats, robust_columns, robust_label,
      p_adjust = p_adjust,
      fun = function(i) {
        s <- samples[[i]]
        n_pairs <- length(s$x)
        h <- n_pairs - 2L * floor(tr * n_pairs)
        if (h < 2L) {
          stop("only ", h, " observation(s) survive trimming ", tr,
               " from each tail of ", n_pairs, " pair(s); 2 are needed.",
               call. = FALSE)
        }
        c(sa_row(n_x = n_pairs, n_y = n_pairs, n_used = n_pairs),
          sa_yuen_paired(s$x, s$y, tr = tr, alternative = alternative,
                         conf_level = conf_level))
      })
  }


  sa_new_comparison(
    analysis  = "two_group_comparison",
    features  = feats,
    design    = list(
      group_lv      = group_lv,
      paired        = paired,
      pairing       = if (!paired) {
        NA_character_
      } else if (is.null(id)) {
        "order"
      } else {
        "id"
      },
      n_dropped     = input$n_dropped,
      unmatched_ids = unmatched
    ),
    parameters = list(
      alternative = alternative,
      conf_level  = conf_level,
      tr          = if (paired) tr else NA_real_,
      fc_mean     = fc_mean,
      p_adjust    = p_adjust
    ),
    effect    = effect,
    tests     = list(
      t_test      = t_result,
      wilcox_test = w_result,
      robust_test = robust_result
    ),
    test_info = list(
      t_test = list(
        id     = if (paired) "paired_t_test" else "welch_t_test",
        label  = t_label,
        paired = paired
      ),
      wilcox_test = list(
        id     = if (paired) "wilcoxon_signed_rank" else "mann_whitney_u",
        label  = w_label,
        paired = paired
      ),
      robust_test = list(
        id     = if (paired) "yuen_paired" else "brunner_munzel",
        label  = robust_label,
        paired = paired
      )
    ),
    subclass  = "sa_two_group"
  )
}
