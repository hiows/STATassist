#' Run every applicable multi-group test at once
#'
#' Compares three or more group levels across any number of numeric features and
#' returns the omnibus tests side by side, each followed by the post-hoc
#' procedure that shares its assumptions. As with [compare_two_groups()],
#' nothing is chosen on the user's behalf: reporting the parametric, the
#' rank-based and the robust result together makes disagreement between them
#' visible, which is the situation where the choice of test actually matters.
#'
#' Which family runs depends on `paired`:
#'
#' | slot           | `paired = FALSE`         | `paired = TRUE`             |
#' |----------------|--------------------------|-----------------------------|
#' | `anova_test`   | One-way ANOVA            | Repeated measures ANOVA     |
#' | `welch_test`   | Welch's ANOVA            | not applicable              |
#' | `robust_test`  | Yuen's trimmed mean ANOVA| not applicable              |
#' | `kruskal_test` | Kruskal-Wallis           | Friedman                    |
#'
#' Each one is followed by its own post-hoc procedure, never by a borrowed one:
#' one-way ANOVA by Tukey's HSD, Welch's ANOVA by Games-Howell, the trimmed mean
#' ANOVA by pairwise Yuen tests, Kruskal-Wallis by Dunn's test, repeated
#' measures ANOVA by pairwise paired t-tests and Friedman by Conover's test. A
#' rank-based omnibus test is therefore never followed by a parametric
#' comparison.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to test.
#'   One output row per entry.
#' @param group Grouping vector with one entry per row of `data`.
#' @param group_lv Character vector of at least three group levels. The first is
#'   the reference the `effect` table is expressed against, and the order fixes
#'   the direction of every post-hoc contrast. Rows belonging to any other level
#'   are dropped.
#' @param id Subject identifier with one entry per row of `data`. Required when
#'   `paired = TRUE` and ignored otherwise.
#' @param paired Logical. If `TRUE`, the levels of `group` are treated as
#'   repeated conditions measured on the subjects named by `id`.
#' @param conf_level Confidence level for the post-hoc intervals.
#' @param tr Trimming proportion for the trimmed mean ANOVA, in `[0, 0.5)`.
#'   Ignored when `paired = TRUE`.
#' @param posthoc Logical. If `FALSE`, no pairwise stage runs and `$posthoc` is
#'   an empty list.
#' @param posthoc_alpha A feature enters the post-hoc stage when its omnibus
#'   `pval_adj` is at or below this value. Set it to 1 to compare every feature
#'   regardless of its omnibus result.
#' @param fc_mean Which centre the fold change divides, `"arith"` for the
#'   arithmetic mean or `"geom"` for the geometric mean. Defaults to `"geom"`
#'   when `input_scale = "log2"` and to `"arith"` otherwise.
#' @param input_scale The scale `data` arrives on, `"raw"` or `"log2"`. On the
#'   log2 scale each observation is raised back through `2^x` before the centres
#'   are taken, so the ratios mean what they do for raw input. This changes the
#'   `effect` table only, never the tests. See [compare_two_groups()] for the
#'   full account.
#' @param p_adjust Multiplicity adjustment applied across `feats` within each
#'   omnibus table, passed to [stats::p.adjust()]. Use `"none"` to disable.
#' @param posthoc_p_adjust Multiplicity adjustment applied across the contrasts
#'   within each feature of a post-hoc table. Ignored for Tukey's HSD and
#'   Games-Howell, whose p-values are already family-wise.
#' @param diagnose Logical. If `TRUE`, the normality, homogeneity of variance
#'   and outlier checks the tests rest on are attached as `$diagnostics`.
#'
#' @return A `sa_comparison` object with the same layout
#'   [compare_two_groups()] returns, plus a populated `posthoc` slot. Its
#'   elements are
#'
#'   \describe{
#'     \item{`schema_version`}{Version of this layout, `"0.2.1"`.}
#'     \item{`analysis`}{`"multi_group_comparison"`.}
#'     \item{`features`}{The feature names, in the row order every table uses.}
#'     \item{`design`}{`group_lv`, `paired`, `pairing`, `n_dropped` and
#'       `unmatched_ids`.}
#'     \item{`parameters`}{The analysis choices as used, plus `n_posthoc`, the
#'       number of features that entered the pairwise stage per test.}
#'     \item{`effect`}{One row per feature, described below.}
#'     \item{`tests`}{The omnibus tables, one row per feature.}
#'     \item{`posthoc`}{One table per omnibus test, one row per feature and pair
#'       of levels, carrying `features`, `contrast`, `group1`, `group2`, `n1`,
#'       `n2`, `estimate`, `stderr`, `statistic`, `df`, `pval`, `pval_adj`,
#'       `lower_conf` and `upper_conf`.}
#'     \item{`test_info`}{Per test, the method `id`, a readable `label` and the
#'       post-hoc procedure that followed it.}
#'     \item{`diagnostics`}{Assumption checks, or `NULL`.}
#'     \item{`metadata`}{`package_version`, `r_version`, `platform` and an
#'       ISO-8601 `timestamp`.}
#'   }
#'
#'   The `effect` table holds `n_used`, `n_groups`, `ref_center` (the centre of
#'   `group_lv[1]`), `extreme_level`, `extreme_center`, `fold_change` and
#'   `log2fc`. `extreme_level` is whichever level sits furthest from the
#'   reference on the log2 scale, and the ratio puts it over the reference, so a
#'   positive `log2fc` means that level is the higher one. Note that this is the
#'   opposite arrangement to [compare_two_groups()], where `group_lv[1]` is the
#'   numerator: with three or more levels the first one reads as a reference
#'   rather than as one side of a contrast. Keeping the column named `log2fc`
#'   means [estimate_significance()] and [draw_volcano_plot()] work on a
#'   multi-group result unchanged.
#'
#' @section Omnibus intervals:
#' The omnibus tables carry `lower_conf` and `upper_conf` as required by the
#' result contract, and both are `NA`. An omnibus test states that the levels
#' are not all alike; it does not state by how much, and there is no single
#' quantity for an interval to be about. The intervals of a multi-group
#' comparison live in `$posthoc`, where each row is one contrast and does have a
#' scale of its own.
#'
#' @section Post-hoc stage:
#' Only features whose omnibus `pval_adj` clears `posthoc_alpha` are compared
#' pairwise, and a feature that did not qualify is absent from the post-hoc
#' table rather than present with `NA`. An absent row means the question was
#' never asked; an `NA` row means it was asked and could not be answered, and
#' the two should not look the same. `parameters$n_posthoc` records how many
#' features entered each stage.
#'
#' Tukey's HSD and Games-Howell control the error rate over the whole set of
#' contrasts through the studentised range, so `pval_adj` equals `pval` for
#' those two and `posthoc_p_adjust` is not applied. Dunn, Conover, pairwise
#' Yuen and pairwise paired t-tests are adjusted across the contrasts within
#' each feature.
#'
#' @section Repeated conditions:
#' A within-subject omnibus test needs a complete rectangle, so `id` is required
#' and subjects missing any condition are dropped whole rather than partially
#' used. The number dropped is reported in `design$unmatched_ids`. Missing values
#' are then handled per feature: a subject with `NA` on one feature is left out
#' of that feature only, which is why `n_used` can differ between features.
#'
#' Row order pairing, which [compare_two_groups()] allows, is deliberately not
#' offered here. With two groups it is at least well defined; with three or more
#' it would also have to assume every condition is stored in the same subject
#' order, and there is no way to notice when it is not.
#'
#' @details
#' Features that cannot be tested do not abort the run. Their row is filled with
#' `NA` and all such features are reported together in a single warning.
#'
#' @seealso [compare_two_groups()] for exactly two levels,
#'   [diagnose_distribution()] for the assumption checks on their own, and
#'   [plot.sa_comparison()] to draw the result.
#'
#' @references
#' Welch, B. L. (1951). On the comparison of several mean values: an alternative
#' approach. *Biometrika*, 38(3-4), 330-336.
#'
#' Kruskal, W. H. and Wallis, W. A. (1952). Use of ranks in one-criterion
#' variance analysis. *Journal of the American Statistical Association*,
#' 47(260), 583-621.
#'
#' Tukey, J. W. (1949). Comparing individual means in the analysis of variance.
#' *Biometrics*, 5(2), 99-114.
#'
#' Games, P. A. and Howell, J. F. (1976). Pairwise multiple comparison
#' procedures with unequal n's and/or variances: a Monte Carlo study.
#' *Journal of Educational Statistics*, 1(2), 113-125.
#'
#' Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics*,
#' 6(3), 241-252.
#'
#' Conover, W. J. (1999). *Practical Nonparametric Statistics*, 3rd edition.
#'
#' @examples
#' ## Independent samples: all three iris species
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#' res <- compare_multiple_groups(
#'   data     = iris,
#'   feats    = feats,
#'   group    = iris$Species,
#'   group_lv = c("setosa", "versicolor", "virginica")
#' )
#' res
#' res$tests$anova_test
#' res$tests$kruskal_test
#'
#' ## The pairwise stage, one row per feature and pair
#' head(res$posthoc$anova_test)
#'
#' ## Setosa is the reference, so a positive log2fc means the most extreme
#' ## other species is the larger one.
#' res$effect
#'
#' ## Repeated conditions: each ChickWeight chick weighed at several times
#' chicks <- ChickWeight[ChickWeight$Time %in% c(0, 6, 12, 18), ]
#' rep_res <- compare_multiple_groups(
#'   data     = data.frame(weight = chicks$weight),
#'   feats    = "weight",
#'   group    = paste0("day", chicks$Time),
#'   group_lv = c("day0", "day6", "day12", "day18"),
#'   id       = chicks$Chick,
#'   paired   = TRUE
#' )
#' rep_res$tests$anova_test   # Mauchly, Greenhouse-Geisser and Huynh-Feldt too
#' rep_res$posthoc$kruskal_test
#'
#' @export
compare_multiple_groups <- function(data,
                                    feats,
                                    group,
                                    group_lv,
                                    id = NULL,
                                    paired = FALSE,
                                    conf_level = 0.95,
                                    tr = 0.2,
                                    posthoc = TRUE,
                                    posthoc_alpha = 0.05,
                                    fc_mean = c("arith", "geom"),
                                    input_scale = c("raw", "log2"),
                                    p_adjust = "BH",
                                    posthoc_p_adjust = "holm",
                                    diagnose = TRUE) {

  input_scale <- match.arg(input_scale)
  fc_mean <- sa_resolve_fc_mean(fc_mean, input_scale, missing(fc_mean))
  sa_check_flag(paired, "paired")
  sa_check_flag(posthoc, "posthoc")
  sa_check_flag(diagnose, "diagnose")
  sa_check_scalar_num(conf_level, "conf_level", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)
  sa_check_scalar_num(tr, "tr", 0, 0.5, upper_open = TRUE)
  sa_check_scalar_num(posthoc_alpha, "posthoc_alpha", 0, 1, lower_open = TRUE)
  sa_check_p_adjust(p_adjust, "p_adjust")
  sa_check_p_adjust(posthoc_p_adjust, "posthoc_p_adjust")

  if (paired && is.null(id)) {
    stop("`paired = TRUE` needs `id` to say which rows belong to the same ",
         "subject. Three or more conditions cannot be matched by row order.",
         call. = FALSE)
  }
  if (!is.null(id) && !paired) {
    warning("`id` is only used to match repeated conditions and is ignored ",
            "when `paired = FALSE`.", call. = FALSE)
  }

  input <- sa_validate_wide_input(data, feats, group, group_lv, id = id,
                                  min_levels = 3L)
  data <- input$data
  feats <- input$feats
  group <- input$group
  id <- input$id
  group_lv <- levels(group)
  n_lv <- length(group_lv)

  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }

  if (paired) {
    aligned <- sa_align_by_subject(id, group, group_lv)
    unmatched <- aligned$unmatched
    if (length(unmatched) > 0L) {
      message("Dropped ", length(unmatched),
              " subject(s) missing at least one condition.")
    }
    # One complete rectangle per feature: a subject with NA anywhere across the
    # conditions is dropped from that feature, since a within-subject test has
    # nothing to compare a partial subject against.
    per_feature <- lapply(feats, function(f) {
      mat <- matrix(data[[f]][aligned$idx], nrow = nrow(aligned$idx),
                    dimnames = list(aligned$subjects, group_lv))
      mat[stats::complete.cases(mat), , drop = FALSE]
    })
  } else {
    unmatched <- character(0)
    per_feature <- lapply(feats, function(f) {
      samples <- lapply(group_lv, function(lv) {
        v <- data[[f]][group == lv]
        v[!is.na(v)]
      })
      stats::setNames(samples, group_lv)
    })
  }
  names(per_feature) <- feats

  effect <- sa_multi_fold_change(per_feature, feats, group_lv, fc_mean, paired,
                                 input_scale)

  specs <- if (paired) {
    sa_multi_specs_repeated(per_feature, conf_level)
  } else {
    sa_multi_specs_independent(per_feature, conf_level, tr)
  }

  tests <- lapply(specs, function(spec) {
    sa_feature_table(feats, spec$columns, spec$label, p_adjust = p_adjust,
                     fun = function(i) spec$omnibus(feats[i]))
  })

  posthoc_tables <- list()
  n_posthoc <- stats::setNames(rep(0L, length(specs)), names(specs))
  if (posthoc) {
    for (nm in names(specs)) {
      spec <- specs[[nm]]
      padj <- tests[[nm]]$pval_adj
      qualified <- feats[!is.na(padj) & padj <= posthoc_alpha]
      n_posthoc[[nm]] <- length(qualified)
      posthoc_tables[[nm]] <- sa_posthoc_table(
        qualified, group_lv, sa_posthoc_columns(), spec$posthoc_label,
        fun = spec$posthoc,
        # A studentised range p-value is already family-wise over the same set
        # of contrasts an adjustment would be applied to, so adjusting it again
        # would correct twice for one comparison.
        p_adjust = if (spec$posthoc_familywise) "none" else posthoc_p_adjust
      )
    }
  }

  diagnostics <- if (diagnose) {
    sa_diagnose_samples(per_feature, feats, group_lv, paired)
  } else {
    NULL
  }

  sa_new_comparison(
    analysis  = "multi_group_comparison",
    features  = feats,
    design    = list(
      group_lv      = group_lv,
      paired        = paired,
      pairing       = if (paired) "id" else NA_character_,
      n_dropped     = input$n_dropped,
      unmatched_ids = unmatched
    ),
    parameters = list(
      alternative      = "two.sided",
      conf_level       = conf_level,
      tr               = if (paired) NA_real_ else tr,
      fc_mean          = fc_mean,
      input_scale      = input_scale,
      p_adjust         = p_adjust,
      posthoc          = posthoc,
      posthoc_alpha    = posthoc_alpha,
      posthoc_p_adjust = posthoc_p_adjust,
      n_posthoc        = n_posthoc
    ),
    effect    = effect,
    tests     = tests,
    posthoc   = posthoc_tables,
    test_info = lapply(specs, function(spec) {
      list(id            = spec$id,
           label         = spec$label,
           paired        = paired,
           posthoc_id    = spec$posthoc_id,
           posthoc_label = spec$posthoc_label)
    }),
    diagnostics = diagnostics,
    subclass  = "sa_multi_group"
  )
}


#' Omnibus and post-hoc pairs for independent levels
#'
#' Each entry carries the omnibus kernel, its column layout and the post-hoc
#' procedure that shares its assumptions. Keeping the two in one object is what
#' makes it impossible to follow a rank-based omnibus test with a parametric
#' comparison by accident.
#'
#' @keywords internal
#' @noRd
sa_multi_specs_independent <- function(per_feature, conf_level, tr) {
  list(
    anova_test = list(
      id = "oneway_anova", label = "One-way ANOVA",
      columns = c("n_used", "n_groups", "f_stat", "df1", "df2", "eta_sq",
                  "omega_sq", "pval", "lower_conf", "upper_conf"),
      omnibus = function(f) sa_oneway_anova(sa_require_groups(per_feature[[f]], 2L)),
      posthoc_id = "tukey_hsd", posthoc_label = "Tukey HSD",
      posthoc_familywise = TRUE,
      posthoc = function(f) sa_tukey(per_feature[[f]], conf_level)
    ),
    welch_test = list(
      id = "welch_anova", label = "Welch's one-way ANOVA",
      columns = c("n_used", "n_groups", "f_stat", "df1", "df2", "eta_sq",
                  "omega_sq", "pval", "lower_conf", "upper_conf"),
      omnibus = function(f) sa_welch_anova(sa_require_groups(per_feature[[f]], 2L)),
      posthoc_id = "games_howell", posthoc_label = "Games-Howell post-hoc test",
      posthoc_familywise = TRUE,
      posthoc = function(f) sa_games_howell(per_feature[[f]], conf_level)
    ),
    robust_test = list(
      id = "yuen_anova", label = "Yuen's trimmed mean one-way ANOVA",
      columns = c("n_used", "n_groups", "f_stat", "df1", "df2",
                  "robust_eta_sq", "pval", "lower_conf", "upper_conf"),
      omnibus = function(f) sa_yuen_anova(sa_require_groups(per_feature[[f]], 2L), tr),
      posthoc_id = "pairwise_yuen", posthoc_label = "Pairwise Yuen tests",
      posthoc_familywise = FALSE,
      posthoc = function(f) sa_pairwise_yuen(per_feature[[f]], tr, conf_level)
    ),
    kruskal_test = list(
      id = "kruskal_wallis", label = "Kruskal-Wallis test",
      columns = c("n_used", "n_groups", "h_stat", "df", "epsilon_sq",
                  "eta_sq_rank", "pval", "lower_conf", "upper_conf"),
      omnibus = function(f) sa_kruskal(sa_require_groups(per_feature[[f]], 1L)),
      posthoc_id = "dunn_test", posthoc_label = "Dunn's post-hoc test",
      posthoc_familywise = FALSE,
      posthoc = function(f) sa_dunn(per_feature[[f]], conf_level)
    )
  )
}


#' Omnibus and post-hoc pairs for repeated conditions
#'
#' @keywords internal
#' @noRd
sa_multi_specs_repeated <- function(per_feature, conf_level) {
  list(
    anova_test = list(
      id = "repeated_measures_anova", label = "Repeated measures ANOVA",
      columns = c("n_used", "n_groups", "f_stat", "df1", "df2",
                  "partial_eta_sq", "gen_eta_sq", "mauchly_w", "mauchly_pval",
                  "gg_eps", "pval_gg", "hf_eps", "pval_hf", "pval",
                  "lower_conf", "upper_conf"),
      omnibus = function(f) sa_rm_anova(per_feature[[f]]),
      posthoc_id = "pairwise_paired_t",
      posthoc_label = "Pairwise paired t-tests",
      posthoc_familywise = FALSE,
      posthoc = function(f) sa_pairwise_paired_t(per_feature[[f]], conf_level)
    ),
    kruskal_test = list(
      id = "friedman_test", label = "Friedman test",
      columns = c("n_used", "n_groups", "chi_sq", "df", "kendalls_w", "pval",
                  "lower_conf", "upper_conf"),
      omnibus = function(f) sa_friedman(per_feature[[f]]),
      posthoc_id = "conover_posthoc", posthoc_label = "Conover post-hoc test",
      posthoc_familywise = FALSE,
      posthoc = function(f) sa_conover(per_feature[[f]], conf_level)
    )
  )
}


#' Reject a feature whose groups are too small before the kernel sees it
#'
#' The kernels raise their own errors, but they see the samples without knowing
#' which level each came from. Checking here means the message names the level.
#'
#' @keywords internal
#' @noRd
sa_require_groups <- function(samples, n_min) {
  sizes <- lengths(samples)
  short <- names(samples)[sizes < n_min]
  if (length(short) > 0L) {
    stop("needs at least ", n_min, " usable observation(s) per group; ",
         paste0(short, " = ", sizes[short], collapse = ", "), ".",
         call. = FALSE)
  }
  samples
}


#' Fold change of the most extreme level against the reference
#'
#' The multi-group counterpart of `sa_fold_change()`. Both reduce a comparison
#' to one signed magnitude per feature so that the same volcano plot path works
#' for either, but with three or more levels there is no single contrast to take
#' the ratio of. The level furthest from the reference on the log2 scale is used,
#' which is the largest change the comparison found.
#'
#' @param per_feature Per feature, either a named list of samples (independent)
#'   or a subjects-by-conditions matrix (repeated).
#' @param feats Feature names, one output row per entry.
#' @param group_lv Group levels, the first being the reference denominator.
#' @param mean_type `"arith"` or `"geom"`.
#' @param paired Whether `per_feature` holds matrices rather than sample lists.
#' @param input_scale `"raw"` or `"log2"`.
#'
#' @return data.frame with `features`, `n_used`, `n_groups`, `ref_center`,
#'   `extreme_level`, `extreme_center`, `fold_change` and `log2fc`.
#'
#' @keywords internal
#' @noRd
sa_multi_fold_change <- function(per_feature, feats, group_lv, mean_type,
                                 paired, input_scale = "raw") {
  label <- paste0(if (mean_type == "arith") "Arithmetic" else "Geometric",
                  " mean fold change")

  out <- sa_feature_table(
    feats,
    c("n_used", "n_groups", "ref_center", "extreme_index", "extreme_center",
      "fold_change", "log2fc"),
    label,
    p_adjust = NULL,
    fun = function(i) {
      samples <- if (paired) {
        mat <- per_feature[[i]]
        stats::setNames(lapply(seq_len(ncol(mat)), function(j) mat[, j]),
                        group_lv)
      } else {
        per_feature[[i]]
      }
      centers <- vapply(seq_along(samples), function(j) {
        sa_fc_center(samples[[j]], group_lv[j], mean_type, input_scale)
      }, numeric(1))

      ratios <- centers / centers[1]
      log_ratios <- suppressWarnings(log2(ratios))
      # A level whose ratio left the domain of log2() cannot be ranked by
      # distance, so it is skipped rather than allowed to win by being NaN.
      rankable <- is.finite(log_ratios)
      rankable[1] <- FALSE
      extreme <- if (any(rankable)) {
        which(rankable)[which.max(abs(log_ratios[rankable]))]
      } else {
        2L
      }

      c(n_used         = if (paired) nrow(per_feature[[i]]) else
                           sum(lengths(samples)),
        n_groups       = length(samples),
        ref_center     = centers[1],
        extreme_index  = extreme,
        extreme_center = centers[extreme],
        fold_change    = ratios[extreme],
        log2fc         = log_ratios[extreme])
    }
  )

  out$extreme_level <- group_lv[out$extreme_index]
  out[c("features", "n_used", "n_groups", "ref_center", "extreme_level",
        "extreme_center", "fold_change", "log2fc")]
}
