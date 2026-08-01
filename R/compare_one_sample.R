#' Compare one sample against a hypothesised value
#'
#' Tests each feature against `mu` and returns a parametric, a rank-based and a
#' proportion result side by side, in the same shape [compare_two_groups()]
#' uses. There is no second group here, so the reference is a number the user
#' supplies rather than a set of observations.
#'
#' \describe{
#'   \item{`t_test`}{One-sample t-test on the mean.}
#'   \item{`wilcox_test`}{One-sample Wilcoxon signed-rank test, with the
#'     Hodges-Lehmann pseudo-median. It tests location without assuming
#'     normality, but it does assume the distribution of `x - mu` is symmetric.}
#'   \item{`prop_test`}{Score test of the proportion of successes against `p`.
#'     It only applies to features that are binary, and features that are not
#'     come back as `NA` rather than being silently coerced.}
#' }
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to test.
#'   One output row per entry.
#' @param mu Hypothesised value for the mean and the pseudo-median. A single
#'   number applied to every feature.
#' @param p Hypothesised proportion for `prop_test`, in `(0, 1)`.
#' @param success The value counted as a success when a feature is binary.
#'   Defaults to 1, so a 0/1 coded column needs nothing further.
#' @param alternative Direction of the alternative hypothesis, one of
#'   `"two.sided"`, `"less"` or `"greater"`. `"greater"` tests whether the
#'   sample exceeds `mu`, and every reported quantity follows that direction.
#' @param conf_level Confidence level for all reported intervals.
#' @param p_adjust Multiplicity adjustment applied across `feats` within each
#'   test table, passed to [stats::p.adjust()]. Use `"none"` to disable.
#' @param diagnose Logical. If `TRUE`, the normality check the t-test rests on
#'   is attached as `$diagnostics`.
#'
#' @return A `sa_comparison` object with the layout described in
#'   [compare_two_groups()], with two differences. `design` carries `mu`,
#'   `p` and `success` instead of `group_lv`, since there are no groups, and
#'   `effect` holds `n_used`, `center` (the sample mean), `mu`, `diff`,
#'   `fold_change` and `log2fc`, the ratio being the sample mean over `mu`.
#'   `posthoc` is empty.
#'
#'   Every test table starts with `n_used` and carries `pval`, `pval_adj`,
#'   `lower_conf` and `upper_conf`. The remaining columns are:
#'
#'   \describe{
#'     \item{`t_test`}{`center`, `mu`, `diff`, `stderr`, `t_stat`, `df` and
#'       `cohens_d`, the mean difference over the sample standard deviation.}
#'     \item{`wilcox_test`}{`hl_shift`, the Hodges-Lehmann pseudo-median, and
#'       `v_stat`.}
#'     \item{`prop_test`}{`n_success`, `proportion`, `p`, `diff`, `chi_sq`,
#'       `df` and `cohens_h`. The interval is a Wilson score interval, which
#'       stays inside `[0, 1]` where a Wald interval does not.}
#'   }
#'
#' @section Fold change against a hypothesised value:
#' `fold_change` is `center / mu` and `log2fc` its base-2 logarithm, so both are
#' undefined when `mu` is zero, which is also its most common value. Both
#' columns are `NA` in that case rather than infinite, a message says so, and
#' [estimate_significance()] will call every feature undecided. Reporting `Inf`
#' would read as an infinitely large increase when what actually happened is
#' that the question has no answer.
#'
#' @details
#' Missing values are dropped per feature, so `n_used` varies between features
#' when the data are incomplete. Features that cannot be tested do not abort the
#' run: their row is filled with `NA` and all such features are reported
#' together in a single warning.
#'
#' @seealso [compare_two_groups()] for two groups and
#'   [compare_multiple_groups()] for three or more.
#'
#' @references
#' Student (1908). The probable error of a mean. *Biometrika*, 6(1), 1-25.
#'
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods.
#' *Biometrics Bulletin*, 1(6), 80-83.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#' statistical inference. *Journal of the American Statistical Association*,
#' 22(158), 209-212.
#'
#' @examples
#' ## Are the iris measurements different from 3 cm?
#' res <- compare_one_sample(
#'   data  = iris,
#'   feats = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   mu    = 3
#' )
#' res
#' res$tests$t_test
#' res$tests$wilcox_test
#' res$effect
#'
#' ## A binary feature reaches the proportion test; a continuous one does not
#' ## and comes back NA with a warning rather than being coerced.
#' cars <- mtcars
#' cars$am <- as.numeric(cars$am)
#' suppressWarnings(
#'   compare_one_sample(cars, c("am", "vs", "mpg"), mu = 0.5, p = 0.5)
#' )$tests$prop_test
#'
#' ## One-sided: is mileage above 20 mpg?
#' compare_one_sample(mtcars, "mpg", mu = 20,
#'                    alternative = "greater")$tests$t_test
#'
#' @export
compare_one_sample <- function(data,
                               feats,
                               mu = 0,
                               p = 0.5,
                               success = 1,
                               alternative = c("two.sided", "less", "greater"),
                               conf_level = 0.95,
                               p_adjust = "BH",
                               diagnose = TRUE) {

  alternative <- match.arg(alternative)
  sa_check_scalar_num(mu, "mu")
  sa_check_scalar_num(p, "p", 0, 1, lower_open = TRUE, upper_open = TRUE)
  sa_check_scalar_num(conf_level, "conf_level", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)
  sa_check_p_adjust(p_adjust, "p_adjust")
  sa_check_flag(diagnose, "diagnose")
  if (!is.numeric(success) || length(success) != 1L || is.na(success)) {
    stop("`success` must be a single non-missing number.", call. = FALSE)
  }

  # One synthetic level, so the same validation that guards the comparison
  # functions guards this one: numeric columns, matching lengths, no duplicates.
  split <- sa_split_for_screening(data, feats, group = NULL, group_lv = NULL)
  data <- split$data

  samples <- lapply(feats, function(f) {
    v <- data[[f]]
    v[!is.na(v)]
  })
  names(samples) <- feats

  effect <- sa_feature_table(
    feats, c("n_used", "center", "mu", "diff", "fold_change", "log2fc"),
    "Fold change against mu", p_adjust = NULL,
    fun = function(i) {
      v <- samples[[i]]
      if (length(v) == 0L) {
        stop("no usable observation left.", call. = FALSE)
      }
      center <- mean(v)
      # mu = 0 is both the usual default and the one value a ratio cannot be
      # taken against. Reporting Inf would suggest an infinitely large increase
      # when what happened is that the question has no answer.
      ratio <- if (mu == 0) NA_real_ else center / mu
      c(n_used      = length(v),
        center      = center,
        mu          = mu,
        diff        = center - mu,
        fold_change = ratio,
        log2fc      = suppressWarnings(log2(ratio)))
    }
  )
  if (mu == 0) {
    message("`mu` is 0, so `fold_change` and `log2fc` are undefined and the ",
            "`effect` table reports them as NA.")
  }

  t_columns <- c("n_used", "center", "mu", "diff", "stderr", "t_stat", "df",
                 "cohens_d", "pval", "lower_conf", "upper_conf")
  t_result <- sa_feature_table(feats, t_columns, "One-sample t-test",
    p_adjust = p_adjust,
    fun = function(i) {
      v <- samples[[i]]
      if (length(v) < 2L) {
        stop("needs at least 2 usable observations, got ", length(v), ".",
             call. = FALSE)
      }
      res <- stats::t.test(v, mu = mu, alternative = alternative,
                           conf.level = conf_level)
      spread <- stats::sd(v)
      sa_row(n_used     = length(v),
             center     = mean(v),
             mu         = mu,
             diff       = mean(v) - mu,
             stderr     = res$stderr,
             t_stat     = res$statistic,
             df         = res$parameter,
             cohens_d   = if (spread > 0) (mean(v) - mu) / spread else NA_real_,
             pval       = res$p.value,
             lower_conf = res$conf.int[1],
             upper_conf = res$conf.int[2])
    })

  w_columns <- c("n_used", "hl_shift", "v_stat", "pval", "lower_conf",
                 "upper_conf")
  w_result <- sa_feature_table(feats, w_columns,
    "One-sample Wilcoxon signed-rank test", p_adjust = p_adjust,
    fun = function(i) {
      v <- samples[[i]]
      if (length(v) < 1L) {
        stop("needs at least 1 usable observation.", call. = FALSE)
      }
      res <- stats::wilcox.test(v, mu = mu, alternative = alternative,
                                conf.int = TRUE, conf.level = conf_level)
      sa_row(n_used     = length(v),
             hl_shift   = res$estimate,
             v_stat     = res$statistic,
             pval       = res$p.value,
             lower_conf = res$conf.int[1],
             upper_conf = res$conf.int[2])
    })

  prop_columns <- c("n_used", "n_success", "proportion", "p", "diff", "chi_sq",
                    "df", "cohens_h", "pval", "lower_conf", "upper_conf")
  prop_result <- sa_feature_table(feats, prop_columns,
    "One-sample proportion test", p_adjust = p_adjust,
    fun = function(i) sa_one_sample_prop(samples[[i]], p, success, alternative,
                                         conf_level))

  sa_new_comparison(
    analysis  = "one_sample_comparison",
    features  = feats,
    design    = list(
      mu        = mu,
      p         = p,
      success   = success,
      paired    = FALSE,
      n_dropped = 0L
    ),
    parameters = list(
      alternative = alternative,
      conf_level  = conf_level,
      p_adjust    = p_adjust
    ),
    effect    = effect,
    tests     = list(
      t_test      = t_result,
      wilcox_test = w_result,
      prop_test   = prop_result
    ),
    test_info = list(
      t_test = list(id = "one_sample_t_test", label = "One-sample t-test",
                    paired = FALSE),
      wilcox_test = list(id = "one_sample_wilcoxon",
                         label = "One-sample Wilcoxon signed-rank test",
                         paired = FALSE),
      prop_test = list(id = "one_sample_proportion",
                       label = "One-sample proportion test", paired = FALSE)
    ),
    diagnostics = if (diagnose) {
      sa_diagnose_samples(lapply(samples, function(v) list(sample = v)),
                          feats, "sample", paired = FALSE)
    } else {
      NULL
    },
    subclass  = "sa_one_sample"
  )
}


#' Score test and Wilson interval for one proportion
#'
#' A feature has to be binary for this to mean anything, and a continuous column
#' silently reduced to "equals `success` or not" would produce a number that
#' looks like a result. It is rejected instead, which turns into an `NA` row and
#' a named warning through `sa_feature_table()`.
#'
#' The interval is Wilson's rather than Wald's. A Wald interval on a proportion
#' near 0 or 1 runs outside `[0, 1]`, which is not a statement the data can
#' make.
#'
#' @param v Numeric vector without missing values.
#' @param p Hypothesised proportion.
#' @param success Value counted as a success.
#' @param alternative Direction of the alternative hypothesis.
#' @param conf_level Confidence level of the interval.
#'
#' @keywords internal
#' @noRd
sa_one_sample_prop <- function(v, p, success, alternative, conf_level) {
  observed <- unique(v)
  if (length(observed) > 2L) {
    stop("the proportion test needs a binary feature, but this one takes ",
         length(observed), " distinct values.", call. = FALSE)
  }
  if (!success %in% observed) {
    stop("the value counted as a success, ", success,
         ", does not occur in this feature.", call. = FALSE)
  }

  n <- length(v)
  n_success <- sum(v == success)
  res <- stats::prop.test(n_success, n, p = p, alternative = alternative,
                          conf.level = conf_level)
  proportion <- n_success / n

  sa_row(n_used     = n,
         n_success  = n_success,
         proportion = proportion,
         p          = p,
         diff       = proportion - p,
         chi_sq     = res$statistic,
         df         = res$parameter,
         # Cohen's h compares proportions on the arcsine scale, where a given
         # difference means the same thing near 0.5 and near the boundaries.
         cohens_h   = 2 * asin(sqrt(proportion)) - 2 * asin(sqrt(p)),
         pval       = res$p.value,
         lower_conf = res$conf.int[1],
         upper_conf = res$conf.int[2])
}
