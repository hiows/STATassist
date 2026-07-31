#' Descriptive summary of several features, optionally split by group
#'
#' Reduces every feature to one row of sample size, central tendency,
#' dispersion, quartiles, outlier fences and distribution shape. With a
#' grouping vector the same row is produced per group level, so the summary
#' lines up with the tests and plots that compare those levels.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to
#'   summarise, in output order.
#' @param group Optional grouping vector with one entry per row of `data`. When
#'   `NULL`, all rows are summarised together and no `group` column is
#'   returned.
#' @param group_lv Group levels to report, in output order. When `NULL`, the
#'   levels present in `group` are used: the factor levels if `group` is a
#'   factor, the sorted unique values otherwise. Rows belonging to any other
#'   level are dropped.
#'
#' @return A data.frame with one row per feature, or per feature and group
#'   level when `group` is supplied. The leading columns are `features` and,
#'   when grouped, `group`. The remaining columns are:
#'
#'   \describe{
#'     \item{`n`, `n_missing`}{Finite values the row is based on, and the
#'       values left out for being `NA`, `NaN` or infinite.}
#'     \item{`mean`, `sd`, `var`, `se`, `cv`}{Mean, standard deviation,
#'       variance, standard error `sd / sqrt(n)` and coefficient of variation
#'       `sd / mean`.}
#'     \item{`min`, `q1`, `median`, `q3`, `max`, `iqr`}{Five-number summary and
#'       the interquartile range. Quartiles come from [stats::quantile()] with
#'       its default type 7.}
#'     \item{`out_lower_bound`, `out_upper_bound`}{Tukey outlier fences
#'       `q1 - 1.5 * iqr` and `q3 + 1.5 * iqr`.}
#'     \item{`mad`}{Median absolute deviation from [stats::mad()], scaled by
#'       the default 1.4826 so that it estimates `sd` for normal data.}
#'     \item{`skewness`, `excess_kurtosis`}{Shape of the distribution, zero for
#'       a normal sample.}
#'   }
#'
#' @details
#' Missing and non-finite values are dropped per feature and per group before
#' anything is computed, so one `Inf` cannot turn a whole row into `Inf` or
#' `NaN`; `n_missing` records how many were left out. A feature with no finite
#' value in a group gives an all-`NA` row rather than aborting the summary.
#'
#' `skewness` and `excess_kurtosis` are the bias-corrected G1 and G2
#' estimators used by SAS and SPSS, the same quantities as
#' `e1071::skewness(type = 2)` and `e1071::kurtosis(type = 2)`. They need three
#' and four observations respectively, and a non-zero spread, and are `NA`
#' otherwise.
#'
#' `cv` is a ratio, so it only reads as relative dispersion when the values are
#' positive. On data that crosses zero the mean shrinks towards it and the
#' ratio explodes without the spread having changed.
#'
#' `out_lower_bound` and `out_upper_bound` are the same fences that
#' [draw_grouped_boxplot()] returns as `lower_bound` and `upper_bound` in
#' `box_summary_stats`. They are where the whiskers may reach, not where they
#' actually end.
#'
#' @seealso [draw_grouped_boxplot()] to see the same quantities as a plot, and
#'   [compare_two_groups()] to test the difference between two levels.
#'
#' @references
#' Joanes, D. N. and Gill, C. A. (1998). Comparing measures of sample skewness
#' and kurtosis. *Journal of the Royal Statistical Society: Series D*, 47(1),
#' 183-189.
#'
#' @examples
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#'
#' ## All rows together
#' summarize_descriptive_stats(iris, feats)
#'
#' ## One row per feature and species
#' by_species <- summarize_descriptive_stats(iris, feats, iris$Species)
#' by_species[by_species$features == "Petal.Length", ]
#'
#' ## Two of the three species, in a chosen order
#' summarize_descriptive_stats(iris, "Petal.Length", iris$Species,
#'                             c("virginica", "versicolor"))
#'
#' @export
summarize_descriptive_stats <- function(data,
                                        feats,
                                        group = NULL,
                                        group_lv = NULL) {

  grouped <- !is.null(group)

  # sa_validate_wide_input() always works in terms of levels, so an ungrouped
  # call is served by a single synthetic level that is dropped again on the way
  # out.
  if (!grouped) {
    if (is.matrix(data)) {
      data <- as.data.frame(data)
    }
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame or a matrix.", call. = FALSE)
    }
    group <- rep("all", nrow(data))
    group_lv <- "all"
  } else if (is.null(group_lv)) {
    group_lv <- if (is.factor(group)) {
      levels(droplevels(group))
    } else {
      sort(unique(as.character(group)))
    }
  }

  input <- sa_validate_wide_input(data, feats, group, group_lv,
                                  min_levels = 1L)
  data <- input$data
  feats <- input$feats
  group <- input$group
  group_lv <- levels(group)

  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }

  # Features vary slowest so that the levels of one feature stay together.
  rows <- unlist(
    lapply(feats, function(f) {
      lapply(group_lv, function(lv) sa_describe_vector(data[[f]][group == lv]))
    }),
    recursive = FALSE
  )

  out <- data.frame(
    features = rep(feats, each = length(group_lv)),
    stringsAsFactors = FALSE
  )
  if (grouped) {
    out$group <- rep(group_lv, times = length(feats))
  }
  out <- cbind(out, as.data.frame(do.call(rbind, rows)))
  rownames(out) <- NULL

  out
}
