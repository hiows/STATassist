#' Fold change between two groups, feature by feature
#'
#' Divides the case group mean by the control group mean for every feature, the
#' quantity a volcano plot puts on its x axis once logged.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data`. One
#'   returned value per entry, in this order.
#' @param group Grouping vector with one entry per row of `data`.
#' @param group_lv Character vector of exactly two group levels. Rows belonging
#'   to any other level are dropped.
#' @param case_label The level of `group_lv` that goes in the numerator. The
#'   remaining level becomes the denominator, so a value above 1 means higher in
#'   `case_label`.
#' @param mean_type Which mean to divide, `"arith"` for the arithmetic mean or
#'   `"geom"` for the geometric mean. The geometric mean requires strictly
#'   positive values.
#'
#' @return Numeric vector of fold changes named by `feats`.
#'
#' @details
#' Missing values are dropped within each group, so the two groups may rest on
#' different numbers of observations for different features.
#'
#' A feature that cannot be summarised does not abort the run. Its entry is
#' `NA` and all such features are reported together in a single warning. This
#' covers a group left empty after dropping `NA` and, under
#' `mean_type = "geom"`, any value at or below zero: `log()` turns those into
#' `NaN` or `-Inf`, and dropping them would quietly return the geometric mean of
#' the positive subset instead.
#'
#' A fold change is a ratio, so it only reads as "n times higher" when both
#' means are positive. Features whose control mean is zero (giving `Inf`) or
#' whose means differ in sign are therefore reported in a `message()`.
#'
#' @seealso [evaluate_significance()] to combine the result with p-values, and
#'   [compare_two_groups()] to obtain those p-values.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#'
#' ## virginica relative to versicolor
#' calculate_fold_change(
#'   data       = iris2,
#'   feats      = feats,
#'   group      = iris2$Species,
#'   group_lv   = c("versicolor", "virginica"),
#'   case_label = "virginica"
#' )
#'
#' ## The default numerator is the first level, which inverts every ratio.
#' calculate_fold_change(iris2, feats, iris2$Species,
#'                       c("versicolor", "virginica"))
#'
#' ## Geometric means are the usual choice for concentration-like data.
#' calculate_fold_change(iris2, feats, iris2$Species,
#'                       c("versicolor", "virginica"),
#'                       mean_type = "geom")
#'
#' @export
calculate_fold_change <- function(data,
                                  feats,
                                  group,
                                  group_lv,
                                  case_label = group_lv[1],
                                  mean_type = c("arith", "geom")) {

  mean_type <- match.arg(mean_type)

  input <- sa_validate_wide_input(data, feats, group, group_lv, n_levels = 2L)
  data <- input$data
  feats <- input$feats
  group <- input$group
  group_lv <- levels(group)

  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }

  if (length(case_label) != 1L || is.na(case_label) ||
      !as.character(case_label) %in% group_lv) {
    stop("`case_label` must be one of the levels in `group_lv` (",
         paste(group_lv, collapse = ", "), "), but is ",
         paste(deparse(case_label), collapse = " "), ".", call. = FALSE)
  }
  case_label <- as.character(case_label)
  control_label <- setdiff(group_lv, case_label)

  # The samples come straight out of the wide columns. Reshaping to long format
  # first and reading a `grp` column back was what made every feature NaN: no
  # such column was ever created, so both groups came back empty.
  idx_case <- which(group == case_label)
  idx_control <- which(group == control_label)

  group_mean <- function(v, side) {
    v <- v[!is.na(v)]
    if (length(v) == 0L) {
      stop("no non-missing observation left in the ", side, " group.",
           call. = FALSE)
    }
    if (mean_type == "arith") {
      return(mean(v))
    }
    if (any(v <= 0)) {
      stop("the geometric mean is undefined for the ", sum(v <= 0),
           " value(s) at or below zero in the ", side,
           " group; use `mean_type = \"arith\"` instead.", call. = FALSE)
    }
    exp(mean(log(v)))
  }

  means <- matrix(NA_real_, nrow = length(feats), ncol = 2L,
                  dimnames = list(feats, c("case", "control")))

  label <- paste0(if (mean_type == "arith") "Arithmetic" else "Geometric",
                  " mean fold change")

  fold_change <- sa_feature_scalar(feats, label, fun = function(i) {
    f <- feats[i]
    case <- group_mean(data[[f]][idx_case], case_label)
    control <- group_mean(data[[f]][idx_control], control_label)
    means[i, ] <<- c(case, control)
    case / control
  })

  zero_control <- feats[!is.na(means[, "control"]) & means[, "control"] == 0]
  if (length(zero_control) > 0L) {
    message("The ", control_label, " mean is zero for ", length(zero_control),
            " feature(s), so the ratio is infinite: ",
            paste(zero_control, collapse = ", "), ".")
  }
  sign_mismatch <- feats[
    !is.na(means[, "case"]) & !is.na(means[, "control"]) &
      means[, "case"] * means[, "control"] < 0
  ]
  if (length(sign_mismatch) > 0L) {
    message("The two group means have opposite signs for ",
            length(sign_mismatch),
            " feature(s), so the ratio does not read as a fold change: ",
            paste(sign_mismatch, collapse = ", "), ".")
  }

  fold_change
}
