#' Flag candidate outliers without removing them
#'
#' Screens each feature, or each feature within each group level, and returns
#' one row per flagged observation. Nothing is deleted and no analysis changes
#' as a result. Which observations belong in a data set is a decision about the
#' experiment rather than about the arithmetic, and the package does not make
#' it on the user's behalf.
#'
#' Three rules are available and they do not agree with each other, which is the
#' point of naming the one used:
#'
#' \describe{
#'   \item{`"iqr"`}{Anything past `Q1 - k * IQR` or `Q3 + k * IQR`. The rule
#'     behind the whiskers of a boxplot. It makes no distributional assumption
#'     and flags a fixed share of any long-tailed sample.}
#'   \item{`"robust_z"`}{Distance from the median in units of the median
#'     absolute deviation. Uses the median and MAD rather than the mean and
#'     standard deviation because one extreme value inflates the standard
#'     deviation enough to hide itself.}
#'   \item{`"grubbs"`}{Tests only the single most extreme observation, against
#'     the null that the sample is normal. It is the only rule that produces a
#'     p-value and the only one that assumes a distribution.}
#' }
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to screen.
#' @param group Optional grouping vector with one entry per row of `data`. When
#'   supplied, each feature is screened within each level separately, which is
#'   what keeps a genuine group difference from being read as a set of outliers.
#' @param group_lv Group levels to keep, in display order. Defaults to the
#'   sorted unique values of `group`. Rows belonging to any other level are
#'   dropped.
#' @param criterion `"iqr"`, `"robust_z"` or `"grubbs"`.
#' @param iqr_multiplier Fence width for `criterion = "iqr"`.
#' @param z_threshold Cut-off for `criterion = "robust_z"`. The 3.5 default is
#'   the Iglewicz and Hoaglin recommendation.
#' @param alpha Significance level for `criterion = "grubbs"`.
#'
#' @return A data.frame with one row per flagged observation and the columns
#'   `features`, `group` (`NA` when no `group` was given), `row` (the row number
#'   in the `data` that was passed in), `value` and `score`. The score is the
#'   quantity the rule thresholded: distance past the nearer quartile in IQR
#'   units for `"iqr"`, the robust z for `"robust_z"`, the Grubbs statistic for
#'   `"grubbs"`. Zero rows means nothing was flagged. The criterion and its
#'   thresholds are attached as attributes.
#'
#' @seealso [diagnose_distribution()], which runs this alongside the normality
#'   and variance checks, and [summarize_descriptive_stats()], whose
#'   `out_lower_bound` and `out_upper_bound` columns are the same IQR fences.
#'
#' @references
#' Iglewicz, B. and Hoaglin, D. C. (1993). *How to Detect and Handle Outliers*.
#' ASQC Quality Press.
#'
#' Grubbs, F. E. (1969). Procedures for detecting outlying observations in
#' samples. *Technometrics*, 11(1), 1-21.
#'
#' @examples
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#'
#' ## Across the whole data set, sepal width is the only measurement with
#' ## observations past the fences.
#' screen_outliers(iris, feats)
#'
#' ## Within species, petal measurements pick up flags that the pooled screen
#' ## hid: the between-species spread was widening the fences.
#' screen_outliers(iris, feats, iris$Species)
#'
#' ## A stricter rule and a distributional one disagree, which is why the
#' ## criterion is named rather than assumed.
#' nrow(screen_outliers(iris, feats, criterion = "robust_z"))
#' screen_outliers(iris, feats, criterion = "grubbs")
#'
#' @export
screen_outliers <- function(data,
                            feats,
                            group = NULL,
                            group_lv = NULL,
                            criterion = c("iqr", "robust_z", "grubbs"),
                            iqr_multiplier = 1.5,
                            z_threshold = 3.5,
                            alpha = 0.05) {

  criterion <- match.arg(criterion)
  sa_check_scalar_num(iqr_multiplier, "iqr_multiplier", 0)
  sa_check_scalar_num(z_threshold, "z_threshold", 0, lower_open = TRUE)
  sa_check_scalar_num(alpha, "alpha", 0, 1, lower_open = TRUE)

  split <- sa_split_for_screening(data, feats, group, group_lv)

  blocks <- lapply(feats, function(f) {
    per_group <- lapply(names(split$rows), function(lv) {
      rows <- split$rows[[lv]]
      v <- split$data[[f]][rows]
      res <- sa_flag_outliers(v, criterion, iqr_multiplier, z_threshold, alpha)
      hit <- which(res$flag)
      data.frame(
        features = rep(f, length(hit)),
        group    = rep(if (split$grouped) lv else NA_character_, length(hit)),
        row      = split$row_id[rows[hit]],
        value    = v[hit],
        score    = res$score[hit],
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, per_group)
  })

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL

  structure(out,
            criterion      = criterion,
            iqr_multiplier = iqr_multiplier,
            z_threshold    = z_threshold,
            alpha          = alpha)
}


#' Row indices to screen, one entry per group level or one for everything
#'
#' Shared by [screen_outliers()] and [diagnose_distribution()] so that an
#' ungrouped call and a grouped one differ in exactly one place.
#'
#' @return List with the validated `data`, a named list of row index vectors in
#'   `rows`, `row_id` translating those back into rows of the `data` the caller
#'   passed in, and `grouped` saying whether `group` was supplied.
#'
#' @keywords internal
#' @noRd
sa_split_for_screening <- function(data, feats, group, group_lv) {
  if (is.null(group)) {
    if (is.matrix(data)) {
      data <- as.data.frame(data)
    }
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame or a matrix.", call. = FALSE)
    }
    if (nrow(data) == 0L) {
      stop("`data` has zero rows.", call. = FALSE)
    }
    sa_check_feat_names(feats)
    unknown <- setdiff(feats, names(data))
    if (length(unknown) > 0L) {
      stop("`feats` not found in `data`: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    non_numeric <- feats[!vapply(data[feats], is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
      stop("`feats` must refer to numeric columns. Not numeric: ",
           paste(non_numeric, collapse = ", "), call. = FALSE)
    }
    return(list(data = data, rows = list(all = seq_len(nrow(data))),
                row_id = seq_len(nrow(data)), grouped = FALSE))
  }

  if (is.null(group_lv)) {
    group_lv <- sort(unique(as.character(group)))
  }
  # min_levels = 1 because a single level is a legitimate thing to screen; it is
  # only the comparison functions that need two or more.
  input <- sa_validate_wide_input(data, feats, group, group_lv, min_levels = 1L)
  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }
  rows <- lapply(levels(input$group), function(lv) which(input$group == lv))
  names(rows) <- levels(input$group)

  # Rows outside `group_lv` were dropped, so a position in the filtered data no
  # longer matches the row the caller can look up. `row_id` translates back.
  list(data = input$data, rows = rows,
       row_id = which(as.character(group) %in% as.character(group_lv)),
       grouped = TRUE)
}
