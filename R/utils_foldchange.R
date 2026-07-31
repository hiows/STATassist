# Fold change between the two group levels of a comparison.
#
# This was `calculate_fold_change()`, an exported function taking `data` and
# deriving its own samples. It now works from the samples the tests were run on,
# which is what keeps the two axes of a volcano plot describing the same
# observations: under `paired = TRUE` the tests keep complete pairs only, and a
# separately computed fold change would quietly average a different set of rows.
#
# The direction is not an argument either. It follows `group_lv`, the same order
# the tests read as `group_lv[1] - group_lv[2]`, so the x and y axes of a volcano
# plot cannot end up pointing opposite ways.

#' Central tendency of one group for the fold change ratio
#'
#' @param v Numeric vector with missing values already removed.
#' @param side Group label, used in the error message.
#' @param mean_type `"arith"` or `"geom"`.
#'
#' @keywords internal
#' @noRd
sa_fc_center <- function(v, side, mean_type) {
  if (length(v) == 0L) {
    stop("no usable observation left in the ", side, " group.", call. = FALSE)
  }
  if (mean_type == "arith") {
    return(mean(v))
  }
  # Dropping the non-positive values instead would silently return the geometric
  # mean of the positive subset, which is a different quantity.
  if (any(v <= 0)) {
    stop("the geometric mean is undefined for the ", sum(v <= 0),
         " value(s) at or below zero in the ", side,
         " group; use `fc_mean = \"arith\"` instead.", call. = FALSE)
  }
  exp(mean(log(v)))
}


#' Fold change table for a two-group comparison
#'
#' @param samples List of `list(x = , y = )` per feature, already reduced to the
#'   observations the tests used.
#' @param feats Feature names, one output row per entry.
#' @param group_lv The two group levels, `group_lv[1]` going in the numerator.
#' @param mean_type `"arith"` or `"geom"`.
#'
#' @return data.frame with `features`, `x_center`, `y_center`, `fold_change` and
#'   `log2fc`.
#'
#' @keywords internal
#' @noRd
sa_fold_change <- function(samples, feats, group_lv, mean_type) {
  label <- paste0(if (mean_type == "arith") "Arithmetic" else "Geometric",
                  " mean fold change")

  out <- sa_feature_table(
    feats,
    c("x_center", "y_center", "fold_change", "log2fc"),
    label,
    p_adjust = NULL,
    fun = function(i) {
      x_center <- sa_fc_center(samples[[i]]$x, group_lv[1], mean_type)
      y_center <- sa_fc_center(samples[[i]]$y, group_lv[2], mean_type)
      fold_change <- x_center / y_center
      c(x_center    = x_center,
        y_center    = y_center,
        fold_change = fold_change,
        # log2() of a negative number warns once per call, which says nothing
        # about which feature caused it. The features are reported below instead.
        log2fc      = suppressWarnings(log2(fold_change)))
    }
  )

  report <- function(mask, what) {
    hit <- out$features[!is.na(mask) & mask]
    if (length(hit) > 0L) {
      message("Fold change: ", what, " for ", length(hit), " feature(s): ",
              paste(hit, collapse = ", "), ".")
    }
  }
  # A ratio only reads as "n times higher" when both centres are positive, so
  # each way of leaving that domain is called out separately.
  report(out$y_center == 0,
         paste0("the ", group_lv[2],
                " centre is zero, so `fold_change` is infinite"))
  report(out$x_center == 0 & out$y_center != 0,
         paste0("the ", group_lv[1],
                " centre is zero, so `log2fc` is -Inf and clears any cutoff"))
  report(out$x_center * out$y_center < 0,
         "the two centres have opposite signs, so `log2fc` is NaN")

  out
}
