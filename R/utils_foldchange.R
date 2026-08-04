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
#
# `input_scale` is the one thing here the tests do not share. They run on the
# values as supplied, which is the point of logging them, while these centres
# are always brought back to the original scale so that a ratio is a ratio. The
# two therefore disagree numerically under `input_scale = "log2"`, and that is
# intended rather than a leak.

#' Resolve `fc_mean` against the scale the data arrived on
#'
#' The default depends on another argument, which a formal default cannot
#' express, so the caller reports whether the user supplied one. On the log2
#' scale the geometric mean is the convention and the only choice that reduces
#' to a difference of means, so it wins by default there.
#'
#' @param fc_mean The argument as received.
#' @param input_scale `"raw"` or `"log2"`, already matched.
#' @param use_default Result of `missing(fc_mean)` in the calling function.
#'
#' @keywords internal
#' @noRd
sa_resolve_fc_mean <- function(fc_mean, input_scale, use_default) {
  if (use_default) {
    return(if (input_scale == "log2") "geom" else "arith")
  }
  match.arg(fc_mean, c("arith", "geom"))
}


#' Central tendency of one group for the fold change ratio
#'
#' @param v Numeric vector with missing values already removed.
#' @param side Group label, used in the error message.
#' @param mean_type `"arith"` or `"geom"`.
#' @param input_scale `"raw"` or `"log2"`.
#'
#' @keywords internal
#' @noRd
sa_fc_center <- function(v, side, mean_type, input_scale = "raw") {
  if (length(v) == 0L) {
    stop("no usable observation left in the ", side, " group.", call. = FALSE)
  }
  # A ratio is only a ratio on the original scale. Undoing the transformation
  # here rather than dividing the log2 values keeps every centre downstream on
  # one scale, so `fold_change == x_center / y_center` still holds.
  if (input_scale == "log2") {
    v <- 2^v
    n_over <- sum(!is.finite(v))
    if (n_over > 0L) {
      stop("2^x overflows to infinity for ", n_over, " value(s) in the ", side,
           " group, so these observations are not on the log2 scale; use ",
           "`input_scale = \"raw\"` instead.", call. = FALSE)
    }
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
#' @param input_scale `"raw"` or `"log2"`.
#'
#' @return data.frame with `features`, `x_center`, `y_center`, `fold_change` and
#'   `log2fc`.
#'
#' @keywords internal
#' @noRd
sa_fold_change <- function(samples, feats, group_lv, mean_type,
                           input_scale = "raw") {
  label <- paste0(if (mean_type == "arith") "Arithmetic" else "Geometric",
                  " mean fold change")

  out <- sa_feature_table(
    feats,
    c("x_center", "y_center", "fold_change", "log2fc"),
    label,
    p_adjust = NULL,
    fun = function(i) {
      x_center <- sa_fc_center(samples[[i]]$x, group_lv[1], mean_type,
                               input_scale)
      y_center <- sa_fc_center(samples[[i]]$y, group_lv[2], mean_type,
                               input_scale)
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
