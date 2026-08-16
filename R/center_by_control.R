#' Centre every feature on the control group
#'
#' Divides, or on the log2 scale subtracts, the centre of the control group out
#' of each feature, so that the control lands on 1 (raw) or 0 (log2) and every
#' other observation reads as its distance from the control rather than as a
#' measurement of its own. The `data` that comes back is the `data` that went in
#' with the `feats` columns replaced, which is what lets the same call be handed
#' straight to a comparison or to [draw_heatmap()].
#'
#' The arguments are the ones [compare_two_groups()] and
#' [compare_multiple_groups()] already take, in the same order and with the same
#' defaults, because the point of the function is that one set of arguments
#' describes both steps:
#'
#' ```
#' centred <- center_by_control(data, feats, group, group_lv,
#'                              input_scale = "log2")
#' compare_two_groups(centred, feats, group, group_lv, input_scale = "log2")
#' ```
#'
#' `fc_mean` is the comparison's own argument rather than a centring method of
#' this function's, and it is resolved the same way: the geometric mean by
#' default on the log2 scale, where it is the convention, and the arithmetic
#' mean otherwise. Passing the same value to both is what keeps the baseline
#' removed here identical to the centre the comparison divides in its `effect`
#' table.
#'
#' Two or more levels are accepted, so the same call serves a two-group and a
#' multi-group design. Only the control level takes part in the baseline; the
#' rest are centred on it.
#'
#' @section What the comparison reports afterwards:
#'
#' The transformation is one constant per feature applied to every row of it, so
#' most of what a comparison reports does not move:
#'
#' \describe{
#'   \item{`fold_change` and `log2fc`}{Unchanged. Both group centres are divided
#'     by the same baseline, so their ratio survives it. This holds for all four
#'     combinations of `fc_mean` and `input_scale`.}
#'   \item{p-values}{Unchanged, in every test family. On the log2 scale the
#'     baseline is subtracted, and the tests are shift invariant; on the raw
#'     scale it is divided, and dividing by a positive constant leaves the t
#'     statistic and every rank untouched.}
#'   \item{The reference centre}{Becomes 1. It is `y_center` in
#'     [compare_two_groups()] and `ref_center` in [compare_multiple_groups()],
#'     and once it is 1 the other centre in the row is the fold change itself.}
#'   \item{`mean_diff`, `hl_shift` and `trim_diff`}{Unchanged on the log2 scale.
#'     On the raw scale they are divided by the baseline along with the data,
#'     since they are differences rather than ratios.}
#' }
#'
#' The reference centre lands on 1 because the observations used here and the
#' ones the comparison centres are the same. Under `paired = TRUE` a comparison
#' keeps complete pairs or complete subjects only, so if any are dropped its
#' control centre is taken on fewer rows than the baseline was and comes back
#' near 1 rather than at it. The ratios are unaffected either way.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to centre.
#'   Columns not named here are returned untouched.
#' @param group Grouping vector with one entry per row of `data`.
#' @param group_lv Character vector of at least two group levels. Unlike in a
#'   comparison, rows belonging to another level are kept: they are centred on
#'   the same baseline and take no part in computing it, which leaves the result
#'   the same length as the `group` the comparison still has to be given. The
#'   comparison drops them itself.
#' @param control_label The level whose centre is divided out. Defaults to
#'   `group_lv[1]`, the reference the order already asked for, so naming it is
#'   only needed to point the baseline at a level that is not first.
#' @param fc_mean Which centre of the control group is removed, `"arith"` for
#'   the arithmetic mean or `"geom"` for the geometric mean. The geometric mean
#'   requires strictly positive values. Defaults to `"geom"` when
#'   `input_scale = "log2"` and to `"arith"` otherwise, exactly as it does in
#'   the comparisons.
#' @param input_scale The scale `data` arrives on, `"raw"` or `"log2"`. This
#'   decides the operation: a ratio on the raw scale, a difference on the log2
#'   scale. The centred data stays on the scale it arrived on, so the same
#'   `input_scale` is what the comparison should be given afterwards.
#'
#' @return The `data` that was passed in, as a data.frame, with the `feats`
#'   columns replaced by their control-relative values. Row count, row order,
#'   row names and every column that is not a feature are left as they arrived.
#'   A feature whose baseline could not be taken comes back as an `NA` column
#'   and is named in a warning.
#'
#' @seealso [compare_two_groups()] and [compare_multiple_groups()], whose
#'   arguments this shares, and [draw_heatmap()], which takes the result with
#'   `scale = "none"` to show departure from the control rather than from each
#'   feature's own mean.
#'
#' @examples
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#' group_lv <- levels(iris$Species)
#'
#' ## setosa is the first level and so the control. Every measurement is now a
#' ## multiple of the setosa mean, and setosa itself centres on 1.
#' centred <- center_by_control(iris, feats, iris$Species, group_lv)
#' colMeans(centred[iris$Species == "setosa", feats])
#'
#' ## The comparison takes the centred data with its arguments unchanged, and
#' ## reports the fold changes it reported before the centring.
#' before <- compare_multiple_groups(iris, feats, iris$Species, group_lv,
#'                                   posthoc = FALSE, diagnose = FALSE)
#' after <- compare_multiple_groups(centred, feats, iris$Species, group_lv,
#'                                  posthoc = FALSE, diagnose = FALSE)
#' cbind(before = before$effect$log2fc, after = after$effect$log2fc)
#'
#' ## What did change is the reference centre, which is now 1 by construction,
#' ## so `extreme_center` reads as the fold change on its own.
#' after$effect[c("ref_center", "extreme_center", "fold_change")]
#'
#' ## On the log2 scale the baseline is subtracted instead, and the control
#' ## centres on 0. Rounded, because the subtraction leaves floating point dust
#' ## rather than a clean zero.
#' d <- iris[iris$Species != "setosa", ]
#' lg <- center_by_control(log2(d[feats]), feats, d$Species,
#'                         c("versicolor", "virginica"),
#'                         input_scale = "log2")
#' round(colMeans(lg[d$Species == "versicolor", ]), 12)
#'
#' ## Pointing the baseline at the other level is one argument rather than a
#' ## rewritten `group_lv`.
#' virg <- center_by_control(iris, feats, iris$Species, group_lv,
#'                           control_label = "virginica")
#' colMeans(virg[iris$Species == "virginica", feats])
#'
#' @export
center_by_control <- function(data,
                              feats,
                              group,
                              group_lv,
                              control_label = group_lv[1],
                              fc_mean = c("arith", "geom"),
                              input_scale = c("raw", "log2")) {

  input_scale <- match.arg(input_scale)
  fc_mean <- sa_resolve_fc_mean(fc_mean, input_scale, missing(fc_mean))

  # Every check a comparison makes on these arguments, made here in the same
  # order and with the same messages. Only the errors and `n_dropped` are taken
  # from the result: its `data` has the rows outside `group_lv` removed, and
  # returning that would leave the output shorter than the `group` vector the
  # caller still has to hand to the comparison.
  input <- sa_validate_wide_input(data, feats, group, group_lv)
  control <- sa_control_first(levels(input$group), control_label)[1]

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  if (input$n_dropped > 0L) {
    message("Kept ", input$n_dropped, " row(s) belonging to a level outside ",
            "`group_lv`. They are centred on the same baseline but take no ",
            "part in it, and the comparison drops them itself.")
  }

  # `which()` rather than a logical mask. A row outside `group_lv` compares NA
  # against the control name, and an NA index would put an NA into the baseline
  # sample and take the whole feature down with it.
  ctrl_rows <- which(as.character(group) == control)

  failures <- character(0)
  for (f in feats) {
    v <- data[[f]]
    baseline <- tryCatch(
      sa_control_baseline(v[ctrl_rows], control, fc_mean, input_scale),
      error = function(e) {
        failures[[f]] <<- conditionMessage(e)
        NA_real_
      }
    )
    data[[f]] <- if (input_scale == "log2") v - baseline else v / baseline
  }

  # One warning naming every feature, rather than one per feature: a scan over
  # hundreds of columns must not be abandoned because one of them has no usable
  # control group.
  if (length(failures) > 0L) {
    warning("The control baseline could not be taken for ", length(failures),
            " of ", length(feats), " feature(s); those columns are all NA:\n",
            paste0("  ", names(failures), ": ", failures, collapse = "\n"),
            call. = FALSE)
  }

  data
}


#' The control centre, on the scale the data arrived on
#'
#' `sa_fc_center()` always reports on the raw scale, since that is the only
#' scale on which a ratio is a ratio, so a log2 input needs its centre brought
#' back before it can be subtracted from log2 values.
#'
#' The non-positive raw centre is rejected rather than reported.
#' `sa_fold_change()` only messages about one, because a ratio against a zero or
#' negative centre is a strange number in a table the reader can see. Here it
#' would divide the feature itself: a zero centre sends every value to infinity,
#' and a negative one reverses the order of all of them, which silently flips
#' the direction of every rank-based test run on the result afterwards.
#'
#' @param v Control group values for one feature, missing values included.
#' @param control The control level name, used in the error messages.
#' @param mean_type `"arith"` or `"geom"`.
#' @param input_scale `"raw"` or `"log2"`.
#'
#' @return The quantity to divide out (raw) or to subtract (log2).
#'
#' @keywords internal
#' @noRd
sa_control_baseline <- function(v, control, mean_type, input_scale) {
  # sa_fc_center() takes the missing values as already gone, the way the
  # comparisons hand it their samples.
  centre <- sa_fc_center(v[!is.na(v)], control, mean_type, input_scale)

  if (input_scale == "log2") {
    baseline <- log2(centre)
    if (!is.finite(baseline)) {
      stop("the ", control, " centre is ", centre, " on the raw scale, which ",
           "has no log2 to subtract.", call. = FALSE)
    }
    return(baseline)
  }

  if (!is.finite(centre) || centre <= 0) {
    consequence <- if (!is.finite(centre)) {
      "send every value to zero"
    } else if (centre == 0) {
      "send every value to infinity"
    } else {
      paste("reverse the order of every value, which would flip the direction",
            "of the rank-based tests run on the result")
    }
    stop("the ", control, " centre is ", centre, ", and dividing by it would ",
         consequence, ". Pass logged values with ",
         "`input_scale = \"log2\"` instead, or restrict `feats` to features ",
         "whose control group is positive.", call. = FALSE)
  }
  centre
}
