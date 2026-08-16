#' Draw a grouped barplot of a descriptive summary
#'
#' Draws one bar per feature and group level, the levels side by side inside
#' each feature's cluster and a legend in a narrow panel on the right. The bar
#' heights are one column of [summarize_descriptive_stats()], read from the same
#' wide input the comparison functions take, so a bar and a row of that table
#' are the same number and neither has to be recomputed to check the other.
#'
#' This is the summary counterpart of [draw_grouped_boxplot()]. A box shows the
#' distribution a group's observations have; a bar shows one number standing for
#' them, which is less of the data and is what a figure wants when the point
#' being made is about a location rather than a spread.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to plot, in
#'   display order along the x axis.
#' @param group Grouping vector with one entry per row of `data`. Required: a
#'   summary of every row together is what [summarize_descriptive_stats()]
#'   returns without one, and it has no clusters to draw.
#' @param group_lv Character vector of at least two group levels, in the order
#'   they should appear inside each cluster. Rows belonging to any other level
#'   are dropped. `NULL` takes the levels from `group` the way
#'   [summarize_descriptive_stats()] does: the factor levels if `group` is a
#'   factor, the sorted unique values otherwise.
#' @param control_label The level to hold as the reference. It moves to the
#'   front of `group_lv` and the rest keep the order they were given, the same
#'   move [compare_two_groups()] and [compare_multiple_groups()] make when the
#'   same argument is passed, so a figure and the analysis it illustrates put
#'   the reference bar in the same place when they are given the same argument.
#'   Read only when `group_lv` is supplied or derived: it re-points a list that
#'   already exists rather than stating the reference for the first time.
#' @param mainbar Which column of [summarize_descriptive_stats()] the bar
#'   heights are. See "What a bar may carry an interval for" for which of them
#'   also take an `errorbar`.
#' @param errorbar What the bars either side of a height are. `"none"`,
#'   `"se"` for one standard error, `"sd"` for one standard deviation, or
#'   `"ci"` for an interval at `conf_level`. Read under `mainbar`, which
#'   decides which of them the height has an answer for.
#' @param conf_level Confidence level of `errorbar = "ci"`, read only under
#'   `mainbar = "mean"`. The median's interval is the notch, whose width is
#'   fixed at the `1.58` that makes it an approximate 95% interval, so there is
#'   no level for it to be stated at.
#' @param gap Blank bar widths inserted between neighbouring clusters.
#' @param lwd Line width of the error bars and of the baseline when the heights
#'   run both ways.
#' @param col Fill colours for the group levels, recycled if short. `NULL`
#'   takes them from [grDevices::hcl.colors()], the palette
#'   [draw_grouped_boxplot()] and [draw_interaction_plot()] already draw group
#'   levels in.
#' @param xlab,ylab,main Axis and title labels. `ylab` left `NULL` names
#'   `mainbar`, since the axis of a bar chart is a particular statistic rather
#'   than the measurement itself; pass `""` for no label at all.
#' @param ylim Numeric length-2 y axis range, or `NULL` to derive one that
#'   covers the bars, their intervals and the zero they are measured from.
#' @param dark If `TRUE`, use a dark background with light text.
#' @param grid_lty,grid_lwd Line type and width of the horizontal grid.
#' @param cex.lab,cex.axis,cex.main,cex.legend Character expansion for axis
#'   labels, axis annotation, the main title and the legend.
#' @param out_statistics If `TRUE`, invisibly return the numbers behind the
#'   bars.
#'
#' @return If `out_statistics = FALSE`, `NULL` invisibly. Otherwise, invisibly,
#'   a data.frame of one row per bar in the order they were drawn, which is the
#'   row order of [summarize_descriptive_stats()]: a feature's levels stay
#'   together, in `group_lv` order.
#'
#'   \describe{
#'     \item{`features`,`group`}{Which bar the row is.}
#'     \item{`n`}{Finite observations the bar was computed from.}
#'     \item{`value`}{The bar height, the `mainbar` column.}
#'     \item{`lower`,`upper`}{The ends of the interval, `NA` under
#'       `errorbar = "none"` and for a bar whose interval was not defined.}
#'   }
#'
#'   Which height and which interval were drawn are attached as the attributes
#'   `"mainbar"` and `"errorbar"`.
#'
#' @section What a bar may carry an interval for:
#' `mainbar` names any of the summary columns, but only two of them are
#' locations that an interval either side says something about, so `errorbar` is
#' read under `mainbar` rather than independently of it.
#'
#' \describe{
#'   \item{`"mean"`}{Takes every bar. `"se"` and `"sd"` are one standard error
#'     and one standard deviation either side, and `"ci"` is Student's interval
#'     at `conf_level`: [stats::qt()] on `n - 1` degrees of freedom times the
#'     standard error.}
#'   \item{`"median"`}{Takes `"ci"` only, the notch interval
#'     `median +/- 1.58 * IQR / sqrt(n)`. That is the same interval
#'     [draw_grouped_boxplot()] returns as `median_confidence_stats`, so a bar
#'     and the notch of the box beside it are the same width on the same data.
#'     `"se"` and `"sd"` describe the observations' spread about their mean and
#'     are refused here rather than drawn around a median they are not about.}
#'   \item{Everything else}{Takes `"none"`. A height that is itself a spread, a
#'     count or a shape has no second quantity for an interval to be about.}
#' }
#'
#' A refused combination is an error rather than a silently dropped interval,
#' because the alternative is a figure that answers a question other than the
#' one it was asked.
#'
#' @details
#' A bar is read against the zero it stands on, so a derived `ylim` always
#' includes zero and puts its headroom on the side the bars run to. Heights that
#' go both ways, which `"skewness"` and a `"mean"` of centred features do, get
#' the baseline drawn as well.
#'
#' A bar whose height is `NA` leaves a blank rather than shifting the ones
#' beside it, which is what makes a group too small for a shape estimate visible
#' instead of silently absent.
#'
#' The function changes graphical parameters and the panel layout, and restores
#' both on exit, so the caller's device is left as it was found.
#'
#' @seealso [summarize_descriptive_stats()] for the table the heights come from,
#'   [draw_grouped_boxplot()] for the observations behind them, and
#'   [compare_two_groups()] or [compare_multiple_groups()] to test the
#'   difference a pair of bars shows.
#'
#' @references
#' McGill, R., Tukey, J. W. and Larsen, W. A. (1978). Variations of box plots.
#' *The American Statistician*, 32(1), 12-16.
#'
#' @examples
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#'
#' ## Mean per species, one standard error either side
#' bars <- draw_grouped_barplot(iris, feats, iris$Species, errorbar = "se")
#' head(bars)
#'
#' ## The heights are the table's own column, so neither has to be checked
#' ## against the other by eye.
#' summ <- summarize_descriptive_stats(iris, feats, iris$Species)
#' all.equal(bars$value, summ$mean)
#'
#' ## Student's interval, at a level of your own
#' draw_grouped_barplot(iris, feats, iris$Species, errorbar = "ci",
#'                      conf_level = 0.99)
#'
#' ## A median takes the notch interval draw_grouped_boxplot() notches with
#' draw_grouped_barplot(iris, feats, iris$Species, mainbar = "median",
#'                      errorbar = "ci")
#'
#' ## A count carries no interval, and `group_lv` picks the levels and their
#' ## order
#' draw_grouped_barplot(iris, feats, iris$Species,
#'                      group_lv = c("virginica", "versicolor"),
#'                      mainbar = "n")
#'
#' ## `control_label` draws the named level first without rewriting `group_lv`
#' draw_grouped_barplot(iris, feats, iris$Species,
#'                      control_label = "setosa")
#'
#' @export
draw_grouped_barplot <- function(data,
                                 feats,
                                 group = NULL,
                                 group_lv = NULL,
                                 control_label = NULL,
                                 mainbar = c("mean", "median", "n", "n_missing",
                                             "sd", "var", "se", "cv", "mad",
                                             "skewness", "excess_kurtosis"),
                                 errorbar = c("none", "se", "sd", "ci"),
                                 conf_level = 0.95,
                                 gap = 1,
                                 lwd = 1.5,
                                 col = NULL,
                                 xlab = NULL,
                                 ylab = NULL,
                                 main = NULL,
                                 ylim = NULL,
                                 dark = FALSE,
                                 grid_lty = 1,
                                 grid_lwd = 0.25,
                                 cex.lab = 1.3,
                                 cex.axis = 1.2,
                                 cex.main = 1.3,
                                 cex.legend = 1.1,
                                 out_statistics = TRUE) {

  mainbar <- match.arg(mainbar)
  errorbar <- match.arg(errorbar)
  sa_bar_check_pair(mainbar, errorbar)

  sa_check_scalar_num(conf_level, "conf_level", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)
  sa_check_scalar_num(gap, "gap", 0)
  sa_check_scalar_num(lwd, "lwd", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.lab, "cex.lab", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.axis, "cex.axis", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.main, "cex.main", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.legend, "cex.legend", 0, lower_open = TRUE)
  sa_check_scalar_num(grid_lwd, "grid_lwd", 0)
  sa_check_flag(dark, "dark")
  sa_check_flag(out_statistics, "out_statistics")
  sa_check_lim(ylim, "ylim")

  input <- sa_bar_input(data, feats, group, group_lv, control_label)
  drawn <- sa_bar_values(input, mainbar, errorbar, conf_level)

  sa_bar_draw(drawn, input, mainbar = mainbar, errorbar = errorbar, gap = gap,
              lwd = lwd, col = col, xlab = xlab, ylab = ylab, main = main,
              ylim = ylim, dark = dark, grid_lty = grid_lty, grid_lwd = grid_lwd,
              cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
              cex.legend = cex.legend)

  if (!out_statistics) {
    return(invisible(NULL))
  }

  attr(drawn, "mainbar") <- mainbar
  attr(drawn, "errorbar") <- errorbar
  invisible(drawn)
}


#' Refuse an interval the height has no answer for
#'
#' Checked before anything is read off `data`, so a call that could only produce
#' a misleading figure fails at the boundary rather than after the summary has
#' been computed.
#'
#' @keywords internal
#' @noRd
sa_bar_check_pair <- function(mainbar, errorbar) {
  if (errorbar == "none" || mainbar == "mean") {
    return(invisible(NULL))
  }

  if (mainbar == "median") {
    if (errorbar == "ci") {
      return(invisible(NULL))
    }
    stop("`errorbar = \"", errorbar, "\"` describes the spread of the ",
         "observations about their mean, so it is not a width to draw either ",
         "side of a median. `mainbar = \"median\"` takes ",
         "errorbar = \"ci\", the notch interval ",
         "median +/- 1.58 * IQR / sqrt(n), or \"none\".", call. = FALSE)
  }

  stop("`mainbar = \"", mainbar, "\"` is itself a spread, a count or a shape, ",
       "so there is no second quantity for an interval either side of it to be ",
       "about. Only \"mean\" and \"median\" take an `errorbar`; this height ",
       "takes errorbar = \"none\".", call. = FALSE)
}


#' Resolve the bars and summarise them
#'
#' The summary is asked for the levels that survived the validation rather than
#' for the ones the caller named, so it has nothing left to drop and does not
#' report the same rows a second time.
#'
#' @return List with `feats`, `lv` in draw order, and `summ`, the descriptive
#'   summary whose rows are the bars.
#'
#' @keywords internal
#' @noRd
sa_bar_input <- function(data, feats, group, group_lv, control_label = NULL) {
  if (is.null(group)) {
    stop("`group` says which bars there are, so it is required: one entry per ",
         "row of `data`. A summary of every row together, with no clusters to ",
         "draw, is what summarize_descriptive_stats() returns without one.",
         call. = FALSE)
  }
  if (is.null(group_lv)) {
    group_lv <- if (is.factor(group)) {
      levels(droplevels(group))
    } else {
      sort(unique(as.character(group)))
    }
  }
  group_lv <- sa_control_first(group_lv, control_label)

  input <- sa_validate_wide_input(data, feats, group, group_lv,
                                 min_levels = 2L)
  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }

  list(
    feats = input$feats,
    lv    = levels(input$group),
    summ  = summarize_descriptive_stats(input$data, input$feats, input$group,
                                       levels(input$group))
  )
}


#' One row per bar: its height and the ends of its interval
#'
#' @keywords internal
#' @noRd
sa_bar_values <- function(input, mainbar, errorbar, conf_level) {
  summ <- input$summ
  interval <- sa_bar_interval(summ, mainbar, errorbar, conf_level)

  data.frame(
    features = summ$features,
    group    = as.character(summ$group),
    n        = summ$n,
    value    = summ[[mainbar]],
    lower    = interval$lower,
    upper    = interval$upper,
    stringsAsFactors = FALSE
  )
}


#' Half-widths either side of the height, per bar
#'
#' Every quantity involved is already a column of the summary, so the interval
#' is arithmetic on the table rather than a second pass over the observations.
#'
#' @return List of `lower` and `upper`, one entry per row of `summ`.
#'
#' @keywords internal
#' @noRd
sa_bar_interval <- function(summ, mainbar, errorbar, conf_level) {
  if (errorbar == "none") {
    return(list(lower = rep(NA_real_, nrow(summ)),
                upper = rep(NA_real_, nrow(summ))))
  }

  half <- switch(
    errorbar,
    se = summ$se,
    sd = summ$sd,
    ci = if (mainbar == "mean") {
      # `pmax()` keeps qt() off zero degrees of freedom, which it answers with a
      # NaN and a warning; a single observation has no interval either way.
      ifelse(
        summ$n > 1L,
        stats::qt(1 - (1 - conf_level) / 2, pmax(summ$n - 1L, 1L)) * summ$se,
        NA_real_
      )
    } else {
      1.58 * summ$iqr / sqrt(summ$n)
    }
  )

  centre <- summ[[mainbar]]
  list(lower = centre - half, upper = centre + half)
}


#' The y range the bars, their intervals and their baseline need
#'
#' @param ylim The caller's range, or `NULL` to derive one.
#'
#' @keywords internal
#' @noRd
sa_bar_span <- function(drawn, mainbar, errorbar, ylim) {
  if (!any(is.finite(drawn$value))) {
    stop("`mainbar = \"", mainbar, "\"` is NA for every feature and group, so ",
         "there is no bar to draw. summarize_descriptive_stats() returns the ",
         "same column: a shape estimate needs three or four observations, and ",
         "every statistic needs one.", call. = FALSE)
  }
  if (!is.null(ylim)) {
    return(ylim)
  }

  # Zero goes in whether or not a bar reaches it, since it is what the height of
  # a bar is measured from.
  vals <- c(0, drawn$value)
  if (errorbar != "none") {
    vals <- c(vals, drawn$lower, drawn$upper)
  }
  span <- range(vals[is.finite(vals)])

  if (diff(span) == 0) {
    # Every bar sits on the baseline, which needs a panel to be a flat line in.
    pad <- max(abs(span), 1) * 0.1
    return(span + c(-pad, pad))
  }

  # Headroom on the side the bars run to and none on the baseline, which a bar
  # has to stand on rather than float above.
  span + c(if (span[1] < 0) -0.04 * diff(span) else 0,
           if (span[2] > 0) 0.04 * diff(span) else 0)
}


#' Draw the clusters of bars and the legend beside them
#'
#' @param drawn One row per bar, as `sa_bar_values()` returns.
#' @param input As `sa_bar_input()` returns.
#'
#' @keywords internal
#' @noRd
sa_bar_draw <- function(drawn, input, mainbar, errorbar, gap, lwd, col, xlab,
                        ylab, main, ylim, dark, grid_lty, grid_lwd, cex.lab,
                        cex.axis, cex.main, cex.legend) {
  feats <- input$feats
  lv <- input$lv
  n_lv <- length(lv)

  # Filled column by column, so a column is a feature's cluster and the matrix
  # reads in the row order of `drawn`: one bar of it is one row of that.
  m <- matrix(drawn$value, nrow = n_lv, ncol = length(feats),
              dimnames = list(lv, feats))

  span <- sa_bar_span(drawn, mainbar, errorbar, ylim)
  cols <- if (is.null(col)) {
    grDevices::hcl.colors(n_lv, "Dark 2")
  } else {
    rep_len(col, n_lv)
  }
  theme <- sa_plot_theme(dark)
  grid_col <- if (dark) "gray80" else "gray40"
  if (is.null(ylab)) {
    ylab <- mainbar
  }

  # Only the parameters this function sets are put back, not a blanket
  # par(no.readonly = TRUE) snapshot, which also carries the absolute sizes
  # `fin`, `pin` and `mai`: restoring those pins the next plot to the size this
  # one happened to be drawn at. `mfrow` comes back because layout() overwrites
  # whatever grid the caller had set up.
  old_par <- graphics::par(c("bg", "fg", "col.axis", "col.lab", "col.main",
                             "mar", "mfrow", "oma"))
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  if (dark) {
    graphics::par(bg = theme$bg, fg = theme$fg, col.axis = theme$fg,
                  col.lab = theme$fg, col.main = theme$fg)
  }
  graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(4, 1))
  graphics::par(mar = c(5.1, 4.1, if (is.null(main)) 2.1 else 4.1, 1.1))

  # Drawn twice with the grid in between so the grid sits behind the bars. One
  # closure keeps the two calls from drifting apart. The annotation is left to
  # the calls below rather than taken from either of these, which would write a
  # title once per pass.
  draw_bars <- function(add) {
    graphics::barplot(m, beside = TRUE, space = c(0, gap), col = cols,
                      border = NA, axes = FALSE, axisnames = FALSE,
                      ylim = span, xpd = FALSE, add = add)
  }

  at <- draw_bars(add = FALSE)
  graphics::grid(
    nx = NA,               # no vertical grid
    ny = NULL,             # horizontal grid at the y axis ticks
    col = grid_col,
    lty = grid_lty,
    lwd = grid_lwd
  )
  draw_bars(add = TRUE)

  if (errorbar != "none") {
    # `at` was filled the way `m` was, so its cells and the rows of `drawn` are
    # the same bars in the same order.
    has_bar <- is.finite(drawn$lower) & is.finite(drawn$upper) &
      drawn$upper > drawn$lower
    if (any(has_bar)) {
      graphics::arrows(as.vector(at)[has_bar], drawn$lower[has_bar],
                       as.vector(at)[has_bar], drawn$upper[has_bar],
                       length = 0.04, angle = 90, code = 3, col = theme$fg,
                       lwd = lwd)
    }
  }

  # A baseline is only worth drawing where it is not already the floor of the
  # panel: bars that run both ways need the zero they are measured from.
  if (span[1] < 0) {
    graphics::abline(h = 0, col = theme$fg, lwd = lwd)
  }

  graphics::axis(2, las = 1, cex.axis = cex.axis)
  graphics::axis(1, at = colMeans(at), labels = feats, tick = FALSE,
                 cex.axis = cex.axis)
  graphics::title(xlab = xlab, ylab = ylab, main = main, cex.lab = cex.lab,
                  cex.main = cex.main)

  graphics::par(mar = c(5, 0, 4, 1))
  graphics::plot.new()
  graphics::legend("center", legend = lv, fill = cols, border = NA, bty = "n",
                   cex = cex.legend, text.col = theme$fg)

  invisible(NULL)
}
