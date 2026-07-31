#' Draw a butterfly histogram of one feature across two groups
#'
#' Draws the distribution of a single feature for two group levels back to
#' back, on a shared set of breaks. The first level runs left from the centre
#' line and the second one right, so the two shapes can be compared bin by bin
#' instead of read off two separate panels.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feat Name of the numeric column in `data` to plot. One feature per
#'   call.
#' @param group Grouping vector with one entry per row of `data`.
#' @param group_lv Character vector of exactly two group levels. The first is
#'   drawn on the left, the second on the right. Rows belonging to any other
#'   level are dropped.
#' @param breaks Histogram break specification, shared by both groups. Accepts
#'   a character rule supported by [graphics::hist()], a single number read as
#'   the approximate bin count, or a strictly increasing vector of break
#'   points.
#' @param scale What the bar length means: `"count"`, `"proportion"` of the
#'   group, or `"density"`. Proportions and densities are computed within each
#'   group, so unequal group sizes still give comparable shapes.
#' @param col Length-2 vector of fill colours, for the left and right group.
#' @param border Colour of the bar borders.
#' @param xlab,ylab,main Axis labels and title. `xlab` defaults to the meaning
#'   of `scale` and `ylab` to `feat`.
#' @param cex.lab,cex.axis,cex.main,cex.legend Character expansion for the axis
#'   labels, the axis annotation, the title and the legend.
#' @param legend.position Position passed to [graphics::legend()], or `NULL` /
#'   `FALSE` to leave the legend out.
#' @param xlim,ylim Numeric length-2 ranges for the bar axis and the value
#'   axis, or `NULL` to derive them from the data. `xlim` is not forced to be
#'   symmetric when it is supplied.
#' @param margin Plot margins in lines, passed to [graphics::par()] as `mar`.
#' @param out_statistics If `TRUE`, invisibly return the numbers behind the
#'   bars.
#' @param ... Additional arguments passed to [graphics::plot.default()].
#'
#' @return If `out_statistics = FALSE`, `NULL` invisibly. Otherwise a list of
#'   two elements, invisibly:
#'
#'   \describe{
#'     \item{`bin_summary_stats`}{One row per bin, with `bin_start`, `bin_end`,
#'       `bin_mid` and one column per group level holding the plotted bar
#'       length. What that length means follows `scale`.}
#'     \item{`group_summary_stats`}{One column per group level with rows `n`
#'       (finite values used), `n_dropped` (missing or non-finite values left
#'       out), `min` and `max`.}
#'     \item{`group_hists`}{One `"histogram"` object per group level, named by
#'       the level, exactly as [graphics::hist()] returns them, so
#'       `plot(res$group_hists[[lv]])` redraws that group on its own. Both
#'       groups are binned on the shared `breaks`, so these are not what
#'       calling [graphics::hist()] on one group alone would give: that would
#'       pick its own breaks from that group only.}
#'   }
#'
#' @details
#' The left group is drawn at negative coordinates so that both distributions
#' share one axis, but the tick labels are the absolute values, so a bar is
#' read the same way on either side.
#'
#' The plot starts from the device background rather than imposing a theme of
#' its own, as [draw_volcano_plot()] does. Only two colours are needed here, so
#' there is no dark variant; that is reserved for [draw_grouped_boxplot()],
#' where three or more group colours have to stay apart.
#'
#' Only the margins are restored on exit, not every graphical parameter, so the
#' coordinate system survives the call and [graphics::abline()],
#' [graphics::points()] and friends can still add to the finished plot.
#'
#' @seealso [draw_grouped_boxplot()] for several features at once, and
#'   [compare_two_groups()] to test the difference the plot shows.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#'
#' draw_butterfly_hist(
#'   data     = iris2,
#'   feat     = "Petal.Length",
#'   group    = iris2$Species,
#'   group_lv = c("versicolor", "virginica"),
#'   breaks   = 12,
#'   scale    = "proportion",
#'   main     = "virginica vs versicolor"
#' )
#'
#' ## Counts, plus the numbers behind the bars
#' res <- draw_butterfly_hist(iris2, "Sepal.Width", iris2$Species,
#'                            c("versicolor", "virginica"))
#' res$bin_summary_stats
#' res$group_summary_stats
#'
#' ## Each group also comes back as a plain histogram object
#' plot(res$group_hists$virginica)
#'
#' @export
draw_butterfly_hist <- function(data,
                                feat,
                                group,
                                group_lv,
                                breaks = "Sturges",
                                scale = c("count", "proportion", "density"),
                                col = c("#4575B4", "#D73027"),
                                border = "white",
                                xlab = NULL,
                                ylab = NULL,
                                main = NULL,
                                cex.lab = 1.3,
                                cex.axis = 1.2,
                                cex.main = 1.3,
                                cex.legend = 1.1,
                                legend.position = "topright",
                                xlim = NULL,
                                ylim = NULL,
                                margin = c(5, 5, 4, 3),
                                out_statistics = TRUE,
                                ...) {

  scale <- match.arg(scale)

  sa_check_scalar_num(cex.lab, "cex.lab", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.axis, "cex.axis", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.main, "cex.main", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.legend, "cex.legend", 0, lower_open = TRUE)
  sa_check_flag(out_statistics, "out_statistics")
  sa_check_margin(margin)
  sa_check_lim(xlim, "xlim")
  sa_check_lim(ylim, "ylim")

  if (length(col) != 2L || anyNA(col)) {
    stop("`col` must contain exactly two colours, one per group level.",
         call. = FALSE)
  }
  if (length(border) != 1L) {
    stop("`border` must be a single colour.", call. = FALSE)
  }
  if (!is.character(feat) || length(feat) != 1L || is.na(feat)) {
    stop("`feat` must be a single column name; this plot shows one feature ",
         "at a time.", call. = FALSE)
  }

  input <- sa_validate_wide_input(data, feat, group, group_lv, n_levels = 2L)
  data <- input$data
  group <- input$group
  group_lv <- levels(group)

  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }

  values <- lapply(group_lv, function(lv) data[[feat]][group == lv])
  names(values) <- group_lv
  n_input <- lengths(values)
  values <- lapply(values, function(v) v[is.finite(v)])
  n_used <- lengths(values)

  empty_lv <- group_lv[n_used == 0L]
  if (length(empty_lv) > 0L) {
    stop("`", feat, "` has no finite value in group level(s): ",
         paste(empty_lv, collapse = ", "), ".", call. = FALSE)
  }
  if (any(n_input > n_used)) {
    message("Left out ", sum(n_input - n_used),
            " missing or non-finite value(s) of `", feat, "`.")
  }

  combined_range <- range(unlist(values, use.names = FALSE))

  if (length(breaks) == 1L && is.character(breaks)) {
    breaks <- graphics::hist(
      unlist(values, use.names = FALSE),
      breaks = breaks,
      plot = FALSE
    )$breaks
  } else if (length(breaks) == 1L && is.numeric(breaks) &&
             is.finite(breaks) && breaks > 0) {
    breaks <- pretty(combined_range, n = breaks)
  } else if (!is.numeric(breaks) || length(breaks) < 2L ||
             any(!is.finite(breaks)) ||
             is.unsorted(breaks, strictly = TRUE)) {
    stop("`breaks` must be a valid histogram rule, a positive bin count, ",
         "or a strictly increasing numeric vector.", call. = FALSE)
  }

  # Both groups go through the same breaks, so bin i means the same interval on
  # either side of the centre line.
  hists <- lapply(group_lv, function(lv) {
    h <- graphics::hist(values[[lv]], breaks = breaks, plot = FALSE,
                        include.lowest = TRUE)
    # hist() takes xname from the deparsed argument, which is the closure's
    # `values[[lv]]` here, so plot() on the returned object would be titled
    # after that instead of the feature.
    h$xname <- paste0(feat, " (", lv, ")")
    h
  })
  names(hists) <- group_lv

  bar_length <- lapply(hists, function(h) {
    switch(
      scale,
      count = h$counts,
      proportion = h$counts / sum(h$counts),
      density = h$density
    )
  })

  bar_max <- max(unlist(bar_length, use.names = FALSE), na.rm = TRUE)
  if (!is.finite(bar_max) || bar_max <= 0) {
    stop("every bin is empty, so there is nothing to draw.", call. = FALSE)
  }

  if (is.null(xlim)) {
    x_at <- pretty(c(-bar_max, bar_max))
    xlim <- range(x_at)
  } else {
    x_at <- pretty(xlim)
  }
  if (is.null(ylim)) {
    ylim <- range(breaks)
  }

  if (is.null(xlab)) {
    xlab <- switch(
      scale,
      count = "Frequency",
      proportion = "Proportion",
      density = "Density"
    )
  }
  if (is.null(ylab)) {
    ylab <- feat
  }

  # Only `mar` is restored, not every settable parameter. A blanket
  # par(no.readonly = TRUE) restore would also put `usr` back, resetting the
  # coordinate system and making it impossible to add anything to the finished
  # plot.
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  graphics::par(mar = margin)

  graphics::plot.default(
    x = NA_real_,
    y = NA_real_,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    xaxt = "n",
    yaxt = "n",
    cex.lab = cex.lab,
    cex.main = cex.main,
    ...
  )

  graphics::rect(
    xleft = -bar_length[[1L]],
    ybottom = breaks[-length(breaks)],
    xright = 0,
    ytop = breaks[-1L],
    col = col[1L],
    border = border
  )
  graphics::rect(
    xleft = 0,
    ybottom = breaks[-length(breaks)],
    xright = bar_length[[2L]],
    ytop = breaks[-1L],
    col = col[2L],
    border = border
  )
  graphics::abline(v = 0)

  x_labels <- abs(x_at)
  if (scale == "proportion") {
    x_labels <- formatC(x_labels, format = "f", digits = 2)
  }
  graphics::axis(side = 1, at = x_at, labels = x_labels, las = 1,
                 cex.axis = cex.axis)
  graphics::axis(side = 2, las = 1, cex.axis = cex.axis)

  if (!is.null(legend.position) && !identical(legend.position, FALSE)) {
    graphics::legend(
      legend.position,
      legend = group_lv,
      fill = col,
      border = border,
      bty = "n",
      cex = cex.legend
    )
  }

  if (!out_statistics) {
    return(invisible(NULL))
  }

  bin_stats <- data.frame(
    bin_start = breaks[-length(breaks)],
    bin_end = breaks[-1L],
    bin_mid = (breaks[-length(breaks)] + breaks[-1L]) / 2
  )
  bin_stats[group_lv] <- bar_length

  group_stats <- as.data.frame(
    vapply(group_lv, function(lv) {
      c(n = n_used[[lv]],
        n_dropped = n_input[[lv]] - n_used[[lv]],
        min = min(values[[lv]]),
        max = max(values[[lv]]))
    }, numeric(4)),
    check.names = FALSE
  )

  invisible(list(
    bin_summary_stats = bin_stats,
    group_summary_stats = group_stats,
    group_hists = hists
  ))
}
