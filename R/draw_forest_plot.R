#' Draw a forest plot of a comparison result
#'
#' Draws the estimate of each feature beside its confidence interval, falling
#' back to a bar of `-log10(pval_adj)` for a table that has no interval to
#' draw, which is every omnibus test.
#'
#' Only the columns the result contract guarantees are read, which is why one
#' function covers every scenario: it never asks whether the object came from
#' two groups, three groups or a single sample, only whether the table it was
#' handed has intervals or p-values, and every table has one or the other.
#' `plot()` on an `sa_comparison` is the same function under the name R users
#' reach for first.
#'
#' Three views are available:
#'
#' \describe{
#'   \item{`"estimate"`}{A forest plot of the effect estimate and its
#'     confidence interval, one row per feature. Available whenever the chosen
#'     table has finite intervals, which the two-group and one-sample scenarios
#'     always do.}
#'   \item{`"posthoc"`}{The same forest plot for the pairwise contrasts of a
#'     multi-group comparison, one row per contrast. Only the first feature of
#'     the post-hoc table is drawn unless `feats` names the ones to draw. An
#'     estimate reads as `group1 - group2`, the direction the row label spells
#'     out, and the reference level is the one subtracted, so a point to the
#'     right of the guide agrees in sign with the `log2fc` a volcano plot of the
#'     same comparison draws.}
#'   \item{`"pvalue"`}{`-log10()` of the p-value per feature with the `alpha`
#'     threshold marked. The fallback when a table has no interval to draw,
#'     which is the case for every omnibus test.}
#' }
#'
#' `type = "auto"`, the default, picks the first of those three that the chosen
#' table can actually support.
#'
#' @param comparison_result A comparison result, as returned by
#'   [compare_two_groups()], [compare_multiple_groups()] or
#'   [compare_one_sample()].
#' @param x The same comparison result, under the name the [plot()] generic
#'   fixes for its first argument.
#' @param test Which test to draw. One of `names(comparison_result$tests)`.
#' @param type `"auto"`, `"estimate"`, `"posthoc"` or `"pvalue"`.
#' @param feats Character vector of features to draw, in display order from the
#'   top of the plot down. `NULL`, the default, draws every feature of the
#'   chosen table, except in the post-hoc view, where it draws the first feature
#'   of the post-hoc table. A feature that has no contrasts because its omnibus
#'   test did not qualify is reported in a `message()` and left out.
#' @param use_adjusted Logical. If `TRUE`, read the `pval_adj` column; if
#'   `FALSE`, the unadjusted `pval`. The colouring, the sorting, the p-value
#'   view and the labels all follow, so the plot always describes the p-value it
#'   actually used.
#' @param alpha Threshold marked on the p-value view and used to colour the
#'   points of the estimate view.
#' @param sort_by `"none"` to keep the feature order of the result, or
#'   `"pvalue"` to draw the most significant rows at the top.
#' @param dark If `TRUE`, use a dark background with light text.
#' @param xlim Numeric length-2 x axis range, or `NULL` to derive it from the
#'   values being drawn. A supplied range is used as given, and the interval
#'   bounds are clamped to it, so an interval running past the range still
#'   reaches the edge of the panel instead of vanishing.
#' @param xlab,main Axis and title labels. Both are derived from the result when
#'   left `NULL`.
#' @param col_signif,col_plain Colours for rows at or below `alpha` and for the
#'   rest.
#' @param cex.axis,cex.lab,cex.main,cex.legend Character expansion for the axis
#'   annotation, the axis label, the title and the legend.
#' @param ... Arguments the `plot()` method passes on. Anything
#'   `draw_forest_plot()` does not name is ignored.
#'
#' @return The plotted data.frame, invisibly, in the row order it was drawn,
#'   with the view that `type = "auto"` resolved to attached as the attribute
#'   `"view"`.
#'
#' @details
#' The function changes graphical parameters and the panel layout, and restores
#' both on exit, so the caller's device is left as it was found. The legend
#' lives in a narrow panel of its own on the right, where it cannot cover the
#' rows it describes.
#'
#' The estimate view marks the null value with a vertical line, at zero for a
#' difference and at one for the Brunner-Munzel relative effect, which is the
#' one quantity in the package whose null is not zero.
#'
#' @seealso [draw_volcano_plot()], which plots effect size against significance
#'   rather than estimate against interval, and [draw_grouped_boxplot()] for the
#'   input data.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#' feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#' res <- compare_two_groups(iris2, feats, iris2$Species,
#'                           c("virginica", "versicolor"))
#'
#' ## Mean differences with their intervals
#' draw_forest_plot(res)
#'
#' ## plot() reaches the same function
#' plot(res)
#'
#' ## The same comparison read through the rank-based test
#' draw_forest_plot(res, test = "wilcox_test", sort_by = "pvalue")
#'
#' ## Two features only, drawn in the order they are named, against the
#' ## unadjusted p-value
#' draw_forest_plot(res, feats = c("Petal.Width", "Petal.Length"),
#'                  use_adjusted = FALSE)
#'
#' ## An axis range of your own, so that two plots can be read against each other
#' draw_forest_plot(res, xlim = c(0, 2))
#'
#' ## A multi-group omnibus table has no interval, so "auto" falls through to
#' ## the pairwise contrasts.
#' multi <- compare_multiple_groups(iris, feats, iris$Species,
#'                                  levels(iris$Species))
#' draw_forest_plot(multi, feats = "Petal.Length")
#'
#' ## Several features at once label each contrast with the feature it belongs to
#' draw_forest_plot(multi, feats = c("Petal.Length", "Sepal.Width"))
#' draw_forest_plot(multi, type = "pvalue")
#'
#' @export
draw_forest_plot <- function(comparison_result,
                             test = names(comparison_result$tests)[1],
                             type = c("auto", "estimate", "posthoc",
                                      "pvalue"),
                             feats = NULL,
                             use_adjusted = TRUE,
                             alpha = 0.05,
                             sort_by = c("none", "pvalue"),
                             dark = FALSE,
                             xlim = NULL,
                             xlab = NULL,
                             main = NULL,
                             col_signif = "#D1495B",
                             col_plain = "#7F8C8D",
                             cex.axis = 0.9,
                             cex.lab = 1.1,
                             cex.main = 1.2,
                             cex.legend = 1,
                             ...) {

  type <- match.arg(type)
  sort_by <- match.arg(sort_by)
  sa_check_flag(dark, "dark")
  sa_check_flag(use_adjusted, "use_adjusted")
  sa_check_scalar_num(alpha, "alpha", 0, 1, lower_open = TRUE)
  sa_check_scalar_num(cex.legend, "cex.legend", 0, lower_open = TRUE)
  sa_check_lim(xlim, "xlim")

  tbl <- sa_pick_test(comparison_result, test, arg = "comparison_result")
  posthoc <- comparison_result$posthoc[[test]]
  p_col <- if (use_adjusted) "pval_adj" else "pval"
  had_posthoc <- !is.null(posthoc) && nrow(posthoc) > 0L
  not_qualified <- character(0)

  if (!is.null(feats)) {
    sa_check_feat_names(feats)
    unknown <- setdiff(feats, tbl$features)
    if (length(unknown) > 0L) {
      stop("`feats` must name features present in the comparison: ",
           paste(tbl$features, collapse = ", "), ". Not found: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    # Selected before `type = "auto"` resolves, so the view is chosen from what
    # will actually be drawn rather than from the whole table. Picking three
    # features whose intervals are all NA falls through to the p-value view
    # instead of drawing an empty panel.
    tbl <- tbl[match(feats, tbl$features), , drop = FALSE]
    if (had_posthoc) {
      not_qualified <- setdiff(feats, posthoc$features)
      posthoc <- posthoc[posthoc$features %in% feats, , drop = FALSE]
      posthoc <- posthoc[order(match(posthoc$features, feats)), , drop = FALSE]
    }
  }

  has_estimate <- !is.null(sa_estimate_column(tbl)) &&
    any(is.finite(tbl$lower_conf) & is.finite(tbl$upper_conf))
  has_posthoc <- !is.null(posthoc) && nrow(posthoc) > 0L

  if (type == "auto") {
    type <- if (has_estimate) {
      "estimate"
    } else if (has_posthoc) {
      "posthoc"
    } else {
      "pvalue"
    }
  }

  if (type == "posthoc") {
    if (!has_posthoc) {
      if (had_posthoc) {
        stop("none of the features named in `feats` has contrasts to draw: ",
             paste(not_qualified, collapse = ", "), ". The pairwise stage ",
             "runs only for the features whose omnibus test qualified.",
             call. = FALSE)
      }
      stop("`comparison_result$posthoc$", test, "` holds no contrasts to ",
           "draw. Only compare_multiple_groups() produces a ",
           "post-hoc stage, and only for features whose omnibus test ",
           "qualified.", call. = FALSE)
    }
    if (length(not_qualified) > 0L) {
      message("No contrasts for ", paste(not_qualified, collapse = ", "),
              "; the omnibus test did not qualify them for the post-hoc ",
              "stage.")
    }
    # Every contrast of every feature at once is a wall of rows that says
    # nothing, so an unselected post-hoc view stays on the first feature.
    drawn <- if (is.null(feats)) {
      posthoc[posthoc$features == posthoc$features[1], , drop = FALSE]
    } else {
      posthoc
    }
    many_feats <- length(unique(drawn$features)) > 1L
    # A contrast label reads the same for every feature, so it only identifies a
    # row once the plot is down to one feature.
    labels <- if (many_feats) {
      paste0(drawn$features, ": ", drawn$contrast)
    } else {
      drawn$contrast
    }
    estimate <- drawn$estimate
    null_value <- 0
    default_xlab <- "estimate (group1 - group2)"
    posthoc_label <- comparison_result$test_info[[test]]$posthoc_label
    default_main <- if (many_feats) {
      posthoc_label
    } else {
      paste0(drawn$features[1], ": ", posthoc_label)
    }
  } else if (type == "estimate") {
    column <- sa_estimate_column(tbl)
    if (is.null(column)) {
      stop("`comparison_result$tests$", test, "` holds no estimate to draw. ",
           "An omnibus test reports that the levels differ, not by how much; ",
           "use type = \"posthoc\" for the contrasts or type = \"pvalue\".",
           call. = FALSE)
    }
    drawn <- tbl
    labels <- drawn$features
    estimate <- drawn[[column]]
    # The relative effect is a probability whose null is 0.5, not a difference
    # whose null is 0. Drawing the guide at 0 would put every point on one side.
    null_value <- if (column == "relative_effect") 0.5 else 0
    default_xlab <- column
    default_main <- comparison_result$test_info[[test]]$label
  } else {
    drawn <- tbl
    labels <- drawn$features
    estimate <- NULL
    null_value <- NULL
    default_xlab <- if (use_adjusted) {
      expression(-log[10] ~ adjusted ~ italic(P))
    } else {
      expression(-log[10] ~ italic(P))
    }
    default_main <- comparison_result$test_info[[test]]$label
  }

  if (sort_by == "pvalue") {
    order_by <- order(drawn[[p_col]], na.last = TRUE)
    drawn <- drawn[order_by, , drop = FALSE]
    labels <- labels[order_by]
    estimate <- estimate[order_by]
  }

  theme <- sa_plot_theme(dark)
  # Only the parameters this function sets are put back, not a blanket
  # par(no.readonly = TRUE) snapshot. That snapshot also carries `fin`, `pin`
  # and `mai`, which are absolute sizes: restoring them pins the figure to the
  # size this plot happened to be drawn at, so the next plot on a device that
  # has since been resized is redrawn small in a corner of it. `mfrow` comes
  # back with the rest, since `layout()` overwrites whatever grid the caller
  # had set up.
  old_par <- graphics::par(c("bg", "fg", "col.axis", "col.lab", "col.main",
                             "mar", "mfrow"))
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  # The legend gets a panel of its own rather than a corner of the plot, where
  # it sat on top of the rows it was describing.
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))

  # A long label earns a wider left margin, but only up to half of the panel it
  # is labelling. Past that the labels are shrunk to fit the margin they get,
  # which keeps every name readable in full instead of squeezing the plot into
  # a strip too narrow to carry an axis.
  wanted <- max(4.1, 0.6 * max(nchar(labels)))
  panel_lines <- graphics::par("din")[1] * 0.8 / graphics::par("csi")
  left <- min(wanted, 0.5 * panel_lines)
  cex_labels <- cex.axis * min(1, left / wanted)

  graphics::par(bg = theme$bg, fg = theme$fg, col.axis = theme$fg,
                col.lab = theme$fg, col.main = theme$fg,
                mar = c(5.1, left, 4.1, 2.1))

  signif <- !is.na(drawn[[p_col]]) & drawn[[p_col]] <= alpha
  cols <- ifelse(signif, col_signif, col_plain)
  # Row 1 at the top, so a sorted plot reads downwards like the table it came
  # from rather than upwards like a default plot.
  at <- rev(seq_len(nrow(drawn)))

  if (type == "pvalue") {
    height <- -log10(drawn[[p_col]])
    height[!is.finite(height)] <- 0
    limits <- if (is.null(xlim)) {
      c(0, max(c(height, -log10(alpha)), na.rm = TRUE) * 1.1)
    } else {
      xlim
    }
    graphics::plot.default(
      height, at, type = "n", yaxt = "n", xaxt = "n", bty = "n",
      xlim = limits,
      ylim = c(0.5, nrow(drawn) + 0.5),
      xlab = if (is.null(xlab)) default_xlab else xlab, ylab = "",
      main = if (is.null(main)) default_main else main,
      cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
    )
    sa_draw_x_axis(cex.axis)
    graphics::abline(v = -log10(alpha), lty = 2, col = theme$guide)
    graphics::segments(0, at, height, at, col = cols, lwd = 6, lend = 1)
  } else {
    span <- range(c(drawn$lower_conf, drawn$upper_conf, estimate, null_value),
                  na.rm = TRUE, finite = TRUE)
    limits <- if (is.null(xlim)) {
      span + c(-1, 1) * diff(span) * 0.08
    } else {
      xlim
    }
    graphics::plot.default(
      estimate, at, type = "n", yaxt = "n", xaxt = "n", bty = "n",
      xlim = limits,
      ylim = c(0.5, nrow(drawn) + 0.5),
      xlab = if (is.null(xlab)) default_xlab else xlab, ylab = "",
      main = if (is.null(main)) default_main else main,
      cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
    )
    sa_draw_x_axis(cex.axis)
    graphics::abline(v = null_value, lty = 2, col = theme$guide)
    # Bounds are clamped to the panel so that a one-sided interval, which runs
    # to infinity on the side it does not test, still draws as a line reaching
    # the edge instead of vanishing. A narrowed `xlim` cuts the same way.
    lower <- pmax(drawn$lower_conf, min(limits))
    upper <- pmin(drawn$upper_conf, max(limits))
    graphics::segments(lower, at, upper, at, col = cols, lwd = 2)
    graphics::points(estimate, at, pch = 18, cex = 1.6, col = cols)
  }

  graphics::axis(2, at = at, labels = labels, las = 1, tick = FALSE,
                 cex.axis = cex_labels)

  # The threshold moves to the legend title so that the two entries can say what
  # they mean rather than repeat the comparison that produced them.
  legend_title <- if (use_adjusted) {
    as.expression(bquote(italic(P)[adj] <= .(alpha)))
  } else {
    as.expression(bquote(italic(P) <= .(alpha)))
  }
  graphics::par(mar = c(5, 0, 4, 1))
  graphics::plot.new()
  graphics::legend("center", title = legend_title,
                   legend = c("Significant", "Not significant"),
                   col = c(col_signif, col_plain), pch = 15, bty = "n",
                   cex = cex.legend, text.col = theme$fg,
                   title.col = theme$fg)

  # The view is carried on the result because `type = "auto"` resolves it here
  # and the caller would otherwise have no way to find out which of the three
  # it got, short of reading the axis label off the device.
  attr(drawn, "view") <- type
  invisible(drawn)
}


#' @rdname draw_forest_plot
#' @export
plot.sa_comparison <- function(x, ...) {
  draw_forest_plot(x, ...)
}


#' Which column of a test table holds the estimate to draw
#'
#' The tests report their estimates under their own names, so the forest plot
#' looks for the first one it recognises rather than requiring every table to
#' agree on a single column name. A table with none of them, which is every
#' omnibus table, gets `NULL` and falls through to the p-value view.
#'
#' @keywords internal
#' @noRd
sa_estimate_column <- function(tbl) {
  candidates <- c("mean_diff", "hl_shift", "trim_diff", "relative_effect",
                  "diff", "estimate")
  hit <- candidates[candidates %in% names(tbl)]
  if (length(hit) == 0L) NULL else hit[1]
}


#' Draw an x axis whose line spans the whole panel
#'
#' The axis line [graphics::plot.default()] draws stops at the outermost tick,
#' so a point or an interval end lying beyond the last labelled value leaves the
#' axis looking broken off in the middle of the panel. The line and the ticks
#' are drawn in two calls instead, the first carrying the line across the panel
#' with no ticks, the second the ticks and labels with no line.
#'
#' @keywords internal
#' @noRd
sa_draw_x_axis <- function(cex.axis) {
  usr <- graphics::par("usr")
  graphics::axis(1, at = usr[1:2], labels = FALSE, lwd.ticks = 0)
  graphics::axis(1, at = graphics::axTicks(1), lwd = 0, lwd.ticks = 1,
                 cex.axis = cex.axis)
  invisible(NULL)
}


#' Foreground, background and guide colours for a plot
#'
#' @keywords internal
#' @noRd
sa_plot_theme <- function(dark) {
  if (dark) {
    list(bg = "#2B2B2B", fg = "white", guide = "gray70")
  } else {
    list(bg = "white", fg = "black", guide = "gray40")
  }
}
