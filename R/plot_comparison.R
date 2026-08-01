#' Plot a comparison result
#'
#' Draws whichever table of a comparison the caller names, using only the
#' columns the result contract guarantees. That is why one method covers every
#' scenario: it never asks whether the object came from two groups, three
#' groups or a single sample, only whether the table it was handed has
#' intervals or p-values, and every table has one or the other.
#'
#' Three views are available:
#'
#' \describe{
#'   \item{`"estimate"`}{A forest plot of the effect estimate and its
#'     confidence interval, one row per feature. Available whenever the chosen
#'     table has finite intervals, which the two-group and one-sample scenarios
#'     always do.}
#'   \item{`"posthoc"`}{The same forest plot for the pairwise contrasts of a
#'     multi-group comparison, one row per contrast. With more than one feature
#'     in the post-hoc table, name the one to draw with `feature`.}
#'   \item{`"pvalue"`}{`-log10(pval_adj)` per feature with the `alpha`
#'     threshold marked. The fallback when a table has no interval to draw,
#'     which is the case for every omnibus test.}
#' }
#'
#' `type = "auto"`, the default, picks the first of those three that the chosen
#' table can actually support.
#'
#' @param x A comparison result, as returned by [compare_two_groups()],
#'   [compare_multiple_groups()] or [compare_one_sample()].
#' @param test Which test to draw. One of `names(x$tests)`.
#' @param type `"auto"`, `"estimate"`, `"posthoc"` or `"pvalue"`.
#' @param feature For `type = "posthoc"`, which feature's contrasts to draw.
#'   Defaults to the first feature present in the post-hoc table.
#' @param alpha Threshold marked on the p-value view and used to colour the
#'   points of the estimate view.
#' @param sort_by `"none"` to keep the feature order of the result, or
#'   `"pvalue"` to draw the most significant rows at the top.
#' @param dark If `TRUE`, use a dark background with light text.
#' @param xlab,main Axis and title labels. Both are derived from the result when
#'   left `NULL`.
#' @param col_signif,col_plain Colours for rows at or below `alpha` and for the
#'   rest.
#' @param cex.axis,cex.lab,cex.main Character expansion for the axis
#'   annotation, the axis label and the title.
#' @param ... Ignored, present for consistency with [plot()].
#'
#' @return The plotted data.frame, invisibly, in the row order it was drawn,
#'   with the view that `type = "auto"` resolved to attached as the attribute
#'   `"view"`.
#'
#' @details
#' The function changes graphical parameters and restores them on exit, so the
#' caller's device is left as it was found.
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
#' plot(res)
#'
#' ## The same comparison read through the rank-based test
#' plot(res, test = "wilcox_test", sort_by = "pvalue")
#'
#' ## A multi-group omnibus table has no interval, so "auto" falls through to
#' ## the pairwise contrasts.
#' multi <- compare_multiple_groups(iris, feats, iris$Species,
#'                                  levels(iris$Species))
#' plot(multi, feature = "Petal.Length")
#' plot(multi, type = "pvalue")
#'
#' @export
plot.sa_comparison <- function(x,
                               test = names(x$tests)[1],
                               type = c("auto", "estimate", "posthoc",
                                        "pvalue"),
                               feature = NULL,
                               alpha = 0.05,
                               sort_by = c("none", "pvalue"),
                               dark = FALSE,
                               xlab = NULL,
                               main = NULL,
                               col_signif = "#D1495B",
                               col_plain = "#7F8C8D",
                               cex.axis = 0.9,
                               cex.lab = 1.1,
                               cex.main = 1.2,
                               ...) {

  type <- match.arg(type)
  sort_by <- match.arg(sort_by)
  sa_check_flag(dark, "dark")
  sa_check_scalar_num(alpha, "alpha", 0, 1, lower_open = TRUE)

  tbl <- sa_pick_test(x, test, arg = "x")
  posthoc <- x$posthoc[[test]]
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
      stop("`x$posthoc$", test, "` is empty, so there are no contrasts to ",
           "draw. Only compare_multiple_groups() produces a post-hoc stage, ",
           "and only for features whose omnibus test qualified.", call. = FALSE)
    }
    if (is.null(feature)) {
      feature <- posthoc$features[1]
    }
    if (!feature %in% posthoc$features) {
      stop("`feature` must name a feature present in the post-hoc table: ",
           paste(unique(posthoc$features), collapse = ", "), ". Got ",
           feature, ".", call. = FALSE)
    }
    drawn <- posthoc[posthoc$features == feature, , drop = FALSE]
    labels <- drawn$contrast
    estimate <- drawn$estimate
    null_value <- 0
    default_xlab <- "estimate (group1 - group2)"
    default_main <- paste0(feature, ": ", x$test_info[[test]]$posthoc_label)
  } else if (type == "estimate") {
    column <- sa_estimate_column(tbl)
    if (is.null(column)) {
      stop("`x$tests$", test, "` holds no estimate to draw. An omnibus test ",
           "reports that the levels differ, not by how much; use ",
           "type = \"posthoc\" for the contrasts or type = \"pvalue\".",
           call. = FALSE)
    }
    drawn <- tbl
    labels <- drawn$features
    estimate <- drawn[[column]]
    # The relative effect is a probability whose null is 0.5, not a difference
    # whose null is 0. Drawing the guide at 0 would put every point on one side.
    null_value <- if (column == "relative_effect") 0.5 else 0
    default_xlab <- column
    default_main <- x$test_info[[test]]$label
  } else {
    drawn <- tbl
    labels <- drawn$features
    estimate <- NULL
    null_value <- NULL
    default_xlab <- expression(-log[10](p[adj]))
    default_main <- x$test_info[[test]]$label
  }

  if (sort_by == "pvalue") {
    order_by <- order(drawn$pval_adj, na.last = TRUE)
    drawn <- drawn[order_by, , drop = FALSE]
    labels <- labels[order_by]
    estimate <- estimate[order_by]
  }

  theme <- sa_plot_theme(dark)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(bg = theme$bg, fg = theme$fg, col.axis = theme$fg,
                col.lab = theme$fg, col.main = theme$fg,
                mar = c(5.1, max(4.1, 0.6 * max(nchar(labels))), 4.1, 2.1))

  signif <- !is.na(drawn$pval_adj) & drawn$pval_adj <= alpha
  cols <- ifelse(signif, col_signif, col_plain)
  # Row 1 at the top, so a sorted plot reads downwards like the table it came
  # from rather than upwards like a default plot.
  at <- rev(seq_len(nrow(drawn)))

  if (type == "pvalue") {
    height <- -log10(drawn$pval_adj)
    height[!is.finite(height)] <- 0
    graphics::plot.default(
      height, at, type = "n", yaxt = "n", bty = "n",
      xlim = c(0, max(c(height, -log10(alpha)), na.rm = TRUE) * 1.1),
      ylim = c(0.5, nrow(drawn) + 0.5),
      xlab = if (is.null(xlab)) default_xlab else xlab, ylab = "",
      main = if (is.null(main)) default_main else main,
      cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
    )
    graphics::abline(v = -log10(alpha), lty = 2, col = theme$guide)
    graphics::segments(0, at, height, at, col = cols, lwd = 6, lend = 1)
  } else {
    span <- range(c(drawn$lower_conf, drawn$upper_conf, estimate, null_value),
                  na.rm = TRUE, finite = TRUE)
    graphics::plot.default(
      estimate, at, type = "n", yaxt = "n", bty = "n",
      xlim = span + c(-1, 1) * diff(span) * 0.08,
      ylim = c(0.5, nrow(drawn) + 0.5),
      xlab = if (is.null(xlab)) default_xlab else xlab, ylab = "",
      main = if (is.null(main)) default_main else main,
      cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
    )
    graphics::abline(v = null_value, lty = 2, col = theme$guide)
    # Bounds are clamped to the panel so that a one-sided interval, which runs
    # to infinity on the side it does not test, still draws as a line reaching
    # the edge instead of vanishing.
    lower <- pmax(drawn$lower_conf, span[1] - diff(span) * 0.08)
    upper <- pmin(drawn$upper_conf, span[2] + diff(span) * 0.08)
    graphics::segments(lower, at, upper, at, col = cols, lwd = 2)
    graphics::points(estimate, at, pch = 18, cex = 1.6, col = cols)
  }

  graphics::axis(2, at = at, labels = labels, las = 1, tick = FALSE,
                 cex.axis = cex.axis)
  graphics::legend("topright", legend = c(paste0("p_adj <= ", alpha),
                                          paste0("p_adj > ", alpha)),
                   col = c(col_signif, col_plain), pch = 15, bty = "n",
                   cex = cex.axis, text.col = theme$fg)

  # The view is carried on the result because `type = "auto"` resolves it here
  # and the caller would otherwise have no way to find out which of the three
  # it got, short of reading the axis label off the device.
  attr(drawn, "view") <- type
  invisible(drawn)
}


#' Which column of a test table holds the estimate to draw
#'
#' The tests report their estimates under their own names, so the plot method
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
