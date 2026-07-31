#' Draw a volcano plot
#'
#' Plots `log2fc` against `-log10(pvalue)` for the table
#' [estimate_significance()] returns, colours the features that clear both
#' cutoffs by direction and labels the strongest of them.
#'
#' @param x The data.frame returned by [estimate_significance()].
#' @param use_adjusted Logical. If `TRUE`, plot and threshold the `adj_pvalue`
#'   column; if `FALSE`, use the unadjusted `pvalue` column. The axis label
#'   follows, so the y axis always describes what was actually plotted.
#' @param log2fc_cutoff,pval_cutoff Cutoffs for calling a feature changed and
#'   significant, drawn as guides. `NULL`, the default, takes the values
#'   [estimate_significance()] recorded on `x`, so the guides agree with its
#'   `is_signif` column. Supply a number to override them. Selecting a subset of
#'   rows keeps the recorded values, but selecting columns drops them, in which
#'   case both must be given.
#' @param anno_feats Logical. If `TRUE`, label the strongest significant
#'   features.
#' @param anno_top How many features to label in each direction, so up to
#'   `2 * anno_top` labels in total.
#' @param cex.anno Character expansion for those labels.
#' @param xlim,ylim Numeric length-2 axis ranges, or `NULL` to derive them from
#'   the data. A supplied range is used as given.
#' @param main Plot title.
#' @param cex.lab,cex.axis,cex.main Character expansion for the axis labels,
#'   the axis annotation and the title.
#' @param margin Plot margins in lines, passed to [graphics::par()] as `mar`.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return `NULL` invisibly.
#'
#' @details
#' The plot margins are restored on exit, but the coordinate system is
#' deliberately left in place, so [graphics::points()], [graphics::text()] and
#' friends can still be used to add to the finished plot.
#'
#' A p-value of exactly zero has no finite `-log10()`. Rather than dropping the
#' most significant points or letting an infinite axis limit blank the plot, the
#' axis is scaled to the largest finite value and those points are drawn at the
#' top of it. Non-finite `log2fc`, which [compare_two_groups()] produces when a
#' group centre is zero, is capped at the edge of the x axis the same way. Either
#' kind of capping is reported in a `message()`, since a capped point no longer
#' sits at its true coordinate.
#'
#' Points are coloured by the same masks that select the labels, so which
#' features are highlighted and which are labelled can never disagree. With the
#' default arguments those masks reproduce the `is_signif` column of `x`. The
#' labels are drawn in a brighter shade than the points on purpose, so that a
#' label stays legible where it overlaps them.
#'
#' @seealso [estimate_significance()], whose output is the only argument this
#'   function needs.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#' res <- compare_two_groups(
#'   data     = iris2,
#'   feats    = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   group    = iris2$Species,
#'   group_lv = c("virginica", "versicolor")
#' )
#' sig <- estimate_significance(res, log2fc_cutoff = 0.1)
#'
#' draw_volcano_plot(sig, main = "virginica vs versicolor")
#'
#' ## Unadjusted p-values, no labels
#' draw_volcano_plot(sig, use_adjusted = FALSE, anno_feats = FALSE)
#'
#' ## A cutoff other than the one the verdict used
#' draw_volcano_plot(sig, log2fc_cutoff = 0.3)
#'
#' @export
draw_volcano_plot <- function(x,
                              use_adjusted = TRUE,
                              log2fc_cutoff = NULL,
                              pval_cutoff = NULL,
                              anno_feats = TRUE,
                              anno_top = 10,
                              cex.anno = 1,
                              xlim = NULL,
                              ylim = NULL,
                              main = NULL,
                              cex.lab = 1.3,
                              cex.axis = 1.2,
                              cex.main = 1.3,
                              margin = c(5, 5, 4, 3),
                              ...) {

  sa_check_flag(use_adjusted, "use_adjusted")
  sa_check_flag(anno_feats, "anno_feats")
  sa_check_scalar_num(anno_top, "anno_top", 0)
  sa_check_scalar_num(cex.anno, "cex.anno", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.lab, "cex.lab", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.axis, "cex.axis", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.main, "cex.main", 0, lower_open = TRUE)
  sa_check_margin(margin)
  sa_check_lim(xlim, "xlim")
  sa_check_lim(ylim, "ylim")

  if (!is.data.frame(x)) {
    stop("`x` must be the data.frame returned by estimate_significance().",
         call. = FALSE)
  }
  p_col <- if (use_adjusted) "adj_pvalue" else "pvalue"
  absent <- setdiff(c("features", "log2fc", p_col), names(x))
  if (length(absent) > 0L) {
    stop("`x` is missing the column(s) ", paste(absent, collapse = ", "),
         ". Pass the table returned by estimate_significance().",
         call. = FALSE)
  }

  # Falling back to the recorded cutoffs is what keeps the guides on the plot and
  # the verdict in the table describing the same rule.
  if (is.null(log2fc_cutoff)) {
    log2fc_cutoff <- attr(x, "log2fc_cutoff")
  }
  if (is.null(pval_cutoff)) {
    pval_cutoff <- attr(x, "pval_cutoff")
  }
  if (is.null(log2fc_cutoff) || is.null(pval_cutoff)) {
    stop("`x` does not carry the cutoffs estimate_significance() records, so ",
         "`log2fc_cutoff` and `pval_cutoff` must be supplied. Selecting ",
         "columns from the table drops them.", call. = FALSE)
  }
  sa_check_scalar_num(log2fc_cutoff, "log2fc_cutoff", 0)
  sa_check_scalar_num(pval_cutoff, "pval_cutoff", 0, 1, lower_open = TRUE)

  feats <- as.character(x$features)
  sa_check_feat_names(feats)
  log2fc <- x$log2fc
  pvalue <- x[[p_col]]
  sa_check_pvalues(pvalue, p_col)

  up_color   <- "#D73027"  # red
  down_color <- "#4575B4"  # blue
  ns_color   <- "grey70"
  # Labels sit on top of the points, so they use a purer and brighter shade than
  # the points do rather than the same one, which would blend in.
  up_label   <- "red"
  down_label <- "blue"

  # A p-value of 0 gives -log10(p) = Inf, which used to make the axis limit
  # infinite and blank the plot. The axis follows the finite values and the
  # infinite ones are drawn at the top of it instead of being discarded.
  neglog_p <- -log10(pvalue)
  y_finite <- neglog_p[is.finite(neglog_p)]
  x_finite <- log2fc[is.finite(log2fc)]
  if (length(y_finite) == 0L || length(x_finite) == 0L) {
    stop("nothing can be plotted: no feature has both a finite `log2fc` and a ",
         "finite -log10(`", p_col, "`).", call. = FALSE)
  }

  y_top <- max(y_finite)
  # A run where every p-value is 1 leaves y_top at 0, which is not a usable
  # axis, so the guide line height sets the floor.
  y_top <- max(y_top, -log10(pval_cutoff), 1)
  x_max <- max(max(abs(x_finite)), log2fc_cutoff)

  if (is.null(ylim)) {
    ylim <- c(0, y_top * 1.1)
  }
  if (is.null(xlim)) {
    xlim <- c(-x_max, x_max) * 1.05
  }
  label_offset <- diff(ylim) * 0.05

  # Capped points are pulled just inside the panel rather than onto its edge, so
  # that a label still fits above them.
  plot_y <- neglog_p
  plot_x <- log2fc
  inf_y <- is.infinite(plot_y)
  inf_x <- is.infinite(plot_x)
  n_capped_y <- sum(inf_y)
  n_capped_x <- sum(inf_x)
  plot_y[inf_y] <- max(ylim) - 2 * label_offset
  plot_x[inf_x] <- ifelse(plot_x[inf_x] > 0,
                          max(xlim) - diff(xlim) * 0.02,
                          min(xlim) + diff(xlim) * 0.02)
  if (n_capped_y > 0L || n_capped_x > 0L) {
    message("Drew ", max(n_capped_y, n_capped_x),
            " point(s) at the edge of the plot",
            if (n_capped_y > 0L) paste0(" (", n_capped_y, " with p = 0)"),
            if (n_capped_x > 0L) {
              paste0(" (", n_capped_x, " with an infinite log2 fold change)")
            },
            "; their true position is off the axis.")
  }

  # One mask per direction, shared by the points and the labels so that a point
  # can never be coloured as changed while its label is left out, or vice versa.
  passes_p <- !is.na(pvalue) & pvalue <= pval_cutoff
  is_up <- passes_p & !is.na(log2fc) & log2fc >= log2fc_cutoff
  is_down <- passes_p & !is.na(log2fc) & log2fc <= -log2fc_cutoff

  # Only `mar` is restored, not every settable parameter. A blanket
  # par(no.readonly = TRUE) restore would also put `usr` back, resetting the
  # coordinate system and making it impossible to add anything to the finished
  # plot.
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  graphics::par(mar = margin)

  y_lab <- if (use_adjusted) {
    expression(-log[10] ~ adjusted ~ italic(P))
  } else {
    expression(-log[10] ~ italic(P))
  }

  graphics::plot(
    plot_x,
    plot_y,
    xlab = expression(log[2] ~ FC),
    ylab = y_lab,
    xlim = xlim,
    ylim = ylim,
    main = main,
    cex.lab = cex.lab,
    cex.axis = cex.axis,
    cex.main = cex.main,
    pch = 16,
    col = ns_color,
    ...
  )
  graphics::abline(
    h = -log10(pval_cutoff),
    v = c(-log2fc_cutoff, log2fc_cutoff),
    col = "green3", lwd = 2, lty = 3
  )

  graphics::points(plot_x[is_up], plot_y[is_up], col = up_color, pch = 16)
  graphics::points(plot_x[is_down], plot_y[is_down], col = down_color, pch = 16)

  if (anno_feats && anno_top >= 1) {
    # Strongest first means smallest p-value, then largest fold change away from
    # zero in the direction concerned.
    pick <- function(mask, decreasing) {
      i <- which(mask)
      i <- i[order(pvalue[i], if (decreasing) -log2fc[i] else log2fc[i])]
      i[seq_len(min(anno_top, length(i)))]
    }
    up_idx <- pick(is_up, TRUE)
    down_idx <- pick(is_down, FALSE)

    graphics::text(
      c(plot_x[up_idx], plot_x[down_idx]),
      c(plot_y[up_idx], plot_y[down_idx]) + label_offset,
      labels = c(feats[up_idx], feats[down_idx]),
      cex = cex.anno,
      col = c(rep(up_label, length(up_idx)), rep(down_label, length(down_idx)))
    )
  }

  invisible(NULL)
}
