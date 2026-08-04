#' Draw a grouped boxplot across several features
#'
#' Draws one cluster of boxes per feature, with the group levels side by side
#' inside each cluster and a legend in a narrow panel on the right. Optionally
#' returns the summary statistics behind the boxes.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to plot,
#'   in display order along the x axis.
#' @param group Grouping vector with one entry per row of `data`.
#' @param group_lv Character vector of at least two group levels, in the order
#'   they should appear inside each cluster. Rows belonging to any other level
#'   are dropped.
#' @param gap Blank box widths inserted between neighbouring feature clusters.
#' @param lwd Line width of the boxes.
#' @param xlab,ylab,main Axis and title labels.
#' @param cex.lab,cex.axis,cex.main,cex.legend Character expansion for axis
#'   labels, axis annotation, the main title and the legend.
#' @param ylim Numeric length-2 y axis range, or `NULL` to let
#'   [graphics::boxplot()] choose.
#' @param dark If `TRUE`, use a dark background with light text.
#' @param grid_lty,grid_lwd Line type and width of the horizontal grid.
#' @param out_statistics If `TRUE`, invisibly return the summary statistics
#'   used by the plot.
#'
#' @return If `out_statistics = FALSE`, `NULL` invisibly. Otherwise a list of
#'   two elements, invisibly:
#'
#'   \describe{
#'     \item{`box_summary_stats`}{One data.frame per feature with rows `min`,
#'       `lower_bound`, `Q1`, `median`, `Q3`, `upper_bound`, `max` and one
#'       column per group level. The bounds are the Tukey whisker fences
#'       `Q1 - 1.5 * IQR` and `Q3 + 1.5 * IQR`, not the drawn whisker ends.}
#'     \item{`median_confidence_stats`}{One data.frame per feature with rows
#'       `n`, `lower_conf`, `upper_conf`, the notch interval
#'       `median +/- 1.58 * IQR / sqrt(n)`. `n` counts non-missing values.}
#'   }
#'
#' @details
#' The function changes graphical parameters and the panel layout, and restores
#' both on exit, so the caller's device is left as it was found.
#'
#' @seealso [compare_two_groups()] to test the same input.
#'
#' @references
#' McGill, R., Tukey, J. W. and Larsen, W. A. (1978). Variations of box plots.
#' *The American Statistician*, 32(1), 12-16.
#'
#' @examples
#' ## All three iris species across the four measurements
#' stats <- draw_grouped_boxplot(
#'   data     = iris,
#'   feats    = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   group    = iris$Species,
#'   group_lv = c("setosa", "versicolor", "virginica"),
#'   ylab     = "cm",
#'   dark     = FALSE
#' )
#' stats$box_summary_stats$Sepal.Length
#' stats$median_confidence_stats$Sepal.Length
#'
#' @export
draw_grouped_boxplot <- function(data,
                                 feats,
                                 group,
                                 group_lv,
                                 gap = 1,
                                 lwd = 1.5,
                                 xlab = NULL,
                                 ylab = NULL,
                                 cex.lab = 1.3,
                                 cex.axis = 1.2,
                                 cex.main = 1.3,
                                 ylim = NULL,
                                 main = NULL,
                                 dark = FALSE,
                                 grid_lty = 1,
                                 grid_lwd = 0.25,
                                 cex.legend = 1.1,
                                 out_statistics = TRUE) {

  sa_check_scalar_num(gap, "gap", 0)
  sa_check_scalar_num(lwd, "lwd", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.lab, "cex.lab", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.axis, "cex.axis", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.main, "cex.main", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.legend, "cex.legend", 0, lower_open = TRUE)
  sa_check_scalar_num(grid_lwd, "grid_lwd", 0)
  sa_check_flag(dark, "dark")
  sa_check_flag(out_statistics, "out_statistics")
  if (!is.null(ylim) && (!is.numeric(ylim) || length(ylim) != 2L)) {
    stop("`ylim` must be NULL or a numeric vector of length 2.", call. = FALSE)
  }

  input <- sa_validate_wide_input(data, feats, group, group_lv,
                                  min_levels = 2L)
  data <- input$data
  feats <- input$feats
  group <- input$group
  group_lv <- levels(group)
  n_feats <- length(feats)
  n_lv <- length(group_lv)

  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }

  # One vector per feature and group level, the group levels consecutive within
  # each feature so the order matches the layout of `at`. Handing boxplot() the
  # list directly replaces a reshape to long format, which was the only reason
  # this package ever needed the `.grp`, `name` and `value` column names.
  boxes <- unlist(
    lapply(feats, function(f) {
      lapply(group_lv, function(lv) data[[f]][group == lv])
    }),
    recursive = FALSE,
    use.names = FALSE
  )

  at <- lapply(seq_len(n_feats), function(i) {
    start <- (i - 1) * (n_lv + gap)
    start + seq_len(n_lv)
  })

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  graphics::layout(
    matrix(c(1, 2), nrow = 1),
    widths = c(4, 1)
  )

  if (dark) {
    graphics::par(
      bg = "#2B2B2B",      # figure background
      fg = "white",        # default foreground
      col.axis = "white",
      col.lab = "white",
      col.main = "white"
    )
    bg_cols <- "#36454F"
    grid_col <- "gray80"
  } else {
    bg_cols <- "white"
    grid_col <- "gray40"
  }

  cols <- grDevices::hcl.colors(n_lv, "Dark 2")

  # Drawn twice with the grid in between so the grid sits behind the boxes. One
  # closure keeps the two calls from drifting apart.
  draw_boxes <- function(add) {
    graphics::boxplot(
      boxes,
      at       = unlist(at),
      xlab     = xlab,
      ylab     = ylab,
      xaxt     = "n",
      cex.lab  = cex.lab,
      cex.axis = cex.axis,
      cex.main = cex.main,
      ylim     = ylim,
      main     = main,
      frame    = FALSE,
      col      = bg_cols,
      border   = cols,
      lwd      = lwd,
      pch      = 16,
      add      = add
    )
  }

  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  draw_boxes(add = FALSE)

  graphics::grid(
    nx = NA,             # no vertical grid
    ny = NULL,           # horizontal grid at the y axis ticks
    col = grid_col,
    lty = grid_lty,
    lwd = grid_lwd
  )

  draw_boxes(add = TRUE)

  graphics::axis(
    side = 1,
    at = vapply(at, mean, numeric(1)),
    labels = feats,
    cex.axis = cex.axis,
    tick = FALSE
  )

  graphics::par(mar = c(5, 0, 4, 1))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = group_lv,
    fill = cols,
    border = "white",
    bty = "n",
    cex = cex.legend
  )

  if (!out_statistics) {
    return(invisible(NULL))
  }

  # One pass per feature/level; the two returned tables are slices of it.
  summaries <- lapply(feats, function(f) {
    vapply(group_lv, function(lv) {
      v <- data[[f]][group == lv]
      v <- v[!is.na(v)]
      n <- length(v)
      if (n == 0L) {
        return(sa_na_row(c("min", "lower_bound", "Q1", "median", "Q3",
                           "upper_bound", "max", "n", "lower_conf",
                           "upper_conf")))
      }
      q <- unname(stats::quantile(v, c(0.25, 0.5, 0.75)))
      iqr <- q[3] - q[1]
      notch <- 1.58 * iqr / sqrt(n)
      c(min         = min(v),
        lower_bound = q[1] - 1.5 * iqr,
        Q1          = q[1],
        median      = q[2],
        Q3          = q[3],
        upper_bound = q[3] + 1.5 * iqr,
        max         = max(v),
        n           = n,
        lower_conf  = q[2] - notch,
        upper_conf  = q[2] + notch)
    }, numeric(10))
  })
  names(summaries) <- feats

  box_rows <- c("min", "lower_bound", "Q1", "median", "Q3", "upper_bound",
                "max")
  conf_rows <- c("n", "lower_conf", "upper_conf")

  invisible(list(
    box_summary_stats = lapply(summaries, function(m) {
      as.data.frame(m[box_rows, , drop = FALSE])
    }),
    median_confidence_stats = lapply(summaries, function(m) {
      as.data.frame(m[conf_rows, , drop = FALSE])
    })
  ))
}
