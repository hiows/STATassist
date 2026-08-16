# Predicted against observed, the picture of a regression evaluation. Two
# lines carry the whole reading: the identity, where a perfect prediction would
# lie, and the calibration line, which is where the predictions actually lie.
# The gap between them is the same thing `$metrics` reports as `calib_slope` and
# `calib_intercept`, so the plot and the table cannot disagree.
#
# Several models are panelled rather than overlaid by default. Overlaying two
# clouds of points that share an axis produces a third cloud that belongs to
# neither, which is the reason the sketch this came from hesitated between
# drawing the points and drawing only the lines. Panelling answers it: each
# model keeps its own points, and the shared limits are what makes the panels
# comparable.


#' Draw one predicted-against-observed panel
#'
#' @keywords internal
#' @noRd
sa_prediction_panel <- function(observed, predicted, row, lim, theme, col,
                                lwd, draw_points, xlab, ylab, main,
                                cex.axis, cex.lab, cex.main) {
  graphics::plot.default(
    observed, predicted, type = "n", bty = "n",
    xlim = lim, ylim = lim, xlab = xlab, ylab = ylab, main = main,
    cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
  )
  # The identity first, so the points and the fitted line sit on top of it.
  graphics::abline(a = 0, b = 1, col = theme$guide, lwd = 2, lty = 3)
  if (draw_points) {
    graphics::points(observed, predicted, pch = 16, col = col)
  }
  # Drawn from the two numbers the table holds rather than from a fresh fit, so
  # that the line and `$metrics` are the same line.
  if (!is.na(row$calib_slope)) {
    graphics::abline(a = row$calib_intercept, b = row$calib_slope,
                     col = col, lwd = lwd)
  }
  invisible(NULL)
}


#' The equation of a calibration line, written the way it reads
#'
#' @keywords internal
#' @noRd
sa_calibration_label <- function(row) {
  if (is.na(row$calib_slope)) {
    return("no calibration line")
  }
  intercept <- row$calib_intercept
  paste0("y = ", sa_fmt_num(row$calib_slope, 3),
         "x ", if (intercept < 0) "- " else "+ ",
         sa_fmt_num(abs(intercept), 3))
}


#' Draw predicted against observed for an evaluated regression
#'
#' The picture of an [evaluate_regression_models()] result: each model's
#' predictions against the outcome they were predicting, with the identity line
#' to read them against and the calibration line to read the identity against.
#'
#' @details
#' Two lines are drawn in every panel. The dotted grey one is the identity,
#' where a prediction that was exactly right would lie. The solid coloured one
#' is `lm(predicted ~ observed)`, taken from the `calib_slope` and
#' `calib_intercept` of `$metrics` rather than fitted again here, so the picture
#' and the table cannot drift apart. A slope under one is a model whose
#' predictions are squeezed towards their own mean, which is the usual shape of
#' a fit scored on rows it has not seen.
#'
#' `type` decides how more than one model is shown. `"panel"` gives each its own
#' square, which is the default past one model because two clouds of points on
#' shared axes make a third cloud that belongs to neither. `"overlay"` puts them
#' all in one square, which is readable with `points = FALSE`, where what is
#' compared is the calibration lines. `"auto"` is `"overlay"` for a single model
#' and `"panel"` beyond that.
#'
#' Both axes span the same range in every panel, taken over every model drawn,
#' so the panels are comparable and the identity line is the diagonal of the
#' square rather than an arbitrary chord.
#'
#' @param performance_result A regression evaluation, as returned by
#'   [evaluate_regression_models()].
#' @param models Which models to draw and in what order, or `NULL` for all of
#'   them in the order the evaluation holds, which puts the baseline first.
#' @param type `"auto"`, `"overlay"` or `"panel"`, as described above.
#' @param panel_nrow Rows of panels under `type = "panel"`. `NULL` draws them in
#'   one row.
#' @param points Whether to draw the predictions themselves. `FALSE` leaves the
#'   calibration lines, which is what makes a crowded overlay readable.
#' @param anno_corr Whether to report each model's correlation beside its name.
#' @param anno_rsq Whether to report each model's held-out R-squared from
#'   `$metrics$r_squared`. This is `1 - SSE/SST` on the scored rows, not
#'   `cor` squared, and the two are annotated separately for the same reason
#'   the table carries both.
#' @param anno_lm Whether to report each model's calibration line as an
#'   equation.
#' @param dark Whether to draw on a dark background.
#' @param lim Range of both axes, or `NULL` to span every value drawn.
#' @param col One colour, or one per drawn model. `NULL` takes them from
#'   `hcl.colors(n, "Dark 2")`.
#' @param lwd Width of the calibration lines.
#' @param xlab,ylab,main Axis and figure labels. `NULL` builds them from the
#'   evaluation.
#' @param cex.axis,cex.lab,cex.main,cex.legend,cex.anno Relative text sizes.
#'   `cex.legend` sizes the model-name legend; `cex.anno` sizes the correlation,
#'   R-squared and calibration labels drawn inside the plot area. `NULL` matches
#'   `cex.legend`.
#'
#' @return The rows of `$metrics` that were drawn, invisibly, carrying the
#'   resolved `type` as a `"view"` attribute.
#'
#' @seealso [evaluate_regression_models()] for the result this draws, and
#'   [draw_roc_curve()] for the classification counterpart.
#'
#' @examples
#' train <- mtcars[1:24, ]
#' test <- mtcars[25:32, ]
#' full <- fit_linear_regression(train, outcome = "mpg",
#'                               predictors = c("wt", "hp", "disp"),
#'                               cv = FALSE)
#' small <- fit_linear_regression(train, outcome = "mpg",
#'                                predictors = "wt", cv = FALSE)
#' res <- evaluate_regression_models(full, list(weight_only = small),
#'                                   newdata = test)
#'
#' draw_prediction_plot(res, anno_lm = TRUE)
#' draw_prediction_plot(res, type = "overlay", points = FALSE,
#'                      anno_corr = TRUE, anno_rsq = TRUE)
#'
#' @export
draw_prediction_plot <- function(performance_result,
                                 models = NULL,
                                 type = c("auto", "overlay", "panel"),
                                 panel_nrow = NULL,
                                 points = TRUE,
                                 anno_corr = FALSE,
                                 anno_rsq = FALSE,
                                 anno_lm = FALSE,
                                 dark = FALSE,
                                 lim = NULL,
                                 col = NULL,
                                 lwd = 2,
                                 xlab = NULL,
                                 ylab = NULL,
                                 main = NULL,
                                 cex.axis = 1.2,
                                 cex.lab = 1.3,
                                 cex.main = 1.3,
                                 cex.legend = 1.1,
                                 cex.anno = NULL) {
  sa_performance_input(performance_result, "regression_performance",
                       "performance_result", "draw_roc_curve()")
  type <- match.arg(type)
  sa_check_flag(points, "points")
  sa_check_flag(anno_corr, "anno_corr")
  sa_check_flag(anno_rsq, "anno_rsq")
  sa_check_flag(anno_lm, "anno_lm")
  sa_check_flag(dark, "dark")
  sa_check_lim(lim, "lim")
  if (is.null(cex.anno)) {
    cex.anno <- cex.legend
  } else {
    sa_check_scalar_num(cex.anno, "cex.anno", 0, lower_open = TRUE)
  }
  if (!is.null(panel_nrow)) {
    panel_nrow <- sa_check_count(panel_nrow, "panel_nrow", 1)
  }

  drawn_models <- sa_performance_models(performance_result, models)
  n_model <- length(drawn_models)
  if (identical(type, "auto")) {
    type <- if (n_model > 1L) "panel" else "overlay"
  }
  cols <- sa_performance_colours(n_model, col)
  theme <- sa_plot_theme(dark)

  metrics <- performance_result$metrics
  metrics <- metrics[match(drawn_models, metrics$model), , drop = FALSE]
  rownames(metrics) <- NULL
  predictions <- performance_result$predictions
  by_model <- lapply(drawn_models, function(nm) {
    predictions[predictions$model == nm, , drop = FALSE]
  })

  # One range over every model drawn, so the panels are comparable and the
  # identity is the diagonal of the square rather than a chord across it.
  span <- if (is.null(lim)) {
    range(c(predictions$observed[predictions$model %in% drawn_models],
            predictions$predicted[predictions$model %in% drawn_models]),
          finite = TRUE)
  } else {
    lim
  }

  label_x <- if (is.null(xlab)) {
    paste0("Observed ", performance_result$design$outcome)
  } else {
    xlab
  }
  label_y <- if (is.null(ylab)) "Predicted" else ylab

  annotation <- function(i) {
    c(if (anno_corr) paste0("Corr = ", sa_fmt_num(metrics$cor[i], 3)),
      if (anno_rsq) paste0("R-sq = ", sa_fmt_num(metrics$r_squared[i], 3)),
      if (anno_lm) sa_calibration_label(metrics[i, ]))
  }

  # Only the parameters this function sets are put back, not a blanket
  # par(no.readonly = TRUE) snapshot. That snapshot also carries `fin`, `pin`
  # and `mai`, which are absolute sizes: restoring them pins the figure to the
  # size this plot happened to be drawn at, so the next plot on a device that
  # has since been resized is redrawn small in a corner of it. `mfrow` comes
  # back with the rest, since `layout()` overwrites whatever grid the caller
  # had set up, and `oma` because a title over several panels lives there.
  old_par <- graphics::par(c("bg", "fg", "col.axis", "col.lab", "col.main",
                             "mar", "mfrow", "oma"))
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)
  graphics::par(bg = theme$bg, fg = theme$fg, col.axis = theme$fg,
                col.lab = theme$fg, col.main = theme$fg)

  if (identical(type, "overlay")) {
    # The legend gets a panel of its own rather than a corner of the square,
    # where it would sit over the region a well calibrated model draws through.
    # It is sized to the model names, which is all it carries: an annotation is
    # a correlation, an R-squared and an equation, and putting those out here too
    # would hand half the figure to the legend, so they go inside beside the
    # lines they describe instead.
    wanted <- (max(nchar(drawn_models)) + 5) * cex.legend *
      graphics::par("cin")[1] + 0.2
    share <- min(0.35, wanted / graphics::par("din")[1])
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1 - share, share))
    graphics::par(mar = c(5.1, 4.6, 4.1, 1.1))
    graphics::plot.default(
      span, span, type = "n", bty = "n", xlim = span, ylim = span,
      xlab = label_x, ylab = label_y,
      main = if (is.null(main)) "Predicted against observed" else main,
      cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
    )
    graphics::abline(a = 0, b = 1, col = theme$guide, lwd = 2, lty = 3)
    for (i in seq_len(n_model)) {
      if (points) {
        graphics::points(by_model[[i]]$observed, by_model[[i]]$predicted,
                         pch = 16, col = cols[i])
      }
      if (!is.na(metrics$calib_slope[i])) {
        graphics::abline(a = metrics$calib_intercept[i],
                         b = metrics$calib_slope[i], col = cols[i], lwd = lwd)
      }
    }

    # One line per model, in that model's colour and in the order the legend
    # beside it lists them, so the two read together without the names being
    # written twice.
    notes <- vapply(seq_len(n_model), function(i) {
      paste(annotation(i), collapse = ", ")
    }, character(1))
    if (any(nzchar(notes))) {
      graphics::legend("bottomright", legend = notes, bty = "n",
                       cex = cex.anno, text.col = cols)
    }

    graphics::par(mar = c(5, 0, 4, 1))
    graphics::plot.new()
    graphics::legend("center", legend = drawn_models, col = cols, lwd = lwd,
                     bty = "n", cex = cex.legend, text.col = theme$fg)
  } else {
    n_row <- min(if (is.null(panel_nrow)) 1L else panel_nrow, n_model)
    n_col <- ceiling(n_model / n_row)
    slots <- c(seq_len(n_model), rep(0L, n_row * n_col - n_model))
    # A title belongs to the figure rather than to a panel of it, so it moves
    # to the outer margin and the panels keep their model names.
    if (!is.null(main)) {
      graphics::par(oma = c(0, 0, 3, 0))
    }
    graphics::layout(matrix(slots, nrow = n_row, byrow = TRUE))
    graphics::par(mar = c(5.1, 4.6, 3.6, 1.6))

    for (i in seq_len(n_model)) {
      sa_prediction_panel(
        by_model[[i]]$observed, by_model[[i]]$predicted, metrics[i, ],
        lim = span, theme = theme, col = cols[i], lwd = lwd,
        draw_points = points, xlab = label_x, ylab = label_y,
        main = drawn_models[i], cex.axis = cex.axis, cex.lab = cex.lab,
        cex.main = cex.main
      )
      note <- annotation(i)
      if (length(note) > 0L) {
        graphics::legend("bottomright", legend = note, bty = "n",
                         cex = cex.anno, text.col = theme$fg)
      }
    }
    if (!is.null(main)) {
      graphics::mtext(main, side = 3, outer = TRUE, line = 0.8,
                      cex = cex.main, font = 2)
    }
  }

  # The view is carried on the result because `type = "auto"` resolves it here
  # and the caller would otherwise have no way to find out which of the two it
  # got, short of counting the panels on the device.
  attr(metrics, "view") <- type
  invisible(metrics)
}
