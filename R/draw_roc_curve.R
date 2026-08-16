# Overlaid ROC curves, the picture of a classification evaluation. Unlike the
# predicted-against-observed scatter, this one overlays rather than panels: a
# curve is a line rather than a cloud, several of them share the unit square
# without obscuring each other, and the whole question the plot is asked is
# which of them is above the others.
#
# The legend goes inside the panel at the bottom right, which is the one corner
# of a ROC plot that no useful curve passes through.


#' Draw the ROC curves of an evaluated classification
#'
#' The picture of an [evaluate_classification_models()] result: one curve per
#' model, all on the rows every model was scored on, with the chance diagonal to
#' read them against.
#'
#' @details
#' The curves are `$curves`, the operating points the evaluation already
#' computed, rather than anything recomputed here, so the picture and the `auc`
#' column of `$metrics` describe the same curve. Consecutive points are joined
#' by straight lines, which is what makes the area under the drawn curve the
#' `auc` beside it: a run of tied predictions cannot be separated by any
#' threshold and the curve crosses it diagonally.
#'
#' Every model is drawn against the same rows, since that is what the evaluation
#' scored them on, so two curves crossing is a statement about the models rather
#' than about which rows each of them managed.
#'
#' @param performance_result A classification evaluation, as returned by
#'   [evaluate_classification_models()].
#' @param models Which models to draw and in what order, or `NULL` for all of
#'   them in the order the evaluation holds, which puts the baseline first.
#' @param anno_auc Whether to add each model's AUC to its legend entry.
#' @param chance Whether to draw the chance diagonal.
#' @param dark Whether to draw on a dark background.
#' @param col One colour, or one per drawn model. `NULL` takes them from
#'   `hcl.colors(n, "Dark 2")`.
#' @param lwd Width of the curves.
#' @param lty One line type, or one per drawn model.
#' @param legend_pos Where to put the legend, or `NULL` for no legend.
#' @param xlab,ylab,main Axis and figure labels. `NULL` builds them from the
#'   evaluation.
#' @param cex.axis,cex.lab,cex.main,cex.legend,cex.anno Relative text sizes.
#'   `cex.legend` sizes the model-name legend when `anno_auc = FALSE`.
#'   `cex.anno` sizes it once each entry carries an AUC; `NULL` matches
#'   `cex.legend`.
#'
#' @return The rows of `$metrics` that were drawn, invisibly.
#'
#' @seealso [evaluate_classification_models()] for the result this draws, and
#'   [draw_prediction_plot()] for the regression counterpart.
#'
#' @examples
#' iris2 <- iris[iris$Species != "setosa", ]
#' iris2$Species <- factor(iris2$Species)
#' train <- iris2[c(1:35, 51:85), ]
#' test <- iris2[c(36:50, 86:100), ]
#'
#' full <- fit_logistic_regression(train, outcome = "Species", cv = FALSE)
#' petal <- fit_logistic_regression(train, outcome = "Species",
#'                                  predictors = "Petal.Width", cv = FALSE)
#' res <- evaluate_classification_models(full, list(petal_only = petal),
#'                                       newdata = test)
#'
#' draw_roc_curve(res, anno_auc = TRUE)
#'
#' @export
draw_roc_curve <- function(performance_result,
                           models = NULL,
                           anno_auc = FALSE,
                           chance = TRUE,
                           dark = FALSE,
                           col = NULL,
                           lwd = 2,
                           lty = 1,
                           legend_pos = "bottomright",
                           xlab = NULL,
                           ylab = NULL,
                           main = NULL,
                           cex.axis = 1.2,
                           cex.lab = 1.3,
                           cex.main = 1.3,
                           cex.legend = 1.1,
                           cex.anno = NULL) {
  sa_performance_input(performance_result, "classification_performance",
                       "performance_result", "draw_prediction_plot()")
  sa_check_flag(anno_auc, "anno_auc")
  sa_check_flag(chance, "chance")
  sa_check_flag(dark, "dark")
  if (is.null(cex.anno)) {
    cex.anno <- cex.legend
  } else {
    sa_check_scalar_num(cex.anno, "cex.anno", 0, lower_open = TRUE)
  }

  drawn_models <- sa_performance_models(performance_result, models)
  n_model <- length(drawn_models)
  cols <- sa_performance_colours(n_model, col)
  if (length(lty) != 1L && length(lty) != n_model) {
    stop("`lty` must hold one line type, or one per drawn model (", n_model,
         ").", call. = FALSE)
  }
  ltys <- rep(lty, length.out = n_model)
  theme <- sa_plot_theme(dark)

  metrics <- performance_result$metrics
  metrics <- metrics[match(drawn_models, metrics$model), , drop = FALSE]
  rownames(metrics) <- NULL

  # Only the parameters this function sets are put back, not a blanket
  # par(no.readonly = TRUE) snapshot. That snapshot also carries `fin`, `pin`
  # and `mai`, which are absolute sizes: restoring them pins the figure to the
  # size this plot happened to be drawn at, so the next plot on a device that
  # has since been resized is redrawn small in a corner of it.
  old_par <- graphics::par(c("bg", "fg", "col.axis", "col.lab", "col.main",
                             "mar"))
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(bg = theme$bg, fg = theme$fg, col.axis = theme$fg,
                col.lab = theme$fg, col.main = theme$fg,
                mar = c(5.1, 4.6, 4.1, 2.1))

  default_main <- paste0("ROC: ", performance_result$design$outcome_lv[2],
                         " against ", performance_result$design$outcome_lv[1])
  graphics::plot.default(
    c(0, 1), c(0, 1), type = "n", bty = "n", xlim = c(0, 1), ylim = c(0, 1),
    xlab = if (is.null(xlab)) "1 - specificity" else xlab,
    ylab = if (is.null(ylab)) "sensitivity" else ylab,
    main = if (is.null(main)) default_main else main,
    cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
  )
  if (chance) {
    graphics::abline(a = 0, b = 1, col = theme$guide, lwd = 2, lty = 3)
  }

  curves <- performance_result$curves
  for (i in seq_len(n_model)) {
    points <- curves[curves$model == drawn_models[i], , drop = FALSE]
    points <- points[order(1 - points$specificity, points$sensitivity), ]
    graphics::lines(1 - points$specificity, points$sensitivity,
                    col = cols[i], lwd = lwd, lty = ltys[i])
  }

  if (!is.null(legend_pos)) {
    entries <- if (anno_auc) {
      paste0(drawn_models, "  (", sa_fmt_num(metrics$auc, 3), ")")
    } else {
      drawn_models
    }
    graphics::legend(legend_pos, legend = entries, col = cols, lwd = lwd,
                     lty = ltys, bty = "n",
                     cex = if (anno_auc) cex.anno else cex.legend,
                     text.col = theme$fg)
  }

  invisible(metrics)
}
