#' Draw an interaction plot of a factorial comparison
#'
#' Joins the cell means of one factor across the levels of another, one line per
#' level of the tracing factor, which is the picture an interaction term is a
#' test of: parallel lines are additive factors and crossing or fanning lines are
#' the interaction the term reports a p-value for.
#'
#' The means come from `comparison_result$cells`, so the plot and the F tests
#' beside it describe the same fit rather than two passes over the data. Where
#' the design holds more factors than the two being drawn, the remaining ones are
#' averaged away, unweighted, which is the same marginal mean the post-hoc stage
#' contrasts: the gap between two points of one line is the `estimate` its
#' contrast reports in `comparison_result$posthoc`.
#'
#' Three views are available:
#'
#' \describe{
#'   \item{`"pairwise"`}{One pair of factors, `trace` against `x`, with every
#'     other factor averaged away. One panel per feature, so several features can
#'     be read side by side.}
#'   \item{`"matrix"`}{Every pair of factors at once, in an upper triangle: the
#'     row is the tracing factor and the column the one on the x axis. One
#'     feature only, since the panels are spent on the factors.}
#'   \item{`"facet"`}{`trace` against `x` again, but with the levels of `facet`
#'     kept apart in panels of their own rather than averaged away, which is what
#'     shows a two-factor interaction that itself depends on a third factor. One
#'     feature only, for the same reason.}
#' }
#'
#' `type = "auto"`, the default, reads the arguments: naming `facet` asks for the
#' facet view, naming `x` or `trace` for the pairwise one, and naming none of
#' them gives the pairwise view for two factors and the matrix for three or more.
#'
#' @param comparison_result A factorial comparison result, as returned by
#'   [compare_factorial_groups()].
#' @param x Name of the factor on the x axis.
#' @param trace Name of the factor one line is drawn per. Unnamed, the factors
#'   are taken in declaration order, `trace` first, which is the order
#'   [draw_grouped_boxplot()] colours a crossed design in.
#' @param facet Name of the factor whose levels are kept in panels of their own
#'   instead of averaged away. Required by `type = "facet"` and refused by the
#'   other two views.
#' @param type `"auto"`, `"pairwise"`, `"matrix"` or `"facet"`.
#' @param feats Character vector of features to draw, in panel order. `NULL`,
#'   the default, draws every feature in the pairwise view and the first feature
#'   in the other two, which is reported in a `message()`.
#' @param errorbar `"none"`, `"se"` for one pooled standard error either side of
#'   each mean, or `"ci"` for a confidence interval at the `conf_level` the
#'   comparison was run at.
#' @param panel_nrow Number of rows to arrange the panels in, or `NULL` for a
#'   default that depends on the view. Not used by the matrix view, whose grid
#'   the factors fix.
#' @param dark If `TRUE`, use a dark background with light text.
#' @param ylim Numeric length-2 y axis range, or `NULL` to derive it from the
#'   means and their bars.
#' @param xlab,ylab,main Axis and title labels. All three are derived from the
#'   result when left `NULL`. A single `xlab` covers every panel of the matrix
#'   view, whose panels are each on a different factor, so leaving it `NULL`
#'   there is what lets each panel name its own.
#' @param col Colours for the levels of the tracing factor, recycled if short.
#'   `NULL` takes them from `grDevices::hcl.colors()`, the palette
#'   [draw_grouped_boxplot()] draws a crossed design in.
#' @param lwd Line width of the traces.
#' @param cex.axis,cex.lab,cex.main,cex.legend Character expansion for the axis
#'   annotation, the axis label, the title and the legend.
#'
#' @return The plotted means, invisibly, one row per point in the order they were
#'   drawn, with columns `panel`, `features`, `x_factor`, `x_level`,
#'   `trace_factor`, `trace_level`, `n_cells`, `mean`, `se`, `lower_conf` and
#'   `upper_conf`. `n_cells` is how many cells were averaged into the point, so
#'   `1` marks a cell mean and anything more a marginal one. The view that
#'   `type = "auto"` resolved to is attached as the attribute `"view"`.
#'
#' @details
#' The function changes graphical parameters and the panel layout, and restores
#' both on exit, so the caller's device is left as it was found.
#'
#' Each level of the tracing factor gets a colour, a plotting symbol and a line
#' type of its own, so the lines stay apart in greyscale and for a reader who
#' cannot tell the two colours apart.
#'
#' @section What an error bar here is not:
#' The bars are built from `cells$se`, which is `sqrt(ms_error / n)` pooled over
#' the whole model, so a bar on a marginal mean over a set of cells has half-width
#' `sqrt(sum(se^2) / length(se)^2)` times one or the `qt()` multiplier. That is
#' the standard error of the mean being drawn, not of the difference between two
#' of them, and the difference is what the interaction is about. Two bars that
#' overlap do not settle the term test, which is in `$terms`, and two points whose
#' gap you want tested are in `$posthoc`.
#'
#' @seealso [compare_factorial_groups()] for the result this reads,
#'   [draw_grouped_boxplot()] for the observations behind the means, and
#'   [draw_volcano_plot()] with `estimate_significance(by = "term")` for which
#'   features have an interaction at all.
#'
#' @examples
#' sim <- simulate_factorial_groups(n_feats = 4, n_per_cell = 8, seed = 1)
#' res <- do.call(compare_factorial_groups, c(sim$args, list(diagnose = FALSE)))
#'
#' ## One panel per feature, the first factor tracing the second
#' draw_interaction_plot(res)
#'
#' ## The same two factors the other way round, with standard error bars
#' draw_interaction_plot(res, x = "treatment", trace = "sex", errorbar = "se")
#'
#' ## Two features only, stacked
#' draw_interaction_plot(res, feats = res$features[1:2], panel_nrow = 2)
#'
#' ## Every pair of factors of a three-factor design, for one feature
#' sim3 <- simulate_factorial_groups(
#'   n_feats = 2, n_per_cell = 8, seed = 2,
#'   factor_lv = list(treatment = c("control", "treat_A"),
#'                    sex       = c("male", "female"),
#'                    site      = c("north", "south", "east"))
#' )
#' res3 <- do.call(compare_factorial_groups,
#'                 c(sim3$args, list(diagnose = FALSE)))
#' draw_interaction_plot(res3, type = "matrix")
#'
#' ## The third factor kept apart in panels instead of averaged away
#' draw_interaction_plot(res3, x = "site", trace = "treatment", facet = "sex",
#'                       errorbar = "ci")
#'
#' @export
draw_interaction_plot <- function(comparison_result,
                                  x = NULL,
                                  trace = NULL,
                                  facet = NULL,
                                  type = c("auto", "pairwise", "matrix",
                                           "facet"),
                                  feats = NULL,
                                  errorbar = c("none", "se", "ci"),
                                  panel_nrow = NULL,
                                  dark = FALSE,
                                  ylim = NULL,
                                  xlab = NULL,
                                  ylab = NULL,
                                  main = NULL,
                                  col = NULL,
                                  lwd = 2,
                                  cex.axis = 1.2,
                                  cex.lab = 1.3,
                                  cex.main = 1.3,
                                  cex.legend = 1.1) {

  type <- match.arg(type)
  errorbar <- match.arg(errorbar)
  sa_check_flag(dark, "dark")
  sa_check_scalar_num(lwd, "lwd", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.legend, "cex.legend", 0, lower_open = TRUE)
  sa_check_lim(ylim, "ylim")
  if (!is.null(panel_nrow)) {
    sa_check_count(panel_nrow, "panel_nrow", 1)
  }

  cells <- sa_inter_cells(comparison_result)
  factor_lv <- comparison_result$design$factor_lv
  roles <- sa_inter_roles(names(factor_lv), x, trace, facet, type)
  type <- roles$type
  drawn_feats <- sa_inter_feats(comparison_result, feats, type)

  # One multiplier per feature, because a confidence interval is read off the
  # error degrees of freedom of that feature's own fit.
  mult <- sa_inter_multiplier(comparison_result, drawn_feats, errorbar)

  panels <- sa_inter_panels(cells, factor_lv, roles, drawn_feats, mult, type)
  drawn <- do.call(rbind, lapply(panels, function(p) p$tbl))
  rownames(drawn) <- NULL

  if (!any(is.finite(drawn$mean))) {
    stop("none of the cells of ", paste(drawn_feats, collapse = ", "),
         " holds a mean to draw, so the model was not fitted for any of them. ",
         "`comparison_result$tests$anova_test` says why.", call. = FALSE)
  }

  if (is.null(ylab)) {
    ylab <- sa_inter_ylab(drawn, comparison_result$parameters$input_scale)
  }
  sa_inter_draw(panels, drawn, factor_lv, roles, type, errorbar, panel_nrow,
                dark, ylim, xlab, ylab, main, col, lwd, cex.axis, cex.lab,
                cex.main, cex.legend)

  # Carried on the result because `type = "auto"` resolves here and the caller
  # would otherwise have no way of finding out which of the three it got.
  attr(drawn, "view") <- type
  invisible(drawn)
}


#' The cell table of a factorial comparison, or the reason there is none
#'
#' @keywords internal
#' @noRd
sa_inter_cells <- function(res) {
  if (!inherits(res, "sa_factorial")) {
    stop("`comparison_result` must be a factorial comparison result, as ",
         "returned by compare_factorial_groups(). An interaction is a ",
         "statement about two crossed factors, and a result over a single ",
         "factor holds no second one to trace against.", call. = FALSE)
  }
  if (is.null(res$cells)) {
    stop("`comparison_result` carries no `$cells` table, so it was produced ",
         "by a version of the package that did not record the cell means. ",
         "Re-run compare_factorial_groups() on the same data.", call. = FALSE)
  }
  res$cells
}


#' Settle which view is being drawn and which factor plays which part
#'
#' @param fac Factor names in declaration order.
#' @param x,trace,facet The arguments as received, any of them `NULL`.
#' @param type The argument as received, already through `match.arg()`.
#'
#' @return List with `type`, `x`, `trace` and `facet`, the last `NULL` outside
#'   the facet view.
#'
#' @keywords internal
#' @noRd
sa_inter_roles <- function(fac, x, trace, facet, type) {
  check_one <- function(nm, arg) {
    if (is.null(nm)) {
      return(invisible(NULL))
    }
    if (!is.character(nm) || length(nm) != 1L || is.na(nm)) {
      stop("`", arg, "` must be a single factor name, or NULL.", call. = FALSE)
    }
    if (!nm %in% fac) {
      stop("`", arg, "` must name one of the factors of the comparison: ",
           paste(fac, collapse = ", "), ". Got ", nm, ".", call. = FALSE)
    }
    invisible(NULL)
  }
  check_one(trace, "trace")
  check_one(x, "x")
  check_one(facet, "facet")

  named <- c(trace = trace, x = x, facet = facet)
  if (length(named) == 0L) {
    named <- character(0)
  }
  if (anyDuplicated(named) > 0L) {
    stop("`x`, `trace` and `facet` must name different factors; ",
         paste(unique(named[duplicated(named)]), collapse = ", "),
         " is named twice. One factor cannot be two axes of the same plot.",
         call. = FALSE)
  }

  if (type == "auto") {
    type <- if (!is.null(facet)) {
      "facet"
    } else if (length(named) > 0L) {
      "pairwise"
    } else if (length(fac) > 2L) {
      # Three factors have three pairs to look at and no reason to prefer one,
      # so the default is to show them all rather than to pick silently.
      "matrix"
    } else {
      "pairwise"
    }
  }

  if (type == "matrix") {
    if (length(named) > 0L) {
      stop("`type = \"matrix\"` draws every pair of factors, so there is ",
           "nothing for ", paste(names(named), collapse = ", "),
           " to choose. Drop it, or use type = \"pairwise\" to draw one pair.",
           call. = FALSE)
    }
    return(list(type = type, x = NULL, trace = NULL, facet = NULL))
  }

  if (type == "facet") {
    if (is.null(facet)) {
      stop("`type = \"facet\"` needs `facet` to name the factor whose levels ",
           "go in panels of their own. Without it the third factor is ",
           "averaged away, which is type = \"pairwise\".", call. = FALSE)
    }
    if (length(fac) < 3L) {
      stop("`type = \"facet\"` needs at least three factors, one for each of ",
           "`x`, `trace` and `facet`, and this comparison holds ", length(fac),
           ": ", paste(fac, collapse = ", "), ".", call. = FALSE)
    }
  } else if (!is.null(facet)) {
    stop("`facet` belongs to type = \"facet\". In the pairwise view every ",
         "factor other than `x` and `trace` is averaged away, so there is no ",
         "panel for ", facet, " to be kept in.", call. = FALSE)
  }

  # Whatever the caller did not name is filled from declaration order, tracing
  # factor first, which is the factor a crossed boxplot colours by.
  free <- setdiff(fac, named)
  for (part in setdiff(c("trace", "x"), names(named))) {
    named[[part]] <- free[1L]
    free <- free[-1L]
  }

  list(type = type, x = named[["x"]], trace = named[["trace"]],
       facet = if (type == "facet") named[["facet"]] else NULL)
}


#' Which features the view has room for
#'
#' @keywords internal
#' @noRd
sa_inter_feats <- function(res, feats, type) {
  one_panel_each <- type == "pairwise"

  if (is.null(feats)) {
    if (one_panel_each) {
      return(res$features)
    }
    # The panels of these two views are spent on the factors, so a second
    # feature has nowhere to go.
    message("Drawing ", res$features[1], ". The ", type, " view spends its ",
            "panels on the factors, so it draws one feature at a time; name ",
            "another in `feats`.")
    return(res$features[1L])
  }

  sa_check_feat_names(feats)
  unknown <- setdiff(feats, res$features)
  if (length(unknown) > 0L) {
    stop("`feats` must name features present in the comparison: ",
         paste(res$features, collapse = ", "), ". Not found: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (!one_panel_each && length(feats) > 1L) {
    stop("`type = \"", type, "\"` draws one feature at a time, and `feats` ",
         "names ", length(feats), ": ", paste(feats, collapse = ", "),
         ". Its panels are spent on the factors; use type = \"pairwise\" for ",
         "a panel per feature.", call. = FALSE)
  }
  feats
}


#' How many standard errors wide a bar is, per feature
#'
#' @keywords internal
#' @noRd
sa_inter_multiplier <- function(res, feats, errorbar) {
  if (errorbar != "ci") {
    # A bar of one standard error needs no distribution behind it, and `NA` for
    # no bar at all leaves the two bound columns absent rather than equal to the
    # mean, which would draw as a bar of zero height.
    return(stats::setNames(rep(if (errorbar == "se") 1 else NA_real_,
                               length(feats)),
                           feats))
  }
  tbl <- res$tests$anova_test
  df <- tbl$df2[match(feats, tbl$features)]
  conf <- res$parameters$conf_level
  stats::setNames(stats::qt(1 - (1 - conf) / 2, df), feats)
}


#' The marginal means of one pair of factors
#'
#' Averaging the cell means without weighting them is what makes the gap between
#' two points of a line the `estimate` the post-hoc stage reports for that
#' contrast: `sa_factorial_tukey()` defines a marginal mean the same way. The
#' pooled `se` of the cells combines the same way its variance does, which is
#' why the cell table stores it pooled rather than within the cell.
#'
#' @param cells The cell rows of one feature, possibly already cut down to one
#'   level of a facet.
#' @param factor_lv Named list of factors and their levels, fixing the order the
#'   points are drawn in.
#' @param x,trace Factor names.
#'
#' @return data.frame of one row per point, x varying fastest so that the rows of
#'   one trace are consecutive.
#'
#' @keywords internal
#' @noRd
sa_inter_marginal <- function(cells, factor_lv, x, trace) {
  out <- expand.grid(x_level = factor_lv[[x]], trace_level = factor_lv[[trace]],
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  out$n_cells <- 0L
  out$mean <- NA_real_
  out$se <- NA_real_

  for (i in seq_len(nrow(out))) {
    at <- cells[[x]] == out$x_level[i] & cells[[trace]] == out$trace_level[i]
    k <- sum(at)
    out$n_cells[i] <- k
    if (k == 0L) next
    # No `na.rm`: a marginal mean over a cell that has none is not a mean of the
    # rest, it is unknown, and a gap in the line says so.
    out$mean[i] <- mean(cells$mean[at])
    out$se[i] <- sqrt(sum(cells$se[at]^2) / k^2)
  }
  out
}


#' Assemble one table per panel, in the order the panels are drawn
#'
#' @return List of panels, each with `label`, `x`, `trace`, `tbl` and, for the
#'   matrix view, the `row` and `col` of the triangle it sits in.
#'
#' @keywords internal
#' @noRd
sa_inter_panels <- function(cells, factor_lv, roles, feats, mult, type) {
  fac <- names(factor_lv)

  one <- function(label, feature, x, trace, keep = NULL) {
    rows <- cells[cells$features == feature, , drop = FALSE]
    if (!is.null(keep)) {
      rows <- rows[rows[[names(keep)]] == keep[[1L]], , drop = FALSE]
    }
    tbl <- sa_inter_marginal(rows, factor_lv, x, trace)
    half <- mult[[feature]] * tbl$se
    out <- data.frame(
      panel        = label,
      features     = feature,
      x_factor     = x,
      x_level      = tbl$x_level,
      trace_factor = trace,
      trace_level  = tbl$trace_level,
      n_cells      = tbl$n_cells,
      mean         = tbl$mean,
      se           = tbl$se,
      lower_conf   = tbl$mean - half,
      upper_conf   = tbl$mean + half,
      stringsAsFactors = FALSE
    )
    list(label = label, x = x, trace = trace, tbl = out)
  }

  if (type == "pairwise") {
    return(lapply(feats, function(f) {
      one(f, f, roles$x, roles$trace)
    }))
  }

  if (type == "facet") {
    return(lapply(factor_lv[[roles$facet]], function(lv) {
      one(paste0(roles$facet, ": ", lv), feats[1L], roles$x, roles$trace,
          keep = stats::setNames(list(lv), roles$facet))
    }))
  }

  # The matrix is an upper triangle of one panel per pair: row `i` traces factor
  # `i` and column `j` puts factor `j + 1` on the x axis, so a pair is drawn
  # once, with the earlier factor tracing the later one.
  out <- list()
  for (i in seq_len(length(fac) - 1L)) {
    for (j in i:(length(fac) - 1L)) {
      panel <- one(paste0(fac[i], " x ", fac[j + 1L]), feats[1L],
                   fac[j + 1L], fac[i])
      panel$row <- i
      panel$col <- j
      out[[length(out) + 1L]] <- panel
    }
  }
  out
}


#' Lay the panels out and draw them
#'
#' @keywords internal
#' @noRd
sa_inter_draw <- function(panels, drawn, factor_lv, roles, type, errorbar,
                          panel_nrow, dark, ylim, xlab, ylab, main, col, lwd,
                          cex.axis, cex.lab, cex.main, cex.legend) {
  theme <- sa_plot_theme(dark)
  n_panel <- length(panels)
  # The panels of the matrix view each trace a different factor, so each carries
  # its own key; the other two views trace one factor throughout and share a
  # legend panel on the right, where it cannot cover the lines it describes.
  shared_legend <- type != "matrix"
  # Panels of one feature are panels of one quantity and have to sit on one
  # scale; panels of different features are different quantities, and a shared
  # range would flatten every one of them.
  free_scale <- type == "pairwise" && n_panel > 1L

  span <- if (is.null(ylim)) {
    sa_inter_span(drawn, errorbar)
  } else {
    NULL
  }

  if (type == "matrix") {
    if (!is.null(panel_nrow)) {
      message("`panel_nrow` is not used by the matrix view, whose grid is ",
              "fixed by the factors.")
    }
    side <- length(factor_lv) - 1L
    grid <- matrix(0L, nrow = side, ncol = side)
    for (p in seq_len(n_panel)) {
      grid[panels[[p]]$row, panels[[p]]$col] <- p
    }
    n_row <- side
    n_col <- side
  } else {
    if (is.null(panel_nrow)) {
      # A panel per feature is many panels, and a single row of them is a strip
      # too thin to read, so those go as near square as they get. Panels sharing
      # a scale read against each other best side by side.
      panel_nrow <- if (free_scale) max(1L, round(sqrt(n_panel))) else 1L
    }
    n_row <- min(panel_nrow, n_panel)
    n_col <- ceiling(n_panel / n_row)
    grid <- matrix(c(seq_len(n_panel), rep(0L, n_row * n_col - n_panel)),
                   nrow = n_row, byrow = TRUE)
  }

  strip <- n_panel > 1L
  default_main <- sa_inter_main(panels, roles, type, strip)
  figure_main <- if (is.null(main)) default_main else main
  outer_main <- strip && !is.null(figure_main)

  # Only the parameters this function sets are put back, not a blanket
  # par(no.readonly = TRUE) snapshot, whose `fin`, `pin` and `mai` are absolute
  # sizes: restoring those pins the figure to the size this plot happened to be
  # drawn at. `mfrow` comes back with the rest because `layout()` overwrites
  # whatever grid the caller had, and `oma` because a title over several panels
  # lives there.
  old_par <- graphics::par(c("bg", "fg", "col.axis", "col.lab", "col.main",
                             "mar", "mfrow", "oma"))
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  if (outer_main) {
    graphics::par(oma = c(0, 0, 3, 0))
  }
  if (shared_legend) {
    graphics::layout(cbind(grid, n_panel + 1L),
                     widths = c(rep(4 / n_col, n_col), 1))
  } else {
    graphics::layout(grid)
  }
  graphics::par(bg = theme$bg, fg = theme$fg, col.axis = theme$fg,
                col.lab = theme$fg, col.main = theme$fg)

  for (p in seq_len(n_panel)) {
    panel <- panels[[p]]
    at_col <- if (type == "matrix") panel$col else (p - 1L) %% n_col + 1L
    # A shared scale is stated once, in the leftmost panel that has one, so the
    # panels can sit against each other; a free one is stated by every panel.
    # The matrix triangle has a gap under its diagonal, so its leftmost drawn
    # panel is the one on it.
    y_annot <- free_scale || at_col == 1L ||
      (type == "matrix" && panel$row == panel$col)

    graphics::par(mar = c(5.1,
                          if (y_annot) 4.6 else 0.5,
                          if (strip) 3.1 else 4.1,
                          if (at_col == n_col) 2.1 else 0.5))

    sa_inter_panel(
      tbl      = panel$tbl,
      x_lv     = factor_lv[[panel$x]],
      trace_lv = factor_lv[[panel$trace]],
      ylim     = if (is.null(ylim) && free_scale) {
        sa_inter_span(panel$tbl, errorbar)
      } else if (is.null(ylim)) {
        span
      } else {
        ylim
      },
      xlab     = if (is.null(xlab)) panel$x else xlab,
      ylab     = if (y_annot) ylab else "",
      main     = if (outer_main) NULL else figure_main,
      col      = sa_inter_cols(col, length(factor_lv[[panel$trace]])),
      lwd      = lwd,
      errorbar = errorbar,
      theme    = theme,
      y_annot  = y_annot,
      key      = if (shared_legend) NULL else panel$trace,
      cex.axis = cex.axis,
      cex.lab  = cex.lab,
      cex.main = cex.main,
      cex.legend = cex.legend
    )

    if (strip) {
      graphics::mtext(panel$label, side = 3, line = 0.5, cex = cex.axis)
    }
  }

  if (shared_legend) {
    trace_lv <- factor_lv[[roles$trace]]
    n_lv <- length(trace_lv)
    graphics::par(mar = c(5, 0, 4, 1))
    graphics::plot.new()
    graphics::legend("center", title = roles$trace, legend = trace_lv,
                     col = sa_inter_cols(col, n_lv), pch = sa_inter_pch(n_lv),
                     lty = sa_inter_lty(n_lv), lwd = lwd, bty = "n",
                     cex = cex.legend, text.col = theme$fg,
                     title.col = theme$fg)
  }

  if (outer_main) {
    # Centred over the panels rather than over the device, whose right fifth is
    # the legend when there is one. `adj` positions the string across the outer
    # margin, so the fraction the panels take is the fraction to centre on.
    graphics::mtext(figure_main, side = 3, outer = TRUE, line = 0.8,
                    cex = cex.main, font = 2, col = theme$fg,
                    adj = if (shared_legend) 0.5 * 4 / 5 else 0.5)
  }

  invisible(NULL)
}


#' Draw one panel of traces
#'
#' @keywords internal
#' @noRd
sa_inter_panel <- function(tbl, x_lv, trace_lv, ylim, xlab, ylab, main, col,
                           lwd, errorbar, theme, y_annot, key, cex.axis,
                           cex.lab, cex.main, cex.legend) {
  n_x <- length(x_lv)
  n_lv <- length(trace_lv)
  pch <- sa_inter_pch(n_lv)
  lty <- sa_inter_lty(n_lv)
  # Room either side of the first and last level, so a point does not sit on the
  # edge of the panel and half of its symbol is not clipped away.
  limits <- c(1 - 0.25, n_x + 0.25)

  graphics::plot.default(
    NA, NA, xlim = limits, ylim = ylim, xaxt = "n", yaxt = "n", bty = "n",
    xlab = xlab, ylab = ylab, main = main,
    cex.lab = cex.lab, cex.main = cex.main
  )
  # Behind the traces, so a line is not cut into dashes by the grid it crosses.
  graphics::grid(nx = NA, ny = NULL, col = theme$guide, lty = 3, lwd = 1)

  # The axis line and the ticks in two calls, the first carrying the line across
  # the whole panel: the line plot.default() draws stops at the outermost tick,
  # which on a categorical axis leaves it short of the padding either side.
  graphics::axis(1, at = limits, labels = FALSE, lwd.ticks = 0)
  graphics::axis(1, at = seq_len(n_x), labels = x_lv, lwd = 0,
                 cex.axis = cex.axis)
  if (y_annot) {
    graphics::axis(2, las = 1, cex.axis = cex.axis)
  }

  for (i in seq_len(n_lv)) {
    at <- tbl$trace_level == trace_lv[i]
    if (!any(at)) next
    point <- tbl[at, , drop = FALSE]
    # In `x_lv` order rather than in the order the rows happen to sit in, so the
    # line joins neighbouring levels.
    point <- point[match(x_lv, point$x_level), , drop = FALSE]

    if (errorbar != "none") {
      has_bar <- is.finite(point$lower_conf) & is.finite(point$upper_conf) &
        point$upper_conf > point$lower_conf
      if (any(has_bar)) {
        graphics::arrows(
          seq_len(n_x)[has_bar], point$lower_conf[has_bar],
          seq_len(n_x)[has_bar], point$upper_conf[has_bar],
          length = 0.04, angle = 90, code = 3, col = col[i], lwd = 1
        )
      }
    }
    graphics::lines(seq_len(n_x), point$mean, col = col[i], lty = lty[i],
                    lwd = lwd)
    graphics::points(seq_len(n_x), point$mean, col = col[i], pch = pch[i],
                     cex = 1.3)
  }

  if (!is.null(key)) {
    graphics::legend("topleft", title = key, legend = trace_lv, col = col,
                     pch = pch, lty = lty, lwd = lwd, bty = "n",
                     cex = cex.legend * 0.85, text.col = theme$fg,
                     title.col = theme$fg)
  }
  invisible(NULL)
}


#' The y range the means and their bars need
#'
#' @keywords internal
#' @noRd
sa_inter_span <- function(tbl, errorbar) {
  vals <- if (errorbar == "none") {
    tbl$mean
  } else {
    c(tbl$mean, tbl$lower_conf, tbl$upper_conf)
  }
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    return(c(0, 1))
  }
  span <- range(vals)
  if (diff(span) == 0) {
    # A design where nothing moved is a flat line, which needs a range to be a
    # line in rather than a panel of zero height.
    return(span + c(-1, 1) * max(abs(span), 1) * 0.1)
  }
  span + c(-1, 1) * diff(span) * 0.08
}


#' What the y axis is measuring
#'
#' @keywords internal
#' @noRd
sa_inter_ylab <- function(drawn, input_scale) {
  # A point built from one cell is a cell mean; a point built from several is a
  # marginal mean, and calling it a cell mean would misname the averaging.
  what <- if (any(drawn$n_cells > 1L)) "marginal mean" else "cell mean"
  # The means are on whatever scale the comparison was handed, and a plot of log2
  # values labelled as if they were raw ones is read off by the wrong factor.
  if (identical(input_scale, "log2")) paste0(what, " (log2)") else what
}


#' The title the view describes itself with
#'
#' @keywords internal
#' @noRd
sa_inter_main <- function(panels, roles, type, strip) {
  if (type == "matrix") {
    return(paste0("Interactions of ", panels[[1L]]$tbl$features[1L]))
  }
  head <- paste0("Interaction of ", roles$trace, " and ", roles$x)
  if (type == "facet") {
    return(paste0(panels[[1L]]$tbl$features[1L], ": ", head, " by ",
                  roles$facet))
  }
  if (strip) {
    return(head)
  }
  paste0(panels[[1L]]$tbl$features[1L], ": ", head)
}


#' Colours, symbols and line types for the levels of a tracing factor
#'
#' Three channels rather than one, so the traces stay apart in greyscale and for
#' a reader who cannot tell the two colours apart.
#'
#' @keywords internal
#' @noRd
sa_inter_cols <- function(col, n) {
  if (is.null(col)) {
    grDevices::hcl.colors(n, "Dark 2")
  } else {
    rep_len(col, n)
  }
}


#' @keywords internal
#' @noRd
sa_inter_pch <- function(n) {
  rep_len(c(16, 15, 17, 18, 8, 4), n)
}


#' @keywords internal
#' @noRd
sa_inter_lty <- function(n) {
  rep_len(c(1, 2, 4, 3, 5, 6), n)
}
