#' Draw a volcano plot
#'
#' Plots `log2fc` against `-log10(pvalue)` for the table
#' [estimate_significance()] returns, colours the features that clear both
#' cutoffs by direction and labels the strongest of them. A term reading of a
#' factorial comparison is drawn as one panel per term in a single figure.
#'
#' @param significance_result The object returned by
#'   [estimate_significance()], whose `significance` element is what is plotted.
#'   With `by = "contrast"` that element holds one table per contrast, so name
#'   the one to draw: `sig$significance[["virginica - setosa"]]`. With
#'   `by = "term"` the whole list is drawn, one panel per term. A bare verdict
#'   data.frame is accepted too.
#' @param terms Which terms get a panel, read only under a term reading. `NULL`,
#'   the default, draws the two main effects of the first two factors and their
#'   interaction, which is every term of a two-factor design and three of the
#'   seven of a three-factor one; the terms left out are named in a `message()`.
#'   Give term labels, as `$terms$terms` spells them, to draw those instead.
#' @param panel_nrow How many rows the panels are laid out in, or `NULL` for a
#'   single row, which is what three panels of a two-factor design want.
#' @param use_adjusted Logical. If `TRUE`, plot and threshold the `adj_pvalue`
#'   column; if `FALSE`, use the unadjusted `pvalue` column. The axis label
#'   follows, so the y axis always describes what was actually plotted.
#' @param log2fc_cutoff,pval_cutoff Cutoffs for calling a feature changed and
#'   significant, drawn as guides. `NULL`, the default, takes the values
#'   [estimate_significance()] recorded on `significance_result`, so the guides
#'   agree with its `is_signif` column. Supply a number to override them.
#'   Selecting a subset of rows keeps the recorded values, but selecting
#'   columns drops them, in which case both must be given.
#' @param anno_feats Logical. If `TRUE`, label the strongest significant
#'   features. A run where no feature clears both cutoffs still draws the plot,
#'   with a `message()` in place of the labels.
#' @param anno_top How many features to label in each direction, so up to
#'   `2 * anno_top` labels in total.
#' @param cex.anno Character expansion for those labels.
#' @param xlim,ylim Numeric length-2 axis ranges, or `NULL` to derive them from
#'   the data. A supplied range is used as given. Across panels a derived range
#'   is derived from all of them at once, so the panels can be compared.
#' @param xlab X axis label, or `NULL` to derive it from what `log2fc` compares
#'   in the comparison behind the verdict.
#' @param main Plot title. With panels it is the title of the figure, written
#'   once above them, and each panel is titled with its term.
#' @param cex.lab,cex.axis,cex.main Character expansion for the axis labels,
#'   the axis annotation and the title.
#' @param margin Plot margins in lines, passed to [graphics::par()] as `mar`.
#'   Applied to every panel.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return `NULL` invisibly.
#'
#' @details
#' The plot margins are restored on exit, but the coordinate system is
#' deliberately left in place, so [graphics::points()], [graphics::text()] and
#' friends can still be used to add to the finished plot. With panels what is
#' left in place is the last panel's coordinate system, and the panel grid itself
#' is undone, `mfrow` included, since a figure of several plots is finished when
#' the last of them is drawn.
#'
#' A p-value of exactly zero has no finite `-log10()`. Rather than dropping the
#' most significant points or letting an infinite axis limit blank the plot, the
#' axis is scaled to the largest finite value and those points are drawn at the
#' top of it. Non-finite `log2fc`, which [compare_two_groups()] produces when a
#' group centre is zero, is capped at the edge of the x axis the same way. Either
#' kind of capping is reported in a `message()`, since a capped point no longer
#' sits at its true coordinate.
#'
#' What `log2fc` compares is not the same question in every scenario, so the x
#' axis label follows the comparison the verdict came from. A multi-group omnibus
#' verdict carries the single fold change its `effect` table holds, the level
#' furthest from the reference rather than any named pair, and the label says so.
#' A factorial omnibus verdict does the same with cells: the reference cell is
#' where every factor sits at its first level, and the label names it. A term
#' table carries an ANOVA component rather than a ratio of two centres, so
#' it is labelled as an effect and not as a fold change; see "The size of a term"
#' in [compare_factorial_groups()] for how to read its magnitude. Every other
#' reading compares two fixed centres and is labelled `log2 FC`. Pass `xlab` to
#' override this.
#'
#' Points are coloured by the same masks that select the labels, so which
#' features are highlighted and which are labelled can never disagree. With the
#' default arguments those masks reproduce the `is_signif` column of the input.
#' The labels are drawn in a brighter shade than the points on purpose, so that
#' a label stays legible where it overlaps them.
#'
#' @section One panel per term:
#' A term reading is drawn as a figure rather than asked to name one table,
#' unlike a contrast reading. A crossed design decomposes into a fixed and small
#' set of terms, three for two factors, and the whole point of the reading is to
#' see which of them a feature responded to, which is a comparison between the
#' panels. The number of pairwise contrasts, by contrast, follows from the level
#' counts and is arbitrary.
#'
#' Every panel is judged by one rule, the cutoffs of the first table, and shares
#' the axes with the others unless `xlim` or `ylim` says otherwise. Both are what
#' makes a point in one panel comparable with a point in another.
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
#' ## A multi-group verdict says on the x axis which levels its log2fc compares
#' multi <- compare_multiple_groups(iris, c("Sepal.Length", "Petal.Length"),
#'                                  iris$Species, levels(iris$Species))
#' draw_volcano_plot(estimate_significance(multi, log2fc_cutoff = 0.1))
#'
#' ## A factorial omnibus verdict names the reference cell on the x axis
#' fact <- compare_factorial_groups(
#'   data    = warpbreaks,
#'   feats   = "breaks",
#'   factors = list(wool = "wool", tension = "tension"),
#'   posthoc = FALSE
#' )
#' draw_volcano_plot(estimate_significance(fact, log2fc_cutoff = 0.1),
#'                   main = "warpbreaks")
#'
#' ## One panel per term of a crossed design: wool, tension, wool:tension
#' draw_volcano_plot(estimate_significance(fact, by = "term", log2fc_cutoff = 0.1),
#'                   main = "warpbreaks, by term")
#'
#' @export
draw_volcano_plot <- function(significance_result,
                              terms = NULL,
                              panel_nrow = NULL,
                              use_adjusted = TRUE,
                              log2fc_cutoff = NULL,
                              pval_cutoff = NULL,
                              anno_feats = TRUE,
                              anno_top = 10,
                              cex.anno = 1,
                              xlim = NULL,
                              ylim = NULL,
                              xlab = NULL,
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
  if (!is.null(panel_nrow)) {
    sa_check_count(panel_nrow, "panel_nrow", 1)
  }
  if (!is.null(terms) &&
        (!is.character(terms) || length(terms) == 0L || anyNA(terms))) {
    stop("`terms` must be NULL or a character vector of term labels.",
         call. = FALSE)
  }

  # The verdict object carries the table beside the scenario name; the table on
  # its own is still accepted, since selecting rows from it produces one.
  if (inherits(significance_result, "sa_significance")) {
    significance_result <- significance_result$significance
  }
  if (!is.data.frame(significance_result)) {
    tables <- sa_volcano_term_tables(significance_result)
    return(sa_volcano_panels(
      tables[sa_volcano_terms(tables, terms)],
      panel_nrow = panel_nrow, use_adjusted = use_adjusted,
      log2fc_cutoff = log2fc_cutoff, pval_cutoff = pval_cutoff,
      anno_feats = anno_feats, anno_top = anno_top, cex.anno = cex.anno,
      xlim = xlim, ylim = ylim, xlab = xlab, main = main, cex.lab = cex.lab,
      cex.axis = cex.axis, cex.main = cex.main, margin = margin, ...
    ))
  }

  sa_volcano_one(significance_result, use_adjusted = use_adjusted,
                 log2fc_cutoff = log2fc_cutoff, pval_cutoff = pval_cutoff,
                 anno_feats = anno_feats, anno_top = anno_top,
                 cex.anno = cex.anno, xlim = xlim, ylim = ylim, xlab = xlab,
                 main = main, cex.lab = cex.lab, cex.axis = cex.axis,
                 cex.main = cex.main, margin = margin, ...)
}


#' The verdict tables of a term reading, or a refusal
#'
#' A term reading is the one list this function draws whole. A contrast reading
#' is sent back to name a table, because a list of contrasts is as long as the
#' level counts make it while the terms of a design are a fixed and small set.
#'
#' @param x The `significance` element, known not to be a data.frame.
#'
#' @return The list, named by the `term` attribute of each table, which is what
#'   the `terms` argument selects on.
#'
#' @keywords internal
#' @noRd
sa_volcano_term_tables <- function(x) {
  if (!is.list(x) || length(x) == 0L ||
        !all(vapply(x, is.data.frame, logical(1)))) {
    stop("`significance_result` must be the object returned by ",
         "estimate_significance().", call. = FALSE)
  }

  labels <- vapply(x, function(tbl) {
    term <- attr(tbl, "term")
    if (is.null(term)) NA_character_ else as.character(term)[1]
  }, character(1))

  if (anyNA(labels)) {
    stop("`significance_result` holds one verdict table per contrast, and a ",
         "volcano plot draws one of them. Name it: ",
         "`sig$significance[[\"", names(x)[1], "\"]]`.", call. = FALSE)
  }

  names(x) <- labels
  x
}


#' Which terms earn a panel
#'
#' @param tables Named term tables, in the order `$terms` lists them.
#' @param terms The argument as received, possibly `NULL`.
#'
#' @keywords internal
#' @noRd
sa_volcano_terms <- function(tables, terms) {
  labels <- names(tables)

  if (!is.null(terms)) {
    unknown <- setdiff(terms, labels)
    if (length(unknown) > 0L) {
      stop("`terms` names term(s) the verdict does not hold: ",
           paste(unknown, collapse = ", "), ". It holds ",
           paste(labels, collapse = ", "), ".", call. = FALSE)
    }
    return(unique(terms))
  }

  # The first two factors are the first two main effects, since `$terms` lists
  # the main effects first and in the order the factors were declared.
  orders <- vapply(tables, function(tbl) as.integer(attr(tbl, "term_order")),
                   integer(1))
  mains <- labels[orders == 1L]
  if (length(mains) < 2L) {
    return(labels)
  }

  pair <- mains[1:2]
  wanted <- intersect(c(pair, paste(pair, collapse = ":")), labels)
  left_out <- setdiff(labels, wanted)
  if (length(left_out) > 0L) {
    message("Drew the terms of ", pair[1], " and ", pair[2],
            ", leaving out ", paste(left_out, collapse = ", "),
            ". Name terms in `terms` to draw those instead.")
  }
  wanted
}


#' One volcano plot per term, on one figure
#'
#' The panels share the rule they are judged by and, unless told otherwise, the
#' axes they are drawn on, which is what makes them comparable rather than three
#' plots that happen to sit together. The drawing itself is `sa_volcano_one()`,
#' the same routine a single plot goes through.
#'
#' @param tables The term tables to draw, already selected and in panel order.
#'
#' @keywords internal
#' @noRd
sa_volcano_panels <- function(tables, panel_nrow, use_adjusted, log2fc_cutoff,
                             pval_cutoff, anno_feats, anno_top, cex.anno, xlim,
                             ylim, xlab, main, cex.lab, cex.axis, cex.main,
                             margin, ...) {
  p_col <- sa_volcano_p_col(use_adjusted)
  for (nm in names(tables)) {
    sa_volcano_check_cols(tables[[nm]], p_col)
  }

  # One rule for the figure, read off the first panel. Every table of a reading
  # came out of one estimate_significance() call, so there is nothing to
  # reconcile; taking them panel by panel would let the guides move between
  # panels drawn on shared axes.
  cutoffs <- sa_volcano_cutoffs(tables[[1]], log2fc_cutoff, pval_cutoff)
  lims <- sa_volcano_lims(tables, p_col, cutoffs$log2fc, cutoffs$pval,
                          xlim, ylim)

  n_panel <- length(tables)
  n_row <- min(if (is.null(panel_nrow)) 1L else panel_nrow, n_panel)
  n_col <- ceiling(n_panel / n_row)
  slots <- c(seq_len(n_panel), rep(0L, n_row * n_col - n_panel))

  # Only the parameters this function sets are put back, not a blanket
  # par(no.readonly = TRUE) snapshot. That snapshot also carries `fin`, `pin`
  # and `mai`, which are absolute sizes: restoring them pins the figure to the
  # size this plot happened to be drawn at, so the next plot on a device that
  # has since been resized is redrawn small in a corner of it. `mfrow` comes
  # back with the rest, since `layout()` overwrites whatever grid the caller had
  # set up, and `oma` because a title over several panels lives there.
  old_par <- graphics::par(c("mar", "mfrow", "oma"))
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  # A title belongs to the figure rather than to a panel of it, so it moves to
  # the outer margin and the panels keep their term names.
  if (!is.null(main)) {
    graphics::par(oma = c(0, 0, 3, 0))
  }
  graphics::layout(matrix(slots, nrow = n_row, byrow = TRUE))

  for (nm in names(tables)) {
    sa_volcano_one(tables[[nm]], use_adjusted = use_adjusted,
                   log2fc_cutoff = cutoffs$log2fc,
                   pval_cutoff = cutoffs$pval, anno_feats = anno_feats,
                   anno_top = anno_top, cex.anno = cex.anno, xlim = lims$xlim,
                   ylim = lims$ylim, xlab = xlab, main = nm, cex.lab = cex.lab,
                   cex.axis = cex.axis, cex.main = cex.main, margin = margin,
                   ...)
  }

  if (!is.null(main)) {
    graphics::mtext(main, side = 3, outer = TRUE, line = 0.8, cex = cex.main,
                    font = 2)
  }

  invisible(NULL)
}


#' Which p-value column a reading is drawn from
#'
#' @keywords internal
#' @noRd
sa_volcano_p_col <- function(use_adjusted) {
  if (use_adjusted) "adj_pvalue" else "pvalue"
}


#' The columns a verdict table has to carry to be plotted
#'
#' @keywords internal
#' @noRd
sa_volcano_check_cols <- function(tbl, p_col) {
  absent <- setdiff(c("features", "log2fc", p_col), names(tbl))
  if (length(absent) > 0L) {
    stop("`significance_result` is missing the column(s) ",
         paste(absent, collapse = ", "),
         ". Pass the table returned by estimate_significance().",
         call. = FALSE)
  }
  invisible(NULL)
}


#' The rule a plot draws its guides for
#'
#' Falling back to the recorded cutoffs is what keeps the guides on the plot and
#' the verdict in the table describing the same rule.
#'
#' @keywords internal
#' @noRd
sa_volcano_cutoffs <- function(tbl, log2fc_cutoff, pval_cutoff) {
  if (is.null(log2fc_cutoff)) {
    log2fc_cutoff <- attr(tbl, "log2fc_cutoff")
  }
  if (is.null(pval_cutoff)) {
    pval_cutoff <- attr(tbl, "pval_cutoff")
  }
  if (is.null(log2fc_cutoff) || is.null(pval_cutoff)) {
    stop("`significance_result` does not carry the cutoffs ",
         "estimate_significance() records, so `log2fc_cutoff` and ",
         "`pval_cutoff` must be supplied. Selecting columns from the table ",
         "drops them.", call. = FALSE)
  }
  sa_check_scalar_num(log2fc_cutoff, "log2fc_cutoff", 0)
  sa_check_scalar_num(pval_cutoff, "pval_cutoff", 0, 1, lower_open = TRUE)
  list(log2fc = log2fc_cutoff, pval = pval_cutoff)
}


#' Axis ranges wide enough for every table they cover
#'
#' Given several tables the range covers all of them, so that panels meant to be
#' compared are drawn on one scale. A supplied range is returned as given, and a
#' figure whose ranges are both supplied is never asked whether it has anything
#' finite to derive them from.
#'
#' @keywords internal
#' @noRd
sa_volcano_lims <- function(tables, p_col, log2fc_cutoff, pval_cutoff, xlim,
                            ylim) {
  if (!is.null(xlim) && !is.null(ylim)) {
    return(list(xlim = xlim, ylim = ylim))
  }

  pull <- function(f) unlist(lapply(tables, f), use.names = FALSE)
  neglog_p <- -log10(pull(function(tbl) tbl[[p_col]]))
  log2fc <- pull(function(tbl) tbl$log2fc)
  y_finite <- neglog_p[is.finite(neglog_p)]
  x_finite <- log2fc[is.finite(log2fc)]
  if (length(y_finite) == 0L || length(x_finite) == 0L) {
    stop("nothing can be plotted: no feature has both a finite `log2fc` and a ",
         "finite -log10(`", p_col, "`).", call. = FALSE)
  }

  # A run where every p-value is 1 leaves the top at 0, which is not a usable
  # axis, so the guide line height sets the floor.
  y_top <- max(max(y_finite), -log10(pval_cutoff), 1)
  x_max <- max(max(abs(x_finite)), log2fc_cutoff)

  list(xlim = if (is.null(xlim)) c(-x_max, x_max) * 1.05 else xlim,
       ylim = if (is.null(ylim)) c(0, y_top * 1.1) else ylim)
}


#' Draw one volcano plot
#'
#' The whole of the drawing, so that a single plot and a panel of a figure cannot
#' come to differ in how a point is placed, coloured or labelled.
#'
#' @param tbl One verdict table.
#' @param log2fc_cutoff,pval_cutoff Possibly `NULL`, in which case the table's
#'   own record of them is used.
#' @param xlim,ylim Possibly `NULL`, in which case they are derived from `tbl`
#'   alone.
#'
#' @keywords internal
#' @noRd
sa_volcano_one <- function(tbl, use_adjusted, log2fc_cutoff, pval_cutoff,
                           anno_feats, anno_top, cex.anno, xlim, ylim, xlab,
                           main, cex.lab, cex.axis, cex.main, margin, ...) {
  p_col <- sa_volcano_p_col(use_adjusted)
  sa_volcano_check_cols(tbl, p_col)

  cutoffs <- sa_volcano_cutoffs(tbl, log2fc_cutoff, pval_cutoff)
  log2fc_cutoff <- cutoffs$log2fc
  pval_cutoff <- cutoffs$pval

  feats <- as.character(tbl$features)
  sa_check_feat_names(feats)
  log2fc <- tbl$log2fc
  pvalue <- tbl[[p_col]]
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
  lims <- sa_volcano_lims(list(tbl), p_col, log2fc_cutoff, pval_cutoff,
                          xlim, ylim)
  xlim <- lims$xlim
  ylim <- lims$ylim
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
  if (is.null(xlab)) {
    xlab <- sa_volcano_xlab(tbl)
  }

  graphics::plot(
    plot_x,
    plot_y,
    xlab = xlab,
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
    idx <- c(up_idx, down_idx)

    # text() treats a zero-length `labels` as an error rather than as nothing to
    # draw, so a run where no feature clears both cutoffs has to stop here.
    if (length(idx) == 0L) {
      message("No feature clears both cutoffs, so nothing was labelled.")
    } else {
      graphics::text(
        plot_x[idx],
        plot_y[idx] + label_offset,
        labels = feats[idx],
        cex = cex.anno,
        col = c(rep(up_label, length(up_idx)),
                rep(down_label, length(down_idx)))
      )
    }
  }

  invisible(NULL)
}


#' The x axis label a verdict table earns
#'
#' Two readings do not compare two centres the caller named, and saying so on the
#' axis is what keeps either from being read as a two-group plot. A multi-group
#' omnibus verdict holds the level furthest from the reference, and which level
#' that is differs per feature. A factorial omnibus verdict does the same with
#' cells rather than levels. A term table holds an ANOVA component, a
#' deviation from what the rest of the model predicts rather than a ratio at all.
#' A contrast table of the same comparison does compare a named pair, so it is
#' left alone.
#'
#' @keywords internal
#' @noRd
sa_volcano_xlab <- function(tbl) {
  term <- attr(tbl, "term")
  if (!is.null(term)) {
    return(as.expression(bquote(log[2] ~ effect ~ .(paste0("(", term, ")")))))
  }
  if (!is.null(attr(tbl, "contrast"))) {
    return(expression(log[2] ~ FC))
  }
  unit <- switch(attr(tbl, "analysis"),
                 multi_group_comparison = "level",
                 factorial_comparison   = "cell",
                 NULL)
  if (is.null(unit)) {
    return(expression(log[2] ~ FC))
  }
  reference <- attr(tbl, "group_lv")[1]
  named <- if (length(reference) == 1L && !is.na(reference)) {
    paste0("(most extreme ", unit, " vs ", reference, ")")
  } else {
    paste0("(most extreme ", unit, " vs reference)")
  }
  as.expression(bquote(log[2] ~ FC ~ .(named)))
}
