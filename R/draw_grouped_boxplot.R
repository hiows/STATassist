#' Draw a grouped boxplot across several features
#'
#' Draws clusters of boxes with the levels of one factor side by side inside
#' each cluster and a legend in a narrow panel on the right. A single factor
#' gives one cluster per feature. A crossed design gives one panel per feature,
#' with the factors after the first along the x axis, so that the crossing sits
#' inside a panel where an interaction can be seen. Optionally returns the
#' summary statistics behind the boxes.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation and one column per feature.
#' @param feats Character vector of numeric column names in `data` to plot,
#'   in display order along the x axis.
#' @param group Grouping vector with one entry per row of `data`, for a single
#'   factor. Leave it `NULL` and name `factors` to draw a crossed design.
#' @param group_lv Character vector of at least two group levels, in the order
#'   they should appear inside each cluster. Rows belonging to any other level
#'   are dropped.
#' @param factors Named list of the crossed factors, each entry either the name
#'   of a column of `data` or a vector with one entry per row of it, exactly as
#'   [compare_factorial_groups()] takes them. There have to be at least two: one
#'   factor is `group`. The **first factor is the primary one**, the factor
#'   whose levels are the coloured boxes inside a cluster and the entries of the
#'   legend; the factors after it are crossed into the other axis.
#' @param factor_lv Named list giving the levels of each factor, with the
#'   reference level first, or `NULL` to take the levels from the data in sorted
#'   order. Naming it fixes the order the boxes, the clusters and the panels are
#'   drawn in and drops the rows belonging to any level it leaves out.
#' @param control_label The level to hold as the reference, one name per factor
#'   it points at, as a named list (`list(treatment = "control", sex = "male")`)
#'   or as a named character vector, exactly as [compare_factorial_groups()]
#'   takes it. The level it names is drawn first in its own factor, so a figure
#'   and the analysis it illustrates put the reference cell in the same place
#'   when they are given the same argument. Read only under `factors`: a single
#'   factor draws in `group_lv` order, which already says where its reference is.
#' @param panel_by Which axis the panels are over, read only under `factors`:
#'   `"feature"` for one panel per feature with the remaining factors along the
#'   x axis, or `"factor"` for one panel per combination of the remaining
#'   factors with the features along the x axis. See "What is drawn where" for
#'   why the first is the default.
#' @param panel_nrow How many rows the panels are laid out in, or `NULL` to let
#'   the arrangement decide: one row under `panel_by = "factor"`, and a grid as
#'   near square as the panel count allows under `"feature"`, where there is one
#'   panel per feature and a row of ten of them cannot be read.
#' @param gap Blank box widths inserted between neighbouring clusters.
#' @param lwd Line width of the boxes.
#' @param xlab,ylab,main Axis and title labels.
#' @param cex.lab,cex.axis,cex.main,cex.legend Character expansion for axis
#'   labels, axis annotation, the main title and the legend.
#' @param ylim Numeric length-2 y axis range, shared by every panel. `NULL`
#'   derives it, and what it derives depends on what a panel is. Panels over the
#'   factors hold the same features and share one range taken from every value
#'   drawn, since panels of the same quantity whose axes differ cannot be read
#'   against each other. Panels over the features hold different quantities,
#'   each with its own baseline, so each keeps the range [graphics::boxplot()]
#'   gives it and is annotated with its own axis.
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
#'       column per box: the group levels under `group`, and the cell labels of
#'       the design under `factors`. The bounds are the Tukey whisker fences
#'       `Q1 - 1.5 * IQR` and `Q3 + 1.5 * IQR`, not the drawn whisker ends.}
#'     \item{`median_confidence_stats`}{One data.frame per feature with rows
#'       `n`, `lower_conf`, `upper_conf`, the notch interval
#'       `median +/- 1.58 * IQR / sqrt(n)`. `n` counts non-missing values.}
#'   }
#'
#' @section What is drawn where:
#' A crossed design has three categorical axes to place and two dimensions to
#' place them in, so one of them has to become the panels. Which one decides
#' whether the picture shows an interaction.
#'
#' `panel_by = "feature"`, the default, gives one panel per feature. Inside it
#' the remaining factors run along the x axis and the primary factor is the
#' coloured boxes within each of those clusters, so **the two factors are side
#' by side in one panel**: an effect of the treatment that reverses between the
#' levels of the other factor is a pattern of colours that visibly flips a
#' couple of centimetres away, which is what an interaction is.
#'
#' `panel_by = "factor"` gives the transpose: one panel per combination of the
#' remaining factors, with the features along the x axis. It puts every feature
#' of one cell together, which is what to ask for when the question is about the
#' features rather than about the crossing. The cost is that the two factors are
#' then split between the legend and the panels, and reading an interaction
#' means comparing a colour profile in one panel against the same in another.
#'
#' Either way the boxes are the same boxes and the returned statistics are
#' identical; only the grouping of them into panels and clusters differs. The
#' cells, their order and their labels come from the same helpers
#' [compare_factorial_groups()] uses, so a box of this plot and a row of that
#' result are the same observations, and the columns of the returned statistics
#' are the cell labels the comparison and [simulate_factorial_groups()] both key
#' on. A cell holding no observation leaves its box blank rather than shifting
#' the ones beside it, and is reported in a message.
#'
#' @details
#' The function changes graphical parameters and the panel layout, and restores
#' both on exit, so the caller's device is left as it was found.
#'
#' @seealso [compare_two_groups()] and [compare_factorial_groups()] to test the
#'   same input.
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
#' ## A crossed design: one panel per feature, the sexes along the x axis and
#' ## the treatment the colours, so an effect that flips with sex is local.
#' sim <- simulate_factorial_groups(n_feats = 6, n_per_cell = 8, seed = 1)
#' cells <- draw_grouped_boxplot(
#'   data      = sim$args$data,
#'   feats     = sim$args$feats,
#'   factors   = sim$args$factors,
#'   factor_lv = sim$args$factor_lv,
#'   ylab      = "log2 abundance"
#' )
#' ## One column per cell, labelled the way the answer key labels them.
#' cells$box_summary_stats$prot_1
#'
#' ## The transpose: one panel per sex, every feature together inside it.
#' draw_grouped_boxplot(
#'   data      = sim$args$data,
#'   feats     = sim$args$feats,
#'   factors   = sim$args$factors,
#'   factor_lv = sim$args$factor_lv,
#'   panel_by  = "factor",
#'   ylab      = "log2 abundance"
#' )
#'
#' ## `control_label` draws the named level first in its own factor, and takes
#' ## the same argument the comparison does, so the two agree on where the
#' ## reference cell sits.
#' draw_grouped_boxplot(
#'   data          = sim$args$data,
#'   feats         = sim$args$feats[1:2],
#'   factors       = sim$args$factors,
#'   control_label = list(sex = "female"),
#'   ylab          = "log2 abundance"
#' )
#'
#' @export
draw_grouped_boxplot <- function(data,
                                 feats,
                                 group = NULL,
                                 group_lv = NULL,
                                 factors = NULL,
                                 factor_lv = NULL,
                                 control_label = NULL,
                                 panel_by = c("feature", "factor"),
                                 panel_nrow = NULL,
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

  panel_by <- match.arg(panel_by)
  if (!is.null(panel_nrow)) {
    sa_check_count(panel_nrow, "panel_nrow", 1)
  }
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

  input <- sa_box_input(data, feats, group, group_lv, factors, factor_lv,
                        control_label)

  sa_box_draw_panels(input, panel_by = panel_by, panel_nrow = panel_nrow,
                     gap = gap, lwd = lwd,
                     xlab = xlab, ylab = ylab, cex.lab = cex.lab,
                     cex.axis = cex.axis, cex.main = cex.main, ylim = ylim,
                     main = main, dark = dark, grid_lty = grid_lty,
                     grid_lwd = grid_lwd, cex.legend = cex.legend)

  if (!out_statistics) {
    return(invisible(NULL))
  }

  # One pass per feature and box; the two returned tables are slices of it.
  summaries <- lapply(input$samples, sa_box_stats)

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


#' Resolve what the boxes are, however the caller said it
#'
#' `group` and `factors` are two ways of saying the same thing, one factor and
#' several, and what comes out is the same either way: a box per feature and per
#' cell of the design, with the primary factor inside a cluster. So the two
#' paths meet here rather than in the drawing code, and a single factor comes
#' out as a design with exactly one group of cells.
#'
#' What this does not decide is the layout. A group of cells is a fact about the
#' design, not a panel: `sa_box_arrange()` turns it into panels one way or the
#' other.
#'
#' @param data,feats,group,group_lv,factors,factor_lv,control_label The
#'   arguments as received.
#'
#' @return List with `feats`, `lv` (the levels inside a cluster, which is what
#'   the legend lists), `samples` (one named list of vectors per feature, whose
#'   names label the boxes and become the columns of the returned statistics),
#'   `groups` (one entry per combination of the factors after the first, holding
#'   its `label` and the `cols` of `samples` it covers, in `lv` order) and
#'   `legend_title`.
#'
#' @keywords internal
#' @noRd
sa_box_input <- function(data, feats, group, group_lv, factors, factor_lv,
                         control_label = NULL) {
  said_group <- !is.null(group) || !is.null(group_lv)
  said_factors <- !is.null(factors) || !is.null(factor_lv)

  if (said_group && said_factors) {
    stop("`group` and `factors` are two ways of saying what the boxes are, so ",
         "a call takes one of them: `group` with `group_lv` for a single ",
         "factor, `factors` with `factor_lv` for a crossed design.",
         call. = FALSE)
  }
  if (!said_group && !said_factors) {
    stop("nothing says what the boxes are. Supply `group` and `group_lv` for ",
         "a single factor, or `factors` for a crossed design.", call. = FALSE)
  }
  if (!is.null(factor_lv) && is.null(factors)) {
    stop("`factor_lv` gives the levels of the factors `factors` holds, which ",
         "was not supplied.", call. = FALSE)
  }
  # A single factor has one list of levels and `group_lv` is it, so a reference
  # named a second way would be a second place for the draw order to be decided.
  if (said_group && !is.null(control_label)) {
    stop("`control_label` names a reference level per factor of a crossed ",
         "design, which `factors` states. A single factor draws in `group_lv` ",
         "order, so put the reference first there.", call. = FALSE)
  }

  if (said_group) {
    sa_box_one_factor(data, feats, group, group_lv)
  } else {
    sa_box_crossed(data, feats, factors, factor_lv, control_label)
  }
}


#' The single-factor case, as a design of one group of cells
#'
#' @keywords internal
#' @noRd
sa_box_one_factor <- function(data, feats, group, group_lv) {
  input <- sa_validate_wide_input(data, feats, group, group_lv,
                                  min_levels = 2L)
  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }
  lv <- levels(input$group)

  # One vector per feature and group level, the group levels consecutive within
  # each feature so the order matches the layout of `at`. Handing boxplot() the
  # list directly replaces a reshape to long format, which was the only reason
  # this package ever needed the `.grp`, `name` and `value` column names.
  samples <- lapply(input$feats, function(f) {
    stats::setNames(lapply(lv, function(l) input$data[[f]][input$group == l]),
                    lv)
  })
  names(samples) <- input$feats

  list(
    feats        = input$feats,
    lv           = lv,
    samples      = samples,
    # No label: with one factor there is nothing a strip would say that the
    # legend does not.
    groups       = list(list(label = NA_character_, cols = seq_along(lv))),
    legend_title = NULL
  )
}


#' The crossed case, one group of cells per combination of the factors after the
#' first
#'
#' The cells come from `sa_fact_layout()`, which `compare_factorial_groups()`
#' calls with the same arguments, so the boxes are the cells the analysis fits
#' and their labels are the ones the answer key uses. What is decided here is
#' only how the cells are dealt out to the groups: the primary factor varies
#' inside a cluster and the rest are read as a mixed-radix number, which is the
#' same arithmetic `sa_fact_cell_index()` numbers the cells with.
#'
#' @keywords internal
#' @noRd
sa_box_crossed <- function(data, feats, factors, factor_lv,
                           control_label = NULL) {
  input <- sa_validate_wide_input(data, feats, group = NULL, group_lv = NULL)
  data <- input$data
  feats <- input$feats

  design <- sa_fact_layout(data, factors, factor_lv, control_label)
  if (design$n_dropped > 0L) {
    message("Dropped ", design$n_dropped,
            " row(s) belonging to a level outside `factor_lv`.")
  }
  if (design$n_empty_cells > 0L) {
    message(design$n_empty_cells, " of ", design$n_cells,
            " cell(s) hold no observation, so their box is left blank: ",
            paste(design$cell_label[design$cell_n == 0L], collapse = ", "), ".")
  }

  primary <- names(design$factor_lv)[1L]
  lv <- design$factor_lv[[primary]]
  rest <- design$factor_lv[-1L]
  group_of_cell <- sa_fact_cell_index(as.matrix(design$cells[names(rest)]),
                                      vapply(rest, length, integer(1)))
  group_label <- sa_fact_cell_labels(rest, sa_fact_grid(rest))

  samples <- lapply(feats, function(f) {
    v <- data[[f]]
    stats::setNames(lapply(design$rows_of_cell, function(at) v[at]),
                    design$cell_label)
  })
  names(samples) <- feats

  groups <- lapply(seq_along(group_label), function(g) {
    at <- which(group_of_cell == g)
    list(label = group_label[g],
         cols = at[order(design$cells[[primary]][at])])
  })

  list(
    feats        = feats,
    lv           = lv,
    samples      = samples,
    groups       = groups,
    legend_title = primary
  )
}


#' Deal the boxes out into panels, one layout or its transpose
#'
#' `sa_box_input()` leaves a `(feature x cell)` matrix of boxes. A panel is a
#' choice of which way to read it: down the features, or across the groups of
#' cells. Both readings hold the same boxes and every panel of either holds the
#' same number of clusters, which is what lets one drawing routine take both.
#'
#' @param input As `sa_box_input()` returns.
#' @param panel_by `"feature"` or `"factor"`.
#'
#' @return List with `panels`, one entry per panel holding its `label` (the
#'   strip), its `cluster_labels` (the x axis annotation) and its `boxes` (one
#'   vector per box, clusters outside and the levels of the primary factor
#'   inside, which is the order `at` is laid out in), and `free_scale`, whether
#'   the panels hold different quantities and so cannot share a y axis.
#'
#' @keywords internal
#' @noRd
sa_box_arrange <- function(input, panel_by) {
  groups <- input$groups
  feats <- input$feats

  # A single factor is one group of cells, so there is nothing to panel over and
  # nothing for `panel_by` to choose between: it reads as "factor" either way.
  if (panel_by == "feature" && length(groups) > 1L) {
    labels <- vapply(groups, function(g) g$label, character(1))
    panels <- lapply(feats, function(f) {
      list(
        label = f,
        cluster_labels = labels,
        boxes = unlist(lapply(groups, function(g) input$samples[[f]][g$cols]),
                       recursive = FALSE, use.names = FALSE)
      )
    })
    # Each feature has its own baseline, so a common axis would flatten all of
    # them; see the `ylim` documentation.
    return(list(panels = panels, free_scale = TRUE))
  }

  panels <- lapply(groups, function(g) {
    list(
      label = g$label,
      cluster_labels = feats,
      boxes = unlist(lapply(input$samples, function(s) s[g$cols]),
                     recursive = FALSE, use.names = FALSE)
    )
  })
  list(panels = panels, free_scale = FALSE)
}


#' Draw the panels and the legend beside them
#'
#' One panel is the whole picture for a single factor, so there is one drawing
#' routine rather than one per mode, and both arrangements of a crossed design
#' reach it as the same `cluster_labels` and `boxes`. The panels are separate
#' plots on a `layout()` grid, so a y axis range they are meant to share has to
#' be settled here for all of them at once.
#'
#' @param input As `sa_box_input()` returns.
#' @param panel_by,panel_nrow The layout arguments as received, `panel_nrow`
#'   possibly `NULL`.
#'
#' @keywords internal
#' @noRd
sa_box_draw_panels <- function(input, panel_by, panel_nrow, gap, lwd, xlab,
                               ylab, cex.lab, cex.axis, cex.main, ylim, main,
                               dark, grid_lty, grid_lwd, cex.legend) {
  arranged <- sa_box_arrange(input, panel_by)
  panels <- arranged$panels
  n_lv <- length(input$lv)
  n_panel <- length(panels)
  n_cluster <- length(panels[[1L]]$cluster_labels)

  at <- lapply(seq_len(n_cluster), function(i) {
    start <- (i - 1) * (n_lv + gap)
    start + seq_len(n_lv)
  })

  # Panels of the same quantity are separate plots, so each would otherwise
  # scale to its own values and the difference between two panels would be
  # invisible. Panels of different quantities are the opposite case: a shared
  # range flattens every one of them, so each keeps its own and says so on its
  # own axis. A lone panel is left to boxplot(), whose range carries padding
  # this one would not.
  free_scale <- arranged$free_scale && is.null(ylim)
  if (is.null(ylim) && !free_scale && n_panel > 1L) {
    vals <- unlist(input$samples, use.names = FALSE)
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) {
      stop("`feats` hold no finite value in any cell, so there is nothing to ",
           "draw.", call. = FALSE)
    }
    ylim <- range(vals)
  }

  strip <- n_panel > 1L
  if (is.null(panel_nrow)) {
    # One panel per feature is many panels, and a single row of them is a strip
    # too thin to read, so the default shape is as near square as it gets.
    panel_nrow <- if (arranged$free_scale) {
      max(1L, round(sqrt(n_panel)))
    } else {
      1L
    }
  }
  n_row <- min(panel_nrow, n_panel)
  n_col <- ceiling(n_panel / n_row)
  slots <- c(seq_len(n_panel), rep(0L, n_row * n_col - n_panel))

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

  # A title belongs to the figure rather than to a panel of it, so with panels
  # it moves to the outer margin and is written once.
  outer_main <- strip && !is.null(main)
  if (outer_main) {
    graphics::par(oma = c(0, 0, 3, 0))
  }

  graphics::layout(
    cbind(matrix(slots, nrow = n_row, byrow = TRUE), n_panel + 1L),
    widths = c(rep(4 / n_col, n_col), 1)
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

  for (p in seq_len(n_panel)) {
    at_col <- (p - 1L) %% n_col + 1L
    boxes <- panels[[p]]$boxes

    # A shared scale is stated once, in the first column, so that the panels can
    # sit against each other; a free one has to be stated by every panel.
    y_annot <- free_scale || at_col == 1L

    graphics::par(mar = c(5.1,
                          if (y_annot) 4.1 else 0.5,
                          if (strip) 3.1 else 4.1,
                          if (at_col == n_col) 2.1 else 0.5))

    # Drawn twice with the grid in between so the grid sits behind the boxes.
    # One closure keeps the two calls from drifting apart.
    draw_boxes <- function(add) {
      graphics::boxplot(
        boxes,
        at       = unlist(at),
        xlab     = xlab,
        ylab     = if (y_annot) ylab else NULL,
        xaxt     = "n",
        yaxt     = if (y_annot) graphics::par("yaxt") else "n",
        cex.lab  = cex.lab,
        cex.axis = cex.axis,
        cex.main = cex.main,
        ylim     = ylim,
        main     = if (outer_main) NULL else main,
        frame    = FALSE,
        col      = bg_cols,
        border   = cols,
        lwd      = lwd,
        pch      = 16,
        add      = add
      )
    }

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
      labels = panels[[p]]$cluster_labels,
      cex.axis = cex.axis,
      tick = FALSE
    )

    if (strip) {
      graphics::mtext(panels[[p]]$label, side = 3, line = 0.5, cex = cex.axis)
    }
  }

  graphics::par(mar = c(5, 0, 4, 1))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = input$lv,
    fill = cols,
    border = "white",
    bty = "n",
    cex = cex.legend,
    title = input$legend_title
  )

  if (outer_main) {
    # Centred over the panels rather than over the device, whose right fifth is
    # the legend. `adj` positions the string across the outer margin, so the
    # fraction the panels take is the fraction to centre on.
    graphics::mtext(main, side = 3, outer = TRUE, line = 0.8, cex = cex.main,
                    font = 2, adj = 0.5 * 4 / 5)
  }

  invisible(NULL)
}


#' The numbers behind one feature's boxes
#'
#' @param values Named list of one vector per box, in the order the columns of
#'   the returned tables are in.
#'
#' @return Matrix of ten named rows, one column per box.
#'
#' @keywords internal
#' @noRd
sa_box_stats <- function(values) {
  vapply(names(values), function(k) {
    v <- values[[k]]
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
}
