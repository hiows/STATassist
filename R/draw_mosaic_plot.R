#' Draw a mosaic plot of a contingency table
#'
#' Splits the x axis by the first variable's marginal shares and each strip by the
#' second variable's conditional shares, so the area of a tile is the cell's share
#' of the table. Three things are drawn on top of that geometry, and each of them
#' answers a question the geometry alone leaves open.
#'
#' The **shading** says which cells made the statistic what it is. It reads the
#' Pearson residual, the quantity that squares and sums to that statistic, at the
#' conventional cuts of 2 and 4. The two colours are the ones
#' [draw_volcano_plot()] already uses for a feature that moved up and one that
#' moved down, so "more than expected" and "less than expected" are the same pair
#' of colours the rest of the package reads as a direction.
#'
#' The **expected line** says what "expected" was. A dashed segment sits at each
#' boundary the tiles of a strip would have had under the null hypothesis, so the
#' departure is the distance between a tile edge and the line beside it rather
#' than something to be inferred by comparing strips by eye. Under independence
#' the lines fall at the same heights in every strip, which is what makes an
#' association visible at a glance; under symmetry they do not, because there the
#' expectation is a cell against its own transpose.
#'
#' The **annotation** says how many observations a tile stands for, which area
#' cannot: a wide short tile and a narrow tall one can hold the same count.
#'
#' @section Which null hypothesis is drawn:
#' The one the result was tested against,
#' `categorical_comparison_result$design$null`. That is the whole point of
#' reading it off the result rather than recomputing it here: a matched design
#' is tested for symmetry, so a mosaic of it shaded by departure from
#' independence would be a picture of a hypothesis nothing in the result has a
#' p-value for. A bare [table()] carries no such hypothesis, which is why one is
#' not accepted here: [compare_categorical_groups()] is what settles the null,
#' the levels and their order, and this function draws what it settled.
#'
#' `residual = "standardized"` is refused under symmetry. The variance correction
#' that residual divides by is derived for a table held against its own margins
#' and `$cells$std_residual` is `NA` there, so the request has no answer rather
#' than a different one.
#'
#' @param categorical_comparison_result A categorical comparison result, as
#'   returned by [compare_categorical_groups()]. The levels that take part, their
#'   order and the null hypothesis the shading is read under are all settled
#'   there, so there is no `category_lv` or `control_label` to restate here.
#' @param shade Logical. If `FALSE`, every tile is drawn in the background
#'   colour and the residual key is omitted.
#' @param residual Which residual the shading reads. `"pearson"` squares and sums
#'   to the test statistic, so it says how the statistic was made up.
#'   `"standardized"` is referred to a standard normal, so it says which cells are
#'   individually surprising, and it is only defined when the null is a statement
#'   about the margins.
#' @param expected_line Logical. If `TRUE`, mark each strip at the tile boundaries
#'   the null hypothesis expects.
#' @param anno_cells What to write on a tile. `"count"` is the observed count,
#'   `"percent"` the tile's share of its own strip, `"both"` puts one over the
#'   other, and `"none"` writes nothing. `"auto"`, the default, writes as much as
#'   the tile has room for, measured against the label rather than against a fixed
#'   fraction of the plot. `TRUE` and `FALSE` are accepted as `"auto"` and
#'   `"none"`.
#' @param gap Gap between neighbouring tiles, as a fraction of the axis, for each
#'   gap rather than for all of them together. Capped so that the gaps never take
#'   more than two fifths of an axis however many levels there are.
#' @param xlab,ylab,main Axis labels and title. `NULL` takes the variable names
#'   from the table, or no title.
#' @param cex.lab,cex.axis,cex.main,cex.legend,cex.anno Character expansion for
#'   the axis labels, the level names, the title, the residual key and the tile
#'   annotation.
#' @param dark Logical. If `TRUE`, draw on a dark background with light
#'   annotation, the same palette [draw_grouped_boxplot()] uses.
#'
#' @return Invisibly, a list of the picture as it was drawn.
#'
#'   \describe{
#'     \item{`cells`}{The cell table with `x1`, `x2`, `y1`, `y2` and `fill` added,
#'       in the order the tiles were painted.}
#'     \item{`widths`}{Marginal share of each strip, named by row level. These are
#'       the widths before the gaps are taken out of them.}
#'     \item{`heights`}{Conditional share of each tile within its strip, as a
#'       matrix of row level by column level.}
#'     \item{`expected_prop`}{The same shares the null hypothesis expects, which
#'       is what the dashed segments were drawn from.}
#'     \item{`empty_levels`}{`row` and `col` levels that hold no observation, so
#'       drew no tile and took no axis label.}
#'     \item{`null`}{Which hypothesis the shading and the segments are about.}
#'     \item{`residual`,`residual_breaks`,`colors`}{The scale the shading read.}
#'   }
#'
#' @details
#' The level names on the y axis are read off the **reference strip**, the first
#' one, which the `control_label` and `category_lv` of
#' [compare_categorical_groups()] decide. No single set
#' of positions can label every strip, since the whole content of a mosaic is that
#' the strips are cut at different heights, so labelling one of them and saying
#' which is the honest version of the choice. The first strip is the one the rest
#' of the package already treats as the reference.
#'
#' Only the `par()` values this function sets are put back. A blanket
#' `par(no.readonly = TRUE)` snapshot also carries `fin`, `pin` and `mai`, which
#' are absolute sizes, and restoring them pins the next plot to the size this
#' one happened to be drawn at.
#'
#' @seealso [compare_categorical_groups()] for the analysis this draws, and
#'   [draw_grouped_boxplot()] for the numeric counterpart.
#'
#' @examples
#' smoking <- data.frame(
#'   smoker = rep(c("y", "n"), each = 60),
#'   grade  = c(rep(c("high", "mid", "low"), c(10, 20, 30)),
#'              rep(c("high", "mid", "low"), c(30, 20, 10)))
#' )
#' res <- compare_categorical_groups(smoking)
#' draw_mosaic_plot(res)
#'
#' ## The dashed segments are what independence expected, so the departure is a
#' ## distance rather than a comparison between strips.
#' drawn <- draw_mosaic_plot(res, anno_cells = "both")
#' drawn$expected_prop
#'
#' ## The levels that take part and the order they sit in are settled by the
#' ## comparison, so a different reference is a different call to it.
#' draw_mosaic_plot(
#'   compare_categorical_groups(smoking, control_label = c(smoker = "y"))
#' )
#'
#' ## A matched design is shaded by departure from symmetry, so the diagonal is
#' ## neutral by construction and only the discordant cells carry colour.
#' before_after <- data.frame(
#'   before = rep(c("pass", "fail"), c(20, 30)),
#'   after  = c(rep(c("pass", "fail"), c(18, 2)), rep(c("pass", "fail"), c(14, 16)))
#' )
#' draw_mosaic_plot(compare_categorical_groups(before_after, paired = TRUE))
#'
#' @export
draw_mosaic_plot <- function(categorical_comparison_result,
                             shade = TRUE,
                             residual = c("pearson", "standardized"),
                             expected_line = TRUE,
                             anno_cells = c("auto", "count", "percent", "both",
                                            "none"),
                             gap = 0.015,
                             xlab = NULL,
                             ylab = NULL,
                             main = NULL,
                             cex.lab = 1.3,
                             cex.axis = 1.2,
                             cex.main = 1.3,
                             cex.legend = 1.1,
                             cex.anno = 1,
                             dark = FALSE) {

  # The argument used to be a flag, and the two flags still say what they said.
  if (is.logical(anno_cells) && length(anno_cells) == 1L && !is.na(anno_cells)) {
    anno_cells <- if (anno_cells) "auto" else "none"
  }
  anno_cells <- match.arg(anno_cells)
  residual <- match.arg(residual)

  sa_check_flag(shade, "shade")
  sa_check_flag(expected_line, "expected_line")
  sa_check_flag(dark, "dark")
  sa_check_scalar_num(gap, "gap", 0, 0.2)
  sa_check_scalar_num(cex.lab, "cex.lab", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.axis, "cex.axis", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.main, "cex.main", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.legend, "cex.legend", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.anno, "cex.anno", 0, lower_open = TRUE)

  input <- sa_mosaic_input(categorical_comparison_result)
  if (identical(residual, "standardized") &&
        identical(input$null, "symmetry")) {
    stop("`residual = \"standardized\"` has no value under symmetry, which is ",
         "the null this result was tested against: the variance correction it ",
         "divides by is derived for a table held against its own margins, so ",
         "`$cells$std_residual` is NA here. Use `residual = \"pearson\"`, whose ",
         "squares sum to McNemar's statistic.", call. = FALSE)
  }

  layout_out <- sa_mosaic_layout(input$cells, gap, input$null)
  key <- sa_mosaic_key(residual)

  if (is.null(xlab)) xlab <- input$row_var
  if (is.null(ylab)) ylab <- input$col_var

  theme <- sa_plot_theme(dark)
  if (dark) {
    tile_border <- "white"
    empty_fill <- "#36454F"
    plain_fill <- "#36454F"
  } else {
    tile_border <- "gray20"
    empty_fill <- "gray92"
    plain_fill <- "white"
  }

  # Measured before the device is split, the way `draw_heatmap()` sizes its key:
  # short band labels ask for less of the device than long ones, and the cap is
  # what stops the key from taking half of a narrow one. Where the cap binds, the
  # key shrinks its own text to the panel it was given rather than overflowing it.
  if (shade) {
    key_frac <- min(0.32, sa_mosaic_key_width(key, cex.legend) /
                      graphics::par("din")[1])
  }

  old_par <- graphics::par(c("bg", "fg", "col.axis", "col.lab", "col.main",
                             "mar", "mfrow", "oma"))
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  if (dark) {
    graphics::par(bg = theme$bg, fg = theme$fg, col.axis = theme$fg,
                  col.lab = theme$fg, col.main = theme$fg)
  }
  if (shade) {
    graphics::layout(matrix(c(1L, 2L), nrow = 1L),
                     widths = c(1 - key_frac, key_frac))
  }

  graphics::par(mar = c(5, 5, if (is.null(main)) 2 else 4, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  tiles <- layout_out$cells
  value <- if (identical(residual, "pearson")) {
    tiles$residual
  } else {
    tiles$std_residual
  }
  tiles$fill <- if (shade) {
    sa_mosaic_fill(value, key, empty_fill)
  } else {
    rep(plain_fill, nrow(tiles))
  }
  text_col <- sa_mosaic_text_col(tiles$fill)

  for (i in seq_len(nrow(tiles))) {
    tile <- tiles[i, ]
    if (!is.finite(tile$x1) || !is.finite(tile$y1) ||
          tile$x2 <= tile$x1 || tile$y2 <= tile$y1) {
      next
    }
    graphics::rect(tile$x1, tile$y1, tile$x2, tile$y2,
                   col = tile$fill, border = tile_border, lwd = 1)
    label <- sa_mosaic_anno(anno_cells, tile$observed, tile$prop_row,
                            tile$x2 - tile$x1, tile$y2 - tile$y1, cex.anno)
    if (!is.null(label)) {
      graphics::text((tile$x1 + tile$x2) / 2, (tile$y1 + tile$y2) / 2,
                     labels = label, cex = cex.anno, col = text_col[i])
    }
  }

  if (expected_line) {
    sa_mosaic_draw_expected(layout_out, theme$guide)
  }

  keep_row <- !layout_out$row_lv %in% layout_out$empty_levels$row
  keep_col <- !layout_out$col_lv %in% layout_out$empty_levels$col
  graphics::axis(1, at = layout_out$x_at[keep_row],
                 labels = layout_out$row_lv[keep_row],
                 cex.axis = cex.axis, tick = FALSE, line = -0.4)
  graphics::axis(2, at = layout_out$y_at[keep_col],
                 labels = layout_out$col_lv[keep_col],
                 cex.axis = cex.axis, tick = FALSE, las = 1, line = -0.4)
  graphics::title(xlab = xlab, ylab = ylab, main = main,
                  cex.lab = cex.lab, cex.main = cex.main)

  if (shade) {
    graphics::par(mar = c(5, 0.4, if (is.null(main)) 2 else 4, 0.6))
    sa_mosaic_draw_key(key, cex.legend, dark)
  }

  invisible(list(
    cells           = tiles,
    widths          = layout_out$widths,
    heights         = layout_out$heights,
    expected_prop   = layout_out$expected_prop,
    empty_levels    = layout_out$empty_levels,
    null            = layout_out$null,
    residual        = residual,
    residual_breaks = key$breaks,
    colors          = key$colors
  ))
}


#' @param x A categorical comparison result, as returned by
#'   [compare_categorical_groups()].
#' @param ... Arguments passed on to [draw_mosaic_plot()].
#'
#' @rdname draw_mosaic_plot
#' @export
plot.sa_categorical <- function(x, ...) {
  draw_mosaic_plot(x, ...)
}


#' The cell table of a categorical comparison, or the reason there is none
#'
#' Carries `null` along with the cells, because the cells alone do not say what
#' their `expected` column is a statement about and the shading has to know.
#'
#' @keywords internal
#' @noRd
sa_mosaic_input <- function(res) {
  if (inherits(res, "sa_categorical")) {
    return(list(
      cells   = res$cells,
      null    = res$design$null,
      row_var = res$design$row_var,
      col_var = res$design$col_var
    ))
  }

  if (inherits(res, "sa_comparison")) {
    stop("`categorical_comparison_result` is a numeric comparison result. ",
         "draw_mosaic_plot() draws a contingency table; ",
         "draw_grouped_boxplot() and draw_volcano_plot() are what read a ",
         "comparison.", call. = FALSE)
  }

  stop("`categorical_comparison_result` must be a categorical comparison ",
       "result, as returned by compare_categorical_groups(). The shading and ",
       "the expected lines are read under the null hypothesis that result was ",
       "tested against, which a bare table or a data.frame does not carry. ",
       "Cross the variables with compare_categorical_groups() first.",
       call. = FALSE)
}


#' Place the tiles of a mosaic in the unit square
#'
#' The first variable takes the x axis, one strip per level, width the marginal
#' share. Each strip is split by the second variable, height the share of that
#' strip. The same arithmetic is run a second time on the expected counts, which
#' is what puts the null hypothesis on the same scale as the tiles and lets it be
#' drawn as a line rather than described in a caption.
#'
#' A level holding no observation has no share to take and no conditional
#' distribution to be cut into, so it draws nothing and is reported instead. Its
#' gap is left in place: a space where a level should have been is the honest
#' picture of a level that was named and never seen.
#'
#' @param cells The cell table, one row per cell.
#' @param gap Fraction of the axis each single gap takes.
#' @param null Which hypothesis `cells$expected` states.
#'
#' @keywords internal
#' @noRd
sa_mosaic_layout <- function(cells, gap, null) {
  row_lv <- unique(cells$row_level)
  col_lv <- unique(cells$col_level)
  n_row <- length(row_lv)
  n_col <- length(col_lv)

  at <- cbind(match(cells$row_level, row_lv), match(cells$col_level, col_lv))
  observed <- matrix(0, n_row, n_col, dimnames = list(row_lv, col_lv))
  expected <- observed
  observed[at] <- cells$observed
  expected[at] <- ifelse(is.finite(cells$expected), cells$expected, 0)

  row_n <- rowSums(observed)
  col_n <- colSums(observed)
  total <- sum(observed)

  widths <- if (total == 0) rep(1 / n_row, n_row) else row_n / total
  heights <- observed / pmax(row_n, 1)
  heights[row_n == 0, ] <- 0

  # Normalised within the strip rather than by `row_n`, because the row sums of
  # the expected table are the row sums of the observed one only under
  # independence. Under symmetry they are the average of the row and column
  # margins, and what is being compared is still one distribution to another.
  exp_row <- rowSums(expected)
  expected_prop <- expected / pmax(exp_row, 1)
  expected_prop[exp_row == 0, ] <- 0

  # `gap` is one gap rather than all of them, so a five-level axis is not five
  # times finer than a two-level one. The cap keeps the tiles from vanishing when
  # a wide gap meets many levels.
  gap_x <- if (n_row > 1L) min(gap, 0.4 / (n_row - 1L)) else 0
  gap_y <- if (n_col > 1L) min(gap, 0.4 / (n_col - 1L)) else 0
  usable_x <- 1 - gap_x * (n_row - 1L)
  usable_y <- 1 - gap_y * (n_col - 1L)

  cum_w <- c(0, cumsum(widths))
  x0 <- cum_w[seq_len(n_row)] * usable_x + (seq_len(n_row) - 1L) * gap_x
  x1 <- x0 + widths * usable_x

  tiles <- cells
  tiles$x1 <- NA_real_
  tiles$x2 <- NA_real_
  tiles$y1 <- NA_real_
  tiles$y2 <- NA_real_

  # One boundary fewer than there are tiles: the top of the last one is the top
  # of the strip under every hypothesis, so there is nothing to mark there.
  expected_y <- matrix(NA_real_, nrow = n_row, ncol = max(n_col - 1L, 0L),
                       dimnames = list(row_lv, col_lv[seq_len(n_col - 1L)]))

  # The reference strip is the first, which the comparison settled through its
  # `control_label` and `category_lv`, unless it is empty and so has no
  # boundaries to label.
  ref <- which(row_n > 0)[1]
  if (is.na(ref)) ref <- 1L
  y_at <- stats::setNames(numeric(n_col), col_lv)

  for (i in seq_len(n_row)) {
    cum_h <- c(0, cumsum(heights[i, ]))
    y0 <- cum_h[seq_len(n_col)] * usable_y + (seq_len(n_col) - 1L) * gap_y
    y1 <- y0 + heights[i, ] * usable_y
    if (i == ref) {
      y_at[] <- (y0 + y1) / 2
    }
    if (n_col > 1L) {
      cum_e <- cumsum(expected_prop[i, ])[seq_len(n_col - 1L)]
      expected_y[i, ] <- cum_e * usable_y + (seq_len(n_col - 1L) - 1L) * gap_y
    }
    for (j in seq_len(n_col)) {
      hit <- at[, 1] == i & at[, 2] == j
      tiles$x1[hit] <- x0[i]
      tiles$x2[hit] <- x1[i]
      tiles$y1[hit] <- y0[j]
      tiles$y2[hit] <- y1[j]
    }
  }

  list(
    cells         = tiles,
    widths        = stats::setNames(as.numeric(widths), row_lv),
    heights       = heights,
    expected_prop = expected_prop,
    expected_y    = expected_y,
    strip_x       = cbind(x0 = x0, x1 = x1),
    empty_levels  = list(row = row_lv[row_n == 0], col = col_lv[col_n == 0]),
    null          = null,
    row_lv        = row_lv,
    col_lv        = col_lv,
    x_at          = (x0 + x1) / 2,
    y_at          = y_at
  )
}


#' Mark each strip where the null hypothesis would have cut it
#'
#' @keywords internal
#' @noRd
sa_mosaic_draw_expected <- function(layout_out, col) {
  expected_y <- layout_out$expected_y
  if (ncol(expected_y) == 0L) {
    return(invisible(NULL))
  }
  strip <- layout_out$strip_x
  empty <- layout_out$row_lv %in% layout_out$empty_levels$row

  for (i in seq_len(nrow(expected_y))) {
    if (empty[i] || strip[i, "x1"] <= strip[i, "x0"]) {
      next
    }
    ok <- is.finite(expected_y[i, ])
    if (!any(ok)) {
      next
    }
    graphics::segments(strip[i, "x0"], expected_y[i, ok],
                       strip[i, "x1"], expected_y[i, ok],
                       col = col, lty = 3, lwd = 1.6)
  }
  invisible(NULL)
}


#' The residual scale a mosaic is shaded on
#'
#' The cuts of 2 and 4 are the ones a residual referred to a standard normal
#' would call surprising and extreme. Pearson residuals grow with the sample, so
#' a large table lights up more tiles than a small one of the same shape; that
#' is the same fact the chi-square statistic itself reports, and it is why the key
#' names which residual it is a scale for rather than only the cuts.
#'
#' @keywords internal
#' @noRd
sa_mosaic_key <- function(residual = "pearson") {
  list(
    residual = residual,
    title    = if (identical(residual, "pearson")) {
      "Pearson residual"
    } else {
      "Std. residual"
    },
    breaks   = c(-Inf, -4, -2, 2, 4, Inf),
    labels   = c("< -4", "-4 to -2", "-2 to 2", "2 to 4", "> 4"),
    colors   = c("#4575B4", "#A6C5DE", "gray88", "#F4A6A0", "#D73027")
  )
}


#' Map residuals onto the mosaic palette
#'
#' @keywords internal
#' @noRd
sa_mosaic_fill <- function(residual, key, empty_fill) {
  out <- rep(empty_fill, length(residual))
  ok <- is.finite(residual)
  if (!any(ok)) {
    return(out)
  }
  bin <- findInterval(residual[ok], key$breaks, rightmost.closed = TRUE,
                      all.inside = TRUE)
  out[ok] <- key$colors[bin]
  out
}


#' Ink that can be read on a given fill
#'
#' Decided from the fill's luminance rather than from the residual that chose it,
#' so an unshaded tile, an empty one and a dark background are all covered by the
#' same rule.
#'
#' @keywords internal
#' @noRd
sa_mosaic_text_col <- function(fill) {
  rgb <- grDevices::col2rgb(fill)
  luminance <- (0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]) / 255
  ifelse(luminance < 0.5, "white", "grey15")
}


#' As much of a tile's numbers as the tile has room for
#'
#' Measured against the label with [graphics::strwidth()] rather than against a
#' fixed fraction of the plot, because a fraction that fits `"7"` does not fit
#' `"1284"` and a mosaic of a large table holds both.
#'
#' @param mode What was asked for. `"auto"` gives up one line at a time.
#' @param observed,prop_row The count and its share of the strip.
#' @param w,h The tile, in user units, which are fractions of the unit square.
#' @param cex Character expansion the label would be drawn at.
#'
#' @return The label, or `NULL` when nothing fits or nothing was asked for.
#'
#' @keywords internal
#' @noRd
sa_mosaic_anno <- function(mode, observed, prop_row, w, h, cex) {
  if (identical(mode, "none")) {
    return(NULL)
  }

  count <- format(observed, trim = TRUE)
  pct <- if (is.finite(prop_row)) {
    paste0(format(round(100 * prop_row), trim = TRUE), "%")
  } else {
    NULL
  }

  wanted <- switch(
    mode,
    count   = count,
    percent = pct,
    both    = if (is.null(pct)) count else paste0(count, "\n", pct),
    auto    = c(if (is.null(pct)) NULL else paste0(count, "\n", pct), count)
  )
  if (length(wanted) == 0L) {
    return(NULL)
  }

  for (label in wanted) {
    if (graphics::strwidth(label, cex = cex) <= 0.92 * w &&
          graphics::strheight(label, cex = cex) <= 0.92 * h) {
      return(label)
    }
  }
  # An explicit request is honoured even where it overflows, since the caller
  # asked for the number rather than for a tidy picture. `"auto"` did not.
  if (identical(mode, "auto")) NULL else wanted[1]
}


#' The fraction of the panel the key leaves free at each side
#'
#' Shared by the two functions below, so the width one of them reserves is the
#' width the other one draws into.
#'
#' @keywords internal
#' @noRd
sa_mosaic_key_pad <- function() {
  c(left = 0.06, right = 0.02)
}


#' How wide a strip the residual key needs
#'
#' Modelled on `sa_key_width()` in `draw_heatmap.R`: the width comes from what
#' goes in the key rather than from a fraction of the device, so a key of short
#' band labels does not reserve room for long ones.
#'
#' What it adds is the panel's own overheads. `layout()` hands out a share of the
#' device, and the plot region inside that share is smaller by the margin and by
#' the padding the key draws with, so a width measured on the text alone reserves
#' a panel the text does not fit in.
#'
#' @keywords internal
#' @noRd
sa_mosaic_key_width <- function(key, cex.legend) {
  char <- graphics::par("cin")[1]
  wid <- function(s, cex, font = 1) {
    max(graphics::strwidth(s, units = "inches", cex = cex, font = font))
  }
  swatch <- 2.2 * char * cex.legend
  bands <- swatch + 0.5 * char * cex.legend + wid(key$labels, cex.legend)
  content <- max(bands, wid(key$title, cex.legend, font = 2))

  pad <- sa_mosaic_key_pad()
  margin <- (0.4 + 0.6) * graphics::par("cin")[2]
  margin + content / (1 - pad[["left"]] - pad[["right"]])
}


#' Draw the residual key in the panel `layout()` already opened
#'
#' Highest band at the top, so the key runs the same way up as the colours do on
#' the tiles: more than expected above, less than expected below.
#'
#' The character expansion is fitted here rather than taken as given. Every width
#' in the key is proportional to it, so the largest one that fits is a division
#' rather than a search, and doing it in the panel is what makes the fit exact:
#' `sa_mosaic_key_width()` asked `layout()` for a share of the device, and this is
#' the panel that share turned out to be.
#'
#' @keywords internal
#' @noRd
sa_mosaic_draw_key <- function(key, cex.legend, dark) {
  n <- length(key$colors)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  lab_col <- if (dark) "white" else "grey15"
  border <- if (dark) "white" else "gray30"

  pad <- sa_mosaic_key_pad()
  unit_char <- graphics::par("cxy")[1]
  unit_wanted <- max(
    2 * unit_char + 0.4 * unit_char + max(graphics::strwidth(key$labels)),
    graphics::strwidth(key$title, font = 2)
  )
  cex <- min(cex.legend, (1 - pad[["left"]] - pad[["right"]]) / unit_wanted)

  char_w <- unit_char * cex
  line_h <- graphics::par("cxy")[2] * cex
  left <- pad[["left"]]
  swatch <- min(0.34, 2 * char_w)
  top <- min(0.96, 1 - 0.4 * line_h)

  graphics::text(left, top, labels = key$title, adj = c(0, 1),
                 cex = cex, font = 2, col = lab_col)

  band_top <- top - 1.7 * line_h
  band_h <- max(min(1.25 * line_h, (band_top - 0.04) / n), 0.02)
  for (i in seq_len(n)) {
    idx <- n - i + 1L
    y1 <- band_top - (i - 1L) * band_h
    y0 <- y1 - band_h
    graphics::rect(left, y0, left + swatch, y1, col = key$colors[idx],
                   border = border)
    graphics::text(left + swatch + 0.4 * char_w, (y0 + y1) / 2,
                   labels = key$labels[idx], adj = c(0, 0.5),
                   cex = cex, col = lab_col)
  }
  invisible(NULL)
}
