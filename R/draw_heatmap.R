#' Draw a clustered heatmap of features by samples
#'
#' Draws one cell per feature and sample, with the sample groups as a coloured
#' annotation strip above the columns and a dendrogram on each axis that was
#' clustered. The input is the wide format the comparison functions take, one row
#' per observation and one column per feature, and it is transposed so that
#' features run down the rows the way an expression heatmap is usually read.
#' [stats::heatmap()] draws the cells, the strip and the trees; this function
#' decides what goes into them and adds the colour key beside them.
#'
#' @param data A data.frame or a matrix in wide format, one row per observation
#'   and one column per feature. This is the same layout
#'   [compare_two_groups()] and [compare_multiple_groups()] take. Either type is
#'   accepted and converted to the numeric matrix the drawing needs; a matrix
#'   with no column names has them made up as `V1`, `V2` and so on, and repeated
#'   row names are kept, since a sample name is a naming choice rather than a
#'   key. There must be at least two features and two samples to draw.
#' @param group Grouping vector with one entry per row of `data`, or `NULL` when
#'   `group_lv` is also `NULL` to draw without a group strip or group legend.
#' @param group_lv Character vector of group levels, in the order they should
#'   appear in the legend, or `NULL` with `group = NULL`. Rows belonging to any
#'   other level are dropped when both are supplied.
#' @param feats Character vector of numeric column names to draw, or `NULL` for
#'   every column of `data`. The order given is the order the rows are drawn in
#'   when they are not clustered.
#' @param scale How to put the features on a comparable scale before drawing:
#'   `"feature"` z-scores each feature across the samples, `"sample"` z-scores
#'   each sample across the features, and `"none"` draws the values as they
#'   arrived. The scaled values are what `matrix` comes back holding.
#' @param zlim Numeric length-2 range the colours span, or `NULL` to derive it
#'   from the values being drawn. Values outside a supplied range are drawn at
#'   the end of the scale rather than left blank, and how many were is reported.
#' @param dist_method Distance behind the clustering. `"correlation"` is
#'   `1 - cor()`, which groups features by the shape of their profile rather
#'   than by how high it sits; the other two are handed to [stats::dist()].
#' @param hclust_method Linkage handed to [stats::hclust()].
#' @param cluster_feats,cluster_samples Whether to cluster and reorder that
#'   axis. `FALSE` keeps the input order.
#' @param feat_labels,sample_labels Labels to draw in place of the column names
#'   of `data` and its row names. `sample_labels` has one entry per row of
#'   `data` and is filtered along with it; `feat_labels` has one entry per
#'   feature in `feats`.
#' @param show_feat_names,show_sample_names Whether to draw those labels. The
#'   labels still exist when they are not drawn, so a hidden axis does not
#'   change anything else about the plot.
#' @param anno If `TRUE`, write each cell value on the cell, rounded to two
#'   decimal places. The numbers are the values in the returned `matrix` (after
#'   `scale`, and not clamped to `zlim`). Missing cells are left blank.
#' @param cex.anno Character expansion for those cell labels, relative to the
#'   smaller of the feature and sample axis label sizes.
#' @param n_colors Number of colours in the blue-white-red ramp.
#' @param main Plot title.
#' @param cex.axis,cex.main,cex.legend Character expansion for the axis labels,
#'   the title and the group legend.
#'
#' @return A list, invisibly:
#'
#'   \describe{
#'     \item{`matrix`}{The scaled matrix as it was drawn: features in rows,
#'       samples in columns, both in the order they appear on the plot, and the
#'       row and column names as they were labelled. Values are not clamped to
#'       `zlim`, which is a decision about colour rather than about the data.}
#'     \item{`feat_order`,`sample_order`}{The permutations the clustering chose,
#'       as indices into `feats` and into the rows of `data` that were kept.}
#'     \item{`feat_hclust`,`sample_hclust`}{The [stats::hclust()] objects behind
#'       the dendrograms, or `NULL` for an axis that was not clustered.}
#'     \item{`zlim`}{The colour range, derived or as supplied.}
#'     \item{`group_colors`}{The colour drawn for each group level, or `NULL`
#'       when no group was supplied.}
#'   }
#'
#' @details
#' Features are z-scored across the samples by default. Without it a single
#' high-abundance feature takes the whole colour range and everything else is
#' left white, since the colour scale is shared by every cell in the plot and
#' features are not measured on a common scale. The plot does not say which
#' `scale` ran: naming it would take a line of text wider than the key it titles,
#' and the numbers beside the key are the scaled values either way. A
#' feature with no variance would divide by zero, so it is only centred and ends
#' up flat at the middle of the scale; how many were is reported in a message.
#'
#' A diverging palette needs a meaningful midpoint. When `zlim` is not given,
#' zero is that midpoint if the values being drawn have both signs, which is
#' always the case after z-scoring, and the range is made symmetric around it.
#' Values of one sign only have no such point, so the range is the range of the
#' values and white falls in the middle of it.
#'
#' Missing values are drawn as grey cells rather than dropped, so a gap is
#' visible as a gap. Clustering leaves them to [stats::dist()], which measures a
#' pair on the values they share. A pair that shares none has no distance at
#' all, and rather than fail, that axis is left in its input order with a
#' message saying so.
#'
#' The clustering is done here rather than left to [stats::heatmap()], which
#' would standardise after clustering and divide a flat feature by zero, and the
#' trees are handed over already built. What comes back is therefore the tree the
#' plot shows. The panel proportions are [stats::heatmap()]'s, including the
#' square cells its `respect = TRUE` layout asks for. That layout keeps the block
#' of panels square and centres it, so a device much wider than it is tall is left
#' with empty space to the right of the plot and its key.
#'
#' [stats::heatmap()] draws no colour key, so this function keeps a strip beside
#' the plot for one and overlays it after: an upright bar with its numbers to the
#' right of it, and the group legend to the right of those, both starting level
#' with the top row of cells. The strip is as wide as those need rather than a
#' fixed fraction of the device, and it is placed against the feature labels
#' rather than at the edge, since the space the centred layout leaves over would
#' otherwise read as a gap between the plot and its key. The panel layout, the
#' margins and the plot region are all restored on exit, so the caller's device
#' is left as it was found.
#'
#' @seealso [draw_grouped_boxplot()], which takes the same wide input.
#'
#' @examples
#' ## Features are z-scored, so the two planted directions separate by colour
#' ## and the clustering of the columns finds the two groups on its own.
#' sim <- simulate_two_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 3)
#' out <- draw_heatmap(
#'   data     = as.matrix(sim$args$data),
#'   group    = sim$args$group,
#'   group_lv = sim$args$group_lv,
#'   show_sample_names = FALSE,
#'   main     = "30 features, 10 of them planted"
#' )
#'
#' ## The clustering is on the result, so what the picture shows can be checked
#' ## rather than eyeballed.
#' head(rownames(out$matrix))
#' out$sample_hclust$method
#'
#' ## Only the planted features, in the order they were named
#' draw_heatmap(
#'   data          = sim$args$data,
#'   group         = sim$args$group,
#'   group_lv      = sim$args$group_lv,
#'   feats         = sim$truth$features[sim$truth$direction != "none"],
#'   cluster_feats = FALSE,
#'   show_sample_names = FALSE
#' )
#'
#' ## Unscaled values, with the four measurements on their own scale
#' draw_heatmap(
#'   data     = iris[1:4],
#'   group    = iris$Species,
#'   group_lv = levels(iris$Species),
#'   scale    = "none",
#'   show_sample_names = FALSE,
#'   main     = "iris, cm"
#' )
#'
#' ## No group strip or group legend when both are omitted
#' draw_heatmap(
#'   data = as.matrix(sim$args$data),
#'   show_sample_names = FALSE
#' )
#'
#' @export
draw_heatmap <- function(data,
                         group = NULL,
                         group_lv = NULL,
                         feats = NULL,
                         scale = c("feature", "sample", "none"),
                         zlim = NULL,
                         dist_method = c("euclidean", "correlation",
                                         "manhattan"),
                         hclust_method = c("average", "complete", "ward.D2"),
                         cluster_feats = TRUE,
                         cluster_samples = TRUE,
                         feat_labels = NULL,
                         sample_labels = NULL,
                         show_feat_names = TRUE,
                         show_sample_names = TRUE,
                         anno = FALSE,
                         cex.anno = 1,
                         n_colors = 101,
                         main = NULL,
                         cex.axis = 0.9,
                         cex.main = 1.5,
                         cex.legend = 1.2) {

  scale <- match.arg(scale)
  dist_method <- match.arg(dist_method)
  hclust_method <- match.arg(hclust_method)
  sa_check_flag(cluster_feats, "cluster_feats")
  sa_check_flag(cluster_samples, "cluster_samples")
  sa_check_flag(show_feat_names, "show_feat_names")
  sa_check_flag(show_sample_names, "show_sample_names")
  sa_check_flag(anno, "anno")
  sa_check_scalar_num(cex.anno, "cex.anno", 0, lower_open = TRUE)
  sa_check_count(n_colors, "n_colors", 3)
  sa_check_scalar_num(cex.axis, "cex.axis", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.main, "cex.main", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.legend, "cex.legend", 0, lower_open = TRUE)
  if (!is.null(zlim)) {
    sa_check_range(zlim, "zlim")
    if (zlim[1] == zlim[2]) {
      stop("`zlim` must have two different ends, but both are ", zlim[1], ".",
           call. = FALSE)
    }
  }
  if (xor(is.null(group), is.null(group_lv))) {
    stop("`group` and `group_lv` must both be supplied or both be `NULL`.",
         call. = FALSE)
  }
  annotate_group <- !is.null(group)

  # A matrix carries its labels in its dimnames, and they are needed after the
  # rows have been filtered. The column names go in as `feats` and the row names
  # ride along as `id`, which is the one argument the validator filters beside
  # the data. Nothing else here needs a pairing key.
  if (is.matrix(data) && is.null(colnames(data))) {
    colnames(data) <- paste0("V", seq_len(ncol(data)))
  }
  if (is.null(feats)) {
    feats <- colnames(data)
  }
  if (!is.null(sample_labels) && length(sample_labels) != NROW(data)) {
    stop("`sample_labels` must have one entry per row of `data`: got ",
         length(sample_labels), " for ", NROW(data), " rows.", call. = FALSE)
  }
  if (is.null(sample_labels)) {
    sample_labels <- rownames(data)
    if (is.null(sample_labels)) {
      sample_labels <- as.character(seq_len(NROW(data)))
    }
  }
  if (is.matrix(data)) {
    # as.data.frame() inside the validator rejects repeated row names, which are
    # a sample naming choice rather than an error. The names are already saved.
    rownames(data) <- NULL
  }

  input <- sa_validate_wide_input(data, feats, group, group_lv,
                                  id = as.character(sample_labels),
                                  min_levels = 1L)
  feats <- input$feats
  group <- input$group
  sample_labels <- input$id

  if (!is.null(feat_labels) && length(feat_labels) != length(feats)) {
    stop("`feat_labels` must have one entry per feature in `feats`: got ",
         length(feat_labels), " for ", length(feats), " feature(s).",
         call. = FALSE)
  }
  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `group_lv`.")
  }

  # Features down the rows and samples across the columns, which is the way an
  # expression heatmap is read and the orientation the group strip needs to sit
  # above the columns it describes.
  m <- t(as.matrix(input$data[feats]))
  rownames(m) <- if (is.null(feat_labels)) feats else as.character(feat_labels)
  colnames(m) <- sample_labels
  n_feats <- nrow(m)
  n_samples <- ncol(m)
  if (n_feats < 2L || n_samples < 2L) {
    stop("`draw_heatmap()` needs at least 2 features and 2 samples to cluster ",
         "and draw, but got ", n_feats, " feature(s) and ", n_samples,
         " sample(s).", call. = FALSE)
  }

  m <- sa_scale_matrix(m, scale)

  feat_hc <- if (cluster_feats) {
    sa_heatmap_hclust(m, dist_method, hclust_method, "feature")
  }
  sample_hc <- if (cluster_samples) {
    sa_heatmap_hclust(t(m), dist_method, hclust_method, "sample")
  }
  feat_order <- if (is.null(feat_hc)) seq_len(n_feats) else feat_hc$order
  sample_order <- if (is.null(sample_hc)) seq_len(n_samples) else sample_hc$order

  m_plot <- m[feat_order, sample_order, drop = FALSE]

  finite <- m[is.finite(m)]
  if (length(finite) == 0L) {
    stop("`data` holds no finite value to draw.", call. = FALSE)
  }
  if (is.null(zlim)) {
    zlim <- if (any(finite > 0) && any(finite < 0)) {
      # Zero is the midpoint the blue-white-red ramp is built around, so the
      # range is made symmetric about it rather than about the data.
      c(-1, 1) * max(abs(finite))
    } else {
      range(finite)
    }
    if (zlim[1] == zlim[2]) {
      zlim <- zlim + c(-0.5, 0.5)
    }
  }

  cols <- grDevices::colorRampPalette(c("blue", "white", "red"))(n_colors)
  breaks <- seq(zlim[1], zlim[2], length.out = n_colors + 1L)

  # image() drops a value outside its breaks, which would punch holes in the
  # plot exactly where the largest effects are. They are drawn at the end of the
  # scale instead, and the count is reported since a clamped cell no longer
  # carries its own value. The clamping is on the input order, since
  # stats::heatmap() reorders what it is handed by the dendrograms it is handed;
  # `m_plot` is the same values in the order they come out.
  z <- m
  clamped <- sum(z < zlim[1] | z > zlim[2], na.rm = TRUE)
  z[!is.na(z) & z < zlim[1]] <- zlim[1]
  z[!is.na(z) & z > zlim[2]] <- zlim[2]
  if (clamped > 0L) {
    message(clamped, " of ", length(z),
            " cell(s) lie outside `zlim` and are drawn at the end of the ",
            "colour scale.")
  }

  group_cols <- if (annotate_group) {
    cols_g <- grDevices::hcl.colors(nlevels(group), "Dark 2")
    stats::setNames(cols_g, levels(group))
  } else {
    NULL
  }

  # The strip on the right that the colour key gets, measured from what goes in
  # it rather than guessed at as a fraction of the device: short group names ask
  # for less of it than long ones. The cap is what stops the key from taking half
  # of a narrow device, and `key_cex` shrinks the text that would then not fit.
  key_wanted <- sa_key_width(zlim, group_cols, cex.axis, cex.legend)
  key_frac <- min(0.4, key_wanted / graphics::par("din")[1])
  key_cex <- min(1, key_frac * graphics::par("din")[1] / key_wanted)

  # Only the parameters this function sets are put back, not a blanket
  # par(no.readonly = TRUE) snapshot. That snapshot also carries `fin`, `pin`
  # and `mai`, which are absolute sizes: restoring them pins the figure to the
  # size this plot happened to be drawn at, so the next plot on a device that
  # has since been resized is redrawn small in a corner of it. `mfrow` comes
  # back with the rest, since `layout()` overwrites whatever grid the caller
  # had set up, and so does `new`, which the key panel sets and nothing clears.
  old_par <- graphics::par(c("mar", "mfrow", "oma", "omd", "cex.main", "new"))
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  # A label earns the margin it needs, up to a third of the figure, and is
  # shrunk to fit past that. Cutting the panel down instead would leave the
  # cells in a strip too narrow to read. The margin is the width the labels
  # actually draw at, since a margin guessed from the number of characters is
  # half again too wide for lowercase text, and every line of it that is not
  # needed reads as a gap between the cells and the key.
  csi <- graphics::par("csi")
  lines_wide <- graphics::par("din")[1] * (1 - key_frac) / csi
  lines_high <- graphics::par("din")[2] / csi
  feat_wanted <- sa_text_lines(rownames(m_plot), cex.axis)
  sample_wanted <- sa_text_lines(colnames(m_plot), cex.axis)
  mar_right <- if (show_feat_names) min(feat_wanted, lines_wide / 3) else 0.4
  mar_bottom <- if (show_sample_names) {
    min(sample_wanted, lines_high / 3)
  } else {
    0.4
  }
  cex_feat <- cex.axis * min(1, mar_right / feat_wanted)
  cex_sample <- cex.axis * min(1, mar_bottom / sample_wanted)
  cex_anno <- min(cex_feat, cex_sample) * cex.anno

  # heatmap() lays its own panels out with layout(), which divides the region
  # inside the outer margins. Narrowing that region is what leaves the strip the
  # colour key is drawn in, since heatmap() draws no key of its own. `omd` and
  # `oma` overwrite one another, so the title is left to heatmap() rather than
  # given an outer margin of its own.
  graphics::par(omd = c(0, 1 - key_frac, 0, 1),
                cex.main = cex.main / 1.5)  # heatmap() draws it at 1.5x

  # Filled in by `add.expr` below, and NULL if that never ran.
  cells_edge <- NULL

  hm_args <- list(
    x = z,
    # A dendrogram goes in already built, so the tree the plot shows is the one
    # this function clustered and reported. NA is the axis that was not
    # clustered: heatmap() then keeps the input order and narrows that panel to
    # nothing rather than needing a layout of its own.
    Rowv = if (is.null(feat_hc)) NA else stats::as.dendrogram(feat_hc),
    Colv = if (is.null(sample_hc)) NA else stats::as.dendrogram(sample_hc),
    # The scaling already ran, above, where a feature with no variance is
    # centred instead of divided by zero. heatmap() would divide.
    scale = "none",
    # First feature at the top, which is where this function has always drawn it
    # and the way an expression heatmap is read.
    revC = TRUE,
    labRow = if (show_feat_names) rownames(z) else rep("", n_feats),
    labCol = if (show_sample_names) colnames(z) else rep("", n_samples),
    cexRow = cex_feat,
    cexCol = cex_sample,
    margins = c(mar_bottom, mar_right),
    # Handed through to image(), which is what makes the colours mean what the
    # key says they mean.
    col = cols,
    breaks = breaks,
    main = main
  )
  hm_args$add.expr <- substitute({
    sa_mark_missing(m_plot)
    if (anno) {
      sa_annotate_heatmap_cells(m_plot, zlim, cex_anno)
    }
    cells_edge <- c(
      right = graphics::grconvertX(graphics::par("usr")[2], "user", "ndc"),
      top = graphics::grconvertY(graphics::par("usr")[4], "user", "ndc")
    )
  }, list(m_plot = m_plot, anno = anno, zlim = zlim, cex_anno = cex_anno))
  if (annotate_group) {
    hm_args$ColSideColors <- unname(group_cols[as.integer(group)])
  }
  do.call(stats::heatmap, hm_args)

  # heatmap() puts the caller's par back on exit, the narrowed region with it,
  # so the key asks for the whole device again and overlays what is beside the
  # cells. It sits against the labels rather than at the edge of the device:
  # everything the centred layout left over would otherwise read as a gap
  # between the plot and its key. The top of the cells is where the key starts,
  # for the same reason.
  key_left <- 1 - key_frac
  key_top <- 1
  if (!is.null(cells_edge)) {
    key_left <- min(key_left, cells_edge[["right"]] +
                      (mar_right + 0.6) * csi / graphics::par("din")[1])
    key_top <- cells_edge[["top"]]
  }
  graphics::par(omd = c(0, 1, 0, 1))
  # No margins, so that the panel's own 0..1 is the device's own 0..1 and the
  # measured top needs no converting.
  graphics::par(fig = c(key_left, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  sa_draw_heatmap_key(cols, zlim, group_cols, cex.axis * key_cex,
                      cex.legend * key_cex, top = key_top)

  invisible(list(
    matrix        = m_plot,
    feat_order    = feat_order,
    sample_order  = sample_order,
    feat_hclust   = feat_hc,
    sample_hclust = sample_hc,
    zlim          = zlim,
    group_colors  = group_cols
  ))
}


#' Centre and standardise a heatmap matrix along one margin
#'
#' A feature with no spread has nothing to standardise by, so it is centred and
#' left flat rather than turned into `NaN` by a division by zero. The same
#' applies to a margin that is entirely missing, whose mean and standard
#' deviation are both `NA` to begin with.
#'
#' @param x Features in rows, samples in columns.
#'
#' @keywords internal
#' @noRd
sa_scale_matrix <- function(x, scale) {
  if (scale == "none") {
    return(x)
  }
  margin <- if (scale == "feature") 1L else 2L
  centre <- apply(x, margin, mean, na.rm = TRUE)
  spread <- apply(x, margin, stats::sd, na.rm = TRUE)
  flat <- !is.finite(spread) | spread == 0
  spread[flat] <- 1
  x <- sweep(x, margin, centre, "-")
  x <- sweep(x, margin, spread, "/")
  n_flat <- sum(flat & is.finite(centre))
  if (n_flat > 0L) {
    message(n_flat, " ", scale, "(s) have no variance to scale by and are ",
            "drawn flat at the middle of the colour scale.")
  }
  x
}


#' Draw the colour key and the group legend beside the cells
#'
#' The two stand side by side, the bar with its numbers on the left and the group
#' legend to the right of them, both starting at the top of the cells. The bar
#' carries no title: what the colours measure is `scale`, which is an argument
#' rather than something the plot can be asked, and a title long enough to say it
#' sets the width of the whole strip.
#'
#' The panel is as tall as the device, so the bar is sized in units of text
#' rather than as a fraction of it: a fraction would grow the bar into a second
#' heatmap on a tall device. It comes out about an inch tall and a fifth of one
#' wide at the default font size. [graphics::par()]`("cxy")` reads a character
#' width and height in the coordinates of the panel it is called in, which is
#' what makes the arithmetic below independent of the size of the device.
#'
#' @param cols The colour ramp, in order, which is the bar itself.
#' @param zlim The range the ramp spans, labelling the bar.
#' @param group_cols One colour per group level, named.
#' @param top Where the top of the cells is, in the coordinates of this panel,
#'   which are the device's own.
#'
#' @keywords internal
#' @noRd
sa_draw_heatmap_key <- function(cols, zlim, group_cols, cex.axis, cex.legend,
                                top = 1) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  char_w <- graphics::par("cxy")[1] * cex.axis
  line_h <- graphics::par("cxy")[2] * cex.axis
  # Five lines tall and a character and a half wide, which is a bar of about an
  # inch by a fifth of one at the default font size: enough to read the ramp off
  # and no more. The height is capped for a device short enough that five lines
  # would be most of it.
  bar_left <- 0.35 * char_w
  bar_right <- bar_left + 1.5 * char_w
  bar_top <- min(top, 1 - 0.3 * line_h)
  bar_bottom <- max(0.05, bar_top - min(0.42, 5 * line_h))

  # The bar is the ramp in the order it was built, so the colours need no breaks
  # to be assigned by: cell i is colour i. It runs upwards, which puts the low
  # end of `zlim` at the bottom.
  edges <- seq(bar_bottom, bar_top, length.out = length(cols) + 1L)
  graphics::rect(bar_left, edges[-length(edges)], bar_right, edges[-1],
                 col = cols, border = NA)

  # The numbers alone, with no tick marks: the bar is the axis.
  at <- sa_key_ticks(zlim)
  at_lab <- format(at, trim = TRUE)
  at_y <- bar_bottom + (bar_top - bar_bottom) * (at - zlim[1]) / diff(zlim)
  graphics::text(bar_right + 0.4 * char_w, at_y, at_lab, adj = c(0, 0.5),
                 cex = cex.axis)

  if (length(group_cols) > 0L) {
    # Placed by coordinate rather than by keyword, since the panel is the height of
    # the device and every keyword would put the legend somewhere in the middle of
    # it. The title is a line of bold text of its own rather than `title = `: the
    # `title.font` that would embolden that one arrived in R 4.2.0 and this package
    # supports 4.1.
    group_x <- bar_right + 0.4 * char_w +
      max(graphics::strwidth(at_lab, cex = cex.axis)) +
      1.2 * graphics::par("cxy")[1] * cex.legend
    # The last word on how wide the legend is belongs to legend() itself, which is
    # asked before anything is drawn: the strip was reserved from an estimate, and
    # a level name long enough to beat it is shrunk rather than drawn off the edge.
    box_w <- graphics::legend(group_x, bar_top, legend = names(group_cols),
                              fill = group_cols, bty = "n", cex = cex.legend,
                              plot = FALSE)$rect$w
    cex.legend <- cex.legend * min(1, (1 - group_x) / box_w)
    line_legend <- graphics::par("cxy")[2] * cex.legend
    graphics::text(group_x + 0.5 * graphics::par("cxy")[1] * cex.legend, bar_top,
                   "group", adj = c(0, 1), cex = cex.legend, font = 2)
    graphics::legend(group_x, bar_top - 1.2 * line_legend,
                     legend = names(group_cols), fill = group_cols,
                     border = "white", bty = "n", cex = cex.legend,
                     xjust = 0, yjust = 1)
  }
  invisible(NULL)
}


#' Paint the cells that hold no value
#'
#' A gap reads as a gap only if something is drawn where the value would have
#' been: [graphics::image()] leaves a missing cell as bare device background.
#' This runs as the `add.expr` of [stats::heatmap()], which evaluates it in the
#' caller's frame once the cells are drawn, so it works in their coordinates.
#'
#' @param x The matrix in the order it was drawn, features in rows.
#'
#' @keywords internal
#' @noRd
sa_mark_missing <- function(x) {
  at <- which(is.na(x), arr.ind = TRUE)
  if (nrow(at) == 0L) {
    return(invisible(NULL))
  }
  # `revC = TRUE` drew the first row at the top, so row i is at n - i + 1.
  y <- nrow(x) - at[, "row"] + 1L
  graphics::rect(at[, "col"] - 0.5, y - 0.5, at[, "col"] + 0.5, y + 0.5,
                 col = "gray92", border = NA)
  invisible(NULL)
}


#' Write rounded values on heatmap cells
#'
#' Runs inside [stats::heatmap()] `add.expr`, in the cell panel's coordinates.
#' Row `i` of `x` is drawn at vertical position `nrow(x) - i + 1` when
#' `revC = TRUE`.
#'
#' @param x Matrix in draw order, features in rows.
#' @param zlim Colour range, used only to pick a readable text colour.
#' @param cex Character expansion for the labels.
#'
#' @keywords internal
#' @noRd
sa_annotate_heatmap_cells <- function(x, zlim, cex) {
  at <- which(is.finite(x), arr.ind = TRUE)
  if (nrow(at) == 0L) {
    return(invisible(NULL))
  }
  nr <- nrow(x)
  span <- diff(zlim)
  for (k in seq_len(nrow(at))) {
    i <- at[k, "row"]
    j <- at[k, "col"]
    v <- x[i, j]
    y <- nr - i + 1L
    lab <- format(round(v, 2), trim = TRUE, scientific = FALSE)
    t <- if (span > 0) (v - zlim[1]) / span else 0.5
    # Mid-range values sit on the pale part of the ramp; the ends are saturated.
    col <- if (t > 0.28 && t < 0.72) "grey15" else "white"
    graphics::text(j, y, labels = lab, cex = cex, col = col)
  }
  invisible(NULL)
}


#' Cluster the rows of a heatmap matrix
#'
#' @param x Objects to cluster in rows.
#' @param axis Named in the message when the clustering cannot run.
#'
#' @return An [stats::hclust()] object, or `NULL` when the distances are not all
#'   defined, in which case the caller keeps the input order.
#'
#' @keywords internal
#' @noRd
sa_heatmap_hclust <- function(x, dist_method, hclust_method, axis) {
  d <- if (dist_method == "correlation") {
    # Rows are the objects, and cor() correlates columns, so it reads the
    # transpose. A row with no variance has no correlation, which the
    # finiteness check below turns into "not clustered" rather than an error.
    r <- suppressWarnings(stats::cor(t(x), use = "pairwise.complete.obs"))
    stats::as.dist(1 - r)
  } else {
    stats::dist(x, method = dist_method)
  }
  if (!all(is.finite(d))) {
    message("Not clustering the ", axis, "s: some distances are undefined, ",
            "which happens when a pair shares no observation or has no ",
            "variance. The input order is kept.")
    return(NULL)
  }
  stats::hclust(d, method = hclust_method)
}


#' Where the colour key is numbered
#'
#' Sizing the key and drawing it have to agree about which numbers go beside the
#' bar, since it is the widest of them that decides how far right the group
#' legend starts.
#'
#' @param zlim The range the ramp spans.
#'
#' @return The tick positions inside `zlim`.
#'
#' @keywords internal
#' @noRd
sa_key_ticks <- function(zlim) {
  at <- pretty(zlim, 3)
  at[at >= zlim[1] & at <= zlim[2]]
}


#' How wide the colour key needs the strip beside the plot to be
#'
#' Measured from what goes in it, in inches, so that the strip can be reserved
#' before anything is drawn: [graphics::strwidth()]`(units = "inches")` answers on
#' a device that holds no plot yet and answers the same whatever size that device
#' is.
#'
#' @param zlim The range the ramp spans.
#' @param group_cols One colour per group level, named.
#'
#' @return A width in inches.
#'
#' @keywords internal
#' @noRd
sa_key_width <- function(zlim, group_cols, cex.axis, cex.legend) {
  char <- graphics::par("cin")[1]
  wid <- function(s, cex, font = 1) {
    max(graphics::strwidth(s, units = "inches", cex = cex, font = font))
  }
  # The bar, a gap, and the widest of its numbers.
  bar <- (1.5 + 0.4) * char * cex.axis +
    wid(format(sa_key_ticks(zlim), trim = TRUE), cex.axis)
  # The swatches and their labels, or the title if that is the longer line. A
  # swatch and the gaps around it come to 3.3 characters, which is measured from
  # what legend() draws rather than derived: the box it lays out is wider than the
  # swatch and the label, and a strip reserved for the narrower guess cuts the
  # labels off at the edge of the device.
  group <- if (length(group_cols) > 0L) {
    max(3.3 * char * cex.legend + wid(names(group_cols), cex.legend),
        wid("group", cex.legend, font = 2))
  } else {
    0
  }
  gap <- if (length(group_cols) > 0L) 1.2 * char * cex.legend else 0
  0.35 * char + bar + gap + group + 0.35 * char
}


#' How many lines of margin a set of axis labels needs
#'
#' @param labels The labels that go in the margin.
#'
#' @return A width in lines of text, the unit [stats::heatmap()]`(margins = )`
#'   takes, with room for the tick and the gap either side of it.
#'
#' @keywords internal
#' @noRd
sa_text_lines <- function(labels, cex) {
  max(graphics::strwidth(labels, units = "inches", cex = cex)) /
    graphics::par("csi") + 0.6
}
