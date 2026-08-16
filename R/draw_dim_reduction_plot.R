# The picture the unsupervised family ends at. A reduction gives every point a
# coordinate and a clustering gives every point a label, and because both read
# their input through `sa_reduce_input()` the two are about the same rows in the
# same order. That is the whole reason this function can exist without asking
# either of them where their rows came from.
#
# Colour and shape are given to two different things on purpose. A clustering is
# what the data was found to say and a `group` is what the caller already knew,
# and the interesting question is almost always whether the two agree. Drawing
# them on one channel would answer it by hiding it: one would have to be left
# out, or the two would have to be crossed into a single set of labels whose
# count is the product of theirs. On two channels the agreement is read directly
# — one colour per shape is a clustering that recovered the groups, and a shape
# split across colours is a group the data does not see as one thing.
#
# Noise is grey rather than a palette colour. `cluster_dbscan()` and
# `cluster_snn()` can leave a point out, and a point left out is the absence of a
# cluster rather than a cluster of its own; giving it a colour beside the others
# would put it in the legend as though it were one.


#' The shapes a group is drawn with, in the order they are handed out
#'
#' Filled first, since a filled marker is the one that reads at a small size, and
#' then hollow. Ten is where it stops: past that the shapes are telling each other
#' apart rather than telling the groups apart.
#'
#' @keywords internal
#' @noRd
sa_scatter_pch <- function() {
  c(16L, 17L, 15L, 18L, 8L, 1L, 2L, 0L, 5L, 6L)
}


#' The two coordinates a reduction is drawn on, and what to call their axes
#'
#' @param reduction_result The caller's first argument, not yet checked.
#' @param dims Which two of the coordinates to draw.
#'
#' @return A list with `x`, `y`, `xlab`, `ylab` and `points`.
#'
#' @keywords internal
#' @noRd
sa_reduction_scatter <- function(reduction_result, dims) {
  # A clustering is the other half of the same pipeline and is easy to reach for
  # first, so it is turned away by name rather than by class.
  if (inherits(reduction_result, "sa_cluster")) {
    stop("`reduction_result` is a clustering, which gives every point a label ",
         "and no coordinate. Reduce the same frame with perform_pca(), ",
         "perform_tsne() or perform_umap() and pass the clustering as ",
         "`cluster_result`.", call. = FALSE)
  }
  if (!inherits(reduction_result, "sa_reduction")) {
    stop("`reduction_result` must be a reduction, as returned by ",
         "perform_pca(), perform_tsne() or perform_umap().", call. = FALSE)
  }

  scores <- reduction_result$scores
  coords <- setdiff(names(scores), "points")

  if (!is.numeric(dims) || length(dims) != 2L) {
    stop("`dims` must be a numeric vector of length 2, naming the two ",
         "coordinates to draw.", call. = FALSE)
  }
  dims <- c(sa_check_count(dims[1], "dims", 1),
            sa_check_count(dims[2], "dims", 1))
  if (dims[1] == dims[2]) {
    stop("`dims` must name two different coordinates, but names ", dims[1],
         " twice.", call. = FALSE)
  }
  if (max(dims) > length(coords)) {
    stop("`dims` asks for coordinate ", max(dims), " of a ",
         reduction_result$analysis, " that has ", length(coords), " (",
         paste(coords, collapse = ", "), "). An embedding has as many as ",
         "`n_dim` asked for.", call. = FALSE)
  }

  drawn <- coords[dims]
  # Only a rotation has a share of the variance to report, and the share is what
  # makes one of its axes wider than another. An embedding's coordinates carry no
  # such number, so their labels are the names and nothing more.
  labels <- drawn
  if (!is.null(reduction_result$variance)) {
    variance <- reduction_result$variance
    at <- match(drawn, variance$component)
    labels <- paste0(drawn, " (", sa_fmt_num(variance$prop_var[at], 3), "%)")
  }

  list(x      = scores[[drawn[1]]],
       y      = scores[[drawn[2]]],
       xlab   = labels[1],
       ylab   = labels[2],
       points = scores$points)
}


#' One colour per cluster, and grey for what did not join one
#'
#' The alignment check is the contract both objects already promise rather than
#' anything new: `sa_cluster` and `sa_reduction` repeat the same `points` in the
#' same order, so a mismatch here means two different frames were reduced and
#' clustered rather than one.
#'
#' @param cluster_result The caller's `cluster_result`, not yet checked.
#' @param points The reduction's point labels.
#' @param col The caller's `col`.
#' @param cluster_lv The caller's `cluster_lv`.
#' @param theme The resolved theme, for the noise colour.
#'
#' @return A list with `cluster`, `col`, `palette`, `levels` and `n_noise`.
#'
#' @keywords internal
#' @noRd
sa_scatter_clusters <- function(cluster_result, points, col, cluster_lv, theme) {
  if (!inherits(cluster_result, "sa_cluster")) {
    stop("`cluster_result` must be a clustering, as returned by ",
         "cluster_hclust(), cluster_kmeans(), cluster_dbscan() or ",
         "cluster_snn().", call. = FALSE)
  }
  if (!identical(cluster_result$assignments$points, points)) {
    stop("`cluster_result` and `reduction_result` describe different points. ",
         "Both read their input through the same function, so this means they ",
         "were given different frames, different `feats`, or one of them the ",
         "sample scale and the other the feature scale.", call. = FALSE)
  }

  cluster <- cluster_result$assignments$cluster
  n <- cluster_result$design$n_clusters
  palette <- if (is.null(col)) {
    grDevices::hcl.colors(max(n, 2L), "Dark 2")[seq_len(n)]
  } else {
    if (length(col) != 1L && length(col) != n) {
      stop("`col` must hold one colour, or one per cluster (", n, ").",
           call. = FALSE)
    }
    rep(col, length.out = n)
  }

  # Indexing with the labels directly would drop the noise points rather than
  # colour them, since `palette[0]` is nothing at all.
  point_col <- character(length(cluster))
  point_col[cluster == 0L] <- theme$guide
  point_col[cluster > 0L] <- palette[cluster[cluster > 0L]]

  levels <- if (is.null(cluster_lv)) {
    paste0("#", seq_len(n))
  } else {
    cluster_lv <- as.character(cluster_lv)
    if (length(cluster_lv) != n) {
      stop("`cluster_lv` must hold one label per cluster (", n, ").",
           call. = FALSE)
    }
    if (anyDuplicated(cluster_lv) > 0L) {
      stop("`cluster_lv` must not repeat a level.", call. = FALSE)
    }
    cluster_lv
  }

  list(cluster = cluster,
       col     = point_col,
       palette = palette,
       levels  = levels,
       n_noise = cluster_result$design$n_noise)
}


#' One shape per group level
#'
#' `group_lv` sets the order the levels are drawn and read in, and nothing else.
#' It does not select rows, which is where this parts company with
#' [draw_grouped_boxplot()]: over there a level left out is a level left out of
#' the analysis, and here the reduction has already placed every point, so a
#' point whose group went unlisted would vanish from a picture it belongs in.
#'
#' @param group The caller's `group`.
#' @param group_lv The caller's `group_lv`.
#' @param points The reduction's point labels.
#' @param design The reduction's `design`, for the message about dropped rows.
#' @param col The caller's `col`, or `NULL` when the clustering takes it.
#' @param pch The caller's `pch`.
#'
#' @return A list with `group`, `pch`, `levels`, `pch_lv`, and `palette` and
#'   `col` when `col` was given.
#'
#' @keywords internal
#' @noRd
sa_scatter_groups <- function(group, group_lv, points, design, col = NULL,
                             pch = NULL) {
  if (length(group) != length(points)) {
    # The likeliest way to get here is to pass the grouping column of the frame
    # that was reduced after the reduction had dropped rows out of it, so that
    # case names itself instead of leaving two numbers to be compared.
    extra <- if (design$n_dropped > 0L && length(group) == design$n_samples) {
      paste0(" The reduction dropped ", design$n_dropped, " of the ",
             design$n_samples, " row(s) it was given, so it holds fewer ",
             "points than `data` had rows. `reduction_result$points` names ",
             "the ones that are left.")
    } else if (identical(design$point_type, "feature")) {
      paste0(" These points are features rather than samples, so `group` ",
             "labels the ", length(points), " feature(s) that were embedded.")
    } else {
      ""
    }
    stop("`group` must hold one label per point, so ", length(points),
         " of them, but holds ", length(group), ".", extra, call. = FALSE)
  }

  group_chr <- as.character(group)
  if (anyNA(group_chr)) {
    stop("`group` must not hold a missing value. Every point the reduction ",
         "placed is drawn, so every one of them needs a shape.", call. = FALSE)
  }

  if (is.null(group_lv)) {
    group_lv <- if (is.factor(group)) levels(group) else sort(unique(group_chr))
    group_lv <- group_lv[group_lv %in% group_chr]
  } else {
    group_lv <- as.character(group_lv)
    if (anyDuplicated(group_lv) > 0L) {
      stop("`group_lv` must not repeat a level.", call. = FALSE)
    }
    unlisted <- setdiff(unique(group_chr), group_lv)
    if (length(unlisted) > 0L) {
      stop("`group_lv` leaves out ", paste(unlisted, collapse = ", "),
           ", which `group` uses. Here it sets the order the levels are drawn ",
           "and read in rather than which rows are kept: every point the ",
           "reduction placed is drawn, so a level left out is a point with no ",
           "shape rather than a point removed.", call. = FALSE)
    }
  }

  shapes <- if (is.null(pch)) {
    sa_scatter_pch()
  } else {
    if (!is.numeric(pch) || any(is.na(pch))) {
      stop("`pch` must be numeric and must not hold a missing value.",
           call. = FALSE)
    }
    if (length(pch) != 1L && length(pch) != length(group_lv)) {
      stop("`pch` must hold one shape, or one per group level (",
           length(group_lv), ").", call. = FALSE)
    }
    rep(pch, length.out = length(group_lv))
  }
  if (is.null(pch) && length(group_lv) > length(shapes)) {
    stop("`group` has ", length(group_lv), " levels and there are ",
         length(shapes), " shapes to tell them apart with. Past that the ",
         "shapes are distinguishing themselves rather than the groups; ",
         "name `pch` with one shape per level, or use `cluster_result` for ",
         "the colours.", call. = FALSE)
  }

  palette <- if (is.null(col)) {
    NULL
  } else {
    if (length(col) != 1L && length(col) != length(group_lv)) {
      stop("`col` must hold one colour, or one per group level (",
           length(group_lv), ").", call. = FALSE)
    }
    rep(col, length.out = length(group_lv))
  }

  at <- match(group_chr, group_lv)
  list(group   = factor(group_chr, levels = group_lv),
       pch     = shapes[at],
       levels  = group_lv,
       pch_lv  = shapes[seq_along(group_lv)],
       palette = palette,
       col     = if (is.null(palette)) NULL else palette[at])
}


#' Draw a reduction as a scatter of its points
#'
#' Plots two coordinates of a [perform_pca()], [perform_tsne()] or
#' [perform_umap()] result against each other. A clustering of the same frame
#' colours the points and a known grouping shapes them, so what the data was
#' found to say and what the caller already knew are read off one picture.
#'
#' @details
#' # Colour and shape say two different things
#'
#' `cluster_result` takes the colours and `group` takes the shapes when both are
#' given, and giving both is the point of the arrangement rather than a conflict
#' to be resolved. With one channel alone, that channel takes the colours when
#' `col` is named; otherwise a lone clustering is coloured and a lone grouping
#' is shaped. The question a clustering usually raises is whether it recovered a
#' grouping that was known all along, and on two channels that is read directly:
#' one colour per shape is a clustering that found the groups, and a shape split
#' across colours is a group the data does not see as one thing. Give neither
#' and the points are drawn in the foreground colour.
#'
#' The clustering has to be of the same points, which the two contracts already
#' promise: `sa_cluster` and `sa_reduction` both read their input through the
#' same function, so a mismatch is refused rather than lined up by position.
#'
#' # Noise
#'
#' [cluster_dbscan()] and [cluster_snn()] can leave a point in no cluster at all.
#' Those points are grey rather than a palette colour, since a point left out is
#' the absence of a cluster and not a cluster of its own, and the legend counts
#' them on a line of their own.
#'
#' # The axes
#'
#' A principal component analysis reports what share of the variance each of its
#' components carries, so its axis labels carry it too: `PC1 (30.2%)` is read
#' from `$variance` rather than recomputed. An embedding has no such number and
#' its axes are labelled with their names alone. `asp = 1` is what makes one unit
#' of the vertical axis the same length as one unit of the horizontal, which is
#' worth setting when the distance between two points is what is being read.
#'
#' @param reduction_result A reduction, as returned by [perform_pca()],
#'   [perform_tsne()] or [perform_umap()].
#' @param group One label per point, or `NULL` for no shape. Points are the rows
#'   the reduction kept, which `reduction_result$points` names, and on the
#'   feature scale they are features rather than samples.
#' @param group_lv The levels of `group` in the order they are drawn and listed
#'   in, or `NULL` for the factor's own order and otherwise alphabetical. Unlike
#'   [draw_grouped_boxplot()]'s, this argument selects no rows: a level `group`
#'   uses and `group_lv` leaves out is an error rather than a point dropped from
#'   the picture.
#' @param cluster_result A clustering of the same points, as returned by
#'   [cluster_hclust()], [cluster_kmeans()], [cluster_dbscan()] or
#'   [cluster_snn()], or `NULL` for no colouring.
#' @param cluster_lv One label per cluster, in the order the clustering numbers
#'   them, or `NULL` for `#1`, `#2`, and so on. Noise is still listed as
#'   `noise (n)` when present.
#' @param dims Which two coordinates to draw, as positions in the score table.
#' @param anno_points Whether to write each point's label beside it.
#' @param dark Whether to draw on a dark background.
#' @param asp Aspect ratio passed to [graphics::plot()], or `NULL` to let the
#'   axes fill the panel independently. `1` makes distances comparable between
#'   the two axes.
#' @param col One colour for every point, one per cluster when `cluster_result`
#'   is given, one per group level when only `group` is given, or `NULL` for
#'   `hcl.colors(n, "Dark 2")` on clusters and the foreground colour elsewhere.
#'   Noise is grey whatever this says.
#' @param pch One shape for every point, one per group level when `group` is
#'   given, or `NULL` for `16` and the default shape sequence respectively.
#' @param cex Size of the plotted points.
#' @param xlim,ylim Axis ranges, or `NULL` to take them from the coordinates.
#' @param xlab,ylab,main Axis and figure labels. `NULL` builds them from the
#'   reduction.
#' @param cex.axis,cex.lab,cex.main,cex.legend,cex.anno Relative text sizes.
#'   `cex.anno` sizes the point labels and `NULL` matches `cex.legend`.
#'
#' @return A data.frame of the points as they were drawn, invisibly: `points`,
#'   `x`, `y`, `col`, `pch`, and `cluster` and `group` when those were given. The
#'   resolved colouring is carried as a `"view"` attribute, one of `"both"`,
#'   `"cluster"`, `"group"` or `"plain"`.
#'
#' @seealso [perform_pca()] for the coordinates and [cluster_kmeans()] for the
#'   labels, and [draw_heatmap()], which shows the same wide input a cell at a
#'   time rather than a point at a time.
#'
#' @examples
#' res <- perform_pca(iris[1:4])
#'
#' ## What was known all along, as shapes. Name `col` to colour the levels
#' ## instead.
#' draw_dim_reduction_plot(res, group = iris$Species)
#'
#' draw_dim_reduction_plot(res, group = iris$Species,
#'                       col = c("#E69F00", "#56B4E9", "#009E73"))
#'
#' ## What the data was found to say, as colours. The clustering was never shown
#' ## the species, so one colour per shape is its own finding.
#' cl <- cluster_kmeans(iris[1:4], n_clust = 3, seed = 1)
#' drawn <- draw_dim_reduction_plot(res, group = iris$Species,
#'                                  cluster_result = cl)
#' table(cluster = drawn$cluster, species = drawn$group)
#'
#' ## Names for the clusters the legend would otherwise number.
#' draw_dim_reduction_plot(res, cluster_result = cl,
#'                       cluster_lv = c("A", "B", "C"))
#'
#' ## A density method can place no point at all, and those are grey.
#' db <- cluster_dbscan(iris[1:4])
#' draw_dim_reduction_plot(res, cluster_result = db, asp = 1)
#'
#' ## The features as the points, labelled, on the third and fourth components.
#' by_feat <- perform_pca(iris[1:4], embedding_scale = "features")
#' draw_dim_reduction_plot(by_feat, dims = c(1, 3), anno_points = TRUE)
#'
#' @export
draw_dim_reduction_plot <- function(reduction_result,
                                    group = NULL,
                                    group_lv = NULL,
                                    cluster_result = NULL,
                                    cluster_lv = NULL,
                                    dims = c(1L, 2L),
                                    anno_points = FALSE,
                                    dark = FALSE,
                                    asp = NULL,
                                    col = NULL,
                                    pch = NULL,
                                    cex = 1.2,
                                    xlim = NULL,
                                    ylim = NULL,
                                    xlab = NULL,
                                    ylab = NULL,
                                    main = NULL,
                                    cex.axis = 1.2,
                                    cex.lab = 1.3,
                                    cex.main = 1.3,
                                    cex.legend = 1.1,
                                    cex.anno = NULL) {
  space <- sa_reduction_scatter(reduction_result, dims)
  sa_check_flag(anno_points, "anno_points")
  sa_check_flag(dark, "dark")
  sa_check_lim(xlim, "xlim")
  sa_check_lim(ylim, "ylim")
  sa_check_scalar_num(cex, "cex", 0, lower_open = TRUE)
  if (!is.null(asp)) {
    sa_check_scalar_num(asp, "asp", 0, lower_open = TRUE)
  }
  if (is.null(cex.anno)) {
    cex.anno <- cex.legend
  } else {
    sa_check_scalar_num(cex.anno, "cex.anno", 0, lower_open = TRUE)
  }
  if (is.null(group) && !is.null(group_lv)) {
    stop("`group_lv` names the levels of `group`, which was not given.",
         call. = FALSE)
  }
  if (is.null(cluster_result) && !is.null(cluster_lv)) {
    stop("`cluster_lv` names the levels of the clustering, which was not given.",
         call. = FALSE)
  }

  theme <- sa_plot_theme(dark)
  points_lab <- space$points

  clusters <- if (is.null(cluster_result)) {
    NULL
  } else {
    sa_scatter_clusters(cluster_result, points_lab, col, cluster_lv, theme)
  }
  groups <- if (is.null(group)) {
    NULL
  } else {
    sa_scatter_groups(group, group_lv, points_lab, reduction_result$design,
                      col = if (is.null(cluster_result)) col else NULL,
                      pch = pch)
  }

  # With no clustering and no group colouring there is nothing for a palette to
  # be one of, so `col` is a single colour for every point and the default is
  # whatever the foreground is on this background.
  point_col <- if (!is.null(clusters)) {
    clusters$col
  } else if (!is.null(groups) && !is.null(groups$col)) {
    groups$col
  } else if (is.null(col)) {
    rep(theme$fg, length(points_lab))
  } else {
    if (length(col) != 1L && length(col) != length(points_lab)) {
      stop("`col` must hold one colour, or one per point (",
           length(points_lab), "), when no `cluster_result` or `group` ",
           "colouring applies.", call. = FALSE)
    }
    rep(col, length.out = length(points_lab))
  }
  point_pch <- if (!is.null(groups)) {
    groups$pch
  } else if (is.null(pch)) {
    rep(16L, length(points_lab))
  } else {
    if (!is.numeric(pch) || any(is.na(pch))) {
      stop("`pch` must be numeric and must not hold a missing value.",
           call. = FALSE)
    }
    if (length(pch) != 1L && length(pch) != length(points_lab)) {
      stop("`pch` must hold one shape, or one per point (",
           length(points_lab), "), when no `group` is given.", call. = FALSE)
    }
    rep(pch, length.out = length(points_lab))
  }

  view <- if (!is.null(clusters) && !is.null(groups)) {
    "both"
  } else if (!is.null(clusters)) {
    "cluster"
  } else if (!is.null(groups)) {
    "group"
  } else {
    "plain"
  }

  span_x <- if (is.null(xlim)) range(space$x, finite = TRUE) else xlim
  span_y <- if (is.null(ylim)) range(space$y, finite = TRUE) else ylim

  # The two legend blocks, built before the layout because their widest line is
  # what decides how much of the figure the legend panel takes.
  cluster_key <- if (is.null(clusters)) {
    NULL
  } else {
    list(labels = c(clusters$levels,
                    if (clusters$n_noise > 0L) {
                      paste0("noise (", clusters$n_noise, ")")
                    }),
         cols   = c(clusters$palette,
                    if (clusters$n_noise > 0L) theme$guide))
  }
  group_key <- if (is.null(groups)) {
    NULL
  } else {
    list(labels = groups$levels,
         cols   = groups$palette,
         pch    = groups$pch_lv)
  }

  legend_text <- c(if (!is.null(cluster_key)) cluster_key$labels,
                   if (!is.null(group_key)) group_key$labels,
                   if (!is.null(cluster_key)) "cluster",
                   if (!is.null(group_key)) "group")

  old_par <- graphics::par(c("bg", "fg", "col.axis", "col.lab", "col.main",
                             "mar", "mfrow"))
  # Only what this function sets is put back. A par(no.readonly = TRUE) snapshot
  # also carries `fin`, `pin` and `mai`, which are absolute sizes, so restoring
  # it pins the figure to the size this plot happened to be drawn at. `mfrow` is
  # in the list because layout() overwrites whatever grid the caller had.
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)
  graphics::par(bg = theme$bg, fg = theme$fg, col.axis = theme$fg,
                col.lab = theme$fg, col.main = theme$fg)

  if (length(legend_text) > 0L) {
    wanted <- (max(nchar(legend_text)) + 5) * cex.legend *
      graphics::par("cin")[1] + 0.2
    share <- min(0.35, wanted / graphics::par("din")[1])
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1 - share, share))
  }
  graphics::par(mar = c(5.1, 4.6, 4.1, 1.1))

  graphics::plot.default(
    space$x, space$y, type = "n", bty = "n",
    xlim = span_x, ylim = span_y, asp = if (is.null(asp)) NA else asp,
    xlab = if (is.null(xlab)) space$xlab else xlab,
    ylab = if (is.null(ylab)) space$ylab else ylab,
    main = if (is.null(main)) reduction_result$engine$label else main,
    cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main
  )
  graphics::points(space$x, space$y, pch = point_pch, col = point_col,
                   cex = cex)
  if (anno_points) {
    graphics::text(space$x, space$y, labels = points_lab, pos = 3,
                   offset = 0.4, cex = cex.anno, col = theme$fg)
  }

  if (length(legend_text) > 0L) {
    graphics::par(mar = c(5, 0, 4, 1))
    graphics::plot.new()
    # Two blocks anchored to the two ends of the panel rather than one after the
    # other, since colour and shape are separate readings and a gap between them
    # is what says so.
    both <- !is.null(cluster_key) && !is.null(group_key)
    if (!is.null(cluster_key)) {
      graphics::legend(if (both) "top" else "center", title = "cluster",
                       legend = cluster_key$labels, col = cluster_key$cols,
                       pch = 15, bty = "n", cex = cex.legend,
                       text.col = theme$fg, title.col = theme$fg)
    }
    if (!is.null(group_key)) {
      graphics::legend(if (both) "bottom" else "center", title = "group",
                       legend = group_key$labels,
                       col = if (is.null(group_key$cols)) {
                         theme$fg
                       } else {
                         group_key$cols
                       },
                       pch = group_key$pch, bty = "n", cex = cex.legend,
                       text.col = theme$fg, title.col = theme$fg)
    }
  }

  drawn <- data.frame(points = points_lab, x = space$x, y = space$y,
                      stringsAsFactors = FALSE)
  if (!is.null(clusters)) {
    drawn$cluster <- clusters$cluster
  }
  if (!is.null(groups)) {
    drawn$group <- groups$group
  }
  drawn$col <- point_col
  drawn$pch <- point_pch
  rownames(drawn) <- NULL

  # Which of the four readings this was is decided from which arguments arrived,
  # so it is carried on the result rather than left to be inferred from which
  # columns are present.
  attr(drawn, "view") <- view
  invisible(drawn)
}


#' @param x A reduction, as returned by [perform_pca()], [perform_tsne()] or
#'   [perform_umap()].
#' @param ... Arguments passed on to [draw_dim_reduction_plot()].
#'
#' @rdname draw_dim_reduction_plot
#' @export
plot.sa_reduction <- function(x, ...) {
  draw_dim_reduction_plot(x, ...)
}
