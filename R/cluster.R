# The second result contract for the unsupervised family. A
# reduction gives every point a coordinate; a clustering gives every point a label.
# Both are answers about the same margin of the same matrix, so this contract keeps
# `points` as its row axis rather than inventing a fifth name: `assignments` is
# aligned with `points` exactly as `scores` is, which is what lets a clustering be
# plotted on top of a reduction of the same rows without either of them being asked
# where its rows came from.
#
# There are four functions rather than one for the reason there are three
# reductions. The four disagree about what a cluster is, and the disagreement is
# the information. `cluster_hclust()` and `cluster_kmeans()` are told how many
# groups to find and will always find that many, so they partition: every point
# lands somewhere and a point in the middle of nowhere lands somewhere anyway.
# `cluster_dbscan()` and `cluster_snn()` are told how dense a group has to be and
# derive the count from that, so they can return two clusters, or nine, or none at
# all, and they can refuse to place a point. A structure all four agree on is a
# different fact from one that only k-means, having been told to find two things,
# found two of.
#
# Which is why the contract has to hold both shapes at once, and it does it with
# one number rather than two slots. `cluster` is an integer per point and `0` is
# noise, which is `dbscan`'s own convention and reads the same way here: a cluster
# number of zero is not a cluster. `clusters` then holds one row per real cluster
# and never a row for zero, so `nrow(clusters)` is the number of groups that were
# found and `design$n_noise` is what did not join one. A partitioning method
# reports `n_noise = 0` and the two shapes read down the same columns.
#
# `silhouette` is here for the same reason `points` is. Coordinates from a t-SNE
# and coordinates from a PCA share no scale and cannot be compared, which is why
# `sa_reduction` promises only the row order; but a silhouette width is a ratio of
# distances in the one matrix all four methods were handed, so it is the number
# that does compare across them. It is computed here rather than left to the
# caller because computing it needs the distance matrix the engine was run on, and
# by the time a result is in the caller's hands that matrix is gone.
#
# `$fit` holds the engine object, the same exception `sa_model` and `sa_reduction`
# make. An `hclust` object is the tree, and the tree is most of what a
# hierarchical clustering has to say — `cutree()` at another height is a question
# the caller should not have to re-run the distance matrix to ask. Drop `$fit` and
# the object writes out as JSON.

#' The clusterings this contract covers
#'
#' `analysis` names the method rather than the family, as it does in `sa_model` and
#' `sa_reduction`, so a result says which of the four it came from without anything
#' having to read `engine`.
#'
#' `"snn"`: shared nearest neighbour clustering on a k-nearest-neighbour graph.
#'
#' @keywords internal
#' @noRd
sa_cluster_analyses <- function() {
  c("hclust", "kmeans", "dbscan", "snn")
}


#' Assemble a clustering result object
#'
#' The checks here guard the contract rather than the user's input, so they fire
#' only on a mistake inside the package and say so. What they are guarding is the
#' promise the object makes: `assignments` is aligned with `points` by position,
#' and `clusters` accounts for every point that is not noise.
#'
#' @param analysis Which clustering this is: `"hclust"`, `"kmeans"`, `"dbscan"` or
#'   `"snn"`.
#' @param points Point labels, the row order every table follows.
#' @param design Named list describing what was clustered.
#' @param parameters Named list of the choices as they were used.
#' @param assignments data.frame of one row per point.
#' @param clusters data.frame of one row per cluster found, noise excluded.
#' @param engine Named list naming what computed the clustering.
#' @param fit The engine object.
#'
#' @keywords internal
#' @noRd
sa_new_cluster <- function(analysis,
                           points,
                           design,
                           parameters,
                           assignments,
                           clusters,
                           engine,
                           fit) {

  if (!analysis %in% sa_cluster_analyses()) {
    stop("internal error: `analysis` must be one of ",
         paste(sa_cluster_analyses(), collapse = ", "), ".", call. = FALSE)
  }
  if (!is.character(points) || length(points) == 0L) {
    stop("internal error: `points` must be a non-empty character vector.",
         call. = FALSE)
  }
  if (!identical(design$point_type, "sample") &&
      !identical(design$point_type, "feature")) {
    stop("internal error: `design$point_type` must be \"sample\" or ",
         "\"feature\".", call. = FALSE)
  }
  if (!is.data.frame(assignments) || !identical(assignments$points, points)) {
    stop("internal error: `assignments` is not a data.frame aligned with ",
         "`points`.", call. = FALSE)
  }
  if (!is.integer(assignments$cluster) || anyNA(assignments$cluster)) {
    stop("internal error: `assignments$cluster` must be an integer vector with ",
         "no missing value; noise is 0.", call. = FALSE)
  }
  if (!is.data.frame(clusters) || any(clusters$cluster == 0L)) {
    stop("internal error: `clusters` must be a data.frame of the clusters that ",
         "were found, which never includes noise.", call. = FALSE)
  }
  # The two tables are one fact counted twice, so they are made to agree here
  # rather than trusted to. Every point is either in a cluster `clusters` has a
  # row for or it is noise, and nothing is both.
  found <- sort(unique(assignments$cluster[assignments$cluster > 0L]))
  if (!identical(clusters$cluster, found)) {
    stop("internal error: `clusters` lists ",
         paste(clusters$cluster, collapse = ", "),
         " but the assignments hold ", paste(found, collapse = ", "), ".",
         call. = FALSE)
  }
  n_noise <- sum(assignments$cluster == 0L)
  if (!identical(design$n_noise, n_noise)) {
    stop("internal error: `design$n_noise` is ", design$n_noise,
         " but ", n_noise, " point(s) were left unassigned.", call. = FALSE)
  }
  if (!identical(design$n_clusters, nrow(clusters))) {
    stop("internal error: `design$n_clusters` is ", design$n_clusters,
         " but `clusters` has ", nrow(clusters), " row(s).", call. = FALSE)
  }
  for (nm in c("package", "method", "label", "overridden")) {
    if (is.null(engine[[nm]])) {
      stop("internal error: `engine` is missing `", nm, "`.", call. = FALSE)
    }
  }

  structure(
    list(analysis    = analysis,
         points      = points,
         design      = design,
         parameters  = parameters,
         assignments = assignments,
         clusters    = clusters,
         engine      = engine,
         fit         = fit,
         metadata    = sa_metadata()),
    class = c("sa_cluster", "sa_result")
  )
}


#' The two tables a clustering reports, built from one vector of labels
#'
#' Every one of the four methods ends with an integer per point and a distance
#' matrix, and everything the result says about the grouping is derived from those
#' two here. Doing it in one place is what keeps `cluster` meaning the same thing
#' in all four: the labels are renumbered from 1 in order of first appearance, so
#' that a k-means run and a DBSCAN run on the same points name their groups by the
#' same rule and neither inherits whatever the engine happened to count from.
#'
#' @param cluster Integer labels as the engine returned them, `0` for noise.
#' @param points Point labels.
#' @param d The [stats::dist()] object the engine was run on, for the silhouette.
#'
#' @return A list with `assignments`, `clusters`, `n_clusters` and `n_noise`.
#'
#' @keywords internal
#' @noRd
sa_cluster_tables <- function(cluster, points, d) {
  cluster <- as.integer(cluster)
  found <- unique(cluster[cluster > 0L])
  # Renumbering by first appearance, not by size: a cluster is not more itself for
  # being large, and first appearance is the one order that does not change when
  # the caller reorders their rows for an unrelated reason.
  renumbered <- integer(length(cluster))
  renumbered[cluster > 0L] <- match(cluster[cluster > 0L], found)

  sil <- sa_silhouette(d, renumbered)

  assignments <- data.frame(
    points     = points,
    cluster    = renumbered,
    silhouette = sil,
    stringsAsFactors = FALSE
  )
  rownames(assignments) <- NULL

  ids <- seq_along(found)
  clusters <- data.frame(
    cluster    = ids,
    size       = as.integer(tabulate(renumbered, nbins = length(found))),
    silhouette = vapply(ids, function(i) {
      mean(sil[renumbered == i], na.rm = TRUE)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  rownames(clusters) <- NULL
  # A single cluster has no other cluster to be far from, so its silhouette is
  # undefined rather than zero, and `mean()` of an all-`NA` group returns `NaN`.
  clusters$silhouette[is.nan(clusters$silhouette)] <- NA_real_

  list(assignments = assignments,
       clusters    = clusters,
       n_clusters  = nrow(clusters),
       n_noise     = sum(renumbered == 0L))
}


#' Print a clustering
#'
#' Summarises what was clustered and what came of it, rather than printing the
#' label of every point. Those are in `x$assignments`, and the engine object is
#' `x$fit`.
#'
#' @param x A clustering, as returned by [cluster_hclust()], [cluster_kmeans()],
#'   [cluster_dbscan()] or [cluster_snn()].
#' @param n Maximum number of clusters to report the size of. The rest are
#'   counted.
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' cluster_kmeans(iris[1:4], n_clust = 3, seed = 1)
#'
#' @export
print.sa_cluster <- function(x, n = 10L, ...) {
  n <- sa_check_count(n, "n", 0)
  design <- x$design
  params <- x$parameters

  cat("<sa_cluster> ", x$analysis, "\n", sep = "")
  cat("  data     : ", design$n_used, " sample(s) x ", design$n_feats,
      " feature(s)",
      if (design$n_dropped > 0L) {
        paste0("  (", design$n_dropped, " incomplete row(s) dropped)")
      },
      "\n", sep = "")
  cat("  points   : ", length(x$points), " ", design$point_type, "(s)\n",
      sep = "")
  cat("  scaling  : ",
      if (params$center && params$scale) {
        "centred and scaled"
      } else if (params$center) {
        "centred"
      } else if (params$scale) {
        "scaled"
      } else {
        "none, values as they arrived"
      },
      "\n", sep = "")

  # The count leads, because for two of the four methods it is the answer rather
  # than something that was asked for.
  cat("  clusters : ", design$n_clusters,
      if (design$n_noise > 0L) {
        paste0("  (", design$n_noise, " point(s) left as noise)")
      },
      "\n", sep = "")

  if (design$n_clusters > 0L) {
    shown <- utils::head(x$clusters, n)
    sa_cat_field("sizes", paste0(
      paste0("#", shown$cluster, " n = ", shown$size,
             ", s = ", sa_fmt_num(shown$silhouette, 3), collapse = "; "),
      if (nrow(shown) < design$n_clusters) {
        paste0("  (", nrow(shown), " of ", design$n_clusters, " shown)")
      }
    ))
    sa_cat_field("silhouette", paste0(
      "mean ", sa_fmt_num(mean(x$assignments$silhouette, na.rm = TRUE), 3),
      " over the ", sum(!is.na(x$assignments$silhouette)), " assigned ",
      design$point_type, "(s), on the ", params$dist_method, " distance"
    ))
  }

  if (identical(x$analysis, "hclust")) {
    cat("  linkage  : ", params$hclust_method, ", cut at k = ", params$n_clust,
        "\n", sep = "")
  }
  if (identical(x$analysis, "kmeans")) {
    cat("  kmeans   : k = ", params$n_clust, ", ", params$n_start,
        " start(s), ", sa_fmt_num(params$tot_withinss, 5), " within-cluster ss",
        if (!is.null(params$seed)) paste0("  (seed = ", params$seed, ")"),
        "\n", sep = "")
  }
  if (identical(x$analysis, "dbscan")) {
    cat("  density  : eps = ", sa_fmt_num(params$eps, 4), " (",
        params$eps_source, "), min_pts = ", params$min_pts, "\n", sep = "")
  }
  if (identical(x$analysis, "snn")) {
    cat("  density  : k = ", params$k, ", eps = ", params$eps,
        " shared neighbour(s), min_pts = ", params$min_pts, "\n", sep = "")
  }
  if (length(design$dropped_feats) > 0L) {
    sa_cat_field("dropped", paste0(
      paste(design$dropped_feats, collapse = ", "), " (no variance)"
    ))
  }

  invisible(x)
}
