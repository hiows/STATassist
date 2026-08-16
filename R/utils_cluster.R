# Everything the four clustering functions share, which is the input side and the
# distance. They differ in what they do with a matrix of points and not at all in
# how they read one, so the same rows are dropped for the same reasons in all four
# and `design` describes the input the same way. That is what makes the four
# comparable, and comparing them is most of the point of having four: a group that
# only k-means finds, having been told to find two things, is a different fact
# from one that DBSCAN found without being told how many to look for.
#
# The input is read through `utils_reduce.R` rather than through a copy of it, so
# a clustering and a reduction of the same frame are about the same rows and the
# assignment can be painted straight onto the scores.
#
# `cluster_scale` is the `embedding_scale` of the reductions under the name that
# fits here, and it means the same thing: which margin of the input becomes a
# point. `center` and `scale` always apply to the features whichever margin that
# is, because that is what scaling a data set means, and standardising the
# transpose instead standardises samples and answers a third question.

#' Read the matrix a clustering is computed on out of the caller's frame
#'
#' The reduction helpers do the work; this puts them in the one order all four
#' clustering functions use, so that none of them can assemble the matrix its own
#' way. What comes back is the point-by-variable matrix the engine is handed,
#' beside everything `design` needs to describe the input it came from.
#'
#' @param data The caller's `data`.
#' @param feats The caller's `feats`, possibly `NULL`.
#' @param cluster_scale Which margin becomes a point, already matched.
#' @param center,scale The caller's two flags, applied to the features.
#' @param fn Name of the calling function, for the message about a table too small
#'   to work with.
#'
#' @return A list with `m`, the numeric matrix of points; `points` and
#'   `point_type`; and the counts `n_samples`, `n_used`, `n_dropped`, `n_feats`
#'   beside `feats` and `dropped_feats`.
#'
#' @keywords internal
#' @noRd
sa_cluster_input <- function(data, feats, cluster_scale, center, scale, fn) {
  input <- sa_reduce_input(data, feats, scale, fn)
  pt <- sa_reduce_points(input$x, input$samples, cluster_scale)
  m <- sa_reduce_embedding_matrix(input$x, cluster_scale, center, scale)
  # The engine object is part of the result, so it should carry labels a reader
  # recognises rather than the positions it was handed.
  rownames(m) <- pt$points

  list(m             = m,
       points        = pt$points,
       point_type    = pt$point_type,
       n_samples     = input$n_samples,
       n_used        = nrow(input$x),
       n_dropped     = input$n_dropped,
       n_feats       = ncol(input$x),
       feats         = colnames(input$x),
       dropped_feats = input$dropped_feats)
}


#' The distance a clustering is computed on
#'
#' `"correlation"` is `1 - cor()`, which groups points by the shape of their
#' profile rather than by how high it sits, and the other two are handed to
#' [stats::dist()]. [draw_heatmap()] reads this too, so the tree it draws and the
#' tree [cluster_hclust()] cuts are the same tree when they are asked for on the
#' same terms.
#'
#' @param x Objects to measure between, in rows.
#' @param dist_method `"euclidean"`, `"manhattan"` or `"correlation"`.
#'
#' @return A [stats::dist()] object.
#'
#' @keywords internal
#' @noRd
sa_cluster_dist <- function(x, dist_method) {
  if (identical(dist_method, "correlation")) {
    # Rows are the objects and cor() correlates columns, so it reads the
    # transpose. A row with no variance has no correlation, which is left to the
    # caller to notice rather than turned into an error here.
    r <- suppressWarnings(stats::cor(t(x), use = "pairwise.complete.obs"))
    stats::as.dist(1 - r)
  } else {
    stats::dist(x, method = dist_method)
  }
}


#' How many clusters to cut to
#'
#' The two partitioning methods are told the count, so the count is theirs to
#' check. One cluster is not a clustering and there cannot be more clusters than
#' points, and both ends are refused here rather than left to the engine, which
#' says it in terms of its own arguments.
#'
#' @param n_clust The caller's `n_clust`.
#' @param n Number of points being clustered.
#' @param point_type What one point is, so the limit is quoted in the caller's
#'   terms.
#'
#' @return `n_clust` as an integer.
#'
#' @keywords internal
#' @noRd
sa_cluster_n_clust <- function(n_clust, n, point_type) {
  n_clust <- sa_check_count(n_clust, "n_clust", 2)
  if (n_clust > n) {
    stop("`n_clust` must not exceed the ", n, " usable ", point_type,
         "(s) being clustered, but is ", n_clust, ". There cannot be more ",
         "groups than there are things to put in them.", call. = FALSE)
  }
  n_clust
}


#' How dense a neighbourhood has to be before DBSCAN calls it a cluster
#'
#' The textbook floor is one more than the number of dimensions: a group in `d`
#' dimensions needs `d + 1` points to be more than a flat piece of one. That rule
#' is unusable on the shape of table this package is usually given, where there are
#' more features than samples and `d + 1` is larger than the whole data set, so it
#' is capped at half the points. A threshold above half can only ever return a
#' single cluster, since two groups that size do not fit.
#'
#' It is floored at 4 rather than at `d + 1` because `d + 1` is 3 on a
#' two-dimensional table, and 3 fragments: Ester et al. (1996), who introduced the
#' method, settled on 4 for two dimensions and found nothing bought by going
#' higher. On the two-blob probe in `Test/cursor_test/2026_08_15/` the difference
#' is four clusters and twelve points of noise against two clusters and none.
#'
#' The derived value is said out loud. It is the whole behaviour of the method and
#' someone clustering 40 features is the last person who would think to check what
#' it came out as.
#'
#' @param min_pts The caller's `min_pts`, possibly `NULL`.
#' @param n Number of points being clustered.
#' @param n_var Number of variables each point is described by.
#' @param point_type What one point is.
#'
#' @return The threshold to run at.
#'
#' @keywords internal
#' @noRd
sa_dbscan_min_pts <- function(min_pts, n, n_var, point_type) {
  if (!is.null(min_pts)) {
    min_pts <- sa_check_count(min_pts, "min_pts", 2)
    if (min_pts > n) {
      stop("`min_pts` must not exceed the ", n, " usable ", point_type,
           "(s) being clustered, but is ", min_pts,
           ". No neighbourhood can hold more points than there are.",
           call. = FALSE)
    }
    return(min_pts)
  }
  derived <- min(n, max(4L, min(as.integer(n_var) + 1L, n %/% 2L)))
  message("Using min_pts = ", derived, ", from the ", n_var,
          " variable(s) describing each ", point_type,
          " and the ", n, " being clustered. Pass `min_pts` to set it.")
  derived
}


#' The neighbourhood shared nearest neighbour clustering is run at
#'
#' Three arguments that only mean anything together, so they are resolved
#' together. `k` is how many neighbours each point keeps, `eps` is how many of them
#' two points must have in common before there is an edge between them, and
#' `min_pts` is how many such edges a point needs before it is a core point. All
#' three are counts of neighbours, which is why `eps` here is nothing like the
#' `eps` of [cluster_dbscan()]: that one is a radius in the units of the data.
#'
#' @param k,eps,min_pts The caller's three, each possibly `NULL`.
#' @param n Number of points being clustered.
#' @param point_type What one point is.
#'
#' @return A list of the three, resolved.
#'
#' @keywords internal
#' @noRd
sa_snn_params <- function(k, eps, min_pts, n, point_type) {
  # Everything that was supplied is checked for being the right kind of thing
  # before anything is derived, so that a call that is going to fail does not
  # first announce a default it never used. The one check that cannot come first
  # is `eps` against `k`, since until `k` is resolved there is nothing to check it
  # against.
  if (!is.null(k)) k <- sa_check_count(k, "k", 2)
  if (!is.null(eps)) eps <- sa_check_count(eps, "eps", 1)
  if (!is.null(min_pts)) min_pts <- sa_check_count(min_pts, "min_pts", 1)

  if (is.null(k)) {
    k <- min(n - 1L, max(3L, as.integer(ceiling(sqrt(n)))))
    message("Using k = ", k, " neighbour(s), from the ", n, " ", point_type,
            "(s) being clustered. Pass `k` to set it.")
  } else if (k > n - 1L) {
    stop("`k` must not exceed one less than the ", n, " usable ", point_type,
         "(s) being clustered, which is ", n - 1L, ", but is ", k, ". A ",
         point_type, " is not its own neighbour.", call. = FALSE)
  }

  if (is.null(eps)) {
    eps <- max(1L, k %/% 2L)
  } else if (eps > k) {
    stop("`eps` must not exceed `k`, which is ", k, ", but is ", eps,
         ". Two points cannot share more neighbours than they each keep. ",
         "Note that this `eps` counts shared neighbours; the one in ",
         "`cluster_dbscan()` is a radius.", call. = FALSE)
  }

  if (is.null(min_pts)) {
    min_pts <- max(2L, k %/% 2L)
  }

  list(k = k, eps = eps, min_pts = min_pts)
}


#' How far a neighbourhood reaches, when the caller has not said
#'
#' DBSCAN's `eps` is the one argument it cannot be given a fixed default for: it is
#' a radius in the units of the data, so any constant would be wrong on the next
#' matrix. It has to be derived from the distances actually present, and what it is
#' derived from is the k-distance curve — every point's distance to its
#' `min_pts - 1`th neighbour — which is the same quantity the manual procedure
#' plots with [dbscan::kNNdistplot()] and reads a knee off by eye.
#'
#' The value taken from it is the 95th percentile rather than the knee. Reading the
#' knee arithmetically, as the point of the sorted curve furthest below the chord
#' joining its ends, is the textbook translation of the manual procedure and it
#' does not survive contact with a curve that rises gradually: on the two-blob
#' probe in `Test/cursor_test/2026_08_15/` it lands at 0.55 where the answer is
#' anywhere from 0.9 up, and splits two clean blobs into three clusters and
#' twenty-eight points of noise.
#'
#' A percentile has no such failure mode and says something the knee never did.
#' `eps` at the 95th percentile is the radius that reaches the `min_pts - 1`th
#' neighbour of all but one point in twenty, so the rule reads as "assume about 5%
#' of the points are noise and set the radius accordingly" — which is a statement a
#' caller can disagree with, rather than a number off a curve. It is still a
#' heuristic, and it is reported as one: `parameters$eps_source` says the value was
#' derived, and a caller who has a radius in mind should pass it.
#'
#' @param m The point-by-variable matrix.
#' @param min_pts The density threshold, whose curve is read.
#'
#' @return The derived `eps`, always positive.
#'
#' @keywords internal
#' @noRd
sa_cluster_eps <- function(m, min_pts) {
  y <- dbscan::kNNdist(m, k = min_pts - 1L)
  eps <- stats::quantile(y[is.finite(y)], 0.95, names = FALSE)

  # Duplicated points give a curve that never leaves zero, and a zero radius would
  # make every point noise. The largest distance there is at least does something.
  if (!is.finite(eps) || eps <= 0) {
    positive <- y[is.finite(y) & y > 0]
    eps <- if (length(positive) > 0L) max(positive) else 1
  }
  eps
}
