# The first of the two that are not told how many groups to find, and the first
# that is allowed to answer "this point is not in one". Both of those follow from
# the same change of question. `cluster_hclust()` and `cluster_kmeans()` are asked
# to divide the points into `n_clust` parts, and a division has no room for a
# leftover; DBSCAN is asked which regions are dense, and a point outside all of
# them is not a small cluster of one, it is noise.
#
# What it costs is that the two arguments it does take are harder to choose than a
# cluster count. `min_pts` is a number of points and can be given a rule; `eps` is
# a radius in the units of the data and cannot, so a constant default would be
# wrong on the next matrix. Rather than refuse to run without it, `eps` is derived
# from the neighbour distances the matrix actually has and reported as derived in
# `parameters$eps_source`. A caller with a radius in mind should pass one; a caller
# without one gets a defensible guess and is told it was a guess.

#' Cluster by finding the dense regions
#'
#' Grows a cluster out of every point that has at least `min_pts` neighbours within
#' `eps` of it, joining the neighbourhoods that overlap, and leaves the points that
#' never fell into one as noise. How many clusters there are is the answer rather
#' than the question, and it can be none.
#'
#' The input is the wide format the comparison functions take: **one row per sample
#' and one column per feature**. Which margin of it becomes the thing being
#' clustered is `cluster_scale`, and `design$point_type` reports the answer.
#'
#' @details
#' # Noise is a result, not a failure
#'
#' A point that never joined a neighbourhood gets cluster `0`, which is
#' [dbscan::dbscan()]'s convention and this contract's. It has no silhouette, since
#' noise is not a cluster to be near or far from, and `clusters` has no row for it;
#' `design$n_noise` is the count. All points being noise is a possible answer and
#' means the density asked for is not present, which is usually `eps` being too
#' small or `min_pts` too large for the data.
#'
#' # Choosing eps
#'
#' `eps` is a radius in the units of the matrix being clustered, which is why there
#' is no default constant that could be right. Left as `NULL` it is derived from
#' the distances actually present: every point's distance to its `min_pts - 1`th
#' neighbour is collected, and `eps` is the 95th percentile of those. The rule
#' therefore reads as "assume about one point in twenty is noise, and set the
#' radius that leaves that many outside", which is a statement worth disagreeing
#' with rather than a number off a curve.
#'
#' The quantity being read is the one [dbscan::kNNdistplot()] plots and the manual
#' procedure reads a knee off by eye. The knee is not what is taken, because
#' locating it arithmetically fails on a curve that rises gradually rather than
#' turning; the probe in `Test/cursor_test/2026_08_15/` has the case. Either way it
#' is a heuristic: `parameters$eps_source` records whether the value was supplied
#' or derived, and a derived one is said out loud when it is used.
#'
#' Note that this `eps` is a distance. [cluster_snn()] has an argument of the same
#' name that counts shared neighbours, because the algorithm it comes from named it
#' that; the two are not comparable.
#'
#' # Which margin is clustered
#'
#' `cluster_scale = "samples"`, the default, puts one point per row of `data`, and
#' `"features"` puts one point per column. See [cluster_hclust()], which documents
#' the same argument and why transposing `data` by hand is a different analysis.
#'
#' @param data A data.frame or a matrix in wide format, one row per sample and one
#'   column per feature.
#' @param feats Column names to cluster on, or `NULL` for every numeric column of
#'   `data`.
#' @param cluster_scale Which margin becomes the points being clustered:
#'   `"samples"`, the default, or `"features"`. `design$point_type` reports it.
#' @param center,scale Whether to centre each feature and divide it by its standard
#'   deviation first. Both always apply to the **columns of `data`**, whatever
#'   `cluster_scale` is. Scaling matters here because `eps` is one radius applied to
#'   all of them at once.
#' @param eps The radius of a neighbourhood, in the units of the matrix being
#'   clustered, or `NULL` to derive it as the radius reaching the `min_pts - 1`th
#'   neighbour of 95% of the points. See the details.
#' @param min_pts How many points must be within `eps` of a point, itself included,
#'   before it can be the core of a cluster. `NULL` derives it from the number of
#'   variables, capped at half the points; the derived value is reported in a
#'   message.
#'
#' @return An object of class `sa_cluster`, a plain list. The slots are
#'   [cluster_hclust()]'s, with `analysis` `"dbscan"`, `fit` the
#'   [dbscan::dbscan()] object, and `design$n_clusters` and `design$n_noise` both
#'   answers rather than arguments.
#'
#' @seealso [cluster_snn()], the other density method here, which measures
#'   closeness by how many neighbours two points share rather than by a radius, and
#'   [cluster_kmeans()], which is told the number of groups and places every point.
#'
#' @examples
#' ## Two well-separated blobs and a couple of strays. The count is not an
#' ## argument, so finding two is the result rather than the setup.
#' set.seed(1)
#' blobs <- rbind(
#'   matrix(rnorm(60, mean = -3), ncol = 2),
#'   matrix(rnorm(60, mean =  3), ncol = 2),
#'   cbind(c(-12, 12), c(12, -12))
#' )
#' colnames(blobs) <- c("x", "y")
#' res <- cluster_dbscan(blobs, scale = FALSE)
#' res
#' table(res$assignments$cluster)
#'
#' plot(blobs, col = res$assignments$cluster + 1L, pch = 16,
#'      main = "cluster 0, in black, is noise")
#'
#' ## On a planted two-group structure, with the radius left to be derived.
#' sim <- simulate_two_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 3)
#' by_samp <- cluster_dbscan(sim$args$data)
#' by_samp$parameters$eps_source
#' table(cluster = by_samp$assignments$cluster, group = sim$args$group)
#'
#' @export
cluster_dbscan <- function(data,
                           feats = NULL,
                           cluster_scale = c("samples", "features"),
                           center = TRUE,
                           scale = TRUE,
                           eps = NULL,
                           min_pts = NULL) {

  cluster_scale <- match.arg(cluster_scale)
  sa_check_flag(center, "center")
  sa_check_flag(scale, "scale")
  # Before the input is read and before `min_pts` is derived, so that a call that
  # is going to fail does not first announce a default it never used.
  if (!is.null(eps)) {
    sa_check_scalar_num(eps, "eps", 0, lower_open = TRUE)
  }

  input <- sa_cluster_input(data, feats, cluster_scale, center, scale,
                            "cluster_dbscan")
  m <- input$m
  min_pts <- sa_dbscan_min_pts(min_pts, nrow(m), ncol(m), input$point_type)

  eps_source <- if (is.null(eps)) "derived" else "supplied"
  if (is.null(eps)) {
    eps <- sa_cluster_eps(m, min_pts)
    message("Using eps = ", format(signif(eps, 4), trim = TRUE),
            ", the radius reaching the ", min_pts - 1L,
            "th neighbour of 95% of the ", input$point_type,
            "(s). Pass `eps` to set it.")
  }

  fit <- dbscan::dbscan(m, eps = eps, minPts = min_pts)

  d <- sa_cluster_dist(m, "euclidean")
  tables <- sa_cluster_tables(fit$cluster, input$points, d)

  if (tables$n_clusters == 0L) {
    message("No ", input$point_type,
            " has min_pts = ", min_pts, " neighbour(s) within eps = ",
            format(signif(eps, 4), trim = TRUE),
            ", so every one of them is noise. A larger `eps` or a smaller ",
            "`min_pts` is what asks for less density.")
  }

  sa_new_cluster(
    analysis = "dbscan",
    points   = input$points,
    design   = list(
      point_type    = input$point_type,
      n_samples     = input$n_samples,
      n_used        = input$n_used,
      n_dropped     = input$n_dropped,
      n_feats       = input$n_feats,
      feats         = input$feats,
      dropped_feats = input$dropped_feats,
      n_clusters    = tables$n_clusters,
      n_noise       = tables$n_noise
    ),
    parameters = list(
      cluster_scale = cluster_scale,
      center        = center,
      scale         = scale,
      eps           = eps,
      eps_source    = eps_source,
      min_pts       = min_pts,
      dist_method   = "euclidean"
    ),
    assignments = tables$assignments,
    clusters    = tables$clusters,
    engine      = list(package = "dbscan", method = "dbscan",
                       label = "Density-based spatial clustering",
                       overridden = character(0)),
    fit         = fit
  )
}
