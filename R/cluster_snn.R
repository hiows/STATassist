# The second density method. Each point keeps its `k` nearest neighbours, two
# points are linked when their neighbour lists overlap, and clusters are the dense
# parts of that graph. A k nearest neighbour classifier is supervised — it labels a
# new point by what its neighbours were labelled — and that is not what this is.
#
# Which is worth having beside `cluster_dbscan()` rather than instead of it,
# because the two disagree about what "close" is in exactly the place a radius
# breaks down. DBSCAN measures one radius everywhere, so it cannot find two
# clusters of different densities at once: an `eps` large enough for the sparse
# one swallows the dense one. Sharing neighbours is a relative measure — the
# neighbours of a point in a sparse region are far away, but they are still its
# neighbours — so this finds both. That is also its known weakness: in high
# dimensions everything's neighbour lists start to overlap, and it will happily
# report structure that is only the curse of dimensionality. Neither is right, and
# a grouping both of them find is worth more than one only one of them does.

#' Cluster by how many neighbours points have in common
#'
#' Builds a graph in which every point keeps its `k` nearest neighbours, links two
#' points that share at least `eps` of them, and grows a cluster out of every point
#' with at least `min_pts` such links. Points that never joined one are left as
#' noise. How many clusters there are is the answer rather than the question.
#'
#' The input is the wide format the comparison functions take: **one row per sample
#' and one column per feature**. Which margin of it becomes the thing being
#' clustered is `cluster_scale`, and `design$point_type` reports the answer.
#'
#' @details
#' # `eps` here counts neighbours, and in cluster_dbscan() it is a distance
#'
#' The two functions have an argument of the same name meaning two different
#' things, which is the algorithms' doing rather than this package's:
#' [dbscan::dbscan()] and [dbscan::sNNclust()] both call their threshold `eps`, and
#' renaming one of them here would make this package's documentation disagree with
#' the engine's. In [cluster_dbscan()] `eps` is a radius in the units of the data
#' and any positive number is meaningful. Here it is a count of shared neighbours,
#' so it is a whole number between 1 and `k`, and a value of `k / 2` says "half of
#' what each of you considers close is the same points".
#'
#' # Why this rather than a radius
#'
#' One radius has to be right everywhere, so DBSCAN cannot find a tight cluster and
#' a loose one in the same call: the `eps` that reaches across the loose one merges
#' the tight one into its surroundings. Shared neighbours are relative, so a point
#' in a sparse region is still close to its own neighbours, and clusters of
#' different densities come out together.
#'
#' The price is dimension. As the number of variables grows, distances concentrate
#' and neighbour lists start to overlap for no reason, so this method will report
#' structure that is an artefact of the width of the table. `clusters$silhouette`
#' is the check worth making, and agreement with [cluster_dbscan()] on the same
#' points is worth more than either alone.
#'
#' # Noise
#'
#' A point that never joined a dense part of the graph gets cluster `0`, has no
#' silhouette and no row in `clusters`; `design$n_noise` is the count. Every point
#' being noise is a possible answer and means the overlap asked for is not there,
#' which is usually `eps` or `min_pts` being too large for the `k` in use.
#'
#' @param data A data.frame or a matrix in wide format, one row per sample and one
#'   column per feature.
#' @param feats Column names to cluster on, or `NULL` for every numeric column of
#'   `data`.
#' @param cluster_scale Which margin becomes the points being clustered:
#'   `"samples"`, the default, or `"features"`. `design$point_type` reports it.
#' @param center,scale Whether to centre each feature and divide it by its standard
#'   deviation first. Both always apply to the **columns of `data`**, whatever
#'   `cluster_scale` is.
#' @param k How many nearest neighbours each point keeps. `NULL` derives it as the
#'   square root of the number of points and reports what it came out as.
#' @param eps How many neighbours two points must **share** before they are linked,
#'   a whole number from 1 to `k`. Not a distance; see the details. `NULL` uses
#'   `k / 2`.
#' @param min_pts How many links a point needs in that graph before it can be the
#'   core of a cluster. `NULL` uses `k / 2`.
#'
#' @return An object of class `sa_cluster`, a plain list. The slots are
#'   [cluster_hclust()]'s, with `analysis` `"snn"`, `fit` the
#'   [dbscan::sNNclust()] object, and `design$n_clusters` and `design$n_noise` both
#'   answers rather than arguments.
#'
#' @seealso [cluster_dbscan()], the other density method here, which measures
#'   closeness with one radius rather than by shared neighbours.
#'
#' @examples
#' ## Two blobs of different spread, which is the case a single radius handles
#' ## badly and shared neighbours handle well.
#' set.seed(1)
#' blobs <- rbind(
#'   matrix(rnorm(80, mean = -6, sd = 0.4), ncol = 2),
#'   matrix(rnorm(80, mean =  6, sd = 2.5), ncol = 2)
#' )
#' colnames(blobs) <- c("x", "y")
#' res <- cluster_snn(blobs, scale = FALSE, k = 10)
#' res
#' table(res$assignments$cluster)
#'
#' plot(blobs, col = res$assignments$cluster + 1L, pch = 16,
#'      main = "cluster 0, in black, is noise")
#'
#' ## On a planted two-group structure, with everything left to be derived.
#' sim <- simulate_two_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 3)
#' by_samp <- cluster_snn(sim$args$data)
#' table(cluster = by_samp$assignments$cluster, group = sim$args$group)
#'
#' @export
cluster_snn <- function(data,
                        feats = NULL,
                        cluster_scale = c("samples", "features"),
                        center = TRUE,
                        scale = TRUE,
                        k = NULL,
                        eps = NULL,
                        min_pts = NULL) {

  cluster_scale <- match.arg(cluster_scale)
  sa_check_flag(center, "center")
  sa_check_flag(scale, "scale")

  input <- sa_cluster_input(data, feats, cluster_scale, center, scale,
                            "cluster_snn")
  m <- input$m
  par <- sa_snn_params(k, eps, min_pts, nrow(m), input$point_type)

  fit <- dbscan::sNNclust(m, k = par$k, eps = par$eps, minPts = par$min_pts)

  d <- sa_cluster_dist(m, "euclidean")
  tables <- sa_cluster_tables(fit$cluster, input$points, d)

  if (tables$n_clusters == 0L) {
    message("No ", input$point_type, " shares eps = ", par$eps,
            " of its k = ", par$k, " neighbour(s) with min_pts = ",
            par$min_pts, " others, so every one of them is noise. A smaller ",
            "`eps` or `min_pts`, or a larger `k`, is what asks for less ",
            "overlap.")
  }

  sa_new_cluster(
    analysis = "snn",
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
      k             = par$k,
      eps           = par$eps,
      min_pts       = par$min_pts,
      dist_method   = "euclidean"
    ),
    assignments = tables$assignments,
    clusters    = tables$clusters,
    engine      = list(package = "dbscan", method = "sNNclust",
                       label = "Shared nearest neighbour clustering",
                       overridden = character(0)),
    fit         = fit
  )
}
