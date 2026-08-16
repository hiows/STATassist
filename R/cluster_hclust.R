# The clustering that keeps its working. The other three return a labelling and
# nothing that says how it was arrived at; this one returns the tree, and the tree
# is a different kind of answer. It holds every cut the data admits, not just the
# one that was asked for, so `n_clust` is a question put to a fitted object rather
# than a parameter of the fit, and `stats::cutree($fit, k = 5)` asks it again
# without recomputing anything.
#
# It is also the only one of the four with a choice of distance. k-means minimises
# squared Euclidean distance by construction and the two density methods measure a
# radius in it, so for those three the distance is the method's and not the
# caller's. Here it is genuinely open, and `"correlation"` is the reason the choice
# is offered: it groups points by the shape of their profile rather than by how
# high it sits, which on a feature scale is the difference between "these genes
# move together" and "these genes are expressed at similar levels".

#' Cluster by building a tree and cutting it
#'
#' Merges the two closest points, then the two closest groups, and keeps going
#' until everything is one group; then cuts the resulting tree so that `n_clust`
#' groups fall out. What comes back is one cluster label per point, beside the
#' [stats::hclust()] tree the labels were cut from.
#'
#' The input is the wide format the comparison functions take: **one row per sample
#' and one column per feature**. Which margin of it becomes the thing being
#' clustered is `cluster_scale`, and `design$point_type` reports the answer.
#'
#' @details
#' # Which margin is clustered
#'
#' `cluster_scale = "samples"`, the default, puts one point per row of `data`, and
#' the question is which samples resemble each other. `"features"` puts one point
#' per column, and the question is which features move together. The second is the
#' one this and [draw_heatmap()] agree on: the features are standardised first and
#' the transpose is clustered as it stands, which is also what
#' `perform_pca(embedding_scale = "features")` reports on.
#'
#' Transposing `data` by hand instead is a third analysis. [scale()] standardises
#' **columns**, so on the transpose it standardises samples, and the picture looks
#' right while answering a different question. `cluster_scale` exists so that the
#' transpose never has to be taken by hand.
#'
#' # Choosing the distance and the linkage
#'
#' `dist_method` decides what "close" means and `hclust_method` decides what
#' "close" means for two groups rather than two points. `"correlation"` is
#' `1 - cor()`, so two points a constant apart are at distance zero: use it when
#' the shape of a profile is the question and its level is not. The default
#' `"average"` linkage is the one that neither chains clusters into strings the way
#' single linkage does nor insists they come out round the way `"ward.D2"` does.
#'
#' The tree is built on exactly the distance [draw_heatmap()] builds its
#' dendrogram from, so a heatmap drawn with the same `dist_method` and
#' `hclust_method` is showing this clustering and not a near relative of it.
#'
#' # What is dropped before anything runs
#'
#' Rows that are not complete and finite across `feats` go before the distance is
#' measured, and `design$n_dropped` reports how many. This is the listwise deletion
#' the rest of the package uses; nothing is imputed. A feature that takes a single
#' value cannot be scaled, so with `scale = TRUE` it is left out with a message and
#' named in `design$dropped_feats`.
#'
#' @param data A data.frame or a matrix in wide format, one row per sample and one
#'   column per feature. This is the same layout [compare_two_groups()] and
#'   [draw_heatmap()] take. Row names are kept as the sample labels, repeated ones
#'   included; rows without a name are labelled by position.
#' @param feats Column names to cluster on, or `NULL` for every numeric column of
#'   `data`. A non-numeric column is left out with a message, so a frame that
#'   carries a grouping column alongside the measurements can be passed as it is.
#' @param cluster_scale Which margin becomes the points being clustered:
#'   `"samples"`, the default, for one point per row of `data`, or `"features"` for
#'   one point per column. `design$point_type` reports it. See the details: this is
#'   not the same as transposing `data` yourself.
#' @param center,scale Whether to centre each feature and divide it by its standard
#'   deviation before measuring anything. Scaling is on by default because features
#'   are not measured on a common scale, and without it the feature with the widest
#'   units decides which points are close. Both always apply to the **columns of
#'   `data`**, whatever `cluster_scale` is.
#' @param n_clust How many groups to cut the tree into. The tree itself does not
#'   depend on this, so `$fit` can be cut again at another `k`.
#' @param dist_method What "close" means. `"correlation"` is `1 - cor()`, which
#'   compares the shape of a profile rather than its level; the other two are
#'   handed to [stats::dist()].
#' @param hclust_method Linkage handed to [stats::hclust()], which is what "close"
#'   means for two groups rather than two points.
#'
#' @return An object of class `sa_cluster`, a plain list.
#'
#'   \describe{
#'     \item{`analysis`}{`"hclust"`.}
#'     \item{`points`}{Labels of the things that were clustered — samples or
#'       features, as `cluster_scale` asked — in the row order `assignments`
#'       follows.}
#'     \item{`design`}{What was clustered: `point_type`, the `feats` kept and any
#'       `dropped_feats`, the counts `n_samples`, `n_used`, `n_dropped` and
#'       `n_feats`, and the two counts of the answer, `n_clusters` and `n_noise`.
#'       `n_noise` is always 0 here: a tree places every point.}
#'     \item{`parameters`}{The choices as they were used.}
#'     \item{`assignments`}{One row per point: `points`, `cluster` and
#'       `silhouette`.}
#'     \item{`clusters`}{One row per cluster: `cluster`, `size` and the mean
#'       `silhouette` of its members.}
#'     \item{`engine`}{What computed the clustering.}
#'     \item{`fit`}{The [stats::hclust()] object. This is the slot that is not
#'       portable; dropping it leaves an object that writes out as JSON.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [cluster_kmeans()], which is told the same `n_clust` and answers
#'   without a tree, and [cluster_dbscan()], which is not told it at all.
#'   [draw_heatmap()] draws the tree this cuts.
#'
#' @examples
#' ## Three species, and the clustering was not told there were three.
#' res <- cluster_hclust(iris[1:4], n_clust = 3)
#' res
#' table(cluster = res$assignments$cluster, species = iris$Species)
#'
#' ## The tree is `$fit`, so another cut costs nothing.
#' table(stats::cutree(res$fit, k = 2))
#'
#' ## On a planted two-group structure, clustering the samples.
#' sim <- simulate_two_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 3)
#' by_samp <- cluster_hclust(sim$args$data, n_clust = 2)
#' table(cluster = by_samp$assignments$cluster, group = sim$args$group)
#'
#' ## `cluster_scale = "features"` asks the other question: which features move
#' ## together. The correlation distance is the one that ignores their levels.
#' by_feat <- cluster_hclust(sim$args$data, cluster_scale = "features",
#'                           n_clust = 3, dist_method = "correlation")
#' by_feat
#' split(by_feat$points, by_feat$assignments$cluster)
#'
#' @export
cluster_hclust <- function(data,
                           feats = NULL,
                           cluster_scale = c("samples", "features"),
                           center = TRUE,
                           scale = TRUE,
                           n_clust = 2,
                           dist_method = c("euclidean", "manhattan",
                                           "correlation"),
                           hclust_method = c("average", "complete", "ward.D2")) {

  cluster_scale <- match.arg(cluster_scale)
  dist_method <- match.arg(dist_method)
  hclust_method <- match.arg(hclust_method)
  sa_check_flag(center, "center")
  sa_check_flag(scale, "scale")

  input <- sa_cluster_input(data, feats, cluster_scale, center, scale,
                            "cluster_hclust")
  n_clust <- sa_cluster_n_clust(n_clust, nrow(input$m), input$point_type)

  d <- sa_cluster_dist(input$m, dist_method)
  if (!all(is.finite(d))) {
    # `draw_heatmap()` falls back to the input order here, because the picture is
    # still a picture without a tree. There is no such fallback when the tree is
    # the whole answer.
    stop("`cluster_hclust()` cannot measure the ", dist_method,
         " distance between every pair of ", input$point_type,
         "s: some of them are undefined, which happens when a ",
         input$point_type, " has no variance across the other margin. ",
         "A different `dist_method` is not governed by this.", call. = FALSE)
  }

  fit <- stats::hclust(d, method = hclust_method)
  tables <- sa_cluster_tables(stats::cutree(fit, k = n_clust), input$points, d)

  sa_new_cluster(
    analysis = "hclust",
    points   = input$points,
    # `design` describes the input, so its counts do not turn with
    # `cluster_scale`: `n_samples` is always rows of `data` and `feats` always the
    # columns kept. `point_type` is what says which of the two became the points.
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
      n_clust       = n_clust,
      dist_method   = dist_method,
      hclust_method = hclust_method
    ),
    assignments = tables$assignments,
    clusters    = tables$clusters,
    engine      = list(package = "stats", method = "hclust",
                       label = "Hierarchical clustering",
                       overridden = character(0)),
    fit         = fit
  )
}
