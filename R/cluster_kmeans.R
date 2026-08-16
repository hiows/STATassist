# The one clustering here that is not deterministic, and the one whose default
# needs the most defending. k-means starts from a random set of centres and
# improves them until they stop moving, which finds a local optimum and not the
# best one; running it once and reporting the answer would make the result depend
# on the state of the random number generator when the caller happened to call it.
# `n_start = 25` is the answer to that — twenty-five starts, best one kept — and it
# is `stats::kmeans()`'s own default of 1 overridden, so it is declared in
# `engine$overridden` the way `perform_tsne()` declares Rtsne's.
#
# The distance is not an argument. k-means does not minimise a distance the caller
# chooses; it minimises the sum of squared Euclidean distances to the centres, and
# a mean is only the centre of its group under that one. Offering `"correlation"`
# here would be offering an algorithm that does not exist.

#' Cluster by moving centres until they stop
#'
#' Places `n_clust` centres, gives every point to its nearest one, moves each
#' centre to the mean of the points it was given, and repeats until nothing moves.
#' What comes back is one cluster label per point, beside the [stats::kmeans()]
#' object holding the centres.
#'
#' The input is the wide format the comparison functions take: **one row per sample
#' and one column per feature**. Which margin of it becomes the thing being
#' clustered is `cluster_scale`, and `design$point_type` reports the answer.
#'
#' @details
#' # Reproducibility
#'
#' The starting centres are random, so two calls can disagree. Two defences are on
#' by default and neither replaces the other. `n_start = 25` runs the whole thing
#' twenty-five times from different starts and keeps the one with the smallest
#' within-cluster sum of squares, which makes a bad local optimum unlikely rather
#' than impossible. `seed` makes the run exactly repeatable, and it restores the
#' random stream afterwards, so seeding this call does not quietly reseed whatever
#' the caller does next.
#'
#' Neither is a guarantee that `n_clust` is the right number. k-means always
#' returns the number of groups it was asked for, and `clusters$silhouette` is
#' where to look for whether they are groups: a cluster near zero is one whose
#' members are about as close to another cluster as to their own.
#'
#' # Which margin is clustered
#'
#' `cluster_scale = "samples"`, the default, puts one point per row of `data`, and
#' `"features"` puts one point per column. See [cluster_hclust()], which documents
#' the same argument and why transposing `data` by hand is a different analysis.
#'
#' # What is dropped before anything runs
#'
#' Rows that are not complete and finite across `feats` go before the centres are
#' placed, and `design$n_dropped` reports how many. A feature that takes a single
#' value cannot be scaled, so with `scale = TRUE` it is left out with a message and
#' named in `design$dropped_feats`.
#'
#' @param data A data.frame or a matrix in wide format, one row per sample and one
#'   column per feature.
#' @param feats Column names to cluster on, or `NULL` for every numeric column of
#'   `data`.
#' @param cluster_scale Which margin becomes the points being clustered:
#'   `"samples"`, the default, or `"features"`. `design$point_type` reports it.
#' @param center,scale Whether to centre each feature and divide it by its standard
#'   deviation first. Both always apply to the **columns of `data`**, whatever
#'   `cluster_scale` is. Scaling matters more here than anywhere else in this
#'   family, since a mean is taken in the units the features arrived in.
#' @param n_clust How many centres to place. Unlike the two density methods, this
#'   is the number of clusters that will come back.
#' @param n_start How many random starts to try, keeping the best. The default of
#'   25 overrides [stats::kmeans()]'s own default of 1.
#' @param iter_max Most iterations one start may take before it is abandoned.
#' @param seed Seed for the starts, or `NULL` to leave the random stream alone. The
#'   stream is restored afterwards either way.
#'
#' @return An object of class `sa_cluster`, a plain list. The slots are
#'   [cluster_hclust()]'s, with `analysis` `"kmeans"`, `fit` the
#'   [stats::kmeans()] object and `design$n_noise` always 0, since every point is
#'   given to a centre.
#'
#' @seealso [cluster_hclust()], which is told the same `n_clust` and returns the
#'   tree it cut, and [cluster_dbscan()], which derives the count from the density
#'   instead and may refuse to place a point.
#'
#' @examples
#' ## Three species, and the clustering was not told there were three.
#' res <- cluster_kmeans(iris[1:4], n_clust = 3, seed = 1)
#' res
#' table(cluster = res$assignments$cluster, species = iris$Species)
#'
#' ## The centres are on `$fit`, in the scaled units the clustering ran in.
#' res$fit$centers
#'
#' ## On a planted two-group structure. The group was never shown to the engine.
#' sim <- simulate_two_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 3)
#' by_samp <- cluster_kmeans(sim$args$data, n_clust = 2, seed = 1)
#' table(cluster = by_samp$assignments$cluster, group = sim$args$group)
#'
#' ## Asking for more groups than there are shows up in the silhouettes rather
#' ## than in an error: k-means returns whatever number it was asked for.
#' cluster_kmeans(sim$args$data, n_clust = 5, seed = 1)$clusters
#'
#' @export
cluster_kmeans <- function(data,
                           feats = NULL,
                           cluster_scale = c("samples", "features"),
                           center = TRUE,
                           scale = TRUE,
                           n_clust = 2,
                           n_start = 25,
                           iter_max = 100,
                           seed = NULL) {

  cluster_scale <- match.arg(cluster_scale)
  sa_check_flag(center, "center")
  sa_check_flag(scale, "scale")
  n_start <- sa_check_count(n_start, "n_start", 1)
  iter_max <- sa_check_count(iter_max, "iter_max", 1)

  input <- sa_cluster_input(data, feats, cluster_scale, center, scale,
                            "cluster_kmeans")
  m <- input$m
  n_clust <- sa_cluster_n_clust(n_clust, nrow(m), input$point_type)

  # `stats::kmeans()` refuses more centres than distinct points, since two centres
  # would have to start in the same place. Said here in terms of the argument the
  # caller passed rather than left to the engine's own wording.
  n_distinct <- nrow(unique(m))
  if (n_clust > n_distinct) {
    stop("`n_clust` is ", n_clust, " but only ", n_distinct, " of the ",
         nrow(m), " ", input$point_type, "(s) being clustered are distinct. ",
         "A centre cannot be placed where another one already is.",
         call. = FALSE)
  }

  restore <- sa_preserve_seed(seed)
  on.exit(restore(), add = TRUE)

  fit <- stats::kmeans(m, centers = n_clust, nstart = n_start,
                       iter.max = iter_max)

  d <- sa_cluster_dist(m, "euclidean")
  tables <- sa_cluster_tables(fit$cluster, input$points, d)

  sa_new_cluster(
    analysis = "kmeans",
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
      n_clust       = n_clust,
      n_start       = n_start,
      iter_max      = iter_max,
      seed          = seed,
      # Not an argument, but the silhouettes were measured on it and the print
      # method says which distance it is quoting.
      dist_method   = "euclidean",
      tot_withinss  = fit$tot.withinss
    ),
    assignments = tables$assignments,
    clusters    = tables$clusters,
    engine      = list(package = "stats", method = "kmeans",
                       label = "k-means clustering",
                       overridden = "nstart = 25"),
    fit         = fit
  )
}
