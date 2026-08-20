# The other embedding, and the one that differs from the two before it in what it
# does to the input first: nothing. `perform_pca()` and `perform_tsne()` standardise
# the features by default because a rotation and a Barnes-Hut gradient both answer
# in euclidean distance and have no other way to be told that a feature measured in
# thousands is not more important than one measured in units. UMAP is handed a
# `metric` instead, and two of its four choices answer that question themselves by
# comparing the shape of a row rather than its size. Standardising first would be
# answering it twice, so the default here is the engine's own and `scale = TRUE` is
# one argument away for the cases that need it.
#
# The default is the `umap` package's own R implementation, `method = "naive"`, so no
# Python installation is involved. `method = "umap-learn"` is the other choice the
# engine offers and needs `reticulate` and a Python environment with `umap-learn`.
#
# `umap` restores the random stream itself, which is worth knowing before comparing
# it with `perform_tsne()`: two seedless UMAP runs from the same state agree, while
# two seedless t-SNE runs do not, because each of those consumes the stream.

#' Embed samples or features with UMAP
#'
#' Places each point in a few dimensions so that the neighbourhood structure of the
#' full feature space survives, by uniform manifold approximation and projection.
#' What comes back is a picture and the coordinates to draw it with: an axis is not a
#' direction the way a principal component is, so there is nothing to read off one
#' but position.
#'
#' The input is the wide format the comparison functions take: **one row per sample
#' and one column per feature**, whichever margin is being embedded. What comes back
#' has one row per point, in one order every table follows, which is what makes an
#' embedding plottable against anything else read from the same frame.
#'
#' @details
#' # Scaling, and why it is off
#'
#' `center` and `scale` are `FALSE` here and `TRUE` in [perform_pca()] and
#' [perform_tsne()]. The difference is `metric`: `"cosine"` and `"pearson"` compare
#' the shape of a row rather than its size, so they have already answered the
#' question standardising would answer, and the `umap` package does not standardise
#' either. The default is therefore the engine's own picture of the data as it
#' arrived.
#'
#' What that leaves the caller is one decision rather than none. With
#' `metric = "euclidean"` or `"manhattan"` on features that are not measured on a
#' common scale, the feature with the widest units decides who is whose neighbour,
#' and `scale = TRUE` is what the picture needs — it is also what makes this
#' embedding comparable with a [perform_pca()] or [perform_tsne()] of the same data,
#' since those two standardise by default.
#'
#' One consequence of the default is that a feature of no variance is kept. It can
#' be kept because nothing divides by its standard deviation; it contributes nothing
#' to any distance either way. With `scale = TRUE` it is left out with a message and
#' named in `design$dropped_feats`.
#'
#' # Choosing the margin
#'
#' `embedding_scale = "features"` embeds the features instead, and unlike
#' [perform_pca()] this really does transpose: UMAP embeds the rows it is handed and
#' has no second answer to read off the same fit. `center` and `scale`, if they are
#' turned on, still apply to the **features** and the transpose is then embedded as
#' it stands, which is what makes this the same margin
#' `perform_pca(embedding_scale = "features")` reports on. `perform_umap(t(data))`
#' with `scale = TRUE` would standardise samples instead, which is a third analysis.
#'
#' A feature margin needs enough features to be worth drawing. `n_neighbors` is read
#' off the number of points, so 8 features force a neighbourhood of 8 and a message
#' says so; at that size the loadings of a PCA are the whole answer.
#'
#' # What is dropped before anything runs
#'
#' `umap` does not accept a missing value, so rows that are not complete and finite
#' across `feats` are dropped before it is called, and `design$n_dropped` reports how
#' many went. This is the listwise deletion the rest of the package uses; nothing is
#' imputed.
#'
#' # Portability
#'
#' Everything but `$fit` is a data.frame, a character vector or a named list, so
#' dropping that one slot leaves an object that writes out as JSON. The `umap`
#' object is the one genuinely unportable thing this package produces: its
#' `config$metric.function` is a function. In a Python transcription this is the
#' separate `umap-learn` package rather than anything in scikit-learn.
#'
#' @param data A data.frame or a matrix in wide format, one row per sample and one
#'   column per feature. Row names are kept as the sample labels, repeated ones
#'   included; rows without a name are labelled by position.
#' @param feats Column names to embed, or `NULL` for every numeric column of `data`.
#'   A non-numeric column is left out with a message.
#' @param embedding_scale Which margin becomes the points of the picture:
#'   `"samples"`, the default, for one point per row of `data`, or `"features"` for
#'   one point per column, in which case the matrix is transposed before it is
#'   embedded. `design$point_type` reports which it was.
#' @param center,scale Whether to centre each feature and divide it by its standard
#'   deviation before embedding. Both are off by default, unlike [perform_pca()] and
#'   [perform_tsne()]; see the details. They always apply to the **columns of
#'   `data`**, whatever `embedding_scale` is.
#' @param n_dim How many dimensions to embed into.
#' @param n_neighbors The neighbourhood size, or `NULL` to read one off the number of
#'   points. Small values follow local detail and large ones the overall shape.
#' @param min_dist How tightly points that belong together are allowed to be packed.
#' @param method Which [umap::umap()] backend to use: `"naive"`, the default, is the
#'   package's pure R implementation and needs no Python; `"umap-learn"` calls the
#'   reference implementation through `reticulate` and requires `umap-learn` installed
#'   in a Python environment the package can see.
#' @param metric Distance neighbours are measured with. `"cosine"` and `"pearson"`
#'   compare the shape of a row rather than its size, which is why the two scaling
#'   flags are off by default.
#' @param seed Seed for the embedding, or `NULL` to use the stream as it stands.
#'   Supplying a seed does not disturb the caller: the previous random number state
#'   is put back when the function returns. It buys less here than it does in
#'   [perform_tsne()], since the `umap` package restores the stream itself, so two
#'   seedless runs from the same state already agree.
#'
#' @return An object of class `sa_reduction`, a plain list.
#'
#'   \describe{
#'     \item{`analysis`}{`"umap"`.}
#'     \item{`points`}{Labels of the things that became points — samples or
#'       features, as `embedding_scale` asked — in the row order `scores` follows.}
#'     \item{`design`}{What was embedded: `point_type`, either `"sample"` or
#'       `"feature"`, the `feats` kept and any `dropped_feats`, and the counts
#'       `n_samples`, `n_used`, `n_dropped` and `n_feats`. The counts describe the
#'       input rather than the picture.}
#'     \item{`parameters`}{The choices as they were used rather than as they were
#'       passed, so a derived `n_neighbors` is the value that was derived.}
#'     \item{`scores`}{Coordinates: `points` beside `UMAP1`, `UMAP2` and so on.}
#'     \item{`engine`}{What computed the embedding, and which of its defaults were
#'       overridden.}
#'     \item{`fit`}{The `umap` object. This is the slot that is not portable;
#'       dropping it leaves an object that writes out as JSON.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [perform_pca()], which answers about the same points in a way that can
#'   say which features moved them, and [perform_tsne()], the other embedding.
#'
#' @examples
#' ## The group was never shown to the embedding, so the split is its to find.
#' sim <- simulate_two_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 3)
#' res <- perform_umap(sim$args$data, seed = 1)
#' head(res$scores)
#'
#' \donttest{
#' plot(res$scores[c("UMAP1", "UMAP2")],
#'      col = as.integer(factor(sim$args$group)), pch = 16)
#' scaled <- perform_umap(sim$args$data, scale = TRUE, seed = 1)
#'
#' ## Feature-scale embedding over planted correlation blocks.
#' cor_mat <- make_block_cor(
#'   n_features = 60,
#'   blocks = list(list(features = 1:20, cor = 0.8),
#'                 list(features = 21:40, cor = 0.6),
#'                 list(features = 41:60, cor = 0.4))
#' )
#' blocks <- simulate_classification(n_pred = 60, cor_mat = cor_mat,
#'                                   seed = 2026)$args$data
#' by_feat <- perform_umap(blocks, feats = paste0("x_", 1:60),
#'                         embedding_scale = "features", scale = TRUE, seed = 1)
#' plot(by_feat$scores[c("UMAP1", "UMAP2")],
#'      col = rep(1:3, each = 20), pch = 16)
#' }
#' @export
perform_umap <- function(data,
                         feats = NULL,
                         embedding_scale = c("samples", "features"),
                         center = FALSE,
                         scale = FALSE,
                         n_dim = 2,
                         n_neighbors = NULL,
                         min_dist = 0.1,
                         method = c("naive", "umap-learn"),
                         metric = c("euclidean", "manhattan", "cosine",
                                    "pearson"),
                         seed = NULL) {

  embedding_scale <- match.arg(embedding_scale)
  method <- match.arg(method)
  metric <- match.arg(metric)
  sa_check_flag(center, "center")
  sa_check_flag(scale, "scale")
  n_dim <- sa_check_count(n_dim, "n_dim", 1)
  sa_check_scalar_num(min_dist, "min_dist", 0, lower_open = TRUE)

  input <- sa_reduce_input(data, feats, scale, "perform_umap")
  x <- input$x
  pt <- sa_reduce_points(x, input$samples, embedding_scale)
  xe <- sa_reduce_embedding_matrix(x, embedding_scale, center, scale)
  n_points <- nrow(xe)

  # Read before the engine is called, so that a rejected neighbourhood is rejected
  # for what it is rather than as whatever the engine makes of it. It is counted in
  # points: on the feature margin it is the features that have neighbours.
  n_neighbors <- sa_umap_neighbors(n_neighbors, n_points, pt$point_type)
  sa_reduce_few_points(n_points, pt$point_type,
                       paste0("n_neighbors = ", n_neighbors))

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  fit <- umap::umap(xe, method = method, n_neighbors = n_neighbors,
                    n_components = n_dim, min_dist = min_dist, metric = metric)

  sa_new_reduction(
    analysis = "umap",
    points   = pt$points,
    design   = list(
      point_type    = pt$point_type,
      n_samples     = input$n_samples,
      n_used        = nrow(x),
      n_dropped     = input$n_dropped,
      n_feats       = ncol(x),
      feats         = colnames(x),
      dropped_feats = input$dropped_feats
    ),
    parameters = list(
      embedding_scale = embedding_scale,
      center          = center,
      scale           = scale,
      n_dim           = n_dim,
      n_neighbors     = n_neighbors,
      min_dist        = min_dist,
      method          = method,
      metric          = metric,
      seed            = seed
    ),
    scores = sa_embedding_frame(fit$layout, pt$points, "UMAP"),
    engine = list(package = "umap", method = "umap",
                  label = "Uniform manifold approximation and projection",
                  overridden = if (identical(method, "naive")) {
                    character(0)
                  } else {
                    paste0('method = "', method, '"')
                  }),
    fit    = fit
  )
}
