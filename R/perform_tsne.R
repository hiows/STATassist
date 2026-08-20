# The stochastic counterpart of `perform_pca()`, and the reason there is more than
# one reduction in the package. A rotation can only draw straight structure. t-SNE
# keeps whichever points were neighbours in the full feature space close in two
# dimensions, so it finds structure that curves — and pays for it: it cannot say
# which feature made the picture, its global distances mean nothing, and a different
# perplexity is a different answer. Read beside a PCA of the same matrix, a cluster
# both of them find is a different fact from one only this method sees.
#
# `Rtsne`'s own preprocessing is turned off. `normalize = TRUE` would undo the
# `center` and `scale` this function was asked for, and `pca = TRUE` would show
# t-SNE the output of a principal component analysis rather than the matrix
# `perform_pca()` sees, which is the one thing that has to stay identical for the
# two to be read side by side. Both overrides are named in `engine$overridden`.
#
# `check_duplicates` is not overridden. Two samples identical across every feature
# make t-SNE refuse, and silently disabling a check the engine's author put there
# is a different thing from keeping the preprocessing we were asked for. `iris` is
# such a data set, which is worth knowing before trying it.

#' Embed samples or features with t-SNE
#'
#' Places each point in two or three dimensions so that the points it was near in
#' the full feature space stay near it, by t-distributed stochastic neighbour
#' embedding. What comes back is a picture and the coordinates to draw it with: an
#' axis is not a direction the way a principal component is, so there is nothing to
#' read off one but position.
#'
#' The input is the wide format the comparison functions take: **one row per sample
#' and one column per feature**, whichever margin is being embedded. What comes back
#' has one row per point, in one order every table follows, which is what makes an
#' embedding plottable against anything else read from the same frame.
#'
#' @details
#' # Scaling, and why it is on
#'
#' Every distance this method measures is a distance across all of `feats` at once,
#' so a feature measured in thousands would decide who is whose neighbour. `center`
#' and `scale` are therefore on by default, as they are for [perform_pca()], and the
#' two functions then see literally the same matrix — which is what lets the two
#' pictures be attributed to the methods rather than to the preprocessing.
#'
#' [perform_umap()] defaults the other way, since UMAP is more often run on
#' coordinates that already mean something, such as the components of a PCA.
#'
#' # Choosing the margin
#'
#' `embedding_scale = "features"` embeds the features instead, and unlike
#' [perform_pca()] this really does transpose: t-SNE embeds the rows it is handed
#' and has no second answer to read off the same fit. The features are standardised
#' first and the transpose is then embedded as it stands, which is what makes this
#' the same margin `perform_pca(embedding_scale = "features")` reports on.
#' Standardising after the transpose would standardise samples, which is what
#' `perform_tsne(t(data))` does and why it is not the same call.
#'
#' A feature margin needs enough features to be worth drawing. `perplexity` is read
#' off the number of points, so 8 features force a perplexity of 2 and a message
#' says so; at that size the loadings of a PCA are the whole answer and an embedding
#' has nothing to add. Sixty features derive a perplexity of 19, and on a simulation
#' with three planted correlation blocks the embedding then recovers them exactly.
#'
#' # What is dropped before anything runs
#'
#' `Rtsne` does not accept a missing value, so rows that are not complete and finite
#' across `feats` are dropped before it is called, and `design$n_dropped` reports how
#' many went. A feature that takes a single value cannot be scaled, so with
#' `scale = TRUE` it is left out with a message and named in `design$dropped_feats`.
#'
#' # Portability
#'
#' Everything but `$fit` is a data.frame, a character vector or a named list, so
#' dropping that one slot leaves an object that writes out as JSON. In a Python
#' transcription this is `sklearn.manifold.TSNE`.
#'
#' @param data A data.frame or a matrix in wide format, one row per sample and one
#'   column per feature. Row names are kept as the sample labels, repeated ones
#'   included; rows without a name are labelled by position.
#' @param feats Column names to embed, or `NULL` for every numeric column of `data`.
#'   A non-numeric column is left out with a message.
#' @param embedding_scale Which margin becomes the points of the picture:
#'   `"samples"`, the default, for one point per row of `data`, or `"features"` for
#'   one point per column, in which case the standardised matrix is transposed
#'   before it is embedded. `design$point_type` reports which it was.
#' @param center,scale Whether to centre each feature and divide it by its standard
#'   deviation before embedding. Both always apply to the **columns of `data`**,
#'   whatever `embedding_scale` is.
#' @param n_dim How many dimensions to embed into. `Rtsne` embeds into at most 3.
#' @param perplexity The neighbourhood size, or `NULL` to read one off the number of
#'   points. It is roughly how many neighbours each point is asked to keep close,
#'   and `Rtsne` requires `3 * perplexity <= n - 1`.
#' @param theta Barnes-Hut approximation. `0` is the exact gradient and slow, and
#'   larger values trade accuracy for speed.
#' @param seed Seed for the embedding, or `NULL` to use the stream as it stands.
#'   Supplying a seed does not disturb the caller: the previous random number state
#'   is put back when the function returns. Without one, two runs on the same data
#'   give two different pictures, which is a property of the method rather than a
#'   fault.
#'
#' @return An object of class `sa_reduction`, a plain list.
#'
#'   \describe{
#'     \item{`analysis`}{`"tsne"`.}
#'     \item{`points`}{Labels of the things that became points — samples or
#'       features, as `embedding_scale` asked — in the row order `scores` follows.}
#'     \item{`design`}{What was embedded: `point_type`, either `"sample"` or
#'       `"feature"`, the `feats` kept and any `dropped_feats`, and the counts
#'       `n_samples`, `n_used`, `n_dropped` and `n_feats`. The counts describe the
#'       input rather than the picture.}
#'     \item{`parameters`}{The choices as they were used rather than as they were
#'       passed, so a derived `perplexity` is the value that was derived.}
#'     \item{`scores`}{Coordinates: `points` beside `tSNE1`, `tSNE2` and so on.}
#'     \item{`engine`}{What computed the embedding, and which of its defaults were
#'       overridden.}
#'     \item{`fit`}{The `Rtsne` object. This is the slot that is not portable;
#'       dropping it leaves an object that writes out as JSON.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [perform_pca()], which answers about the same points in a way that can
#'   say which features moved them, and [perform_umap()], the other embedding.
#'
#' @examples
#' ## The group was never shown to the embedding, so the split is its to find.
#' sim <- simulate_two_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 3)
#' res <- perform_tsne(sim$args$data, seed = 1)
#' head(res$scores)
#'
#' \donttest{
#' plot(res$scores[c("tSNE1", "tSNE2")],
#'      col = as.integer(factor(sim$args$group)), pch = 16)
#' pca <- perform_pca(sim$args$data)
#' table(group = sim$args$group, side = pca$scores$PC1 > 0)
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
#' by_feat <- perform_tsne(blocks, feats = paste0("x_", 1:60),
#'                         embedding_scale = "features", seed = 1)
#' plot(by_feat$scores[c("tSNE1", "tSNE2")],
#'      col = rep(1:3, each = 20), pch = 16)
#' }
#'
#' @export
perform_tsne <- function(data,
                         feats = NULL,
                         embedding_scale = c("samples", "features"),
                         center = TRUE,
                         scale = TRUE,
                         n_dim = 2,
                         perplexity = NULL,
                         theta = 0.5,
                         seed = NULL) {

  embedding_scale <- match.arg(embedding_scale)
  sa_check_flag(center, "center")
  sa_check_flag(scale, "scale")
  n_dim <- sa_check_count(n_dim, "n_dim", 1)
  sa_check_scalar_num(theta, "theta", 0, 1)
  if (n_dim > 3L) {
    stop("`n_dim` must be at most 3, but is ", n_dim, ". Rtsne embeds into 1, 2 ",
         "or 3 dimensions; `perform_umap()` embeds into more.", call. = FALSE)
  }

  input <- sa_reduce_input(data, feats, scale, "perform_tsne")
  x <- input$x
  pt <- sa_reduce_points(x, input$samples, embedding_scale)
  xe <- sa_reduce_embedding_matrix(x, embedding_scale, center, scale)
  n_points <- nrow(xe)

  # Read before the engine is called, so that a rejected neighbourhood is rejected
  # for what it is rather than as whatever the engine makes of it. It is counted in
  # points: on the feature margin it is the features that have neighbours.
  perplexity <- sa_tsne_perplexity(perplexity, n_points, pt$point_type)
  sa_reduce_few_points(n_points, pt$point_type,
                       paste0("perplexity = ", perplexity))

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  fit <- Rtsne::Rtsne(xe, dims = n_dim, perplexity = perplexity, theta = theta,
                      normalize = FALSE, pca = FALSE, verbose = FALSE)

  sa_new_reduction(
    analysis = "tsne",
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
      perplexity      = perplexity,
      theta           = theta,
      seed            = seed
    ),
    scores = sa_embedding_frame(fit$Y, pt$points, "tSNE"),
    engine = list(package = "Rtsne", method = "Rtsne",
                  label = "t-distributed stochastic neighbour embedding",
                  overridden = c("normalize = FALSE", "pca = FALSE")),
    fit    = fit
  )
}
