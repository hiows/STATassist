# The first unsupervised function of the machine learning family, and the first one
# with no outcome at all. There is nothing to predict and nothing to score, so the
# result is `sa_reduction` and its row axis is `points`, whichever margin the caller
# asked for.
#
# This is the only one of the three that answers on both margins from one fit. A
# principal component analysis is a singular value decomposition, so `$x` holds the
# samples and `$rotation` holds the features and the two are the same fit read from
# either end. That is why `embedding_scale` does not transpose anything here: the
# matrix goes in as it arrived and the argument only decides which side of the
# decomposition becomes `$scores`. `perform_tsne()` and `perform_umap()` have no
# such dual, and the transpose they are handed is what makes the three answer about
# the same margin.
#
# The feature scale reports `$rotation` rescaled from unit length to
# variance-weighted length, `%*% diag(sdev * sqrt(n - 1))`, rather than the raw
# loadings. Unit-length loadings are the same map on a scale nothing else in the
# result uses; rescaled, they are the coordinates a feature would have if the
# features had been the rows, and they can be read on the same axis proportions
# `$variance` reports.

#' Reduce many features to a few components
#'
#' Rotates a wide table of samples and features onto the axes that carry the most
#' variance, so that a few coordinates stand in for all of them. The components are
#' ordered by how much of the variance they carry, `$variance` says how much that
#' is, and `$loadings` says which features built each one.
#'
#' The input is the wide format the comparison functions take: **one row per sample
#' and one column per feature**. What comes back has one row per point, in one order
#' every table follows, which is what makes a reduction plottable against anything
#' else read from the same frame.
#'
#' @details
#' # Which slot holds what
#'
#' `$fit` is the [stats::prcomp()] object itself, always fitted to the samples,
#' since that is the fit both margins are read from. `$scores` is the margin
#' `embedding_scale` asked for and is what a scatter plot of the points is drawn
#' from. On the default sample scale the two agree: `$scores` is `$fit$x` beside a
#' `points` column, and `$loadings` is `$fit$rotation`.
#'
#' # Choosing the margin, rather than transposing
#'
#' `embedding_scale = "features"` is how a map of the features is asked for, and
#' hand-transposing the input is not the same thing. [stats::prcomp()] centres and
#' scales **columns**: on a sample-by-feature matrix that standardises features,
#' which is what this method needs, and on the transpose it standardises samples
#' instead. So `perform_pca(t(data))` runs a third analysis that is neither of the
#' two on offer. On the 8-feature simulation in the examples its first component
#' correlates 0.82 with the loadings, where `embedding_scale = "features"`
#' correlates 1.
#'
#' What the argument does instead is leave the matrix alone. One decomposition
#' already answers on both margins, so the features are the rows of `$rotation` and
#' nothing has to be recomputed. What `$scores` reports is those loadings rescaled
#' to variance-weighted length, `$rotation %*% diag(sdev * sqrt(n - 1))`: the same
#' map, on the scale a coordinate has rather than the unit length a direction has.
#' `$loadings` is then the margin that was not embedded, one row per sample, and
#' `$variance` does not change at all, so the axis labels a plot carries are the
#' same on both scales.
#'
#' Transposing by hand is right only when the rows of `data` really are features,
#' which is how an expression matrix usually arrives. Then the transpose puts the
#' data into this function's layout and `embedding_scale` chooses from there.
#'
#' # What is dropped before anything runs
#'
#' [stats::prcomp()] does not accept a missing value, so rows that are not complete
#' and finite across `feats` are dropped before it is called, and `design$n_dropped`
#' reports how many went. This is the listwise deletion the rest of the package
#' uses; nothing is imputed.
#'
#' A feature that takes a single value cannot be scaled — the division is by zero —
#' so with `scale = TRUE` it is left out with a message and named in
#' `design$dropped_feats`. With `scale = FALSE` it stays and becomes a component of
#' no variance.
#'
#' # Portability
#'
#' Everything but `$fit` is a data.frame, a character vector or a named list, so
#' dropping that one slot leaves an object that writes out as JSON. In a Python
#' transcription this is `sklearn.decomposition.PCA`.
#'
#' @param data A data.frame or a matrix in wide format, one row per sample and one
#'   column per feature. This is the same layout [compare_two_groups()] and
#'   [draw_heatmap()] take. Row names are kept as the sample labels, repeated ones
#'   included, since a sample name is a naming choice rather than a key; rows
#'   without a name are labelled by position.
#' @param feats Column names to reduce, or `NULL` for every numeric column of
#'   `data`. A non-numeric column is left out with a message, so a frame that
#'   carries a grouping column alongside the measurements can be passed as it is.
#' @param embedding_scale Which margin becomes the points of the picture:
#'   `"samples"`, the default, for one point per row of `data`, or `"features"` for
#'   one point per column. `design$point_type` reports it. See the details: this is
#'   not the same as transposing `data` yourself.
#' @param center,scale Whether to centre each feature and divide it by its standard
#'   deviation before rotating. Scaling is on by default because features are not
#'   measured on a common scale, and without it the feature with the widest units
#'   decides where the first component points. Both always apply to the **columns of
#'   `data`**, whatever `embedding_scale` is.
#'
#' @return An object of class `sa_reduction`, a plain list.
#'
#'   \describe{
#'     \item{`analysis`}{`"pca"`.}
#'     \item{`points`}{Labels of the things that became points — samples or
#'       features, as `embedding_scale` asked — in the row order `scores` follows.
#'       This is what `features` is to a comparison result.}
#'     \item{`design`}{What was reduced: `point_type`, either `"sample"` or
#'       `"feature"`, the `feats` kept and any `dropped_feats`, and the counts
#'       `n_samples`, `n_used`, `n_dropped` and `n_feats`. The counts describe the
#'       input rather than the picture, so `n_samples` is always rows of `data` and
#'       `n_feats` always kept columns, whichever margin became the points.}
#'     \item{`parameters`}{`embedding_scale` as resolved, and the two flags.}
#'     \item{`variance`}{One row per component: `component`, `sdev`, `prop_var` as
#'       a percentage and `cum_var`. Every component is here, since a share of the
#'       variance is only a share if the rest of it is there to be a share of.}
#'     \item{`loadings`}{The margin that was not embedded: one row per variable,
#'       `variables` beside one column per component.}
#'     \item{`scores`}{Coordinates: `points` beside one column per component.}
#'     \item{`engine`}{What computed the rotation.}
#'     \item{`fit`}{The [stats::prcomp()] object, always fitted to the samples.
#'       This is the slot that is not portable; dropping it leaves an object that
#'       writes out as JSON.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [perform_tsne()] and [perform_umap()], which answer the same question
#'   in a way that can follow structure a rotation cannot, and [draw_heatmap()],
#'   which shows the same wide input a cell at a time instead of a point at a time.
#'
#' @examples
#' ## The sample coordinates are `$scores`, and `$variance` is what labels the axes.
#' res <- perform_pca(iris[1:4])
#' res
#' head(res$scores)
#' res$variance
#'
#' plot(res$scores[c("PC1", "PC2")],
#'      xlab = paste0("PC1 (", round(res$variance$prop_var[1], 2), "%)"),
#'      ylab = paste0("PC2 (", round(res$variance$prop_var[2], 2), "%)"),
#'      col = as.integer(iris$Species), pch = 16)
#'
#' ## `$loadings` is the other axis: which features built each component.
#' res$loadings
#'
#' ## On data with a planted two-group structure. The group was never shown to the
#' ## rotation, so the split is its to find.
#' sim <- simulate_two_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 3)
#' by_samp <- perform_pca(sim$args$data)
#' table(group = sim$args$group, side = by_samp$scores$PC1 > 0)
#'
#' ## `embedding_scale = "features"` turns the features into the points. Here the
#' ## correlation blocks are planted, and the rotation was not told about them.
#' cor_mat <- make_block_cor(
#'   n_features = 8,
#'   blocks = list(list(features = 1:2, cor = 0.8),
#'                 list(features = 3:5, cor = 0.5),
#'                 list(features = 7:8, cor = 0.9))
#' )
#' blocks <- simulate_classification(cor_mat = cor_mat, seed = 2026)$args$data
#' by_feat <- perform_pca(blocks, feats = paste0("x_", 1:8),
#'                        embedding_scale = "features")
#' by_feat
#' by_feat$points
#'
#' plot(by_feat$scores[c("PC1", "PC2")], type = "n",
#'      xlab = paste0("PC1 (", round(by_feat$variance$prop_var[1], 2), "%)"),
#'      ylab = paste0("PC2 (", round(by_feat$variance$prop_var[2], 2), "%)"))
#' text(by_feat$scores[c("PC1", "PC2")], labels = by_feat$points, font = 2)
#'
#' ## The same fit read from the other end: the loadings of the sample scale,
#' ## rescaled from unit length to variance-weighted length, on the same axis
#' ## proportions.
#' same <- perform_pca(blocks, feats = paste0("x_", 1:8))
#' all.equal(by_feat$variance, same$variance)
#' round(diag(cor(by_feat$scores[-1], same$loadings[-1])), 6)
#'
#' @export
perform_pca <- function(data,
                        feats = NULL,
                        embedding_scale = c("samples", "features"),
                        center = TRUE,
                        scale = TRUE) {

  embedding_scale <- match.arg(embedding_scale)
  sa_check_flag(center, "center")
  sa_check_flag(scale, "scale")

  input <- sa_reduce_input(data, feats, scale, "perform_pca")
  x <- input$x
  pt <- sa_reduce_points(x, input$samples, embedding_scale)

  # The unscaled columns go in with the two flags rather than a matrix `scale()`
  # has already been over: `prcomp()` builds the same matrix out of the same call,
  # so the numbers are identical, and the object then carries the centre and the
  # scale, which is what lets `predict()` on it project a sample it has not seen.
  fit <- stats::prcomp(x, center = center, scale. = scale)

  if (identical(embedding_scale, "features")) {
    # One decomposition, read from the other end. `diag()` is given its size
    # explicitly because a single component would otherwise be read as a
    # dimension rather than as a value.
    weights <- fit$sdev * sqrt(nrow(x) - 1)
    coords <- fit$rotation %*% diag(weights, nrow = length(weights))
    colnames(coords) <- colnames(fit$rotation)
    other <- fit$x
  } else {
    coords <- fit$x
    other <- fit$rotation
  }

  v <- fit$sdev^2
  variance <- data.frame(
    component = colnames(fit$rotation),
    sdev      = fit$sdev,
    prop_var  = v / sum(v) * 100,
    cum_var   = cumsum(v) / sum(v) * 100,
    stringsAsFactors = FALSE
  )
  loadings <- data.frame(variables = rownames(other), stringsAsFactors = FALSE)
  loadings <- cbind(loadings, as.data.frame(other))
  rownames(loadings) <- NULL

  sa_new_reduction(
    analysis = "pca",
    points   = pt$points,
    # `design` describes the input, so its counts do not turn with
    # `embedding_scale`: `n_samples` is always rows of `data` and `feats` always
    # the columns kept. `point_type` is what says which of the two became the
    # points.
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
      scale           = scale
    ),
    variance = variance,
    loadings = loadings,
    scores   = sa_embedding_frame(coords, pt$points, "PC"),
    engine   = list(package = "stats", method = "prcomp",
                    label = "Principal component analysis",
                    overridden = character(0)),
    fit      = fit
  )
}
