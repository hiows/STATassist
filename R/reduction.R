# The result contract shared by the three reductions, and the counterpart of
# `result.R` for the unsupervised family. A comparison result is organised around a
# feature axis and a model around a term axis. A reduction repeats neither: what
# every table here repeats is `points`. A point is whichever margin was embedded —
# a sample by default, a feature when asked — and `design$point_type` says which,
# because a reduction is the one result in this package whose row axis is the
# caller's choice rather than the question's.
#
# There are three functions rather than one for the reason the comparison wrappers
# run several tests in one call and this family does not: a comparison's tests
# answer the same question and can be read down one table, while these three
# answer in coordinates that share no scale. `perform_pca()` is a rotation, so it
# is reproducible, invertible and readable as "which features moved this sample",
# but it can only ever draw straight structure. `perform_tsne()` and
# `perform_umap()` find structure that curves and cannot say which feature made it,
# and each pays for that with an arbitrary global scale. Nothing can be joined
# across the three except the points, and the points are what this contract fixes:
# every one of them lists its coordinates in the same row order, so two of the
# three can be plotted against each other and against anything else read from the
# same rows.
#
# One rule is broken here, in the same place `sa_model` breaks it. `$fit` holds the
# engine object, because a `prcomp` object is what anyone who has done a PCA in R
# already knows how to read and hiding it behind our own tables would make the
# familiar call the long way round. Only the `umap` object is genuinely unportable
# — its `config$metric.function` is a function — but the slot is declared as the
# exception for all three so that the rule stays one sentence: drop `$fit` and the
# object writes out as JSON.

#' The reductions this contract covers
#'
#' `analysis` names the method rather than the family, as `"linear_regression"` and
#' `"svm"` do in `sa_model`, so a result says which of the three it came from
#' without anything having to read `engine`.
#'
#' @keywords internal
#' @noRd
sa_reduction_analyses <- function() {
  c("pca", "tsne", "umap")
}


#' Assemble a dimensionality reduction result object
#'
#' The checks here guard the contract rather than the user's input, so they fire
#' only on a mistake inside the package and say so. What they are guarding is the
#' one promise the object makes: `scores` is aligned with `points` by position, so
#' a reduction can be plotted against anything else read from the same rows.
#'
#' @param analysis Which reduction this is: `"pca"`, `"tsne"` or `"umap"`.
#' @param points Point labels, the row order every table follows.
#' @param design Named list describing what was reduced.
#' @param parameters Named list of the choices as they were used.
#' @param scores data.frame of coordinates, one row per point.
#' @param variance,loadings The two PCA tables, or `NULL` for an embedding.
#' @param engine Named list naming what computed the reduction.
#' @param fit The engine object.
#'
#' @keywords internal
#' @noRd
sa_new_reduction <- function(analysis,
                             points,
                             design,
                             parameters,
                             scores,
                             variance = NULL,
                             loadings = NULL,
                             engine,
                             fit) {

  if (!analysis %in% sa_reduction_analyses()) {
    stop("internal error: `analysis` must be one of ",
         paste(sa_reduction_analyses(), collapse = ", "), ".", call. = FALSE)
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
  if (!is.data.frame(scores) || !identical(scores$points, points)) {
    stop("internal error: `scores` is not a data.frame aligned with `points`.",
         call. = FALSE)
  }
  for (nm in c("package", "method", "label", "overridden")) {
    if (is.null(engine[[nm]])) {
      stop("internal error: `engine` is missing `", nm, "`.", call. = FALSE)
    }
  }
  # The two PCA tables are the rotation's other side, so either both are here or
  # neither is, and which it is follows from the method rather than from what the
  # engine happened to return.
  is_pca <- identical(analysis, "pca")
  if (is_pca != !is.null(variance) || is_pca != !is.null(loadings)) {
    stop("internal error: `variance` and `loadings` are present exactly when ",
         "the analysis is a principal component analysis.", call. = FALSE)
  }
  # `loadings` is the margin that was not embedded, so what it should have one row
  # of depends on which margin that is.
  n_vars <- if (identical(design$point_type, "feature")) {
    design$n_used
  } else {
    design$n_feats
  }
  if (!is.null(loadings) && nrow(loadings) != n_vars) {
    stop("internal error: `loadings` has ", nrow(loadings), " row(s) for ",
         n_vars, " variable(s).", call. = FALSE)
  }

  structure(
    c(
      list(analysis   = analysis,
           points     = points,
           design     = design,
           parameters = parameters),
      # A slot is present only when it has something to say, so an embedding does
      # not offer an empty version of the rotation's tables.
      if (!is.null(variance)) list(variance = variance),
      if (!is.null(loadings)) list(loadings = loadings),
      list(scores   = scores,
           engine   = engine,
           fit      = fit,
           metadata = sa_metadata())
    ),
    class = c("sa_reduction", "sa_result")
  )
}


#' Print a dimensionality reduction
#'
#' Summarises what was reduced and how far, rather than printing the coordinates.
#' Those are in `x$scores`, and the engine object is `x$fit`.
#'
#' @param x A reduction, as returned by [perform_pca()], [perform_tsne()] or
#'   [perform_umap()].
#' @param n Maximum number of components to report the variance of. The rest are
#'   counted. Ignored by an embedding, which has no components.
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' perform_pca(iris[1:4])
#'
#' @export
print.sa_reduction <- function(x, n = 3L, ...) {
  n <- sa_check_count(n, "n", 0)
  design <- x$design
  params <- x$parameters

  cat("<sa_reduction> ", x$analysis, "\n", sep = "")
  # `data` says what came in and `points` says which side of it was embedded, so
  # the two lines together read as the whole choice rather than as a repetition.
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

  if (!is.null(x$variance)) {
    shown <- utils::head(x$variance, n)
    if (nrow(shown) > 0L) {
      sa_cat_field("variance", paste0(
        paste0(shown$component, " ", sa_fmt_num(shown$prop_var, 4), "%",
               collapse = ", "),
        "  (", nrow(shown), " of ", nrow(x$variance), " component(s), ",
        sa_fmt_num(shown$cum_var[nrow(shown)], 4), "% cumulative)"
      ))
    }
  }
  if (identical(x$analysis, "tsne")) {
    cat("  tsne     : ", params$n_dim, " dimension(s), perplexity = ",
        sa_fmt_num(params$perplexity, 3), ", theta = ",
        sa_fmt_num(params$theta, 3),
        if (!is.null(params$seed)) paste0("  (seed = ", params$seed, ")"),
        "\n", sep = "")
  }
  if (identical(x$analysis, "umap")) {
    cat("  umap     : ", params$n_dim, " dimension(s), method = ",
        params$method, ", n_neighbors = ",
        params$n_neighbors, ", min_dist = ", sa_fmt_num(params$min_dist, 3),
        ", ", params$metric,
        if (!is.null(params$seed)) paste0("  (seed = ", params$seed, ")"),
        "\n", sep = "")
  }
  if (length(design$dropped_feats) > 0L) {
    sa_cat_field("dropped", paste0(
      paste(design$dropped_feats, collapse = ", "), " (no variance)"
    ))
  }

  invisible(x)
}
