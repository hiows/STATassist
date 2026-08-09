# Everything the three reduction functions share, which is the whole of the input
# side. `perform_pca()`, `perform_tsne()` and `perform_umap()` differ in what they
# do with a matrix and not at all in how they read one, so a caller who moves from
# one to the next finds the same rows dropped for the same reasons and the same
# `design` describing them. That is what makes the three comparable: a cluster
# that only one of them finds is a fact about the method, and it can only be read
# that way if the three were handed the same numbers.
#
# `embedding_scale` lives here for the same reason, and it is not a preprocessing
# option. It names which margin of the input becomes a point, and each function
# answers it in the way its engine allows. A principal component analysis is a
# singular value decomposition, so one fit answers on both margins at once and the
# matrix never has to turn: `$x` holds the samples and `$rotation` holds the
# features. t-SNE and UMAP have no such dual and only ever embed the rows they
# were handed, so on the feature scale it is the transpose that goes in.
#
# Which is why the argument exists rather than a note telling callers to transpose.
# `center` and `scale` always apply to the features, whichever margin is being
# embedded, because that is what scaling a data set means; standardising the
# transpose instead standardises samples and answers a third question that looks
# exactly like an answer to this one.

#' Read the matrix a reduction is computed on out of the caller's frame
#'
#' Everything that decides which rows and which columns reach an engine happens
#' here, so that the three reductions cannot be given three different matrices.
#'
#' @param data The caller's `data`.
#' @param feats The caller's `feats`, possibly `NULL`.
#' @param scale Whether the features will be divided by their standard deviation,
#'   which is what makes a feature of no variance impossible rather than useless.
#' @param fn Name of the calling function, so that the message about a table too
#'   small to reduce names the function the caller actually called.
#'
#' @return A list with `x`, the numeric matrix of usable rows and kept features;
#'   `samples`, its row labels; `n_samples`, the rows that arrived; `n_dropped`,
#'   those that could not be used; and `dropped_feats`.
#'
#' @keywords internal
#' @noRd
sa_reduce_input <- function(data, feats, scale, fn) {
  # A matrix carries its labels in its dimnames and the row names are needed after
  # the rows have been filtered, so they ride through the validator as `id`, the
  # one argument it filters beside the data. This is what `draw_heatmap()` does
  # with the same input.
  if (is.matrix(data) && is.null(colnames(data))) {
    colnames(data) <- paste0("V", seq_len(ncol(data)))
  }
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("`data` must be a data.frame or a matrix.", call. = FALSE)
  }
  sample_labels <- rownames(data)
  if (is.null(sample_labels)) {
    sample_labels <- as.character(seq_len(NROW(data)))
  }
  if (is.matrix(data)) {
    # as.data.frame() inside the validator rejects repeated row names, which are a
    # sample naming choice rather than an error. The names are already saved.
    rownames(data) <- NULL
  }

  if (is.null(feats)) {
    numeric_col <- vapply(as.data.frame(data), is.numeric, logical(1))
    if (!any(numeric_col)) {
      stop("`data` holds no numeric column, so there is nothing to reduce.",
           call. = FALSE)
    }
    if (!all(numeric_col)) {
      # Not an error: a wide frame usually carries the grouping column beside the
      # measurements, and asking the caller to strip it would be asking them to
      # write `feats` out by hand.
      message("Left out ", sum(!numeric_col), " non-numeric column(s): ",
              paste(names(numeric_col)[!numeric_col], collapse = ", "), ".")
    }
    feats <- names(numeric_col)[numeric_col]
  }

  # One synthetic level, as an ungrouped summarize_descriptive_stats() call does:
  # the validator is written in terms of group levels and a reduction has none.
  input <- sa_validate_wide_input(data, feats,
                                  group = rep("all", NROW(data)),
                                  group_lv = "all",
                                  id = sample_labels,
                                  min_levels = 1L)
  x <- as.matrix(input$data[input$feats])
  rownames(x) <- NULL

  # No engine here takes a hole, so the rows go before any of them is called. An
  # infinite value is dropped with the missing ones: it survives
  # `complete.cases()` and then turns the whole column into `NaN` on the way
  # through `scale()`.
  usable <- rowSums(!is.finite(x)) == 0L
  n_dropped <- sum(!usable)
  if (n_dropped > 0L) {
    message("Dropped ", n_dropped, " row(s) that are not complete and finite ",
            "across the feature(s) being reduced.")
  }
  x <- x[usable, , drop = FALSE]
  samples <- input$id[usable]
  # Put the labels back on the matrix the engine is given, so that the fitted
  # object carries them too and can be read on its own terms. They were only off
  # it to get past `as.data.frame()`, which refuses a repeated one.
  rownames(x) <- samples

  dropped_feats <- character(0)
  if (scale && nrow(x) > 1L) {
    # A feature that never moves would be divided by zero. It carries nothing
    # either way, so it is named and left out rather than allowed to turn the
    # whole matrix into `NaN`.
    flat <- apply(x, 2, function(v) !is.finite(stats::sd(v)) || stats::sd(v) == 0)
    if (any(flat)) {
      dropped_feats <- colnames(x)[flat]
      message("Left out ", length(dropped_feats),
              " feature(s) of no variance, which `scale = TRUE` cannot rescale: ",
              paste(dropped_feats, collapse = ", "), ".")
      x <- x[, !flat, drop = FALSE]
    }
  }

  if (nrow(x) < 2L || ncol(x) < 2L) {
    stop("`", fn, "()` needs at least 2 samples and 2 features to reduce, but ",
         "got ", nrow(x), " usable sample(s) and ", ncol(x),
         " usable feature(s).", call. = FALSE)
  }

  list(x = x,
       samples = samples,
       n_samples = NROW(data),
       n_dropped = n_dropped,
       dropped_feats = dropped_feats)
}


#' What one point is, and what it is called
#'
#' The row axis of a reduction is the one axis in this package the caller chooses,
#' so it is resolved in one place and every function reports it the same way.
#'
#' @param x The matrix `sa_reduce_input()` returned.
#' @param samples Its row labels.
#' @param embedding_scale The caller's `embedding_scale`, already matched.
#'
#' @return A list with `points` and `point_type`.
#'
#' @keywords internal
#' @noRd
sa_reduce_points <- function(x, samples, embedding_scale) {
  if (identical(embedding_scale, "features")) {
    list(points = colnames(x), point_type = "feature")
  } else {
    list(points = samples, point_type = "sample")
  }
}


#' The matrix an embedding engine is handed
#'
#' For [perform_tsne()] and [perform_umap()], which embed the rows they are given
#' and nothing else. The features are standardised first and the transpose is then
#' embedded as it stands, which is the definition under which the feature scale
#' agrees with what a principal component analysis says about the same features.
#' Standardising after the transpose would standardise samples.
#'
#' [perform_pca()] does not call this. `prcomp()` builds the same matrix out of the
#' same call to [scale()] and then remembers the centre and the scale, which is
#' what lets `predict()` on it project a row it has not seen.
#'
#' @param x The matrix `sa_reduce_input()` returned, one row per sample.
#' @param embedding_scale Which margin becomes a point.
#' @param center,scale The caller's two flags, applied to the features.
#'
#' @return A numeric matrix, one row per point.
#'
#' @keywords internal
#' @noRd
sa_reduce_embedding_matrix <- function(x, embedding_scale, center, scale) {
  xs <- base::scale(x, center = center, scale = scale)
  if (identical(embedding_scale, "features")) t(xs) else xs
}


#' Say out loud how small the neighbourhood came out
#'
#' Not a warning: the run goes through and its output is a picture like any other.
#' It is said because the derived neighbourhood is the whole behaviour of these two
#' methods, and someone embedding eight features is the last person who would think
#' to check what it came out as.
#'
#' @param n_points Number of points being embedded.
#' @param point_type What one point is.
#' @param size The neighbourhood as it was derived, already formatted.
#'
#' @keywords internal
#' @noRd
sa_reduce_few_points <- function(n_points, point_type, size) {
  if (n_points >= 16L) {
    return(invisible(NULL))
  }
  message("Only ", n_points, " ", point_type, "(s) to embed (", size,
          "). This method describes a neighbourhood, and below about 16 points ",
          "there is not much of one to describe. `perform_pca()` is not ",
          "governed by one.")
  invisible(NULL)
}


#' The neighbourhood size t-SNE is run at
#'
#' `Rtsne` refuses `n - 1 < 3 * perplexity`, since a perplexity is a number of
#' neighbours and there have to be that many rows to be neighbours with. A value
#' that was asked for and cannot be honoured is an error naming the limit, the way
#' an `mtry` above the predictor count is in [fit_rf()]. A value this function
#' derived and cannot honour is an error too, and says instead that there are too
#' few points for the method at all: when the whole call is one method there is
#' nothing left for a skip to protect, which is what this used to be for.
#'
#' @param perplexity The caller's `perplexity`, possibly `NULL`.
#' @param n Number of points being embedded.
#' @param point_type What one point is, `"sample"` or `"feature"`, so that the
#'   limit is quoted in terms of the margin the caller chose.
#'
#' @return The perplexity to run at.
#'
#' @keywords internal
#' @noRd
sa_tsne_perplexity <- function(perplexity, n, point_type = "sample") {
  upper <- (n - 1) / 3
  if (is.null(perplexity)) {
    derived <- min(30, floor(upper))
    if (derived < 1) {
      stop("`perform_tsne()` cannot embed ", n, " ", point_type,
           "(s): they admit no perplexity of 1 or more, since Rtsne requires ",
           "3 * perplexity <= n - 1. `perform_pca()` has no neighbourhood and ",
           "is not limited this way.", call. = FALSE)
    }
    return(derived)
  }
  sa_check_scalar_num(perplexity, "perplexity", 1)
  if (perplexity > upper) {
    stop("`perplexity` must not exceed (n - 1) / 3, which is ",
         format(upper, digits = 4), " for the ", n, " usable ", point_type,
         "(s), but is ", perplexity,
         ". t-SNE keeps a neighbourhood rather than a pair, and there are not ",
         "that many points to fill one.", call. = FALSE)
  }
  perplexity
}


#' The neighbourhood size UMAP is run at
#'
#' Two is the smallest neighbourhood there is and the number of points is the
#' largest, and `umap` rejects both ends itself. The same split as
#' `sa_tsne_perplexity()`: an impossible request and an impossible derived value
#' are both errors, and they say different things about whose mistake it was.
#'
#' @param n_neighbors The caller's `n_neighbors`, possibly `NULL`.
#' @param n Number of points being embedded.
#' @param point_type What one point is, `"sample"` or `"feature"`.
#'
#' @return The neighbourhood size to run at.
#'
#' @keywords internal
#' @noRd
sa_umap_neighbors <- function(n_neighbors, n, point_type = "sample") {
  if (is.null(n_neighbors)) {
    derived <- min(15L, n)
    if (derived < 2L) {
      stop("`perform_umap()` cannot embed ", n, " ", point_type,
           "(s): they admit no neighbourhood of 2 or more. `perform_pca()` has ",
           "no neighbourhood and is not limited this way.", call. = FALSE)
    }
    return(derived)
  }
  n_neighbors <- sa_check_count(n_neighbors, "n_neighbors", 2)
  if (n_neighbors > n) {
    stop("`n_neighbors` must not exceed the ", n, " usable ", point_type,
         "(s) being embedded, but is ", n_neighbors, ". A ", point_type,
         " cannot have more neighbours than there are ", point_type, "s.",
         call. = FALSE)
  }
  n_neighbors
}


#' Put an embedding beside the points it describes
#'
#' The engines do not agree on what a coordinate matrix looks like: `prcomp()`
#' names its columns and keeps row names, `Rtsne` does neither. The row a
#' coordinate belongs to is its position, as it is everywhere else in this
#' package, so the labels come from `points` rather than from whatever the engine
#' happened to carry through.
#'
#' @param m Coordinate matrix, one row per point.
#' @param points Point labels.
#' @param prefix What to call the columns when the engine did not name them.
#'
#' @return A data.frame of `points` beside one column per dimension.
#'
#' @keywords internal
#' @noRd
sa_embedding_frame <- function(m, points, prefix) {
  colnames(m) <- if (is.null(colnames(m))) {
    paste0(prefix, seq_len(ncol(m)))
  } else {
    colnames(m)
  }
  rownames(m) <- NULL
  out <- cbind(data.frame(points = points, stringsAsFactors = FALSE),
               as.data.frame(m))
  rownames(out) <- NULL
  out
}
