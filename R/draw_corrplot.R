# The correlation matrix drawn as a heatmap, which is what a corrplot is. There
# is no second drawing engine here: draw_heatmap() already owns the cells, the
# diverging ramp, the colour key and the dendrograms, and this function is the
# three decisions a correlation matrix needs that a feature-by-sample matrix does
# not.
#
# The first is that nothing may be standardised. A coefficient is already on a
# common scale, and z-scoring it would replace the number the reader came for.
# The second is that both axes hold the same features, so they need one order
# rather than two: a symmetric matrix clustered twice can come back with its rows
# and columns in different orders, and the diagonal then wanders off the diagonal.
# The clustering is therefore done once here and draw_heatmap() is handed a matrix
# that is already in that order. The third is that a cell may be blanked for
# having no evidence behind it, which has to happen after the clustering, so that
# what is drawn and what the tree was built from stay the same matrix.

#' Draw a correlation matrix, with the cells that failed the test left blank
#'
#' A corrplot: one cell per pair of features, coloured by the coefficient on a
#' fixed -1 to 1 scale, with both axes in one order so that the diagonal runs
#' corner to corner and blocks of features that move together sit next to each
#' other. Given the p-values as well, the pairs that did not clear `sig_level`
#' are drawn as blank cells, so that what is left coloured is what there is
#' evidence for.
#'
#' [draw_heatmap()] draws it. This function decides what it is handed: the matrix
#' unscaled, the colour range fixed at the range a correlation can take, one
#' clustering shared by the two axes, and the blanking applied afterwards.
#'
#' @param cor_matrix The result of [summarize_association_stats()], or a
#'   correlation matrix on its own. A matrix must be square and symmetric, with
#'   the features as its dimnames, and hold at least two of them.
#' @param method Which coefficient to draw when `cor_matrix` is a result, naming
#'   one of the slots it holds. `NULL` takes the first method it was computed
#'   with. It must be `NULL` when `cor_matrix` is a matrix, there being only one
#'   matrix to draw.
#' @param pvalue Matrix of p-values laid out like `cor_matrix`, used only when
#'   `cor_matrix` is a matrix. When `cor_matrix` is a result the p-values come
#'   from the same slot as the coefficients, so this must be left `NULL`.
#' @param use_adjusted Whether to read `adj_pvalue` rather than `pvalue` out of
#'   the result. Ignored when `pvalue` is supplied directly.
#' @param sig_level Largest p-value a cell may have and still be drawn. Cells
#'   above it are blanked. The diagonal is never blanked, a feature not being
#'   tested against itself.
#' @param cluster Whether to reorder the features by clustering them. `FALSE`
#'   keeps the order they arrive in.
#' @param hclust_method Linkage handed to [stats::hclust()].
#' @param zlim Numeric length-2 range the colours span. The default is the range
#'   a correlation can take, so that the same colour means the same strength from
#'   one plot to the next.
#' @param anno Whether to write each coefficient on its cell, rounded to two
#'   decimal places by [draw_heatmap()]. Blanked cells stay blank.
#' @param main Plot title.
#' @param cex.anno Character expansion for those cell labels, relative to the
#'   smaller of the feature axis label sizes.
#' @param cex.axis,cex.main,cex.legend Character expansion for the axis labels,
#'   title and colour key, passed to [draw_heatmap()].
#' @param ... Passed to [draw_heatmap()]. The arguments this function decides
#'   (`data`, `group`, `group_lv`, `scale`, `cluster_feats`, `cluster_samples`,
#'   `anno`, `cex.anno`, `cex.axis`, `cex.main`, `cex.legend`) cannot be given
#'   again.
#'
#' @return A list, invisibly: everything [draw_heatmap()] returns, and beside it
#'
#'   \describe{
#'     \item{`corr`}{The matrix as it was drawn, in the drawn order and with the
#'       blanked cells `NA`.}
#'     \item{`pvalue`}{The p-values in that same order, or `NULL`.}
#'     \item{`order`}{The permutation of the input the clustering chose.}
#'     \item{`hclust`}{The [stats::hclust()] object behind it, or `NULL` when the
#'       features were not clustered.}
#'     \item{`n_masked`}{How many cells were blanked.}
#'   }
#'
#' @details
#' The distance the clustering runs on is `1 - cor()`, the same one
#' [cluster_hclust()] and [draw_heatmap()] mean by `dist_method =
#' "correlation"`, so a corrplot and a heatmap of the same features group them
#' the same way. It is computed once and both axes are permuted by it. A matrix
#' holding an `NA`, which is what a feature with no variance leaves behind, has
#' no distance for that feature, and rather than fail the features are left in
#' their input order with a message saying so.
#'
#' Blanking happens after the clustering rather than before it. A cell removed
#' for its p-value would otherwise change the tree, and the picture would no
#' longer be the matrix the reader is being shown.
#'
#' A cell whose p-value is `NA`, a pair that could not be tested, is left as it
#' arrived rather than blanked: there is no evidence against it either, and the
#' coefficient beside it is usually already `NA`.
#'
#' The matrix is symmetric, so the transpose [draw_heatmap()] takes on the way in
#' leaves it unchanged and the features come out on both axes.
#'
#' @seealso [summarize_association_stats()], which produces what this draws, and
#'   [draw_heatmap()], which draws it.
#'
#' @examples
#' feats <- c("mpg", "cyl", "disp", "hp", "drat", "wt")
#' res <- summarize_association_stats(mtcars, feats, methods = "pearson")
#'
#' ## Every pair, with the features clustered so that the blocks sit together
#' drawn <- draw_corrplot(res, main = "mtcars")
#' drawn$order
#'
#' ## Only the pairs that cleared the adjustment at 1%
#' draw_corrplot(res, sig_level = 0.01, cex.anno = 0.8, cex.axis = 0.7)
#'
#' ## A bare matrix, in the order it arrives and with nothing to blank
#' draw_corrplot(stats::cor(mtcars[feats]), cluster = FALSE)
#'
#' @export
draw_corrplot <- function(cor_matrix,
                          method = NULL,
                          pvalue = NULL,
                          use_adjusted = TRUE,
                          sig_level = 0.05,
                          cluster = TRUE,
                          hclust_method = c("average", "complete", "ward.D2"),
                          zlim = c(-1, 1),
                          anno = TRUE,
                          main = NULL,
                          cex.anno = 1,
                          cex.axis = 0.9,
                          cex.main = 1.5,
                          cex.legend = 1.2,
                          ...) {

  hclust_method <- match.arg(hclust_method)
  sa_check_flag(use_adjusted, "use_adjusted")
  sa_check_flag(cluster, "cluster")
  sa_check_flag(anno, "anno")
  sa_check_scalar_num(sig_level, "sig_level", 0, 1, lower_open = TRUE)
  sa_check_scalar_num(cex.anno, "cex.anno", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.axis, "cex.axis", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.main, "cex.main", 0, lower_open = TRUE)
  sa_check_scalar_num(cex.legend, "cex.legend", 0, lower_open = TRUE)

  read <- sa_corrplot_input(cor_matrix, method, pvalue, use_adjusted)
  corr <- read$corr
  pvalue <- read$pvalue
  feats <- colnames(corr)
  p <- ncol(corr)

  ord <- seq_len(p)
  hc <- NULL
  if (cluster) {
    d <- stats::as.dist(1 - corr)
    if (anyNA(d)) {
      message("Some pair of features has no correlation to measure a distance ",
              "between, so the features are drawn in the order they arrived.")
    } else {
      hc <- stats::hclust(d, method = hclust_method)
      ord <- hc$order
    }
  }

  drawn <- corr[ord, ord, drop = FALSE]
  pv <- if (is.null(pvalue)) NULL else pvalue[ord, ord, drop = FALSE]

  n_masked <- 0L
  if (!is.null(pv)) {
    blank <- !is.na(pv) & pv > sig_level
    diag(blank) <- FALSE
    n_masked <- sum(blank)
    drawn[blank] <- NA_real_
  }

  out <- draw_heatmap(data = drawn,
                      group = NULL,
                      group_lv = NULL,
                      scale = "none",
                      zlim = zlim,
                      cluster_feats = FALSE,
                      cluster_samples = FALSE,
                      anno = anno,
                      cex.anno = cex.anno,
                      main = main,
                      cex.axis = cex.axis,
                      cex.main = cex.main,
                      cex.legend = cex.legend,
                      ...)

  out$corr <- drawn
  out$pvalue <- pv
  out$order <- ord
  out$hclust <- hc
  out$n_masked <- n_masked
  # feats is what the order indexes into, which is not recoverable from the
  # drawn matrix once the permutation has been applied to it.
  out$feats <- feats

  invisible(out)
}


#' Read the coefficient matrix and its p-values out of either kind of input
#'
#' `draw_corrplot()` takes the result of [summarize_association_stats()] so that
#' the two functions meet in one line, and a bare matrix so that a correlation
#' computed some other way can still be drawn. Which one arrived decides where
#' the p-values come from, and the two ways of naming them cannot both be used.
#'
#' @param cor_matrix The caller's `cor_matrix`.
#' @param method The caller's `method`.
#' @param pvalue The caller's `pvalue`.
#' @param use_adjusted The caller's `use_adjusted`.
#'
#' @return A list with the validated `corr` and `pvalue`, the latter possibly
#'   `NULL`.
#'
#' @keywords internal
#' @noRd
sa_corrplot_input <- function(cor_matrix, method, pvalue, use_adjusted) {
  is_result <- is.list(cor_matrix) && !is.data.frame(cor_matrix) &&
    !is.null(cor_matrix$design$methods)

  if (is_result) {
    if (!is.null(pvalue)) {
      stop("`pvalue` cannot be given for a `summarize_association_stats()` ",
           "result: the p-values come from the same slot as the coefficients. ",
           "Use `use_adjusted` to choose between them.", call. = FALSE)
    }
    methods <- cor_matrix$design$methods
    if (is.null(method)) {
      method <- methods[1]
    } else if (!is.character(method) || length(method) != 1L ||
               !method %in% methods) {
      stop("`method` must name one of the methods `cor_matrix` holds: ",
           paste(methods, collapse = ", "), ".", call. = FALSE)
    }
    slot <- cor_matrix[[method]]
    corr <- slot$corr
    pvalue <- if (use_adjusted) slot$adj_pvalue else slot$pvalue
  } else {
    if (!is.null(method)) {
      stop("`method` names a slot of a `summarize_association_stats()` result, ",
           "and `cor_matrix` is a matrix. Leave it NULL.", call. = FALSE)
    }
    corr <- cor_matrix
  }

  if (!is.matrix(corr) || !is.numeric(corr)) {
    stop("`cor_matrix` must be a numeric correlation matrix or the result of ",
         "`summarize_association_stats()`.", call. = FALSE)
  }
  if (nrow(corr) != ncol(corr)) {
    stop("`cor_matrix` must be square, but is ", nrow(corr), " by ",
         ncol(corr), ".", call. = FALSE)
  }
  if (ncol(corr) < 2L) {
    stop("`draw_corrplot()` needs at least 2 features to draw, but got ",
         ncol(corr), ".", call. = FALSE)
  }
  if (!isTRUE(all.equal(corr, t(corr), check.attributes = FALSE))) {
    stop("`cor_matrix` must be symmetric: a correlation between two features ",
         "is one number, so the two cells that hold it must agree.",
         call. = FALSE)
  }
  finite <- corr[is.finite(corr)]
  if (length(finite) > 0L && (min(finite) < -1 || max(finite) > 1)) {
    stop("`cor_matrix` holds value(s) outside [-1, 1], so it is not a matrix ",
         "of correlations.", call. = FALSE)
  }

  # A matrix from cor() always carries its dimnames; one assembled by hand may
  # not, and draw_heatmap() would then invent a name for one axis only.
  if (is.null(colnames(corr))) {
    colnames(corr) <- paste0("V", seq_len(ncol(corr)))
  }
  rownames(corr) <- colnames(corr)

  if (!is.null(pvalue)) {
    if (!is.matrix(pvalue) || !is.numeric(pvalue) ||
        !identical(dim(pvalue), dim(corr))) {
      stop("`pvalue` must be a numeric matrix laid out like `cor_matrix`: ",
           nrow(corr), " by ", ncol(corr), ".", call. = FALSE)
    }
    if (!is.null(colnames(pvalue)) &&
        !identical(colnames(pvalue), colnames(corr))) {
      stop("`pvalue` must name the same features as `cor_matrix`, in the same ",
           "order.", call. = FALSE)
    }
    dimnames(pvalue) <- dimnames(corr)
  }

  list(corr = corr, pvalue = pvalue)
}
