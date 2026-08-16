# The kernel behind summarize_association_stats(). A correlation screen asks one
# question of every pair of features, so the four matrices come out of a single
# pass over the upper triangle and are mirrored into the lower one, rather than
# testing each pair twice.
#
# The diagonal is set rather than estimated. A feature's correlation with itself
# is a property of the matrix and not a test, so `corr` carries 1 there even for a
# feature with no variance to correlate, while `pvalue` and `adj_pvalue` carry NA.
# That convention is what lets draw_corrplot() mask a cell by its p-value without
# having to make an exception of the diagonal.

#' The p-value of one correlation test, or `NA` when there is no test to run
#'
#' [stats::cor.test()] refuses a pair it cannot test, fewer than three complete
#' observations or a vector with no variance, and the refusal is an error rather
#' than a p-value. A screen over every pair cannot stop at the first such pair, so
#' the refusal becomes `NA` and the rest of the matrix is still computed.
#'
#' The tie warning the exact branches of Spearman's and Kendall's tests emit is
#' expected on real data and is not passed on. What comes back when it fires is
#' the normal approximation, which is the p-value the engine itself falls back to.
#'
#' @param u,v The two columns, with non-finite values already `NA`.
#' @param method `"pearson"`, `"spearman"` or `"kendall"`.
#'
#' @return The p-value as a length-1 numeric, or `NA_real_`.
#'
#' @keywords internal
#' @noRd
sa_cor_test_pvalue <- function(u, v, method) {
  out <- tryCatch(
    suppressWarnings(stats::cor.test(u, v, method = method)),
    error = function(e) NULL
  )
  if (is.null(out) || !is.finite(out$p.value)) NA_real_ else out$p.value
}


#' How many observations each pair of features shares
#'
#' [base::crossprod()] on the indicator of what is present counts, for every pair
#' of columns, the rows where both are. The diagonal is then the count for one
#' column on its own, which is what the pair of a feature with itself would have.
#'
#' @param x Numeric matrix, features in columns, non-finite values already `NA`.
#'
#' @return An integer matrix, features by features.
#'
#' @keywords internal
#' @noRd
sa_pairwise_n <- function(x) {
  present <- !is.na(x)
  storage.mode(present) <- "double"
  n <- crossprod(present)
  storage.mode(n) <- "integer"
  dimnames(n) <- list(colnames(x), colnames(x))
  n
}


#' The four matrices one correlation method produces
#'
#' @param x Numeric matrix, features in columns, non-finite values already `NA`
#'   and any listwise deletion already done.
#' @param method `"pearson"`, `"spearman"` or `"kendall"`.
#' @param adj_type Multiplicity adjustment, already checked.
#'
#' @return A list of the four features-by-features matrices `corr`, `pvalue`,
#'   `adj_pvalue` and `n`.
#'
#' @details
#' The family the adjustment runs over is the pairs that produced a p-value, not
#' every cell of the triangle. A pair [stats::cor.test()] refused is not a test
#' that was performed, and counting it would shrink the others for a comparison
#' that never happened.
#'
#' @keywords internal
#' @noRd
sa_association_matrices <- function(x, method, adj_type) {
  feats <- colnames(x)
  p <- ncol(x)
  empty <- matrix(NA_real_, p, p, dimnames = list(feats, feats))

  # A pair with no variance on either side has no correlation, which cor() warns
  # about once per pair. The NA it returns is the answer here.
  corr <- suppressWarnings(
    stats::cor(x, method = method, use = "pairwise.complete.obs")
  )
  diag(corr) <- 1

  pvalue <- empty
  for (j in seq_len(p - 1L)) {
    for (k in seq.int(j + 1L, p)) {
      pv <- sa_cor_test_pvalue(x[, j], x[, k], method)
      pvalue[j, k] <- pv
      pvalue[k, j] <- pv
    }
  }

  adj_pvalue <- empty
  tested <- upper.tri(pvalue) & !is.na(pvalue)
  adj_pvalue[tested] <- stats::p.adjust(pvalue[tested], method = adj_type)
  lower <- lower.tri(adj_pvalue)
  adj_pvalue[lower] <- t(adj_pvalue)[lower]

  list(corr       = corr,
       pvalue     = pvalue,
       adj_pvalue = adj_pvalue,
       n          = sa_pairwise_n(x))
}
