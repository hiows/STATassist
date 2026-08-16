#' Correlation between every pair of features, with all three coefficients
#'
#' Reduces a set of features to the association between each pair of them, as a
#' square matrix per quantity: the coefficient, its p-value, the p-value adjusted
#' across the pairs, and how many observations the pair shared. Pearson,
#' Spearman and Kendall are reported side by side on the same pairs, the way the
#' comparison functions report a parametric, a rank-based and a robust test side
#' by side, so that a linear coefficient and a monotonic one disagreeing is
#' visible rather than a matter of which call was made.
#'
#' This is a screen rather than a test of one hypothesis. It is the companion of
#' [summarize_descriptive_stats()], which reduces one feature at a time to a row
#' of its own; this reduces a pair at a time. Neither returns a `sa_comparison`.
#'
#' @param data A data.frame (or matrix) in wide format, one row per observation
#'   and one column per feature.
#' @param feats Character vector of numeric column names in `data` to correlate,
#'   in output order, or `NULL` for every numeric column of `data`. At least two
#'   are needed.
#' @param methods Which coefficients to compute, drawn from `"pearson"`,
#'   `"spearman"` and `"kendall"`. Each one named gets a slot of the result, in
#'   the order given; the ones not named have no slot at all.
#' @param adj_type Multiplicity adjustment applied across the pairs, any method
#'   from [stats::p.adjust.methods].
#' @param use How a missing value is handled. `"pairwise.complete.obs"` reads
#'   each pair on the observations that pair shares, and `"complete.obs"` drops
#'   any row with a missing value in any feature first, so that every pair is
#'   read on one set of rows.
#'
#' @return A plain list. One slot per entry of `methods`, named after it and
#'   holding four features-by-features matrices with the same dimnames:
#'
#'   \describe{
#'     \item{`corr`}{The coefficient, from [stats::cor()].}
#'     \item{`pvalue`}{The two-sided p-value from [stats::cor.test()].}
#'     \item{`adj_pvalue`}{`pvalue` adjusted by `adj_type` across the pairs.}
#'     \item{`n`}{Observations the pair shared.}
#'   }
#'
#'   Beside those, `design` records `feats`, `n_obs`, `methods`, `adj_type` and
#'   `use`.
#'
#' @details
#' Every matrix is symmetric. The upper triangle is computed and mirrored, so a
#' pair is tested once rather than twice.
#'
#' The diagonal is a property of the matrix rather than an estimate: `corr` is 1
#' and `pvalue` and `adj_pvalue` are `NA`, since a feature is not tested against
#' itself. `n` on the diagonal is the number of observations that feature has.
#'
#' The adjustment runs over the pairs that produced a p-value. A pair
#' [stats::cor.test()] refused, having fewer than three shared observations or no
#' variance on one side, is not a test that was performed, and counting it into
#' the family would shrink the other pairs for a comparison that never happened.
#' Such a pair comes back `NA` in all three of `corr`, `pvalue` and `adj_pvalue`
#' rather than aborting the screen, and features with no variance at all are
#' named in a message.
#'
#' Missing and non-finite values are treated alike, an `Inf` being as much "no
#' value to correlate" as an `NA`, and `n` counts what was left.
#'
#' The cost is one [stats::cor.test()] per pair per method: thirty features are
#' 435 pairs, and 1305 tests with all three methods asked for. Kendall's is the
#' slowest of the three, so naming `methods` is worth doing on a wide frame.
#'
#' @seealso [draw_corrplot()] to draw a method's `corr` with the non-significant
#'   cells left blank, and [summarize_descriptive_stats()] for the one-feature-at-a-time
#'   summary this is the pairwise counterpart of.
#'
#' @examples
#' feats <- c("mpg", "disp", "hp", "wt")
#'
#' ## All three coefficients on the same pairs
#' res <- summarize_association_stats(mtcars, feats)
#' round(res$pearson$corr, 3)
#' round(res$pearson$adj_pvalue, 4)
#'
#' ## A linear coefficient and a monotonic one on the same pair
#' c(pearson  = res$pearson$corr["mpg", "disp"],
#'   spearman = res$spearman$corr["mpg", "disp"])
#'
#' ## One method only
#' rho <- summarize_association_stats(mtcars, feats, methods = "spearman")
#' round(rho$spearman$corr, 3)
#'
#' @export
summarize_association_stats <- function(data,
                                        feats = NULL,
                                        methods = c("pearson", "spearman",
                                                    "kendall"),
                                        adj_type = "BH",
                                        use = c("pairwise.complete.obs",
                                                "complete.obs")) {

  use <- match.arg(use)
  sa_check_p_adjust(adj_type, "adj_type")

  known <- c("pearson", "spearman", "kendall")
  if (!is.character(methods) || length(methods) == 0L || anyNA(methods)) {
    stop("`methods` must be a non-empty character vector drawn from: ",
         paste(known, collapse = ", "), ".", call. = FALSE)
  }
  unknown <- setdiff(methods, known)
  if (length(unknown) > 0L) {
    stop("`methods` must be drawn from: ", paste(known, collapse = ", "),
         ". Not recognised: ", paste(unknown, collapse = ", "), ".",
         call. = FALSE)
  }
  dup_methods <- unique(methods[duplicated(methods)])
  if (length(dup_methods) > 0L) {
    stop("`methods` contains duplicated names: ",
         paste(dup_methods, collapse = ", "), call. = FALSE)
  }

  # A matrix keeps its column names as the feature names; repeated row names are
  # a sample naming choice and nothing here reads them, so they are dropped
  # rather than allowed to fail the conversion.
  if (is.matrix(data)) {
    rownames(data) <- NULL
    data <- as.data.frame(data)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or a matrix.", call. = FALSE)
  }
  if (is.null(feats)) {
    feats <- names(data)[vapply(data, is.numeric, logical(1))]
    if (length(feats) == 0L) {
      stop("`data` holds no numeric column to correlate.", call. = FALSE)
    }
  }

  input <- sa_validate_wide_input(data, feats, group = NULL, group_lv = NULL)
  data <- input$data
  feats <- input$feats

  if (length(feats) < 2L) {
    stop("`summarize_association_stats()` needs at least 2 features to ",
         "correlate, but got ", length(feats), ".", call. = FALSE)
  }

  x <- as.matrix(data[feats])
  # An Inf is as much "no value to correlate" as an NA, and letting one through
  # would turn every coefficient that feature has into NaN.
  x[!is.finite(x)] <- NA
  colnames(x) <- feats

  if (use == "complete.obs") {
    keep <- stats::complete.cases(x)
    if (!any(keep)) {
      stop("`use = \"complete.obs\"` leaves no row: every row has a missing ",
           "value in at least one of `feats`.", call. = FALSE)
    }
    if (any(!keep)) {
      message("Dropped ", sum(!keep),
              " row(s) with a missing value, as `use = \"complete.obs\"` asks.")
    }
    x <- x[keep, , drop = FALSE]
  }

  spread <- apply(x, 2, stats::sd, na.rm = TRUE)
  flat <- !is.finite(spread) | spread == 0
  if (any(flat)) {
    message(sum(flat), " feature(s) have no variance to correlate and come ",
            "back as NA: ", paste(feats[flat], collapse = ", "), ".")
  }

  out <- lapply(methods, function(m) sa_association_matrices(x, m, adj_type))
  names(out) <- methods

  out$design <- list(feats    = feats,
                     n_obs    = nrow(x),
                     methods  = methods,
                     adj_type = adj_type,
                     use      = use)
  out
}
