# The correlation matrix the supervised learning simulators take. It is a
# separate function rather than an argument of them because a matrix written out
# by hand is unreadable past about four predictors, while the structure that
# matters is almost always blocks: a handful of predictors that measure nearly the
# same thing, and the rest unrelated to them.
#
# It is checked here rather than at the draw. This is where the caller wrote the
# blocks down, so it is where a message about them is worth anything.

#' Build a block correlation matrix
#'
#' Assembles a correlation matrix out of groups of predictors that correlate with
#' each other, which is the structure worth simulating when the question is how a
#' model behaves under collinearity. Every predictor inside a block correlates
#' with every other one in it at the same value, and everything outside every
#' block correlates at `default_cor`.
#'
#' The matrix is what [simulate_regression()] and [simulate_classification()] take
#' as `cor_mat`. Its point there is that a null predictor correlated with a planted
#' one is not distinguishable from the planted one by the data alone, so its
#' estimated coefficient is drawn away from the zero it really has. That is the
#' single most common reason a coefficient table names the wrong predictor, and it
#' cannot be shown at all with independent predictors.
#'
#' @details
#' Blocks may not overlap. A predictor in two blocks would have two correlations
#' with the same partner and only one of them could be written down, so the second
#' block would silently win. Nest the smaller correlation as `default_cor` and
#' name one block instead.
#'
#' The result is checked for positive definiteness, which is the property that
#' separates a matrix of correlations from a matrix of numbers between -1 and 1.
#' Three predictors each correlated 0.9 with the other two are not jointly
#' possible, and the moment to find that out is here rather than three arguments
#' later inside a simulator.
#'
#' @param n_features Number of predictors the matrix describes, so its size.
#' @param blocks List of blocks, each a list with `features`, the indices in the
#'   block, and `cor`, the correlation they share. `list()` gives a matrix with
#'   `default_cor` everywhere off the diagonal.
#' @param default_cor Correlation between any two predictors that are not in a
#'   block together. The default of `0` leaves them independent.
#'
#' @return A symmetric `n_features` by `n_features` matrix with 1 on the diagonal.
#'
#' @seealso [simulate_regression()] and [simulate_classification()], which take the
#'   result as `cor_mat`.
#'
#' @examples
#' ## Six predictors: the first two nearly interchangeable, the next three
#' ## moderately related, and the sixth on its own.
#' cor_mat <- make_block_cor(
#'   n_features = 6,
#'   blocks = list(
#'     list(features = 1:2, cor = 0.8),
#'     list(features = 3:5, cor = 0.5)
#'   )
#' )
#' round(cor_mat, 2)
#'
#' ## Handed to a simulator, it is what makes a null predictor look like a
#' ## planted one. `x_2` has a coefficient of exactly zero and correlates 0.8
#' ## with a predictor that does not.
#' sim <- simulate_regression(n_pred = 6, beta = c(2, 0, 0, 0, 0, 0),
#'                            cor_mat = cor_mat, seed = 1)
#' sim$truth[1:2, c("predictors", "role", "beta", "max_cor_signal")]
#'
#' ## Four predictors cannot all disagree with each other, so a request for it is
#' ## refused here rather than met later by a draw that quietly ignores it.
#' try(make_block_cor(4, default_cor = -0.5))
#'
#' @export
make_block_cor <- function(n_features,
                           blocks = list(),
                           default_cor = 0) {

  n_features <- sa_check_count(n_features, "n_features", 1)
  sa_check_scalar_num(default_cor, "default_cor", -1, 1)
  if (!is.list(blocks)) {
    stop("`blocks` must be a list of blocks, each a list with `features` and ",
         "`cor`.", call. = FALSE)
  }

  cor_mat <- matrix(default_cor, nrow = n_features, ncol = n_features)
  diag(cor_mat) <- 1

  claimed <- integer(0)
  for (k in seq_along(blocks)) {
    block <- blocks[[k]]
    label <- paste0("blocks[[", k, "]]")
    if (!is.list(block) || is.null(block$features) || is.null(block$cor)) {
      stop("`", label, "` must be a list with `features` and `cor`.",
           call. = FALSE)
    }

    idx <- block$features
    if (!is.numeric(idx) || length(idx) < 2L || anyNA(idx) ||
          any(idx != trunc(idx)) || anyDuplicated(idx) > 0L) {
      stop("`", label, "$features` must be at least two distinct whole ",
           "numbers, the indices of the predictors in the block.",
           call. = FALSE)
    }
    idx <- as.integer(idx)
    if (any(idx < 1L) || any(idx > n_features)) {
      stop("`", label, "$features` indexes predictor(s) outside the ",
           n_features, " that `n_features` asks for: ",
           paste(sort(unique(idx[idx < 1L | idx > n_features])),
                 collapse = ", "), ".", call. = FALSE)
    }
    # A predictor in two blocks would need two correlations with the same
    # partner, and only the later one would survive being written down.
    overlap <- intersect(idx, claimed)
    if (length(overlap) > 0L) {
      stop("`", label, "$features` overlaps an earlier block at predictor(s) ",
           paste(sort(overlap), collapse = ", "), ". A predictor can only ",
           "carry one within-block correlation, so nest the smaller one as ",
           "`default_cor` instead.", call. = FALSE)
    }
    claimed <- c(claimed, idx)

    sa_check_scalar_num(block$cor, paste0(label, "$cor"), -1, 1)
    cor_mat[idx, idx] <- block$cor
    diag(cor_mat)[idx] <- 1
  }

  # Rejected here rather than by `chol()` inside a simulator, where the message
  # would be about a matrix the caller never wrote out.
  if (is.null(sa_sim_chol(cor_mat))) {
    stop("these blocks do not describe a possible correlation matrix: it is ",
         "not positive definite, so no data has these correlations. A block ",
         "asking for near-perfect agreement among three or more predictors is ",
         "the usual cause, and a `default_cor` of the opposite sign to the ",
         "blocks is the other.", call. = FALSE)
  }
  cor_mat
}
