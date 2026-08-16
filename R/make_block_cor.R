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
#' block correlates at `default_cor`. A block that names `against` is split in
#' two instead: each side agrees within itself and disagrees with the other side.
#'
#' The matrix is what [simulate_regression()] and [simulate_classification()] take
#' as `cor_mat`. Its point there is that a null predictor correlated with a planted
#' one is not distinguishable from the planted one by the data alone, so its
#' estimated coefficient is drawn away from the zero it really has. That is the
#' single most common reason a coefficient table names the wrong predictor, and it
#' cannot be shown at all with independent predictors.
#'
#' @details
#' Each block is its own `list()`. A single `list()` that names `features` twice
#' is one block and not two, because `$` reads the first of a repeated name and
#' nothing else, so the second pair of values would be dropped without a word.
#' Repeated and unknown names in a block are refused for that reason.
#'
#' Blocks may not overlap. A predictor in two blocks would have two correlations
#' with the same partner and only one of them could be written down, so the second
#' block would silently win. Nest the smaller correlation as `default_cor` and
#' name one block instead, or, when the second group is what the first moves
#' against, name it as `against` in one block.
#'
#' How strong a negative correlation one block can hold depends on how many
#' predictors are in it. A block of `k` predictors sharing one value is positive
#' definite only above `-1/(k - 1)`: two predictors may disagree at -0.9, three at
#' no more than -0.5, four at no more than -0.333. Three predictors cannot all
#' disagree strongly, since whichever way the third moves it agrees with one of
#' the first two.
#'
#' `against` is how a strong negative correlation is written down instead.
#' `list(features = 1:3, cor = 0.9, against = 4:6)` puts 0.9 among the first
#' three, 0.9 among the last three and -0.9 between the two sides. Splitting a
#' block by which way its predictors move leaves it positive definite for any
#' `cor` below 1 whatever its size, since it is then one factor with a sign per
#' predictor rather than a demand that everything disagree at once.
#'
#' The result is checked for positive definiteness, which is the property that
#' separates a matrix of correlations from a matrix of numbers between -1 and 1.
#' The blocks and `default_cor` are checked one at a time first, so that a value
#' no block of that size could hold is named as such, and the eigenvalue of the
#' assembled matrix is what reports a `default_cor` the blocks cannot sit inside.
#'
#' @param n_features Number of predictors the matrix describes, so its size.
#' @param blocks List of blocks, each a list with `features`, the indices in the
#'   block, and `cor`, the correlation they share. A block may also name
#'   `against`, further indices that correlate at `cor` among themselves and at
#'   `-cor` with the ones in `features`. `list()` gives a matrix with
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
#' ## Three predictors moving one way and three the other. `against` carries the
#' ## sign, so -0.9 between the sides is available at any block size.
#' opposed <- make_block_cor(
#'   n_features = 6,
#'   blocks = list(list(features = 1:3, cor = 0.9, against = 4:6))
#' )
#' round(opposed, 2)
#'
#' ## Three predictors cannot all disagree at -0.6, and the limit for a block of
#' ## three is named rather than left to the matrix.
#' try(make_block_cor(6, list(list(features = 1:3, cor = -0.6))))
#'
#' ## Each block needs its own `list()`. One `list()` naming `features` twice is
#' ## a single block whose second value R would never read.
#' try(make_block_cor(6, list(list(features = 1:3, cor = 0.9,
#'                                 features = 4:6, cor = -0.4))))
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
  # With no block to break the pattern the whole matrix holds one value off the
  # diagonal, so the bound on it is exact and belongs beside the argument. With
  # blocks it is the combination that can fail, and the eigenvalue at the end is
  # what says so.
  if (length(blocks) == 0L && n_features > 1L) {
    sa_block_shared_bound(default_cor, "default_cor", n_features,
                          paste(n_features, "predictors"))
  }

  cor_mat <- matrix(default_cor, nrow = n_features, ncol = n_features)
  diag(cor_mat) <- 1

  claimed <- integer(0)
  for (k in seq_along(blocks)) {
    block <- blocks[[k]]
    label <- paste0("blocks[[", k, "]]")
    if (!is.list(block)) {
      stop("`", label, "` must be a list with `features` and `cor`.",
           call. = FALSE)
    }
    sa_block_names(block, label)
    if (is.null(block$features) || is.null(block$cor)) {
      stop("`", label, "` must be a list with `features` and `cor`.",
           call. = FALSE)
    }

    # A side of one predictor is meaningful only against another side, so the two
    # indices are enough for a split block and two for `features` otherwise.
    two_sided <- !is.null(block$against)
    feat <- sa_block_index(block$features, paste0(label, "$features"),
                           n_features, if (two_sided) 1L else 2L)
    agn <- if (two_sided) {
      sa_block_index(block$against, paste0(label, "$against"), n_features, 1L)
    } else {
      integer(0)
    }
    on_both <- intersect(feat, agn)
    if (length(on_both) > 0L) {
      stop("`", label, "` names predictor(s) ",
           paste(sort(on_both), collapse = ", "), " in both `features` and ",
           "`against`, and a predictor cannot move against itself.",
           call. = FALSE)
    }

    idx <- c(feat, agn)
    signs <- rep(c(1, -1), c(length(feat), length(agn)))

    # A predictor in two blocks would need two correlations with the same
    # partner, and only the later one would survive being written down. Asked
    # before the value, since how negative the value may be is a question about
    # the size of the block, and the block is not settled until this holds.
    overlap <- intersect(idx, claimed)
    if (length(overlap) > 0L) {
      stop("`", label, "` overlaps an earlier block at predictor(s) ",
           paste(sort(overlap), collapse = ", "), ". A predictor can only ",
           "carry one within-block correlation, so nest the smaller one as ",
           "`default_cor`, or, when this block is what the earlier one moves ",
           "against, name it as `against` in that block instead.",
           call. = FALSE)
    }
    claimed <- c(claimed, idx)

    cor_label <- paste0(label, "$cor")
    sa_check_scalar_num(block$cor, cor_label, -1, 1)
    if (two_sided) {
      # `against` is what carries the sign here. A negative `cor` would make each
      # side disagree within itself and agree across, which is the limit this
      # argument exists to lift, written backwards.
      if (block$cor <= 0) {
        stop("`", cor_label, "` must be above 0 when `against` is given, but ",
             "is ", block$cor, ". Each side agrees at `cor` and disagrees with ",
             "the other side at -`cor`, so it is `against` that makes a ",
             "correlation negative.", call. = FALSE)
      }
      if (block$cor >= 1) {
        stop("`", cor_label, "` of ", block$cor, " puts each side of the block ",
             "at perfect agreement, which is one variable repeated rather than ",
             "several, so the matrix is singular rather than a correlation ",
             "matrix. Use a value below 1.", call. = FALSE)
      }
    } else {
      sa_block_shared_bound(block$cor, cor_label, length(idx),
                            paste("the", length(idx),
                                  "predictors of the block"))
    }

    cor_mat[idx, idx] <- block$cor * outer(signs, signs)
    diag(cor_mat)[idx] <- 1
  }

  # Rejected here rather than by `chol()` inside a simulator, where the message
  # would be about a matrix the caller never wrote out. Every block holds on its
  # own by this point, so what is left to catch is how they meet.
  if (is.null(sa_sim_chol(cor_mat))) {
    min_eigen <- min(eigen(cor_mat, symmetric = TRUE,
                           only.values = TRUE)$values)
    stop("these blocks do not describe a possible correlation matrix: its ",
         "smallest eigenvalue is ", signif(min_eigen, 3), ", where a ",
         "correlation matrix has none below 0, so no data has these ",
         "correlations. Every block holds on its own, so what does not is ",
         "`default_cor` beside them, most often a value of the opposite sign ",
         "to the blocks.", call. = FALSE)
  }
  cor_mat
}


#' Check the names a block holds
#'
#' The trap this exists for is `list(features = 1:3, cor = 0.9, features = 4:6,
#' cor = -0.4)`, which reads as two blocks and is one. `$` returns the first of a
#' repeated name, so the second pair would be dropped and the matrix would come
#' back holding a third of what was asked for, with nothing said.
#'
#' @keywords internal
#' @noRd
sa_block_names <- function(block, label) {
  if (length(block) == 0L) {
    return(invisible(block))
  }
  nm <- names(block)
  if (is.null(nm) || anyNA(nm) || !all(nzchar(nm))) {
    stop("`", label, "` must name every element it holds: `features`, `cor`, ",
         "and `against` when its predictors do not all move the same way.",
         call. = FALSE)
  }
  repeated <- unique(nm[duplicated(nm)])
  if (length(repeated) > 0L) {
    stop("`", label, "` names ",
         paste0("`", repeated, "`", collapse = " and "), " more than once. ",
         "`$` reads the first of a repeated name and nothing else, so every ",
         "later value would be dropped without a word. Several blocks are ",
         "several `list()`s: blocks = list(list(features = 1:3, cor = 0.9), ",
         "list(features = 4:6, cor = -0.4)).", call. = FALSE)
  }
  unknown <- setdiff(nm, c("features", "cor", "against"))
  if (length(unknown) > 0L) {
    stop("`", label, "` holds ", paste0("`", unknown, "`", collapse = ", "),
         ", which a block has no use for. A block is `features` and `cor`, and ",
         "`against` when its predictors do not all move the same way.",
         call. = FALSE)
  }
  invisible(block)
}


#' Check one side of a block's predictor indices
#'
#' @param min_len Indices a side must hold. Two for a block that is not split,
#'   since a correlation needs a pair, and one for each side of a split block,
#'   since the pair is then across the two sides.
#'
#' @return The indices as integers.
#'
#' @keywords internal
#' @noRd
sa_block_index <- function(idx, label, n_features, min_len) {
  ok <- is.numeric(idx) && length(idx) >= min_len && !anyNA(idx) &&
    all(idx == trunc(idx)) && anyDuplicated(idx) == 0L
  if (!ok) {
    detail <- if (min_len == 1L) {
      paste("one or more distinct whole numbers, the indices of the predictors",
            "on that side of the block.")
    } else {
      paste("at least two distinct whole numbers, the indices of the",
            "predictors in the block.")
    }
    stop("`", label, "` must be ", detail, call. = FALSE)
  }
  idx <- as.integer(idx)
  outside <- sort(unique(idx[idx < 1L | idx > n_features]))
  if (length(outside) > 0L) {
    stop("`", label, "` indexes predictor(s) outside the ", n_features,
         " that `n_features` asks for: ", paste(outside, collapse = ", "), ".",
         call. = FALSE)
  }
  idx
}


#' Check that one shared correlation could hold among `k` predictors
#'
#' A `k` by `k` matrix with 1 on the diagonal and one value everywhere else has
#' eigenvalues `1 - value` and `1 + (k - 1) * value`, so it is a correlation
#' matrix exactly for `value` in `(-1/(k - 1), 1)`. Both ends are worth a
#' sentence of their own: the top is one predictor written twice, and the bottom
#' is the reason `against` exists.
#'
#' @param among Noun phrase for the predictors sharing the value, since the same
#'   bound is what `default_cor` and an unsplit block are each held to.
#'
#' @keywords internal
#' @noRd
sa_block_shared_bound <- function(value, label, k, among) {
  if (value >= 1) {
    stop("`", label, "` of ", value, " puts ", among, " at perfect agreement, ",
         "which is one variable repeated rather than several, so the matrix is ",
         "singular rather than a correlation matrix. Use a value below 1.",
         call. = FALSE)
  }
  bound <- -1 / (k - 1)
  if (value <= bound) {
    stop("`", label, "` of ", value, " is not possible among ", among,
         ": one value shared by every pair holds only above ", signif(bound, 3),
         ", since they cannot all disagree with each other at once. Name the ",
         "ones that move the other way as `against` in a block instead, which ",
         "carries the sign and has no such limit.", call. = FALSE)
  }
  invisible(value)
}
