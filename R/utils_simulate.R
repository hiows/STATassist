# Internal helpers shared by the supervised learning simulators. The two of them
# differ only in what they do with the linear predictor: one adds noise to it and
# one runs it through a logistic link. Everything before that is the same
# question twice over — which coefficients were planted, what the predictors look
# like, which rows belong to which subject — so it is answered once here, the way
# `utils_model.R` answers the input resolution once for both fitting functions.
#
# The correlated draw is done with `chol()` rather than `MASS::mvrnorm()`. It is
# three lines either way, and doing it here keeps `Imports` at the one package
# the machine learning family already needed. It is also the stricter of the two:
# a correlation matrix that no data could have is rejected instead of quietly
# being projected onto one that could.
#
# Every argument is checked before any of them is drawn on. A simulator that
# consumed the random stream and then rejected its own arguments would give a
# different data set to the seed that fixed it, depending on which call failed
# first.

#' Check a numeric argument that may be given once or once per predictor
#'
#' The alternative would be to draw the means and standard deviations from a
#' range, the way the expression simulators draw theirs. It is not done here
#' because `beta` is a coefficient per unit of its predictor, so the size of the
#' effect that is planted depends on the spread of the column it is planted on. A
#' spread that was drawn at random would leave the signal-to-noise ratio of the
#' simulation unreadable from its arguments.
#'
#' @keywords internal
#' @noRd
sa_sim_recycle <- function(x, n, arg, lower = -Inf) {
  if (!is.numeric(x) || !length(x) %in% c(1L, n) || !all(is.finite(x))) {
    stop("`", arg, "` must be a finite numeric vector of length 1 or ", n,
         ", the number of numeric predictors.", call. = FALSE)
  }
  if (any(x < lower)) {
    stop("`", arg, "` must not go below ", lower, ".", call. = FALSE)
  }
  rep_len(as.numeric(x), n)
}


#' Check a correlation matrix and return the factor the draw needs
#'
#' The three properties checked first are the ones that make a matrix a
#' correlation matrix, and the fourth is the one that makes it a possible one. A
#' symmetric matrix with a unit diagonal can still describe no data at all: ask
#' for three predictors each correlated 0.9 with the other two and no set of
#' vectors satisfies it. `chol()` is what finds that out, so its factor is
#' returned rather than recomputed at the draw.
#'
#' @param cor_mat The matrix as received, or `NULL` for independence.
#' @param n_pred Number of numeric predictors it has to describe.
#' @param arg Argument name to name in the error.
#'
#' @return The upper triangular Cholesky factor.
#'
#' @keywords internal
#' @noRd
sa_sim_cor_root <- function(cor_mat, n_pred, arg = "cor_mat") {
  if (is.null(cor_mat)) {
    return(diag(n_pred))
  }
  if (!is.matrix(cor_mat) || !is.numeric(cor_mat) ||
        !identical(dim(cor_mat), c(n_pred, n_pred))) {
    stop("`", arg, "` must be a numeric ", n_pred, " x ", n_pred,
         " matrix, one row and column per numeric predictor.", call. = FALSE)
  }
  if (!all(is.finite(cor_mat))) {
    stop("`", arg, "` must not contain missing or non-finite values.",
         call. = FALSE)
  }
  if (!isSymmetric(unname(cor_mat))) {
    stop("`", arg, "` must be symmetric.", call. = FALSE)
  }
  if (!all(diag(cor_mat) == 1)) {
    stop("`", arg, "` must have 1 on its diagonal, since a variable is ",
         "perfectly correlated with itself.", call. = FALSE)
  }
  if (any(abs(cor_mat) > 1)) {
    stop("`", arg, "` holds correlation(s) outside [-1, 1].", call. = FALSE)
  }

  root <- sa_sim_chol(cor_mat)
  if (is.null(root)) {
    stop("`", arg, "` is not positive definite, so no data has these ",
         "correlations. Build it with make_block_cor(), which says which of ",
         "the blocks cannot hold.", call. = FALSE)
  }
  root
}


#' Factorise a correlation matrix, or say that it cannot be
#'
#' A symmetric matrix with a unit diagonal and no entry outside -1 and 1 can
#' still describe no data at all: ask for three predictors each correlated 0.9
#' with the other two and no set of vectors satisfies it. `chol()` is what finds
#' that out, and the factor it produces is also what the draw needs, so the two
#' questions are answered by one call.
#'
#' The failure comes back as `NULL` rather than as an error, because the sentence
#' worth printing depends on whether the caller wrote the matrix down or asked
#' for blocks that add up to it.
#'
#' @keywords internal
#' @noRd
sa_sim_chol <- function(cor_mat) {
  tryCatch(chol(unname(cor_mat)), error = function(e) NULL)
}


#' Draw correlated normal predictors
#'
#' @param n Rows to draw.
#' @param value_mean,value_sd One entry per predictor.
#' @param root Cholesky factor of the correlation matrix.
#'
#' @return An `n` by `length(value_mean)` matrix.
#'
#' @keywords internal
#' @noRd
sa_sim_mvnorm <- function(n, value_mean, value_sd, root) {
  n_pred <- length(value_mean)
  z <- matrix(stats::rnorm(n * n_pred), nrow = n, ncol = n_pred)
  # The columns are scaled after the rotation rather than by assembling a
  # covariance matrix, since `diag()` of a single standard deviation is the
  # identity of that size rather than a 1 x 1 matrix holding it.
  out <- z %*% root
  out * rep(value_sd, each = n) + rep(value_mean, each = n)
}


#' Settle how many numeric predictors there are, and how their coefficients come
#'
#' Two ways of saying what the coefficients are, and they are not alternatives
#' that need reconciling: one plants a known number of them and leaves the rest at
#' exactly zero, and the other states all of them. Naming both is refused rather
#' than resolved, the same way `simulate_multiple_groups()` refuses a `group_lv`
#' and an `n_treat` that count differently.
#'
#' Nothing is drawn here. The counts and the lengths are settled, and
#' `sa_sim_plant_beta()` is what turns them into numbers once every other
#' argument has been checked too.
#'
#' @param explicit Names of the arguments the caller actually supplied, which is
#'   how a default is told apart from a value that was asked for.
#'
#' @return List with `n_pred`, `n_pos`, `n_neg`, the `beta` that was supplied or
#'   `NULL`, and `value_mean` / `value_sd` recycled to `n_pred`.
#'
#' @keywords internal
#' @noRd
sa_sim_pred_spec <- function(n_pred, beta, n_pos, n_neg, beta_range,
                             value_mean, value_sd, explicit) {
  if (is.null(beta)) {
    n_pred <- sa_check_count(n_pred, "n_pred", 1)
    n_pos <- sa_check_count(n_pos, "n_pos")
    n_neg <- sa_check_count(n_neg, "n_neg")
    if (n_pos + n_neg > n_pred) {
      stop("`n_pos` + `n_neg` is ", n_pos + n_neg, ", which is more ",
           "coefficients than the ", n_pred, " numeric predictor(s) that ",
           "`n_pred` asks for.", call. = FALSE)
    }
  } else {
    clash <- intersect(explicit, c("n_pos", "n_neg"))
    if (length(clash) > 0L) {
      stop("`beta` states every coefficient, so there is nothing left for ",
           paste0("`", clash, "`", collapse = " and "), " to plant. Drop one ",
           "of the two.", call. = FALSE)
    }
    if (!is.numeric(beta) || length(beta) == 0L || !all(is.finite(beta))) {
      stop("`beta` must be a finite numeric vector, one coefficient per ",
           "numeric predictor and no intercept among them.", call. = FALSE)
    }
    # `beta` holds one entry per predictor, which makes its length the number of
    # them, exactly as the length of `n_treat` is the number of treatment groups
    # in `simulate_multiple_groups()`.
    if ("n_pred" %in% explicit && length(beta) != n_pred) {
      stop("`n_pred` asks for ", n_pred, " numeric predictor(s) but `beta` ",
           "gives ", length(beta), " coefficient(s). The intercept is not one ",
           "of them: it is `intercept` for a regression and `event_rate` for a ",
           "classification.", call. = FALSE)
    }
    n_pred <- length(beta)
    n_pos <- 0L
    n_neg <- 0L
  }

  list(
    n_pred     = n_pred,
    n_pos      = n_pos,
    n_neg      = n_neg,
    beta       = beta,
    value_mean = sa_sim_recycle(value_mean, n_pred, "value_mean"),
    value_sd   = sa_sim_recycle(value_sd, n_pred, "value_sd", 0)
  )
}


#' Turn the settled counts into coefficients
#'
#' A planted coefficient lands on a predictor drawn at random, but how many are
#' positive and how many are negative is a function of the arguments alone.
#' Drawing the signs instead would make those two counts move with the seed, which
#' is the kind of thing about a simulation that should not have to be looked up.
#'
#' @return List with `beta` and `direction`, both of length `n_pred`.
#'
#' @keywords internal
#' @noRd
sa_sim_plant_beta <- function(spec, beta_range) {
  if (!is.null(spec$beta)) {
    coefs <- as.numeric(spec$beta)
    return(list(
      beta = coefs,
      direction = ifelse(coefs > 0, "up", ifelse(coefs < 0, "down", "none"))
    ))
  }

  coefs <- numeric(spec$n_pred)
  direction <- rep("none", spec$n_pred)
  if (spec$n_pos + spec$n_neg > 0L) {
    # Head and tail of one shuffled draw, as in the expression simulators. Taking
    # the negative set as the complement of the positive one would return every
    # predictor rather than none when the positive set is empty.
    picked <- sample.int(spec$n_pred, spec$n_pos + spec$n_neg)
    pos_idx <- utils::head(picked, spec$n_pos)
    neg_idx <- utils::tail(picked, spec$n_neg)
    direction[pos_idx] <- "up"
    direction[neg_idx] <- "down"
    coefs[pos_idx] <- stats::runif(spec$n_pos, beta_range[1], beta_range[2])
    coefs[neg_idx] <- -stats::runif(spec$n_neg, beta_range[1], beta_range[2])
  }
  list(beta = coefs, direction = direction)
}


#' Work out how many subjects there are and how many rows each one carries
#'
#' `n_per_subject` carries one row count per subject, which makes its length the
#' number of subjects, the same rule `n_treat` follows in
#' `simulate_multiple_groups()`. `n_samples` says the same thing from the other
#' side, so the two are settled together rather than in two passes that could
#' disagree, and a single count has an obvious number of subjects to be spread
#' over as soon as `n_samples` says how many rows there are in all.
#'
#' @param use_default_n Result of `missing(n_samples)` in the calling function.
#'
#' @return List with `sizes`, one row count per subject or `NULL` when there are
#'   no subjects, and the resulting `n_samples`.
#'
#' @keywords internal
#' @noRd
sa_sim_subject_sizes <- function(n_samples, n_per_subject, use_default_n) {
  if (is.null(n_per_subject)) {
    return(list(sizes = NULL,
                n_samples = sa_check_count(n_samples, "n_samples", 2)))
  }
  if (!is.numeric(n_per_subject) || length(n_per_subject) == 0L) {
    stop("`n_per_subject` must be one or more row counts, one per subject, or ",
         "NULL for one row per subject.", call. = FALSE)
  }

  if (length(n_per_subject) == 1L) {
    n_samples <- sa_check_count(n_samples, "n_samples", 2)
    per <- sa_check_count(n_per_subject, "n_per_subject", 1)
    if (n_samples %% per != 0L) {
      stop("`n_per_subject` = ", per, " does not divide the ", n_samples,
           " row(s) `n_samples` asks for. Pass a row count per subject, such ",
           "as `n_per_subject = rep(", per, ", ", n_samples %/% per, ")`.",
           call. = FALSE)
    }
    sizes <- rep(per, n_samples %/% per)
  } else {
    sizes <- vapply(seq_along(n_per_subject), function(k) {
      sa_check_count(n_per_subject[k], paste0("n_per_subject[", k, "]"), 1)
    }, integer(1))
    total <- sum(sizes)
    # The counts already say how many rows there are, so a default `n_samples`
    # has nothing to add. One that was asked for and disagrees is not guessed at.
    if (!use_default_n && sa_check_count(n_samples, "n_samples", 2) != total) {
      stop("`n_per_subject` gives ", length(sizes), " subject(s) holding ",
           total, " row(s) in all, but `n_samples` asks for ", n_samples,
           ". Drop one of the two.", call. = FALSE)
    }
    n_samples <- total
  }

  if (length(sizes) < 2L) {
    stop("`n_per_subject` describes ", length(sizes), " subject(s), and a ",
         "split taken over subjects needs at least 2.", call. = FALSE)
  }
  list(sizes = as.integer(sizes), n_samples = as.integer(n_samples))
}


#' Hand out factor levels in counts that do not depend on the seed
#'
#' A permutation of a balanced vector rather than a draw with replacement, so that
#' which unit gets which level is random while how many of each there are is not.
#' `rep_len()` leaves the counts differing by at most one.
#'
#' @keywords internal
#' @noRd
sa_sim_balanced_levels <- function(n, levels) {
  factor(levels[sample(rep_len(seq_along(levels), n))], levels = levels)
}


#' Plant an offset on every factor level beyond the reference
#'
#' The reference level carries no offset, because it is what the intercept absorbs
#' and what the other levels are contrasts against. The magnitudes are drawn from
#' the same range the numeric coefficients use, and the signs alternate rather
#' than being drawn, for the reason the counts of positive and negative
#' coefficients are not drawn either.
#'
#' @keywords internal
#' @noRd
sa_sim_factor_offsets <- function(factor_lv, beta_range) {
  k <- length(factor_lv)
  offsets <- stats::setNames(numeric(k), factor_lv)
  signs <- rep_len(c(1, -1), k - 1L)
  offsets[-1L] <- signs * stats::runif(k - 1L, beta_range[1], beta_range[2])
  offsets
}


#' Punch holes in the predictors after the outcome has been computed
#'
#' The outcome is generated from the complete predictors and the holes are made
#' afterwards, which is what missing completely at random means: the value was
#' there and doing its work, and it is the record of it that is gone. Computing
#' the outcome from the holed frame instead would make the missingness part of the
#' truth rather than something the analysis has to survive.
#'
#' Only the numeric predictors are holed. A hole in the factor predictor would
#' make it unusable as the stratifier of a split, and a hole in a constant
#' predictor would stop it being constant, since `unique()` counts `NA` as one of
#' the values a column takes.
#'
#' The number of cells is a function of `p_missing` and the size of the frame;
#' only which cells they are is drawn.
#'
#' @keywords internal
#' @noRd
sa_sim_mask_missing <- function(x, p_missing) {
  if (p_missing == 0 || ncol(x) == 0L) {
    return(x)
  }
  n_row <- nrow(x)
  n_na <- round(p_missing * n_row * ncol(x))
  if (n_na == 0L) {
    return(x)
  }

  at <- sample.int(n_row * ncol(x), n_na)
  rows <- ((at - 1L) %% n_row) + 1L
  cols <- ((at - 1L) %/% n_row) + 1L
  for (j in unique(cols)) {
    x[[j]][rows[cols == j]] <- NA
  }
  x
}


#' Find the intercept that gives the requested event rate
#'
#' Solved on the linear predictor that was actually drawn rather than on its
#' expectation, so the rate the data comes out with is the rate that was asked for
#' up to the Bernoulli draw itself. `plogis()` is increasing in the intercept, so
#' the root is unique and the only work is finding a bracket that contains it.
#'
#' @param eta Linear predictor without an intercept, one entry per row.
#' @param event_rate Target proportion of events, strictly between 0 and 1.
#'
#' @keywords internal
#' @noRd
sa_sim_solve_intercept <- function(eta, event_rate) {
  gap <- function(a) mean(stats::plogis(a + eta)) - event_rate

  lower <- -1
  upper <- 1
  while (gap(lower) > 0 && lower > -1e4) lower <- lower * 2
  while (gap(upper) < 0 && upper < 1e4) upper <- upper * 2
  if (gap(lower) > 0 || gap(upper) < 0) {
    stop("no intercept gives an event rate of ", event_rate,
         " on these predictors. A rate this far from a half needs a smaller ",
         "`beta_range` or a smaller `subject_sd`.", call. = FALSE)
  }
  stats::uniroot(gap, c(lower, upper), tol = .Machine$double.eps^0.5)$root
}


#' Build the predictors, the subjects and the linear predictor they imply
#'
#' Everything the two simulators share. What comes back is a design and not yet a
#' data set: the outcome is the one thing each of them adds for itself, since a
#' continuous outcome is this linear predictor plus noise and a class is a draw
#' from the logistic function of it.
#'
#' The subject offset is inside `eta` rather than added to the outcome afterwards,
#' so that in a classification it moves the probability of the class rather than
#' the class itself. Either way it is the between-subject variation that makes a
#' row-wise split leak: a subject seen in training is partly known before its test
#' rows are read, and `subject_sd` is how much of it is.
#'
#' Factor predictors are drawn per subject rather than per row when there are
#' subjects. A subject attribute is what a factor predictor usually is in a
#' repeated-measures design, and it is also the only thing in the frame that can
#' stratify a split taken over subjects, since a stratifier has to be constant
#' within a unit.
#'
#' The numeric predictors are split between a subject level and a row level by
#' `subject_share`, which is what makes two rows of one subject resemble each
#' other rather than merely share an outcome offset. Without that resemblance a
#' row-wise split has nothing to give away: a model cannot recognise a subject it
#' was trained on if its rows look like anyone else's. The two parts are drawn
#' through the same correlation factor and their variances add to `value_sd^2`, so
#' `subject_share` moves the intraclass correlation of a column without moving its
#' distribution.
#'
#' @return List with the predictor frame `x`, the `predictors` in it and the three
#'   kinds of them by name, the planted `beta` and `direction`, `eta` without an
#'   intercept, the `subject` labels or `NULL`, `subject_offset`, `n_samples`,
#'   `sizes`, and the per-predictor and per-term truth tables.
#'
#' @keywords internal
#' @noRd
sa_sim_supervised_design <- function(n_samples,
                                     n_pred,
                                     beta,
                                     n_pos,
                                     n_neg,
                                     beta_range,
                                     value_mean,
                                     value_sd,
                                     cor_mat,
                                     n_factor_pred,
                                     factor_lv,
                                     n_constant_pred,
                                     p_missing,
                                     n_per_subject,
                                     subject_sd,
                                     subject_share,
                                     pred_prefix,
                                     explicit,
                                     use_default_n) {

  spec <- sa_sim_pred_spec(n_pred, beta, n_pos, n_neg, beta_range, value_mean,
                           value_sd, explicit)
  n_pred <- spec$n_pred
  # Checked whether or not it planted the coefficients: the factor offsets are
  # drawn from it either way, so a rejected value must not depend on `beta`.
  sa_check_range(beta_range, "beta_range", 0)
  root <- sa_sim_cor_root(cor_mat, n_pred)

  n_factor_pred <- sa_check_count(n_factor_pred, "n_factor_pred")
  n_constant_pred <- sa_check_count(n_constant_pred, "n_constant_pred")
  sa_check_scalar_num(p_missing, "p_missing", 0, 1, upper_open = TRUE)
  sa_check_scalar_num(subject_sd, "subject_sd", 0)
  sa_check_scalar_num(subject_share, "subject_share", 0, 1)
  if (n_factor_pred > 0L &&
        (!is.character(factor_lv) || length(factor_lv) < 2L ||
           anyNA(factor_lv) || anyDuplicated(factor_lv) > 0L)) {
    stop("`factor_lv` must be at least two distinct non-missing level names, ",
         "the first being the reference.", call. = FALSE)
  }
  if (!is.character(pred_prefix) || length(pred_prefix) != 1L ||
        is.na(pred_prefix) || !nzchar(pred_prefix)) {
    stop("`pred_prefix` must be a single non-empty string.", call. = FALSE)
  }

  subjects <- sa_sim_subject_sizes(n_samples, n_per_subject, use_default_n)
  sizes <- subjects$sizes
  n_samples <- subjects$n_samples
  n_unit <- if (is.null(sizes)) n_samples else length(sizes)

  numeric_pred <- paste0(pred_prefix, "_", seq_len(n_pred))
  factor_pred <- if (n_factor_pred > 0L) {
    paste0(pred_prefix, "_cat_", seq_len(n_factor_pred))
  } else {
    character(0)
  }
  constant_pred <- if (n_constant_pred > 0L) {
    paste0(pred_prefix, "_const_", seq_len(n_constant_pred))
  } else {
    character(0)
  }

  # Nothing above this line has drawn anything.
  planted <- sa_sim_plant_beta(spec, beta_range)
  values <- if (is.null(sizes)) {
    sa_sim_mvnorm(n_samples, spec$value_mean, spec$value_sd, root)
  } else {
    between <- sa_sim_mvnorm(length(sizes), spec$value_mean,
                             spec$value_sd * sqrt(subject_share), root)
    within <- sa_sim_mvnorm(n_samples, rep(0, n_pred),
                            spec$value_sd * sqrt(1 - subject_share), root)
    between[rep(seq_along(sizes), times = sizes), , drop = FALSE] + within
  }
  x <- as.data.frame(values)
  names(x) <- numeric_pred
  eta <- drop(values %*% planted$beta)

  offsets <- stats::setNames(vector("list", n_factor_pred), factor_pred)
  for (k in seq_len(n_factor_pred)) {
    level <- sa_sim_balanced_levels(n_unit, factor_lv)
    if (!is.null(sizes)) {
      level <- rep(level, times = sizes)
    }
    offsets[[k]] <- sa_sim_factor_offsets(factor_lv, beta_range)
    x[[factor_pred[k]]] <- level
    eta <- eta + offsets[[k]][as.character(level)]
  }

  for (nm in constant_pred) {
    x[[nm]] <- rep(1, n_samples)
  }

  subject <- NULL
  subject_offset <- numeric(n_samples)
  if (!is.null(sizes)) {
    subject <- rep(paste0("subject_", seq_along(sizes)), times = sizes)
    subject_offset <- rep(stats::rnorm(length(sizes), 0, subject_sd),
                          times = sizes)
  }
  eta <- unname(eta + subject_offset)

  x[numeric_pred] <- sa_sim_mask_missing(x[numeric_pred], p_missing)
  rownames(x) <- NULL

  list(
    x              = x,
    predictors     = c(numeric_pred, factor_pred, constant_pred),
    numeric_pred   = numeric_pred,
    factor_pred    = factor_pred,
    constant_pred  = constant_pred,
    beta           = planted$beta,
    direction      = planted$direction,
    offsets        = offsets,
    eta            = eta,
    subject        = subject,
    subject_offset = subject_offset,
    n_samples      = n_samples,
    sizes          = sizes,
    truth          = sa_sim_truth_pred(planted, spec, numeric_pred,
                                       factor_pred, constant_pred, cor_mat),
    truth_term     = sa_sim_truth_term(planted$beta, numeric_pred, offsets)
  )
}


#' Per predictor answer, and what it is up against
#'
#' `max_cor_signal` is the reason `make_block_cor()` exists. A null predictor
#' correlated with a planted one is the case where a coefficient of exactly zero
#' is estimated well away from zero, and looking the correlation up accounts for
#' that rather than leaving it as a false positive with no explanation.
#'
#' A factor predictor has one offset per level rather than one coefficient, so its
#' `beta` is missing here and its answer is in `truth_term`.
#'
#' @keywords internal
#' @noRd
sa_sim_truth_pred <- function(planted, spec, numeric_pred, factor_pred,
                              constant_pred, cor_mat) {
  n_pred <- spec$n_pred
  signal <- which(planted$beta != 0)
  cors <- if (is.null(cor_mat)) diag(n_pred) else abs(unname(cor_mat))
  max_cor <- vapply(seq_len(n_pred), function(i) {
    others <- setdiff(signal, i)
    if (length(others) == 0L) 0 else max(cors[i, others])
  }, numeric(1))

  n_other <- length(factor_pred) + length(constant_pred)
  data.frame(
    predictors     = c(numeric_pred, factor_pred, constant_pred),
    role           = c(ifelse(planted$beta == 0, "null", "signal"),
                       rep("factor", length(factor_pred)),
                       rep("constant", length(constant_pred))),
    beta           = c(planted$beta,
                       rep(NA_real_, length(factor_pred)),
                       rep(0, length(constant_pred))),
    direction      = c(planted$direction,
                       rep(NA_character_, length(factor_pred)),
                       rep("none", length(constant_pred))),
    value_mean     = c(spec$value_mean, rep(NA_real_, n_other)),
    value_sd       = c(spec$value_sd, rep(NA_real_, n_other)),
    max_cor_signal = c(max_cor, rep(NA_real_, n_other)),
    stringsAsFactors = FALSE
  )
}


#' Per model term answer, in the row order the coefficient table follows
#'
#' The predictors that were passed in and the terms that come back are not the
#' same list: a factor with `k` levels becomes `k - 1` terms named after the level
#' each stands for, and a predictor that takes one value becomes no term at all.
#' So the table that scores `coefficients` is built on the term axis rather than
#' reindexed from the predictor axis afterwards, which is why `truth_contrast`
#' exists beside `truth` in `simulate_multiple_groups()`.
#'
#' The intercept row is added by `sa_sim_add_intercept()`, since a regression and a
#' classification arrive at their intercept by different routes.
#'
#' @keywords internal
#' @noRd
sa_sim_truth_term <- function(beta, numeric_pred, offsets) {
  terms <- numeric_pred
  values <- beta
  predictors <- numeric_pred

  for (nm in names(offsets)) {
    lv <- names(offsets[[nm]])
    # `lm()` names a dummy column by pasting the level onto the column name, so
    # the term is predicted here rather than read back off a fit.
    terms <- c(terms, paste0(nm, lv[-1L]))
    values <- c(values, unname(offsets[[nm]][-1L]))
    predictors <- c(predictors, rep(nm, length(lv) - 1L))
  }

  data.frame(
    terms      = terms,
    predictors = predictors,
    beta       = values,
    stringsAsFactors = FALSE
  )
}


#' Put the intercept at the top of the term answer
#'
#' @keywords internal
#' @noRd
sa_sim_add_intercept <- function(truth_term, intercept) {
  out <- rbind(
    data.frame(terms = "(Intercept)", predictors = NA_character_,
               beta = intercept, stringsAsFactors = FALSE),
    truth_term
  )
  rownames(out) <- NULL
  out
}


#' What a split of this data set should be told to preserve
#'
#' A stratifier has to be constant within a sampling unit, since a unit goes to
#' one side of the split as a whole. That rules the outcome out for a regression
#' measured repeatedly, because it varies from row to row within a subject, and
#' the subject-level factor predictor is what is left. A classification has no
#' such problem: a subject is a case or a control as a whole.
#'
#' @param stratify_outcome Whether the outcome is constant within a subject.
#'
#' @keywords internal
#' @noRd
sa_sim_split_args <- function(data, design, stratify_outcome) {
  stratified <- if (stratify_outcome || is.null(design$subject)) {
    "y"
  } else if (length(design$factor_pred) > 0L) {
    design$factor_pred[1L]
  } else {
    NULL
  }

  list(
    data       = data,
    stratified = stratified,
    id         = if (is.null(design$subject)) NULL else "subject"
  )
}
