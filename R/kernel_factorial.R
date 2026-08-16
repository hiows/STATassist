# The kernels of a fully crossed between-subject analysis, written to the same
# rule as `kernel_anova.R`: numbers in, a named numeric vector or a matrix out,
# no fitted model kept anywhere.
#
# Two things are different here. A factorial analysis answers on two axes at
# once, the whole model and the individual terms, and both come out of one call
# for the reason `sa_oneway_anova()` and `sa_tukey()` share a mean square error:
# the same sums of squares computed twice in two code paths is how the two ends
# of a result object come to disagree.
#
# And the arithmetic runs on the cell means rather than on the observations. A
# fully crossed model gives every row of a cell the same predictor values, so the
# residual sum of squares of any sub-model is the within-cell sum of squares plus
# the weighted residual sum of squares of the cell means, with the cell counts as
# weights. Every sum of squares here is a difference of two such residuals, so
# the within-cell part cancels and what is left is a regression of `n_cells`
# numbers rather than of `n_used` of them. It is exact rather than approximate:
# the two formulations have the same normal equations.


#' Sum-to-zero model matrix of the cells of a crossed design
#'
#' One row per cell and one column per model coefficient, the intercept first and
#' then the terms in `sa_fact_terms()` order. Sum-to-zero (`contr.sum`) coding is
#' what makes a term's columns orthogonal to the terms it does not contain in a
#' balanced design, which is why Type I, II and III agree there and why a Type III
#' sum of squares means what the unweighted marginal means say it does.
#'
#' The columns of an interaction are the elementwise products of the columns of
#' the factors it is over, which is what `model.matrix()` builds for a `:` term.
#' Built here instead so that the column order is the term order and no formula
#' has to be parsed at run time.
#'
#' @param factor_lv Named list of factors and their levels.
#' @param cells Grid of level indices, as `sa_fact_grid()` returns.
#'
#' @return List with the matrix `x`, the term index `assign` of each of its
#'   columns (`0` for the intercept) and the `terms` themselves.
#'
#' @keywords internal
#' @noRd
sa_fact_cell_matrix <- function(factor_lv, cells) {
  terms <- sa_fact_terms(names(factor_lv))
  codes <- lapply(factor_lv, function(lv) unname(stats::contr.sum(length(lv))))
  n_cells <- nrow(cells)

  x <- matrix(1, nrow = n_cells, ncol = 1L)
  assign <- 0L
  for (k in seq_along(terms)) {
    block <- matrix(1, nrow = n_cells, ncol = 1L)
    for (f in terms[[k]]) {
      cf <- codes[[f]][cells[[f]], , drop = FALSE]
      block <- block[, rep(seq_len(ncol(block)), times = ncol(cf)),
                     drop = FALSE] *
        cf[, rep(seq_len(ncol(cf)), each = ncol(block)), drop = FALSE]
    }
    x <- cbind(x, block)
    assign <- c(assign, rep(k, ncol(block)))
  }

  list(x = x, assign = assign, terms = terms)
}


#' The two models whose difference is a term's sum of squares
#'
#' Every sum of squares in an ANOVA table is a model comparison, and the three
#' types differ only in which pair of models is compared. Naming the pairs once
#' means the three types share one arithmetic path and cannot come to disagree
#' about anything except what they are meant to disagree about.
#'
#' \describe{
#'   \item{III}{Everything else stays in and the term comes out, so the sum of
#'     squares is what this term explains that no other term can. Unweighted by
#'     cell size, which is what makes it the type to score a simulation against:
#'     the planted main effect is a statement about the levels, not about how
#'     many observations happened to land in each.}
#'   \item{II}{The term is added to the model holding every term that does not
#'     contain it, so a main effect is adjusted for the other main effects but not
#'     for the interaction it is part of.}
#'   \item{I}{Sequential. Each term is added to the ones before it, so the sums of
#'     squares add up to the between-cell sum of squares exactly and the answer
#'     depends on the order the factors were declared in. This is what
#'     [stats::aov()] reports, which is what makes it the type an external check
#'     can be matched against on unbalanced data.}
#' }
#'
#' @param terms List of character vectors, one per model term.
#' @param assign Term index of each column of the model matrix.
#' @param ss_type `"III"`, `"II"` or `"I"`.
#'
#' @return List of one entry per term, each holding the `base` and `full` column
#'   indices of the two models compared.
#'
#' @keywords internal
#' @noRd
sa_fact_ss_plan <- function(terms, assign, ss_type) {
  n_terms <- length(terms)
  cols_of <- lapply(seq_len(n_terms), function(k) which(assign == k))
  intercept <- which(assign == 0L)

  lapply(seq_len(n_terms), function(k) {
    keep <- switch(
      ss_type,
      III = setdiff(seq_len(n_terms), k),
      # A term contains this one when its factors include all of them, which is
      # also true of the term itself, so it drops out without a second condition.
      II  = which(vapply(terms, function(u) !all(terms[[k]] %in% u),
                         logical(1))),
      I   = seq_len(k - 1L),
      stop("internal error: unknown `ss_type` `", ss_type, "`.", call. = FALSE)
    )
    base <- c(intercept, unlist(cols_of[keep], use.names = FALSE))
    list(base = sort(base), full = sort(c(base, cols_of[[k]])))
  })
}


#' Everything about the model that does not depend on the data
#'
#' The design matrix and the model comparisons are the same for every feature, so
#' they are settled once and handed to the kernel rather than rebuilt per feature.
#' With hundreds of features that is the difference between building one matrix
#' and building hundreds of identical ones.
#'
#' @param factor_lv Named list of factors and their levels.
#' @param cells Grid of level indices.
#' @param ss_type `"III"`, `"II"` or `"I"`.
#'
#' @return List with `x`, `assign`, `terms`, `labels`, `orders` and `ss`.
#'
#' @keywords internal
#' @noRd
sa_factorial_plan <- function(factor_lv, cells, ss_type) {
  mat <- sa_fact_cell_matrix(factor_lv, cells)
  list(
    x       = mat$x,
    assign  = mat$assign,
    terms   = mat$terms,
    labels  = sa_fact_term_labels(mat$terms),
    orders  = lengths(mat$terms),
    ss      = sa_fact_ss_plan(mat$terms, mat$assign, ss_type)
  )
}


#' Factorial analysis of variance, whole model and term by term
#'
#' The whole-model test is the one-way ANOVA that treats the cells as groups,
#' which is the same test as the F test of a fully crossed model: the crossed
#' model is the cell means model written in another basis, so it fits the same
#' values and leaves the same residuals. That is what lets one feature-wise row
#' stand for "this feature responds to the design" while the term rows say which
#' part of it responds.
#'
#' The term sums of squares are model comparisons on the weighted cell means, as
#' described at the top of this file, so `df` comes from the ranks of the two
#' matrices rather than from a formula and stays right if a design ever arrives
#' that is not of full rank.
#'
#' @param samples List of numeric vectors, one per cell in cell order, no missing
#'   values.
#' @param plan The list `sa_factorial_plan()` returns.
#'
#' @return List with `model`, the named numeric row of the whole-model test;
#'   `terms`, a matrix with one row per term; and `means`, `n`, `ms_error` and
#'   `df_error`, which the post-hoc stage reuses so that it is scaled by the same
#'   mean square error the F tests were.
#'
#' @keywords internal
#' @noRd
sa_factorial_anova <- function(samples, plan) {
  n <- lengths(samples)
  if (any(n == 0L)) {
    stop("cell(s) with no usable observation, which leaves a crossed model ",
         "with nothing to estimate there: ",
         paste(names(samples)[n == 0L], collapse = ", "), ".", call. = FALSE)
  }

  # The whole-model row is delegated rather than rewritten, so a factorial result
  # and a one-way result over the same cells cannot report different F values.
  model <- sa_oneway_anova(samples)
  names(model)[names(model) == "n_groups"] <- "n_cells"

  means <- vapply(samples, mean, numeric(1))
  ss_within <- sum(vapply(seq_along(samples), function(c) {
    sum((samples[[c]] - means[c])^2)
  }, numeric(1)))
  df_error <- model[["df2"]]
  ms_error <- ss_within / df_error
  grand <- sum(n * means) / sum(n)
  ss_total <- ss_within + sum(n * (means - grand)^2)

  # One weighted least squares problem per distinct sub-model. Several terms ask
  # for the same one under Type I and II, so the answers are kept.
  sw <- sqrt(n)
  xw <- plan$x * sw
  yw <- means * sw
  seen <- new.env(parent = emptyenv())
  fit <- function(cols) {
    key <- paste(cols, collapse = ",")
    got <- seen[[key]]
    if (is.null(got)) {
      q <- qr(xw[, cols, drop = FALSE])
      got <- list(rss = sum(qr.resid(q, yw)^2), rank = q$rank)
      assign(key, got, envir = seen)
    }
    got
  }

  columns <- c("n_used", "df", "ss", "ms", "f_stat", "df_error", "eta_sq",
               "partial_eta_sq", "pval")
  rows <- lapply(plan$ss, function(pair) {
    base <- fit(pair$base)
    full <- fit(pair$full)
    df <- full$rank - base$rank
    # Subtracting two residual sums of squares of nearly equal size can land a
    # hair below zero on a term that explains nothing at all.
    ss <- max(base$rss - full$rss, 0)
    if (df < 1L) {
      return(sa_row(n_used = sum(n), df = 0, ss = ss, ms = NA_real_,
                    f_stat = NA_real_, df_error = df_error, eta_sq = NA_real_,
                    partial_eta_sq = NA_real_, pval = NA_real_))
    }
    ms <- ss / df
    f_stat <- ms / ms_error
    sa_row(n_used         = sum(n),
           df             = df,
           ss             = ss,
           ms             = ms,
           f_stat         = f_stat,
           df_error       = df_error,
           eta_sq         = ss / ss_total,
           partial_eta_sq = ss / (ss + ss_within),
           pval           = stats::pf(f_stat, df, df_error, lower.tail = FALSE))
  })

  list(
    model    = model,
    terms    = matrix(unlist(rows, use.names = FALSE), nrow = length(rows),
                      byrow = TRUE, dimnames = list(plan$labels, columns)),
    means    = means,
    n        = n,
    ms_error = ms_error,
    df_error = df_error
  )
}


#' Tukey-Kramer comparisons of marginal means and of simple effects
#'
#' The post-hoc stage of a factorial model. Every contrast is scaled by the mean
#' square error of the whole model, which is what makes these comparisons
#' consistent with the F tests they follow rather than a second analysis of the
#' same data.
#'
#' A marginal mean is the **unweighted** mean of the cell means, not the mean of
#' the observations, so a level's mean is not pulled towards whichever combination
#' of the other factors happened to be sampled most. That is also the quantity
#' `truth_contrast` records, being `rowMeans()` over the cells of a level, so an
#' estimate here and a planted delta there are the same quantity and can be
#' compared without a correction. The weights being `1/m`, the variance of a
#' difference is `MSE * (sum(w^2 / n_c) + sum(w^2 / n_c))`, the Kramer form,
#' which reduces to Tukey's own when the cells are equal in size.
#'
#' The family is one block of contrasts rather than the whole table: the
#' studentised range is over the number of levels of the factor being compared,
#' so a marginal comparison of a four-level factor and a simple comparison of the
#' same factor within one stratum are judged against the same distribution, and
#' the p-values are family-wise within each block without further adjustment.
#'
#' @param fit The list `sa_factorial_anova()` returns.
#' @param skeleton The list `sa_fact_contrast_skeleton()` returns.
#' @param nmeans Number of means the studentised range spans, one per row of the
#'   skeleton.
#' @param rows Skeleton rows to compute, in the order they are wanted.
#' @param conf_level Confidence level for the reported intervals.
#'
#' @return Matrix of `sa_posthoc_columns()`, one row per entry of `rows`.
#'
#' @references
#' Tukey, J. W. (1949). Comparing individual means in the analysis of variance.
#' *Biometrics*, 5(2), 99-114.
#'
#' Kramer, C. Y. (1956). Extension of multiple range tests to group means with
#' unequal numbers of replications. *Biometrics*, 12(3), 307-310.
#'
#' @keywords internal
#' @noRd
sa_factorial_tukey <- function(fit, skeleton, nmeans, rows, conf_level = 0.95) {
  means <- fit$means
  n <- fit$n
  mse <- fit$ms_error
  df <- fit$df_error
  if (mse <= 0) {
    stop("the mean square error of the model is zero, so no contrast can be ",
         "scaled.", call. = FALSE)
  }

  out <- lapply(rows, function(k) {
    s1 <- skeleton$sel1[[k]]
    s2 <- skeleton$sel2[[k]]
    estimate <- mean(means[s1]) - mean(means[s2])
    variance <- mse * (sum((1 / length(s1))^2 / n[s1]) +
                         sum((1 / length(s2))^2 / n[s2]))
    # The studentised range is the range of the means over the standard error of
    # one of them, so the divisor carries a 1/2 that a t statistic does not.
    stderr <- sqrt(variance / 2)
    q_stat <- estimate / stderr
    sa_row(n1         = sum(n[s1]),
           n2         = sum(n[s2]),
           estimate   = estimate,
           stderr     = stderr,
           statistic  = q_stat,
           df         = df,
           pval       = stats::ptukey(abs(q_stat), nmeans[k], df,
                                      lower.tail = FALSE),
           lower_conf = estimate - stats::qtukey(conf_level, nmeans[k], df) *
             stderr,
           upper_conf = estimate + stats::qtukey(conf_level, nmeans[k], df) *
             stderr)[sa_posthoc_columns()]
  })

  matrix(unlist(out, use.names = FALSE), nrow = length(rows), byrow = TRUE,
         dimnames = list(NULL, sa_posthoc_columns()))
}
