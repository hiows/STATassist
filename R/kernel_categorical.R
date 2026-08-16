# Contingency table kernels, written to the same rule as `kernel_anova.R`: plain
# input in, a named numeric vector out, no fitted object kept anywhere and
# nothing said to the user.
#
# That last part is a rule and not an accident. A kernel is a function of its
# arguments, so anything it would want to tell the caller -- that a correction was
# applied, that a branch was taken -- is a number it can return instead, and the
# scenario function is the one place that decides what is worth a `message()`.
# Every fact this file could have printed is in the vector it hands back.
#
# Two things differ from the numeric kernels. There is one table rather than one
# sample per feature, so a kernel here is called once and an unusable table is an
# error the caller sees rather than an NA row in a long scan. And the engines warn
# about the expected counts, which `$diagnostics` already reports as a number; the
# warning is muffled at the call site rather than passed on, because the answer to
# it is the exact test that is already sitting in the same result.
#
# None of these tests reports an interval for the association itself. A test of
# independence says the two variables are not independent; it does not say by how
# much, and the measures that do say it live in `$association`. The `lower_conf`
# and `upper_conf` columns are therefore present and `NA` except where the test
# itself defines an interval, which only Fisher's exact test on a 2 x 2 table
# does.


#' Pearson's chi-square test of independence
#'
#' @param counts Two-dimensional table of counts.
#' @param correct Whether to apply Yates' continuity correction. [stats::chisq.test()]
#'   applies it to 2 x 2 tables only, whatever it is told, so it changes nothing
#'   on a larger table.
#' @param simulate_p_value,n_resamples Monte Carlo p-value for a sparse table, and
#'   how many tables to draw for it.
#'
#' @return Named numeric vector: `n_used`, `statistic`, `df`, `pval`,
#'   `lower_conf`, `upper_conf`.
#'
#' @references
#' Pearson, K. (1900). On the criterion that a given system of deviations from
#' the probable in the case of a correlated system of variables is such that it
#' can be reasonably supposed to have arisen from random sampling.
#' *Philosophical Magazine*, 50(302), 157-175.
#'
#' @keywords internal
#' @noRd
sa_chisq <- function(counts, correct = TRUE, simulate_p_value = FALSE,
                     n_resamples = 9999) {
  if (any(rowSums(counts) == 0) || any(colSums(counts) == 0)) {
    stop("a row or column holding no observation leaves an expected count of ",
         "zero, which the chi-square statistic divides by. Drop the level from ",
         "`category_lv`.", call. = FALSE)
  }

  # The only warning this engine raises is the one about the approximation, and
  # `$diagnostics` states the same fact as the smallest expected count.
  fit <- suppressWarnings(
    stats::chisq.test(counts, correct = correct,
                      simulate.p.value = simulate_p_value, B = n_resamples)
  )

  c(n_used     = sum(counts),
    statistic  = as.numeric(fit$statistic),
    # A simulated p-value is not referred to a chi-square distribution, so there
    # are no degrees of freedom to report rather than zero of them.
    df         = if (simulate_p_value) NA_real_ else as.numeric(fit$parameter),
    pval       = as.numeric(fit$p.value),
    lower_conf = NA_real_,
    upper_conf = NA_real_)
}


#' Fisher's exact test on a contingency table
#'
#' Conditions on both margins and reads the p-value off the hypergeometric
#' distribution, so it needs no expected count to be large. There is no statistic
#' referred to a null distribution, which is why `statistic` and `df` are `NA`
#' rather than zero.
#'
#' On a 2 x 2 table the test defines an odds ratio of its own, the conditional
#' maximum likelihood estimate, and an interval for it. That is reported here as
#' `odds_ratio_cond` and is **not** the sample odds ratio `$association` carries:
#' the conditional estimate is pulled towards 1 and the two differ most where the
#' counts are smallest, which is the situation this test exists for.
#'
#' The enumeration has a size limit, and reaching it is not an error in the data.
#' The network algorithm walks every table with the observed margins, and on a
#' large r x c table there are more of them than the workspace holds. So the
#' p-value comes back `NA` with `enumerated = 0` rather than as a condition,
#' because the chi-square test standing beside it in the same result was computed
#' and losing it would be the more expensive failure. The caller is what says so
#' and points at `simulate_p_value = TRUE`.
#'
#' @param counts Two-dimensional table of counts.
#' @param conf_level Confidence level for the conditional odds ratio interval.
#' @param simulate_p_value,n_resamples Monte Carlo variant for a large r x c
#'   table, where the exact enumeration is infeasible.
#'
#' @return Named numeric vector: `n_used`, `statistic`, `df`, `pval`,
#'   `lower_conf`, `upper_conf`, `odds_ratio_cond`, `enumerated`. `enumerated` is
#'   whether the test could be computed at all, which the caller reports and does
#'   not carry into the test table.
#'
#' @references
#' Fisher, R. A. (1935). The logic of inductive inference. *Journal of the Royal
#' Statistical Society*, 98(1), 39-82.
#'
#' @keywords internal
#' @noRd
sa_fisher <- function(counts, conf_level = 0.95, simulate_p_value = FALSE,
                      n_resamples = 9999) {
  square <- identical(dim(counts), c(2L, 2L))

  fit <- tryCatch(
    if (square) {
      stats::fisher.test(counts, conf.level = conf_level)
    } else {
      stats::fisher.test(counts, simulate.p.value = simulate_p_value,
                         B = n_resamples)
    },
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(c(n_used = sum(counts), statistic = NA_real_, df = NA_real_,
             pval = NA_real_, lower_conf = NA_real_, upper_conf = NA_real_,
             odds_ratio_cond = NA_real_, enumerated = 0))
  }

  c(n_used          = sum(counts),
    statistic       = NA_real_,
    df              = NA_real_,
    pval            = as.numeric(fit$p.value),
    lower_conf      = if (square) as.numeric(fit$conf.int[1]) else NA_real_,
    upper_conf      = if (square) as.numeric(fit$conf.int[2]) else NA_real_,
    odds_ratio_cond = if (square) as.numeric(fit$estimate) else NA_real_,
    enumerated      = 1)
}


#' McNemar's test of symmetry
#'
#' Reads only the two discordant cells. A pair that answered the same way under
#' both conditions carries no information about which condition raises the
#' response, so the concordant cells do not enter the statistic and the number of
#' pairs does not either, beyond fixing how many discordant ones there could
#' have been.
#'
#' The uncorrected statistic is `(b - c)^2 / (b + c)`, which is exactly the sum of
#' squared Pearson residuals of the table against the symmetry expectation
#' `sa_expected_symmetry()` builds. The cell table of the result and the p-value
#' here are the same arithmetic read two ways.
#'
#' @param counts A 2 x 2 table crossing the two conditions.
#' @param correct Whether to apply the continuity correction to the chi-square
#'   approximation. Ignored by the exact branch, which needs none.
#' @param exact `TRUE` for the exact binomial test on the discordant pairs,
#'   `FALSE` for the chi-square approximation, or `NULL` to take the exact test
#'   when there are fewer than 25 discordant pairs, which is the rule
#'   `assumptions.yaml` records as `discordant_pair_count`.
#'
#' @return Named numeric vector: `n_used`, `n_discordant`, `exact_used`,
#'   `statistic`, `df`, `pval`, `lower_conf`, `upper_conf`. `exact_used` is how
#'   the branch taken under `exact = NULL` reaches the caller, which records it as
#'   `parameters$exact`; it is a setting rather than a finding, so it does not
#'   travel on into the test table.
#'
#' @references
#' McNemar, Q. (1947). Note on the sampling error of the difference between
#' correlated proportions or percentages. *Psychometrika*, 12(2), 153-157.
#'
#' @keywords internal
#' @noRd
sa_mcnemar <- function(counts, correct = TRUE, exact = NULL) {
  b <- counts[1, 2]
  c_ <- counts[2, 1]
  n_discordant <- b + c_

  if (n_discordant == 0L) {
    stop("every pair answered the same way under both conditions, so there is ",
         "no discordance for McNemar's test to be about.", call. = FALSE)
  }

  use_exact <- if (is.null(exact)) n_discordant < 25L else exact
  if (use_exact) {
    pval <- stats::binom.test(b, n_discordant, p = 0.5)$p.value
    statistic <- NA_real_
    df <- NA_real_
  } else {
    adjust <- if (correct) 1 else 0
    statistic <- max(abs(b - c_) - adjust, 0)^2 / n_discordant
    df <- 1
    pval <- stats::pchisq(statistic, df = 1, lower.tail = FALSE)
  }

  c(n_used       = sum(counts),
    n_discordant = n_discordant,
    exact_used   = as.numeric(use_exact),
    statistic    = statistic,
    df           = df,
    pval         = pval,
    lower_conf   = NA_real_,
    upper_conf   = NA_real_)
}


#' Cochran's Q test for three or more repeated binary conditions
#'
#' The extension of McNemar's test past two conditions. Written out rather than
#' taken from an engine because the statistic is also what Kendall's W is built
#' from, and computing it twice from two code paths is how the two end up
#' disagreeing.
#'
#' A subject who answered the same way under every condition contributes nothing
#' to the numerator, which is the same fact that makes the concordant cells drop
#' out of McNemar's test. Such subjects are kept in the denominator, where the
#' formula puts them.
#'
#' @param responses Subjects-by-conditions matrix of 0 and 1, no missing values.
#'
#' @return Named numeric vector: `n_used`, `n_conditions`, `statistic`, `df`,
#'   `pval`, `lower_conf`, `upper_conf`.
#'
#' @references
#' Cochran, W. G. (1950). The comparison of percentages in matched samples.
#' *Biometrika*, 37(3-4), 256-266.
#'
#' @keywords internal
#' @noRd
sa_cochran_q <- function(responses) {
  k <- ncol(responses)
  n <- nrow(responses)
  col_n <- colSums(responses)
  row_n <- rowSums(responses)
  total <- sum(responses)

  denominator <- k * total - sum(row_n^2)
  if (denominator <= 0) {
    stop("every subject answered the same way under every condition, leaving ",
         "Cochran's Q undefined: there is no within-subject variation for it ",
         "to compare across conditions.", call. = FALSE)
  }

  q_stat <- k * (k - 1) * sum((col_n - total / k)^2) / denominator

  c(n_used       = n,
    n_conditions = k,
    statistic    = q_stat,
    df           = k - 1,
    pval         = stats::pchisq(q_stat, df = k - 1, lower.tail = FALSE),
    lower_conf   = NA_real_,
    upper_conf   = NA_real_)
}


#' Association measures for an independent r x c table
#'
#' An effect-size builder rather than a test kernel, so it returns a data.frame in
#' the shape `$association` takes, the way `sa_multi_fold_change()` returns the
#' shape `$effect` takes.
#'
#' Every measure is built from the **uncorrected** chi-square statistic, whichever
#' way `correct` was set for the test. Yates' correction is about the tail
#' probability of a discrete statistic referred to a continuous distribution; it
#' is not about how far the table sits from independence, and letting it into the
#' effect size would make the reported strength of an association depend on a
#' choice made about its p-value.
#'
#' @param counts Two-dimensional table of counts.
#' @param conf_level Confidence level for the odds ratio interval.
#'
#' @return data.frame with `measure`, `estimate`, `lower_conf`, `upper_conf`.
#'
#' @references
#' Cramer, H. (1946). *Mathematical Methods of Statistics*. Princeton University
#' Press.
#'
#' Agresti, A. (2002). *Categorical Data Analysis*, 2nd ed. Wiley.
#'
#' @keywords internal
#' @noRd
sa_assoc_measures <- function(counts, conf_level = 0.95) {
  n <- sum(counts)
  expected <- sa_expected_independence(counts)
  chi_sq <- sum((counts - expected)^2 / expected)
  min_df <- min(nrow(counts), ncol(counts)) - 1L

  rows <- list(
    sa_assoc_row("cramers_v", sqrt(chi_sq / (n * min_df))),
    sa_assoc_row("contingency_coefficient", sqrt(chi_sq / (chi_sq + n)))
  )

  if (identical(dim(counts), c(2L, 2L))) {
    rows <- c(rows,
              list(sa_assoc_row("phi_coefficient", sa_phi(counts)),
                   sa_odds_ratio(counts, conf_level)))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Association measures for a matched 2 x 2 table
#'
#' All three read the discordant cells, which is where the whole of a matched
#' comparison lives. `b` is the pairs that answered `category_lv[2]` under the
#' second condition and `category_lv[1]` under the first, and `c` the reverse.
#'
#' When every discordant pair moved the same way the odds ratio is unbounded, so
#' its estimate and one of its limits are `NA`. The other two measures are bounded
#' by construction and still carry the finding: `cohens_g` reaches its extreme of
#' 0.5 and `risk_difference_paired` reports how large a share of the pairs moved.
#'
#' @inheritParams sa_assoc_measures
#'
#' @keywords internal
#' @noRd
sa_assoc_measures_paired <- function(counts, conf_level = 0.95) {
  n <- sum(counts)
  b <- counts[1, 2]
  c_ <- counts[2, 1]
  n_discordant <- b + c_

  if (n_discordant == 0L) {
    share <- c(NA_real_, NA_real_, NA_real_)
  } else {
    # Clopper-Pearson on the split of the discordant pairs, which is exactly the
    # quantity the exact branch of McNemar's test is about, so the interval and
    # the p-value cannot disagree about whether 0.5 is plausible.
    ci <- stats::binom.test(b, n_discordant, p = 0.5,
                            conf.level = conf_level)$conf.int
    share <- c(b / n_discordant, ci[1], ci[2])
  }

  # The paired odds ratio is a monotone transform of that same share, so its
  # interval is the transformed one rather than a second calculation.
  paired_or <- share / (1 - share)

  diff <- (b - c_) / n
  se <- sqrt(b + c_ - (b - c_)^2 / n) / n
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  out <- rbind(
    sa_assoc_row("odds_ratio_paired", paired_or[1], paired_or[2], paired_or[3]),
    sa_assoc_row("risk_difference_paired", diff,
                 max(diff - z * se, -1), min(diff + z * se, 1)),
    sa_assoc_row("cohens_g", share[1] - 0.5, share[2] - 0.5, share[3] - 0.5)
  )
  rownames(out) <- NULL
  out
}


#' Association measure for three or more repeated binary conditions
#'
#' Kendall's W rescales Cochran's Q by the largest value it could have taken for
#' this many subjects and conditions, which is what turns a statistic that grows
#' with the sample into a measure that does not.
#'
#' @param q_stat Cochran's Q, as returned by `sa_cochran_q()`.
#' @param n_subjects,k Subjects and conditions the statistic was computed on.
#'
#' @keywords internal
#' @noRd
sa_assoc_measures_repeated <- function(q_stat, n_subjects, k) {
  out <- sa_assoc_row("kendalls_w", q_stat / (n_subjects * (k - 1)))
  rownames(out) <- NULL
  out
}


#' One row of an association table
#'
#' @keywords internal
#' @noRd
sa_assoc_row <- function(measure, estimate, lower_conf = NA_real_,
                         upper_conf = NA_real_) {
  data.frame(
    measure    = measure,
    estimate   = sa_finite_or_na(as.numeric(estimate)),
    lower_conf = sa_finite_or_na(as.numeric(lower_conf)),
    upper_conf = sa_finite_or_na(as.numeric(upper_conf)),
    stringsAsFactors = FALSE
  )
}


#' The signed phi coefficient of a 2 x 2 table
#'
#' Its magnitude is `sqrt(chi_sq / n)` on the uncorrected statistic, so phi says
#' nothing Cramer's V does not on a 2 x 2 table, where the two are equal in size.
#' What it adds is the sign, and the sign is the finding: phi above zero means the
#' second level of each variable occurs with the second level of the other more
#' often than independence would give, which is the same direction the odds ratio
#' reads above 1.
#'
#' @keywords internal
#' @noRd
sa_phi <- function(counts) {
  a <- counts[1, 1]
  b <- counts[1, 2]
  c_ <- counts[2, 1]
  d <- counts[2, 2]
  denominator <- sqrt(prod(c(a + b, c_ + d, a + c_, b + d)))
  if (denominator == 0) {
    return(NA_real_)
  }
  (a * d - b * c_) / denominator
}


#' Whether a 2 x 2 table holds a cell the odds ratio cannot be read off
#'
#' The scenario asks this rather than reading it off the measure, so that the one
#' message about the Haldane-Anscombe correction is raised where every other
#' message to the user is raised and `sa_odds_ratio()` stays a function of its
#' arguments.
#'
#' @keywords internal
#' @noRd
sa_has_zero_cell <- function(counts) {
  identical(dim(counts), c(2L, 2L)) && any(counts == 0)
}


#' The sample odds ratio of a 2 x 2 table, with a Wald interval
#'
#' Read against the first level of each variable, which `control_label` fixes. An
#' odds ratio above 1 says the two second levels go together, and pointing
#' `control_label` at the other level of either variable inverts it.
#'
#' The interval is built on the log scale and exponentiated, so it is asymmetric
#' about the estimate on the ratio scale and cannot reach below zero. A zero cell
#' leaves the log undefined, so the Haldane-Anscombe correction of half an
#' observation per cell is applied: the estimate is then a shrunken one rather
#' than an infinite one. `sa_has_zero_cell()` is how the scenario knows to say so.
#'
#' @keywords internal
#' @noRd
sa_odds_ratio <- function(counts, conf_level = 0.95) {
  cells <- c(counts[1, 1], counts[1, 2], counts[2, 1], counts[2, 2])
  if (any(cells == 0)) {
    cells <- cells + 0.5
  }

  log_or <- log(cells[1] * cells[4] / (cells[2] * cells[3]))
  se <- sqrt(sum(1 / cells))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  sa_assoc_row("odds_ratio", exp(log_or), exp(log_or - z * se),
               exp(log_or + z * se))
}
