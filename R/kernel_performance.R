# Performance kernels for a two-class outcome, written out rather than taken
# from a package.
#
# `pROC` would cover the first three of these and none of the last three, which
# is what decided it: IDI and NRI have no implementation in `pROC` at all, and
# DeLong's test has no counterpart in `scikit-learn` or `scipy` for the Python
# port to call. Depending on `pROC` would therefore buy the easy half and leave
# the hard half to be written anyway, in two languages, against two different
# sets of defaults. Written here once, the formula is the spec, and `pROC` earns
# its place in `Suggests` as the oracle the tests check these against.
#
# Every function takes plain numeric vectors and returns a scalar or a named
# numeric vector, so the port is a transcription rather than a redesign.
#
# `response` is always 0/1 with 1 for the event, which is `outcome_lv[2]`: the
# class `predict(model, type = "response")` reports the probability of. The
# callers convert once and these never see a label.


#' Check that a kernel was handed a usable pair of vectors
#'
#' @keywords internal
#' @noRd
sa_check_response <- function(response, predictor) {
  if (length(response) != length(predictor)) {
    stop("internal error: `response` and `predictor` differ in length.",
         call. = FALSE)
  }
  if (!all(response %in% c(0, 1))) {
    stop("internal error: `response` must be 0/1 with 1 for the event.",
         call. = FALSE)
  }
  n_event <- sum(response == 1)
  if (n_event == 0L || n_event == length(response)) {
    stop("the scored rows hold a single class, so there is nothing to ",
         "discriminate between. Both classes have to be present in `answer`.",
         call. = FALSE)
  }
  invisible(NULL)
}


#' Operating points of a ROC curve
#'
#' One row per distinct predicted value, plus the point above all of them, so a
#' curve of `k` distinct predictions has `k + 1` points running from `(0, 1)` to
#' `(1, 0)`. A row is positive when its prediction is greater than or equal to
#' the threshold.
#'
#' Ties are one point rather than several. Rows sharing a predicted value cannot
#' be separated by any threshold, so the curve steps diagonally through them,
#' which is the same thing counting them as half a concordant pair does in
#' `sa_auc()`.
#'
#' @param response 0/1, 1 for the event.
#' @param predictor Predicted probability of the event.
#'
#' @return data.frame of `threshold`, `sensitivity` and `specificity`.
#'
#' @keywords internal
#' @noRd
sa_roc_points <- function(response, predictor) {
  sa_check_response(response, predictor)
  n_event <- sum(response == 1)
  n_other <- length(response) - n_event

  at <- order(predictor, decreasing = TRUE)
  sorted <- predictor[at]
  hit <- cumsum(response[at])
  miss <- cumsum(1 - response[at])
  # The last index of each run of equal predictions is the only one a threshold
  # can stop at, since every row of the run crosses together.
  last <- !duplicated(sorted, fromLast = TRUE)

  data.frame(
    threshold   = c(Inf, sorted[last]),
    sensitivity = c(0, hit[last]) / n_event,
    specificity = 1 - c(0, miss[last]) / n_other,
    stringsAsFactors = FALSE
  )
}


#' Placement values, the per-row terms an AUC is the mean of
#'
#' DeLong's structural components. For an event, the share of non-events it
#' outranks, counting a tie as half; for a non-event, the share of events that
#' outrank it. Both average to the AUC, and their variances are what its
#' standard error is built from, which is what makes the statistic a mean of
#' independent terms rather than a U statistic to be approximated.
#'
#' Computed from ranks rather than from the `n_event * n_other` comparisons.
#' The rank of an event within the pooled sample, less its rank within the
#' events alone, is the number of non-events below it plus half the number tied
#' with it, which is that row's placement value times `n_other`.
#'
#' @return A list of `event` and `other`, the placement values of each class.
#'
#' @keywords internal
#' @noRd
sa_placement_values <- function(response, predictor) {
  is_event <- response == 1
  x <- predictor[is_event]
  y <- predictor[!is_event]
  n_event <- length(x)
  n_other <- length(y)

  pooled <- rank(c(x, y))
  event <- (pooled[seq_len(n_event)] - rank(x)) / n_other
  other <- 1 - (pooled[n_event + seq_len(n_other)] - rank(y)) / n_event

  list(event = event, other = other)
}


#' Area under the ROC curve
#'
#' The Mann-Whitney statistic scaled to a probability: the chance that a
#' randomly drawn event is ranked above a randomly drawn non-event, with a tie
#' counted as half. This is the same quantity `compare_two_groups()` reports as
#' the `relative_effect` of its Brunner-Munzel test, read here as a statement
#' about a classifier rather than about two samples.
#'
#' @keywords internal
#' @noRd
sa_auc <- function(response, predictor) {
  sa_check_response(response, predictor)
  is_event <- response == 1
  n_event <- sum(is_event)
  n_other <- length(response) - n_event
  ranks <- rank(predictor)
  (sum(ranks[is_event]) - n_event * (n_event + 1) / 2) / (n_event * n_other)
}


#' Area under the ROC curve with DeLong's standard error
#'
#' The variance of a mean of placement values, one variance per class, which is
#' the non-parametric standard error of an AUC and the one the interval beside
#' it in `$metrics` is built from.
#'
#' @return Named numeric of `auc` and `se`.
#'
#' @keywords internal
#' @noRd
sa_auc_delong <- function(response, predictor) {
  sa_check_response(response, predictor)
  placement <- sa_placement_values(response, predictor)
  n_event <- length(placement$event)
  n_other <- length(placement$other)

  auc <- mean(placement$event)
  # A single row of a class leaves its variance undefined, and an AUC resting
  # on one observation has no spread to report rather than a spread of zero.
  variance <- if (n_event > 1L && n_other > 1L) {
    stats::var(placement$event) / n_event +
      stats::var(placement$other) / n_other
  } else {
    NA_real_
  }

  c(auc = auc, se = sqrt(variance))
}


#' DeLong's test for two AUCs measured on the same rows
#'
#' Paired, because both models ranked the same rows and their placement values
#' therefore covary. Ignoring that covariance would treat two models that agree
#' almost everywhere as two independent estimates and make every difference
#' between them look more surprising than it is.
#'
#' @param response 0/1, 1 for the event.
#' @param predictor_1,predictor_2 The two sets of predicted probabilities, in
#'   the direction `predictor_1 - predictor_2`.
#'
#' @return Named numeric of `delta`, `se`, `statistic` and `pval`.
#'
#' @keywords internal
#' @noRd
sa_delong_test <- function(response, predictor_1, predictor_2) {
  sa_check_response(response, predictor_1)
  sa_check_response(response, predictor_2)

  first <- sa_placement_values(response, predictor_1)
  second <- sa_placement_values(response, predictor_2)
  n_event <- length(first$event)
  n_other <- length(first$other)

  delta <- mean(first$event) - mean(second$event)

  if (n_event < 2L || n_other < 2L) {
    return(c(delta = delta, se = NA_real_, statistic = NA_real_,
             pval = NA_real_))
  }

  s_event <- stats::cov(cbind(first$event, second$event))
  s_other <- stats::cov(cbind(first$other, second$other))
  s <- s_event / n_event + s_other / n_other
  variance <- s[1, 1] + s[2, 2] - 2 * s[1, 2]

  # Two models that rank every row identically differ by exactly zero with a
  # standard error of exactly zero, and the ratio of the two is not a number.
  # Reporting 1 would claim a test was run against a distribution that has no
  # spread, so the columns say there is nothing here instead.
  if (!is.finite(variance) || variance <= 0) {
    return(c(delta = delta, se = 0, statistic = NA_real_, pval = NA_real_))
  }

  se <- sqrt(variance)
  statistic <- delta / se
  c(delta     = delta,
    se        = se,
    statistic = statistic,
    pval      = 2 * stats::pnorm(-abs(statistic)))
}


#' Integrated discrimination improvement
#'
#' How much further apart the two classes' predicted probabilities moved. The
#' mean predicted probability rises among events and falls among non-events for
#' a model that discriminates better, and the IDI is the sum of those two
#' movements, so it is on the probability scale rather than on the scale of a
#' rank.
#'
#' It answers something a difference of AUCs cannot. An AUC sees only the order
#' of the predictions, so a new model that pushes every event's probability up
#' by a tenth without reordering anything leaves it untouched and moves this.
#'
#' The standard error is the one Pencina's paired form gives: the two class-wise
#' mean changes are independent, so their variances add.
#'
#' @return Named numeric of `idi`, `se`, `statistic` and `pval`.
#'
#' @references
#' Pencina, M. J., D'Agostino, R. B., D'Agostino, R. B. and Vasan, R. S. (2008).
#' Evaluating the added predictive ability of a new marker: from area under the
#' ROC curve to reclassification and beyond. *Statistics in Medicine*, 27(2),
#' 157-172.
#'
#' @keywords internal
#' @noRd
sa_idi <- function(response, predictor_old, predictor_new) {
  sa_check_response(response, predictor_old)
  sa_check_response(response, predictor_new)

  is_event <- response == 1
  moved_event <- (predictor_new - predictor_old)[is_event]
  moved_other <- (predictor_new - predictor_old)[!is_event]
  n_event <- length(moved_event)
  n_other <- length(moved_other)

  idi <- mean(moved_event) - mean(moved_other)

  if (n_event < 2L || n_other < 2L) {
    return(c(idi = idi, se = NA_real_, statistic = NA_real_, pval = NA_real_))
  }

  variance <- stats::var(moved_event) / n_event +
    stats::var(moved_other) / n_other
  if (!is.finite(variance) || variance <= 0) {
    return(c(idi = idi, se = 0, statistic = NA_real_, pval = NA_real_))
  }

  se <- sqrt(variance)
  statistic <- idi / se
  c(idi       = idi,
    se        = se,
    statistic = statistic,
    pval      = 2 * stats::pnorm(-abs(statistic)))
}


#' Continuous net reclassification improvement
#'
#' How often a probability moved the right way, which is the third question and
#' the coarsest. Only the direction of each row's change is counted, so it is
#' unmoved by how large the changes were and answers where the IDI does not: a
#' new model that helps most rows a little and hurts a few a great deal has a
#' positive NRI and can have a negative IDI.
#'
#' Category-free, so no risk strata are named. A stratified NRI depends on cut
#' points that are a clinical convention rather than a property of the data, and
#' there is no default for them this package could pick.
#'
#' The two class-wise components are reported beside the total, since they are
#' what it is made of and they routinely point opposite ways.
#'
#' @details
#' The standard error is the non-null one, the variance of the difference of two
#' proportions with the `(p_up - p_down)^2` term kept. Pencina's published test
#' drops that term, which is correct under the null the test is against, where
#' the two proportions are equal and it vanishes. Keeping it is what makes the
#' interval and the p-value here come from one standard error rather than two,
#' at the cost of a p-value very slightly different from the one
#' `Hmisc::improveProb()` reports.
#'
#' @return Named numeric of `nri`, `nri_event`, `nri_other`, `se`, `statistic`
#'   and `pval`.
#'
#' @references
#' Pencina, M. J., D'Agostino, R. B. and Steyerberg, E. W. (2011). Extensions of
#' net reclassification improvement calculations to measure usefulness of new
#' biomarkers. *Statistics in Medicine*, 30(1), 11-21.
#'
#' @keywords internal
#' @noRd
sa_nri <- function(response, predictor_old, predictor_new) {
  sa_check_response(response, predictor_old)
  sa_check_response(response, predictor_new)

  is_event <- response == 1
  moved <- predictor_new - predictor_old
  moved_event <- moved[is_event]
  moved_other <- moved[!is_event]
  n_event <- length(moved_event)
  n_other <- length(moved_other)

  up_event <- mean(moved_event > 0)
  down_event <- mean(moved_event < 0)
  up_other <- mean(moved_other > 0)
  down_other <- mean(moved_other < 0)

  # An event whose probability rose was reclassified towards the truth; a
  # non-event whose probability rose was reclassified away from it, which is why
  # the second component reads the other way round.
  nri_event <- up_event - down_event
  nri_other <- down_other - up_other
  nri <- nri_event + nri_other

  variance <- (up_event + down_event - nri_event^2) / n_event +
    (up_other + down_other - nri_other^2) / n_other

  if (!is.finite(variance) || variance <= 0) {
    return(c(nri = nri, nri_event = nri_event, nri_other = nri_other,
             se = 0, statistic = NA_real_, pval = NA_real_))
  }

  se <- sqrt(variance)
  statistic <- nri / se
  c(nri       = nri,
    nri_event = nri_event,
    nri_other = nri_other,
    se        = se,
    statistic = statistic,
    pval      = 2 * stats::pnorm(-abs(statistic)))
}


#' Brier score
#'
#' Mean squared distance between the predicted probability and the outcome, so
#' it is to a classification what `rmse` squared is to a regression: a model
#' that ranks perfectly but predicts every event at 0.6 is scored here and not
#' by an AUC.
#'
#' @keywords internal
#' @noRd
sa_brier <- function(response, predictor) {
  sa_check_response(response, predictor)
  mean((predictor - response)^2)
}


#' Accuracy, sensitivity and specificity at one stated threshold
#'
#' A row is called an event when its predicted probability is greater than or
#' equal to `threshold`, the same direction `sa_roc_points()` steps in.
#'
#' @keywords internal
#' @noRd
sa_threshold_scores <- function(response, predictor, threshold) {
  sa_check_response(response, predictor)
  called <- as.numeric(predictor >= threshold)
  is_event <- response == 1
  c(accuracy    = mean(called == response),
    sensitivity = mean(called[is_event] == 1),
    specificity = mean(called[!is_event] == 0))
}
