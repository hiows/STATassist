# The classification counterpart of `simulate_regression()`. The two are
# deliberately near-identical below the roxygen: the same predictors, the same
# planted coefficients, the same subjects. What is specific to this one is that the
# outcome is a draw rather than a sum, and the two consequences of that.
#
# The first is the intercept. In a regression it is a number the caller picks and
# nothing depends on it; here it is what decides how many events there are, and
# how many events there are is what decides whether a split has to be stratified.
# So it is not asked for. `event_rate` is asked for and the intercept is solved
# from it.
#
# The second is that a class is drawn once per subject rather than once per row.
# A subject is a case or a control as a whole, which is the shape of the real
# design this stands in for and also the reason a row-wise split of it is
# worthless: the same subject's other rows carry its label.

#' Simulate a two-class outcome whose coefficients are known
#'
#' Generates a two-class outcome from a logistic model of its predictors and
#' returns the coefficients that were planted alongside the data, so that a fitted
#' model can be scored against what was actually there. The same design arguments
#' [simulate_regression()] takes are available, together with the class balance,
#' which is the reason a split of a classification has to be stratified.
#'
#' A classification differs from a regression in what can go wrong with it, and the
#' defaults are set so that all of it can be seen. Classes are imbalanced, so an
#' unstratified split can hand a fold too few events to fit on. A subject is a case
#' or a control as a whole, so a split that does not respect `id` scores the model
#' on rows whose label it already holds. And a predictor with a coefficient of
#' exactly zero is null in the strict sense, so an odds ratio away from 1 on one is
#' a mistake by definition rather than by judgement.
#'
#' @details
#' The linear predictor is
#' `intercept + sum(beta * x) + factor offsets + subject offset`, the class
#' probability is its logistic transform, and the class is a Bernoulli draw from
#' that probability. There is no noise argument: the draw is the noise, which is
#' why a logistic regression recovers less from the same number of rows than a
#' linear one does.
#'
#' The intercept is not an argument. It is solved for so that the mean class
#' probability over the rows that were actually drawn equals `event_rate`, which
#' means the balance of the data is what was asked for rather than whatever the
#' coefficients happened to imply. `truth_model$intercept` reports the value it
#' took and `truth_model$achieved_event_rate` the proportion the Bernoulli draw
#' then produced.
#'
#' `outcome_lv` fixes the direction by the rule the rest of the package follows:
#' the first level is the reference, so a planted positive coefficient raises the
#' chance of `outcome_lv[2]` and its odds ratio comes back above 1. It is carried
#' in `args` rather than left out, because [fit_logistic_regression()] sorts the
#' classes when it is not told them, and `sort(c("case", "control"))` puts `case`
#' first, which would report the odds of the wrong class.
#'
#' @inheritParams simulate_regression
#' @param n_pred Number of numeric predictors. Columns are named `x_1` upwards, or
#'   whatever `pred_prefix` asks for.
#' @param event_rate Proportion of rows in `outcome_lv[2]`, the class being
#'   modelled. The default is deliberately away from a half, since a balanced
#'   outcome makes `split_data()`'s stratification look unnecessary.
#' @param outcome_lv The two class labels, the reference first, so that the
#'   coefficients describe the odds of the second one.
#' @param n_per_subject Rows measured on each subject, one entry per subject, so
#'   that its length is how many subjects there are. A single number is spread over
#'   `n_samples` rows and must divide them. `NULL`, the default, gives no `subject`
#'   column. With subjects the class is drawn once per subject from the mean of its
#'   rows' probabilities, so a subject is a case or a control as a whole and the
#'   outcome can still stratify a split taken over subjects.
#'
#' @return A list of six elements, the same shape [simulate_regression()] returns,
#'   with these differences:
#'
#'   \describe{
#'     \item{`args`}{Also carries `outcome_lv`, since the direction of every
#'       coefficient depends on it and the default would sort the labels the other
#'       way round.}
#'     \item{`split_args`}{`stratified` is always the outcome. Unlike a continuous
#'       one it is constant within a subject, so it stratifies a split over
#'       subjects as readily as one over rows.}
#'     \item{`truth_model`}{`intercept` as solved, the `event_rate` asked for and
#'       the `achieved_event_rate` the draw produced, `signal_var`, `subject_var`,
#'       `n_samples`, `n_subject` and `subject_sd`. There is no `r_squared`: the
#'       outcome is a draw, so no share of its variance is recoverable in that
#'       sense.}
#'     \item{`truth_row`}{`prob`, the class probability of the row, and
#'       `draw_prob`, the probability the Bernoulli draw actually used, which is the
#'       subject's mean when there are subjects and `prob` itself when there are
#'       not.}
#'   }
#'
#' @section What the defaults are tuned for:
#' The same design [simulate_regression()] uses — eight numeric predictors, four of
#' them planted, one categorical predictor and 200 rows — at an event rate of 0.3,
#' so that the two simulators differ in the outcome and in nothing else. Averaged
#' over twenty seeds that recovers 0.83 of the planted coefficients at
#' `p <= 0.05`, the difficulty the rest of the simulators in this package are tuned
#' to. Coefficients of 0.5 to 2 on the log odds scale are odds ratios of roughly
#' 1.6 to 7, which is what it takes for a class label to carry as much as a number
#' does: the Bernoulli draw is this model's noise and there is no argument to turn
#' it down.
#'
#' @seealso [fit_logistic_regression()], which consumes `args` directly,
#'   [simulate_regression()] for a continuous outcome, [split_data()], which
#'   consumes `split_args`, and [make_block_cor()] for `cor_mat`.
#'
#' @examples
#' sim <- simulate_classification(seed = 1)
#' table(sim$args$data$y)
#'
#' ## The names in `args` are fit_logistic_regression()'s own, `outcome_lv`
#' ## included, so the fit is one call away and points the way it was planted.
#' fit <- do.call(fit_logistic_regression, c(sim$args, cv = FALSE))
#'
#' ## A planted positive coefficient raises the chance of the second level, so
#' ## its odds ratio is above 1. `truth_term` is in the row order the table uses.
#' scored <- merge(fit$coefficients, sim$truth_term, by = "terms")
#' scored[, c("terms", "beta", "estimate", "odds_ratio", "pval")]
#' table(planted = scored$beta != 0, called = scored$pval <= 0.05)
#'
#' ## The intercept was not asked for, it was solved for: the balance of the data
#' ## is the balance that was requested.
#' c(asked = sim$truth_model$event_rate,
#'   drawn = sim$truth_model$achieved_event_rate)
#'
#' ## Two samples per subject, and a subject is a case or a control as a whole.
#' ## The outcome can therefore still be the stratifier of a split over subjects.
#' rep_sim <- simulate_classification(n_per_subject = rep(2, 100), seed = 2)
#' sp <- do.call(split_data, c(rep_sim$split_args, seed = 1))
#' sp
#'
#' @export
simulate_classification <- function(n_samples = 200,
                                    n_pred = 8,
                                    n_pos = round(0.25 * n_pred),
                                    n_neg = round(0.25 * n_pred),
                                    beta = NULL,
                                    beta_range = c(0.5, 2),
                                    event_rate = 0.3,
                                    outcome_lv = c("control", "case"),
                                    value_mean = 0,
                                    value_sd = 1,
                                    cor_mat = NULL,
                                    n_factor_pred = 1,
                                    factor_lv = c("low", "mid", "high"),
                                    n_constant_pred = 0,
                                    p_missing = 0,
                                    n_per_subject = NULL,
                                    subject_sd = 1,
                                    subject_share = 0.5,
                                    pred_prefix = "x",
                                    seed = NULL) {

  sa_check_scalar_num(event_rate, "event_rate", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)
  if (!is.character(outcome_lv) || length(outcome_lv) != 2L ||
        anyNA(outcome_lv) || outcome_lv[1] == outcome_lv[2]) {
    stop("`outcome_lv` must be two distinct non-missing class labels, the ",
         "reference first.", call. = FALSE)
  }
  explicit <- c(if (!missing(n_pred)) "n_pred", if (!missing(n_pos)) "n_pos",
                if (!missing(n_neg)) "n_neg")

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  design <- sa_sim_supervised_design(
    n_samples, n_pred, beta, n_pos, n_neg, beta_range, value_mean, value_sd,
    cor_mat, n_factor_pred, factor_lv, n_constant_pred, p_missing,
    n_per_subject, subject_sd, subject_share, pred_prefix, explicit,
    missing(n_samples)
  )

  intercept <- sa_sim_solve_intercept(design$eta, event_rate)
  eta <- intercept + design$eta
  prob <- stats::plogis(eta)

  if (is.null(design$sizes)) {
    draw_prob <- prob
    event <- stats::rbinom(design$n_samples, 1L, prob)
  } else {
    # One draw per subject, from the mean of the probabilities of its rows. A
    # draw per row would let a subject be a case in one sample and a control in
    # the next, which is not the design `id` exists to protect.
    key <- factor(rep(seq_along(design$sizes), times = design$sizes),
                  levels = seq_along(design$sizes))
    unit_prob <- vapply(split(prob, key), mean, numeric(1))
    draw_prob <- rep(unit_prob, times = design$sizes)
    event <- rep(stats::rbinom(length(unit_prob), 1L, unit_prob),
                 times = design$sizes)
  }

  data <- data.frame(y = outcome_lv[1L + event], stringsAsFactors = FALSE)
  data[design$predictors] <- design$x
  if (!is.null(design$subject)) {
    data$subject <- design$subject
  }

  list(
    args = list(
      data       = data,
      outcome    = "y",
      predictors = design$predictors,
      outcome_lv = outcome_lv
    ),
    split_args = sa_sim_split_args(data, design, stratify_outcome = TRUE),
    truth      = design$truth,
    truth_term = sa_sim_add_intercept(design$truth_term, intercept),
    truth_model = list(
      intercept           = intercept,
      event_rate          = event_rate,
      achieved_event_rate = mean(event),
      signal_var          = stats::var(design$eta - design$subject_offset),
      subject_var         = if (is.null(design$sizes)) {
        0
      } else {
        stats::var(design$subject_offset)
      },
      n_samples           = design$n_samples,
      n_subject           = if (is.null(design$sizes)) {
        NA_integer_
      } else {
        length(design$sizes)
      },
      subject_sd          = if (is.null(design$sizes)) NA_real_ else subject_sd
    ),
    truth_row = data.frame(
      subject        = if (is.null(design$subject)) {
        rep(NA_character_, design$n_samples)
      } else {
        design$subject
      },
      subject_offset = design$subject_offset,
      eta            = eta,
      prob           = prob,
      draw_prob      = unname(draw_prob),
      stringsAsFactors = FALSE
    )
  )
}
