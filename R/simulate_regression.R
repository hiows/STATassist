# The counterpart of `simulate_two_groups()` for the supervised learning family.
# The expression simulators plant an effect in some features and leave the rest
# strictly null, so that both kinds of mistake a comparison can make are defined.
# A regression is scored the same way, on the axis it has: a coefficient of
# exactly zero is a predictor that a p-value below the cutoff is wrong about.
#
# What is new here is that the two mistakes are not independent. A null predictor
# correlated with a planted one carries real information about the outcome, so its
# estimate is pulled away from the zero it truly has, and no amount of data fixes
# it. `cor_mat` is how that is put into the data on purpose rather than met by
# accident.

#' Simulate a regression whose coefficients are known
#'
#' Generates a continuous outcome from a linear combination of predictors and
#' returns the coefficients that were planted alongside the data, so that a fitted
#' model can be scored against what was actually there. Everything a model has to
#' survive can be asked for: predictors that correlate, a categorical predictor, a
#' predictor that takes one value, missing cells, and repeated measurements of the
#' same subject.
#'
#' The point of the exercise is the gap between the coefficient table and the
#' truth, and the three things that open it. A planted coefficient can be too
#' small for the noise to let it through, a null predictor correlated with a
#' planted one is estimated away from zero however much data there is, and a
#' subject measured repeatedly makes a row-wise split score the model on rows it
#' half knows already. The defaults are set so that a run recovers most but not
#' all of what was planted, because a simulation that recovers everything teaches
#' none of this.
#'
#' @details
#' Each row is drawn as
#' `y = intercept + sum(beta * x) + factor offsets + subject offset + noise`,
#' with the numeric predictors drawn from a multivariate normal whose correlations
#' are `cor_mat` and the noise normal with standard deviation `noise_sd`. Because
#' the outcome is built from the coefficients and nothing else, a predictor whose
#' coefficient is zero is null in the strict sense, and a p-value below the cutoff
#' on one is a false positive by definition.
#'
#' Which predictors carry a planted coefficient is drawn at random, but how many
#' carry a positive one and how many a negative one is not: `n_pos` and `n_neg` are
#' counts, so they do not move with the seed. `beta` states every coefficient
#' instead, in which case nothing is planted and its length is how many numeric
#' predictors there are.
#'
#' @param n_samples Rows to generate. Ignored when `n_per_subject` gives a row
#'   count per subject, since those already say how many rows there are.
#' @param n_pred Number of numeric predictors. Columns are named `x_1` upwards, or
#'   whatever `pred_prefix` asks for.
#' @param n_pos,n_neg How many numeric predictors are given a positive and a
#'   negative coefficient. Their sum cannot exceed `n_pred`, and every other
#'   numeric predictor is left with a coefficient of exactly zero. The defaults are
#'   a fraction of `n_pred` rather than a fixed count, so that asking for fewer
#'   predictors plants fewer coefficients instead of failing.
#' @param beta The coefficients themselves, one per numeric predictor and no
#'   intercept among them, or `NULL` to plant `n_pos` and `n_neg` of them. Its
#'   length is then how many numeric predictors there are, so `n_pred` need not be
#'   given as well; naming both and disagreeing is an error rather than a guess.
#'   Supplying it together with `n_pos` or `n_neg` is refused, since the two are
#'   different ways of saying the same thing.
#' @param beta_range Range the magnitude of a planted coefficient is drawn from.
#'   The offsets of the factor predictors are drawn from it too, so it is read
#'   whether or not `beta` was given.
#' @param intercept The intercept. It is not part of `beta`.
#' @param value_mean,value_sd Mean and standard deviation of each numeric
#'   predictor, given once for all of them or once each. A coefficient is a change
#'   in the outcome per unit of its predictor, so these fix what `beta_range`
#'   means; that is why they are given rather than drawn from a range as the
#'   expression simulators draw their spreads.
#' @param noise_sd Standard deviation of the residual noise. This and `beta_range`
#'   together decide how much of the outcome is recoverable at all;
#'   `truth_model$r_squared` reports how much, on the design that was drawn.
#' @param cor_mat Correlation matrix of the numeric predictors, as built by
#'   [make_block_cor()], or `NULL` to leave them independent. A null predictor
#'   correlated with a planted one is the case a coefficient table gets wrong no
#'   matter how many rows it is given.
#' @param n_factor_pred Number of categorical predictors, named `x_cat_1` upwards.
#'   Each becomes `length(factor_lv) - 1` terms in the model rather than one, which
#'   is why `truth_term` exists beside `truth`. Levels are handed out in balanced
#'   counts, and each level beyond the first carries a planted offset.
#' @param factor_lv Levels of each categorical predictor, the first being the
#'   reference that carries no offset.
#' @param n_constant_pred Number of predictors that take a single value, named
#'   `x_const_1` upwards. They cannot contribute, so [fit_linear_regression()]
#'   leaves them out and names them in `design$dropped_predictors`. The default
#'   plants none.
#' @param p_missing Proportion of numeric predictor cells to replace with `NA`,
#'   drawn after the outcome has been computed from the complete values. The rows
#'   they fall in are the ones a model drops before its folds are laid out, and
#'   `design$n_dropped` counts them. Only the numeric predictors are holed, since a
#'   hole in the categorical one would stop it being able to stratify a split and a
#'   hole in a constant one would stop it being constant.
#' @param n_per_subject Rows measured on each subject, one entry per subject, so
#'   that its length is how many subjects there are. A single number is spread over
#'   `n_samples` rows and must divide them. `NULL`, the default, gives no `subject`
#'   column and one row per sampling unit. This is what [split_data()]'s `id`
#'   argument exists for: a subject partly seen in training is partly known before
#'   its test rows are read.
#' @param subject_sd Standard deviation of the per-subject offset on the outcome,
#'   which is variation no predictor accounts for. Ignored without
#'   `n_per_subject`.
#' @param subject_share Share of each numeric predictor's variance that lies
#'   between subjects rather than within one, so its intraclass correlation. This
#'   is what makes two rows of one subject resemble each other, and it is that
#'   resemblance a row-wise split gives away: at `0` the rows of a subject are
#'   independent draws that happen to share an outcome offset, and no model could
#'   tell one subject's rows from another's. How much a model gains from the
#'   resemblance depends on the model, and a regression with a fixed set of
#'   coefficients gains little, so the argument is what makes the data set the
#'   shape `id` is for rather than a way to inflate a score. The distribution of
#'   each column is the same whatever it is set to. Ignored without
#'   `n_per_subject`.
#' @param pred_prefix Prefix for the generated predictor names. `"x"` gives `x_1`,
#'   `x_cat_1` and `x_const_1`.
#' @param seed Seed for the draw, or `NULL` to use the stream as it stands.
#'   Supplying one does not disturb the caller: the previous random number state is
#'   put back when the function returns.
#'
#' @return A list of six elements.
#'
#'   \describe{
#'     \item{`args`}{`data`, `outcome` and `predictors`, named after the arguments
#'       of [fit_linear_regression()] so that
#'       `do.call(fit_linear_regression, sim$args)` fits the model. `predictors`
#'       is given explicitly rather than left to its `NULL` default, which would
#'       take the `subject` column as a predictor and let the model fit on which
#'       subject a row came from.}
#'     \item{`split_args`}{`data`, `stratified` and `id`, named after the arguments
#'       of [split_data()]. The outcome is the stratifier when there are no
#'       subjects; with subjects it varies within a subject and so cannot stratify
#'       a split taken over them, and the first categorical predictor, which is
#'       drawn per subject, is used instead.}
#'     \item{`truth`}{One row per predictor, in the column order of `data`, holding
#'       `predictors`, `role` (`"signal"`, `"null"`, `"factor"` or `"constant"`),
#'       `beta`, `direction`, `value_mean`, `value_sd` and `max_cor_signal`, the
#'       largest correlation this predictor has with a planted one. The last is
#'       what accounts for a null predictor that came back significant.}
#'     \item{`truth_term`}{One row per model term, in the row order
#'       `coefficients` follows, holding `terms`, the `predictors` each term came
#'       from, and `beta`. This is the table that scores the coefficients, since a
#'       categorical predictor is several terms and a constant one is none.}
#'     \item{`truth_model`}{The model as a whole: `intercept`, `noise_sd`,
#'       `signal_var`, `subject_var` and `r_squared`, the share of the variance of
#'       the outcome the predictors account for, which is what `fit_stats$r_squared`
#'       estimates. Also `n_samples`, `n_subject` and `subject_sd`.}
#'     \item{`truth_row`}{One row per observation, holding `subject`,
#'       `subject_offset`, `eta` (the whole linear predictor, intercept included)
#'       and `noise`, so that `y` is exactly `eta + noise`.}
#'   }
#'
#' @section What the defaults are tuned for:
#' Eight numeric predictors, four of them planted, one categorical predictor and
#' 200 rows, with `noise_sd = 3` against coefficients of 0.5 to 2 on unit-variance
#' predictors. Averaged over twenty seeds that leaves 47% of the variance of the
#' outcome accounted for by the predictors and recovers 0.83 of the planted
#' coefficients at `p <= 0.05`, which is the same difficulty
#' [simulate_two_groups()] is tuned to. Lowering `noise_sd` to 2 recovers 0.94 and
#' raising it to 5 costs a third of them.
#'
#' The rate at which a null predictor is called is 0.05 and stays there, because a
#' coefficient table applies no multiplicity adjustment across its terms. That is
#' a property of the model rather than of these defaults, and it is one of the
#' things a simulation with strictly null predictors is for: with eight predictors
#' and a cutoff of 0.05, a table that names a predictor it should not is the
#' expected outcome of roughly every other run.
#'
#' @seealso [fit_linear_regression()], which consumes `args` directly,
#'   [simulate_classification()] for a two-class outcome, [split_data()], which
#'   consumes `split_args`, and [make_block_cor()] for `cor_mat`.
#'
#' @examples
#' sim <- simulate_regression(seed = 1)
#' sim$truth[, c("predictors", "role", "beta")]
#'
#' ## The names in `args` are fit_linear_regression()'s own, so the fit is one
#' ## call away.
#' fit <- do.call(fit_linear_regression, c(sim$args, cv = FALSE))
#'
#' ## Scored against what was planted. `truth_term` is already in the row order
#' ## the coefficient table follows, so the two line up without matching.
#' scored <- merge(fit$coefficients, sim$truth_term, by = "terms")
#' scored[, c("terms", "beta", "estimate", "pval")]
#'
#' ## Both kinds of mistake are defined, because a null predictor's coefficient
#' ## is exactly zero rather than merely small.
#' table(planted = scored$beta != 0, called = scored$pval <= 0.05)
#'
#' ## The model as a whole is scored too: r_squared estimates the share of the
#' ## variance the predictors actually carry.
#' c(planted = sim$truth_model$r_squared, fitted = fit$fit_stats$r_squared)
#'
#' ## Correlated predictors are the reason a coefficient table names the wrong
#' ## one. `x_2` is null and correlates 0.9 with `x_1`, which is not.
#' pair <- simulate_regression(
#'   n_pred = 4, beta = c(2, 0, 0, 0),
#'   cor_mat = make_block_cor(4, list(list(features = 1:2, cor = 0.9))),
#'   seed = 2
#' )
#' pair$truth[, c("predictors", "role", "beta", "max_cor_signal")]
#'
#' ## Three measurements per subject, so a split has to be taken over subjects.
#' rep_sim <- simulate_regression(n_per_subject = rep(3, 40), seed = 3)
#' sp <- do.call(split_data, c(rep_sim$split_args, seed = 1))
#' intersect(sp$datasets[[1]]$train_data$subject,
#'           sp$datasets[[1]]$test_data$subject)
#'
#' @export
simulate_regression <- function(n_samples = 200,
                                n_pred = 8,
                                n_pos = round(0.25 * n_pred),
                                n_neg = round(0.25 * n_pred),
                                beta = NULL,
                                beta_range = c(0.5, 2),
                                intercept = 0,
                                value_mean = 0,
                                value_sd = 1,
                                noise_sd = 3,
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

  sa_check_scalar_num(intercept, "intercept")
  sa_check_scalar_num(noise_sd, "noise_sd", 0)
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

  eta <- intercept + design$eta
  noise <- stats::rnorm(design$n_samples, 0, noise_sd)

  data <- data.frame(y = eta + noise)
  data[design$predictors] <- design$x
  if (!is.null(design$subject)) {
    data$subject <- design$subject
  }

  # The subject offset is variance the predictors cannot account for, so it sits
  # with the noise in the denominator rather than with the signal. Counting it as
  # signal would make `r_squared` a number no model could reach.
  signal_var <- stats::var(design$eta - design$subject_offset)
  subject_var <- if (is.null(design$sizes)) {
    0
  } else {
    stats::var(design$subject_offset)
  }

  list(
    args = list(
      data       = data,
      outcome    = "y",
      predictors = design$predictors
    ),
    split_args = sa_sim_split_args(data, design, stratify_outcome = FALSE),
    truth      = design$truth,
    truth_term = sa_sim_add_intercept(design$truth_term, intercept),
    truth_model = list(
      intercept   = intercept,
      noise_sd    = noise_sd,
      signal_var  = signal_var,
      subject_var = subject_var,
      r_squared   = signal_var / (signal_var + subject_var + noise_sd^2),
      n_samples   = design$n_samples,
      n_subject   = if (is.null(design$sizes)) {
        NA_integer_
      } else {
        length(design$sizes)
      },
      subject_sd  = if (is.null(design$sizes)) NA_real_ else subject_sd
    ),
    truth_row = data.frame(
      subject        = if (is.null(design$subject)) {
        rep(NA_character_, design$n_samples)
      } else {
        design$subject
      },
      subject_offset = design$subject_offset,
      eta            = eta,
      noise          = noise,
      stringsAsFactors = FALSE
    )
  )
}
