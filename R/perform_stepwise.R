# The second function here that searches rather than fits, and the first that
# searches without holding anything out. `perform_rfe()` scores every subset size
# on rows that did not choose it, which is what a resampled score is for. This one
# has no resampled score at all: what it compares candidates by is a penalised
# likelihood, the fit's own log likelihood with a charge levied per parameter, so a
# predictor has to earn its coefficient before the criterion will keep it.
#
# That is why there is no `cv` argument and no `seed`. Nothing is resampled, so
# nothing is random and the same rows give the same path every time. It is also
# why the criterion is not a performance claim: it is computed on the rows the
# model was fitted to, and the model was chosen because it scored best on them, so
# the number that chose the model cannot also be an honest estimate of it.
#
# The charge per parameter is the whole difference between the two criteria. AIC
# levies 2 and BIC levies `log(n)`, so past seven observations BIC charges more and
# keeps fewer predictors. Both are computed at every step whichever one is
# searching, so a path chosen by one can be read against the other.
#
# What this does not have to do is fold dummy columns back onto the columns they
# came from, which is most of the work inside `perform_rfe()`. `step()` moves whole
# terms of the formula, and a factor is one term however many dummy columns it
# becomes, so what enters and leaves the model is always a column of the input and
# `$selected` is a set of names a `fit_*()` call takes as it stands.

#' Stepwise feature selection by information criteria
#'
#' Runs a stepwise search: the model is refitted with each candidate term taken out
#' or put in, the move that lowers AIC or BIC the most is taken, and the search
#' stops when no single move lowers it any further. What comes back is the
#' predictors of the model it stopped at, the path it walked to get there, and what
#' each candidate is worth to that model, so a set of four predictors can be read
#' against what the other six would have cost.
#'
#' The input is the wide format the model functions take, **one row per observation
#' with one column as the outcome**, and is normally the training half of a
#' [split_data()] result. Selecting on data a model is later scored on is how a
#' selection flatters itself; see the details.
#'
#' @details
#' # What the criterion charges
#'
#' Both criteria are the fit's own log likelihood with a charge levied per
#' parameter, and they differ only in the size of the charge: `"AIC"` levies 2 and
#' `"BIC"` levies `log(n)`, so past seven observations BIC charges more and keeps
#' fewer predictors. Smaller is better for both, which is why `maximize` is `FALSE`
#' and is reported rather than asked for. The charge as it was used is
#' `parameters$k`.
#'
#' Whichever one searches, `profile` reports both at every step, so a path chosen
#' by AIC can be read against what BIC would have said about the same models. They
#' are [stats::AIC()]'s and [stats::BIC()]'s own numbers, the ones
#' [fit_linear_regression()] reports in `fit_stats`, so a step of the path and a
#' fitted model are the same number. [stats::step()] itself searches on
#' [stats::extractAIC()]'s scale, which for a linear model differs from
#' [stats::AIC()]'s by a constant that is the same for every model on `n` rows;
#' the path is therefore ordered identically and only the printed values differ.
#'
#' A criterion is comparable across models only when they were fitted to the same
#' rows, which is what the listwise deletion below is for. It is not comparable
#' across outcomes, across transformations of an outcome, or between a linear and a
#' logistic fit, so the numbers in `profile` say which of these models to prefer
#' and nothing more.
#'
#' # What the search is not
#'
#' The criterion is computed on the rows the model was fitted to, and the model was
#' kept because it scored best on them. Three things follow, and none of them is a
#' fault of the implementation:
#'
#' \describe{
#'   \item{The criterion is not a validation}{Score the selection on the test half
#'     of [split_data()], which the search never saw. That is the only number here
#'     that is a claim about new rows, and it is why this function has no `cv`
#'     argument: cross-validating a penalised likelihood computed on the training
#'     rows would resample the wrong quantity.}
#'   \item{The p-values of the selected model are no longer honest}{A coefficient
#'     kept because it was significant enough to survive is being tested against a
#'     null it was already screened on, so the p-value [fit_linear_regression()]
#'     reports for it on the same rows reads smaller than it is. Refit on held-out
#'     rows, or read the fit as a description of these rows rather than as a test.}
#'   \item{The path is greedy}{One move is taken at a time, so a pair of predictors
#'     that is worth keeping only together can be dropped one at a time and never
#'     come back. `direction = "both"` reconsiders a dropped term at every later
#'     step, which is what it is for; nothing short of fitting every subset removes
#'     the problem entirely.}
#' }
#'
#' [fit_elastic_net()] answers the same question without a path at all, by
#' shrinking a coefficient to exactly zero, and [perform_rfe()] answers it with a
#' resampled score rather than with a penalty.
#'
#' # Where the search starts and which way it moves
#'
#' `"backward"` starts at every candidate and only drops. `"forward"` starts at the
#' intercept and only adds. `"both"` starts at every candidate and may add a term
#' back after dropping it, so its path can visit the same size twice; that is the
#' one direction whose `profile` is not a ladder.
#'
#' Every direction is bounded by the same two models: the intercept alone below and
#' all of `predictors` above. A search that walks back to the intercept has kept
#' nothing, which is an answer — no candidate pays for itself at this charge — but
#' not one this contract can carry, since `selected` would be empty. It is an error
#' saying so rather than a result with a hole in it.
#'
#' # Which model does the searching
#'
#' `model` names what is fitted at every step, and the outcome has to agree with
#' it: `"linear"` is a continuous outcome and `"logistic"` a two-class one. A
#' disagreement is an error naming the model that would have fitted, rather than a
#' silently different analysis. There is no `"rf"` here, unlike [perform_rfe()]: a
#' criterion is a likelihood with a charge against its parameter count, and a
#' forest has neither.
#'
#' `outcome_lv` is read as it is everywhere else in this package: the first level
#' is the reference, so `outcome_lv = c("control", "case")` searches for the
#' predictors of `case`, and `control_label` names that first level on its own. For
#' the search itself the direction changes nothing, since the likelihood of a model
#' is the same number whichever class is called the reference; what the levels
#' decide is which class a later [fit_logistic_regression()] on `$selected` reads.
#'
#' # What is dropped before anything runs
#'
#' The same listwise deletion the model functions use: rows missing the outcome or
#' any candidate go first, and `design$n_dropped` says how many went. Here it is
#' also what makes the search legible, since a criterion compares models only when
#' they were fitted to the same rows, and leaving the deletion to the engine would
#' have each step fitted to whatever its own terms happened to be complete on. A
#' candidate that takes a single value is left out with a message, since a column
#' with nothing in it cannot earn a parameter.
#'
#' # Portability
#'
#' Everything but `$fit` is a scalar, a character vector, a named list or a
#' data.frame, so dropping that one slot leaves an object that writes out as JSON.
#'
#' @param data A data.frame (or matrix) in wide format, one row per observation.
#'   Typically the training half of a [split_data()] result.
#' @param outcome The outcome, either the name of a column of `data` or a vector
#'   with one entry per row.
#' @param predictors Candidate column names, or `NULL` for every column of `data`
#'   except the outcome. Numeric, logical, factor and character columns are all
#'   accepted; a factor is one candidate however many dummy columns it becomes.
#' @param outcome_lv For a two-class outcome, the two classes with the reference
#'   first. `NULL` sorts them, which puts `"control"` before `"treated"` and `0`
#'   before `1`. Naming it is also what tells this function that a numeric column
#'   of zeroes and ones is two classes rather than two numbers.
#' @param control_label The reference class on its own, for when the other one
#'   needs no saying. Defaults to `outcome_lv[1]`, so a call that names one of the
#'   two names the reference either way; naming both and disagreeing is an error.
#' @param model What is fitted at every step: `"linear"` for a continuous outcome
#'   or `"logistic"` for a two-class one.
#' @param criterion Which criterion the moves are judged by, `"AIC"` or `"BIC"`.
#'   The search is the same one either way; what changes is the charge per
#'   parameter, 2 against `log(n)`.
#' @param direction `"backward"` starts at every candidate and only drops,
#'   `"forward"` starts at the intercept and only adds, and `"both"` starts at
#'   every candidate and may put a dropped term back.
#'
#' @return An object of class `sa_selection`, a plain list.
#'
#'   \describe{
#'     \item{`analysis`}{`"stepwise"`.}
#'     \item{`candidates`}{The predictors that were offered, most important first.
#'       This is the row order `ranking` follows, and what `terms` is to a model.}
#'     \item{`design`}{What the search saw: the `outcome` and its `outcome_type`,
#'       `outcome_lv` with `n_events` and `event_rate` for a classification, the
#'       row counts `n_obs`, `n_used` and `n_dropped`, the `predictors` in the
#'       order they arrived, any `dropped_predictors`, and `predictor_lv` for those
#'       that are factors.}
#'     \item{`parameters`}{The choices as they were used: `model`, `criterion`,
#'       `maximize`, which is `FALSE` because a smaller criterion is a better
#'       model, `k`, the charge per parameter it levied, and `direction`.}
#'     \item{`selected`}{The predictors of the model the search stopped at, most
#'       important first. These are the names to hand to `predictors =` in a
#'       `fit_*()` call.}
#'     \item{`ranking`}{One row per candidate: `candidates`, the `estimate` it was
#'       ranked by, its `rank`, and whether it was `selected`. The estimate is what
#'       leaving that one predictor out of the selected model costs the criterion,
#'       so it is positive for a predictor worth keeping and negative for one worth
#'       leaving out.}
#'     \item{`profile`}{One row per step of the path: `n_vars`, both `AIC` and
#'       `BIC` of the model at that step, the `step` that was taken to reach it,
#'       and `chosen`, which is `TRUE` on the last row.}
#'     \item{`resampling`}{`NULL`. Nothing was resampled.}
#'     \item{`engine`}{What ran the search, including `importance`, the name of
#'       what `ranking$estimate` measures.}
#'     \item{`fit`}{The [stats::step()] result, which is the selected `lm` or `glm`
#'       with the path attached as `$anova`. This is the slot that is not portable;
#'       dropping it leaves an object that writes out as JSON.}
#'     \item{`metadata`}{Package version, R version, platform and timestamp.}
#'   }
#'
#' @seealso [split_data()], which defines the rows a selection should be run on,
#'   [perform_rfe()], which chooses by a resampled score rather than by a penalty,
#'   [fit_elastic_net()], which chooses by shrinking a coefficient to exactly zero,
#'   and [fit_linear_regression()] and [fit_logistic_regression()] for fitting the
#'   predictors it kept.
#'
#' @examples
#' ## Six candidates and one continuous outcome. The path says what was dropped and
#' ## when, and `$selected` is where the search stopped.
#' res <- perform_stepwise(mtcars, outcome = "mpg",
#'                     predictors = c("wt", "hp", "disp", "qsec", "drat", "carb"))
#' res
#' res$selected
#' res$profile
#'
#' ## The same search at a heavier charge per parameter keeps fewer predictors.
#' perform_stepwise(mtcars, outcome = "mpg",
#'              predictors = c("wt", "hp", "disp", "qsec", "drat", "carb"),
#'              criterion = "BIC")$selected
#'
#' ## The selection is a set of column names, so it goes straight back into a fit.
#' fit_linear_regression(mtcars, outcome = "mpg", predictors = res$selected,
#'                       cv = FALSE)$coefficients
#'
#' ## A two-class outcome, with the reference named on its own.
#' iris2 <- iris[iris$Species != "setosa", ]
#' perform_stepwise(iris2, outcome = "Species", control_label = "versicolor",
#'              model = "logistic")
#'
#' @export
perform_stepwise <- function(data,
                         outcome,
                         predictors = NULL,
                         outcome_lv = NULL,
                         control_label = outcome_lv[1],
                         model = c("linear", "logistic"),
                         criterion = c("AIC", "BIC"),
                         direction = c("backward", "both", "forward")) {

  model <- match.arg(model)
  criterion <- match.arg(criterion)
  direction <- match.arg(direction)

  input <- sa_resolve_model_input(data, outcome, predictors)

  # `model` is the first say and the outcome is the second, which is the rule
  # `perform_rfe()` follows and for the same reason: the name of a model is already
  # an answer to the question the outcome would have been asked.
  classify <- model == "logistic" || !is.null(outcome_lv) ||
    !is.null(control_label) || !is.numeric(input$y)

  if (model == "linear" && classify) {
    stop("`model = \"linear\"` searches for the predictors of a number, and ",
         "`outcome` is a set of class labels. Use `model = \"logistic\"` for a ",
         "two-class outcome.", call. = FALSE)
  }
  if (model == "logistic" && is.numeric(input$y) &&
        length(unique(input$y)) > 2L) {
    stop("`model = \"logistic\"` classifies two classes, and `outcome` is a ",
         "numeric column taking ", length(unique(input$y)),
         " values. Use `model = \"linear\"` for a continuous outcome.",
         call. = FALSE)
  }
  if (!classify && length(unique(input$y)) == 2L) {
    message("`outcome` is numeric and takes two values, so it was searched as ",
            "a regression. Pass `control_label`, or a factor column, with ",
            "`model = \"logistic\"` to treat it as a classification.")
  }

  if (classify) {
    y <- sa_outcome_levels(input$y, outcome_lv, control_label,
                           model = "a stepwise selection")
    outcome_lv <- levels(y)
  } else {
    if (!all(is.finite(input$y))) {
      stop("`outcome` holds non-finite value(s), which a model fitted at each ",
           "step has no likelihood for.", call. = FALSE)
    }
    y <- input$y
  }

  k <- if (criterion == "BIC") log(input$n_used) else 2
  label <- sa_search_label(model, classify)

  # One grouped message for the whole procedure rather than one per stage. A
  # condition of the data — a logistic regression whose predictors separate its
  # classes perfectly, say — is raised by every model the path fits and again by
  # the ones refitted afterwards to price each candidate, and it is the same
  # condition every time.
  found <- sa_quiet_engine(
    {
      search <- sa_step_search(sa_search_frame(input$x, y), input$predictors,
                               model, direction, k)
      kept <- attr(stats::terms(search$fit), "term.labels")
      list(fit     = search$fit,
           kept    = kept,
           ranking = sa_step_ranking(search, input$predictors, kept, k,
                                     criterion),
           profile = sa_step_profile(search$fit))
    },
    label
  )

  if (length(found$kept) == 0L) {
    stop("the search walked back to the intercept: at a charge of ",
         sa_fmt_num(k, 3), " per parameter, none of the ",
         length(input$predictors), " candidate(s) lowers ", criterion,
         " by more than it costs. That is an answer rather than a failure, but ",
         "not one this function can return, since `$selected` would be empty. ",
         if (criterion == "BIC") {
           "`criterion = \"AIC\"` charges 2 per parameter instead. "
         },
         "Read it as none of these predictors being worth a coefficient on ",
         "these rows.", call. = FALSE)
  }
  ranking <- found$ranking

  design <- list(
    outcome      = input$outcome,
    outcome_type = if (classify) "two classes" else "continuous"
  )
  if (classify) {
    n_events <- sum(y == outcome_lv[2])
    design <- c(design, list(outcome_lv = outcome_lv,
                             n_events   = n_events,
                             event_rate = n_events / input$n_used))
  }

  sa_new_selection(
    analysis   = "stepwise",
    candidates = ranking$candidates,
    design     = c(design, list(
      n_obs              = input$n_obs,
      n_used             = input$n_used,
      n_dropped          = input$n_dropped,
      predictors         = input$predictors,
      dropped_predictors = input$dropped_predictors
    ), sa_design_lv(input$predictor_lv)),
    # `maximize` is recorded and not asked for: a criterion is a cost, so a caller
    # who could say otherwise could ask for the worst model by accident. This is
    # the rule `sa_rfe_metric()` follows for a resampled metric.
    parameters = list(model     = model,
                      criterion = criterion,
                      maximize  = FALSE,
                      k         = k,
                      direction = direction),
    # Most important first, which is the ranking's order rather than the formula's.
    # `sa_new_selection()` checks the two agree, so this is the one place the
    # selected set is spelled out.
    selected   = ranking$candidates[ranking$selected],
    ranking    = ranking,
    profile    = found$profile,
    resampling = NULL,
    engine     = list(
      package    = "stats",
      method     = "step",
      label      = label,
      metrics    = c("AIC", "BIC"),
      importance = paste(criterion, "increase when the predictor is left out")
    ),
    fit = found$fit
  )
}


#' Walk the path, from the model the direction starts at
#'
#' The formulas are built rather than written, for the reason `sa_rfe_formula()`
#' is: `.outcome ~ .` written out is a symbol that exists only inside a model
#' frame, which is the one thing `R CMD check` cannot tell from a typo.
#'
#' Their environment is set deliberately and is not decoration. [stats::step()]
#' moves by [stats::update()], which re-evaluates the fit's `data` argument in the
#' environment of its formula, so a formula carrying the environment
#' [stats::reformulate()] happened to be called from would send the refit looking
#' for `frame` in the caller's frame and fail there. Pointing all three at this
#' function's own frame, where `frame` is an argument, is what makes every step of
#' the path fit the same rows as the first.
#'
#' `steps` is left at [stats::step()]'s own 1000. The path is bounded by the
#' candidate count in one direction and does not run far past it in the other, so a
#' cap low enough to bind would stop the search somewhere it had not finished.
#'
#' @param frame Model frame with the outcome as `.outcome`.
#' @param predictors The candidates, which are the terms of the upper model.
#' @param model Which model is fitted at each step, already resolved.
#' @param direction Which moves are allowed, already resolved.
#' @param k Charge per parameter: 2 for AIC, `log(n)` for BIC.
#'
#' @return List with the `fit` [stats::step()] returned and the `upper` formula,
#'   which `sa_step_ranking()` needs to ask what an excluded candidate is worth.
#'
#' @keywords internal
#' @noRd
sa_step_search <- function(frame, predictors, model, direction, k) {
  upper <- stats::reformulate(predictors)
  full <- stats::reformulate(predictors, response = ".outcome")
  lower <- stats::reformulate("1", response = ".outcome")
  environment(upper) <- environment()
  environment(full) <- environment()
  environment(lower) <- environment()

  start <- if (direction == "forward") lower else full
  fitted <- switch(
    model,
    linear = stats::lm(start, data = frame),
    logistic = stats::glm(start, data = frame, family = stats::binomial()),
    stop("internal error: unhandled `model` ", model, ".", call. = FALSE)
  )
  fit <- stats::step(fitted, scope = list(lower = lower, upper = full),
                     direction = direction, trace = FALSE, k = k,
                     keep = sa_step_keep)

  list(fit = fit, upper = upper)
}


#' What is recorded about every model the path visits
#'
#' [stats::step()] calls this once per accepted model, which is once per row of the
#' `$anova` table it returns, and that is the only way to the term count of a model
#' the search has already moved on from. Reading it off the `Step` labels instead
#' would mean parsing `"- carb"` back into a name, which a term with a space in it
#' breaks.
#'
#' Both criteria are recorded whichever one is searching, since they are two
#' charges against the same likelihood and the fit that has one has the other.
#'
#' @param fit The model at this step.
#' @param aic Its criterion on [stats::extractAIC()]'s scale, which `step()` passes
#'   in and which is not what is kept; see the note on scales in `perform_stepwise()`.
#'
#' @return A named list, which [stats::step()] arranges into a matrix with one
#'   column per step.
#'
#' @keywords internal
#' @noRd
sa_step_keep <- function(fit, aic) {
  list(n_vars = length(attr(stats::terms(fit), "term.labels")),
       AIC    = stats::AIC(fit),
       BIC    = stats::BIC(fit))
}


#' The profile table, one row per step of the path
#'
#' `n_vars` counts predictors rather than coefficients, so a factor counts once,
#' and it is the field name the contract shares with `perform_rfe()`, where the row
#' axis is a subset size instead of a step. `chosen` is the last row: the search
#' stops when no move improves the criterion, so where it stopped is what it chose.
#'
#' @param fit The [stats::step()] result.
#'
#' @return data.frame with `n_vars`, `AIC`, `BIC`, `step` and `chosen`.
#'
#' @keywords internal
#' @noRd
sa_step_profile <- function(fit) {
  path <- fit$keep
  steps <- as.character(fit$anova$Step)
  if (is.null(path) || ncol(path) != length(steps)) {
    stop("internal error: `step()` recorded ", if (is.null(path)) 0L else
      ncol(path), " model(s) and ", length(steps), " path row(s).",
      call. = FALSE)
  }

  out <- data.frame(
    n_vars = as.integer(unlist(path["n_vars", ])),
    AIC    = as.numeric(unlist(path["AIC", ])),
    BIC    = as.numeric(unlist(path["BIC", ])),
    # The first row is the model the search started at, which no move reached.
    step   = trimws(steps),
    stringsAsFactors = FALSE
  )
  out$chosen <- seq_len(nrow(out)) == nrow(out)
  rownames(out) <- NULL
  out
}


#' The ranking table, one row per candidate
#'
#' One number for both groups of candidates, which is what makes the table
#' readable as a single column: what the criterion would be with this predictor
#' left out of the selected model, minus what it is with it in. A predictor the
#' search kept is worth the rise that dropping it would cause, so it is positive; a
#' predictor the search left out would raise the criterion by being added, so its
#' number is the same difference and comes out negative. The sign is the search's
#' own decision about it and the size is by how much, which puts the selection at
#' the top of the table with no separate sort.
#'
#' [stats::drop1()] and [stats::add1()] report those models on
#' [stats::extractAIC()]'s scale rather than [stats::AIC()]'s. Both are differences
#' taken within one table against its own `<none>` row, and the two scales differ
#' by a constant across models on the same rows, so the difference is one number
#' either way and comparable with the columns of `profile`.
#'
#' Ties are broken by name, which is what makes the ranking of an otherwise
#' symmetric set of candidates reproducible, and a candidate neither table could
#' answer for sorts last on an `NA` rather than on a fabricated zero.
#'
#' @param search The list `sa_step_search()` returned.
#' @param predictors The candidates, in the order they arrived.
#' @param kept The terms of the selected model.
#' @param k Charge per parameter, the same one the search moved by.
#' @param criterion Which criterion that charge belongs to, for the error a table
#'   with no criterion column would otherwise raise further in.
#'
#' @return data.frame with `candidates`, `estimate`, `rank` and `selected`.
#'
#' @keywords internal
#' @noRd
sa_step_ranking <- function(search, predictors, kept, k, criterion) {
  estimate <- stats::setNames(rep(NA_real_, length(predictors)), predictors)

  dropped <- sa_step_criterion(stats::drop1(search$fit, k = k), criterion)
  estimate[kept] <- dropped[kept] - dropped[["<none>"]]

  outside <- setdiff(predictors, kept)
  if (length(outside) > 0L) {
    # The scope is the upper model rather than the excluded terms alone, because
    # `add.scope()` reads it as the model to grow towards and requires the current
    # terms to be inside it.
    added <- sa_step_criterion(
      stats::add1(search$fit, scope = search$upper, k = k), criterion
    )
    estimate[outside] <- added[["<none>"]] - added[outside]
  }

  candidates <- sort(predictors)
  at <- order(estimate[candidates], decreasing = TRUE, na.last = TRUE)

  out <- data.frame(
    candidates = candidates[at],
    estimate   = unname(estimate[candidates][at]),
    rank       = seq_along(candidates),
    selected   = candidates[at] %in% kept,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}


#' The criterion column of a `drop1()` or `add1()` table
#'
#' Named `AIC` whatever `k` was, since the column holds the criterion the charge
#' produces and `k` is the charge. What is read out of it are differences against
#' the `<none>` row, so the label is the only thing about it that mentions AIC.
#'
#' @param tab A [stats::drop1()] or [stats::add1()] table.
#' @param criterion Which criterion was asked for, for the message.
#'
#' @return Named numeric vector, one entry per row of `tab`.
#'
#' @keywords internal
#' @noRd
sa_step_criterion <- function(tab, criterion) {
  if (is.null(tab[["AIC"]])) {
    stop("internal error: the single-term table holds no criterion column to ",
         "read ", criterion, " out of; it has ",
         paste(names(tab), collapse = ", "), ".", call. = FALSE)
  }
  stats::setNames(as.numeric(tab[["AIC"]]), rownames(tab))
}
