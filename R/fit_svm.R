# The second model here with no coefficients, and it arrives at that from the
# other side than `fit_rf()` does. A forest has too many numbers per predictor to
# report one; a radial-kernel machine has none at all. What it holds is a set of
# support vectors and a weight for each, and those are points in the data rather
# than directions in the predictors: the surface they describe lives in a space the
# kernel never builds, so there is nothing to write in the column that would say
# what one unit of a predictor is worth.
#
# What can still be asked of any fitted model is what it would lose without a
# given predictor, and that is what `estimate` carries here. The measurement is
# the same idea as a forest's permutation importance and is not the same number:
# a forest permutes within the rows each tree left out of its bootstrap sample, and
# a machine has no such rows, so this is measured on the rows it was fitted to.
#
# The other departure is the design matrix. `kernlab` measures a distance between
# rows, so it takes a numeric matrix and `caret` hands it `as.matrix()` of whatever
# it was given — the same trap `fit_elastic_net()` avoids, and avoided the same
# way, by dummy coding here so that a factor arrives as its levels rather than as
# its integer codes.

#' Fit a support vector machine
#'
#' Fits a support vector machine with a radial basis kernel on a set of
#' predictors, choosing the cost and the kernel width by cross-validation, and
#' scores the chosen machine on the same folds. Which outcome it is decides which
#' machine is fitted: a numeric outcome is a support vector regression and a
#' two-class one a classification, as in [fit_elastic_net()] and [fit_rf()], so one
#' function covers what [fit_linear_regression()] and [fit_logistic_regression()]
#' cover between them.
#'
#' What comes back in place of a coefficient table is an importance table: one row
#' per term, in descending order of how much the machine lost when that term's
#' values were shuffled among the rows. It says which terms carried the fit and not
#' which way they pushed it, since a kernel machine bends its surface in different
#' directions in different regions of the data.
#'
#' @details
#' \describe{
#'   \item{`estimate` is permutation importance, and there is no standard error
#'     beside it}{The fitted machine is scored on the rows it was fitted to, then
#'     scored again with one term's values shuffled among those rows, and
#'     `estimate` is the mean loss over `n_permute` shuffles. The metric is the one
#'     the resampling tuned on: the rise in `RMSE` for a regression and the fall in
#'     `Accuracy` for a classification, so a larger value is a term the machine
#'     needed either way. It is not an estimate of anything a distribution is
#'     defined over, so the inference columns are absent from the table rather than
#'     present and `NA`, exactly as they are for [fit_elastic_net()] and
#'     [fit_rf()], and `is.null(fit$coefficients$pval)` tells the two kinds of
#'     table apart.}
#'   \item{The importance is measured in sample, and a forest's is not}{This is the
#'     one number here that is weaker than [fit_rf()]'s counterpart. Every tree of
#'     a forest leaves about a third of the rows out of its bootstrap sample, so
#'     the permutation can be scored on rows that tree had not seen; a machine is
#'     fitted on all of them at once and has no such rows. A term the machine
#'     fitted noise with therefore earns some importance here that it would not
#'     earn on held-out rows. `performance` is the number that was measured on rows
#'     the procedure had not seen, and it describes the whole machine rather than
#'     one term of it.}
#'   \item{A negative importance is a value, not a gap}{Shuffling a term that
#'     carried nothing can leave the machine very slightly better than it was, and
#'     the mean loss then comes out below zero. It reads as it should: this term
#'     did no better than its own permutation.}
#'   \item{Terms are not predictors}{A factor or character predictor with `k`
#'     levels becomes `k - 1` terms named after the levels, and the coding is the
#'     one [stats::lm()] would have used, so the tables can be read side by side.
#'     Each of those terms is shuffled on its own, so one level of a factor can
#'     matter to the machine and another not. [fit_rf()] is the model where the
#'     terms are the predictors themselves, because a tree splits a factor
#'     directly.}
#'   \item{The predictors are centred and scaled}{A radial kernel reads one
#'     distance between two rows, so a predictor measured in millimetres would
#'     dominate the same predictor in metres and the kernel width would mean
#'     something different along each axis. Every term is centred and scaled before
#'     the machine sees it, and [predict.sa_model()] applies the same centre and
#'     scale to new rows, so this is not something to undo afterwards. It is also
#'     why `sigma` is a width on the standardised scale rather than on the scale
#'     the columns arrived in.}
#'   \item{Cross-validation chooses `C` and `sigma`}{Both are chosen by the
#'     resampled metric, which is `RMSE` for a regression and `Accuracy` for a
#'     classification, and the final machine is then fitted on all usable rows at
#'     the pair that won, as in [fit_elastic_net()]. `performance` therefore holds
#'     one row per candidate and `parameters` the pair that was fitted. With
#'     `cv = FALSE` there is nothing to choose between, so `C` and `sigma` must
#'     name exactly one candidate.}
#' }
#'
#' `sigma = NULL` is a width read off the data rather than a grid, since the
#' distances between the rows are what fixes the scale a kernel is comparing them
#' on. [kernlab::sigest()] reports three widths, from the 10th, 50th and 90th
#' percentile of those distances, and the middle one is used. It samples pairs of
#' rows to do it, so this is one of the things `seed` fixes.
#'
#' Rows with a missing value in the outcome or in any predictor are dropped before
#' the folds are drawn rather than inside each fit, and `design$n_dropped` reports
#' how many went. A predictor that takes a single value is left out with a message.
#'
#' One predictor is enough. A distance can be measured along one axis as readily as
#' along ten, so a single column is a machine here rather than the error
#' [fit_elastic_net()] raises for want of a budget to divide.
#'
#' For a two-class outcome the direction is `outcome_lv`, read as
#' [fit_logistic_regression()] reads it: the first level is the reference, so
#' `sensitivity` is the share of `outcome_lv[2]` the machine found, `specificity`
#' the share of `outcome_lv[1]` it left alone, and
#' `predict(model, newdata, type = "response")` is the probability of
#' `outcome_lv[2]`. Those probabilities are Platt's, a logistic curve fitted to the
#' decision values by an internal cross-validation of `kernlab`'s own, which is a
#' second thing `seed` fixes and a reason two machines fitted without one predict
#' slightly different probabilities from the same decisions. The importance table
#' has no direction to report, which is why there is no `odds_ratio` column beside
#' it. `control_label` names the reference on its own and is enough by itself to
#' make a column of zeroes and ones a classification; naming it alongside an
#' `outcome_lv` that puts the other class first is an error, as it is in
#' [fit_logistic_regression()].
#'
#' @param data A data.frame (or matrix) in wide format, one row per observation.
#'   Typically the training half of a [split_data()] result.
#' @param outcome The outcome, either the name of a column of `data` or a vector
#'   with one entry per row. A numeric outcome is fitted as a regression, and a
#'   factor, character or logical one as a two-class classification.
#' @param predictors Column names to fit on, or `NULL` for every column of `data`
#'   except the outcome. Numeric, logical, factor and character columns are
#'   accepted; a column that takes a single value is left out with a message,
#'   since it cannot contribute.
#' @param outcome_lv The two classes, reference first, so that the reported
#'   sensitivity and the predicted probability are about the second one. Supplying
#'   it also states that the outcome is to be classified, which is the one thing a
#'   numeric column holding two values cannot say on its own. `NULL` sorts the
#'   classes of an outcome that is not numeric, and leaves a numeric one to the
#'   regression.
#' @param control_label The reference class on its own, for when the other one
#'   needs no saying. Defaults to `outcome_lv[1]`, so a call that names one of the
#'   two names the reference either way; naming both and disagreeing is an error.
#' @param C Cost of violating the margin, one value to score or several to choose
#'   between. Larger fits the training rows harder and generalises less.
#' @param sigma Width of the radial kernel on the standardised scale, one value or
#'   several. Larger makes the kernel narrower, so each support vector reaches
#'   fewer rows. `NULL` is the middle of the three widths [kernlab::sigest()]
#'   reports for these rows.
#' @param n_permute Shuffles per term behind each importance. More is steadier and
#'   costs one prediction of every row per shuffle.
#' @param cv Whether to cross-validate. `FALSE` fits the single candidate the grid
#'   names, once, and reports no resampled performance.
#' @param cv_method Resampling scheme: `"repeated_kfold"` for `n_repeat` runs of
#'   `n_fold`-fold cross-validation, `"kfold"` for one, or `"loocv"` for
#'   leave-one-out.
#' @param n_fold Folds per run, used by `"repeated_kfold"` and `"kfold"`.
#' @param n_repeat Number of runs, used by `"repeated_kfold"`.
#' @param seed Seed for the kernel width, the class probabilities, the shuffles
#'   behind the importance and the fold assignment, or `NULL` to use the stream as
#'   it stands. Supplying one does not disturb the caller: the previous random
#'   number state is put back when the function returns.
#'
#' @return An object of class `sa_model`, the same eleven elements
#'   [fit_linear_regression()] returns, with these differences:
#'
#'   \describe{
#'     \item{`analysis`}{`"svm"`.}
#'     \item{`terms`}{The model terms, in descending order of importance.}
#'     \item{`design`}{Holds `outcome_lv`, `n_events` and `event_rate` for a
#'       two-class outcome, as [fit_logistic_regression()] does, and neither for a
#'       continuous one.}
#'     \item{`parameters`}{Holds `kernel`, the `C` and `sigma` that were fitted
#'       rather than the grids that were searched, `n_candidates`, how many pairs
#'       were scored, and `n_permute`. The grids themselves are the rows of
#'       `performance`.}
#'     \item{`coefficients`}{`estimate`, the permutation importance. There is no
#'       intercept row, no `odds_ratio` and no inference column.}
#'     \item{`fit_stats`}{Measured on the rows the machine was fitted to:
#'       `r_squared`, `rmse` and `mae` for a regression, and `accuracy`, `error`,
#'       `kappa`, `sensitivity` and `specificity` for a classification, each with
#'       `n_support_vector` and `support_vector_rate`, the rows the machine kept as
#'       support vectors and their share of the rows it was fitted to. These
#'       describe the fit rather than what it would do next; `performance` is the
#'       one that was measured on rows the machine had not seen.}
#'     \item{`performance`}{One row per candidate, the chosen one being the row
#'       that matches `parameters$C` and `parameters$sigma`, or `NULL` when
#'       `cv = FALSE`.}
#'   }
#'
#' @seealso [fit_rf()] for the other model that answers with importance rather
#'   than with coefficients, [fit_linear_regression()],
#'   [fit_logistic_regression()] and [fit_elastic_net()] for the models that do
#'   report an effect per predictor, [split_data()], which defines the rows this is
#'   fitted on, [coef.sa_model()] for the importance table, and
#'   [predict.sa_model()] for predicting the rows it was not fitted on.
#'
#' @examples
#' ## Fitted on every row, scored on rows each fold had not seen. The table is in
#' ## descending order, so its first rows are the terms the machine leaned on.
#' fit <- fit_svm(mtcars, outcome = "mpg", predictors = c("wt", "hp", "disp"),
#'                C = 1, cv_method = "kfold", n_fold = 5, seed = 1)
#' fit
#' fit$coefficients
#' fit$fit_stats
#'
#' ## `C` and `sigma` are chosen by the resampling: a single pair is scored,
#' ## several are compared.
#' tuned <- fit_svm(mtcars, outcome = "mpg", predictors = c("wt", "hp", "disp"),
#'                  C = c(0.5, 2, 8), cv_method = "kfold", n_fold = 5, seed = 1)
#' tuned$parameters[c("C", "sigma", "n_candidates")]
#' tuned$performance[c("C", "RMSE")]
#'
#' ## A two-class outcome is a classification, and `outcome_lv` fixes which class
#' ## the sensitivity and `type = "response"` are about.
#' iris2 <- iris[iris$Species != "setosa", ]
#' clf <- fit_svm(iris2, outcome = "Species",
#'                outcome_lv = c("versicolor", "virginica"), C = 1, cv = FALSE,
#'                seed = 1)
#' clf$fit_stats[c("accuracy", "sensitivity", "specificity")]
#' head(predict(clf, newdata = iris2, type = "response"))
#'
#' ## Scored against an answer that is known: a simulated regression plants some
#' ## coefficients at exactly zero, so the null terms are the ones the importance
#' ## ought to leave at the bottom of the table.
#' sim <- simulate_regression(seed = 1)
#' svm <- do.call(fit_svm, c(sim$args, list(C = 1, cv = FALSE, seed = 1)))
#' scored <- merge(svm$coefficients, sim$truth_term, by = "terms")
#' tapply(scored$estimate, scored$beta != 0, mean)
#'
#' @export
fit_svm <- function(data,
                    outcome,
                    predictors = NULL,
                    outcome_lv = NULL,
                    control_label = outcome_lv[1],
                    C = 2^seq(-5, 10, by = 2),
                    sigma = NULL,
                    n_permute = 10,
                    cv = TRUE,
                    cv_method = c("repeated_kfold", "kfold", "loocv"),
                    n_fold = 5,
                    n_repeat = 5,
                    seed = NULL) {

  cv_method <- match.arg(cv_method)
  n_permute <- sa_check_count(n_permute, "n_permute", 1)

  input <- sa_resolve_model_input(data, outcome, predictors)

  # One argument, two machines, and the outcome decides which, the same way it
  # does in `fit_elastic_net()` and `fit_rf()`. A numeric column is a regression,
  # since that is what a number usually is, and anything else is a set of class
  # labels. `outcome_lv` overrules the guess, and so does `control_label` on its
  # own, which is the only way to say that a column of zeroes and ones is two
  # classes rather than two numbers.
  classify <- !is.null(outcome_lv) || !is.null(control_label) ||
    !is.numeric(input$y)
  if (!classify && length(unique(input$y)) == 2L) {
    message("`outcome` is numeric and takes two values, so it was fitted as a ",
            "regression. Pass `outcome_lv` or `control_label`, or a factor ",
            "column, to model it as a classification.")
  }

  if (classify) {
    y <- sa_outcome_levels(input$y, outcome_lv, control_label,
                           model = "a support vector machine")
    outcome_lv <- levels(y)
  } else {
    if (!all(is.finite(input$y))) {
      stop("`outcome` holds non-finite value(s), which the loss of a support ",
           "vector regression has no residual for.", call. = FALSE)
    }
    y <- input$y
  }

  # Dummy coded here rather than left to the engine, for the reason it is in
  # `fit_elastic_net()`: what `caret` hands `kernlab` is `as.matrix()` of the
  # frame, which turns a factor into its integer codes and fits without
  # complaint. See `sa_design_matrix()`.
  x <- sa_design_matrix(input$x)
  ctrl <- sa_train_control(cv, cv_method, n_fold, n_repeat, input$n_used)

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  sigma <- sa_svm_sigma(sigma, x)
  grid <- sa_svm_grid(C, sigma, cv)

  # `prob.model` reaches `kernlab::ksvm()` through `caret`'s `...`, and passing it
  # is what makes `type = "response"` and `type = "prob"` answerable at all on a
  # classification: without it the machine reports a class and a decision value and
  # no probability. It is `classify` rather than `TRUE` because `caret` reads
  # `prob.model` out of `trainControl(classProbs = )` when the argument is absent,
  # and this package does not ask the caller for that; a regression then gets the
  # `FALSE` that is `ksvm`'s own default for it.
  fit <- sa_fit_engine(
    caret::train(x = x, y = y, method = "svmRadial", trControl = ctrl$control,
                 tuneGrid = grid, preProcess = c("center", "scale"),
                 prob.model = classify),
    "Support vector machine"
  )
  model <- fit$finalModel

  coefs <- sa_svm_importance(fit, x, y, classify, n_permute)
  fitted_value <- stats::predict(fit, newdata = x)
  n_support_vector <- kernlab::nSV(model)
  fit_stats <- c(
    if (classify) {
      sa_svm_class_stats(fitted_value, y, outcome_lv)
    } else {
      sa_svm_reg_stats(fitted_value, y)
    },
    list(n_support_vector    = n_support_vector,
         support_vector_rate = n_support_vector / input$n_used)
  )
  fit_stats <- lapply(fit_stats, function(v) unname(as.numeric(v)))

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

  sa_new_model(
    analysis = "svm",
    terms    = coefs$terms,
    design   = c(design, list(
      n_obs              = input$n_obs,
      n_used             = input$n_used,
      n_dropped          = input$n_dropped,
      predictors         = input$predictors,
      dropped_predictors = input$dropped_predictors
    ), sa_design_lv(input$predictor_lv)),
    # The pair that ran, not the grids that were asked for: `performance` holds
    # every candidate, so recording the grids here as well would say the same
    # thing twice and leave two places for it to be wrong.
    parameters = list(
      kernel       = "radial",
      C            = unname(fit$bestTune$C),
      sigma        = unname(fit$bestTune$sigma),
      n_candidates = nrow(grid),
      n_permute    = n_permute,
      cv           = cv,
      cv_method    = ctrl$cv_method,
      n_fold       = ctrl$n_fold,
      n_repeat     = ctrl$n_repeat,
      seed         = seed
    ),
    coefficients = coefs,
    fit_stats    = fit_stats,
    performance  = sa_model_frame(fit$results),
    resampling   = sa_model_frame(fit$resample),
    engine       = list(
      package = "caret",
      method  = "svmRadial",
      kernel  = "radial",
      label   = paste("Support vector machine",
                      if (classify) "classification" else "regression",
                      "(radial basis kernel)"),
      metrics = fit$perfNames,
      # A machine fitted on a design matrix reads one by position, so the columns
      # it was given have to be rebuilt in that order to predict new rows.
      # `caret` records no `coefnames` for a fit it was handed a matrix, which is
      # why the result keeps them itself; `predict.sa_model()` is what reads this.
      x_names = colnames(x)
    ),
    fit = fit
  )
}


#' The kernel width, read off the rows when the caller names none
#'
#' A radial kernel compares two rows by one distance, so the width worth starting
#' from is a property of the distances the data holds rather than a constant.
#' [kernlab::sigest()] reports the widths matching the 10th, 50th and 90th
#' percentile of the squared distances between pairs of rows, and the middle one is
#' taken: the narrower one fits a surface that turns inside the bulk of the data
#' and the wider one is close to a straight one.
#'
#' `scaled = TRUE` because the machine sees standardised columns. Measuring the
#' distances on the raw scale would give a width for a space the kernel is never
#' asked about, which for predictors in different units is not a rescaling of the
#' right answer but a different one.
#'
#' It samples pairs rather than reading all of them, so this is a random quantity
#' and two calls without a seed choose slightly different widths.
#'
#' @param sigma What the caller passed, or `NULL`.
#' @param x The design matrix, already dummy coded.
#'
#' @return `sigma` unchanged, or the single width read off `x`.
#'
#' @keywords internal
#' @noRd
sa_svm_sigma <- function(sigma, x) {
  if (!is.null(sigma)) {
    return(sigma)
  }
  estimated <- kernlab::sigest(x, scaled = TRUE)
  # Named "50%" by `kernlab`, and a name is not something the result should carry
  # into `parameters` from an engine's internals.
  unname(estimated[2L])
}


#' The importance table a machine answers with instead of coefficients
#'
#' What a fitted model can always be asked is what it would lose without a given
#' term, and shuffling that term's values among the rows is how the question is put
#' to a model that has no coefficient to remove. The machine is not refitted, so
#' this is what the surface it already found is worth per term rather than what a
#' machine fitted without that term would have found.
#'
#' The metric is the one the resampling tuned on, so that a term's importance and
#' the machine's own score are in the same units and `performance` can be read
#' beside this table. A regression is scored by `RMSE` and a classification by
#' `Accuracy`, and the sign is turned so that a large `estimate` is a term that was
#' needed in both cases.
#'
#' Descending rather than in the order the columns arrived. `terms` is the row
#' order of every table in the result, and for this model the order worth reading
#' first is the most important term; the order the columns were read in is
#' `engine$x_names`, which keeps it.
#'
#' @param fit The fitted engine object, predicted through rather than refitted.
#' @param x The design matrix the machine was fitted to, on the scale it was
#'   handed: `caret` centres and scales inside its own `predict()`, so a shuffled
#'   column has to go in the same way the original did.
#' @param observed The outcome as the fit received it.
#' @param classify Whether the outcome was classified, which is what decides which
#'   metric the loss is measured in.
#' @param n_permute Shuffles per term, averaged.
#'
#' @return data.frame with `terms` and `estimate`.
#'
#' @keywords internal
#' @noRd
sa_svm_importance <- function(fit, x, observed, classify, n_permute) {
  score <- if (classify) {
    function(prediction) {
      mean(as.character(prediction) == as.character(observed))
    }
  } else {
    function(prediction) sqrt(mean((observed - prediction)^2))
  }

  as_fitted <- score(stats::predict(fit, newdata = x))
  loss <- vapply(colnames(x), function(nm) {
    shuffled <- vapply(seq_len(n_permute), function(i) {
      permuted <- x
      permuted[, nm] <- sample(permuted[, nm])
      score(stats::predict(fit, newdata = permuted))
    }, numeric(1))
    # Accuracy falls when a term mattered and RMSE rises, so one of the two is
    # subtracted the other way round and both come out positive for a term the
    # machine needed.
    if (classify) as_fitted - mean(shuffled) else mean(shuffled) - as_fitted
  }, numeric(1))

  out <- data.frame(
    terms    = colnames(x),
    estimate = unname(loss),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$estimate, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}


#' Goodness of fit of a machine on the rows it was fitted to
#'
#' The same arithmetic `fit_elastic_net()` reports, and for the same reason: a
#' `ksvm` object holds no summary of itself, only the support vectors it kept and
#' the weight on each. These describe the surface that was found rather than what
#' it would do next, which is what `performance` is for. [fit_rf()] is the one
#' model in this family whose in-sample counterpart is honest, because a forest has
#' rows every tree left out.
#'
#' @param fitted_value The machine's prediction for the rows it was fitted to.
#' @param observed The outcome as the fit received it.
#'
#' @return Named list of scalars.
#'
#' @keywords internal
#' @noRd
sa_svm_reg_stats <- function(fitted_value, observed) {
  residual <- observed - fitted_value
  sse <- sum(residual^2)
  sst <- sum((observed - mean(observed))^2)

  list(
    r_squared = if (sst > 0) 1 - sse / sst else NA_real_,
    rmse      = sqrt(mean(residual^2)),
    mae       = mean(abs(residual))
  )
}


#' @rdname sa_svm_reg_stats
#' @param outcome_lv The two classes, reference first, so that the sensitivity is
#'   the one of `outcome_lv[2]` and the whole result reads in one direction.
#' @keywords internal
#' @noRd
sa_svm_class_stats <- function(fitted_value, observed, outcome_lv) {
  called <- as.character(fitted_value) == outcome_lv[2]
  positive <- as.character(observed) == outcome_lv[2]

  accuracy <- mean(called == positive)
  # Cohen's kappa on the two-by-two table: the agreement above what the observed
  # and the predicted labels would reach from their marginals alone.
  expected <- mean(called) * mean(positive) + mean(!called) * mean(!positive)

  list(
    accuracy    = accuracy,
    error       = 1 - accuracy,
    kappa       = if (expected < 1) {
      (accuracy - expected) / (1 - expected)
    } else {
      NA_real_
    },
    sensitivity = if (any(positive)) mean(called[positive]) else NA_real_,
    specificity = if (any(!positive)) mean(!called[!positive]) else NA_real_
  )
}
