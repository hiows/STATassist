# Internal helpers shared by the model fitting functions, and by the searches that
# fit a model at every step. A model function's own body should be the part that is
# specific to that model: which engine runs and what its summary means. Everything
# before that is the same question every time — which columns are the predictors,
# which rows can be used, and how the resampling scheme was described — so it is
# answered once here.
#
# The resampling scheme is the reason this file exists rather than each function
# building its own `trainControl()`. The four schemes were already written out
# twice in the draft and would have been written out five more times by the
# elastic net, random forest and gradient boosting functions.

#' Resolve the resampling arguments into the scheme caret names
#'
#' The three schemes differ in which arguments they consume, and the ones they do
#' not consume are silently ignored by the control object they end up in. LOOCV
#' has no fold count and no repeats and plain k-fold has no repeats, so reporting
#' `n_fold = 5` for a LOOCV run would be a plausible-looking record of something
#' that never happened. What the scheme actually used comes back beside caret's
#' own name for it, and it is those values that a result records.
#'
#' This is the part `sa_train_control()` and `sa_rfe_control()` share. The two
#' build different control objects — one resamples a fixed set of predictors and
#' the other resamples the elimination itself — but the question of what
#' `cv_method = "kfold", n_fold = 5` means is the same one, so it is answered
#' once.
#'
#' @param cv_method Scheme name, already resolved by [match.arg()].
#' @param n_fold,n_repeat Fold count and number of repeats, as requested.
#' @param n_obs Rows available, used to reject more folds than observations.
#'
#' @return List with caret's `method` string and the `cv_method` / `n_fold` /
#'   `n_repeat` the scheme actually used, `NA` where it uses none of them.
#'
#' @keywords internal
#' @noRd
sa_resample_scheme <- function(cv_method, n_fold, n_repeat, n_obs) {
  n_fold <- sa_check_count(n_fold, "n_fold", 2)
  n_repeat <- sa_check_count(n_repeat, "n_repeat", 1)

  if (cv_method != "loocv" && n_fold > n_obs) {
    stop("`n_fold` = ", n_fold, " exceeds the ", n_obs,
         " usable observation(s), so a fold would be empty. Lower `n_fold` or ",
         "use `cv_method = \"loocv\"`.", call. = FALSE)
  }

  switch(
    cv_method,
    repeated_kfold = list(method = "repeatedcv", cv_method = cv_method,
                          n_fold = n_fold, n_repeat = n_repeat),
    kfold = list(method = "cv", cv_method = cv_method,
                 n_fold = n_fold, n_repeat = NA_integer_),
    loocv = list(method = "LOOCV", cv_method = cv_method,
                 n_fold = NA_integer_, n_repeat = NA_integer_),
    stop("internal error: unhandled `cv_method` ", cv_method, ".",
         call. = FALSE)
  )
}


#' Turn the resampling arguments into a caret control object
#'
#' `cv = FALSE` is the fourth scheme, the one with no folds and no repeats at
#' all, and it is this function's own rather than `sa_resample_scheme()`'s: a
#' model can be fitted once and reported without a resampled score, which is not
#' a thing an elimination can be.
#'
#' Every argument is validated whether or not the chosen scheme reads it. A
#' rejected value would otherwise depend on `cv_method`, which is the kind of
#' conditional strictness that is impossible to guess from the outside.
#'
#' @param cv Whether to resample at all.
#' @param cv_method Scheme name, already resolved by [match.arg()].
#' @param n_fold,n_repeat Fold count and number of repeats, as requested.
#' @param n_obs Rows available, used to reject more folds than observations.
#'
#' @return List with `control`, and the `cv_method` / `n_fold` / `n_repeat` the
#'   scheme actually used, `NA` where it uses none of them.
#'
#' @keywords internal
#' @noRd
sa_train_control <- function(cv, cv_method, n_fold, n_repeat, n_obs) {
  sa_check_flag(cv, "cv")
  n_fold <- sa_check_count(n_fold, "n_fold", 2)
  n_repeat <- sa_check_count(n_repeat, "n_repeat", 1)

  if (!cv) {
    return(list(control = caret::trainControl(method = "none"),
                cv_method = NA_character_, n_fold = NA_integer_,
                n_repeat = NA_integer_))
  }

  scheme <- sa_resample_scheme(cv_method, n_fold, n_repeat, n_obs)
  control <- switch(
    scheme$method,
    repeatedcv = caret::trainControl(method = "repeatedcv",
                                     number = scheme$n_fold,
                                     repeats = scheme$n_repeat),
    cv = caret::trainControl(method = "cv", number = scheme$n_fold),
    LOOCV = caret::trainControl(method = "LOOCV")
  )
  list(control = control, cv_method = scheme$cv_method,
       n_fold = scheme$n_fold, n_repeat = scheme$n_repeat)
}


#' Turn the resampling arguments into a caret elimination control object
#'
#' The same three schemes as `sa_train_control()`, wrapped around
#' [caret::rfeControl()] instead, which resamples the elimination rather than one
#' fit: the ranking is recomputed inside every fold, so a predictor that looks
#' useful only on the rows it was ranked on is scored on rows it was not.
#'
#' `rerank` is left at `FALSE`, caret's default. Re-ranking after every drop
#' refits the model once per remaining predictor per fold, which is a different
#' and far slower procedure than the one this function is named after.
#'
#' @param funcs The `functions` list the elimination fits and ranks with.
#' @param cv_method Scheme name, already resolved by [match.arg()].
#' @param n_fold,n_repeat Fold count and number of repeats, as requested.
#' @param n_obs Rows available, used to reject more folds than observations.
#'
#' @return List with `control`, and the `cv_method` / `n_fold` / `n_repeat` the
#'   scheme actually used, `NA` where it uses none of them.
#'
#' @keywords internal
#' @noRd
sa_rfe_control <- function(funcs, cv_method, n_fold, n_repeat, n_obs) {
  scheme <- sa_resample_scheme(cv_method, n_fold, n_repeat, n_obs)
  control <- switch(
    scheme$method,
    repeatedcv = caret::rfeControl(functions = funcs, method = "repeatedcv",
                                   number = scheme$n_fold,
                                   repeats = scheme$n_repeat),
    cv = caret::rfeControl(functions = funcs, method = "cv",
                           number = scheme$n_fold),
    LOOCV = caret::rfeControl(functions = funcs, method = "LOOCV")
  )
  list(control = control, cv_method = scheme$cv_method,
       n_fold = scheme$n_fold, n_repeat = scheme$n_repeat)
}


#' The model frame a search fits on
#'
#' The outcome column is called `.outcome`, which is `caret`'s own name for it, so
#' that a predictor called `y` is a predictor rather than a collision with the
#' formula. Both searches use the name: `perform_rfe()` because that is what its
#' `functions` list is handed, and `perform_stepwise()` so that the two build the same
#' frame out of the same predictors.
#'
#' @param x Predictor columns, as a data.frame or a matrix.
#' @param y The outcome of those rows, or `NULL` when the frame is being built to
#'   predict on.
#'
#' @return data.frame of `x`, with `.outcome` appended when `y` was given.
#'
#' @keywords internal
#' @noRd
sa_search_frame <- function(x, y = NULL) {
  out <- if (is.data.frame(x)) x else as.data.frame(x, stringsAsFactors = TRUE)
  if (!is.null(y)) {
    out[[".outcome"]] <- y
  }
  out
}


#' What to call the model inside a search, in a message and in `engine$label`
#'
#' @param model Which model the search fits, already resolved.
#' @param classify Whether the outcome is being classified, which is the only
#'   thing that distinguishes the two forests.
#'
#' @keywords internal
#' @noRd
sa_search_label <- function(model, classify) {
  switch(
    model,
    linear = "Linear regression",
    logistic = "Binomial logistic regression",
    rf = paste("Random forest", if (classify) "classification" else "regression")
  )
}


#' Turn the penalty arguments into a caret tuning grid
#'
#' `penalty` is a name for a corner of the same model: the elastic net penalty is
#' a mixture of the L1 and L2 ones, and `alpha` is the mixing weight, so a lasso
#' is `alpha = 1` and a ridge is `alpha = 0`. Naming the corner rather than the
#' number means the two cases cannot be asked for wrongly, and the grid the
#' resampling searched is what the result records.
#'
#' `alpha` is validated even when `penalty` fixes it, for the same reason
#' `sa_train_control()` validates the fold count of a scheme with no folds.
#'
#' The one-row rule is `caret`'s. With `trainControl(method = "none")` there is no
#' resampling to choose between candidates, so `caret` aborts on a grid with more
#' than one row. Its message names `tuneGrid`, which is an argument this package
#' does not have, so the condition is caught here and said in terms of the
#' arguments that were passed.
#'
#' @param penalty Penalty name, already resolved by [match.arg()].
#' @param alpha,lambda The mixing weights and penalty sizes to score.
#' @param cv Whether anything is being resampled.
#'
#' @return data.frame with one row per candidate, columns `alpha` and `lambda`.
#'
#' @keywords internal
#' @noRd
sa_enet_grid <- function(penalty, alpha, lambda, cv) {
  sa_check_num_vector(alpha, "alpha", 0, 1)
  sa_check_num_vector(lambda, "lambda", 0)

  alpha <- switch(
    penalty,
    lasso = 1,
    ridge = 0,
    elastic_net = unique(alpha),
    stop("internal error: unhandled `penalty` ", penalty, ".", call. = FALSE)
  )
  grid <- expand.grid(alpha = alpha, lambda = unique(lambda))

  if (!cv && nrow(grid) > 1L) {
    stop("`cv = FALSE` fits one model, so the grid must hold one candidate, ",
         "but `alpha` and `lambda` give ", nrow(grid), ". Name a single ",
         "`lambda`", if (penalty == "elastic_net") " and a single `alpha`",
         ", or leave `cv = TRUE` so that the resampling can choose.",
         call. = FALSE)
  }
  grid
}


#' Turn the forest arguments into a caret tuning grid
#'
#' `mtry` is the one argument of a random forest that [caret::train()] tunes, so
#' it is the one that becomes a grid; `ntree` and `nodesize` are passed through to
#' the engine and are the same for every candidate.
#'
#' `NULL` resolves to the rule of thumb rather than to a grid, which is where this
#' departs from `sa_enet_grid()`. A penalty has no default size, so the elastic net
#' has to search for one and `cv = FALSE` is refused unless the caller names a
#' single `lambda`. `mtry` does have a default worth fitting — the square root of
#' the predictor count for a classification and a third of it for a regression —
#' so `fit_rf(data, outcome, cv = FALSE)` is a complete call.
#'
#' A value above the predictor count is refused rather than passed on.
#' `randomForest()` resets it to the valid range with a warning and fits, so what
#' comes back would be a forest at a different `mtry` from the one the result
#' records.
#'
#' @param mtry Predictors to offer at each split, one value or several, or `NULL`
#'   for the rule of thumb.
#' @param p Predictors the model has. A forest splits a factor directly rather
#'   than on dummy columns, so this is the column count and not a term count.
#' @param classify Whether the outcome is being classified, which is what decides
#'   which rule of thumb applies.
#' @param cv Whether anything is being resampled.
#'
#' @return data.frame with one row per candidate, column `mtry`.
#'
#' @keywords internal
#' @noRd
sa_rf_grid <- function(mtry, p, classify, cv) {
  if (is.null(mtry)) {
    mtry <- if (classify) floor(sqrt(p)) else floor(p / 3)
    mtry <- max(1L, mtry)
  } else {
    sa_check_num_vector(mtry, "mtry", 1)
    fractional <- unique(mtry[mtry != trunc(mtry)])
    if (length(fractional) > 0L) {
      stop("`mtry` counts predictors, so it must hold whole numbers, but holds ",
           paste(fractional, collapse = ", "), ".", call. = FALSE)
    }
    above <- unique(mtry[mtry > p])
    if (length(above) > 0L) {
      stop("`mtry` cannot exceed the ", p, " predictor(s) the model has, but ",
           "holds ", paste(above, collapse = ", "),
           ". `randomForest()` would reset it to the valid range and fit a ",
           "forest at a different `mtry` than the one reported.", call. = FALSE)
    }
  }
  grid <- expand.grid(mtry = unique(as.integer(mtry)))

  if (!cv && nrow(grid) > 1L) {
    stop("`cv = FALSE` fits one forest, so the grid must hold one candidate, ",
         "but `mtry` gives ", nrow(grid), ". Name a single `mtry`, or leave ",
         "`cv = TRUE` so that the resampling can choose.", call. = FALSE)
  }
  grid
}


#' Turn the machine's arguments into a caret tuning grid
#'
#' Both arguments of a radial-kernel machine are tuned, so both become the grid,
#' which puts this between the other two: the elastic net searches for a penalty
#' that has no default, and a forest's `mtry` has one worth fitting. `sigma` does
#' have a default, but it is read off the data by `sa_svm_sigma()` before this is
#' called, so what arrives here is always a number.
#'
#' Zero is rejected by name rather than by the bound. `sa_check_num_vector()`
#' bounds inclusively, and neither argument is answerable at zero: a machine at
#' `C = 0` pays nothing for violating its margin and fits a flat surface, and a
#' kernel at `sigma = 0` reports the same distance between every pair of rows.
#' `kernlab` takes both and returns a fit, so the refusal has to be here.
#'
#' @param C,sigma The costs and kernel widths to score.
#' @param cv Whether anything is being resampled.
#'
#' @return data.frame with one row per candidate, columns `sigma` and `C`.
#'
#' @keywords internal
#' @noRd
sa_svm_grid <- function(C, sigma, cv) {
  specs <- list(
    list(value = C, arg = "C",
         reason = paste("a machine that pays nothing for violating its margin",
                        "fits a flat surface")),
    list(value = sigma, arg = "sigma",
         reason = paste("a kernel of no width reports the same distance between",
                        "every pair of rows"))
  )
  for (spec in specs) {
    sa_check_num_vector(spec$value, spec$arg, 0)
    if (any(spec$value == 0)) {
      stop("`", spec$arg, "` must be above 0, but holds 0: ", spec$reason,
           ", and `kernlab` fits it without complaint.", call. = FALSE)
    }
  }
  grid <- expand.grid(sigma = unique(sigma), C = unique(C))

  if (!cv && nrow(grid) > 1L) {
    stop("`cv = FALSE` fits one machine, so the grid must hold one candidate, ",
         "but `C` and `sigma` give ", nrow(grid), ". Name a single `C` and a ",
         "single `sigma`, or leave `cv = TRUE` so that the resampling can ",
         "choose.", call. = FALSE)
  }
  grid
}


#' Resolve the outcome and the predictors out of one data frame
#'
#' The model functions take the same wide data frame the comparison functions
#' take, which here is the training half [split_data()] handed back. What differs
#' is that one of its columns is the outcome and the rest are candidates for
#' predictors, so the two have to be told apart before anything is fitted.
#'
#' Rows with a missing value anywhere in the model are dropped here rather than
#' inside the engine. Left to `lm()`, deletion happens once per fold on whatever
#' each fold happens to hold, so the folds would be scored on different subsets
#' of the data and the resampled numbers would not be comparable. Dropping first
#' means every fold sees the same rows, and how many went is reported.
#'
#' A predictor that takes one value cannot contribute, and leaving it in makes
#' the engine return an `NA` coefficient that reads like a failure rather than
#' like a column with nothing in it. It is dropped with a message instead, the
#' same way an axis with no usable distance is left unclustered in
#' [draw_heatmap()] rather than aborting the plot.
#'
#' @param data Wide data.frame (or matrix), one row per observation.
#' @param outcome Column name of the outcome, or a vector with one entry per row.
#' @param predictors Column names to use, or `NULL` for every column except the
#'   outcome.
#'
#' @return List with the predictor frame `x`, the outcome vector `y`, the
#'   `outcome` label, the `predictors` kept and those `dropped_predictors`, the
#'   `predictor_lv` of every predictor that has levels, and the row counts
#'   `n_obs`, `n_used`, `n_dropped`.
#'
#' @keywords internal
#' @noRd
sa_resolve_model_input <- function(data, outcome, predictors = NULL) {
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or a matrix.", call. = FALSE)
  }
  n_obs <- nrow(data)
  if (n_obs == 0L) {
    stop("`data` has zero rows.", call. = FALSE)
  }

  # NA is allowed through the resolver here: an unusable outcome is a row the
  # listwise deletion below removes and counts, not a call that cannot proceed.
  resolved <- sa_resolve_row_vector(outcome, "outcome", data, allow_na = TRUE)
  if (is.null(resolved$value)) {
    stop("`outcome` must name a column of `data` or hold one entry per row ",
         "of it.", call. = FALSE)
  }
  y <- resolved$value
  outcome_label <- resolved$label

  if (is.null(predictors)) {
    # An outcome passed as a vector leaves no column to exclude, so every column
    # is a candidate. Passed as a name, that column is the one thing a predictor
    # must not be.
    predictors <- setdiff(names(data), outcome_label)
  }
  if (!is.character(predictors) || length(predictors) == 0L ||
        anyNA(predictors)) {
    stop("`predictors` must be a non-empty character vector of column names, ",
         "or NULL for every column except `outcome`.", call. = FALSE)
  }
  dup <- unique(predictors[duplicated(predictors)])
  if (length(dup) > 0L) {
    stop("`predictors` contains duplicated names: ",
         paste(dup, collapse = ", "), ".", call. = FALSE)
  }
  unknown <- setdiff(predictors, names(data))
  if (length(unknown) > 0L) {
    stop("`predictors` not found in `data`: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (outcome_label %in% predictors) {
    stop("`predictors` contains the outcome column `", outcome_label,
         "`, which would let the model predict from the answer.", call. = FALSE)
  }

  x <- data[predictors]
  unsupported <- predictors[!vapply(
    x,
    function(v) is.numeric(v) || is.logical(v) || is.factor(v) ||
      is.character(v),
    logical(1)
  )]
  if (length(unsupported) > 0L) {
    stop("`predictors` must be numeric, logical, factor or character columns. ",
         "Not usable: ", paste(unsupported, collapse = ", "), ".",
         call. = FALSE)
  }

  keep <- stats::complete.cases(x) & !is.na(y)
  n_used <- sum(keep)
  if (n_used < 2L) {
    stop("only ", n_used, " row(s) of `data` are complete across `outcome` ",
         "and `predictors`; at least 2 are needed.", call. = FALSE)
  }
  x <- x[keep, , drop = FALSE]
  y <- y[keep]
  rownames(x) <- NULL

  # Character columns become factors here rather than inside the engine, so that
  # the levels are fixed before the folds are drawn. Unused levels left over from
  # the row filtering are dropped for the same reason a constant column is: an
  # all-zero dummy column is not a predictor.
  for (nm in predictors) {
    if (is.character(x[[nm]])) {
      x[[nm]] <- factor(x[[nm]])
    } else if (is.factor(x[[nm]])) {
      x[[nm]] <- droplevels(x[[nm]])
    }
  }

  constant <- predictors[vapply(x, function(v) length(unique(v)) < 2L,
                               logical(1))]
  if (length(constant) > 0L) {
    x <- x[setdiff(predictors, constant)]
    message("predictor(s) with a single value cannot contribute and were left ",
            "out: ", paste(constant, collapse = ", "), ".")
  }
  predictors <- setdiff(predictors, constant)
  if (length(predictors) == 0L) {
    stop("every predictor takes a single value over the usable rows, so there ",
         "is nothing to fit.", call. = FALSE)
  }

  # The levels are settled above and nowhere else, so this is where they are
  # recorded. `predict.sa_model()` has to code a factor in `newdata` the way the
  # fit coded it, and by then the engine object holds only the names of the dummy
  # columns that came out.
  is_factor <- vapply(x, is.factor, logical(1))

  list(
    x                  = x,
    y                  = y,
    outcome            = outcome_label,
    predictors         = predictors,
    dropped_predictors = constant,
    predictor_lv       = lapply(x[is_factor], levels),
    n_obs              = n_obs,
    n_used             = n_used,
    n_dropped          = n_obs - n_used
  )
}


#' The levels slot of `design`, left out when nothing has levels
#'
#' `design` reports what the model saw, and a set of models sees no factor at
#' all. An empty named list there would be a slot that says "these are the
#' levels" about nothing, which is the same reason `outcome_lv` is absent from a
#' regression rather than present and empty.
#'
#' @param predictor_lv Named list of levels, possibly empty.
#'
#' @return A one-element list to splice into `design`, or `NULL`.
#'
#' @keywords internal
#' @noRd
sa_design_lv <- function(predictor_lv) {
  if (length(predictor_lv) == 0L) {
    return(NULL)
  }
  list(predictor_lv = predictor_lv)
}


#' Dummy code the predictor frame into the matrix a penalized engine needs
#'
#' `lm()` and `glm()` are handed the predictor frame as it is, because `caret`
#' builds a formula for them and the model frame does the dummy coding. `glmnet`
#' takes a numeric matrix, and what `caret` does to a data.frame before passing it
#' on is `Matrix::as.matrix()`, which turns a factor into its integer codes. The
#' result fits without complaint and is wrong: a three-level factor arrives as one
#' evenly spaced numeric predictor rather than as two dummies, so the model
#' assumes an order and a spacing between the levels that nobody stated.
#'
#' The coding is [stats::model.matrix()]'s, which is the same one `lm()` would
#' have applied, so a `k`-level factor becomes the same `k - 1` terms under the
#' same names and the two models' coefficient tables can be read side by side.
#' The intercept column is dropped because `glmnet` fits an unpenalized intercept
#' of its own.
#'
#' Predicting on new rows codes them here as well, which is what `xlev` and
#' `want` are for, since a matrix built by any other route is a matrix that could
#' disagree with the one the model was fitted to. `xlev` fixes the levels, so a
#' level that no row of `newdata` happens to take still gets its column of
#' zeroes, and `want` fixes the order by name, since a design matrix is read by
#' position once it reaches `glmnet`. Rows are kept whatever they hold: a missing
#' cell has to reach the caller as a missing prediction rather than as a row that
#' quietly went absent, so the model frame passes `NA` through instead of
#' dropping it.
#'
#' @param x Predictor frame from `sa_resolve_model_input()`, already reduced to
#'   the usable rows and with its character columns turned into factors, or the
#'   same columns of the rows to predict.
#' @param xlev Named list of factor levels to code against, as
#'   `design$predictor_lv` holds them, or `NULL` to read the levels off `x`.
#' @param want Column names the result must have, in the order it must have
#'   them, or `NULL` to keep what the coding produced.
#'
#' @return Numeric matrix with one column per model term, no intercept.
#'
#' @keywords internal
#' @noRd
sa_design_matrix <- function(x, xlev = NULL, want = NULL) {
  mf <- stats::model.frame(~ ., data = x, xlev = xlev,
                           na.action = stats::na.pass)
  mm <- stats::model.matrix(~ ., data = mf)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

  if (!is.null(want)) {
    absent <- setdiff(want, colnames(mm))
    if (length(absent) > 0L) {
      stop("internal error: the coding of `newdata` is missing term(s) the ",
           "model has: ", paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    mm <- mm[, want, drop = FALSE]
  }
  rownames(mm) <- NULL
  mm
}


#' Reduce new rows to the predictors the model was fitted on
#'
#' What a model needs of `newdata` is its own predictors, coded the way they were
#' coded when it was fitted, and nothing else. Columns it never saw are ignored
#' rather than refused, since the rows to predict usually arrive as the other half
#' of the same data frame the fit was given, outcome column and all. A predictor
#' that is not there at all is an error naming it: the alternative is a model
#' predicting from fewer columns than it has coefficients for, which is either a
#' failure further in or a silently wrong number.
#'
#' The levels are put back rather than read off `newdata`, and that is the whole
#' reason this exists. A factor in a held-out half may be missing a level, or
#' carry its levels in another order, and either would code to a different matrix
#' from the one the model was fitted to. A level the fit never saw is an error
#' instead, since there is no coefficient to apply to it.
#'
#' @param newdata Rows to predict, a data.frame or matrix.
#' @param design The `design` element of the model, which names the predictors
#'   and holds `predictor_lv`.
#'
#' @return data.frame of the predictor columns, in the model's own order, with
#'   every factor at the levels the fit used.
#'
#' @keywords internal
#' @noRd
sa_predict_frame <- function(newdata, design) {
  if (is.matrix(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data.frame or a matrix.", call. = FALSE)
  }
  if (nrow(newdata) == 0L) {
    stop("`newdata` has zero rows.", call. = FALSE)
  }

  predictors <- design$predictors
  absent <- setdiff(predictors, names(newdata))
  if (length(absent) > 0L) {
    stop("`newdata` is missing predictor column(s) the model was fitted on: ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }

  x <- newdata[predictors]
  for (nm in predictors) {
    v <- x[[nm]]
    lv <- design$predictor_lv[[nm]]

    if (is.null(lv)) {
      if (!is.numeric(v) && !is.logical(v)) {
        stop("`", nm, "` was a numeric predictor when the model was fitted, ",
             "and `newdata` holds it as ", class(v)[1L], ". The coding of a ",
             "column cannot change between fitting and predicting.",
             call. = FALSE)
      }
      next
    }

    if (is.factor(v) || is.logical(v) || is.numeric(v)) {
      v <- as.character(v)
    }
    if (!is.character(v)) {
      stop("`", nm, "` was a factor when the model was fitted, and `newdata` ",
           "holds it as ", class(v)[1L], ", which has no levels to match.",
           call. = FALSE)
    }
    unseen <- setdiff(unique(v[!is.na(v)]), lv)
    if (length(unseen) > 0L) {
      stop("`newdata` holds level(s) of `", nm, "` the model was not fitted ",
           "on, so there is no coefficient for them: ",
           paste(sort(unseen), collapse = ", "), ". Fitted on: ",
           paste(lv, collapse = ", "), ".", call. = FALSE)
    }
    x[[nm]] <- factor(v, levels = lv)
  }
  rownames(x) <- NULL
  x
}


#' Put a binary outcome in the order that fixes the direction
#'
#' The first level is the reference at every group count in this package, so it
#' is the level a fold change divides by and the one a post-hoc contrast
#' subtracts. A logistic regression follows the same rule: with
#' `outcome_lv = c("control", "case")` every coefficient is the change in the log
#' odds of `case`, and its odds ratio is above 1 for a predictor that raises the
#' chance of `case`.
#'
#' That happens to be what the engine does on its own. [stats::glm()] models the
#' probability of the *last* level of a factor, so ordering the levels
#' reference-first is enough and nothing has to be reversed afterwards. The draft
#' this replaces relevelled the class of interest to the front, which is the
#' order [caret::train()] wants for its class-specific metrics but the opposite of
#' what `glm()` reads, so the reported odds belonged to the other class.
#'
#' A third level is an error rather than a dropped set of rows. Two levels are
#' what this model is, so silently fitting a different subset of the data than
#' was passed in would answer a question nobody asked.
#'
#' `control_label` names the same level `outcome_lv[1]` names, and exists because
#' most calls have nothing to say about the other one. Sorting already puts
#' `"control"` before `"treated"` and `0` before `1`, so what a caller usually
#' wants is to correct the one case where the sort is wrong, and writing both
#' levels out to move one of them is more than the correction is worth.
#'
#' Naming the reference twice and disagreeing is an error rather than a
#' precedence rule. Either of the two arguments is a complete answer on its own,
#' so a call that carries both and contradicts itself has no reading that is more
#' likely than the other.
#'
#' @param y Outcome vector, already reduced to the usable rows.
#' @param outcome_lv The two levels, reference first, or `NULL` to sort them.
#' @param control_label The reference level on its own, or `NULL` to leave the
#'   order as `outcome_lv` or the sort gave it.
#' @param model What to call the model in the message a third class raises, since
#'   more than one model reads an outcome this way.
#'
#' @return `y` as a factor whose levels are `outcome_lv`.
#'
#' @keywords internal
#' @noRd
sa_outcome_levels <- function(y,
                              outcome_lv = NULL,
                              control_label = NULL,
                              model = "a logistic regression") {
  named <- !is.null(outcome_lv)
  if (is.factor(y)) {
    y <- as.character(y)
  }
  if (is.logical(y) || is.numeric(y)) {
    y <- as.character(y)
  }
  if (!is.character(y)) {
    stop("`outcome` must be a factor, character, logical or numeric column ",
         "holding two classes.", call. = FALSE)
  }
  present <- unique(y)
  if (length(present) < 2L) {
    stop("`outcome` takes a single value over the usable rows, so there is ",
         "nothing to classify.", call. = FALSE)
  }

  if (is.null(outcome_lv)) {
    if (length(present) > 2L) {
      stop("`outcome` holds ", length(present), " classes, but ", model,
           " models two: ", paste(sort(present), collapse = ", "),
           ". Name the two to model with `outcome_lv`, or reduce `data` to ",
           "them first.", call. = FALSE)
    }
    outcome_lv <- sort(present)
  } else {
    outcome_lv <- as.character(outcome_lv)
    if (length(outcome_lv) != 2L || anyNA(outcome_lv) ||
          anyDuplicated(outcome_lv) > 0L) {
      stop("`outcome_lv` must be two distinct level names, the reference ",
           "first.", call. = FALSE)
    }
    absent <- setdiff(outcome_lv, present)
    if (length(absent) > 0L) {
      stop("`outcome_lv` level(s) absent from `outcome`: ",
           paste(absent, collapse = ", "), ". Present: ",
           paste(sort(present), collapse = ", "), ".", call. = FALSE)
    }
  }

  # A named pair that leaves classes out would fit the model on a subset of the
  # rows that were passed in, which is a different data set from the one the call
  # describes. Reducing `data` first says so in the caller's own code.
  extra <- setdiff(present, outcome_lv)
  if (length(extra) > 0L) {
    stop("`outcome` holds ", length(present), " classes and `outcome_lv` names ",
         "two of them, so ", length(extra), " would be silently left out: ",
         paste(sort(extra), collapse = ", "),
         ". Reduce `data` to the two classes first.", call. = FALSE)
  }

  if (!is.null(control_label)) {
    if (length(control_label) != 1L || is.na(control_label)) {
      stop("`control_label` must be a single level name, the one to hold as ",
           "the reference.", call. = FALSE)
    }
    control_label <- as.character(control_label)
    if (!control_label %in% outcome_lv) {
      stop("`control_label` names a class `outcome` does not hold: ",
           control_label, ". Present: ",
           paste(sort(present), collapse = ", "), ".", call. = FALSE)
    }
    if (named && !identical(control_label, outcome_lv[1])) {
      stop("`control_label` names ", control_label, " as the reference and ",
           "`outcome_lv` puts ", outcome_lv[1], " first, so the two disagree ",
           "about which class the other one is compared against. Pass one of ",
           "them.", call. = FALSE)
    }
    outcome_lv <- c(control_label, setdiff(outcome_lv, control_label))
  }

  factor(y, levels = outcome_lv)
}


#' The confidence interval that agrees with the standard error beside it
#'
#' Built from the summary table rather than taken from [stats::confint()], for two
#' reasons. `confint()` on a rank-deficient fit indexes the coefficients against a
#' shorter variance matrix and fails, while a term that could not be estimated
#' should come back as an `NA` row like any other unanswerable question. And on a
#' `glm` the default `confint()` is a profile likelihood interval, which is a
#' better interval but a different quantity from the Wald standard error and z
#' value in the same row, so the three numbers would not agree.
#'
#' @param coef_matrix `summary()$coefficients`, estimates in column 1 and
#'   standard errors in column 2.
#' @param conf_level Two-sided confidence level.
#' @param df Residual degrees of freedom for a t interval, or `NULL` for the
#'   normal approximation.
#'
#' @return Two-column matrix of limits, rows named as `coef_matrix` is.
#'
#' @keywords internal
#' @noRd
sa_wald_interval <- function(coef_matrix, conf_level, df = NULL) {
  estimate <- coef_matrix[, 1L]
  stderr <- coef_matrix[, 2L]
  upper_tail <- (1 - conf_level) / 2
  crit <- if (is.null(df) || !is.finite(df)) {
    stats::qnorm(1 - upper_tail)
  } else {
    stats::qt(1 - upper_tail, df)
  }

  out <- cbind(estimate - crit * stderr, estimate + crit * stderr)
  dimnames(out) <- list(rownames(coef_matrix), c("lower_conf", "upper_conf"))
  out
}


#' Run the engine, and report its warnings once rather than per fold
#'
#' Cross-validation fits the same model many times, so a condition of the data
#' such as a perfectly separated logistic regression is raised once per fold. The
#' warnings are collected and re-emitted as one grouped message with a count, the
#' same way `sa_feature_table()` handles a warning raised once per feature. They
#' are not discarded: a model that did not converge has to say so.
#'
#' @param expr The fitting call, evaluated here.
#' @param label Human readable model name used in the message.
#'
#' @return Whatever `expr` evaluated to.
#'
#' @keywords internal
#' @noRd
sa_quiet_engine <- function(expr, label) {
  caught <- character(0)
  fit <- withCallingHandlers(
    expr,
    warning = function(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(caught) > 0L) {
    grouped <- table(caught)
    message(label, ": engine note(s) while fitting:\n",
            paste0("  [", grouped, " time(s)] ", names(grouped),
                   collapse = "\n"))
  }
  fit
}


#' Run the engine and class the fit so that the usual questions reach it
#'
#' The `sa_fit` class goes on in front of the engine's own, which is what lets
#' `coef()` and `summary()` reach the fitted model inside a [caret::train()]
#' object; see `coef.sa_fit()`. It is prepended rather than substituted, so every
#' method the engine already carries — `predict()`, `fitted()`, `residuals()`,
#' `print()` — is still found by inheritance. Note that the class names the fit
#' this returns, while the function names the act of producing it.
#'
#' Only a fit gets the class. `perform_rfe()` groups its warnings the same way
#' but leaves the object it stores alone, since an elimination has no
#' `$finalModel` for the two methods to reach and `coef()` on it would answer with
#' a message about a fit that is not there.
#'
#' @param expr The fitting call, evaluated here.
#' @param label Human readable model name used in the message.
#'
#' @keywords internal
#' @noRd
sa_fit_engine <- function(expr, label) {
  fit <- sa_quiet_engine(expr, label)
  structure(fit, class = c("sa_fit", class(fit)))
}


#' Assemble the coefficient table of a fitted linear or generalised linear model
#'
#' The row set comes from `coef()` rather than from the summary, because a term
#' the engine could not estimate is absent from `summary()$coefficients` while
#' `coef()` keeps it as `NA`. Dropping it would make the table quietly shorter
#' than the model it describes; kept as an `NA` row it says that the term was in
#' the model and could not be estimated, which is the same distinction the
#' feature tables make.
#'
#' The interval is supplied by the caller rather than computed here, since the
#' two models do not refer their statistic to the same distribution: `lm()` reads
#' a t on its residual degrees of freedom and `glm()` a Wald z.
#'
#' @param model A fitted `lm` or `glm`.
#' @param interval Two-column matrix of confidence limits, rows named by term.
#' @param df Residual degrees of freedom the statistic was referred to, or `NA`
#'   for a statistic that has none.
#'
#' @return data.frame with the columns of `sa_model_coef_columns()` and those of
#'   `sa_model_inference_columns()`, since these are the models that have both.
#'
#' @keywords internal
#' @noRd
sa_coef_table <- function(model, interval, df) {
  estimate <- stats::coef(model)
  terms <- names(estimate)
  summ <- summary(model)$coefficients
  at <- match(terms, rownames(summ))
  limits <- interval[match(terms, rownames(interval)), , drop = FALSE]

  out <- data.frame(
    terms      = terms,
    estimate   = as.numeric(estimate),
    stderr     = as.numeric(summ[at, 2L]),
    statistic  = as.numeric(summ[at, 3L]),
    df         = rep(as.numeric(df), length(terms)),
    pval       = as.numeric(summ[at, 4L]),
    lower_conf = as.numeric(limits[, 1L]),
    upper_conf = as.numeric(limits[, 2L]),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}


#' Reduce a caret result table to the rows and columns worth keeping
#'
#' `caret::train()` reports one row per hyperparameter combination it scored, and
#' with `method = "none"` it scored none, so the table comes back with its
#' columns and no rows. An empty table reads as a result that is missing
#' something rather than as one for which the question does not arise, which is
#' the same reason `posthoc` is absent from a two-group comparison instead of
#' present and empty.
#'
#' @keywords internal
#' @noRd
sa_model_frame <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(NULL)
  }
  rownames(df) <- NULL
  df
}
