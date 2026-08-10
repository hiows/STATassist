# The result contract of the feature selection family, and the fourth row axis in
# this package. A comparison repeats `features`, a model repeats `terms` and a
# reduction repeats `points`. What a selection repeats is `candidates`: the
# predictors that were offered, in the order the search ranked them rather than
# the order they arrived, so the first row is the one the procedure would keep
# last.
#
# A selection is not a model and is not stored as one. A model answers what each
# predictor is worth given all the others, and every table it holds is about the
# fit that is there. A selection answers a question one level up — how many
# predictors are worth keeping, and which — so it holds two tables that a model
# has no counterpart for: `ranking`, one row per candidate, and `profile`, one row
# per model the search compared. The size is the answer, and a table with one row
# per candidate model is the only place it can be read.
#
# What `profile` repeats therefore depends on how the search moved. `perform_rfe()`
# scores one model per subset size, so its rows are sizes and each appears once;
# `perform_stepwise()` walks a path one term at a time, so its rows are steps and
# `direction = "both"` can visit a size twice. What the two have in common is the
# `n_vars` column and the single `chosen` row, which is all a reader needs to say
# how large the answer is and what it was compared against.
#
# `resampling` is what tells the two apart at a glance. An elimination has no score
# without holding rows out, so it always resamples; a criterion is a penalised
# likelihood computed on the rows the model was fitted to, so it never does, and
# the slot is `NULL`.
#
# `$fit` is the engine object, the same exception `sa_model` and `sa_reduction`
# make and for the same reason: everything else here is a scalar, a character
# vector, a named list or a data.frame, so dropping that one slot leaves an object
# that writes out as JSON.

#' The selections this contract covers
#'
#' `analysis` names the procedure rather than the family, the way
#' `sa_reduction_analyses()` names the three reductions, so that a result says
#' which search produced it without anything having to read `engine`.
#'
#' @keywords internal
#' @noRd
sa_selection_analyses <- function() {
  c("rfe", "stepwise")
}


#' Assemble a feature selection result object
#'
#' The checks here guard the contract rather than the user's input, so they fire
#' only on a mistake inside the package and say so. What they guard is the promise
#' the object makes: `ranking` is aligned with `candidates` by position, `selected`
#' is a subset of them, and exactly one row of `profile` is the size that was
#' chosen, so the three tables cannot describe different searches.
#'
#' @param analysis Which search this is, `"rfe"` or `"stepwise"`.
#' @param candidates Candidate predictor names, the row order `ranking` follows.
#' @param design Named list describing the data the search saw.
#' @param parameters Named list of the choices as they were used.
#' @param selected The candidates the search kept.
#' @param ranking data.frame of one row per candidate.
#' @param profile data.frame of one row per model the search compared: one subset
#'   size for an elimination, one step of the path for a stepwise search.
#' @param resampling data.frame of one row per resample at the chosen size, or
#'   `NULL` for a search that held nothing out.
#' @param engine Named list naming what ran the search.
#' @param fit The engine object.
#'
#' @keywords internal
#' @noRd
sa_new_selection <- function(analysis,
                             candidates,
                             design,
                             parameters,
                             selected,
                             ranking,
                             profile,
                             resampling = NULL,
                             engine,
                             fit) {

  if (!analysis %in% sa_selection_analyses()) {
    stop("internal error: `analysis` must be one of ",
         paste(sa_selection_analyses(), collapse = ", "), ".", call. = FALSE)
  }
  if (!is.character(candidates) || length(candidates) == 0L) {
    stop("internal error: `candidates` must be a non-empty character vector.",
         call. = FALSE)
  }
  if (!is.data.frame(ranking) || !identical(ranking$candidates, candidates)) {
    stop("internal error: `ranking` is not a data.frame aligned with ",
         "`candidates`.", call. = FALSE)
  }
  # A selection that kept nothing is not a selection, and one that kept something
  # it was never offered is a table that cannot be read against `ranking`.
  if (!is.character(selected) || length(selected) == 0L) {
    stop("internal error: `selected` must be a non-empty character vector.",
         call. = FALSE)
  }
  unknown <- setdiff(selected, candidates)
  if (length(unknown) > 0L) {
    stop("internal error: `selected` holds name(s) that are not candidates: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (!identical(ranking$selected, candidates %in% selected)) {
    stop("internal error: `ranking$selected` disagrees with `selected`.",
         call. = FALSE)
  }
  if (!is.data.frame(profile) || nrow(profile) == 0L ||
        is.null(profile$n_vars)) {
    stop("internal error: `profile` must be a non-empty data.frame with an ",
         "`n_vars` column.", call. = FALSE)
  }
  if (sum(profile$chosen) != 1L) {
    stop("internal error: exactly one row of `profile` is the size that was ",
         "chosen, but ", sum(profile$chosen), " are marked.", call. = FALSE)
  }
  if (profile$n_vars[profile$chosen] != length(selected)) {
    stop("internal error: the chosen row of `profile` is a subset of ",
         profile$n_vars[profile$chosen], " variable(s) and `selected` holds ",
         length(selected), ".", call. = FALSE)
  }
  # `importance` is required here and not in `sa_model`'s engine, because a
  # ranking is only readable once it says what it ranked by.
  for (nm in c("package", "method", "label", "metrics", "importance")) {
    if (is.null(engine[[nm]])) {
      stop("internal error: `engine` is missing `", nm, "`.", call. = FALSE)
    }
  }
  if (!is.null(resampling) && !is.data.frame(resampling)) {
    stop("internal error: `resampling` must be a data.frame or NULL.",
         call. = FALSE)
  }

  structure(
    list(
      analysis   = analysis,
      candidates = candidates,
      design     = design,
      parameters = parameters,
      selected   = selected,
      ranking    = ranking,
      profile    = profile,
      resampling = resampling,
      engine     = engine,
      fit        = fit,
      metadata   = sa_metadata()
    ),
    class = c("sa_selection", "sa_result")
  )
}


#' Print a feature selection
#'
#' Summarises what was searched and what was kept, rather than printing every
#' table. The score of every model the search compared is in `x$profile`, which is
#' what says whether the answer won by much or by nothing, and the engine object is
#' `x$fit`.
#'
#' @param x A selection, as returned by [perform_rfe()] or [perform_stepwise()].
#' @param n Maximum number of candidates to show. The rest are counted.
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' perform_rfe(mtcars, outcome = "mpg", predictors = c("wt", "hp", "disp"),
#'             cv_method = "kfold", n_fold = 3, seed = 1)
#'
#' ## The same method, on a search that walked a path rather than a ladder of
#' ## subset sizes and that held no rows out to score it.
#' perform_stepwise(mtcars, outcome = "mpg", predictors = c("wt", "hp", "disp"))
#'
#' @export
print.sa_selection <- function(x, n = 10L, ...) {
  n <- sa_check_count(n, "n", 0)
  design <- x$design
  params <- x$parameters
  # `metric` for a search that resamples and `criterion` for one that charges a
  # penalty. They are the same field under the name each search calls its own
  # number by, and either one names the column of `profile` the answer was chosen
  # by, which is all this function needs of it.
  chose_by <- if (is.null(params$metric)) params$criterion else params$metric

  cat("<sa_selection> ", x$analysis, "\n", sep = "")
  cat("  outcome  : ", design$outcome, "  (", design$outcome_type, ")\n",
      sep = "")
  if (!is.null(design$outcome_lv)) {
    cat("             modelling ", design$outcome_lv[2], " against ",
        design$outcome_lv[1], ", ", design$n_events, " of ", design$n_used,
        " row(s)\n", sep = "")
  }
  cat("  rows     : ", design$n_used, " used",
      if (design$n_dropped > 0L) {
        paste0("  (", design$n_dropped, " incomplete row(s) dropped)")
      },
      "\n", sep = "")
  # What was compared is read off `profile` rather than out of `parameters`, since
  # the table of the models the search visited is the record of what was asked for.
  # An elimination scores one model per subset size and a stepwise search takes one
  # move at a time, so the same table is counted two ways.
  sa_cat_field("search", paste0(
    x$engine$label, " over ", length(x$candidates), " candidate(s), ",
    switch(
      x$analysis,
      rfe = paste0("size(s) ", paste(x$profile$n_vars, collapse = ", ")),
      stepwise = paste0(nrow(x$profile) - 1L, " step(s)"),
      paste0(nrow(x$profile), " model(s) compared")
    )
  ))
  if (identical(x$analysis, "stepwise")) {
    cat("  settings : ", params$direction, " search, ", chose_by, " ",
        if (params$maximize) "maximised" else "minimised",
        " at ", sa_fmt_num(params$k, 3), " per parameter\n", sep = "")
  } else {
    cat("  settings : ", params$cv_method,
        if (!is.na(params$n_fold)) paste0(", ", params$n_fold, " fold(s)"),
        if (!is.na(params$n_repeat)) paste0(" x ", params$n_repeat,
                                            " repeat(s)"),
        ", ", chose_by, " ",
        if (params$maximize) "maximised" else "minimised",
        "\n", sep = "")
  }

  # The score of the model that was kept is the one number that says whether the
  # search answered anything, so it is printed on the line that reports the answer.
  # Three significant digits show all of a resampled metric, which is a small
  # number. A criterion is on the scale of the row count while the differences that
  # decided its path are of order one, so it needs the digits to show them.
  best <- x$profile[x$profile$chosen, , drop = FALSE]
  digits <- if (identical(x$analysis, "stepwise")) 6L else 3L
  sd_col <- paste0(chose_by, "SD")
  cat("  selected : ", length(x$selected), " of ", length(x$candidates),
      "  (", chose_by, " = ", sa_fmt_num(best[[chose_by]], digits),
      if (!is.null(best[[sd_col]])) {
        paste0(" (SD ", sa_fmt_num(best[[sd_col]], 2), ")")
      },
      if (!is.null(x$resampling)) {
        paste0(" over ", nrow(x$resampling), " resample(s)")
      },
      ")\n", sep = "")

  cat("\n  ranking  (", x$engine$importance, ")\n", sep = "")
  shown <- utils::head(x$ranking, n)
  width <- if (nrow(shown) > 0L) max(nchar(shown$candidates)) else 0L
  for (i in seq_len(nrow(shown))) {
    row <- shown[i, ]
    cat("    ", formatC(row$candidates, width = -width), "  ",
        formatC(sa_fmt_num(row$estimate), width = 10),
        if (isTRUE(row$selected)) "  selected" else "  dropped",
        "\n", sep = "")
  }
  if (nrow(x$ranking) > nrow(shown)) {
    cat("    ... and ", nrow(x$ranking) - nrow(shown),
        " more candidate(s) in $ranking\n", sep = "")
  }

  if (length(design$dropped_predictors) > 0L) {
    cat("\n  dropped  : ", paste(design$dropped_predictors, collapse = ", "),
        " (single valued)\n", sep = "")
  }

  invisible(x)
}
