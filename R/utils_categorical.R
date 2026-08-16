# Resolving categorical input, laying it out as a contingency table, and saying
# what each cell of that table was expected to hold.
#
# `sa_validate_wide_input()` cannot serve here: it requires every named column to
# be numeric, which is exactly the requirement this scenario inverts. What it can
# share is everything after the columns are read -- the level ordering rules, the
# `control_label` reordering, the row dropping and the messages that report it --
# so those come from `sa_control_first()` and `sa_fact_control_first()` rather
# than being written a third time.
#
# The expected counts are built by one function per null hypothesis rather than
# by one function with a branch, because "expected" is not a property of a table.
# It is a property of a table and a claim about it, and the three claims this
# scenario tests make three different numbers out of the same counts.


#' Resolve the categorical variables and the levels of each
#'
#' Everything the analysis needs to know about where an observation sits, settled
#' in one pass, so that the levels, the table and the row selection cannot be
#' derived from each other twice in different orders. This is the counterpart of
#' `sa_fact_layout()` for a design whose variables are the answers rather than
#' the strata.
#'
#' Two kinds of unusable row are counted apart rather than together. A row naming
#' a level `category_lv` leaves out was measured and excluded; a row missing a
#' value was not measured. Both leave the table, and `match()` would make them
#' indistinguishable, so they are found before it is used.
#'
#' A matched design is the one place the level sets are not resolved per
#' variable. Repeated measurements of one thing share a level set by definition,
#' and a square table is what McNemar's test needs, so the levels are unified
#' across the conditions and `control_label` is a single name rather than one per
#' variable.
#'
#' @param data Wide data.frame or matrix, one row per observation.
#' @param category_lv Named list of levels per variable, or `NULL` to take every
#'   column of `data` with its levels in sorted order.
#' @param control_label The level to hold first, one name per variable for an
#'   independent design and a single name for a matched one, or `NULL`.
#' @param paired Whether the columns are repeated conditions on the same rows.
#' @param max_levels How many levels a variable may take before it is refused as
#'   a category. See `sa_check_level_count()`.
#'
#' @return List with the row-filtered `data` holding only the used columns as
#'   factors, `variables` in order, `category_lv`, `n_used`, `n_dropped` and
#'   `n_incomplete`.
#'
#' @keywords internal
#' @noRd
sa_validate_categorical_input <- function(data,
                                         category_lv,
                                         control_label,
                                         paired,
                                         max_levels = 20L) {
  if (is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or a matrix.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("`data` has zero rows.", call. = FALSE)
  }

  named_lv <- !is.null(category_lv)
  if (named_lv) {
    if (!is.list(category_lv) || length(category_lv) == 0L ||
          is.null(names(category_lv)) || anyNA(names(category_lv)) ||
          !all(nzchar(names(category_lv))) ||
          anyDuplicated(names(category_lv)) > 0L) {
      stop("`category_lv` must be a named list, one entry per categorical ",
           "variable, holding that variable's levels: ",
           "list(cat_1 = c(\"y\", \"n\"), cat_2 = c(\"high\", \"low\")).",
           call. = FALSE)
    }
    variables <- names(category_lv)
    unknown <- setdiff(variables, names(data))
    if (length(unknown) > 0L) {
      stop("`category_lv` names column(s) absent from `data`: ",
           paste(unknown, collapse = ", "), ". Present: ",
           paste(names(data), collapse = ", "), ".", call. = FALSE)
    }
  } else {
    variables <- names(data)
  }

  if (length(variables) < 2L) {
    stop("a categorical comparison crosses at least two variables, and ",
         if (named_lv) "`category_lv` names " else "`data` holds ",
         length(variables), ".", call. = FALSE)
  }

  not_atomic <- variables[!vapply(data[variables], is.atomic, logical(1))]
  if (length(not_atomic) > 0L) {
    stop("a categorical variable must be an atomic column. Not atomic: ",
         paste(not_atomic, collapse = ", "), ".", call. = FALSE)
  }
  # A factor arrives as its labels and a 0/1 coded variable as "0" and "1", so
  # both are categorical here and neither keeps a storage mode the table would
  # have to remember.
  values <- lapply(variables, function(nm) as.character(data[[nm]]))
  names(values) <- variables

  observed <- lapply(values, function(v) sort(unique(v[!is.na(v)])))
  names(observed) <- variables
  empty <- variables[lengths(observed) == 0L]
  if (length(empty) > 0L) {
    stop("variable(s) holding no non-missing value: ",
         paste(empty, collapse = ", "), ".", call. = FALSE)
  }

  if (named_lv) {
    for (nm in variables) {
      lv <- category_lv[[nm]]
      if (!is.character(lv) && !is.factor(lv) && !is.atomic(lv)) {
        stop("`category_lv$", nm, "` must hold that variable's levels.",
             call. = FALSE)
      }
      lv <- as.character(lv)
      if (length(lv) < 2L || anyNA(lv) || anyDuplicated(lv) > 0L) {
        stop("`category_lv$", nm, "` must hold at least two distinct ",
             "non-missing levels.", call. = FALSE)
      }
      absent <- lv[!lv %in% observed[[nm]]]
      if (length(absent) > 0L) {
        stop("`category_lv$", nm, "` level(s) absent from `data$", nm, "`: ",
             paste(absent, collapse = ", "), ".", call. = FALSE)
      }
      category_lv[[nm]] <- lv
    }
  } else {
    category_lv <- observed
  }

  # Checked on the resolved levels rather than on the column, so naming three
  # levels of a variable that happens to take fifty is a way through rather than
  # a second thing to argue with.
  sa_check_level_count(category_lv, max_levels, named_lv)

  # The levels are settled here and nowhere later, so the table, the cell labels
  # and the row selection are all built from the order the reference ended up in
  # rather than corrected afterwards.
  if (paired) {
    category_lv <- sa_categorical_shared_lv(category_lv, observed, named_lv,
                                            control_label)
  } else {
    category_lv <- sa_fact_control_first(category_lv, control_label,
                                         if (named_lv) "category_lv" else "data")
  }

  incomplete <- Reduce(`|`, lapply(values, is.na))
  # `%in%` reads NA as absent, so the missing rows are masked out first and the
  # two counts stay separate.
  outside <- !incomplete &
    Reduce(`|`, lapply(variables, function(nm) {
      !values[[nm]] %in% category_lv[[nm]]
    }))
  keep <- !incomplete & !outside

  if (sum(keep) < 2L) {
    stop("only ", sum(keep), " row(s) hold a level of every variable, which is ",
         "not a table. Dropped: ", sum(incomplete), " for a missing value and ",
         sum(outside), " for a level outside `category_lv`.", call. = FALSE)
  }

  out <- data[keep, variables, drop = FALSE]
  for (nm in variables) {
    out[[nm]] <- factor(values[[nm]][keep], levels = category_lv[[nm]])
  }
  rownames(out) <- NULL

  list(
    data         = out,
    variables    = variables,
    category_lv  = category_lv,
    n_used       = sum(keep),
    n_dropped    = sum(outside),
    n_incomplete = sum(incomplete)
  )
}


#' Refuse a variable that is a measurement rather than a category
#'
#' `as.character()` turns any column into labels, which is what lets a factor, a
#' logical and a 0/1 code all read as categorical here. It also turns a
#' continuous measurement into one label per observation, and a table with a cell
#' per observation is not a thing a test of association has anything to say
#' about: every expected count is a fraction, Fisher's enumeration does not
#' finish, and the answer that comes back is arithmetic on noise.
#'
#' So there is a ceiling, and it is an argument rather than a constant, because a
#' genuinely many-levelled category exists and the caller is the one who knows
#' whether theirs is one.
#'
#' @keywords internal
#' @noRd
sa_check_level_count <- function(category_lv, max_levels, named_lv) {
  max_levels <- sa_check_count(max_levels, "max_levels", 2)

  over <- names(category_lv)[lengths(category_lv) > max_levels]
  if (length(over) == 0L) {
    return(invisible(NULL))
  }

  counted <- paste0(over, " (", lengths(category_lv)[over], ")",
                    collapse = ", ")
  stop("variable(s) taking more levels than the ", max_levels,
       " `max_levels` allows: ", counted,
       ". A measurement read as a category makes a table with a cell per ",
       "observation, which no test of association is about. Bin the variable ",
       "first, or name the levels to keep in `category_lv`",
       if (!named_lv) " (which also drops the rows at the rest)" else "",
       ". Raise `max_levels` if the variable really has this many.",
       call. = FALSE)
}


#' Unify the level sets of a matched design
#'
#' Repeated measurements of one thing take one set of levels, so a matched design
#' has one to resolve rather than one per condition. Taking the union when the
#' levels come from the data is what keeps a condition in which nobody answered
#' `"n"` from silently producing a table that is not square.
#'
#' @keywords internal
#' @noRd
sa_categorical_shared_lv <- function(category_lv, observed, named_lv,
                                     control_label) {
  variables <- names(category_lv)

  if (named_lv) {
    first <- category_lv[[1]]
    disagree <- variables[!vapply(category_lv, function(lv) {
      setequal(lv, first)
    }, logical(1))]
    if (length(disagree) > 0L) {
      stop("`paired = TRUE` reads the columns as repeated measurements of one ",
           "thing, so every variable takes the same levels and the table is ",
           "square. `category_lv` disagrees at: ",
           paste(disagree, collapse = ", "), ".", call. = FALSE)
    }
    shared <- first
  } else {
    shared <- sort(unique(unlist(observed, use.names = FALSE)))
  }

  if (!is.null(control_label)) {
    if (is.list(control_label)) {
      control_label <- unlist(control_label, use.names = TRUE)
    }
    if (length(control_label) != 1L) {
      stop("a matched design has one level set shared by every condition, so ",
           "`control_label` is a single level name rather than one per ",
           "variable. Got ", length(control_label), ".", call. = FALSE)
    }
    shared <- sa_control_first(shared, unname(control_label),
                               lv_arg = if (named_lv) "category_lv" else "data")
  }

  stats::setNames(rep(list(shared), length(variables)), variables)
}


#' Cross two variables into a contingency table
#'
#' @param data Validated data whose columns are factors at their levels.
#' @param variables The two variable names, row first.
#'
#' @keywords internal
#' @noRd
sa_categorical_counts <- function(data, variables) {
  out <- table(data[[variables[1]]], data[[variables[2]]])
  dimnames(out) <- stats::setNames(dimnames(out), variables)
  out
}


#' Summarise repeated binary conditions as a condition-by-response table
#'
#' Cochran's Q is asked of a subjects-by-conditions matrix, which has no
#' two-variable cross-classification to tabulate. What it is asked *about* is
#' whether the conditions share a marginal response rate, and that question has a
#' table: one row per condition and one column per level of the shared response.
#' It is the table to plot for the same reason -- it holds the rates the test
#' compares -- and it is not the paired table McNemar's test is read from, which
#' crosses two conditions against each other and only exists when there are two.
#'
#' @keywords internal
#' @noRd
sa_categorical_condition_counts <- function(data, variables, levels) {
  out <- t(vapply(variables, function(nm) {
    as.integer(table(data[[nm]]))
  }, integer(length(levels))))
  dimnames(out) <- list(condition = variables, response = levels)
  as.table(out)
}


#' The counts a table was expected to hold if the two variables were independent
#'
#' The product of the margins over the total. This is what the chi-square test of
#' independence and Fisher's exact test are both about, and it is also what the
#' condition-by-response table of a repeated design is held against: there the
#' row margins are fixed at the number of subjects, so the same arithmetic says
#' that every condition shows the pooled response rate, which is marginal
#' homogeneity rather than independence. The formula does not distinguish the two
#' claims; `design$null` is what records which one was made.
#'
#' @param counts Two-dimensional table or matrix of counts.
#'
#' @keywords internal
#' @noRd
sa_expected_independence <- function(counts) {
  out <- outer(rowSums(counts), colSums(counts)) / sum(counts)
  dimnames(out) <- dimnames(counts)
  out
}


#' The counts a square table was expected to hold if it were symmetric
#'
#' The average of each cell and its transpose. Two things follow, and both are
#' the whole content of a matched comparison.
#'
#' The diagonal is expected at exactly what it holds, so it carries no residual.
#' That is not an approximation: a pair that answered the same way under both
#' conditions says nothing about which condition raises the response, which is
#' the same fact that drops the concordant cells out of McNemar's statistic.
#'
#' And the residuals of a discordant pair square and sum to that statistic. For a
#' 2 x 2 table each of `b` and `c` is expected at `(b + c) / 2`, so the two
#' Pearson residuals square to `(b - c)^2 / (2 * (b + c))` each and to
#' `(b - c)^2 / (b + c)` together, which is McNemar's uncorrected chi-square. The
#' cell table and the p-value beside it are therefore about the same hypothesis.
#'
#' @param counts A square two-dimensional table or matrix of counts.
#'
#' @keywords internal
#' @noRd
sa_expected_symmetry <- function(counts) {
  if (nrow(counts) != ncol(counts)) {
    stop("internal error: symmetry is a claim about a square table, and this ",
         "one is ", paste(dim(counts), collapse = " x "), ".", call. = FALSE)
  }
  out <- (counts + t(counts)) / 2
  dimnames(out) <- dimnames(counts)
  out
}


#' Expand a contingency table into the canonical cell table
#'
#' A matrix vectorises down its columns and [expand.grid()] varies its first
#' argument fastest, so the labels and the counts line up without either being
#' matched by name.
#'
#' @param counts A two-dimensional [table()] or matrix with dimnames.
#' @param null Which null hypothesis the expected counts state. One of
#'   `sa_categorical_nulls()`.
#'
#' @return data.frame with the columns of `sa_categorical_cell_columns()`.
#'
#' @keywords internal
#' @noRd
sa_categorical_cells <- function(counts, null = "independence") {
  if (!isTRUE(null %in% sa_categorical_nulls())) {
    stop("internal error: `null` must name one of ",
         paste(sa_categorical_nulls(), collapse = ", "), ".", call. = FALSE)
  }

  n <- sum(counts)
  row_n <- rowSums(counts)
  col_n <- colSums(counts)

  expected <- if (identical(null, "symmetry")) {
    sa_expected_symmetry(counts)
  } else {
    sa_expected_independence(counts)
  }

  # Dividing by an empty margin is not a residual of zero, it is a residual that
  # does not exist, so the non-finite results are kept as NA rather than passed
  # on as NaN.
  residual <- sa_finite_or_na((counts - expected) / sqrt(expected))

  # The variance correction the standardized residual divides by is derived for a
  # two-way table held against its own margins. Under symmetry the comparison is
  # a cell against its transpose and that correction has no counterpart, so the
  # column is NA there rather than a number that looks referable to a standard
  # normal and is not.
  std_residual <- if (identical(null, "symmetry")) {
    array(NA_real_, dim = dim(counts))
  } else {
    sa_finite_or_na(
      (counts - expected) /
        sqrt(expected * outer(1 - row_n / n, 1 - col_n / n))
    )
  }

  grid <- expand.grid(row_level = rownames(counts),
                      col_level = colnames(counts),
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)

  out <- data.frame(
    row_level    = grid$row_level,
    col_level    = grid$col_level,
    observed     = as.numeric(counts),
    expected     = as.numeric(expected),
    residual     = as.numeric(residual),
    std_residual = as.numeric(std_residual),
    prop_total   = as.numeric(counts) / n,
    prop_row     = sa_finite_or_na(as.numeric(counts / row_n)),
    prop_col     = sa_finite_or_na(as.numeric(t(t(counts) / col_n))),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}


#' Replace non-finite values with NA, keeping the shape
#'
#' @keywords internal
#' @noRd
sa_finite_or_na <- function(x) {
  x[!is.finite(x)] <- NA_real_
  x
}


#' The approximation rule this design rests on
#'
#' Reported rather than enforced, which is the same choice
#' [diagnose_distribution()] makes about normality: a failed check does not change
#' which tests run, it changes which of them deserves the weight. Here it is a
#' short walk from the check to the answer, because the test that does not need
#' the approximation is already in the result beside it.
#'
#' Each design rests on a different rule, so the three builders below report
#' different numbers, and every one of them ends in the same two fields:
#' `approx_ok`, and a `note` that says in one sentence what to read instead when
#' it is `FALSE`. That is what lets `print.sa_categorical()` report the check
#' without knowing which design produced it. Each `rule` is an id
#' `Configuration/registry/assumptions.yaml` records.
#'
#' The independent rule is `expected_count_min`: every expected count at least 5,
#' or at most a fifth of the cells below 5 with none below 1.
#'
#' @keywords internal
#' @noRd
sa_diagnose_expected <- function(cells) {
  expected <- cells$expected
  n_cells <- length(expected)
  n_lt5 <- sum(expected < 5)
  n_lt1 <- sum(expected < 1)
  ok <- n_lt5 == 0L || (n_lt1 == 0L && n_lt5 / n_cells <= 0.2)

  list(
    rule              = "expected_count_min",
    min_expected      = min(expected),
    n_cells           = n_cells,
    n_expected_lt5    = n_lt5,
    prop_expected_lt5 = n_lt5 / n_cells,
    approx_ok         = ok,
    note              = paste0(
      "The smallest expected count is ", sa_fmt_est(min(expected)), " and ",
      n_lt5, " of ", n_cells, " cell(s) fall below 5, so the chi-square ",
      "approximation is doubtful here. Read `$tests$fisher_test`, which needs no ",
      "approximation, or set `simulate_p_value = TRUE`."
    )
  )
}


#' The discordant pair rule McNemar's approximation rests on
#'
#' `assumptions.yaml` records it as `discordant_pair_count`: at least 25
#' discordant pairs for the chi-square approximation. The exact branch is what
#' runs below that under the default `exact = NULL`, so this reports whether the
#' approximation *would* have been sound rather than whether the answer is.
#'
#' @keywords internal
#' @noRd
sa_diagnose_discordance <- function(n_discordant) {
  ok <- n_discordant >= 25L

  list(
    rule         = "discordant_pair_count",
    n_discordant = n_discordant,
    approx_ok    = ok,
    note         = paste0(
      "Only ", n_discordant, " discordant pair(s) carry the comparison, which ",
      "is below the 25 the chi-square approximation asks for. The exact ",
      "binomial branch is the one to read; `parameters$exact` says whether it ",
      "ran."
    )
  )
}


#' The sample size rule Cochran's Q rests on
#'
#' Q is referred to a chi-square distribution on `k - 1` degrees of freedom, and
#' the usual rule of thumb for that approximation is that the number of subjects
#' times the number of conditions reaches 24. `assumptions.yaml` records it as
#' `sample_size_repeated`. There is no exact test in the result to fall back on,
#' so the note says what the number means rather than where to look instead.
#'
#' @keywords internal
#' @noRd
sa_diagnose_repeated <- function(n_subjects, k) {
  cells <- n_subjects * k
  ok <- cells >= 24L

  list(
    rule         = "sample_size_repeated",
    n_subjects   = n_subjects,
    n_conditions = k,
    n_cells      = cells,
    approx_ok    = ok,
    note         = paste0(
      n_subjects, " subject(s) over ", k, " condition(s) is ", cells,
      " observation(s), below the 24 the chi-square approximation for Q asks ",
      "for, so read its p-value as an indication rather than as a rate."
    )
  )
}
