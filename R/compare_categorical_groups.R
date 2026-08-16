#' Test a contingency table with every applicable test at once
#'
#' Crosses two categorical variables into a contingency table and returns an
#' asymptotic and an exact test of it side by side, together with the measures of
#' how strong the association is. Nothing is chosen on the user's behalf:
#' reporting both makes disagreement between them visible, which is the situation
#' where the choice between an approximation and an exact enumeration actually
#' matters.
#'
#' This is the one scenario in the package with no feature axis. Every other
#' comparison asks its question once per numeric column and returns a table with
#' one row per column; a contingency table is asked about as a whole, so the
#' result carries the table itself in `$cells` and one row per test in `$tests`.
#' That is also why the result is not an `sa_comparison` and does not go through
#' [estimate_significance()]: a volcano plot needs a signed effect per feature and
#' this scenario has no feature axis to carry one.
#' [estimate_categorical_significance()] is the counterpart that reads it, one
#' verdict per cell. See "Why this is not a comparison result".
#'
#' @section One test, three questions:
#' The chi-square statistic answers what look like three different questions, and
#' the arithmetic is the same in each.
#'
#' \describe{
#'   \item{Independence}{Both variables were measured on one sample: is there any
#'     association between them? This is the reading the default label uses.}
#'   \item{Homogeneity}{One variable says which sample a row came from and the
#'     other is what was measured: do the samples share the same distribution
#'     over the categories? The sampling scheme differs, the null hypothesis is
#'     stated differently, and the test is the same one.}
#'   \item{Goodness of fit}{One variable against a set of expected proportions.
#'     This one is genuinely different, being about a single margin rather than a
#'     cross-classification, and it is not implemented here.}
#' }
#'
#' Which of the first two a caller means is a fact about their sampling scheme,
#' not about their data, and the data cannot be inspected to find out. So this
#' function does not ask, and both readings are the same call.
#'
#' @section Which tests run:
#' Decided by `paired` and by how many variables `category_lv` names, so there is
#' no argument naming a test.
#'
#' | design | null hypothesis | tests |
#' |---|---|---|
#' | independent, two variables | independence | Chi-square test of independence, Fisher's exact test |
#' | matched, two binary conditions | symmetry | McNemar's test |
#' | matched, three or more binary conditions | marginal homogeneity | Cochran's Q test |
#'
#' A matched design reads the columns as repeated measurements of one thing on the
#' same row, so pairing is by row and there is no `id` argument. It also requires
#' the levels to be binary: the tests of symmetry that generalise McNemar's test
#' past two levels, Bowker's and Stuart-Maxwell's, are not implemented, and a
#' three-level matched design is an error naming them rather than a table quietly
#' collapsed to two levels.
#'
#' @section What the cells are expected to hold:
#' `design$null` names the hypothesis the whole result is about, and
#' `$cells$expected` is read under it. This matters because a contingency table
#' can be held against more than one null and the three designs above hold it
#' against three different ones.
#'
#' Under **independence** a cell is expected at the product of its margins over
#' the total. Under **symmetry** it is expected at the average of it and its
#' transpose, so the diagonal is expected at exactly what it holds and carries no
#' residual: a pair that answered the same way under both conditions says nothing
#' about which condition raises the response. Under **marginal homogeneity** every
#' condition is expected to show the pooled response rate, which on a
#' condition-by-response table is the same arithmetic as independence and a
#' different claim about the world.
#'
#' The residuals follow the same hypothesis, which is what keeps them and the
#' p-value beside them talking about one thing. In a matched 2 x 2 the squared
#' Pearson residuals of the two discordant cells sum to `(b - c)^2 / (b + c)`,
#' which is McNemar's uncorrected statistic exactly; the cell table and the test
#' are the same arithmetic read two ways. [draw_mosaic_plot()] shades on those
#' residuals, so the picture is about the hypothesis the result tested rather than
#' about independence in every case.
#'
#' `std_residual` is the exception. Its variance correction is derived for a table
#' held against its own margins and has no counterpart under symmetry, so it is
#' `NA` there rather than a number that looks referable to a standard normal and
#' is not.
#'
#' @section Why this is not a comparison result:
#' `sa_new_comparison()` holds every table in a comparison to one row per feature,
#' and [estimate_significance()] reads one `log2fc` per row off that axis. Neither
#' exists here. There is one table, so there is one p-value and no multiplicity to
#' adjust across, which is why no `pval_adj` column is carried. And an association
#' has no sign: `cramers_v` says how far the table sits from independence but not
#' in which direction, because past a 2 x 2 table there is no single direction to
#' name. On a 2 x 2 table there is one, and `odds_ratio` reports it.
#'
#' The result is therefore an `sa_categorical`, which keeps the vocabulary of a
#' comparison -- `design`, `parameters`, `tests` beside `test_info`, `metadata` --
#' without claiming a contract it cannot meet. [draw_mosaic_plot()] is what reads
#' it, in the place [draw_volcano_plot()] holds for the numeric scenarios.
#'
#' Every sentence above is about the table. A **cell** does have both axes:
#' `observed / expected` is a signed departure defined on a table of any shape,
#' and `std_residual` beside it is referred to a standard normal, so the cells of
#' one table are a family to adjust across.
#' [estimate_categorical_significance()] is what reads that axis, and it is a
#' function of `$cells` alone, so no slot is carried for it.
#'
#' @section Direction:
#' Set once, by the order of the levels in `category_lv`, and `control_label` is
#' the second way of stating it. Every quantity that has a direction follows it.
#'
#' `odds_ratio` is above 1 when the **second** level of each variable occurs with
#' the second level of the other more often than independence would give, and
#' `phi_coefficient` is above zero in the same situation. Pointing
#' `control_label` at the other level of either variable inverts both, and
#' pointing it at the other level of both leaves them where they were. In a
#' matched design the second level is the response the conditions are compared on,
#' so `odds_ratio_paired` is above 1 when the later condition raises it.
#'
#' @param data A data.frame (or matrix) whose columns are the categorical
#'   variables. Unlike the other comparison scenarios there is no `feats`
#'   argument: the columns are what is tested rather than what is measured, so a
#'   factor, a character vector, a logical or a 0/1 coded column all read as
#'   categorical and are taken as their labels.
#' @param category_lv Named list giving the levels of each variable, with the
#'   reference level first, or `NULL` to take every column of `data` with its
#'   levels in sorted order. Naming it does three things sorting cannot: it picks
#'   which columns take part, it fixes which level is the reference the odds ratio
#'   is read against, and it drops the rows belonging to any level it leaves out.
#'   Every level it names has to be present in `data`.
#' @param control_label The level to hold as the reference. For an independent
#'   design this is one name per variable it points at, as a named list
#'   (`list(smoker = "n", grade = "low")`) or as a named character vector
#'   (`c(smoker = "n")`); the two shapes say the same thing, and a variable it
#'   says nothing about is left as it arrived. A matched design has one level set
#'   shared by every condition, so there it is a single level name. Defaults to
#'   `NULL`, the first level of every variable, which under `category_lv = NULL`
#'   is the first in sorted order.
#' @param paired Logical. If `TRUE`, the columns are read as repeated conditions
#'   measured on the same row rather than as different variables. See "Which tests
#'   run".
#' @param conf_level Confidence level for the association intervals and for the
#'   conditional odds ratio of Fisher's exact test.
#' @param correct Whether to apply the continuity correction to the chi-square
#'   approximation. [stats::chisq.test()] applies it to 2 x 2 tables only, so it
#'   changes nothing on a larger one, and the exact branch of McNemar's test needs
#'   none.
#' @param exact Read only by McNemar's test. `TRUE` for the exact binomial test on
#'   the discordant pairs, `FALSE` for the chi-square approximation, or `NULL` to
#'   take the exact test when there are fewer than 25 discordant pairs.
#'   `parameters$exact` records which one ran.
#' @param simulate_p_value Read only by the tests of an independent design.
#'   Replaces the chi-square approximation with a Monte Carlo p-value, and is the
#'   way to get an answer out of Fisher's test on a large r x c table that cannot
#'   be enumerated.
#' @param n_resamples How many tables the Monte Carlo p-value is taken over.
#' @param max_levels How many levels a variable may take before it is refused as a
#'   category. `as.character()` turns a continuous measurement into one label per
#'   observation, which makes a table with a cell per row and no test of
#'   association to run on it; this is the ceiling that catches that rather than
#'   letting it through. Checked against the levels actually used, so naming three
#'   levels of a fifty-valued column in `category_lv` is a way through. Raise it
#'   for a variable that genuinely has this many.
#' @param seed Seed for the Monte Carlo p-value, restored on exit, so a simulated
#'   p-value is reproducible without the caller's random stream being disturbed.
#' @param diagnose Logical. If `TRUE`, the rule the reported approximation rests
#'   on is attached as `$diagnostics`.
#'
#' @return An `sa_categorical` object: a plain list, so it survives being written
#'   out as JSON and read back in another language, with an S3 class on top that
#'   supplies [print()], [plot()] and [as.table()]. Its elements are
#'
#'   \describe{
#'     \item{`analysis`}{`"categorical_comparison"`.}
#'     \item{`variables`}{The variable names, in the order `category_lv` fixed.}
#'     \item{`design`}{`category_lv`, `null`, `paired`, `pairing`, `dim`,
#'       `row_var`, `col_var`, `n_used`, `n_dropped` and `n_incomplete`. `null` is
#'       `"independence"`, `"symmetry"` or `"marginal_homogeneity"`, and it is what
#'       `expected` and the residuals are read under.}
#'     \item{`parameters`}{The analysis choices as used, so `exact` says which
#'       branch of McNemar's test ran rather than what was passed.}
#'     \item{`cells`}{One row per cell of the table: `row_level`, `col_level`,
#'       `observed`, `expected`, `residual`, `std_residual`, `prop_total`,
#'       `prop_row` and `prop_col`. This is the canonical form of the table, and
#'       what [draw_mosaic_plot()] reads. `as.table()` folds it back into a
#'       [table()], which is the shape to read it in rather than a second copy to
#'       keep.}
#'     \item{`tests`}{One one-row data.frame per test, named `chisq_test`,
#'       `fisher_test`, `mcnemar_test` or `cochran_q`, each carrying `n_used`,
#'       `statistic`, `df`, `pval`, `lower_conf` and `upper_conf`. There is no
#'       `pval_adj`.}
#'     \item{`test_info`}{The method `id`, a readable `label` and whether the test
#'       is a matched one, per element of `tests`.}
#'     \item{`association`}{One row per measure: `measure`, `estimate`,
#'       `lower_conf`, `upper_conf`. Which measures are defined depends on the
#'       design and on the size of the table.}
#'     \item{`diagnostics`}{The approximation rule this design rests on, or
#'       `NULL`.}
#'     \item{`metadata`}{`package_version`, `r_version`, `platform` and an
#'       ISO-8601 `timestamp`.}
#'   }
#'
#' @section Which measures are defined:
#' An independent table always reports `cramers_v` and
#' `contingency_coefficient`, and a 2 x 2 one adds `phi_coefficient` and
#' `odds_ratio`. A matched 2 x 2 reports `odds_ratio_paired`,
#' `risk_difference_paired` and `cohens_g`, all three read off the discordant
#' cells. Three or more matched conditions report `kendalls_w`.
#'
#' Only `odds_ratio` and the paired measures carry an interval. The others are
#' functions of the chi-square statistic whose sampling distribution has no
#' closed-form interval, so their `lower_conf` and `upper_conf` are `NA`: the
#' columns exist because every result table in the package carries them, not
#' because every number in them is finite.
#'
#' Every measure is built from the **uncorrected** chi-square statistic, whatever
#' `correct` was set to. Yates' correction is about referring a discrete statistic
#' to a continuous distribution, which is a statement about a p-value; letting it
#' into an effect size would make the reported strength of an association depend
#' on a choice made about its tail probability.
#'
#' @details
#' A row is dropped for one of two reasons, and the two are counted apart. A row
#' naming a level `category_lv` leaves out was measured and excluded, and appears
#' as `design$n_dropped`. A row missing a value in any variable was not measured,
#' and appears as `design$n_incomplete`; a table needs the whole row, so the
#' deletion is listwise.
#'
#' The warning [stats::chisq.test()] raises about small expected counts is not
#' passed on. `$diagnostics` states the same fact as a number, and the exact test
#' that does not need the approximation is already in the same result, so a
#' warning would be a less precise version of what is sitting next to it.
#'
#' Fisher's exact test enumerates every table with the observed margins, and on a
#' large r x c one there are more of those than the algorithm's workspace holds.
#' That is a limit of the enumeration rather than a fault in the data, so
#' `$tests$fisher_test$pval` comes back `NA` with a message saying so, instead of
#' the whole call failing and taking the chi-square result with it.
#' `simulate_p_value = TRUE` is the way to get an answer there.
#'
#' @seealso [simulate_categorical_groups()] for a table whose association is
#'   known, [estimate_categorical_significance()] to reduce the result to a
#'   verdict per cell, [draw_mosaic_plot()] to draw it, and
#'   [compare_multiple_groups()] for the case where the measurement is numeric and
#'   only the grouping is categorical.
#'
#' @references
#' Pearson, K. (1900). On the criterion that a given system of deviations from the
#' probable in the case of a correlated system of variables is such that it can be
#' reasonably supposed to have arisen from random sampling. *Philosophical
#' Magazine*, 50(302), 157-175.
#'
#' Fisher, R. A. (1935). The logic of inductive inference. *Journal of the Royal
#' Statistical Society*, 98(1), 39-82.
#'
#' McNemar, Q. (1947). Note on the sampling error of the difference between
#' correlated proportions or percentages. *Psychometrika*, 12(2), 153-157.
#'
#' Cochran, W. G. (1950). The comparison of percentages in matched samples.
#' *Biometrika*, 37(3-4), 256-266.
#'
#' Agresti, A. (2002). *Categorical Data Analysis*, 2nd ed. Wiley.
#'
#' @examples
#' ## Both columns are categorical, so `data` is the whole input.
#' smoking <- data.frame(
#'   smoker = rep(c("y", "n"), each = 60),
#'   grade  = c(rep(c("high", "mid", "low"), c(10, 20, 30)),
#'              rep(c("high", "mid", "low"), c(30, 20, 10)))
#' )
#' res <- compare_categorical_groups(smoking)
#' res
#'
#' ## The table itself, and the cells that made the statistic what it is.
#' as.table(res)
#' res$cells[c("row_level", "col_level", "observed", "expected", "std_residual")]
#'
#' ## How strong the association is, which no test reports.
#' res$association
#'
#' ## `category_lv` picks the levels and their order, and drops the rest. Two
#' ## levels of `grade` make a 2 x 2 table, which is where the odds ratio exists.
#' two_by_two <- compare_categorical_groups(
#'   smoking,
#'   category_lv = list(smoker = c("n", "y"), grade = c("low", "high"))
#' )
#' subset(two_by_two$association, measure == "odds_ratio")
#'
#' ## `control_label` restates the reference without the levels being retyped,
#' ## and pointing it at one variable inverts the odds ratio.
#' flipped <- compare_categorical_groups(
#'   smoking,
#'   category_lv   = list(smoker = c("n", "y"), grade = c("low", "high")),
#'   control_label = c(smoker = "y")
#' )
#' subset(flipped$association, measure == "odds_ratio")
#'
#' ## A matched design: the columns are the same question asked twice, so pairing
#' ## is by row, the null is symmetry rather than independence, and McNemar's test
#' ## reads only the discordant pairs.
#' before_after <- data.frame(
#'   before = rep(c("pass", "fail"), c(20, 30)),
#'   after  = c(rep(c("pass", "fail"), c(18, 2)), rep(c("pass", "fail"), c(14, 16)))
#' )
#' matched <- compare_categorical_groups(before_after, paired = TRUE)
#' matched$design$null
#'
#' ## Which is what the cell table is read under: the diagonal is expected at
#' ## exactly what it holds, so only the discordant cells carry a residual.
#' matched$cells[c("row_level", "col_level", "observed", "expected", "residual")]
#'
#' ## A simulated table hands over exactly the arguments this function takes.
#' sim <- simulate_categorical_groups(n_samples = 400, assoc = 0.4, seed = 1)
#' fit <- do.call(compare_categorical_groups, sim$args)
#' cbind(planted = sim$truth$cramers_v,
#'       estimated = fit$association$estimate[1])
#'
#' @export
compare_categorical_groups <- function(data,
                                       category_lv = NULL,
                                       control_label = NULL,
                                       paired = FALSE,
                                       conf_level = 0.95,
                                       correct = TRUE,
                                       exact = NULL,
                                       simulate_p_value = FALSE,
                                       n_resamples = 9999,
                                       max_levels = 20L,
                                       seed = NULL,
                                       diagnose = TRUE) {

  sa_check_flag(paired, "paired")
  sa_check_flag(correct, "correct")
  sa_check_flag(simulate_p_value, "simulate_p_value")
  sa_check_flag(diagnose, "diagnose")
  sa_check_scalar_num(conf_level, "conf_level", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)
  n_resamples <- sa_check_count(n_resamples, "n_resamples", 199)
  max_levels <- sa_check_count(max_levels, "max_levels", 2)
  if (!is.null(exact)) {
    sa_check_flag(exact, "exact")
  }
  if (!is.null(exact) && !paired) {
    warning("`exact` is only read by McNemar's test and is ignored by an ",
            "independent design. `simulate_p_value` is what replaces the ",
            "chi-square approximation there.", call. = FALSE)
  }
  if (simulate_p_value && paired) {
    warning("`simulate_p_value` is only read by the tests of an independent ",
            "design and is ignored by a matched one. `exact` is what chooses ",
            "the exact branch of McNemar's test.", call. = FALSE)
  }

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  input <- sa_validate_categorical_input(data, category_lv, control_label,
                                         paired, max_levels)
  if (input$n_dropped > 0L) {
    message("Dropped ", input$n_dropped,
            " row(s) belonging to a level outside `category_lv`.")
  }
  if (input$n_incomplete > 0L) {
    message("Dropped ", input$n_incomplete,
            " row(s) missing a value in one of the variables; a contingency ",
            "table needs the whole row.")
  }

  layout <- if (paired) {
    sa_categorical_matched(input, correct, exact, conf_level)
  } else {
    sa_categorical_independent(input, correct, simulate_p_value, n_resamples,
                               conf_level)
  }

  # Raised here rather than inside `sa_odds_ratio()`, so that every message this
  # scenario can produce is raised in one place and the kernel stays a function
  # of its arguments.
  if (sa_has_zero_cell(layout$counts)) {
    message("A zero cell leaves the odds ratio undefined, so the ",
            "Haldane-Anscombe correction of 0.5 per cell was applied to it. ",
            "The tests read the table as it is.")
  }
  if (isFALSE(layout$enumerated)) {
    message("Fisher's exact test could not enumerate a ",
            paste(dim(layout$counts), collapse = " x "), " table of ",
            sum(layout$counts), " observation(s), so its p-value is NA. The ",
            "chi-square test in the same result was computed; set ",
            "`simulate_p_value = TRUE` for the Monte Carlo variant of both.")
  }

  diagnostics <- if (diagnose) layout$diagnostics else NULL
  if (diagnose && !layout$diagnostics$approx_ok) {
    message(layout$diagnostics$note)
  }

  sa_new_categorical(
    analysis  = "categorical_comparison",
    variables = input$variables,
    design    = list(
      category_lv  = input$category_lv,
      null         = layout$null,
      paired       = paired,
      pairing      = if (paired) "row" else NA_character_,
      dim          = dim(layout$counts),
      row_var      = layout$row_var,
      col_var      = layout$col_var,
      n_used       = input$n_used,
      n_dropped    = input$n_dropped,
      n_incomplete = input$n_incomplete
    ),
    parameters = list(
      conf_level       = conf_level,
      correct          = correct,
      exact            = layout$exact_used,
      simulate_p_value = simulate_p_value,
      n_resamples      = if (simulate_p_value) n_resamples else NA_integer_,
      max_levels       = max_levels,
      seed             = if (is.null(seed)) NA_real_ else as.numeric(seed)
    ),
    cells       = layout$cells,
    tests       = layout$tests,
    test_info   = layout$test_info,
    association = layout$association,
    diagnostics = diagnostics
  )
}


#' Run the independent design: chi-square beside Fisher's exact test
#'
#' @param input The validated input.
#' @param correct,simulate_p_value,n_resamples,conf_level The arguments as
#'   received.
#'
#' @return List with `counts`, `cells`, `null`, `tests`, `test_info`,
#'   `association`, `row_var`, `col_var`, `exact_used` and `diagnostics`.
#'
#' @keywords internal
#' @noRd
sa_categorical_independent <- function(input, correct, simulate_p_value,
                                       n_resamples, conf_level) {
  variables <- input$variables
  if (length(variables) != 2L) {
    stop("an independent categorical comparison crosses exactly two variables, ",
         "and ", length(variables), " were given: ",
         paste(variables, collapse = ", "),
         ". Name the two in `category_lv`, or set `paired = TRUE` if the ",
         "columns are repeated measurements of one thing.", call. = FALSE)
  }

  counts <- sa_categorical_counts(input$data, variables)
  cells <- sa_categorical_cells(counts, "independence")

  fisher <- sa_fisher(counts, conf_level, simulate_p_value, n_resamples)
  tests <- list(
    chisq_test = sa_categorical_row(
      sa_chisq(counts, correct, simulate_p_value, n_resamples)
    ),
    # `enumerated` is a fact about whether the test ran rather than a finding
    # about the table, and an NA p-value already says it in the table itself.
    fisher_test = sa_categorical_row(fisher, drop = "enumerated")
  )

  list(
    counts    = counts,
    cells     = cells,
    null      = "independence",
    enumerated = as.logical(fisher[["enumerated"]]),
    tests     = tests,
    test_info = list(
      chisq_test = list(
        id     = "chisq_independence",
        label  = paste0("Chi-square test of independence",
                        if (simulate_p_value) {
                          paste0(" (Monte Carlo, ", n_resamples, " resamples)")
                        } else if (correct && identical(dim(counts), c(2L, 2L))) {
                          " (Yates' continuity correction)"
                        } else {
                          ""
                        }),
        paired = FALSE
      ),
      fisher_test = list(
        id     = "fisher_exact",
        label  = "Fisher's exact test",
        paired = FALSE
      )
    ),
    association = sa_assoc_measures(counts, conf_level),
    row_var     = variables[1],
    col_var     = variables[2],
    exact_used  = NA,
    diagnostics = sa_diagnose_expected(cells)
  )
}


#' Run the matched design: McNemar's test or Cochran's Q
#'
#' The two branches read different tables against different nulls, and neither
#' difference is incidental. McNemar's test is about a square table crossing two
#' conditions against each other, held against symmetry, where the discordant
#' cells are the whole of the evidence. Cochran's Q has no such table to be about,
#' since three conditions do not cross into two dimensions; what it compares is
#' the response rate of each condition, so the table is one row per condition and
#' the null is that those rates agree.
#'
#' @param input The validated input.
#' @param correct,exact,conf_level The arguments as received.
#'
#' @keywords internal
#' @noRd
sa_categorical_matched <- function(input, correct, exact, conf_level) {
  variables <- input$variables
  levels <- input$category_lv[[1]]

  if (length(levels) != 2L) {
    stop("`paired = TRUE` needs binary conditions, and the shared level set ",
         "holds ", length(levels), ": ", paste(levels, collapse = ", "),
         ". The tests of symmetry that generalise McNemar's test past two ",
         "levels, Bowker's and Stuart-Maxwell's, are not implemented in this ",
         "version. Name two levels in `category_lv` to reduce the design ",
         "rather than have it collapsed here.", call. = FALSE)
  }

  if (length(variables) == 2L) {
    counts <- sa_categorical_counts(input$data, variables)
    row <- sa_mcnemar(counts, correct, exact)
    exact_used <- as.logical(row[["exact_used"]])

    return(list(
      counts    = counts,
      cells     = sa_categorical_cells(counts, "symmetry"),
      null      = "symmetry",
      # `exact_used` is a setting rather than a finding, and `parameters$exact`
      # is where a setting is recorded, so it does not go on into the table.
      tests     = list(mcnemar_test = sa_categorical_row(row,
                                                         drop = "exact_used")),
      test_info = list(mcnemar_test = list(
        id     = "mcnemar_test",
        label  = if (exact_used) {
          "McNemar's test (exact binomial on the discordant pairs)"
        } else if (correct) {
          "McNemar's test (continuity corrected)"
        } else {
          "McNemar's test"
        },
        paired = TRUE
      )),
      association = sa_assoc_measures_paired(counts, conf_level),
      row_var     = variables[1],
      col_var     = variables[2],
      exact_used  = exact_used,
      diagnostics = sa_diagnose_discordance(row[["n_discordant"]])
    ))
  }

  # The second level is the response, which is the direction rule every other
  # scenario reads: `category_lv[[1]][1]` is the reference and what is counted is
  # departure from it.
  responses <- vapply(variables, function(nm) {
    as.integer(input$data[[nm]] == levels[2])
  }, integer(nrow(input$data)))
  responses <- matrix(responses, nrow = nrow(input$data),
                      dimnames = list(NULL, variables))

  row <- sa_cochran_q(responses)
  counts <- sa_categorical_condition_counts(input$data, variables, levels)

  list(
    counts    = counts,
    cells     = sa_categorical_cells(counts, "marginal_homogeneity"),
    null      = "marginal_homogeneity",
    tests     = list(cochran_q = sa_categorical_row(row)),
    test_info = list(cochran_q = list(
      id     = "cochran_q",
      label  = paste0("Cochran's Q test over ", length(variables),
                      " repeated condition(s)"),
      paired = TRUE
    )),
    association = sa_assoc_measures_repeated(row[["statistic"]],
                                             row[["n_used"]],
                                             length(variables)),
    row_var     = "condition",
    col_var     = "response",
    exact_used  = NA,
    diagnostics = sa_diagnose_repeated(row[["n_used"]], length(variables))
  )
}


#' Turn a kernel row into the one-row table `$tests` holds
#'
#' @param row The named numeric vector a kernel returned.
#' @param drop Names the kernel reported for the scenario to record elsewhere,
#'   which are not findings about the table and so do not belong in it.
#'
#' @keywords internal
#' @noRd
sa_categorical_row <- function(row, drop = character(0)) {
  row <- row[setdiff(names(row), drop)]
  out <- as.data.frame(as.list(row), stringsAsFactors = FALSE)
  absent <- setdiff(sa_categorical_test_columns(), names(out))
  if (length(absent) > 0L) {
    stop("internal error: kernel row is missing column(s): ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }
  rownames(out) <- NULL
  out[c(sa_categorical_test_columns(),
        setdiff(names(out), sa_categorical_test_columns()))]
}
