#' Analyse a crossed-factor design as one model
#'
#' The factorial counterpart of [compare_multiple_groups()], and the one place in
#' the package where a scenario function fits **one** model rather than running
#' every applicable test side by side. Two crossed factors are a two-way ANOVA,
#' three are a three-way ANOVA and more than three are a factorial ANOVA, and
#' those are three names for the same fully crossed linear model rather than three
#' procedures to choose between. Which name applies follows from
#' `length(factor_lv)`, so there is no argument for it; the answer is reported in
#' `design$anova_type` and in the readable label of the test.
#'
#' What multiplies instead is the number of answers. One factor asks whether its
#' levels are alike; crossing a second one asks that of each factor and asks
#' whether either one's effect depends on the other, and those questions fail
#' separately. A single p-value cannot carry them, so the result has an axis the
#' sibling functions do not: `$terms`, one row per feature and model term.
#'
#' @param data A data.frame (or matrix) in wide format, one row per observation
#'   and one column per feature.
#' @param feats Character vector of numeric column names in `data` to test. One
#'   row of `$tests` and `$effect` per entry, and one row of `$terms` per entry
#'   and term.
#' @param factors Named list of the crossed factors, each entry either the name of
#'   a column of `data` or a vector with one entry per row of it. Its length is how
#'   many factors are crossed, and there have to be at least two: one factor is
#'   [compare_multiple_groups()]. The **first factor is the primary one** and the
#'   order of the list is the order the terms are listed in, which is also the
#'   order Type I sums of squares are taken in.
#' @param factor_lv Named list giving the levels of each factor, with the
#'   reference level first, or `NULL` to take the levels from the data in sorted
#'   order. Naming it does two things sorting cannot: it fixes which level is the
#'   reference that `effect` and the cell labels are read against, and it drops
#'   the rows belonging to any level it leaves out. Every level it names has to be
#'   present in `factors`. `control_label` fixes the reference on its own, for
#'   when that is all there is to say.
#' @param control_label The level to hold as the reference, one name per factor
#'   it points at, as a named list (`list(treatment = "control", sex = "male")`)
#'   or as a named character vector (`c(treatment = "control")`). The two shapes
#'   say the same thing. The level it names moves to the front of its own factor
#'   and the rest keep the order they were given, so the reference cell is the one
#'   where every factor sits at the level named here and every fold change is
#'   read against it. A factor it says nothing about is left as it arrived, which
#'   is what makes pointing one factor of three a sentence rather than a rewrite
#'   of all three. Defaults to `NULL`, the first level of every factor, which
#'   under `factor_lv = NULL` is the first in sorted order.
#' @param within Names of the factors measured within subjects. **Not implemented
#'   in this version**: a non-empty value is an error rather than a design that is
#'   silently analysed as though the repeated measurements were independent. The
#'   argument exists so that the error names what is missing, which an
#'   `unused argument (within)` from [do.call()] would not.
#' @param id Subject identifier with one entry per row of `data`. Used only by a
#'   within-subject design, so it is ignored with a warning here.
#' @param conf_level Confidence level for the post-hoc intervals.
#' @param ss_type Which sums of squares the term tests are built from, `"III"`,
#'   `"II"` or `"I"`. They agree on a balanced design and differ on an unbalanced
#'   one; see "Unbalanced cells" below.
#' @param posthoc Logical. If `FALSE`, no pairwise stage runs and the result
#'   carries no `$posthoc` slot.
#' @param posthoc_alpha A contrast is computed when the term it belongs to has a
#'   `pval_adj` at or below this value. Set it to 1 to compare every level of
#'   every factor regardless of the term tests.
#' @param posthoc_scope Which contrasts the pairwise stage covers: `"marginal"`
#'   for the comparisons a main effect is about, `"simple"` for the same
#'   comparisons inside one combination of the other factors, or `"both"`.
#' @param fc_mean Which centre the fold change divides, `"arith"` for the
#'   arithmetic mean or `"geom"` for the geometric mean. Defaults to `"geom"`
#'   when `input_scale = "log2"` and to `"arith"` otherwise.
#' @param input_scale The scale `data` arrives on, `"raw"` or `"log2"`. On the
#'   log2 scale each observation is raised back through `2^x` before the centres
#'   are taken, so the ratios mean what they do for raw input. This changes the
#'   `effect` table only, never the tests. See [compare_two_groups()] for the
#'   full account.
#' @param p_adjust Multiplicity adjustment passed to [stats::p.adjust()], applied
#'   along the feature axis. Use `"none"` to disable. In `$terms` it is applied
#'   within each term rather than over the whole table, for the reason given under
#'   "Two axes" below.
#' @param diagnose Logical. If `TRUE`, the normality, homogeneity of variance and
#'   outlier checks the model rests on are attached as `$diagnostics`, taken over
#'   the cells.
#'
#' @return A `sa_comparison` object with the layout [compare_multiple_groups()]
#'   returns, minus `$pairwise` and plus `$terms` and `$cells`. Its elements are
#'
#'   \describe{
#'     \item{`analysis`}{`"factorial_comparison"`.}
#'     \item{`features`}{The feature names, in the row order every table uses.}
#'     \item{`design`}{`factor_lv`, `anova_type`, `n_factors`, `group_lv` (the
#'       cell labels), `cell_n`, `n_empty_cells`, `paired`, `pairing`, `within`,
#'       `n_dropped` and `unmatched_ids`.}
#'     \item{`parameters`}{The analysis choices as used, plus `n_posthoc`, the
#'       number of features that entered the pairwise stage.}
#'     \item{`effect`}{One row per feature: `n_used`, `n_cells`, `ref_center`,
#'       `extreme_cell`, `extreme_center`, `fold_change` and `log2fc`. The
#'       reference cell is the one where every factor sits at its first level,
#'       which `factor_lv` or `control_label` states, and
#'       `extreme_cell` is whichever cell sits furthest from it on the log2 scale,
#'       the same definition `simulate_factorial_groups()` records in
#'       `truth$extreme_cell`.}
#'     \item{`tests`}{One table, `anova_test`, holding the whole-model F test with
#'       one row per feature.}
#'     \item{`terms`}{One row per feature and model term, holding `features`,
#'       `terms`, `term_order`, `n_used`, `df`, `ss`, `ms`, `f_stat`, `df_error`,
#'       `eta_sq`, `partial_eta_sq`, `log2_effect`, `pval` and `pval_adj`.}
#'     \item{`cells`}{One row per feature and cell of the grid, holding
#'       `features`, one column per factor named after it and holding the level
#'       name, `cell` (the dot-joined label, as in `design$group_lv`), `n`,
#'       `mean`, `sd` and `se`. `mean` is the arithmetic cell mean the model was
#'       fitted on and `se` is `sqrt(ms_error / n)`, pooled over the whole model,
#'       so the variance of a marginal mean over a set of cells is recoverable
#'       from it. This is the table [draw_interaction_plot()] reads, and the one
#'       part of the grid the rest of the result cannot be read backwards into.}
#'     \item{`posthoc`}{`anova_test`, one row per feature and contrast, carrying
#'       the post-hoc contract columns plus `factor` and `stratum`. Absent when
#'       `posthoc = FALSE`.}
#'     \item{`test_info`}{The method `id`, a readable `label` naming the ANOVA
#'       that ran, and the post-hoc procedure that followed it.}
#'     \item{`diagnostics`}{Assumption checks over the cells, or `NULL`.}
#'     \item{`metadata`}{`package_version`, `r_version`, `platform` and an
#'       ISO-8601 `timestamp`.}
#'   }
#'
#' @section Two axes:
#' The whole-model test in `$tests$anova_test` asks whether a feature responds to
#' the design at all. That is one question per feature, and it is exactly the
#' one-way ANOVA that treats the cells as groups: a fully crossed model is the
#' cell means model written in another basis, so the two fit the same values and
#' leave the same residuals. Keeping it there is what lets
#' [estimate_significance()], [print()], [draw_forest_plot()] and
#' [draw_volcano_plot()] read a factorial result without knowing that it is one.
#'
#' Which **part** of the design a feature responds to is the question a crossed
#' model was fitted to answer, and it has one answer per term. Those live in
#' `$terms`, whose `terms` and `term_order` columns are the ones
#' `simulate_factorial_groups()` writes into `truth_term`, so a result and an
#' answer key merge on `c("features", "terms")` with neither side renamed.
#'
#' `pval_adj` is adjusted within each term, across features. A term is one family:
#' asking of five hundred features whether the treatment matters is five hundred
#' instances of one question, while asking whether the treatment matters and
#' whether it depends on sex are two questions and pooling them would correct each
#' for the other's multiplicity.
#'
#' There is no `$pairwise` slot. Its keys are contrast labels, and a factorial
#' design has a two-dimensional contrast axis, factor by stratum, that a flat
#' list keyed by label cannot hold without a naming convention. Until there is
#' one, `estimate_significance(by = "contrast")` reports that the result has no
#' pairwise stage, which is the true statement.
#'
#' @section The size of a term:
#' `eta_sq` and `partial_eta_sq` say how much of the variance a term accounts for,
#' and neither has a sign, so neither can be read as a direction. `log2_effect`
#' can: it is the largest ANOVA component of the term, with its sign, taken by
#' decomposing `log2()` of the same cell centres `effect` is built from. It is
#' what [estimate_significance()] puts on the effect axis in a term reading, and
#' so what [draw_volcano_plot()] plots a term panel against when passed the
#' verdict of `estimate_significance(fact, by = "term")`.
#'
#' A component is a deviation from what the rest of the model already predicts,
#' **not** a difference between two levels. A two-level factor whose levels differ
#' by one log2 unit has components `-0.5` and `+0.5`, so `|log2_effect|` is `0.5`
#' where the marginal fold change is 2. Exact and near ties in absolute value
#' (within the factorial tolerance) keep the earlier, reference cell's sign; do
#' not read that sign as a treatment-versus-control direction. The definition is
#' kept because it is the one `simulate_factorial_groups()` records in
#' `truth_term$max_abs_delta`, which makes the column scorable against an answer
#' key; the consequence is that cutoffs meant for a fold change,
#' `log2fc_cutoff = 1` among them, are stricter here than they look.
#'
#' Components do not depend on which cell is the reference, since subtracting one
#' cell from every cell cancels out of the deviation. That is the one place the
#' term axis is simpler than `effect$log2fc`, which is read against the reference
#' cell by definition.
#'
#' @section Unbalanced cells:
#' With equal cell sizes the three types of sums of squares are identical and
#' `ss_type` changes nothing. They part company when the cells are unequal,
#' because the factors are then no longer orthogonal and some of the variation
#' can be attributed to more than one of them.
#'
#' The default is **Type III**, each term adjusted for every other term. It is the
#' type whose main effects are statements about the levels rather than about how
#' many observations happened to land in each, which is what makes it agree with
#' the unweighted marginal means the post-hoc stage compares and with the
#' decomposition `simulate_factorial_groups()` plants. **Type II** adjusts a term
#' for every term that does not contain it, leaving a main effect unadjusted for
#' the interaction it is part of. **Type I** is sequential and therefore depends
#' on the order `factors` was written in; it is what [stats::aov()] reports, which
#' makes it the type to ask for when checking these numbers against an external
#' implementation on unbalanced data.
#'
#' A cell holding no usable observation leaves the crossed model with nothing to
#' estimate there, so that feature's rows are `NA` and the reason is reported in
#' one warning. `design$n_empty_cells` counts the cells that hold no rows at all.
#'
#' @section Marginal contrasts and simple effects:
#' A marginal contrast compares two levels of one factor with the other factors
#' averaged away, which is the comparison a main effect is a statement about. A
#' simple effect compares the same two levels inside one combination of the other
#' factors, and it is the only one of the two that means anything when an
#' interaction is real: if the treatment helps males and harms females, the
#' average of the two is a number that describes nobody. `stratum` tells them
#' apart, `NA` for a marginal contrast and the levels of the other factors joined
#' by a dot for a simple one, and `factor` names the factor being compared. Both
#' columns are the ones `truth_contrast` carries.
#'
#' The marginal mean is the **unweighted** mean of the cell means rather than the
#' mean of the observations, so a level is not pulled towards whichever
#' combination of the other factors was sampled most heavily. Every contrast is
#' scaled by the mean square error of the whole model, which is what makes the
#' pairwise stage consistent with the term tests instead of a second analysis of
#' the same data, and judged against the studentised range over the number of
#' levels of its own factor. Those p-values are already family-wise within a
#' block, so there is no `posthoc_p_adjust` argument to apply and
#' `parameters$posthoc_p_adjust` is `NA`.
#'
#' Which contrasts are computed is decided term by term. A marginal contrast of a
#' factor runs when that factor's main effect cleared `posthoc_alpha`, and a
#' simple effect runs when the interaction of the factor with the factors held
#' fixed cleared it. Gating everything on the whole-model test instead would
#' compare the levels of a factor that the model says has no effect, on the
#' strength of a different factor that has one.
#'
#' @details
#' Features that cannot be tested do not abort the run. Their rows are filled with
#' `NA` and all such features are reported together in a single warning.
#'
#' @seealso [compare_multiple_groups()] for a single factor,
#'   [simulate_factorial_groups()] for data whose answer per term is known, and
#'   [draw_forest_plot()] to draw the result.
#'
#' @references
#' Fisher, R. A. (1935). *The Design of Experiments*.
#'
#' Yates, F. (1934). The analysis of multiple classifications with unequal numbers
#' in the different classes. *Journal of the American Statistical Association*,
#' 29(185), 51-66.
#'
#' Tukey, J. W. (1949). Comparing individual means in the analysis of variance.
#' *Biometrics*, 5(2), 99-114.
#'
#' Kramer, C. Y. (1956). Extension of multiple range tests to group means with
#' unequal numbers of replications. *Biometrics*, 12(3), 307-310.
#'
#' @examples
#' ## Two crossed factors, so a two-way ANOVA. The factors name columns of `data`.
#' res <- compare_factorial_groups(
#'   data    = warpbreaks,
#'   feats   = "breaks",
#'   factors = list(wool = "wool", tension = "tension")
#' )
#' res
#'
#' ## The whole-model test, one row per feature: does this feature respond at all.
#' res$tests$anova_test
#'
#' ## The answer per term, which is what the design was crossed to get.
#' res$terms
#'
#' ## Every fold change is read against the cell where each factor sits at its
#' ## first level. `control_label` moves one of those levels without the rest
#' ## having to be listed, so the reference cell becomes A.M rather than A.L.
#' medium <- compare_factorial_groups(
#'   data          = warpbreaks,
#'   feats         = "breaks",
#'   factors       = list(wool = "wool", tension = "tension"),
#'   control_label = list(tension = "M"),
#'   posthoc       = FALSE
#' )
#' medium$design$group_lv[1]
#' medium$effect[c("ref_center", "extreme_cell", "log2fc")]
#'
#' ## The whole-model volcano names the reference cell on the x axis.
#' draw_volcano_plot(estimate_significance(res, log2fc_cutoff = 0.1))
#'
#' ## Marginal contrasts and simple effects, told apart by `stratum`.
#' subset(res$posthoc$anova_test, factor == "tension" & is.na(stratum))
#'
#' ## A simulated design hands over exactly the arguments this function takes,
#' ## and `truth_term` is the answer key for `$terms`.
#' sim <- simulate_factorial_groups(n_feats = 12, n_per_cell = 8, seed = 1)
#' fac <- do.call(compare_factorial_groups, sim$args)
#' scored <- merge(fac$terms, sim$truth_term, by = c("features", "terms"))
#' with(scored, table(called = pval_adj <= 0.05, planted = is_effect, terms))
#'
#' @export
compare_factorial_groups <- function(data,
                                     feats,
                                     factors,
                                     factor_lv = NULL,
                                     control_label = NULL,
                                     within = NULL,
                                     id = NULL,
                                     conf_level = 0.95,
                                     ss_type = c("III", "II", "I"),
                                     posthoc = TRUE,
                                     posthoc_alpha = 0.05,
                                     posthoc_scope = c("both", "marginal",
                                                       "simple"),
                                     fc_mean = c("arith", "geom"),
                                     input_scale = c("raw", "log2"),
                                     p_adjust = "BH",
                                     diagnose = TRUE) {

  input_scale <- match.arg(input_scale)
  fc_mean <- sa_resolve_fc_mean(fc_mean, input_scale, missing(fc_mean))
  ss_type <- match.arg(ss_type)
  posthoc_scope <- match.arg(posthoc_scope)
  sa_check_flag(posthoc, "posthoc")
  sa_check_flag(diagnose, "diagnose")
  sa_check_scalar_num(conf_level, "conf_level", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)
  sa_check_scalar_num(posthoc_alpha, "posthoc_alpha", 0, 1, lower_open = TRUE)
  sa_check_p_adjust(p_adjust, "p_adjust")

  if (length(within) > 0L) {
    stop("`within` names factor(s) measured within subjects, which this ",
         "version cannot analyse: ", paste(within, collapse = ", "),
         ". A mixed factorial model needs a second error stratum, which is ",
         "not implemented yet. Analyse the between-subject factors on one ",
         "level of the within one, or use compare_multiple_groups(paired = ",
         "TRUE) on the repeated factor alone.", call. = FALSE)
  }
  if (!is.null(id)) {
    warning("`id` is only used to match repeated measurements and is ignored ",
            "by a between-subject factorial analysis.", call. = FALSE)
  }

  input <- sa_validate_wide_input(data, feats, group = NULL, group_lv = NULL)
  data <- input$data
  feats <- input$feats

  design <- sa_fact_layout(data, factors, factor_lv, control_label)
  if (design$n_dropped > 0L) {
    message("Dropped ", design$n_dropped,
            " row(s) belonging to a level outside `factor_lv`.")
  }
  if (design$n_empty_cells > 0L) {
    message(design$n_empty_cells, " of ", design$n_cells,
            " cell(s) hold no observation, so no crossed model can be fitted ",
            "on them: ",
            paste(design$cell_label[design$cell_n == 0L], collapse = ", "), ".")
  }

  per_feature <- lapply(feats, function(f) {
    v <- data[[f]]
    stats::setNames(lapply(design$rows_of_cell, function(at) {
      x <- v[at]
      x[!is.na(x)]
    }), design$cell_label)
  })
  names(per_feature) <- feats

  centers <- sa_group_centers(per_feature, feats, design$cell_label, fc_mean,
                             paired = FALSE, input_scale = input_scale)
  effect <- sa_fact_effect(centers, feats, design$cell_label, fc_mean)

  plan <- sa_factorial_plan(design$factor_lv, design$cells, ss_type)
  label <- sa_fact_anova_label(design$anova_type, ss_type)

  # One fit per feature, both axes of it. The failures are held rather than
  # raised so that they can be re-raised inside `sa_feature_table()`, where they
  # become an NA row and one grouped warning for the whole scan.
  fits <- vector("list", length(feats))
  errors <- rep(NA_character_, length(feats))
  for (i in seq_along(feats)) {
    # `fits[i] <- list(...)` rather than `fits[[i]] <- ...`: assigning NULL with
    # `[[<-` deletes the element instead of emptying it, which shortens the list
    # under a feature that could not be fitted and leaves every later feature's
    # fit at the wrong index.
    fits[i] <- list(tryCatch(
      sa_factorial_anova(per_feature[[i]], plan),
      error = function(e) {
        errors[i] <<- conditionMessage(e)
        NULL
      }
    ))
  }

  tests <- list(anova_test = sa_feature_table(
    feats,
    c("n_used", "n_cells", "f_stat", "df1", "df2", "eta_sq", "omega_sq",
      "pval", "lower_conf", "upper_conf"),
    label,
    p_adjust = p_adjust,
    fun = function(i) {
      if (!is.na(errors[i])) {
        stop(errors[i], call. = FALSE)
      }
      fits[[i]]$model
    }
  ))

  terms_tbl <- sa_fact_term_table(feats, plan, fits, p_adjust, centers,
                                  design$cells)
  cells_tbl <- sa_fact_cell_table(feats, per_feature, fits, design)

  posthoc_tables <- list()
  n_posthoc <- 0L
  if (posthoc) {
    stage <- sa_fact_posthoc_stage(feats, fits, terms_tbl, design,
                                   posthoc_scope, posthoc_alpha, conf_level)
    posthoc_tables$anova_test <- stage$table
    n_posthoc <- stage$n_posthoc
  }

  diagnostics <- if (diagnose) {
    sa_diagnose_samples(per_feature, feats, design$cell_label, paired = FALSE)
  } else {
    NULL
  }

  sa_new_comparison(
    analysis  = "factorial_comparison",
    features  = feats,
    design    = list(
      factor_lv     = design$factor_lv,
      anova_type    = design$anova_type,
      n_factors     = length(design$factor_lv),
      group_lv      = design$cell_label,
      cell_n        = design$cell_n,
      n_empty_cells = design$n_empty_cells,
      paired        = FALSE,
      pairing       = NA_character_,
      within        = character(0),
      n_dropped     = design$n_dropped,
      unmatched_ids = character(0)
    ),
    parameters = list(
      alternative      = "two.sided",
      conf_level       = conf_level,
      ss_type          = ss_type,
      fc_mean          = fc_mean,
      input_scale      = input_scale,
      p_adjust         = p_adjust,
      posthoc          = posthoc,
      posthoc_alpha    = posthoc_alpha,
      # Every contrast is judged against the studentised range of its own block,
      # so its p-value is family-wise already and there is no adjustment for this
      # field to record.
      posthoc_p_adjust = NA_character_,
      posthoc_scope    = posthoc_scope,
      n_posthoc        = n_posthoc
    ),
    effect    = effect,
    tests     = tests,
    terms     = terms_tbl,
    cells     = cells_tbl,
    posthoc   = posthoc_tables,
    test_info = list(anova_test = list(
      id            = "factorial_anova",
      label         = label,
      paired        = FALSE,
      posthoc_id    = "factorial_tukey",
      posthoc_label = "Tukey HSD on marginal means and simple effects"
    )),
    diagnostics = diagnostics,
    subclass  = "sa_factorial"
  )
}


#' Resolve the crossed factors and lay the observations out in cells
#'
#' Everything the analysis needs to know about where an observation sits, settled
#' in one pass so that the levels, the cells and the row selections cannot be
#' derived from each other twice in different orders. The cell numbering and the
#' labels come from `utils_factorial.R`, which is what the simulator counts with.
#'
#' @param data Wide data.frame, already validated.
#' @param factors,factor_lv,control_label The arguments as received.
#'
#' @return List with `factor_lv` in declaration order, the `cells` grid,
#'   `n_cells`, `cell_label`, `cell_n`, `rows_of_cell`, `n_empty_cells`,
#'   `n_dropped` and `anova_type`.
#'
#' @keywords internal
#' @noRd
sa_fact_layout <- function(data, factors, factor_lv, control_label = NULL) {
  if (!is.list(factors) || length(factors) < 2L || is.null(names(factors)) ||
        anyNA(names(factors)) || !all(nzchar(names(factors))) ||
        anyDuplicated(names(factors)) > 0L) {
    stop("`factors` must be a named list of at least two crossed factors, each ",
         "entry a column name of `data` or one value per row of it. Use ",
         "compare_multiple_groups() for a single factor.", call. = FALSE)
  }

  # The cell table names a column after each factor and holds its level there,
  # which is what lets a reader subset on a factor without parsing a cell label.
  # A factor sharing a name with one of the statistics columns would silently
  # overwrite it, so the collision is refused here rather than discovered in the
  # result.
  taken <- intersect(names(factors), sa_cell_table_columns())
  if (length(taken) > 0L) {
    stop("`factors` may not name a factor after a column of the cell table: ",
         paste(taken, collapse = ", "), ". Reserved: ",
         paste(sa_cell_table_columns(), collapse = ", "),
         ". Rename the factor.", call. = FALSE)
  }

  values <- lapply(names(factors), function(nm) {
    # A row that does not say which cell it is in cannot be used, but it can be
    # dropped, which is what happens to a row outside `factor_lv` too.
    as.character(sa_resolve_row_vector(factors[[nm]], paste0("factors$", nm),
                                       data, allow_na = TRUE)$value)
  })
  names(values) <- names(factors)

  named_lv <- !is.null(factor_lv)
  if (!named_lv) {
    factor_lv <- lapply(values, function(v) sort(unique(v[!is.na(v)])))
  } else {
    if (!is.list(factor_lv) || is.null(names(factor_lv)) ||
          !setequal(names(factor_lv), names(factors))) {
      stop("`factor_lv` must be a named list holding the levels of every ",
           "factor `factors` names, or NULL to take them from the data. ",
           "`factors` has: ", paste(names(factors), collapse = ", "), ".",
           call. = FALSE)
    }
    # Declaration order is the term order and the Type I order, and it is
    # `factor_lv` that states it when both lists are given.
    values <- values[names(factor_lv)]
  }

  for (nm in names(factor_lv)) {
    lv <- as.character(factor_lv[[nm]])
    if (length(lv) < 2L || anyNA(lv) || !all(nzchar(lv)) ||
          anyDuplicated(lv) > 0L) {
      stop("`factor_lv$", nm, "` must be at least two distinct non-empty ",
           "level names, the first being the reference.", call. = FALSE)
    }
    absent <- lv[!lv %in% values[[nm]]]
    if (length(absent) > 0L) {
      stop("`factor_lv$", nm, "` level(s) absent from `factors$", nm, "`: ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    factor_lv[[nm]] <- lv
  }

  # After the levels are settled and before anything is counted from them, so
  # that the grid, the cell labels and the row selections are all built from the
  # order the reference ended up in rather than corrected afterwards.
  factor_lv <- sa_fact_control_first(factor_lv, control_label,
                                     if (named_lv) "factor_lv" else "factors")

  level_idx <- matrix(0L, nrow = nrow(data), ncol = length(factor_lv),
                      dimnames = list(NULL, names(factor_lv)))
  for (nm in names(factor_lv)) {
    level_idx[, nm] <- match(values[[nm]], factor_lv[[nm]])
  }
  # `match()` leaves NA for a value the levels do not hold, so a row naming a
  # level outside `factor_lv` and a row naming none are dropped together.
  keep <- stats::complete.cases(level_idx)

  cells <- sa_fact_grid(factor_lv)
  n_cells <- nrow(cells)
  cell_idx <- sa_fact_cell_index(level_idx[keep, , drop = FALSE],
                                 vapply(factor_lv, length, integer(1)))
  rows_of_cell <- split(which(keep), factor(cell_idx, levels = seq_len(n_cells)))

  list(
    factor_lv     = factor_lv,
    cells         = cells,
    n_cells       = n_cells,
    cell_label    = sa_fact_cell_labels(factor_lv, cells),
    cell_n        = lengths(rows_of_cell, use.names = FALSE),
    rows_of_cell  = unname(rows_of_cell),
    n_empty_cells = sum(lengths(rows_of_cell) == 0L),
    n_dropped     = sum(!keep),
    anova_type    = sa_fact_anova_type(length(factor_lv))
  )
}


#' Which ANOVA the number of factors makes this
#'
#' Two factors are a two-way ANOVA and three are a three-way one, and past that
#' the ordinal stops being said out loud. All three are the same fully crossed
#' model, so this names what happened rather than choosing it.
#'
#' @keywords internal
#' @noRd
sa_fact_anova_type <- function(n_factors) {
  switch(as.character(n_factors),
         "2" = "two_way",
         "3" = "three_way",
         "factorial")
}


#' The readable name of the analysis that ran
#'
#' @keywords internal
#' @noRd
sa_fact_anova_label <- function(anova_type, ss_type) {
  paste0(switch(anova_type,
                two_way   = "Two-way ANOVA",
                three_way = "Three-way ANOVA",
                "Factorial ANOVA"),
         " (Type ", ss_type, " sums of squares)")
}


#' Fold change of the most extreme cell against the reference cell
#'
#' `sa_multi_fold_change()` with the cells as its levels. The quantity is the same
#' one a multi-group comparison reports, so it is computed by the same code and
#' only the two columns that name what a level is get renamed: a cell is a
#' combination of levels rather than a level.
#'
#' @keywords internal
#' @noRd
sa_fact_effect <- function(centers, feats, cell_label, mean_type) {
  out <- sa_multi_fold_change(centers, feats, cell_label, mean_type)
  names(out)[names(out) == "n_groups"] <- "n_cells"
  names(out)[names(out) == "extreme_level"] <- "extreme_cell"
  out[c("features", "n_used", "n_cells", "ref_center", "extreme_cell",
        "extreme_center", "fold_change", "log2fc")]
}


#' Stack the cell means into one feature by cell table
#'
#' The one thing the rest of the result cannot be read backwards into. `$effect`
#' reduces the grid to the reference cell and the most extreme one, and a term
#' test reduces it to a p-value; neither says which way the lines of an
#' interaction run, which is what [draw_interaction_plot()] needs.
#'
#' The means are recomputed from the samples rather than taken from `fits`, so
#' that a feature whose model could not be fitted still reports what its cells
#' held. They are the same numbers either way: `sa_factorial_anova()` takes the
#' arithmetic mean of these very vectors. `se` is the one column that does need
#' the fit, because it is pooled over the whole model, and it is `NA` for a
#' feature that has no fit to pool over.
#'
#' @param feats Feature names, fixing the block order.
#' @param per_feature One list of cell samples per feature, in cell order, with
#'   the missing values already dropped.
#' @param fits One `sa_factorial_anova()` result per feature, `NULL` where it
#'   failed.
#' @param design The layout `sa_fact_layout()` returned.
#'
#' @keywords internal
#' @noRd
sa_fact_cell_table <- function(feats, per_feature, fits, design) {
  n_cells <- design$n_cells
  n_feats <- length(feats)

  out <- data.frame(features = rep(feats, each = n_cells),
                    stringsAsFactors = FALSE)
  # One column per factor, holding the level name rather than its index, in
  # declaration order. `design$cells` counts the first factor fastest, which is
  # the order `cell_label` and `cell_n` are already in.
  for (nm in names(design$factor_lv)) {
    out[[nm]] <- rep(design$factor_lv[[nm]][design$cells[[nm]]],
                     times = n_feats)
  }
  out$cell <- rep(design$cell_label, times = n_feats)

  # Counted per feature rather than taken from `design$cell_n`, which counts
  # rows: a cell of ten rows holding three missing values contributes seven
  # observations to this feature and ten to one with no gaps.
  out$n <- as.integer(unlist(lapply(per_feature, lengths), use.names = FALSE))
  out$mean <- unlist(lapply(per_feature, function(cells) {
    vapply(cells, function(x) if (length(x) == 0L) NA_real_ else mean(x),
           numeric(1))
  }), use.names = FALSE)
  out$sd <- unlist(lapply(per_feature, function(cells) {
    vapply(cells, function(x) if (length(x) < 2L) NA_real_ else stats::sd(x),
           numeric(1))
  }), use.names = FALSE)
  out$se <- unlist(lapply(seq_len(n_feats), function(i) {
    if (is.null(fits[[i]])) {
      rep(NA_real_, n_cells)
    } else {
      sqrt(fits[[i]]$ms_error / lengths(per_feature[[i]]))
    }
  }), use.names = FALSE)

  rownames(out) <- NULL
  out
}


#' Stack the per-feature term results into one feature by term table
#'
#' The counterpart of `sa_feature_table()` for the term axis. A feature whose
#' model could not be fitted is present with `NA` in every term, since the term
#' rows of a feature exist as questions whether or not they were answerable, and
#' the reason it failed is reported once by the whole-model table rather than
#' twice.
#'
#' @param feats Feature names, fixing the block order.
#' @param plan The list `sa_factorial_plan()` returns, holding the term labels.
#' @param fits One `sa_factorial_anova()` result per feature, `NULL` where it
#'   failed.
#' @param p_adjust Method passed to [stats::p.adjust()], applied across features
#'   within each term.
#' @param centers The list `sa_group_centers()` returns, whose matrix the term
#'   effect sizes are decomposed out of.
#' @param cells Grid of level indices, as `sa_fact_grid()` returns.
#'
#' @keywords internal
#' @noRd
sa_fact_term_table <- function(feats, plan, fits, p_adjust, centers, cells) {
  columns <- c("n_used", "df", "ss", "ms", "f_stat", "df_error", "eta_sq",
               "partial_eta_sq", "pval")
  n_terms <- length(plan$labels)

  blocks <- lapply(fits, function(fit) {
    if (is.null(fit)) {
      matrix(NA_real_, nrow = n_terms, ncol = length(columns),
             dimnames = list(NULL, columns))
    } else {
      fit$terms[, columns, drop = FALSE]
    }
  })

  out <- data.frame(
    features   = rep(feats, each = n_terms),
    terms      = rep(plan$labels, times = length(feats)),
    term_order = rep(plan$orders, times = length(feats)),
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(do.call(rbind, blocks)))
  rownames(out) <- NULL

  # The effect size of a term, on the log2 scale the centres put it on. Taken
  # from the same centre matrix `effect` is, so the two are the same numbers read
  # at two grains, and a negative or zero centre leaves it NaN or -Inf there as
  # it does there. Components are deviations from what the other terms predict,
  # so they do not depend on which cell is the reference the way `effect$log2fc`
  # does: subtracting the reference from every cell cancels.
  # One column per feature, so reading it down its columns is the feature-major
  # order the rows are already in.
  log2_center <- suppressWarnings(log2(centers$centers))
  out$log2_effect <- as.vector(vapply(seq_along(feats), function(i) {
    sa_fact_term_effect(log2_center[i, ], cells, plan$terms)
  }, numeric(n_terms)))

  # One term is one family. The same question asked of every feature is what a
  # feature-axis adjustment is for; the main effect and the interaction are two
  # different questions, and pooling them would correct each for the other's
  # multiplicity.
  out$pval_adj <- NA_real_
  for (nm in plan$labels) {
    at <- out$terms == nm
    out$pval_adj[at] <- stats::p.adjust(out$pval[at], method = p_adjust)
  }

  out[c("features", "terms", "term_order", "n_used", "df", "ss", "ms", "f_stat",
        "df_error", "eta_sq", "partial_eta_sq", "log2_effect", "pval",
        "pval_adj")]
}


#' Run the contrasts each feature's term tests earned it
#'
#' The pairwise stage of a factorial model is ragged on both axes: a feature has
#' contrasts only for the terms it cleared, so the marginal comparisons of one
#' factor can be present while another factor's are absent for the same feature.
#' That is why this is assembled here rather than by `sa_posthoc_table()`, whose
#' table is one fixed set of pairs per qualifying feature.
#'
#' @param feats Feature names.
#' @param fits One `sa_factorial_anova()` result per feature, `NULL` where it
#'   failed.
#' @param terms_tbl The term table, whose `pval_adj` decides what runs.
#' @param design The layout `sa_fact_layout()` returned.
#' @param scope `"both"`, `"marginal"` or `"simple"`.
#' @param alpha Threshold on the gating term's `pval_adj`.
#' @param conf_level Confidence level for the intervals.
#'
#' @return List with `table`, holding the post-hoc contract columns plus `factor`
#'   and `stratum`, and `n_posthoc`, the number of features that contributed a
#'   row.
#'
#' @keywords internal
#' @noRd
sa_fact_posthoc_stage <- function(feats, fits, terms_tbl, design, scope, alpha,
                                  conf_level) {
  factor_lv <- design$factor_lv
  skeleton <- sa_fact_contrast_skeleton(design)
  rows <- skeleton$table
  columns <- c("features", "factor", "stratum",
               setdiff(sa_posthoc_table_columns(), "features"))

  wanted <- switch(scope,
                   both     = rep(TRUE, nrow(rows)),
                   marginal = is.na(rows$stratum),
                   simple   = !is.na(rows$stratum))
  candidates <- which(wanted)

  # The term that licenses a contrast: a factor's own main effect for a marginal
  # comparison, and its interaction with everything held fixed for a simple one.
  # In `factor_lv` order, which is the order the term labels are written in.
  gate <- vapply(candidates, function(k) {
    used <- if (is.na(rows$stratum[k])) {
      rows$factor[k]
    } else {
      names(factor_lv)
    }
    paste(names(factor_lv)[names(factor_lv) %in% used], collapse = ":")
  }, character(1))
  nmeans <- vapply(rows$factor, function(f) length(factor_lv[[f]]), integer(1))

  blocks <- list()
  failures <- character(0)
  for (i in seq_along(feats)) {
    if (is.null(fits[[i]])) {
      next
    }
    padj <- terms_tbl$pval_adj[terms_tbl$features == feats[i]]
    names(padj) <- terms_tbl$terms[terms_tbl$features == feats[i]]
    earned <- padj[gate]
    take <- candidates[!is.na(earned) & earned <= alpha]
    if (length(take) == 0L) {
      next
    }
    mat <- tryCatch(
      sa_factorial_tukey(fits[[i]], skeleton, nmeans, take, conf_level),
      error = function(e) {
        failures[[feats[i]]] <<- conditionMessage(e)
        matrix(NA_real_, nrow = length(take),
               ncol = length(sa_posthoc_columns()),
               dimnames = list(NULL, sa_posthoc_columns()))
      }
    )
    stats_df <- as.data.frame(mat)
    # The studentised range controls the error rate over the block of contrasts
    # this one sits in, so adjusting it again would correct twice for one
    # comparison.
    stats_df$pval_adj <- stats_df$pval
    blocks[[length(blocks) + 1L]] <- data.frame(
      features = rep(feats[i], length(take)),
      rows[take, c("factor", "stratum", "contrast", "group1", "group2"),
           drop = FALSE],
      stats_df,
      stringsAsFactors = FALSE
    )
  }

  if (length(failures) > 0L) {
    warning("Tukey HSD on marginal means and simple effects could not be ",
            "computed for ", length(failures), " of ", length(feats),
            " feature(s); those rows are NA:\n",
            paste0("  ", names(failures), ": ", failures, collapse = "\n"),
            call. = FALSE)
  }

  if (length(blocks) == 0L) {
    empty <- data.frame(features = character(0), factor = character(0),
                        stratum = character(0), contrast = character(0),
                        group1 = character(0), group2 = character(0),
                        stringsAsFactors = FALSE)
    for (nm in setdiff(columns, names(empty))) {
      empty[[nm]] <- numeric(0)
    }
    return(list(table = empty[columns], n_posthoc = 0L))
  }

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  list(table = out[columns], n_posthoc = length(blocks))
}
