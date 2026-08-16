#' Simulate a crossed-factor experiment whose answer is known
#'
#' The factorial counterpart of [simulate_multiple_groups()]. Crosses any number
#' of factors, lets each one be measured between subjects or within them, and
#' returns the planted answer alongside the data so that a two-way or an n-way
#' analysis can be scored against what was actually put there.
#'
#' One factor asks one question: are the levels alike. Crossing a second one asks
#' three, and they fail separately and for different reasons. Each factor has a
#' main effect, the pair has an interaction, and a design that is read as though
#' the second factor were not there answers none of them. So the effect is
#' planted in five shapes rather than one, chosen to make the three questions come
#' apart: a shape whose main effects are real and whose interaction is not, and a
#' shape whose interaction is real and whose main effects are exactly zero, are
#' both here, and no single test tells them apart.
#'
#' @param n_feats Number of features to generate. Columns are named `prot_1`
#'   upwards, or whatever `feat_prefix` asks for.
#' @param factor_lv Named list of factors, each entry the levels of one factor
#'   with the reference level first. Its length is how many factors are crossed,
#'   so there is no separate argument for that, and there have to be at least
#'   two: one factor is [simulate_multiple_groups()].
#'
#'   The **first factor is the primary one**, the treatment the experiment is
#'   about, and the ones after it are the other factors the effect may or may not
#'   depend on. The shapes below are written in those terms. The cell in which
#'   every factor sits at its reference level is the reference cell, and it is
#'   what the planted effects are measured against.
#' @param within Names of the factors measured within subjects, or `NULL` for a
#'   design that is entirely between them. A subject belongs to one combination
#'   of the between factors and is measured under every combination of the within
#'   ones, so naming no factor gives a factorial ANOVA, naming all of them gives
#'   a fully repeated design, and naming some of them gives a mixed one. When any
#'   factor is named, `args` gains `id` and `within`.
#' @param n_per_cell Observations per cell. When there are within factors this is
#'   also the number of subjects per combination of the between factors, because
#'   each of those subjects contributes exactly one observation to each cell it
#'   is measured in. One number spreads over every cell; a vector carries one
#'   size per combination of the between factors, which makes its length how many
#'   of those there are, the rule `n_treat` follows in
#'   [simulate_multiple_groups()]. The within factors are fully crossed with
#'   every subject, so they cannot hold different sizes.
#' @param n_up,n_down How many features are moved up and down. Their sum cannot
#'   exceed `n_feats`, and every other feature is left with a true effect of
#'   exactly zero in every cell and every term. The defaults are a fraction of
#'   `n_feats` rather than a fixed count, so that asking for fewer features
#'   plants fewer effects instead of failing.
#' @param term_mix Named vector of relative weights over the five shapes of
#'   effect described under "The five shapes" below, which decide **which terms**
#'   of the model an effect is planted in. Set a weight to zero to leave that
#'   shape out. The planted features are split between the shapes by the largest
#'   remainder method rather than drawn at random, so the counts are exactly what
#'   the weights ask for and do not move with the seed.
#' @param pattern_mix Named vector of relative weights over `"all"`,
#'   `"gradient"` and `"single"`, which decide **how a factor spreads its main
#'   effect** over its non-reference levels, exactly as in
#'   [simulate_multiple_groups()]. The two mixes are on different axes and are
#'   crossed at random: `term_mix` says which terms move and `pattern_mix` says
#'   what the profile along a factor looks like. Split by largest remainder too,
#'   so both sets of counts are a function of the arguments alone.
#' @param expr_range Range the baseline log2 abundance of each feature is drawn
#'   from. Every cell shares the baseline, which is what makes an unplanted
#'   feature null in every term at once.
#' @param ref_sd,cell_sd Ranges the per-feature standard deviation of the
#'   reference cell and of every other cell are drawn from. Every cell draws its
#'   own, so the design is heteroscedastic and unbalanced variance is something
#'   the analysis has to survive rather than something it is spared. Pass the
#'   same range twice for equal variances. The defaults are the two ranges
#'   [simulate_multiple_groups()] uses.
#'
#'   Keeping them apart costs something that is worth knowing about. An `aov()`
#'   interaction test on a `"main_only"` feature, whose interaction is exactly
#'   zero, is called about eight times in a hundred rather than five, because the
#'   cells whose means differ are also the cells whose spread differs. Passing
#'   the same range twice brings it back to about three. The anticonservatism is
#'   the test's rather than the simulation's, and finding it is the sort of thing
#'   a simulation with a known answer is for.
#' @param deg_log2fc Range the magnitude of the planted effect is drawn from, on
#'   the log2 scale. One magnitude is drawn per planted feature and the shape
#'   decides how it is spread over the terms.
#' @param interaction_scale Size of the interaction relative to the main effect
#'   under the `"interaction"` shape, as a fraction of the magnitude drawn from
#'   `deg_log2fc`. It has no bearing on `"crossover"`, where the interaction
#'   carries the whole effect because there is nothing else to carry it. Must be
#'   above zero: a scale of zero would leave a feature whose `pattern` says
#'   `"interaction"` with no interaction in it.
#'
#'   The default puts the `"interaction"` shape between the other two things an
#'   interaction row can be. An `aov()` on the defaults finds its interaction
#'   about three times in five, against about nine in ten for the `"crossover"`
#'   shape and about one in twenty where no interaction was planted, so the row
#'   is neither trivially recovered nor indistinguishable from noise. Raising it
#'   to 1 takes the shape to nine in ten as well and the two stop differing.
#' @param subject_sd Range the per-feature subject standard deviation is drawn
#'   from. A subject's offset is drawn once per feature and reused across every
#'   condition it is measured under, which is what a within-subject test exists
#'   to remove. Ignored when `within` is `NULL`.
#' @param feat_prefix Prefix for the generated feature names. `"prot"` gives
#'   `prot_1`, `prot_2` and so on.
#' @param seed Seed for the draw, or `NULL` to use the stream as it stands.
#'   Supplying one does not disturb the caller: the previous random number state
#'   is put back when the function returns.
#'
#' @return A list of five elements.
#'
#'   \describe{
#'     \item{`args`}{`data`, `feats`, `factors`, `factor_lv` and `input_scale`,
#'       and under a within design also `within` and `id`. `factors` is a named
#'       list holding one vector per factor, each as long as `data` has rows, and
#'       `factor_lv` is the level order of each, exactly as it was passed in.
#'       These are the names the factorial comparison will take, so it will be
#'       one `do.call()` away when it exists.}
#'     \item{`truth`}{One row per feature, aligned with `feats`, holding
#'       `features`, `pattern`, `spread`, `direction`, `partner`,
#'       `extreme_cell`, `extreme_tied`, `log2fc`, `baseline` and `sd_subject`.}
#'     \item{`truth_term`}{One row per feature and model term, every main effect
#'       and every interaction of every order, holding `features`, `terms`,
#'       `term_order`, `is_within`, `max_abs_delta` and `is_effect`. This is the
#'       table that scores an ANOVA table row by row, and the one that has no
#'       counterpart in [simulate_multiple_groups()].}
#'     \item{`truth_cell`}{One row per feature and cell, holding `features`, one
#'       column per factor, then `is_ref`, `delta`, `center`, `sd` and `n`. A
#'       feature the analysis missed can be looked up here rather than guessed
#'       at: a large `sd` explains a miss that the effect size alone does not.}
#'     \item{`truth_contrast`}{One row per feature and pair of levels, in the row
#'       order and direction a post-hoc table uses, holding `features`,
#'       `factor`, `stratum`, `contrast`, `group1`, `group2`, `delta` and
#'       `is_diff`. A `stratum` of `NA` is the marginal contrast, averaged over
#'       the other factors; anything else names the combination of the other
#'       factors the contrast was taken inside, which is the simple effect.}
#'   }
#'
#' @section The five shapes:
#' Each planted feature is given a magnitude `d` drawn from `deg_log2fc`,
#' positive for an up feature and negative for a down one, a shape from
#' `term_mix`, and for every shape but the first a partner factor drawn at random
#' from the factors after the primary one. The shape decides which terms of the
#' model end up carrying `d`.
#'
#' \describe{
#'   \item{`"main_only"`}{The primary factor moves and nothing else does. Every
#'     other main effect and every interaction is exactly zero. This is
#'     [simulate_multiple_groups()] inside a factorial frame, and the case a
#'     two-way analysis should answer with one row.}
#'   \item{`"additive"`}{The primary factor and the partner each move, and their
#'     interaction is exactly zero. The cell means are the sum of the two, so the
#'     profiles are parallel, and an interaction reported as significant here is
#'     a false positive by construction.}
#'   \item{`"interaction"`}{The primary factor moves and the size of its effect
#'     depends on the level of the partner. The primary main effect and the
#'     interaction are both real; the partner's own main effect is left at
#'     exactly zero, so the two terms that should be called are the only two
#'     there are.}
#'   \item{`"crossover"`}{Pure interaction. The primary factor rises at the
#'     partner's reference level and falls at the others, by amounts arranged so
#'     that **both main effects are exactly zero** while the cells plainly
#'     differ. This is the shape a main-effect test has to miss and an
#'     interaction test has to catch, which is the reason a factorial design is
#'     analysed as one.}
#'   \item{`"nuisance_only"`}{The partner moves and the primary factor is exactly
#'     zero. Read as one factor, the primary factor looks null with inflated
#'     within-group spread; read as a factorial design, the spread is accounted
#'     for and belongs to a term of its own.}
#' }
#'
#' A feature that was not planted has a delta of exactly zero in every cell and a
#' component of exactly zero in every term. Both kinds of mistake are therefore
#' defined for every row of `truth_term`: a term called significant on a zero
#' component is a false positive, and a non-zero component that was not called is
#' a miss.
#'
#' @section What the defaults leave recoverable:
#' The defaults are set so that an analysis gets most of the answer and not all
#' of it, which is the band between a simulation that is trivially recovered and
#' one that looks broken. An `aov(y ~ treatment * sex)` on the default 4 by 2
#' design finds the treatment main effect about four times in five, which is the
#' rate [simulate_multiple_groups()] was tuned to, the two-level factor's main
#' effect rather more often than that, and the interaction between the two,
#' depending on which shape planted it. Every term that was not planted is called
#' at about a twentieth, which is what makes a false positive rate readable off
#' this simulation at all.
#'
#' The rates are a function of `n_per_cell` as much as of the effect sizes, and
#' they fall away quickly below the default: at eight per cell the
#' `"interaction"` shape's interaction is no longer distinguishable from a null
#' term. A design small enough to be quick is not a design these numbers describe.
#'
#' @section How the effect is planted, and how it is reported:
#' The effect is built in the space an ANOVA decomposes into: a main effect is a
#' profile along its own factor that sums to zero, and an interaction sums to zero
#' along each of its factors. The components are added up and the value at the
#' reference cell is then subtracted from the whole array, which leaves the
#' reference cell at exactly zero delta without touching any term, since a
#' constant belongs to the grand mean alone.
#'
#' That is why the two tables read differently and both are right.
#' `truth_cell$delta` is the shift of a cell from the reference cell, the way
#' `truth_group$delta` is in [simulate_multiple_groups()], and for a
#' `"main_only"` feature the cell at level `j` of the primary factor carries
#' exactly what the `pattern_mix` profile put there. `truth_term$max_abs_delta`
#' is the largest component of the ANOVA effect itself, which is the quantity
#' that is exactly zero for a term that was not planted. Components smaller than
#' `1e-8` in absolute value are recorded as exactly zero: they are the rounding
#' left over from averaging, and a term left holding `3e-17` would score every
#' row of an ANOVA table against the wrong answer.
#'
#' @section Directions:
#' `direction` is the sign of `d`, and it is the sign of the primary factor's
#' effect at the reference level of the partner. `truth$log2fc` is the delta of
#' whichever cell sits furthest from the reference cell, so for every shape but
#' `"crossover"` an up feature is positive there. Under `"crossover"` the primary
#' factor moves in opposite directions at different levels of the partner, so
#' which cell is furthest, and its sign, follow from the shape rather than from
#' `direction`.
#'
#' `truth_contrast$delta` is `group1 - group2` with `group1` the later level of
#' the factor, which is the direction and the row order a post-hoc table uses. It
#' comes from `sa_level_pairs()`, the same helper the post-hoc tables are built
#' from, so the two cannot drift apart.
#'
#' `extreme_cell` is the levels of that cell joined by a dot, so it reads back
#' against `truth_cell` without a lookup. When more than one cell is equally far
#' from the reference, which is what the `"all"` profile does on purpose, it
#' records the first of them and `extreme_tied` is `TRUE`: the flag that says to
#' score the magnitude rather than the name of the cell. It is `TRUE` with
#' `extreme_cell` missing for an unplanted feature, whose cells are all zero and
#' none of which is furthest.
#'
#' @section Within and between:
#' A subject belongs to one combination of the between factors and is measured
#' under every combination of the within ones, so no subject is dropped and the
#' within-subject rectangle is complete. Each subject gets an offset per feature,
#' drawn once and added to all of its rows, which is the between-subject
#' variation a within-subject test removes. The residual standard deviation still
#' differs from cell to cell, so sphericity does not hold and the corrections a
#' repeated measures analysis reports have something to report.
#'
#' @seealso [simulate_multiple_groups()] for the one-factor case, whose
#'   `pattern_mix` shapes this reuses, and [simulate_two_groups()] for two
#'   groups.
#'
#' @examples
#' ## A 4 x 2 design: four treatments crossed with sex, both between subjects.
#' sim <- simulate_factorial_groups(n_feats = 20, n_up = 5, n_down = 5,
#'                                  n_per_cell = 6, seed = 1)
#' table(pattern = sim$truth$pattern, direction = sim$truth$direction)
#'
#' ## The answer per term. A "crossover" feature has an interaction and no main
#' ## effect at all, which is what makes the shape worth planting.
#' cross <- sim$truth$features[sim$truth$pattern == "crossover"][1]
#' subset(sim$truth_term, features == cross,
#'        select = c("terms", "max_abs_delta", "is_effect"))
#'
#' ## Scored against a two-way ANOVA. The term names line up with `truth_term`,
#' ## so the two tables merge without being renamed.
#' long <- data.frame(y = sim$args$data[[cross]], sim$args$factors)
#' summary(stats::aov(y ~ treatment * sex, data = long))
#'
#' ## An "additive" feature has both main effects and no interaction, so it is
#' ## the shape that makes an interaction call a false positive.
#' add <- sim$truth$features[sim$truth$pattern == "additive"][1]
#' subset(sim$truth_term, features == add, select = c("terms", "is_effect"))
#'
#' ## A missed feature is looked up rather than guessed at. The row per cell
#' ## carries the spread the cell was given, which is one of the reasons.
#' head(subset(sim$truth_cell, features == cross), 4)
#'
#' ## Marginal contrasts and simple effects are both in `truth_contrast`. A
#' ## stratum of NA is the marginal one.
#' head(subset(sim$truth_contrast, features == cross & factor == "treatment"), 3)
#'
#' ## Time measured on the same subjects, treatment and sex between them: a
#' ## mixed design. `args` then carries `id` and `within`.
#' mixed <- simulate_factorial_groups(
#'   n_feats = 10, n_up = 2, n_down = 2, n_per_cell = 4,
#'   factor_lv = list(treatment = c("control", "treat_A", "treat_B"),
#'                    sex       = c("male", "female"),
#'                    time      = c("T0", "T1", "T2")),
#'   within = "time", seed = 1
#' )
#' names(mixed$args)
#' ## Every subject under every time point, so the rectangle is complete.
#' table(table(mixed$args$id))
#'
#' ## Which terms are tested in the within-subject error stratum is recorded, so
#' ## the answer table knows which half of a mixed ANOVA to score.
#' unique(mixed$truth_term[c("terms", "term_order", "is_within")])
#'
#' @export
simulate_factorial_groups <- function(n_feats = 100,
                                      factor_lv = list(
                                        treatment = c("control", "treat_A",
                                                      "treat_B", "treat_C"),
                                        sex       = c("male", "female")
                                      ),
                                      within = NULL,
                                      n_per_cell = 20,
                                      n_up = round(0.15 * n_feats),
                                      n_down = round(0.15 * n_feats),
                                      term_mix = c(main_only = 1, additive = 1,
                                                   interaction = 1,
                                                   crossover = 1,
                                                   nuisance_only = 1),
                                      pattern_mix = c(all = 1, gradient = 1,
                                                      single = 1),
                                      expr_range = c(2, 12),
                                      ref_sd = c(1.2, 2.4),
                                      cell_sd = c(1.8, 3.2),
                                      deg_log2fc = c(1, 2.5),
                                      interaction_scale = 0.8,
                                      subject_sd = c(2, 4),
                                      feat_prefix = "prot",
                                      seed = NULL) {

  design <- sa_fact_design(factor_lv, within, n_per_cell)
  factor_lv <- design$factor_lv
  fac_names <- names(factor_lv)
  n_cells <- design$n_cells
  has_within <- length(design$within) > 0L

  # `n_up` and `n_down` default to a fraction of `n_feats`, so their promises are
  # forced after this line and see the checked value rather than whatever was
  # passed in.
  n_feats <- sa_check_count(n_feats, "n_feats", 1)
  n_up <- sa_check_count(n_up, "n_up")
  n_down <- sa_check_count(n_down, "n_down")
  if (n_up + n_down > n_feats) {
    stop("`n_up` + `n_down` is ", n_up + n_down, ", which is more features ",
         "than the ", n_feats, " that `n_feats` asks for.", call. = FALSE)
  }
  shapes <- sa_sim_pattern_mix(term_mix, sa_fact_shapes(), "term_mix")
  mix <- sa_sim_pattern_mix(pattern_mix)
  sa_check_range(expr_range, "expr_range")
  sa_check_range(ref_sd, "ref_sd", 0)
  sa_check_range(cell_sd, "cell_sd", 0)
  sa_check_range(deg_log2fc, "deg_log2fc", 0)
  sa_check_scalar_num(interaction_scale, "interaction_scale", 0,
                      lower_open = TRUE)
  sa_check_range(subject_sd, "subject_sd", 0)
  if (!is.character(feat_prefix) || length(feat_prefix) != 1L ||
      is.na(feat_prefix) || !nzchar(feat_prefix)) {
    stop("`feat_prefix` must be a single non-empty string.", call. = FALSE)
  }

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  feats <- paste0(feat_prefix, "_", seq_len(n_feats))
  baseline <- stats::runif(n_feats, expr_range[1], expr_range[2])

  sd_mat <- matrix(0, nrow = n_feats, ncol = n_cells,
                   dimnames = list(feats, design$cell_label))
  sd_mat[, design$ref_cell] <- stats::runif(n_feats, ref_sd[1], ref_sd[2])
  for (j in setdiff(seq_len(n_cells), design$ref_cell)) {
    sd_mat[, j] <- stats::runif(n_feats, cell_sd[1], cell_sd[2])
  }
  sd_subject <- if (has_within) {
    stats::runif(n_feats, subject_sd[1], subject_sd[2])
  } else {
    rep(NA_real_, n_feats)
  }

  delta <- matrix(0, nrow = n_feats, ncol = n_cells,
                  dimnames = list(feats, design$cell_label))
  direction <- rep("none", n_feats)
  pattern <- rep("none", n_feats)
  spread <- rep("none", n_feats)
  partner <- rep(NA_character_, n_feats)

  if (n_up + n_down > 0L) {
    # Taken from the head and the tail of one shuffled draw. Selecting the down
    # set as the complement of the up set with `-idx` would return everything
    # rather than nothing when the up set is empty.
    picked <- sample.int(n_feats, n_up + n_down)
    up_idx <- utils::head(picked, n_up)
    down_idx <- utils::tail(picked, n_down)
    direction[up_idx] <- "up"
    direction[down_idx] <- "down"

    plant_idx <- c(up_idx, down_idx)
    plant_mag <- c(stats::runif(n_up, deg_log2fc[1], deg_log2fc[2]),
                   -stats::runif(n_down, deg_log2fc[1], deg_log2fc[2]))
    # Each direction is split between the shapes on its own, so a mix holds
    # within the up set and within the down set rather than only in total.
    plant_shape <- c(rep(names(shapes), sa_sim_allocate(n_up, shapes)),
                     rep(names(shapes), sa_sim_allocate(n_down, shapes)))
    # The two mixes are handed out in blocks over the same features, so without
    # a shuffle the term shape and the level profile would arrive in lockstep
    # and their crossing would never be covered. Shuffling the profiles leaves
    # both sets of counts exact and makes only the pairing random.
    plant_spread <- c(sa_fact_shuffle(rep(names(mix),
                                          sa_sim_allocate(n_up, mix))),
                      sa_fact_shuffle(rep(names(mix),
                                          sa_sim_allocate(n_down, mix))))

    for (k in seq_along(plant_idx)) {
      i <- plant_idx[k]
      mate <- sa_fact_partner(plant_shape[k], fac_names)
      eff <- sa_fact_plant(plant_mag[k], plant_shape[k], plant_spread[k],
                           mate, factor_lv, design$cells, interaction_scale)
      # The reference cell is put at exactly zero. A constant shift belongs to
      # the grand mean, so no term of the decomposition moves with it.
      delta[i, ] <- eff - eff[design$ref_cell]
      pattern[i] <- plant_shape[k]
      spread[i] <- plant_spread[k]
      partner[i] <- mate
    }
  }

  center <- baseline + delta
  # One offset per subject and feature, drawn before the rows and added to all
  # of the ones that subject owns. Drawing it per row would make it noise rather
  # than a subject effect, and the within-subject tests would have nothing to
  # gain over the independent ones.
  offsets <- if (has_within) {
    vapply(seq_len(n_feats), function(i) {
      stats::rnorm(design$n_units, 0, sd_subject[i])
    }, numeric(design$n_units))
  } else {
    NULL
  }

  cell_idx <- design$cell_idx
  values <- vapply(seq_len(n_feats), function(i) {
    v <- stats::rnorm(design$n_rows, mean = center[i, cell_idx],
                      sd = sd_mat[i, cell_idx])
    if (has_within) v + offsets[design$subject_idx, i] else v
  }, numeric(design$n_rows))

  data <- as.data.frame(values)
  names(data) <- feats
  rownames(data) <- NULL

  args <- list(
    data      = data,
    feats     = feats,
    factors   = design$factors,
    factor_lv = factor_lv
  )
  if (has_within) {
    args$within <- design$within
    args$id <- design$subject
  }
  args$input_scale <- "log2"

  list(
    args           = args,
    truth          = sa_fact_truth(feats, delta, design, pattern, spread,
                                   direction, partner, baseline, sd_subject),
    truth_term     = sa_fact_truth_term(feats, delta, design,
                                        pattern != "none"),
    truth_cell     = sa_fact_truth_cell(feats, delta, center, sd_mat, design),
    truth_contrast = sa_fact_truth_contrast(feats, delta, design)
  )
}


#' The shapes an effect can take across the terms of a factorial model
#'
#' Named in one place so that the weights argument, the validation and the
#' dispatch in `sa_fact_plant()` cannot come to disagree about what exists.
#'
#' @keywords internal
#' @noRd
sa_fact_shapes <- function() {
  c("main_only", "additive", "interaction", "crossover", "nuisance_only")
}


#' Permute a vector, leaving a vector of one alone
#'
#' `sample()` of a length-one vector reads its argument as a count, so the guard
#' is not decoration.
#'
#' @keywords internal
#' @noRd
sa_fact_shuffle <- function(x) {
  if (length(x) > 1L) x[sample.int(length(x))] else x
}


#' Work out the cells, the rows and the subjects the design implies
#'
#' The factors, which of them are repeated and how big a cell is all come from
#' three arguments that have to agree, so they are settled together rather than
#' in passes that could disagree. `factor_lv` names the factors and their levels,
#' which makes its length the number of crossed factors; `within` names a subset
#' of them; and `n_per_cell` carries one size per combination of the factors that
#' are left, which makes its length say how many of those there are.
#'
#' Nothing is drawn here. What comes back is the grid of cells, the row each
#' observation sits in and the subject that owns it.
#'
#' @param factor_lv,within,n_per_cell The arguments as received.
#'
#' @return List with the checked `factor_lv`, `within` and `between` in
#'   `factor_lv` order, the `cells` grid of level indices, `n_cells`,
#'   `cell_label`, `ref_cell`, `cell_n`, the model `terms`, the per-row
#'   `cell_idx`, `factors`, `subject`, `subject_idx`, `n_units` and `n_rows`.
#'
#' @keywords internal
#' @noRd
sa_fact_design <- function(factor_lv, within, n_per_cell) {
  if (!is.list(factor_lv) || length(factor_lv) < 2L ||
      is.null(names(factor_lv)) || anyNA(names(factor_lv)) ||
      !all(nzchar(names(factor_lv))) ||
      anyDuplicated(names(factor_lv)) > 0L) {
    stop("`factor_lv` must be a named list of at least two crossed factors, ",
         "each entry the levels of one factor with the reference level first. ",
         "Use simulate_multiple_groups() for a single factor.", call. = FALSE)
  }
  # `truth_cell` gives each factor a column of its own beside these, so a factor
  # named after one of them would produce a table with two columns of one name.
  reserved <- intersect(names(factor_lv),
                        c("features", "is_ref", "delta", "center", "sd", "n"))
  if (length(reserved) > 0L) {
    stop("`factor_lv` names factor(s) that the answer tables already use as ",
         "columns: ", paste(reserved, collapse = ", "), ".", call. = FALSE)
  }
  for (nm in names(factor_lv)) {
    lv <- factor_lv[[nm]]
    if (!is.character(lv) || length(lv) < 2L || anyNA(lv) ||
        !all(nzchar(lv)) || anyDuplicated(lv) > 0L) {
      stop("`factor_lv$", nm, "` must be at least two distinct non-empty ",
           "level names, the first being the reference.", call. = FALSE)
    }
  }

  if (is.null(within)) {
    within <- character(0)
  }
  if (!is.character(within) || anyNA(within) || anyDuplicated(within) > 0L) {
    stop("`within` must be NULL, or the distinct names of the factors ",
         "measured within subjects.", call. = FALSE)
  }
  unknown <- setdiff(within, names(factor_lv))
  if (length(unknown) > 0L) {
    stop("`within` names factor(s) that `factor_lv` does not hold: ",
         paste(unknown, collapse = ", "), ". Known factors are: ",
         paste(names(factor_lv), collapse = ", "), ".", call. = FALSE)
  }
  # Put back into `factor_lv` order, so that every table built from either list
  # is in the order the factors were declared in rather than the order they were
  # named in here.
  within <- names(factor_lv)[names(factor_lv) %in% within]
  between <- setdiff(names(factor_lv), within)

  dims <- vapply(factor_lv, length, integer(1))
  cells <- sa_fact_grid(factor_lv)
  n_cells <- nrow(cells)
  cell_label <- sa_fact_cell_labels(factor_lv, cells)

  between_cells <- sa_fact_grid(factor_lv[between])
  n_between <- nrow(between_cells)
  within_cells <- sa_fact_grid(factor_lv[within])
  n_within <- nrow(within_cells)

  if (!is.numeric(n_per_cell) || !length(n_per_cell) %in% c(1L, n_between)) {
    stop("`n_per_cell` must be one size, or one size per combination of the ",
         "between-subject factors, of which this design has ", n_between,
         ". The within-subject factors are crossed with every subject, so ",
         "they cannot hold sizes of their own.", call. = FALSE)
  }
  sizes <- if (length(n_per_cell) == 1L) {
    rep(sa_check_count(n_per_cell, "n_per_cell", 2), n_between)
  } else {
    vapply(seq_along(n_per_cell), function(k) {
      sa_check_count(n_per_cell[k], paste0("n_per_cell[", k, "]"), 2)
    }, integer(1))
  }

  # A unit is a subject when there are within factors and an observation when
  # there are not, which is the same statement either way: a unit sits in one
  # combination of the between factors and contributes one row to each
  # combination of the within ones.
  unit_between <- rep(seq_len(n_between), times = sizes)
  n_units <- length(unit_between)
  n_rows <- n_units * n_within
  subject_idx <- rep(seq_len(n_units), each = n_within)
  within_row <- rep(seq_len(n_within), times = n_units)

  level_idx <- matrix(0L, nrow = n_rows, ncol = length(factor_lv),
                      dimnames = list(NULL, names(factor_lv)))
  for (f in between) {
    level_idx[, f] <- between_cells[[f]][unit_between[subject_idx]]
  }
  for (f in within) {
    level_idx[, f] <- within_cells[[f]][within_row]
  }

  cell_idx <- sa_fact_cell_index(level_idx, dims)
  # Which between-subject combination each cell belongs to, so that a cell can
  # say how many observations it holds without counting rows.
  cell_between <- sa_fact_cell_index(as.matrix(cells[between]), dims[between])

  factors <- lapply(names(factor_lv), function(f) {
    factor_lv[[f]][level_idx[, f]]
  })
  names(factors) <- names(factor_lv)

  list(
    factor_lv   = factor_lv,
    within      = within,
    between     = between,
    cells       = cells,
    n_cells     = n_cells,
    cell_label  = cell_label,
    ref_cell    = 1L,
    cell_n      = sizes[cell_between],
    terms       = sa_fact_terms(names(factor_lv)),
    cell_idx    = cell_idx,
    factors     = factors,
    subject     = if (length(within) > 0L) {
      paste0("subject_", subject_idx)
    } else {
      NULL
    },
    subject_idx = subject_idx,
    n_units     = n_units,
    n_rows      = n_rows
  )
}


#' Pick the factor an effect leans on besides the primary one
#'
#' Drawn at random from the factors after the first, the way the `"single"`
#' profile picks the level it moves, and reported in `truth$partner` so that it
#' is looked up rather than inferred. `"main_only"` needs no partner and takes
#' none, which is also why it draws nothing: the stream is a function of the
#' shapes that were handed out.
#'
#' @keywords internal
#' @noRd
sa_fact_partner <- function(shape, fac_names) {
  if (shape == "main_only") {
    return(NA_character_)
  }
  others <- fac_names[-1L]
  others[sample.int(length(others), 1L)]
}


#' A profile along one factor, centred so that it is a main effect
#'
#' The uncentred profile is the one `simulate_multiple_groups()` plants: zero at
#' the reference level and the magnitude spread over the others by the shape.
#' Centring it turns it into the main effect an ANOVA would report, and since the
#' whole array is shifted to put the reference cell at zero afterwards, the
#' uncentred profile is what `truth_cell$delta` comes back holding anyway.
#'
#' @keywords internal
#' @noRd
sa_fact_profile <- function(d, spread, k) {
  raw <- c(0, sa_sim_pattern_delta(d, spread, k - 1L))
  raw - mean(raw)
}


#' Signs that turn a profile into a pure interaction
#'
#' A main effect along the primary factor becomes an interaction with the partner
#' by being multiplied by a different number at each partner level. For the result
#' to be a *pure* interaction, both of its main effects have to vanish: the
#' primary one vanishes when these numbers average to zero, and the partner one
#' when the profile they multiply is centred, which `sa_fact_profile()` sees to.
#'
#' Alternating signs are the arrangement that makes the effect reverse rather than
#' merely differ, which is what a crossover is. They are centred so that they
#' average to exactly zero at an odd number of levels too, and then scaled so that
#' the largest of them is one, which keeps the size of the interaction the size
#' that was asked for rather than a multiple of it that depends on how many levels
#' the partner has.
#'
#' @keywords internal
#' @noRd
sa_fact_flip <- function(k) {
  s <- rep_len(c(1, -1), k)
  s <- s - mean(s)
  s / max(abs(s))
}


#' Spread one magnitude over the terms according to its shape
#'
#' @param d Signed magnitude of the effect, on the log2 scale.
#' @param shape One of `sa_fact_shapes()`.
#' @param spread `"all"`, `"gradient"` or `"single"`, the profile along a factor.
#' @param mate Partner factor name, or `NA` for `"main_only"`.
#' @param cells Grid of level indices, one row per cell.
#'
#' @return Numeric vector of one effect per cell, in `cells` order, summing to
#'   zero along every factor it does not use.
#'
#' @keywords internal
#' @noRd
sa_fact_plant <- function(d, shape, spread, mate, factor_lv, cells,
                          interaction_scale) {
  primary <- names(factor_lv)[1L]
  profile <- function(f) {
    sa_fact_profile(d, spread, length(factor_lv[[f]]))
  }
  main <- function(f, v) v[cells[[f]]]
  # The primary profile, reversed level by level along the partner. Both main
  # effects of this are exactly zero, so it is interaction and nothing else.
  crossed <- function(v) {
    outer(v, sa_fact_flip(length(factor_lv[[mate]])))[
      cbind(cells[[primary]], cells[[mate]])
    ]
  }

  switch(
    shape,
    main_only     = main(primary, profile(primary)),
    additive      = main(primary, profile(primary)) + main(mate, profile(mate)),
    nuisance_only = main(mate, profile(mate)),
    interaction   = {
      p <- profile(primary)
      main(primary, p) + interaction_scale * crossed(p)
    },
    crossover     = crossed(profile(primary)),
    stop("internal error: unknown effect shape `", shape, "`.", call. = FALSE)
  )
}


#' Feature-level answer
#'
#' @param delta Features by cells, the reference cell being zero throughout.
#'
#' @keywords internal
#' @noRd
sa_fact_truth <- function(feats, delta, design, pattern, spread, direction,
                          partner, baseline, sd_subject) {
  abs_delta <- abs(delta)
  largest <- apply(abs_delta, 1, max)
  # A shape that leaves several cells equally far from the reference leaves no
  # single cell furthest, so the tie is reported rather than broken silently and
  # scored on.
  tied <- rowSums(abs_delta == largest) > 1L
  which_max <- max.col(abs_delta, ties.method = "first")

  extreme_cell <- design$cell_label[which_max]
  extreme_cell[largest == 0] <- NA_character_

  data.frame(
    features     = feats,
    pattern      = pattern,
    spread       = spread,
    direction    = direction,
    partner      = partner,
    extreme_cell = extreme_cell,
    extreme_tied = tied,
    log2fc       = delta[cbind(seq_along(feats), which_max)],
    baseline     = baseline,
    sd_subject   = sd_subject,
    stringsAsFactors = FALSE
  )
}


#' Per feature and term answer, in the row order an ANOVA table follows
#'
#' @param planted One logical per feature. An unplanted feature is zero in every
#'   cell, so its every component is zero too and decomposing it would be work
#'   done to arrive back at the matrix it started as.
#'
#' @keywords internal
#' @noRd
sa_fact_truth_term <- function(feats, delta, design, planted) {
  terms <- design$terms
  n_terms <- length(terms)
  labels <- sa_fact_term_labels(terms)
  orders <- lengths(terms)
  is_within <- vapply(terms, function(t) any(t %in% design$within), logical(1))

  comp <- matrix(0, nrow = length(feats), ncol = n_terms)
  for (i in which(planted)) {
    for (k in seq_len(n_terms)) {
      comp[i, k] <- max(abs(sa_fact_component(delta[i, ], design$cells,
                                              terms[[k]])))
    }
  }
  flat <- as.vector(t(comp))

  data.frame(
    features      = rep(feats, each = n_terms),
    terms         = rep(labels, times = length(feats)),
    term_order    = rep(orders, times = length(feats)),
    is_within     = rep(is_within, times = length(feats)),
    max_abs_delta = flat,
    is_effect     = flat > 0,
    stringsAsFactors = FALSE
  )
}


#' Per feature and cell answer
#'
#' @keywords internal
#' @noRd
sa_fact_truth_cell <- function(feats, delta, center, sd_mat, design) {
  n_cells <- design$n_cells
  n_feats <- length(feats)

  out <- data.frame(features = rep(feats, each = n_cells),
                    stringsAsFactors = FALSE)
  for (f in names(design$factor_lv)) {
    out[[f]] <- rep(design$factor_lv[[f]][design$cells[[f]]], times = n_feats)
  }
  out$is_ref <- rep(seq_len(n_cells) == design$ref_cell, times = n_feats)
  out$delta <- as.vector(t(delta))
  out$center <- as.vector(t(center))
  out$sd <- as.vector(t(sd_mat))
  out$n <- rep(design$cell_n, times = n_feats)
  out
}


#' Per feature and pair answer, in a post-hoc table's own direction
#'
#' The rows and the cells each of them averages come from
#' `sa_fact_contrast_skeleton()`, which `compare_factorial_groups()` reads too,
#' so the row order and the `group1 - group2` direction cannot drift apart from
#' the tables this is meant to score.
#'
#' @keywords internal
#' @noRd
sa_fact_truth_contrast <- function(feats, delta, design) {
  skel <- sa_fact_contrast_skeleton(design)
  n_rows <- nrow(skel$table)

  mat <- matrix(0, nrow = length(feats), ncol = n_rows)
  for (k in seq_len(n_rows)) {
    mat[, k] <- rowMeans(delta[, skel$sel1[[k]], drop = FALSE]) -
      rowMeans(delta[, skel$sel2[[k]], drop = FALSE])
  }
  # Averaging a set of cells divides, so a contrast that is zero by construction
  # can come back a rounding away from it. `is_diff` has to be exact.
  mat[abs(mat) < sa_fact_tol()] <- 0
  flat <- as.vector(t(mat))

  out <- cbind(
    data.frame(features = rep(feats, each = n_rows), stringsAsFactors = FALSE),
    skel$table[rep(seq_len(n_rows), times = length(feats)), , drop = FALSE],
    data.frame(delta = flat, is_diff = flat != 0, stringsAsFactors = FALSE)
  )
  rownames(out) <- NULL
  out
}
