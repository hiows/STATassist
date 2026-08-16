#' Simulate a contingency table whose association is known
#'
#' Generates two or more categorical variables with a chosen amount of
#' association planted between them, and returns the planted answer alongside the
#' data. Every quantity [compare_categorical_groups()] estimates can then be
#' checked against what was actually put there, which is what a real data set can
#' never offer.
#'
#' The point of the exercise is the gap. A table drawn from a distribution is not
#' that distribution: the observed Cramer's V carries a sampling error of its own,
#' and it is biased upwards, because every departure from independence counts
#' towards it whether it was planted or drawn. So a table simulated at
#' `assoc = 0` does not come back with an estimated association of zero, and
#' reading that number as the answer rather than the p-value beside it is the
#' mistake the simulator exists to make visible.
#'
#' @section How the association is planted:
#' The margins are fixed first, and then mass is moved between the cells in a way
#' that leaves them exactly where they were. That is what `assoc` scales: a
#' perturbation matrix whose every row and column sums to zero, added to the
#' independent joint distribution.
#'
#' Keeping the margins fixed is what makes the null hypothesis the only thing
#' that moves. A test of independence compares the table to the product of its own
#' observed margins, so a construction that shifted the margins too would plant an
#' association and change what it is being measured against at the same time.
#'
#' `assoc` is the fraction of the largest step of that shape the table can take
#' before a cell reaches zero, so it runs from 0 to 1 whatever the margins are.
#' `assoc = 0` is the product of the margins **exactly**, so it is null in the
#' strict sense and is what a type I error rate can be measured on.
#' `assoc = 1` puts a structural zero in the table, which is the one setting
#' where the exact test and the approximation part company sharply.
#'
#' @section The shape of the association:
#' `pattern` picks which cells the mass moves between. All three keep the margins
#' fixed, so all three are the same kind of departure at different places.
#'
#' \describe{
#'   \item{`"corner"`}{Mass moves into the two cells where both variables sit at
#'     their first level or both at their second, and out of the two where they
#'     disagree. On a 2 x 2 table this is the whole of the association, and the
#'     odds ratio it plants is above 1.}
#'   \item{`"single"`}{One level of each variable is special and the rest are
#'     alike: mass moves into the cell where both first levels meet and is taken
#'     evenly from the others. This is the shape a test of homogeneity is usually
#'     looking for, one category behaving differently from the pack.}
#'   \item{`"gradient"`}{A monotone association, planted by a centred linear ramp
#'     over the levels of each variable in the order `category_lv` gives them.
#'     This is the shape an ordered variable such as `high` / `mid` / `low`
#'     actually takes, and the one a chi-square test is least efficient at
#'     finding, since it spends degrees of freedom on departures that are not
#'     there.}
#' }
#'
#' @section A matched design:
#' `paired = TRUE` generates repeated binary measurements on the same subject
#' rather than a cross-classification, so `assoc`, `pattern` and `margins` do not
#' apply and `discordance` takes over. Every subject starts at either level with
#' equal probability, and then moves between consecutive conditions: a subject at
#' the first level takes the second with probability `discordance[1]`, and one at
#' the second takes the first with probability `discordance[2]`.
#'
#' Two things follow from that, and both are what the matched tests are about.
#' Over two conditions the planted paired odds ratio is exactly
#' `discordance[1] / discordance[2]`, because the equal start makes the two
#' discordant cells half of each transition probability. Over three or more, the
#' response rate climbs from condition to condition whenever the first
#' probability exceeds the second, and that climb is what Cochran's Q is testing
#' for. Equal probabilities are the strict null of both: the transitions cancel,
#' the rate stays at one half, and nothing is planted.
#'
#' @section Which arguments each design reads:
#' A cross-classification reads `margins`, `assoc` and `pattern`. A matched design
#' reads `discordance`. Passing one design an argument belonging to the other is a
#' warning naming both lists rather than a value silently ignored, because the two
#' designs plant different things and a call that names the wrong knob is a call
#' that expected the other design.
#'
#' @param n_samples Number of observations, meaning rows for a cross-classified
#'   design and subjects for a matched one.
#' @param category_lv Named list giving the levels of each variable, passed
#'   straight through to [compare_categorical_groups()], so its names are the
#'   column names of the generated data and its order fixes the reference level of
#'   each variable. A matched design needs one binary level set shared by every
#'   condition, and the default changes to reflect that when `paired = TRUE`.
#' @param margins Named list of the marginal probabilities of each variable, in
#'   the order its levels are given, or `NULL` for a uniform margin. Each entry is
#'   normalised to sum to one, so relative weights are enough.
#' @param assoc How much association to plant, from 0 for exact independence to 1
#'   for the largest margin-preserving departure the margins allow.
#' @param pattern Which cells the association moves mass between. See "The shape
#'   of the association".
#' @param paired Logical. If `TRUE`, the columns are repeated binary measurements
#'   on the same subject rather than different variables.
#' @param discordance The two transition probabilities of a matched design, from
#'   the first level to the second and back. Their ratio is the planted paired
#'   odds ratio, and their being equal is the strict null.
#' @param seed Seed for the draw, restored on exit, so the caller's random stream
#'   is left as it was found.
#'
#' @return A plain list, the shape every simulator in the package returns:
#'
#'   \describe{
#'     \item{`args`}{`data`, `category_lv` and `paired`, named after the arguments
#'       of [compare_categorical_groups()] so that
#'       `do.call(compare_categorical_groups, sim$args)` runs the analysis this
#'       data was made for.}
#'     \item{`truth`}{One row summarising what was planted: `n_samples`,
#'       `pattern`, and the population value of every association measure the
#'       comparison reports and this design defines. A cross-classification adds
#'       `assoc` and `cramers_v`, plus `phi_coefficient` and `odds_ratio` on a
#'       2 x 2 table. A matched design adds the transition probabilities, and then
#'       either the three paired measures or the sequence of response rates.
#'       These are the values the estimates should approach as `n_samples` grows,
#'       and they carry no sampling error.}
#'     \item{`truth_cell`}{One row per cell: `row_level`, `col_level`,
#'       `p_independent`, `p_planted`, `lift` and `expected_n`. It merges with
#'       `$cells` of the comparison on `c("row_level", "col_level")` with neither
#'       side renamed, so a residual can be read against the departure that
#'       produced it. A matched pair of conditions adds `p_symmetric` and
#'       `expected_symmetry_n`, which are what `$cells$expected` holds there:
#'       symmetry is the null that design is tested against, so it is the one the
#'       planted table has to be scored on. For three or more conditions the cells
#'       are condition by response, and there marginal homogeneity and
#'       independence are the same arithmetic, so `p_independent` already is the
#'       null.}
#'   }
#'
#' @seealso [compare_categorical_groups()] for the analysis this feeds, and
#'   [simulate_two_groups()] for the numeric counterpart.
#'
#' @examples
#' ## A planted association, recovered by both tests.
#' sim <- simulate_categorical_groups(n_samples = 300, assoc = 0.4, seed = 1)
#' fit <- do.call(compare_categorical_groups, sim$args)
#' fit$tests$chisq_test$pval
#'
#' ## The planted strength beside the estimated one. The estimate is high, which
#' ## is what an association measure built from a sum of squares always is.
#' cbind(planted = sim$truth$cramers_v,
#'       estimated = fit$association$estimate[1])
#'
#' ## The residual of each cell against the lift that was planted there.
#' scored <- merge(fit$cells, sim$truth_cell,
#'                 by = c("row_level", "col_level"))
#' scored[c("row_level", "col_level", "lift", "std_residual")]
#'
#' ## `assoc = 0` is the product of the margins exactly, so it is null in the
#' ## strict sense: the lift of every cell is 1.
#' null <- simulate_categorical_groups(n_samples = 300, assoc = 0, seed = 2)
#' range(null$truth_cell$lift)
#'
#' ## A monotone association over ordered levels, which is the shape a chi-square
#' ## test spends degrees of freedom it does not need on. A 3 x 3 table of this
#' ## many observations is past what Fisher's enumeration can walk, so the
#' ## simulated variant is what answers there.
#' ordered <- simulate_categorical_groups(
#'   n_samples   = 300,
#'   category_lv = list(dose = c("low", "mid", "high"),
#'                      response = c("none", "partial", "full")),
#'   assoc       = 0.5,
#'   pattern     = "gradient",
#'   seed        = 3
#' )
#' fit_ordered <- do.call(compare_categorical_groups,
#'                        c(ordered$args, list(simulate_p_value = TRUE)))
#' as.table(fit_ordered)
#'
#' ## A matched design plants the paired odds ratio as a ratio of transitions,
#' ## and the cells carry what symmetry expects rather than what independence
#' ## would, since symmetry is the null McNemar's test is about.
#' matched <- simulate_categorical_groups(n_samples = 200, paired = TRUE,
#'                                        discordance = c(0.3, 0.1), seed = 4)
#' matched$truth$odds_ratio_paired
#' matched$truth_cell[c("row_level", "col_level", "p_symmetric",
#'                      "expected_symmetry_n")]
#'
#' @export
simulate_categorical_groups <- function(n_samples = 200,
                                        category_lv = NULL,
                                        margins = NULL,
                                        assoc = 0.3,
                                        pattern = c("corner", "single",
                                                    "gradient"),
                                        paired = FALSE,
                                        discordance = c(0.25, 0.10),
                                        seed = NULL) {

  # Recorded before anything is reassigned, since `match.arg()` on `pattern`
  # would make `missing()` say it was supplied whether it was or not.
  given <- c(margins = !is.null(margins), assoc = !missing(assoc),
             pattern = !missing(pattern), discordance = !missing(discordance))

  pattern <- match.arg(pattern)
  sa_check_flag(paired, "paired")
  n_samples <- sa_check_count(n_samples, "n_samples", 2)
  sa_check_scalar_num(assoc, "assoc", 0, 1)

  if (is.null(category_lv)) {
    category_lv <- if (paired) {
      list(before = c("fail", "pass"), after = c("fail", "pass"))
    } else {
      list(cat_1 = c("y", "n"), cat_2 = c("high", "mid", "low"))
    }
  }
  category_lv <- sa_sim_cat_levels(category_lv)

  sa_sim_cat_check_given(given, paired)

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  if (paired) {
    sa_sim_cat_matched(n_samples, category_lv, discordance)
  } else {
    sa_sim_cat_crossed(n_samples, category_lv, margins, assoc, pattern)
  }
}


#' Warn about an argument belonging to the other design
#'
#' Both lists are named rather than only the offending one. An argument passed to
#' the wrong design is usually a call that meant the other design, so knowing what
#' this one does read is the part that resolves it.
#'
#' @param given Named logical, which optional arguments the caller supplied.
#' @param paired Which design is running.
#'
#' @keywords internal
#' @noRd
sa_sim_cat_check_given <- function(given, paired) {
  reads <- if (paired) "discordance" else c("margins", "assoc", "pattern")
  ignored <- setdiff(names(given)[given], reads)
  if (length(ignored) == 0L) {
    return(invisible(NULL))
  }

  warning("a ", if (paired) "matched" else "cross-classified",
          " design reads ", paste(reads, collapse = ", "), ", so the ",
          "value(s) given for ", paste(ignored, collapse = ", "),
          " were ignored. Set `paired = ", if (paired) "FALSE" else "TRUE",
          "` for the design those belong to.", call. = FALSE)
  invisible(NULL)
}


#' Check the level sets a simulated table is built on
#'
#' @keywords internal
#' @noRd
sa_sim_cat_levels <- function(category_lv) {
  if (!is.list(category_lv) || length(category_lv) < 2L ||
        is.null(names(category_lv)) || anyNA(names(category_lv)) ||
        !all(nzchar(names(category_lv))) ||
        anyDuplicated(names(category_lv)) > 0L) {
    stop("`category_lv` must be a named list of at least two variables, each ",
         "entry holding that variable's levels.", call. = FALSE)
  }
  lapply(stats::setNames(names(category_lv), names(category_lv)), function(nm) {
    lv <- as.character(category_lv[[nm]])
    if (length(lv) < 2L || anyNA(lv) || anyDuplicated(lv) > 0L) {
      stop("`category_lv$", nm, "` must hold at least two distinct ",
           "non-missing levels.", call. = FALSE)
    }
    lv
  })
}


#' Draw a cross-classification with a planted, margin-preserving association
#'
#' @keywords internal
#' @noRd
sa_sim_cat_crossed <- function(n_samples, category_lv, margins, assoc, pattern) {
  if (length(category_lv) != 2L) {
    stop("a cross-classified simulation plants an association between exactly ",
         "two variables, and `category_lv` names ", length(category_lv),
         ". Set `paired = TRUE` for repeated measurements of one thing.",
         call. = FALSE)
  }
  variables <- names(category_lv)
  row_lv <- category_lv[[1]]
  col_lv <- category_lv[[2]]

  probs <- sa_sim_cat_margins(margins, category_lv)
  independent <- outer(probs[[1]], probs[[2]])
  dimnames(independent) <- list(row_lv, col_lv)

  planted <- sa_sim_cat_perturb(independent, assoc, pattern)

  # `sample()` over the flattened joint keeps the two draws in step, which
  # drawing each variable in turn could not: an association is a statement about
  # the pair, so the pair is what is drawn.
  drawn <- sample(seq_along(planted), size = n_samples, replace = TRUE,
                  prob = as.numeric(planted))
  row_idx <- (drawn - 1L) %% length(row_lv) + 1L
  col_idx <- (drawn - 1L) %/% length(row_lv) + 1L

  data <- data.frame(row_lv[row_idx], col_lv[col_idx],
                     stringsAsFactors = FALSE)
  names(data) <- variables
  sa_sim_cat_check_drawn(data, category_lv)

  list(
    args = list(data = data, category_lv = category_lv, paired = FALSE),
    truth = sa_sim_cat_truth(planted, independent, n_samples, assoc, pattern),
    truth_cell = sa_sim_cat_truth_cell(planted, independent, n_samples)
  )
}


#' Whether every named level was actually drawn
#'
#' A level the draw happened to miss leaves an empty row or column, which
#' [compare_categorical_groups()] refuses rather than test, so the simulator says
#' which one it was here instead of the analysis failing on data it was handed.
#'
#' @keywords internal
#' @noRd
sa_sim_cat_check_drawn <- function(data, category_lv) {
  absent <- unlist(lapply(names(category_lv), function(nm) {
    missed <- setdiff(category_lv[[nm]], unique(data[[nm]]))
    if (length(missed) == 0L) NULL else paste0(nm, ": ",
                                               paste(missed, collapse = ", "))
  }))
  if (length(absent) > 0L) {
    warning("no row was drawn at level(s) ", paste(absent, collapse = "; "),
            ", which leaves an empty row or column that no test of ",
            "independence can be run on. Raise `n_samples`, or give the level ",
            "more weight in `margins`.", call. = FALSE)
  }
  invisible(data)
}


#' Normalise the marginal probabilities of each variable
#'
#' @keywords internal
#' @noRd
sa_sim_cat_margins <- function(margins, category_lv) {
  if (is.null(margins)) {
    return(lapply(category_lv, function(lv) rep(1 / length(lv), length(lv))))
  }
  if (!is.list(margins) || is.null(names(margins)) ||
        !setequal(names(margins), names(category_lv))) {
    stop("`margins` must be a named list holding one entry per variable of ",
         "`category_lv`: ", paste(names(category_lv), collapse = ", "), ".",
         call. = FALSE)
  }
  lapply(stats::setNames(names(category_lv), names(category_lv)), function(nm) {
    m <- margins[[nm]]
    sa_check_num_vector(m, paste0("margins$", nm), 0)
    if (length(m) != length(category_lv[[nm]])) {
      stop("`margins$", nm, "` must hold one weight per level of `category_lv$",
           nm, "`: got ", length(m), " for ", length(category_lv[[nm]]),
           " level(s).", call. = FALSE)
    }
    if (sum(m) <= 0 || any(m == 0)) {
      stop("`margins$", nm, "` must give every level a positive weight; a level ",
           "with none is a level to leave out of `category_lv`.", call. = FALSE)
    }
    m / sum(m)
  })
}


#' Move mass between the cells without moving the margins
#'
#' The perturbation is the outer product of two vectors that each sum to zero,
#' which is what makes every row and column sum of it zero as well: the row sums
#' are one vector scaled by the total of the other, and that total is zero. One
#' rank is enough for every pattern here, and it is what keeps the step a single
#' number to scale.
#'
#' @keywords internal
#' @noRd
sa_sim_cat_perturb <- function(independent, assoc, pattern) {
  if (assoc == 0) {
    return(independent)
  }
  r <- nrow(independent)
  c_ <- ncol(independent)

  zero_sum <- function(k) {
    switch(pattern,
           corner   = c(1, -1, rep(0, k - 2L)),
           single   = c(1, rep(-1 / (k - 1L), k - 1L)),
           gradient = {
             ramp <- seq(1, -1, length.out = k)
             ramp - mean(ramp)
           })
  }
  step <- outer(zero_sum(r), zero_sum(c_))

  # The largest step this shape can take is the one that first empties a cell.
  # Scaling by it is what gives `assoc` the same meaning whatever the margins
  # are, rather than a meaning that has to be found by trial for each table.
  losing <- step < 0
  max_step <- min(independent[losing] / -step[losing])

  out <- independent + assoc * max_step * step
  # A cell that reaches exactly zero at `assoc = 1` can land a hair below it
  # through floating point, and a negative probability is not something to hand
  # to `sample()`.
  out[out < 0] <- 0
  dimnames(out) <- dimnames(independent)
  out
}


#' The population value of every measure this design defines
#'
#' Computed from the planted distribution rather than from the drawn table, so it
#' carries no sampling error and is the number an estimate should approach.
#'
#' The mean square contingency `sum((p - p_ind)^2 / p_ind)` is what a chi-square
#' statistic divided by `n` estimates, so Cramer's V follows from it without `n`
#' entering anywhere. The contingency coefficient does not, since its denominator
#' holds `n` itself, so it is not reported here.
#'
#' @keywords internal
#' @noRd
sa_sim_cat_truth <- function(planted, independent, n_samples, assoc, pattern) {
  phi_sq <- sum((planted - independent)^2 / independent)
  min_df <- min(dim(planted)) - 1L

  out <- data.frame(
    n_samples = n_samples,
    pattern   = pattern,
    assoc     = assoc,
    cramers_v = sqrt(phi_sq / min_df),
    stringsAsFactors = FALSE
  )

  if (identical(dim(planted), c(2L, 2L))) {
    cross <- planted[1, 1] * planted[2, 2] - planted[1, 2] * planted[2, 1]
    out$phi_coefficient <- cross /
      sqrt(prod(c(rowSums(planted), colSums(planted))))
    out$odds_ratio <- sa_finite_or_na(
      planted[1, 1] * planted[2, 2] / (planted[1, 2] * planted[2, 1])
    )
  }
  rownames(out) <- NULL
  out
}


#' The planted distribution, cell by cell
#'
#' Keyed on `row_level` and `col_level`, the columns `$cells` of a comparison
#' carries, so the two merge with neither side renamed.
#'
#' @param planted The joint distribution the data was drawn from.
#' @param independent The product of its margins, which is what a test of
#'   independence holds it against.
#' @param total How many observations the table counts, which is the number of
#'   rows for a cross-classification and the number of subjects times the number
#'   of conditions for a repeated one. `p_planted` is a share of the whole table
#'   either way, so it is the only thing that turns it back into a count.
#'
#' @keywords internal
#' @noRd
sa_sim_cat_truth_cell <- function(planted, independent, total) {
  grid <- expand.grid(row_level = rownames(planted),
                      col_level = colnames(planted),
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  out <- data.frame(
    row_level     = grid$row_level,
    col_level     = grid$col_level,
    p_independent = as.numeric(independent),
    p_planted     = as.numeric(planted),
    lift          = sa_finite_or_na(as.numeric(planted / independent)),
    expected_n    = total * as.numeric(planted),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}


#' Draw repeated binary measurements with a planted transition
#'
#' Every subject starts at either level with equal probability, which is what
#' makes the two discordant cells half of their transition probabilities and the
#' planted paired odds ratio their ratio. The chain then runs across the
#' conditions, so the response rate of condition `j + 1` follows from that of
#' condition `j` and nothing else, and the whole sequence of rates is known before
#' a single subject is drawn.
#'
#' @keywords internal
#' @noRd
sa_sim_cat_matched <- function(n_samples, category_lv, discordance) {
  if (!is.numeric(discordance) || length(discordance) != 2L) {
    stop("`discordance` must be two transition probabilities, from the first ",
         "level to the second and back.", call. = FALSE)
  }
  sa_check_num_vector(discordance, "discordance", 0, 1)

  levels <- category_lv[[1]]
  if (!all(vapply(category_lv, function(lv) identical(lv, levels),
                  logical(1)))) {
    stop("a matched simulation measures one thing repeatedly, so every entry ",
         "of `category_lv` holds the same levels in the same order. They ",
         "differ.", call. = FALSE)
  }
  if (length(levels) != 2L) {
    stop("a matched simulation needs binary conditions, and `category_lv` ",
         "holds ", length(levels), " level(s). McNemar's test and Cochran's Q ",
         "are both about a binary response.", call. = FALSE)
  }

  variables <- names(category_lv)
  k <- length(variables)
  move_up <- discordance[1]
  move_down <- discordance[2]

  # The response, here and in the comparison, is the second level: the first is
  # the reference every other scenario reads against.
  state <- matrix(NA, nrow = n_samples, ncol = k,
                  dimnames = list(NULL, variables))
  state[, 1] <- stats::runif(n_samples) < 0.5
  for (j in seq_len(k - 1L)) {
    move <- stats::runif(n_samples) < ifelse(state[, j], move_down, move_up)
    state[, j + 1L] <- xor(state[, j], move)
  }

  rates <- numeric(k)
  rates[1] <- 0.5
  for (j in seq_len(k - 1L)) {
    rates[j + 1L] <- rates[j] * (1 - move_down) + (1 - rates[j]) * move_up
  }

  data <- as.data.frame(
    lapply(seq_len(k), function(j) levels[state[, j] + 1L]),
    stringsAsFactors = FALSE
  )
  names(data) <- variables

  list(
    args = list(data = data, category_lv = category_lv, paired = TRUE),
    truth = sa_sim_cat_truth_matched(n_samples, k, discordance, rates),
    truth_cell = sa_sim_cat_truth_cell_matched(n_samples, k, variables, levels,
                                               discordance, rates)
  )
}


#' What a matched simulation planted
#'
#' @keywords internal
#' @noRd
sa_sim_cat_truth_matched <- function(n_samples, k, discordance, rates) {
  move_up <- discordance[1]
  move_down <- discordance[2]
  b <- 0.5 * move_up
  c_ <- 0.5 * move_down

  out <- data.frame(
    n_samples    = n_samples,
    pattern      = "transition",
    n_conditions = k,
    move_up      = move_up,
    move_down    = move_down,
    stringsAsFactors = FALSE
  )

  if (k == 2L) {
    out$odds_ratio_paired <- sa_finite_or_na(b / c_)
    out$risk_difference_paired <- b - c_
    out$cohens_g <- sa_finite_or_na(b / (b + c_) - 0.5)
  } else {
    # Three or more conditions have no single pair to be about, so what was
    # planted is the climb in the response rate, which is what Q is testing.
    out$rate_first <- rates[1]
    out$rate_last <- rates[k]
    out$rate_range <- max(rates) - min(rates)
  }
  rownames(out) <- NULL
  out
}


#' The planted distribution of a matched design, cell by cell
#'
#' Over two conditions this is the paired square table McNemar's test reads, and
#' the null it is read against is symmetry, so `p_symmetric` is added: it is the
#' average of each cell and its transpose, which is the share `$cells$expected`
#' holds there. The diagonal is symmetric by construction, so its `p_symmetric`
#' equals its `p_planted` and nothing was planted there to find.
#'
#' Over three or more conditions there is no such table, so it is the condition by
#' response table Cochran's Q is about, which is the one the comparison carries as
#' well. Marginal homogeneity and independence are the same arithmetic on that
#' table, so `p_independent` already is the null and no column is added.
#'
#' @keywords internal
#' @noRd
sa_sim_cat_truth_cell_matched <- function(n_samples, k, variables, levels,
                                          discordance, rates) {
  if (k == 2L) {
    move_up <- discordance[1]
    move_down <- discordance[2]
    planted <- matrix(
      c(0.5 * (1 - move_up), 0.5 * move_down,
        0.5 * move_up,       0.5 * (1 - move_down)),
      nrow = 2L, dimnames = list(levels, levels)
    )
    independent <- outer(rowSums(planted), colSums(planted))
    dimnames(independent) <- dimnames(planted)

    out <- sa_sim_cat_truth_cell(planted, independent, n_samples)
    symmetric <- (planted + t(planted)) / 2
    out$p_symmetric <- as.numeric(symmetric)
    out$expected_symmetry_n <- n_samples * as.numeric(symmetric)
    return(out)
  }

  planted <- cbind(1 - rates, rates) / k
  dimnames(planted) <- list(variables, levels)
  independent <- outer(rowSums(planted), colSums(planted))
  dimnames(independent) <- dimnames(planted)
  # Every subject is measured under every condition, so the table this counts is
  # subjects times conditions rather than subjects.
  sa_sim_cat_truth_cell(planted, independent, n_samples * k)
}
