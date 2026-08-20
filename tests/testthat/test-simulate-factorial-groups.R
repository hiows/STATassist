# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# The value of a simulator is that the answer is known, so what is pinned here is
# mostly that the answer it reports is the answer it planted. Beyond what the
# one-factor simulator has to promise, this one plants effects across the terms of
# a crossed model, so the term each shape claims to move, and just as importantly
# the terms it claims to leave alone, are pinned exactly.

small <- function(...) {
  given <- list(...)
  base <- list(n_feats = 12, n_up = 3, n_down = 3, n_per_cell = 4)
  do.call(simulate_factorial_groups,
          c(given, base[setdiff(names(base), names(given))]))
}

three <- function(...) {
  small(factor_lv = list(treatment = c("control", "treat_A", "treat_B"),
                         sex       = c("male", "female"),
                         time      = c("T0", "T1", "T2")),
        ...)
}

cells_of <- function(sim, feature) {
  sim$truth_cell[sim$truth_cell$features == feature, ]
}

# Which terms an answer table says moved, for one feature.
terms_of <- function(sim, feature) {
  rows <- sim$truth_term[sim$truth_term$features == feature, ]
  rows$terms[rows$is_effect]
}


test_that("the returned arguments name the factors and their levels", {
  sim <- small(seed = 1)

  expect_named(sim, c("args", "truth", "truth_term", "truth_cell",
                      "truth_contrast"))
  expect_named(sim$args, c("data", "feats", "factors", "factor_lv",
                           "input_scale"))
  expect_identical(sim$args$input_scale, "log2")
  expect_identical(names(sim$args$data), sim$args$feats)

  # `factors` says which level every row sits at and `factor_lv` says the order
  # those levels are in, which is the pair `group` and `group_lv` make in
  # simulate_multiple_groups(). One is as long as the data, the other is not.
  expect_named(sim$args$factors, c("treatment", "sex"))
  expect_true(all(vapply(sim$args$factors, length, integer(1)) ==
                    nrow(sim$args$data)))
  expect_identical(sim$args$factor_lv$treatment,
                   c("control", "treat_A", "treat_B", "treat_C"))
  expect_identical(sim$args$factor_lv$sex, c("male", "female"))
  # The reference cell is the first level of every factor.
  expect_identical(sim$truth_cell$treatment[1], "control")
  expect_identical(sim$truth_cell$sex[1], "male")
  expect_true(sim$truth_cell$is_ref[1])
  expect_equal(sum(sim$truth_cell$is_ref), 12L)
})

test_that("a within factor carries its pairing in the arguments", {
  sim <- three(seed = 1, within = "time")

  expect_named(sim$args, c("data", "feats", "factors", "factor_lv", "within",
                           "id", "input_scale"))
  expect_identical(sim$args$within, "time")
  # Treatment and sex are between subjects, so there are 3 x 2 combinations of
  # them holding four subjects each, and every subject is seen at every time.
  expect_equal(nrow(sim$args$data), 3L * 2L * 4L * 3L)
  expect_length(unique(sim$args$id), 24L)
  expect_equal(unname(table(sim$args$id)), rep(3L, 24L), ignore_attr = TRUE)
  # A subject belongs to one combination of the between factors and stays there.
  by_subject <- split(sim$args$factors$treatment, sim$args$id)
  expect_true(all(vapply(by_subject, function(x) length(unique(x)), integer(1))
                  == 1L))
  # And is seen once at each level of the within factor.
  by_time <- split(sim$args$factors$time, sim$args$id)
  expect_true(all(vapply(by_time, function(x) length(unique(x)), integer(1))
                  == 3L))
})

test_that("the answer tables are aligned with the features, cells and terms", {
  sim <- small(seed = 1)

  expect_identical(sim$truth$features, sim$args$feats)
  expect_named(sim$truth, c("features", "pattern", "spread", "direction",
                            "partner", "extreme_cell", "extreme_tied",
                            "log2fc", "baseline", "sd_subject"))
  expect_named(sim$truth_term, c("features", "terms", "term_order",
                                 "is_within", "max_abs_delta", "is_effect"))
  expect_named(sim$truth_cell, c("features", "treatment", "sex", "is_ref",
                                 "delta", "center", "sd", "n"))
  expect_named(sim$truth_contrast, c("features", "factor", "stratum",
                                     "contrast", "group1", "group2", "delta",
                                     "is_diff"))

  # Eight cells and three terms in a 4 x 2 design.
  expect_equal(nrow(sim$truth_cell), 12L * 8L)
  expect_equal(nrow(sim$truth_term), 12L * 3L)
  expect_identical(unique(sim$truth_term$terms),
                   c("treatment", "sex", "treatment:sex"))
  # Main effects come before the interaction, which is the order an ANOVA table
  # lists them in.
  expect_equal(sim$truth_term$term_order[1:3], c(1L, 1L, 2L))
})

test_that("every cell gets the size it was given", {
  balanced <- small(seed = 1, n_per_cell = 5)
  expect_equal(nrow(balanced$args$data), 8L * 5L)
  expect_true(all(balanced$truth_cell$n == 5L))
  expect_equal(unname(table(balanced$args$factors$treatment)),
               rep(10L, 4L), ignore_attr = TRUE)

  # One size per combination of the between factors, so an unbalanced design is
  # written rather than assembled out of one recycled number. The grid runs with
  # the first factor varying fastest, so these are the four treatments at male
  # followed by the four at female.
  uneven <- small(seed = 1, n_per_cell = c(5, 6, 7, 8, 9, 10, 11, 12))
  expect_equal(nrow(uneven$args$data), sum(5:12))
  expect_equal(uneven$truth_cell$n[1:8], 5:12)
  expect_equal(unname(table(uneven$args$factors$sex)[c("male", "female")]),
               c(sum(5:8), sum(9:12)), ignore_attr = TRUE)
})

test_that("only the between factors carry sizes of their own", {
  # With time within, a size is one per treatment-by-sex combination and each of
  # them is a count of subjects rather than of rows.
  sim <- three(seed = 1, within = "time", n_per_cell = c(2, 3, 4, 5, 6, 7))
  expect_length(unique(sim$args$id), sum(2:7))
  expect_equal(nrow(sim$args$data), sum(2:7) * 3L)
  # Every cell of the full grid holds as many observations as its between-subject
  # combination holds subjects, because each subject contributes one to each.
  first <- cells_of(sim, sim$args$feats[1])
  expect_equal(first$n[first$time == "T0"], as.integer(2:7))
  expect_equal(first$n[first$time == "T2"], as.integer(2:7))

  expect_error(three(within = "time", n_per_cell = rep(4, 18)),
               "one size per combination")
})

test_that("exactly the requested number of features is planted", {
  sim <- small(n_feats = 40, n_up = 7, n_down = 11, seed = 3)

  expect_equal(sum(sim$truth$direction == "up"), 7L)
  expect_equal(sum(sim$truth$direction == "down"), 11L)
  expect_equal(sum(sim$truth$direction == "none"), 22L)
  expect_equal(sum(sim$truth$pattern == "none"), 22L)
})

test_that("both mixes are handed out in the proportions asked for", {
  # Split by largest remainder rather than by lot, so the counts are a function
  # of the arguments and can be pinned exactly. Ten up features over five equal
  # weights is two each.
  even <- small(n_feats = 40, n_up = 10, n_down = 10, seed = 1)
  expect_equal(unname(table(even$truth$pattern)[sa_fact_shapes()]),
               rep(4L, 5L), ignore_attr = TRUE)

  # The remainder goes to the earlier shapes, in the order the weights are in.
  odd <- small(n_feats = 40, n_up = 7, n_down = 0, seed = 1)
  expect_equal(unname(table(odd$truth$pattern)[sa_fact_shapes()]),
               c(2L, 2L, 1L, 1L, 1L), ignore_attr = TRUE)

  weighted <- small(n_feats = 40, n_up = 8, n_down = 0, seed = 1,
                    term_mix = c(crossover = 3, additive = 1))
  expect_equal(sum(weighted$truth$pattern == "crossover"), 6L)
  expect_equal(sum(weighted$truth$pattern == "additive"), 2L)
  expect_equal(sum(weighted$truth$pattern == "main_only"), 0L)

  # `pattern_mix` is on the other axis and is split the same way, so the two are
  # crossed without either set of counts moving with the seed.
  expect_equal(unname(table(even$truth$spread)[c("all", "gradient",
                                                 "single")]),
               c(8L, 6L, 6L), ignore_attr = TRUE)
  narrow <- small(n_feats = 40, n_up = 10, n_down = 10, seed = 1,
                  pattern_mix = c(all = 1, gradient = 0, single = 0))
  expect_true(all(narrow$truth$spread[narrow$truth$pattern != "none"] ==
                    "all"))
})

test_that("each shape moves exactly the terms it says it does", {
  sim <- small(n_feats = 40, n_up = 10, n_down = 10, seed = 2)
  moved <- function(shape) {
    lapply(sim$truth$features[sim$truth$pattern == shape], terms_of,
           sim = sim)
  }
  all_are <- function(shape, expected) {
    got <- moved(shape)
    expect_gt(length(got), 0L)
    for (g in got) expect_setequal(g, expected)
  }

  all_are("main_only", "treatment")
  all_are("additive", c("treatment", "sex"))
  all_are("interaction", c("treatment", "treatment:sex"))
  all_are("nuisance_only", "sex")
  # The shape the whole design exists for: the cells differ and neither main
  # effect does. A main-effect test has to miss this and an interaction test has
  # to catch it, so "roughly zero" would not do.
  all_are("crossover", "treatment:sex")
  cross <- sim$truth$features[sim$truth$pattern == "crossover"]
  mains <- subset(sim$truth_term, features %in% cross & term_order == 1)
  expect_true(all(mains$max_abs_delta == 0))
  expect_gt(min(subset(sim$truth_term,
                       features %in% cross & term_order == 2)$max_abs_delta), 0)
})

test_that("a shape leans on one partner factor and says which", {
  sim <- three(n_feats = 40, n_up = 10, n_down = 10, seed = 2)

  # The primary factor needs no partner, so it has none.
  expect_true(all(is.na(sim$truth$partner[sim$truth$pattern == "main_only"])))
  leaning <- sim$truth$pattern %in% c("additive", "interaction", "crossover",
                                      "nuisance_only")
  expect_true(all(sim$truth$partner[leaning] %in% c("sex", "time")))
  expect_true(all(is.na(sim$truth$partner[sim$truth$pattern == "none"])))

  # Whichever partner was drawn is the one the moved terms are over, and the
  # third factor is left out of it entirely.
  for (f in sim$truth$features[leaning]) {
    mate <- sim$truth$partner[sim$truth$features == f]
    got <- terms_of(sim, f)
    expect_true(all(got %in% c("treatment", mate,
                               paste("treatment", mate, sep = ":"))))
  }
})

test_that("the level profile is the one the one-factor simulator plants", {
  sim <- small(n_feats = 60, n_up = 15, n_down = 15, seed = 5,
               term_mix = c(main_only = 1))

  for (f in sim$truth$features[sim$truth$spread == "all"]) {
    cells <- cells_of(sim, f)
    d <- cells$delta[cells$sex == "male"]
    expect_equal(d[1], 0)
    expect_equal(length(unique(d[-1])), 1L)
    expect_true(d[2] != 0)
    # A main effect does not depend on the other factor, so the profile repeats.
    expect_equal(cells$delta[cells$sex == "female"], d)
  }
  for (f in sim$truth$features[sim$truth$spread == "gradient"]) {
    cells <- cells_of(sim, f)
    d <- abs(cells$delta[cells$sex == "male"])
    expect_equal(d[1], 0)
    expect_true(all(diff(d) > 0))
  }
  for (f in sim$truth$features[sim$truth$spread == "single"]) {
    cells <- cells_of(sim, f)
    d <- cells$delta[cells$sex == "male"]
    expect_equal(sum(d != 0), 1L)
  }
})

test_that("an unplanted feature is null in the strict sense", {
  sim <- small(seed = 4)
  # Not "small", exactly zero, and in every cell and every term rather than on
  # average. Everything downstream that scores a false positive rate depends on
  # there being no true effect here at all.
  null_feats <- sim$truth$features[sim$truth$direction == "none"]

  expect_true(all(sim$truth$log2fc[sim$truth$direction == "none"] == 0))
  expect_true(all(sim$truth_cell$delta[
    sim$truth_cell$features %in% null_feats] == 0))
  expect_true(all(sim$truth_term$max_abs_delta[
    sim$truth_term$features %in% null_feats] == 0))
  expect_true(!any(sim$truth_term$is_effect[
    sim$truth_term$features %in% null_feats]))
  expect_true(all(sim$truth_contrast$delta[
    sim$truth_contrast$features %in% null_feats] == 0))
  expect_true(!any(sim$truth_contrast$is_diff[
    sim$truth_contrast$features %in% null_feats]))
})

test_that("the reference cell is the zero everything else is measured from", {
  sim <- small(n_feats = 20, n_up = 5, n_down = 5, seed = 1)

  expect_true(all(sim$truth_cell$delta[sim$truth_cell$is_ref] == 0))
  # Which makes the centre of the reference cell the baseline exactly.
  expect_equal(sim$truth_cell$center[sim$truth_cell$is_ref],
               sim$truth$baseline)
})

test_that("the planted magnitudes stay inside deg_log2fc", {
  # Read on the shape whose extreme cell carries the whole magnitude and nothing
  # else. The other shapes add terms together, so their extreme cell is a sum of
  # magnitudes rather than one of them.
  sim <- small(n_feats = 30, n_up = 6, n_down = 6, deg_log2fc = c(1.5, 2),
               term_mix = c(main_only = 1), seed = 5)
  planted <- abs(sim$truth$log2fc[sim$truth$direction != "none"])

  expect_true(all(planted >= 1.5))
  expect_true(all(planted <= 2))
})

test_that("log2fc points the way direction does, crossover aside", {
  sim <- small(n_feats = 40, n_up = 10, n_down = 10, seed = 2)
  # Under crossover the primary factor rises at one level of the partner and
  # falls at another, so which cell is furthest from the reference is a fact
  # about the shape rather than about the sign of the magnitude.
  ordinary <- sim$truth$pattern %in% c("main_only", "additive", "interaction",
                                       "nuisance_only")

  expect_gt(sum(ordinary), 0L)
  expect_true(all(sim$truth$log2fc[ordinary & sim$truth$direction == "up"] > 0))
  expect_true(all(sim$truth$log2fc[ordinary &
                                     sim$truth$direction == "down"] < 0))
  # And it is the delta of the cell the answer names.
  for (f in sim$truth$features[!sim$truth$extreme_tied]) {
    row <- sim$truth[sim$truth$features == f, ]
    cells <- cells_of(sim, f)
    expect_equal(row$log2fc, cells$delta[which.max(abs(cells$delta))])
  }
})

test_that("extreme_cell names a cell only when one stands out", {
  sim <- small(n_feats = 30, n_up = 8, n_down = 8, seed = 5)

  expect_true(all(sim$truth$extreme_tied[sim$truth$pattern == "none"]))
  expect_true(all(is.na(sim$truth$extreme_cell[sim$truth$pattern == "none"])))
  expect_true(all(!is.na(sim$truth$extreme_cell[sim$truth$pattern != "none"])))
  # The label is the levels of the cell, joined, so it can be read back against
  # truth_cell without a lookup table.
  named <- sim$truth[!sim$truth$extreme_tied, ]
  for (i in seq_len(nrow(named))) {
    cells <- cells_of(sim, named$features[i])
    at <- which.max(abs(cells$delta))
    expect_identical(named$extreme_cell[i],
                     paste(cells$treatment[at], cells$sex[at], sep = "."))
  }
})

test_that("truth_contrast follows the post-hoc pair order and direction", {
  sim <- small(seed = 1)
  pairs <- sa_level_pairs(sim$args$factor_lv$treatment)
  first <- sim$truth_contrast[sim$truth_contrast$features == "prot_1", ]

  # The marginal contrasts of the first factor lead, in the pair order the
  # post-hoc tables use.
  expect_identical(first$contrast[1:6], pairs$contrast)
  expect_identical(first$group1[1:6], pairs$group1)
  expect_true(all(is.na(first$stratum[1:6])))
  # Then one block per combination of the other factors, which is the simple
  # effect the marginal one may be hiding.
  expect_identical(unique(first$stratum[first$factor == "treatment"]),
                   c(NA, "male", "female"))
  expect_identical(unique(first$stratum[first$factor == "sex"]),
                   c(NA, "control", "treat_A", "treat_B", "treat_C"))
  # Six pairs of treatments marginally and in each of two strata, plus one pair
  # of sexes marginally and in each of four strata.
  expect_equal(nrow(first), 6L * 3L + 1L * 5L)

  # A post-hoc estimate is `group1 - group2` with the reference being subtracted,
  # so the answer for a feature the treatments raised is positive here just as
  # its log2fc is. Getting this backwards would score every call wrong.
  up <- small(seed = 1, n_up = 6, n_down = 0,
              term_mix = c(main_only = 1))
  f <- up$truth$features[up$truth$direction == "up"][1]
  con <- up$truth_contrast[up$truth_contrast$features == f, ]
  expect_gt(con$delta[con$contrast == "treat_C - control" &
                        is.na(con$stratum)], 0)
  expect_identical(con$is_diff, con$delta != 0)
})

test_that("a marginal contrast and a simple effect part company", {
  # Under crossover the primary factor moves in opposite directions at the two
  # levels of the partner, so it has simple effects everywhere and no marginal
  # effect anywhere. This is the pair of statements the two kinds of row exist
  # to keep apart.
  sim <- small(n_feats = 20, n_up = 5, n_down = 5, seed = 2,
               term_mix = c(crossover = 1))
  f <- sim$truth$features[sim$truth$pattern == "crossover"][1]
  con <- sim$truth_contrast[sim$truth_contrast$features == f &
                              sim$truth_contrast$factor == "treatment", ]

  expect_true(all(con$delta[is.na(con$stratum)] == 0))
  expect_true(!any(con$is_diff[is.na(con$stratum)]))
  expect_true(any(con$is_diff[!is.na(con$stratum)]))
  # And they reverse: what rises at male falls at female.
  male <- con$delta[con$stratum == "male" & con$contrast == "treat_A - control"]
  female <- con$delta[con$stratum == "female" &
                        con$contrast == "treat_A - control"]
  expect_equal(male, -female)
})

test_that("a within factor is recorded on the terms it is part of", {
  sim <- three(seed = 1, within = "time")
  flags <- unique(sim$truth_term[c("terms", "is_within")])

  expect_true(flags$is_within[flags$terms == "time"])
  expect_true(flags$is_within[flags$terms == "treatment:time"])
  expect_true(flags$is_within[flags$terms == "treatment:sex:time"])
  expect_false(flags$is_within[flags$terms == "treatment"])
  expect_false(flags$is_within[flags$terms == "treatment:sex"])
  # Nothing is within when nothing was named.
  expect_true(!any(three(seed = 1)$truth_term$is_within))

  # The subject offset is real, so the answer records the spread it was drawn
  # with rather than leaving it to be guessed at.
  expect_true(all(!is.na(sim$truth$sd_subject)))
  expect_true(all(is.na(small(seed = 1)$truth$sd_subject)))
})

test_that("the subject offset is a subject effect and not noise", {
  # One offset per subject and feature, reused across every condition. If it were
  # drawn per row the within-subject correlation would be gone and a repeated
  # design would have nothing a between-subject one lacks.
  sim <- small(n_feats = 1, n_up = 0, n_down = 0, n_per_cell = 30,
               factor_lv = list(treatment = c("control", "treat_A"),
                                time = c("T0", "T1")),
               within = "time", subject_sd = c(8, 8), cell_sd = c(0.2, 0.2),
               ref_sd = c(0.2, 0.2), seed = 1)
  y <- sim$args$data[[1]]
  within_subject <- stats::var(tapply(y, sim$args$id, diff))
  between_subject <- stats::var(tapply(y, sim$args$id, mean))

  expect_lt(within_subject, between_subject / 10)
})

test_that("no feature is planted when none is asked for", {
  # sample() on an empty selection is the classic place for this to return
  # everything instead of nothing.
  sim <- small(seed = 1, n_up = 0, n_down = 0)
  expect_true(all(sim$truth$direction == "none"))
  expect_true(all(sim$truth$pattern == "none"))
  expect_true(all(sim$truth$spread == "none"))
  expect_true(all(sim$truth_cell$delta == 0))
  expect_true(!any(sim$truth_term$is_effect))

  up_only <- small(seed = 1, n_up = 2, n_down = 0)
  expect_equal(sum(up_only$truth$direction == "up"), 2L)
  expect_equal(sum(up_only$truth$direction == "down"), 0L)
})

test_that("n_up and n_down follow n_feats instead of being fixed counts", {
  # Called without n_up or n_down, so the defaults are the thing under test. A
  # fixed count would leave a ten feature simulation with nowhere to plant.
  planted <- function(n_feats) {
    sim <- simulate_factorial_groups(n_feats = n_feats, n_per_cell = 2,
                                     seed = 1)
    c(sum(sim$truth$direction == "up"), sum(sim$truth$direction == "down"))
  }

  expect_equal(planted(10), c(2L, 2L))
  expect_equal(planted(20), c(3L, 3L))
  expect_equal(planted(100), c(15L, 15L))
})

test_that("feat_prefix names the features", {
  sim <- small(seed = 1, feat_prefix = "gene")
  expect_identical(sim$args$feats[1:2], c("gene_1", "gene_2"))
  expect_identical(sim$truth$features, sim$args$feats)
})

test_that("a seed makes the draw reproducible without stealing the stream", {
  expect_equal(small(seed = 42), small(seed = 42))
  expect_false(isTRUE(all.equal(small(seed = 1)$args$data,
                                small(seed = 2)$args$data)))

  # Seeding inside the function must not reseed the caller, or every draw after
  # a simulation would silently repeat.
  set.seed(99)
  before <- stats::runif(3)
  set.seed(99)
  invisible(small(seed = 7))
  expect_equal(stats::runif(3), before)
})

test_that("without a seed the draw follows the caller's stream", {
  set.seed(11)
  first <- small()
  set.seed(11)
  expect_equal(small(), first)
})

test_that("the factors can be crossed any number of ways", {
  sim <- three(seed = 1)

  # Three factors give three main effects, three pairs and one triple.
  expect_equal(nrow(sim$truth_term), 12L * 7L)
  expect_identical(unique(sim$truth_term$terms),
                   c("treatment", "sex", "time", "treatment:sex",
                     "treatment:time", "sex:time", "treatment:sex:time"))
  expect_equal(nrow(sim$truth_cell), 12L * 3L * 2L * 3L)
  expect_equal(nrow(sim$args$data), 3L * 2L * 3L * 4L)

  # And two are the fewest there can be: one factor is a job for
  # simulate_multiple_groups().
  expect_error(small(factor_lv = list(treatment = c("a", "b", "c"))),
               "at least two crossed factors")
})

test_that("impossible or malformed arguments are rejected", {
  expect_error(small(n_feats = 10, n_up = 6, n_down = 6), "more features")
  expect_error(small(n_feats = 2.5), "whole number")
  expect_error(small(n_per_cell = 1), "must be in")
  expect_error(small(n_per_cell = c(4, 1, 4, 4, 4, 4, 4, 4)), "must be in")
  expect_error(small(n_per_cell = c(4, 4)), "one size per combination")
  expect_error(small(factor_lv = c("a", "b")), "named list")
  expect_error(small(factor_lv = list(a = c("x", "y"), c("p", "q"))),
               "named list")
  expect_error(small(factor_lv = list(a = c("x", "y"), a = c("p", "q"))),
               "named list")
  expect_error(small(factor_lv = list(treatment = c("a", "b"), sex = "male")),
               "at least two distinct")
  expect_error(small(factor_lv = list(treatment = c("a", "a"),
                                      sex = c("m", "f"))),
               "at least two distinct")
  # The answer tables give each factor a column beside their own, so a factor
  # named after one of those would make a table with two columns of one name.
  expect_error(small(factor_lv = list(treatment = c("a", "b"),
                                      delta = c("m", "f"))),
               "already use as columns")
  expect_error(small(expr_range = c(12, 2)), "increasing")
  expect_error(small(cell_sd = c(-1, 2)), "must not go below")
  expect_error(small(deg_log2fc = 1), "length 2")
  expect_error(small(interaction_scale = 0), "must be in")
  expect_error(small(feat_prefix = ""), "non-empty")
})

test_that("a malformed within is rejected", {
  expect_error(small(within = "nope"), "does not hold")
  expect_error(small(within = c("sex", "sex")), "distinct")
  expect_error(small(within = 1), "must be NULL")
  # Naming every factor is a fully repeated design, not a mistake.
  expect_silent(small(within = c("treatment", "sex"), seed = 1))
})

test_that("a malformed term_mix is rejected", {
  expect_error(small(term_mix = c(1, 1)), "named numeric vector")
  expect_error(small(term_mix = c(main_only = 1, nope = 1)), "unknown shape")
  expect_error(small(term_mix = c(main_only = -1, additive = 1)),
               "not be negative")
  expect_error(small(term_mix = c(main_only = 0, crossover = 0)),
               "at least one positive weight")
  # The one-factor mix is checked the same way and still names its own shapes.
  expect_error(small(pattern_mix = c(all = 1, nope = 1)), "unknown shape")
  expect_error(small(pattern_mix = c(main_only = 1)), "unknown shape")
})
