# The value of a simulator is that the answer is known, so what is pinned here
# is that the answer it reports is the answer it planted: the strict nullity of
# assoc = 0, the merge keys against a comparison, and that args is a call.
#
# A matched design has a second null to plant against, since McNemar's test is
# about symmetry rather than independence, so `truth_cell` carries the symmetric
# share as well and the tests below check that it is the share the comparison
# expects.

test_that("the returned arguments are the ones compare_categorical_groups takes", {
  sim <- simulate_categorical_groups(n_samples = 80, seed = 1)

  expect_named(sim, c("args", "truth", "truth_cell"))
  expect_true(all(names(sim$args) %in% names(formals(compare_categorical_groups))))
  expect_identical(names(sim$args$data), names(sim$args$category_lv))
  expect_false(sim$args$paired)
  expect_identical(sim$truth$n_samples, 80L)
})

test_that("assoc = 0 is null in the strict sense", {
  sim <- simulate_categorical_groups(n_samples = 200, assoc = 0, seed = 2)

  expect_equal(sim$truth$assoc, 0)
  expect_equal(sim$truth$cramers_v, 0)
  expect_true(all(sim$truth_cell$lift == 1))
  expect_equal(sim$truth_cell$p_planted, sim$truth_cell$p_independent)
})

test_that("a planted association is recovered by both tests", {
  sim <- simulate_categorical_groups(n_samples = 400, assoc = 0.5, seed = 3)
  fit <- do.call(compare_categorical_groups,
                 c(sim$args, list(diagnose = FALSE)))

  expect_lt(fit$tests$chisq_test$pval, 0.05)
  expect_lt(fit$tests$fisher_test$pval, 0.05)
  scored <- merge(fit$cells, sim$truth_cell, by = c("row_level", "col_level"))
  expect_equal(nrow(scored), nrow(fit$cells))
})

test_that("the planted expectation is the one the comparison holds the table to", {
  sim <- simulate_categorical_groups(n_samples = 4000, assoc = 0, seed = 7)
  fit <- do.call(compare_categorical_groups,
                 c(sim$args, list(diagnose = FALSE)))
  scored <- merge(fit$cells, sim$truth_cell, by = c("row_level", "col_level"))

  expect_identical(fit$design$null, "independence")
  # A large table drawn at exact independence, so the expected counts the
  # comparison computes from the observed margins land near the planted ones.
  expect_equal(scored$expected, scored$expected_n, tolerance = 0.05)
})

test_that("the seed is restored and a seeded call is reproducible", {
  stream_before <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
  a <- simulate_categorical_groups(n_samples = 60, seed = 11)
  stream_after <- get(".Random.seed", envir = globalenv())
  b <- simulate_categorical_groups(n_samples = 60, seed = 11)

  expect_identical(a$args$data, b$args$data)
  if (!is.null(stream_before)) {
    expect_identical(stream_after, stream_before)
  }
})

test_that("a matched simulation plants the paired odds ratio as a ratio", {
  sim <- simulate_categorical_groups(n_samples = 200, paired = TRUE,
                                     discordance = c(0.3, 0.1), seed = 4)

  expect_true(sim$args$paired)
  expect_equal(sim$truth$odds_ratio_paired, 0.3 / 0.1)
  expect_equal(sim$truth$risk_difference_paired, 0.5 * 0.3 - 0.5 * 0.1)
  fit <- do.call(compare_categorical_groups,
                 c(sim$args, list(diagnose = FALSE)))
  expect_named(fit$tests, "mcnemar_test")
})

test_that("a matched simulation plants the symmetric share as well", {
  sim <- simulate_categorical_groups(n_samples = 200, paired = TRUE,
                                     discordance = c(0.3, 0.1), seed = 4)

  expect_true(all(c("p_symmetric", "expected_symmetry_n") %in%
                    names(sim$truth_cell)))
  expect_equal(sum(sim$truth_cell$p_symmetric), 1)
  expect_equal(sim$truth_cell$expected_symmetry_n,
               200 * sim$truth_cell$p_symmetric)
  # Symmetry says nothing about the concordant pairs, so the diagonal is already
  # where the null puts it and only the discordant cells hold a departure.
  diagonal <- sim$truth_cell$row_level == sim$truth_cell$col_level
  expect_equal(sim$truth_cell$p_symmetric[diagonal],
               sim$truth_cell$p_planted[diagonal])
  expect_false(isTRUE(all.equal(sim$truth_cell$p_symmetric[!diagonal],
                                sim$truth_cell$p_planted[!diagonal])))
})

test_that("the symmetric truth is recovered by the cells of a matched fit", {
  sim <- simulate_categorical_groups(n_samples = 6000, paired = TRUE,
                                     discordance = c(0.3, 0.1), seed = 8)
  fit <- do.call(compare_categorical_groups,
                 c(sim$args, list(diagnose = FALSE)))
  scored <- merge(fit$cells, sim$truth_cell, by = c("row_level", "col_level"))

  expect_identical(fit$design$null, "symmetry")
  expect_equal(nrow(scored), nrow(fit$cells))
  # An absolute margin rather than a relative one: the smallest planted share is
  # 0.10, and a share that small carries a sampling error of a few per cent of
  # itself even at 6000 draws.
  expect_lt(max(abs(scored$expected / 6000 - scored$p_symmetric)), 0.02)
})

test_that("equal discordance is the strict null of a matched design", {
  sim <- simulate_categorical_groups(n_samples = 200, paired = TRUE,
                                     discordance = c(0.2, 0.2), seed = 5)
  expect_equal(sim$truth$odds_ratio_paired, 1)
  expect_equal(sim$truth$risk_difference_paired, 0)
  expect_equal(sim$truth$cohens_g, 0)
  # Which is symmetry exactly, so the planted table is its own transpose and
  # every cell already sits where the null puts it.
  expect_equal(sim$truth_cell$p_symmetric, sim$truth_cell$p_planted)
})

test_that("three matched conditions plant a climb in the response rate", {
  sim <- simulate_categorical_groups(
    n_samples = 150,
    category_lv = list(t1 = c("n", "y"), t2 = c("n", "y"), t3 = c("n", "y")),
    paired = TRUE,
    discordance = c(0.35, 0.10),
    seed = 6
  )
  expect_gt(sim$truth$rate_last, sim$truth$rate_first)
  fit <- do.call(compare_categorical_groups,
                 c(sim$args, list(diagnose = FALSE)))
  expect_named(fit$tests, "cochran_q")
  expect_identical(fit$design$null, "marginal_homogeneity")
})

test_that("a repeated truth counts the whole table rather than the subjects", {
  sim <- simulate_categorical_groups(
    n_samples = 150,
    category_lv = list(t1 = c("n", "y"), t2 = c("n", "y"), t3 = c("n", "y")),
    paired = TRUE,
    discordance = c(0.35, 0.10),
    seed = 6
  )
  fit <- do.call(compare_categorical_groups,
                 c(sim$args, list(diagnose = FALSE)))

  # Every subject is measured under every condition, so 150 subjects make 450
  # observations, and the planted counts have to be on the same scale as the
  # cells they are merged against.
  expect_equal(sum(sim$truth_cell$expected_n), 150 * 3)
  expect_equal(sum(fit$cells$observed), 150 * 3)
  # And marginal homogeneity is the same arithmetic as independence there, so
  # the null needs no column of its own.
  expect_false("p_symmetric" %in% names(sim$truth_cell))
  scored <- merge(fit$cells, sim$truth_cell, by = c("row_level", "col_level"))
  expect_equal(nrow(scored), nrow(fit$cells))
})

test_that("an argument of the other design is a warning naming both", {
  expect_warning(
    simulate_categorical_groups(n_samples = 60, paired = TRUE, assoc = 0.4,
                                seed = 9),
    "matched design reads discordance"
  )
  expect_warning(
    simulate_categorical_groups(n_samples = 60, paired = TRUE, assoc = 0.4,
                                seed = 9),
    "given for assoc were ignored"
  )
  expect_warning(
    simulate_categorical_groups(n_samples = 60, discordance = c(0.3, 0.2),
                                seed = 9),
    "cross-classified design reads margins, assoc, pattern"
  )
  # An argument the design does read is silent, `pattern` included, which
  # `match.arg()` would otherwise make look supplied.
  expect_silent(simulate_categorical_groups(n_samples = 60, assoc = 0.4,
                                            pattern = "gradient", seed = 9))
  expect_silent(simulate_categorical_groups(n_samples = 60, paired = TRUE,
                                            discordance = c(0.3, 0.2), seed = 9))
})
