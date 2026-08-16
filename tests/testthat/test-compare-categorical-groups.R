# A contingency table is one question, so what is pinned here is the layout of
# that question -- which tests run, which measures exist, the direction of the
# odds ratio -- and agreement with the engines the kernels wrap.
#
# The other thing pinned here is that `design$null` is load bearing. `expected`
# is not a property of a table but of a table and a claim about it, so the tests
# below check that each design names its own claim and that the residuals are
# read under it. The matched 2 x 2 case has an exact identity to hold: the
# squared Pearson residuals of the discordant cells sum to McNemar's uncorrected
# statistic, which is what makes the cell table and the p-value beside it
# demonstrably about one hypothesis.

test_that("an independent design runs chi-square beside Fisher", {
  res <- sa_categorical_fixture()

  expect_s3_class(res, c("sa_categorical", "sa_result"), exact = TRUE)
  expect_false(inherits(res, "sa_comparison"))
  expect_identical(res$analysis, "categorical_comparison")
  expect_identical(res$variables, c("smoker", "grade"))
  expect_named(res$tests, c("chisq_test", "fisher_test"))
  expect_setequal(names(res$tests), names(res$test_info))
  expect_equal(nrow(res$tests$chisq_test), 1L)
  expect_true(all(sa_categorical_test_columns() %in% names(res$tests$chisq_test)))
  expect_false("pval_adj" %in% names(res$tests$chisq_test))
})

test_that("each design names the null hypothesis it tested", {
  expect_identical(sa_categorical_fixture()$design$null, "independence")
  expect_identical(sa_categorical_matched_fixture()$design$null, "symmetry")
  expect_identical(sa_categorical_repeated_fixture()$design$null,
                   "marginal_homogeneity")
  expect_true(all(c("independence", "symmetry", "marginal_homogeneity") %in%
                    sa_categorical_nulls()))
})

test_that("the contract refuses a null it does not define", {
  res <- sa_categorical_fixture()
  design <- res$design
  design$null <- "agreement"
  expect_error(
    sa_new_categorical("categorical_comparison", res$variables, design,
                       res$parameters, res$cells, res$tests, res$test_info,
                       res$association),
    "design\\$null"
  )
})

test_that("the cell table is the long form of the table", {
  res <- sa_categorical_fixture()
  counts <- as.table(res)

  expect_true(all(sa_categorical_cell_columns() %in% names(res$cells)))
  expect_equal(nrow(res$cells), prod(res$design$dim))
  expect_equal(sum(res$cells$observed), res$design$n_used)
  expect_equal(sum(res$cells$prop_total), 1, tolerance = 1e-12)
  expect_identical(as.numeric(counts), res$cells$observed)
})

test_that("the table is a method rather than a slot", {
  res <- sa_categorical_fixture()

  expect_null(res$counts)
  counts <- as.table(res)
  expect_s3_class(counts, "table")
  expect_identical(dim(counts), res$design$dim)
  expect_identical(names(dimnames(counts)),
                   c(res$design$row_var, res$design$col_var))
  # The table the tests were actually run on, so it feeds the engine directly
  # and comes back with the statistic the object already carries.
  theirs <- suppressWarnings(stats::chisq.test(counts, correct = TRUE))
  expect_equal(res$tests$chisq_test$statistic, unname(theirs$statistic))
})

test_that("the chi-square statistic matches stats::chisq.test", {
  res <- sa_categorical_fixture()
  theirs <- suppressWarnings(stats::chisq.test(as.table(res), correct = TRUE))

  expect_equal(res$tests$chisq_test$statistic, unname(theirs$statistic))
  expect_equal(res$tests$chisq_test$df, unname(theirs$parameter))
  expect_equal(res$tests$chisq_test$pval, unname(theirs$p.value))
})

test_that("Fisher's p-value matches stats::fisher.test", {
  res <- sa_categorical_fixture()
  theirs <- stats::fisher.test(as.table(res))

  expect_equal(res$tests$fisher_test$pval, unname(theirs$p.value))
  expect_true(is.na(res$tests$fisher_test$statistic))
})

test_that("independence expects the product of the margins", {
  res <- sa_categorical_fixture()
  counts <- as.table(res)
  theirs <- suppressWarnings(stats::chisq.test(counts, correct = TRUE))

  expect_equal(res$cells$expected, as.numeric(theirs$expected))
  expect_equal(res$cells$residual, as.numeric(theirs$residuals))
  expect_equal(res$cells$std_residual, as.numeric(theirs$stdres))
  # Pearson residuals square and sum to the uncorrected statistic, which is what
  # makes them a statement about how it was made up.
  expect_equal(sum(res$cells$residual^2),
               sum((counts - theirs$expected)^2 / theirs$expected))
})

test_that("a 2 x 2 table reports phi and an odds ratio", {
  res <- compare_categorical_groups(
    sa_categorical_frame(),
    category_lv = list(smoker = c("n", "y"), grade = c("high", "low")),
    diagnose = FALSE
  )

  expect_identical(res$association$measure,
                   c("cramers_v", "contingency_coefficient",
                     "phi_coefficient", "odds_ratio"))
  or <- res$association$estimate[res$association$measure == "odds_ratio"]
  expect_gt(or, 1)
  expect_true(all(is.finite(res$association$lower_conf[
    res$association$measure == "odds_ratio"])))
})

test_that("a larger table leaves phi and the odds ratio out", {
  res <- sa_categorical_fixture()
  expect_identical(res$association$measure,
                   c("cramers_v", "contingency_coefficient"))
})

test_that("control_label inverts the odds ratio of one variable", {
  lv <- list(smoker = c("n", "y"), grade = c("low", "high"))
  a <- compare_categorical_groups(sa_categorical_frame(), category_lv = lv,
                                  diagnose = FALSE)
  b <- compare_categorical_groups(sa_categorical_frame(), category_lv = lv,
                                  control_label = c(smoker = "y"),
                                  diagnose = FALSE)

  or_a <- a$association$estimate[a$association$measure == "odds_ratio"]
  or_b <- b$association$estimate[b$association$measure == "odds_ratio"]
  expect_equal(or_a, 1 / or_b, tolerance = 1e-10)
  expect_identical(a$design$category_lv$smoker, c("n", "y"))
  expect_identical(b$design$category_lv$smoker, c("y", "n"))
})

test_that("category_lv drops levels it leaves out", {
  res <- compare_categorical_groups(
    sa_categorical_frame(),
    category_lv = list(smoker = c("n", "y"), grade = c("low", "high")),
    diagnose = FALSE
  )
  expect_equal(res$design$n_dropped, 40L)
  expect_identical(res$design$dim, c(2L, 2L))
})

test_that("a zero cell is reported by the scenario, not by the kernel", {
  # A table on the diagonal: the odds ratio needs the Haldane-Anscombe
  # correction, and who says so is the point. No kernel in the package speaks.
  d <- data.frame(u = rep(c("a", "b"), each = 10),
                  v = rep(c("x", "y"), each = 10))
  expect_message(compare_categorical_groups(d, diagnose = FALSE),
                 "Haldane-Anscombe")
  quiet <- suppressMessages(compare_categorical_groups(d, diagnose = FALSE))
  or <- quiet$association$estimate[quiet$association$measure == "odds_ratio"]
  expect_true(is.finite(or))
  # The correction touched the measure and not the table.
  expect_identical(as.numeric(as.table(quiet)), c(10, 0, 0, 10))
})

test_that("max_levels refuses a measurement read as a category", {
  set.seed(4)
  d <- data.frame(a = rnorm(50), b = rnorm(50))
  expect_error(compare_categorical_groups(d), "max_levels")
  expect_error(compare_categorical_groups(d), "a \\(50\\)")
  expect_error(compare_categorical_groups(d, max_levels = 1), "max_levels")
})

test_that("max_levels is checked against the levels actually used", {
  set.seed(5)
  d <- data.frame(
    dose = c(rep(c("low", "high"), each = 20), round(rnorm(20), 6)),
    resp = rep(c("n", "y"), 30),
    stringsAsFactors = FALSE
  )
  # Naming two levels of the many-valued column is a way through, and the rows
  # at the rest are dropped rather than counted as levels.
  res <- suppressMessages(compare_categorical_groups(
    d,
    category_lv = list(dose = c("low", "high"), resp = c("n", "y")),
    diagnose = FALSE
  ))
  expect_identical(res$design$dim, c(2L, 2L))
  expect_equal(res$design$n_dropped, 20L)
  expect_identical(res$parameters$max_levels, 20L)
})

test_that("an unenumerable Fisher is an NA rather than a failed call", {
  # A 3 x 3 table of this many observations has more tables with the observed
  # margins than the network algorithm's workspace holds. Losing the chi-square
  # result over that would be the more expensive failure.
  sim <- simulate_categorical_groups(
    n_samples   = 300,
    category_lv = list(dose = c("low", "mid", "high"),
                       response = c("none", "partial", "full")),
    assoc       = 0.5,
    seed        = 3
  )
  expect_message(do.call(compare_categorical_groups,
                         c(sim$args, list(diagnose = FALSE))),
                 "could not enumerate")
  res <- suppressMessages(do.call(compare_categorical_groups,
                                  c(sim$args, list(diagnose = FALSE))))

  expect_true(is.na(res$tests$fisher_test$pval))
  expect_true(is.finite(res$tests$chisq_test$pval))
  expect_false("enumerated" %in% names(res$tests$fisher_test))

  # And the Monte Carlo variant is what answers there.
  simulated <- do.call(compare_categorical_groups,
                       c(sim$args, list(simulate_p_value = TRUE,
                                        diagnose = FALSE)))
  expect_true(is.finite(simulated$tests$fisher_test$pval))
})

test_that("three variables without paired is an error", {
  d <- data.frame(a = c("x", "y"), b = c("p", "q"), c = c("u", "v"))
  expect_error(compare_categorical_groups(d),
               "exactly two variables")
})

test_that("a matched 2 x 2 runs McNemar's test", {
  res <- compare_categorical_groups(sa_categorical_matched_frame(),
                                    paired = TRUE, exact = FALSE,
                                    correct = FALSE, diagnose = FALSE)

  expect_named(res$tests, "mcnemar_test")
  theirs <- stats::mcnemar.test(as.table(res), correct = FALSE)
  expect_equal(res$tests$mcnemar_test$statistic, unname(theirs$statistic))
  expect_equal(res$tests$mcnemar_test$pval, unname(theirs$p.value))
  expect_identical(res$association$measure,
                   c("odds_ratio_paired", "risk_difference_paired", "cohens_g"))
})

test_that("exact_used is a setting rather than a column of the test table", {
  res <- sa_categorical_matched_fixture()

  expect_false("exact_used" %in% names(res$tests$mcnemar_test))
  expect_true("n_discordant" %in% names(res$tests$mcnemar_test))
  expect_true(res$parameters$exact)
  expect_false(compare_categorical_groups(sa_categorical_matched_frame(),
                                          paired = TRUE, exact = FALSE,
                                          diagnose = FALSE)$parameters$exact)
})

test_that("McNemar's exact branch matches binom.test on the discordant cells", {
  res <- sa_categorical_matched_fixture(exact = TRUE)
  counts <- as.table(res)
  b <- counts[1, 2]
  n_disc <- counts[1, 2] + counts[2, 1]

  expect_equal(res$tests$mcnemar_test$pval,
               stats::binom.test(b, n_disc, p = 0.5)$p.value)
  expect_true(res$parameters$exact)
})

test_that("symmetry expects a cell at the average of it and its transpose", {
  res <- sa_categorical_matched_fixture()
  counts <- as.table(res)
  expected <- (counts + t(counts)) / 2

  expect_equal(res$cells$expected, as.numeric(expected))
  # The diagonal is expected at exactly what it holds, so it carries no residual
  # at all. That is the concordant pairs dropping out of the comparison, which is
  # the same fact McNemar's statistic rests on.
  diagonal <- res$cells$row_level == res$cells$col_level
  expect_identical(res$cells$residual[diagonal], c(0, 0))
})

test_that("the discordant residuals square to McNemar's statistic exactly", {
  res <- sa_categorical_matched_fixture()
  counts <- as.table(res)
  b <- counts[1, 2]
  c_ <- counts[2, 1]
  off <- res$cells$residual[res$cells$row_level != res$cells$col_level]

  expect_equal(sum(off^2), (b - c_)^2 / (b + c_))
  # Which is the uncorrected statistic the engine reports, so the cell table and
  # the test are one piece of arithmetic read two ways.
  theirs <- stats::mcnemar.test(counts, correct = FALSE)
  expect_equal(sum(off^2), unname(theirs$statistic))
})

test_that("the standardized residual is NA under symmetry", {
  res <- sa_categorical_matched_fixture()

  expect_true(all(is.na(res$cells$std_residual)))
  # And present where its variance correction is defined.
  expect_false(any(is.na(sa_categorical_fixture()$cells$std_residual)))
})

test_that("three binary conditions run Cochran's Q", {
  res <- sa_categorical_repeated_fixture()

  expect_named(res$tests, "cochran_q")
  expect_identical(res$association$measure, "kendalls_w")
  expect_identical(res$design$row_var, "condition")
  expect_identical(res$design$col_var, "response")
  expect_identical(rownames(as.table(res)), c("t1", "t2", "t3"))
  # Every subject is measured under every condition, so the table counts more
  # than the rows the design used.
  expect_equal(sum(res$cells$observed), res$design$n_used * 3L)
})

test_that("Cochran's Q matches the formula the registry records", {
  d <- data.frame(
    t1 = c("n", "n", "n", "y", "y", "y", "n", "y"),
    t2 = c("n", "y", "y", "n", "y", "y", "y", "n"),
    t3 = c("y", "n", "y", "n", "n", "y", "y", "y")
  )
  res <- compare_categorical_groups(
    d, category_lv = list(t1 = c("n", "y"), t2 = c("n", "y"), t3 = c("n", "y")),
    paired = TRUE, diagnose = FALSE
  )
  mat <- matrix(c(0, 0, 0, 1, 1, 1, 0, 1,
                  0, 1, 1, 0, 1, 1, 1, 0,
                  1, 0, 1, 0, 0, 1, 1, 1),
                nrow = 8L)
  k <- 3
  col_n <- colSums(mat)
  row_n <- rowSums(mat)
  total <- sum(mat)
  q_stat <- k * (k - 1) * sum((col_n - total / k)^2) / (k * total - sum(row_n^2))
  expect_equal(res$tests$cochran_q$statistic, q_stat, tolerance = 1e-12)
  expect_equal(res$tests$cochran_q$pval,
               stats::pchisq(q_stat, df = 2, lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("marginal homogeneity expects every condition at the pooled rate", {
  res <- sa_categorical_repeated_fixture()
  counts <- as.table(res)
  pooled <- colSums(counts) / sum(counts)

  # One expectation per condition, and it is the same one for all of them: that
  # is what the null says and what the mosaic draws as a line.
  by_row <- matrix(res$cells$expected, nrow = 3L)
  for (i in seq_len(3L)) {
    expect_equal(by_row[i, ] / sum(by_row[i, ]), unname(pooled))
  }
})

test_that("a matched design with more than two levels is refused", {
  d <- data.frame(a = c("low", "mid", "high"),
                  b = c("mid", "high", "low"))
  expect_error(compare_categorical_groups(d, paired = TRUE),
               "binary conditions")
})

test_that("each design attaches the approximation rule the registry records", {
  expect_identical(
    compare_categorical_groups(sa_categorical_frame())$diagnostics$rule,
    "expected_count_min"
  )
  expect_identical(
    suppressMessages(compare_categorical_groups(
      sa_categorical_matched_frame(), paired = TRUE))$diagnostics$rule,
    "discordant_pair_count"
  )
  expect_identical(
    compare_categorical_groups(sa_categorical_repeated_frame(),
                               paired = TRUE)$diagnostics$rule,
    "sample_size_repeated"
  )
})

test_that("print names the null hypothesis the result is about", {
  res <- sa_categorical_fixture()
  expect_output(print(res), "categorical_comparison")
  expect_output(print(res), "smoker")
  expect_output(print(res), "null")
  expect_output(print(res), "independence")
  expect_output(print(res), "chisq_test")
  expect_output(print(res), "cramers_v")
  expect_invisible(print(res))

  expect_output(print(sa_categorical_matched_fixture()), "symmetry")
  expect_output(print(sa_categorical_repeated_fixture()),
                "marginal homogeneity")
})
