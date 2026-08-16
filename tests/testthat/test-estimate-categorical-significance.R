# The point of the cell reading is that it does not change shape with the table,
# so most of what follows is run against a 2 x 2, a 2 x 3 and something larger at
# once rather than against whichever one happened to be convenient.

sa_cat_2x2 <- function() {
  # Naming two of the three grades is what makes the table 2 x 2, and dropping
  # the third is worth a message that is not this file's subject.
  suppressMessages(compare_categorical_groups(
    sa_categorical_frame(), diagnose = FALSE,
    category_lv = list(smoker = c("n", "y"), grade = c("low", "high"))
  ))
}

sa_cat_4x3 <- function() {
  set.seed(4)
  # Fisher's exact test cannot enumerate a table this size and says so, which is
  # likewise not what is under test here.
  suppressMessages(compare_categorical_groups(
    data.frame(site  = sample(c("a", "b", "c", "d"), 300, replace = TRUE),
               grade = sample(c("low", "mid", "high"), 300, replace = TRUE),
               stringsAsFactors = FALSE),
    diagnose = FALSE
  ))
}

sa_cat_shapes <- function() {
  list(`2x2` = sa_cat_2x2(),
       `2x3` = sa_categorical_fixture(),
       `4x3` = sa_cat_4x3())
}

sa_cell_columns <- c("row_level", "col_level", "observed", "expected", "lift",
                     "log2_lift", "std_residual", "pvalue", "adj_pvalue",
                     "is_signif")


test_that("a cell reading has the same shape whatever size the table is", {
  for (nm in names(sa_cat_shapes())) {
    res <- sa_cat_shapes()[[nm]]
    sig <- estimate_categorical_significance(res)

    expect_s3_class(sig, "sa_categorical_significance")
    expect_identical(sig$analysis_type, "categorical_comparison", info = nm)
    expect_identical(names(sig$significance), sa_cell_columns, info = nm)
    expect_identical(nrow(sig$significance), nrow(res$cells), info = nm)
    expect_type(sig$significance$is_signif, "logical")
  }
})


test_that("the cell key is the one $cells and truth_cell use", {
  res <- sa_categorical_fixture()
  sig <- estimate_categorical_significance(res)

  merged <- merge(sig$significance, res$cells,
                  by = c("row_level", "col_level"))
  expect_identical(nrow(merged), nrow(res$cells))
  expect_equal(merged$observed.x, merged$observed.y)
})


test_that("pvalue is exactly the two-sided normal tail of std_residual", {
  for (nm in names(sa_cat_shapes())) {
    tbl <- estimate_categorical_significance(
      sa_cat_shapes()[[nm]])$significance
    expect_equal(tbl$pvalue, 2 * stats::pnorm(-abs(tbl$std_residual)),
                 info = nm)
  }
})


test_that("lift is observed over expected and log2_lift is its log2", {
  for (nm in names(sa_cat_shapes())) {
    tbl <- estimate_categorical_significance(
      sa_cat_shapes()[[nm]])$significance
    expect_equal(tbl$lift, tbl$observed / tbl$expected, info = nm)
    expect_equal(tbl$log2_lift, log2(tbl$lift), info = nm)
  }
})


test_that("an unobserved cell is an infinite shortfall and clears any cutoff", {
  # A level of `grade` that only one group ever reaches, so the other group's
  # cell holds nothing while its expected count is comfortably positive.
  empty <- compare_categorical_groups(
    data.frame(smoker = rep(c("y", "n"), each = 40),
               grade  = c(rep(c("high", "low"), c(20, 20)),
                          rep(c("high", "mid"), c(20, 20))),
               stringsAsFactors = FALSE),
    diagnose = FALSE
  )
  tbl <- estimate_categorical_significance(empty)$significance
  gone <- tbl[tbl$observed == 0, ]

  expect_gt(nrow(gone), 0L)
  expect_true(all(gone$lift == 0))
  expect_true(all(gone$log2_lift == -Inf))
  expect_true(all(is.finite(gone$expected) & gone$expected > 0))
  # -Inf clears every finite magnitude cutoff, which is the point: the cell is
  # as far below what was expected as a ratio can go.
  expect_true(all(abs(gone$log2_lift) >= 99))
})


test_that("adj_pvalue is adjusted across the cells and follows adj_type", {
  res <- sa_cat_4x3()

  bh <- estimate_categorical_significance(res)$significance
  expect_equal(bh$adj_pvalue, stats::p.adjust(bh$pvalue, "BH"))

  holm <- estimate_categorical_significance(res, adj_type = "holm")$significance
  expect_equal(holm$adj_pvalue, stats::p.adjust(holm$pvalue, "holm"))
  expect_false(isTRUE(all.equal(bh$adj_pvalue, holm$adj_pvalue)))

  none <- estimate_categorical_significance(res, adj_type = "none")$significance
  expect_equal(none$adj_pvalue, none$pvalue)

  # The family is this table's cells and nothing else, so the adjustment is
  # over as many p-values as the table has cells.
  expect_identical(nrow(bh), nrow(res$cells))
})


test_that("is_signif is three-valued and NA means undecided", {
  res <- sa_categorical_fixture()
  tbl <- estimate_categorical_significance(res, log2_lift_cutoff = 0.1,
                                           pval_cutoff = 0.05)$significance

  expect_true(all(tbl$is_signif %in% c(TRUE, FALSE, NA)))
  expect_identical(
    tbl$is_signif,
    abs(tbl$log2_lift) >= 0.1 & tbl$adj_pvalue <= 0.05
  )

  # A cell whose expected count is zero has no lift to take, and the verdict
  # says so rather than deciding against it.
  undecided <- estimate_categorical_significance(res)$significance
  undecided$log2_lift[1] <- NA_real_
  expect_true(is.na(abs(undecided$log2_lift[1]) >= 1 &
                      undecided$adj_pvalue[1] <= 0.05))
})


test_that("a looser magnitude cutoff can only add cells, never remove them", {
  res <- sa_cat_4x3()
  strict <- estimate_categorical_significance(res, log2_lift_cutoff = 1)
  loose <- estimate_categorical_significance(res, log2_lift_cutoff = 0.2)

  expect_true(all(which(strict$significance$is_signif %in% TRUE) %in%
                    which(loose$significance$is_signif %in% TRUE)))
})


test_that("the planted lift is recovered from the cell reading", {
  sim <- simulate_categorical_groups(n_samples = 4000, assoc = 0.4, seed = 11)
  fit <- do.call(compare_categorical_groups,
                 c(sim$args, list(diagnose = FALSE)))
  scored <- merge(estimate_categorical_significance(fit)$significance,
                  sim$truth_cell, by = c("row_level", "col_level"))

  expect_identical(nrow(scored), nrow(sim$truth_cell))
  # `lift.x` is estimated and `lift.y` planted. Four thousand rows is enough for
  # the two to agree to within a tenth on every cell.
  expect_equal(scored$lift.x, scored$lift.y, tolerance = 0.1)
})


test_that("a table reading is one row with the association beside the p", {
  sig <- estimate_categorical_significance(sa_categorical_fixture(),
                                           by = "table")

  expect_identical(nrow(sig$significance), 1L)
  expect_identical(names(sig$significance),
                   c("measure", "estimate", "lower_conf", "upper_conf",
                     "pvalue", "is_signif"))
  # One table is one question, so there is no family and no adjusted column.
  expect_false("adj_pvalue" %in% names(sig$significance))
})


test_that("measure = auto takes the one each design defines", {
  auto <- function(res) {
    estimate_categorical_significance(res, by = "table")$significance$measure
  }

  expect_identical(auto(sa_cat_2x2()), "odds_ratio")
  expect_identical(auto(sa_categorical_fixture()), "cramers_v")
  expect_identical(auto(sa_cat_4x3()), "cramers_v")
  expect_identical(auto(sa_categorical_matched_fixture()), "odds_ratio_paired")
  expect_identical(auto(sa_categorical_repeated_fixture()), "kendalls_w")
})


test_that("a table reading reports the measure it was asked for", {
  sig <- estimate_categorical_significance(
    sa_categorical_fixture(), by = "table",
    measure = "contingency_coefficient"
  )
  res <- sa_categorical_fixture()
  row <- res$association[
    res$association$measure == "contingency_coefficient", ]

  expect_identical(sig$significance$measure, "contingency_coefficient")
  expect_equal(sig$significance$estimate, row$estimate)
})


test_that("effect_cutoff is read on the measure's own scale", {
  two <- sa_cat_2x2()
  or <- estimate_categorical_significance(two, by = "table")$significance$estimate

  # A ratio centred at 1: a cutoff of c is cleared from either side.
  expect_true(or < 1)
  expect_true(estimate_categorical_significance(
    two, by = "table", effect_cutoff = 1 / or)$significance$is_signif)
  expect_false(estimate_categorical_significance(
    two, by = "table", effect_cutoff = 100)$significance$is_signif)

  # A magnitude measure is compared against the cutoff directly.
  v <- estimate_categorical_significance(
    sa_categorical_fixture(), by = "table")$significance$estimate
  expect_true(estimate_categorical_significance(
    sa_categorical_fixture(), by = "table",
    effect_cutoff = v - 0.01)$significance$is_signif)
  expect_false(estimate_categorical_significance(
    sa_categorical_fixture(), by = "table",
    effect_cutoff = v + 0.01)$significance$is_signif)

  # A signed one is compared against its magnitude, so a fall below zero is as
  # large a departure as the same rise above it.
  g <- estimate_categorical_significance(
    sa_categorical_matched_fixture(), by = "table",
    measure = "cohens_g")$significance
  expect_true(estimate_categorical_significance(
    sa_categorical_matched_fixture(), by = "table", measure = "cohens_g",
    effect_cutoff = abs(g$estimate) - 0.01)$significance$is_signif)
})


test_that("no effect_cutoff means the p-value alone decides", {
  sig <- estimate_categorical_significance(sa_categorical_fixture(),
                                           by = "table")
  expect_identical(sig$significance$is_signif,
                   sig$significance$pvalue <= 0.05)
  expect_null(attr(sig$significance, "effect_cutoff"))
})


test_that("a table reading follows the test it was pointed at", {
  res <- sa_categorical_fixture()
  for (nm in names(res$tests)) {
    sig <- estimate_categorical_significance(res, by = "table", test = nm)
    expect_equal(sig$significance$pvalue, res$tests[[nm]]$pval, info = nm)
    expect_identical(attr(sig$significance, "test"), nm)
  }
})


test_that("the verdict says which rule produced it", {
  cell <- estimate_categorical_significance(sa_categorical_fixture(),
                                            log2_lift_cutoff = 0.5,
                                            pval_cutoff = 0.1,
                                            adj_type = "holm")$significance

  expect_identical(attr(cell, "by"), "cell")
  expect_identical(attr(cell, "null"), "independence")
  expect_identical(attr(cell, "log2_lift_cutoff"), 0.5)
  expect_identical(attr(cell, "pval_cutoff"), 0.1)
  expect_identical(attr(cell, "adj_type"), "holm")
  expect_identical(attr(cell, "table_dim"), c(2L, 3L))

  table_read <- estimate_categorical_significance(sa_categorical_fixture(),
                                                  by = "table",
                                                  effect_cutoff = 0.2)$significance
  expect_identical(attr(table_read, "by"), "table")
  expect_identical(attr(table_read, "measure"), "cramers_v")
  expect_identical(attr(table_read, "effect_cutoff"), 0.2)
  expect_type(attr(table_read, "test_label"), "character")
})


test_that("a matched pair of conditions has no cell reading", {
  expect_error(
    estimate_categorical_significance(sa_categorical_matched_fixture()),
    "symmetry"
  )
  expect_error(
    estimate_categorical_significance(sa_categorical_matched_fixture()),
    "std_residual"
  )
  # And the refusal names the reading that does work on it.
  expect_error(
    estimate_categorical_significance(sa_categorical_matched_fixture()),
    'by = "table"'
  )
  expect_s3_class(
    estimate_categorical_significance(sa_categorical_matched_fixture(),
                                      by = "table"),
    "sa_categorical_significance"
  )
})


test_that("three or more matched conditions do have a cell reading", {
  res <- sa_categorical_repeated_fixture()
  expect_identical(res$design$null, "marginal_homogeneity")

  tbl <- estimate_categorical_significance(res)$significance
  expect_identical(names(tbl), sa_cell_columns)
  # The condition-by-response table is arithmetically an independent one, so
  # every residual exists rather than being NA the way a symmetry null leaves
  # them.
  expect_true(all(is.finite(tbl$std_residual)))
  expect_true(all(is.finite(tbl$pvalue)))
})


test_that("a numeric comparison is refused and pointed at its own function", {
  expect_error(
    estimate_categorical_significance(sa_two_group_fixture()),
    "estimate_significance"
  )
  expect_error(
    estimate_categorical_significance(sa_factorial_fixture()),
    "estimate_significance"
  )
  expect_error(
    estimate_categorical_significance(data.frame(x = 1)),
    "compare_categorical_groups"
  )
})


test_that("a measure or test this design has no value for is refused", {
  res <- sa_categorical_fixture()
  # A 2 x 3 table has no odds ratio, and the refusal lists what it does have.
  expect_error(
    estimate_categorical_significance(res, by = "table",
                                      measure = "odds_ratio"),
    "cramers_v"
  )
  expect_error(
    estimate_categorical_significance(res, by = "table", test = "mcnemar_test"),
    "chisq_test"
  )
})


test_that("a ratio cutoff below 1 is refused rather than admitting everything", {
  expect_error(
    estimate_categorical_significance(sa_cat_2x2(), by = "table",
                                      effect_cutoff = 0.5),
    "at least 1"
  )
  # The same number is a perfectly ordinary cutoff on a magnitude scale.
  expect_s3_class(
    estimate_categorical_significance(sa_categorical_fixture(), by = "table",
                                      effect_cutoff = 0.5),
    "sa_categorical_significance"
  )
})


test_that("a setting the chosen reading does not read is called out", {
  res <- sa_categorical_fixture()

  expect_warning(estimate_categorical_significance(res, test = "fisher_test"),
                 "not read by")
  expect_warning(estimate_categorical_significance(res, measure = "cramers_v"),
                 "`measure`")
  expect_warning(estimate_categorical_significance(res, effect_cutoff = 0.3),
                 "`effect_cutoff`")
  expect_warning(
    estimate_categorical_significance(res, by = "table", adj_type = "holm"),
    "`adj_type`"
  )
  expect_warning(
    estimate_categorical_significance(res, by = "table",
                                      log2_lift_cutoff = 2),
    "`log2_lift_cutoff`"
  )
  # Passing nothing extra says nothing, and the cutoffs each reading does read
  # are silent on both paths.
  expect_silent(estimate_categorical_significance(res, log2_lift_cutoff = 2,
                                                  adj_type = "holm"))
  expect_silent(estimate_categorical_significance(res, by = "table",
                                                  test = "fisher_test",
                                                  effect_cutoff = 0.3))
})


test_that("an out-of-range cutoff is refused", {
  res <- sa_categorical_fixture()
  expect_error(estimate_categorical_significance(res, log2_lift_cutoff = -1),
               "log2_lift_cutoff")
  expect_error(estimate_categorical_significance(res, pval_cutoff = 2),
               "pval_cutoff")
  expect_error(estimate_categorical_significance(res, adj_type = "nope"),
               "adj_type")
  expect_error(estimate_categorical_significance(res, by = "cells"))
})


test_that("the verdict is not an sa_significance and volcano refuses it", {
  sig <- estimate_categorical_significance(sa_categorical_fixture())

  expect_false(inherits(sig, "sa_significance"))
  expect_s3_class(sig, "sa_result")
  expect_error(draw_volcano_plot(sig))
})


test_that("print reports the rule rather than the table", {
  cell <- estimate_categorical_significance(sa_categorical_fixture())
  out <- capture.output(print(cell))

  expect_true(any(grepl("sa_categorical_significance", out)))
  expect_true(any(grepl("reading", out)))
  expect_true(any(grepl("2 x 3 table", out)))
  expect_true(any(grepl("log2_lift", out)))
  expect_true(any(grepl("cell\\(s\\) significant", out)))

  table_out <- capture.output(print(
    estimate_categorical_significance(sa_categorical_fixture(), by = "table")))
  expect_true(any(grepl("cramers_v", table_out)))
  expect_true(any(grepl("chisq_test", table_out)))

  capture.output(expect_invisible(print(cell)))
})
