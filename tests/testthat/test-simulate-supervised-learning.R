# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# A simulator is judged on whether the answer it hands back is the answer that is
# in the data. So the properties pinned here are mostly identities rather than
# tolerances: the outcome is exactly the linear predictor plus the noise, a null
# coefficient is exactly zero, and the term names the truth predicts are the term
# names a fit produces. The few tolerances left are the ones that have to be
# tolerances, where an estimate is being compared with what it estimates.

test_that("both simulators return the same six slots", {
  reg <- simulate_regression(n_samples = 60, n_pred = 4, seed = 1)
  cls <- simulate_classification(n_samples = 60, n_pred = 4, seed = 1)

  expect_named(reg, c("args", "split_args", "truth", "truth_term",
                      "truth_model", "truth_row"))
  expect_named(cls, names(reg))

  # `args` is named after the fitting function it feeds, which is what makes
  # do.call() the whole of the connection between the two.
  expect_named(reg$args, c("data", "outcome", "predictors"))
  expect_named(cls$args, c("data", "outcome", "predictors", "outcome_lv"))
  expect_named(reg$split_args, c("data", "stratified", "id"))

  expect_identical(nrow(reg$truth), length(reg$args$predictors))
  expect_identical(nrow(reg$truth_row), 60L)
  expect_identical(reg$args$outcome, "y")
  # The subject column is never a predictor: with `predictors = NULL` a model
  # would fit on which subject a row came from.
  expect_false("y" %in% reg$args$predictors)
})


test_that("the outcome is exactly what was planted plus the noise", {
  sim <- simulate_regression(n_samples = 80, seed = 2)

  # The identity that makes an unplanted predictor null in the strict sense: the
  # outcome has no source other than the coefficients and the noise.
  expect_equal(sim$args$data$y, sim$truth_row$eta + sim$truth_row$noise)

  null_beta <- sim$truth$beta[sim$truth$role == "null"]
  expect_true(length(null_beta) > 0L)
  expect_true(all(null_beta == 0))
  expect_identical(sim$truth$direction[sim$truth$role == "null"],
                   rep("none", length(null_beta)))
})


test_that("how many coefficients are planted does not depend on the seed", {
  counts <- vapply(1:5, function(s) {
    truth <- simulate_regression(n_samples = 40, n_pred = 10, n_pos = 3,
                                 n_neg = 2, seed = s)$truth
    c(up = sum(truth$direction %in% "up"),
      down = sum(truth$direction %in% "down"))
  }, numeric(2))

  # Which predictors carry them is drawn; how many do is a function of the
  # arguments, so it never has to be looked up.
  expect_true(all(counts["up", ] == 3L))
  expect_true(all(counts["down", ] == 2L))

  # And the sign of a planted coefficient agrees with the direction beside it.
  truth <- simulate_regression(n_samples = 40, n_pred = 10, n_pos = 3,
                               n_neg = 2, seed = 1)$truth
  expect_true(all(truth$beta[truth$direction %in% "up"] > 0))
  expect_true(all(truth$beta[truth$direction %in% "down"] < 0))
})


test_that("`beta` states the coefficients and its length says how many", {
  sim <- simulate_regression(n_samples = 40, beta = c(1.5, 0, -2),
                             n_factor_pred = 0, seed = 3)

  expect_identical(sim$args$predictors, c("x_1", "x_2", "x_3"))
  expect_identical(sim$truth$beta, c(1.5, 0, -2))
  expect_identical(sim$truth$role, c("signal", "null", "signal"))
  expect_identical(sim$truth$direction, c("up", "none", "down"))

  # Two arguments that could disagree are settled together rather than one of
  # them silently winning.
  expect_error(simulate_regression(n_pred = 4, beta = c(1, 2)),
               "`n_pred` asks for 4")
  expect_error(simulate_regression(beta = c(1, 2), n_pos = 1),
               "nothing left for `n_pos` to plant")
  expect_error(simulate_regression(beta = c(1, 2), n_pos = 1, n_neg = 1),
               "`n_pos` and `n_neg`")
})


test_that("the coefficient table is aligned with `truth_term`", {
  sim <- simulate_regression(n_samples = 150, seed = 4)
  fit <- do.call(fit_linear_regression, c(sim$args, cv = FALSE))

  # A factor predictor is several terms and the intercept is a term that is no
  # predictor, so the truth is built on the term axis rather than reindexed from
  # the predictor axis afterwards.
  expect_identical(fit$terms, sim$truth_term$terms)
  expect_identical(fit$coefficients$terms, sim$truth_term$terms)
  expect_identical(sim$truth_term$terms[1], "(Intercept)")
  expect_identical(sim$truth_term$beta[1], sim$truth_model$intercept)

  numeric_at <- match(sim$truth$predictors[sim$truth$role %in%
                                             c("signal", "null")],
                      sim$truth_term$terms)
  expect_equal(sim$truth_term$beta[numeric_at],
               sim$truth$beta[sim$truth$role %in% c("signal", "null")])
})


test_that("a fit recovers what was planted, sign and size", {
  sim <- simulate_regression(n_samples = 400, noise_sd = 1, seed = 5)
  fit <- do.call(fit_linear_regression, c(sim$args, cv = FALSE))

  planted <- sim$truth_term$beta != 0
  planted[1] <- FALSE
  expect_identical(sign(fit$coefficients$estimate[planted]),
                   sign(sim$truth_term$beta[planted]))
  expect_equal(fit$coefficients$estimate[planted],
               sim$truth_term$beta[planted], tolerance = 0.2)
  expect_true(all(fit$coefficients$pval[planted] < 0.01))

  # The model as a whole is scored too, against the share of the variance the
  # predictors were actually given.
  expect_equal(fit$fit_stats$r_squared, sim$truth_model$r_squared,
               tolerance = 0.05)
})


test_that("the intercept is solved for the requested event rate", {
  sim <- simulate_classification(n_samples = 400, event_rate = 0.2, seed = 6)

  expect_equal(sim$truth_row$prob, stats::plogis(sim$truth_row$eta))
  # Solved on the linear predictor that was drawn, so the mean probability is the
  # rate that was asked for exactly and only the Bernoulli draw moves off it.
  expect_equal(mean(sim$truth_row$prob), 0.2)
  expect_lt(abs(sim$truth_model$achieved_event_rate - 0.2), 0.05)
  expect_equal(mean(sim$args$data$y == "case"),
               sim$truth_model$achieved_event_rate)

  # A rate no intercept reaches is said to be one rather than approached.
  expect_error(simulate_classification(n_samples = 40, event_rate = 1),
               "`event_rate` must be in")
})


test_that("`outcome_lv` points the coefficients the way they were planted", {
  sim <- simulate_classification(n_samples = 400, beta = c(2, -2),
                                 n_factor_pred = 0, seed = 7)
  fit <- do.call(fit_logistic_regression, c(sim$args, cv = FALSE))

  expect_identical(sim$args$outcome_lv, c("control", "case"))
  expect_identical(fit$design$outcome_lv, c("control", "case"))
  # A planted positive coefficient raises the chance of the second level, which
  # is the same rule `group_lv` follows in a comparison.
  expect_gt(fit$coefficients$odds_ratio[fit$terms == "x_1"], 1)
  expect_lt(fit$coefficients$odds_ratio[fit$terms == "x_2"], 1)

  # Carried in `args` rather than left out: sorting the labels would put `case`
  # first and report the odds of the class that was not planted.
  reversed <- do.call(fit_logistic_regression,
                      c(sim$args[c("data", "outcome", "predictors")],
                        cv = FALSE))
  expect_identical(reversed$design$outcome_lv, c("case", "control"))
  expect_equal(reversed$coefficients$estimate[reversed$terms == "x_1"],
               -fit$coefficients$estimate[fit$terms == "x_1"])
})


test_that("the answer is not among the predictors", {
  cls <- simulate_classification(n_samples = 60, seed = 8)

  # The draft returned the class probability as a column of the data, which
  # `predictors = NULL` would have read as a predictor and separated the classes
  # perfectly on.
  expect_false("prob" %in% names(cls$args$data))
  expect_false("eta" %in% names(cls$args$data))
  expect_true(all(c("prob", "eta", "draw_prob") %in% names(cls$truth_row)))
  expect_identical(names(cls$args$data), c("y", cls$args$predictors))
})


test_that("a factor predictor becomes one term per level beyond the first", {
  sim <- simulate_regression(n_samples = 150, n_pred = 2, n_factor_pred = 2,
                             factor_lv = c("a", "b", "c"), seed = 9)
  fit <- do.call(fit_linear_regression, c(sim$args, cv = FALSE))

  expect_identical(sim$truth$role[sim$truth$predictors == "x_cat_1"], "factor")
  # One offset per level rather than one coefficient, so the per-predictor answer
  # points at the per-term one instead of holding a number it does not have.
  expect_true(is.na(sim$truth$beta[sim$truth$predictors == "x_cat_1"]))
  expect_identical(
    sim$truth_term$terms,
    c("(Intercept)", "x_1", "x_2", "x_cat_1b", "x_cat_1c", "x_cat_2b",
      "x_cat_2c")
  )
  expect_identical(fit$terms, sim$truth_term$terms)
  expect_identical(sim$truth_term$predictors,
                   c(NA, "x_1", "x_2", "x_cat_1", "x_cat_1", "x_cat_2",
                     "x_cat_2"))

  # Levels are handed out in balanced counts, so how many rows carry each is a
  # function of the arguments and only which rows is drawn.
  expect_equal(unname(c(table(sim$args$data$x_cat_1))), c(50L, 50L, 50L))
  expect_identical(levels(sim$args$data$x_cat_1), c("a", "b", "c"))
})


test_that("a constant predictor is a term the model does not get", {
  sim <- simulate_regression(n_samples = 60, n_pred = 3, n_factor_pred = 0,
                             n_constant_pred = 2, seed = 10)

  expect_identical(sim$truth$predictors[4:5], c("x_const_1", "x_const_2"))
  expect_identical(sim$truth$role[4:5], c("constant", "constant"))
  # It is in `truth` because it was generated, and absent from `truth_term`
  # because it cannot become one.
  expect_false(any(grepl("const", sim$truth_term$terms)))

  expect_message(fit <- do.call(fit_linear_regression, c(sim$args, cv = FALSE)),
                 "single value")
  expect_identical(fit$design$dropped_predictors,
                   c("x_const_1", "x_const_2"))
  expect_identical(fit$terms, sim$truth_term$terms)
})


test_that("`p_missing` holes the numeric predictors and nothing else", {
  sim <- simulate_regression(n_samples = 200, n_pred = 5, n_constant_pred = 1,
                             p_missing = 0.05, seed = 11)
  data <- sim$args$data
  numeric_pred <- paste0("x_", 1:5)

  # The count is a function of the argument and the size of the frame; only which
  # cells is drawn.
  expect_identical(sum(is.na(data[numeric_pred])), 50L)
  expect_false(anyNA(data$y))
  expect_false(anyNA(data$x_cat_1))
  expect_false(anyNA(data$x_const_1))

  # The outcome was computed from the complete values, so a hole is a record that
  # is gone rather than a value that was never there.
  expect_equal(data$y, sim$truth_row$eta + sim$truth_row$noise)

  incomplete <- sum(!stats::complete.cases(data[sim$args$predictors]))
  expect_gt(incomplete, 0L)
  suppressMessages(
    fit <- do.call(fit_linear_regression, c(sim$args, cv = FALSE))
  )
  expect_identical(fit$design$n_dropped, incomplete)
  expect_identical(fit$design$n_used, 200L - incomplete)
})


test_that("`n_per_subject` folds the rows into subjects", {
  sim <- simulate_regression(n_per_subject = c(4, 3, 3, 2), seed = 12)

  expect_identical(sim$truth_model$n_samples, 12L)
  expect_identical(sim$truth_model$n_subject, 4L)
  expect_identical(unname(c(table(sim$args$data$subject))),
                   c(4L, 3L, 3L, 2L))
  expect_identical(sim$truth_row$subject, sim$args$data$subject)
  # The offset is drawn once per subject and reused, which is what makes a
  # subject worth knowing and so what makes a row-wise split leak.
  expect_identical(
    unname(lengths(lapply(split(sim$truth_row$subject_offset,
                               sim$truth_row$subject), unique))),
    rep(1L, 4L)
  )
  expect_identical(sim$truth_model$subject_sd, 1)

  # A single count is spread over the rows, and one that does not divide them is
  # rejected rather than rounded.
  expect_identical(
    simulate_regression(n_samples = 60, n_per_subject = 3,
                        seed = 1)$truth_model$n_subject, 20L
  )
  expect_error(simulate_regression(n_samples = 200, n_per_subject = 3),
               "does not divide")
  expect_error(simulate_regression(n_samples = 40, n_per_subject = c(3, 3)),
               "Drop one of the two")
  expect_error(simulate_regression(n_per_subject = 200), "at least 2")
})


test_that("`subject_share` is the intraclass correlation of the predictors", {
  sim <- simulate_regression(n_pred = 3, n_factor_pred = 0,
                             n_per_subject = rep(4, 500), subject_share = 0.8,
                             value_sd = 2, seed = 16)
  data <- sim$args$data

  icc <- vapply(paste0("x_", 1:3), function(nm) {
    between <- stats::var(tapply(data[[nm]], data$subject, mean))
    between / stats::var(data[[nm]])
  }, numeric(1))
  expect_lt(max(abs(icc - 0.8)), 0.06)

  # The two parts add up to the spread that was asked for, so moving the share
  # moves how alike a subject's rows are and not what the column looks like.
  expect_lt(max(abs(vapply(paste0("x_", 1:3),
                           function(nm) stats::sd(data[[nm]]), numeric(1)) - 2)),
            0.15)

  # At zero a subject's rows are independent draws, so nothing but the outcome
  # offset tells them apart from any other subject's.
  flat <- simulate_regression(n_pred = 3, n_factor_pred = 0,
                              n_per_subject = rep(4, 500), subject_share = 0,
                              seed = 16)$args$data
  expect_lt(stats::var(tapply(flat$x_1, flat$subject, mean)) /
              stats::var(flat$x_1), 0.35)

  expect_error(simulate_regression(subject_share = 1.5),
               "`subject_share` must be in")
})


test_that("a subject is a case or a control as a whole", {
  sim <- simulate_classification(n_per_subject = rep(3, 60), seed = 13)
  data <- sim$args$data

  # Drawn once per subject rather than once per row: a subject that changed class
  # between its samples is not the design `id` exists to protect.
  per_subject <- tapply(data$y, data$subject, function(v) length(unique(v)))
  expect_true(all(per_subject == 1L))
  expect_equal(sim$truth_row$draw_prob,
               unname(rep(vapply(split(sim$truth_row$prob,
                                       factor(data$subject,
                                              levels = unique(data$subject))),
                                 mean, numeric(1)), each = 3L)))
})


test_that("`split_args` names something a split over units can stratify on", {
  flat <- simulate_regression(n_samples = 60, seed = 14)
  reg <- simulate_regression(n_per_subject = rep(3, 40), seed = 14)
  cls <- simulate_classification(n_per_subject = rep(3, 40), seed = 14)

  expect_identical(flat$split_args$stratified, "y")
  expect_null(flat$split_args$id)
  # A continuous outcome varies within a subject, so it cannot be the stratum of
  # a unit that goes to one side of the split as a whole. The categorical
  # predictor is drawn per subject and can.
  expect_identical(reg$split_args$stratified, "x_cat_1")
  expect_identical(reg$split_args$id, "subject")
  expect_identical(cls$split_args$stratified, "y")

  expect_null(simulate_regression(n_per_subject = rep(3, 40),
                                  n_factor_pred = 0,
                                  seed = 14)$split_args$stratified)

  for (sim in list(flat, reg, cls)) {
    sp <- do.call(split_data, c(sim$split_args, seed = 1))
    train <- sp$datasets[[1]]$train_data
    test <- sp$datasets[[1]]$test_data
    if (!is.null(sim$split_args$id)) {
      expect_length(intersect(train$subject, test$subject), 0L)
    }
    expect_identical(nrow(train) + nrow(test), nrow(sim$args$data))
  }
})


test_that("`cor_mat` is the correlation the predictors come out with", {
  cor_mat <- make_block_cor(4, list(list(features = 1:2, cor = 0.8)))
  sim <- simulate_regression(n_samples = 4000, n_pred = 4,
                             beta = c(2, 0, 0, 0), cor_mat = cor_mat,
                             n_factor_pred = 0, seed = 15)

  # Compared cell by cell rather than with a relative tolerance, since the cells
  # that matter most are the zeros.
  observed <- stats::cor(sim$args$data[sim$args$predictors])
  expect_lt(max(abs(unname(observed) - cor_mat)), 0.05)

  # The reason the argument exists: a null predictor correlated with a planted
  # one accounts for a coefficient estimated well away from the zero it has.
  expect_identical(sim$truth$max_cor_signal, c(0, 0.8, 0, 0))
  expect_identical(sim$truth$role, c("signal", "null", "null", "null"))
  fit <- do.call(fit_linear_regression, c(sim$args, cv = FALSE))
  expect_gt(abs(stats::cor(sim$args$data$y, sim$args$data$x_2)), 0.3)
  expect_lt(abs(fit$coefficients$estimate[fit$terms == "x_2"]), 0.2)
})


test_that("a correlation matrix that describes no data is refused", {
  expect_error(simulate_regression(n_pred = 3, cor_mat = diag(2)),
               "must be a numeric 3 x 3 matrix")
  expect_error(
    simulate_regression(n_pred = 2, cor_mat = matrix(c(1, 0.5, 0.2, 1), 2)),
    "must be symmetric"
  )
  expect_error(
    simulate_regression(n_pred = 2, cor_mat = matrix(c(2, 0.5, 0.5, 2), 2)),
    "must have 1 on its diagonal"
  )
  expect_error(
    simulate_regression(n_pred = 2, cor_mat = matrix(c(1, 1.5, 1.5, 1), 2)),
    "outside"
  )
  # Symmetric, unit diagonal, every entry a valid correlation, and still no data
  # has it.
  expect_error(
    simulate_regression(n_pred = 3,
                        cor_mat = matrix(c(1, -0.8, -0.8,
                                           -0.8, 1, -0.8,
                                           -0.8, -0.8, 1), 3)),
    "not positive definite"
  )
})


test_that("make_block_cor assembles what the blocks ask for", {
  cor_mat <- make_block_cor(
    6,
    blocks = list(list(features = 1:2, cor = 0.8),
                  list(features = 3:5, cor = 0.5)),
    default_cor = 0
  )

  expect_identical(dim(cor_mat), c(6L, 6L))
  expect_true(isSymmetric(cor_mat))
  expect_identical(diag(cor_mat), rep(1, 6L))
  expect_identical(cor_mat[1, 2], 0.8)
  expect_identical(cor_mat[3, 4], 0.5)
  expect_identical(cor_mat[4, 5], 0.5)
  expect_identical(cor_mat[1, 3], 0)
  expect_identical(cor_mat[6, 1], 0)

  # Everything outside every block, which is what `default_cor` is for.
  spread <- make_block_cor(3, list(list(features = 1:2, cor = 0.6)),
                           default_cor = 0.2)
  expect_identical(spread[1, 3], 0.2)
  expect_identical(spread[1, 2], 0.6)
  expect_identical(diag(make_block_cor(3, default_cor = 0.4)), rep(1, 3L))
})


test_that("`against` splits a block into two sides that disagree", {
  split <- make_block_cor(
    8, list(list(features = 1:3, cor = 0.9, against = 4:6))
  )

  # Each side agrees with itself at `cor` and disagrees with the other at -`cor`.
  expect_identical(split[1, 2], 0.9)
  expect_identical(split[4, 5], 0.9)
  expect_identical(split[1, 4], -0.9)
  expect_identical(split[3, 6], -0.9)
  expect_identical(split[1, 7], 0)
  expect_identical(diag(split), rep(1, 8L))

  # The point of the argument: a block of six sharing one value could hold no
  # correlation below -0.2, while the split one holds -0.9 at any size, since its
  # smallest eigenvalue is 1 - cor whatever the block covers.
  expect_equal(min(eigen(split, symmetric = TRUE, only.values = TRUE)$values),
               0.1)
  wide <- make_block_cor(
    10, list(list(features = 1:5, cor = 0.95, against = 6:10))
  )
  expect_gt(min(eigen(wide, symmetric = TRUE, only.values = TRUE)$values), 0)

  # A side of one predictor is a predictor that moves against a block.
  lone <- make_block_cor(4, list(list(features = 1, cor = 0.8, against = 2:4)))
  expect_identical(lone[1, 2], -0.8)
  expect_identical(lone[2, 3], 0.8)

  # And the predictors come out with it, which is the whole claim.
  sim <- simulate_regression(n_samples = 4000, n_pred = 6,
                             beta = c(2, 0, 0, 0, 0, 0),
                             cor_mat = split[1:6, 1:6], n_factor_pred = 0,
                             seed = 16)
  drawn <- stats::cor(sim$args$data[paste0("x_", 1:6)])
  expect_equal(drawn[1, 4], -0.9, tolerance = 0.02)
  expect_equal(drawn[4, 5], 0.9, tolerance = 0.02)
})


test_that("make_block_cor refuses what it cannot write down", {
  # A predictor in two blocks would need two correlations with the same partner,
  # and only the later block's would survive.
  expect_error(
    make_block_cor(4, list(list(features = 1:2, cor = 0.5),
                           list(features = 2:3, cor = 0.3))),
    "overlaps an earlier block at predictor\\(s\\) 2"
  )
  expect_error(make_block_cor(3, list(list(features = 2:4, cor = 0.5))),
               "outside the 3")
  expect_error(make_block_cor(3, list(list(features = 1, cor = 0.5))),
               "at least two distinct whole numbers")
  expect_error(make_block_cor(3, list(list(features = 1:2))),
               "must be a list with `features` and `cor`")
  expect_error(make_block_cor(3, list(list(features = 1:2, cor = 1.2))),
               "`blocks\\[\\[1\\]\\]\\$cor` must be in")
  expect_error(make_block_cor(3, default_cor = -2), "`default_cor` must be in")
  expect_error(make_block_cor(0), "`n_features` must be in")

  # One `list()` naming `features` twice is one block and not two, and `$` reads
  # the first of a repeated name, so the rest used to be dropped without a word.
  expect_error(
    make_block_cor(10, list(list(features = 1:3, cor = 0.9,
                                 features = 4:6, cor = -0.4))),
    "names `features` and `cor` more than once"
  )
  expect_error(make_block_cor(3, list(list(features = 1:2, cor = 0.5,
                                           blocks = 3))),
               "holds `blocks`, which a block has no use for")
  expect_error(make_block_cor(3, list(list(1:2, 0.5))),
               "must name every element it holds")

  # A value no block of that size could hold is named as such, with the bound it
  # would have to clear, rather than left to the assembled matrix.
  expect_error(make_block_cor(6, list(list(features = 1:3, cor = -0.6))),
               "holds only above -0.5")
  expect_error(make_block_cor(6, list(list(features = 1:4, cor = -0.4))),
               "holds only above -0.333")
  expect_error(make_block_cor(4, default_cor = -0.5),
               "`default_cor` of -0.5 is not possible among 4 predictors")
  expect_error(make_block_cor(2, list(list(features = 1:2, cor = 1))),
               "at perfect agreement")

  # `against` is what carries the sign, so it takes a positive `cor` and its two
  # sides are two sides.
  expect_error(
    make_block_cor(6, list(list(features = 1:3, cor = -0.9, against = 4:6))),
    "must be above 0 when `against` is given"
  )
  expect_error(
    make_block_cor(6, list(list(features = 1:3, cor = 1, against = 4:6))),
    "puts each side of the block at perfect agreement"
  )
  expect_error(
    make_block_cor(6, list(list(features = 1:3, cor = 0.9, against = 3:5))),
    "predictor\\(s\\) 3 in both `features` and `against`"
  )
  expect_error(
    make_block_cor(6, list(list(features = 1:3, cor = 0.9, against = 4:7))),
    "`blocks\\[\\[1\\]\\]\\$against` indexes predictor\\(s\\) outside the 6"
  )
  # Both sides are claimed, so a later block may not reuse either of them.
  expect_error(
    make_block_cor(8, list(list(features = 1:3, cor = 0.9, against = 4:6),
                           list(features = 6:7, cor = 0.5))),
    "overlaps an earlier block at predictor\\(s\\) 6"
  )

  # What is left for the assembled matrix to catch is how the blocks meet
  # `default_cor`, every block holding on its own by then.
  expect_error(
    make_block_cor(8, list(list(features = 1:3, cor = 0.9, against = 4:6)),
                   default_cor = 0.2),
    "smallest eigenvalue is -0\\.235"
  )
})


test_that("the arguments are checked before anything is drawn", {
  expect_error(simulate_regression(n_pred = 0), "`n_pred` must be in")
  expect_error(simulate_regression(n_pred = 4, n_pos = 3, n_neg = 3),
               "more coefficients than the 4")
  expect_error(simulate_regression(beta_range = c(2, 1)),
               "`beta_range` must be increasing")
  expect_error(simulate_regression(value_sd = c(1, 2)),
               "length 1 or 8")
  expect_error(simulate_regression(value_sd = -1), "must not go below 0")
  expect_error(simulate_regression(noise_sd = -1), "`noise_sd` must be in")
  expect_error(simulate_regression(p_missing = 1), "`p_missing` must be in")
  expect_error(simulate_regression(factor_lv = "only"),
               "at least two distinct")
  expect_error(simulate_regression(pred_prefix = ""),
               "single non-empty string")
  expect_error(simulate_classification(outcome_lv = c("a", "a")),
               "two distinct non-missing class labels")

  # The stream is untouched by a call that rejected its arguments, so a seed
  # cannot mean one data set on one call and another on the next.
  set.seed(21)
  before <- stats::runif(3)
  set.seed(21)
  expect_error(simulate_regression(n_pred = 4, n_pos = 5))
  expect_equal(stats::runif(3), before)
})


test_that("a seed makes the draw reproducible without stealing the stream", {
  expect_equal(simulate_regression(n_samples = 40, seed = 42)$args$data,
               simulate_regression(n_samples = 40, seed = 42)$args$data)
  expect_false(isTRUE(all.equal(
    simulate_regression(n_samples = 40, seed = 1)$args$data,
    simulate_regression(n_samples = 40, seed = 2)$args$data
  )))
  expect_equal(simulate_classification(n_samples = 40, seed = 42)$args$data,
               simulate_classification(n_samples = 40, seed = 42)$args$data)

  set.seed(99)
  before <- stats::runif(3)
  set.seed(99)
  invisible(simulate_regression(n_samples = 40, seed = 7))
  invisible(simulate_classification(n_samples = 40, seed = 7))
  expect_equal(stats::runif(3), before)
})


test_that("without a seed the draw follows the caller's stream", {
  set.seed(11)
  first <- simulate_regression(n_samples = 40)$args$data
  set.seed(11)
  expect_equal(simulate_regression(n_samples = 40)$args$data, first)
})
