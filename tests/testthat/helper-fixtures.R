# Fixtures shared by the test files. Built from the datasets that ship with R so
# that a failure can be reproduced without loading anything, and seeded where
# random so that a failure is the same failure twice.

sa_feats <- function() {
  c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
}

sa_two_group_fixture <- function() {
  d <- iris[iris$Species != "setosa", ]
  compare_two_groups(d, sa_feats(), d$Species, c("versicolor", "virginica"),
                     diagnose = FALSE)
}

sa_multi_group_fixture <- function(...) {
  compare_multiple_groups(iris, sa_feats(), iris$Species,
                          levels(iris$Species), diagnose = FALSE, ...)
}

# A balanced within-subject design with a real condition effect, so that the
# omnibus tests have something to find and the post-hoc stage runs.
sa_repeated_frame <- function(n = 14L, effects = c(0, 0.8, 1.4, 0.5)) {
  set.seed(20260801)
  k <- length(effects)
  subject <- rnorm(n, 0, 2)
  values <- as.vector(outer(subject, effects, "+")) + rnorm(n * k, 0, 1)
  data.frame(
    value = values,
    cond  = rep(paste0("t", seq_len(k)), each = n),
    subj  = rep(paste0("s", seq_len(n)), times = k),
    stringsAsFactors = FALSE
  )
}

sa_repeated_fixture <- function(...) {
  d <- sa_repeated_frame()
  compare_multiple_groups(d["value"], "value", d$cond,
                          sort(unique(d$cond)), id = d$subj, paired = TRUE,
                          diagnose = FALSE, ...)
}

sa_one_sample_fixture <- function(...) {
  suppressWarnings(compare_one_sample(iris, sa_feats(), mu = 3,
                                      diagnose = FALSE, ...))
}

# warpbreaks is the balanced two-by-three factorial that ships with R, and being
# balanced is what makes the three types of sums of squares agree, so a test of
# the numbers can be written against aov() without naming a type.
sa_factorial_fixture <- function(...) {
  compare_factorial_groups(warpbreaks, "breaks",
                           list(wool = "wool", tension = "tension"),
                           diagnose = FALSE, ...)
}

# Two blocks of samples and two blocks of features, in the wide layout the
# comparison functions take: one row per sample, one column per feature. Every
# feature follows the sample blocks, so clustering has to recover both splits and
# a test of the drawn order does not rest on which of two near-equal distances
# happened to merge first.
sa_heatmap_matrix <- function() {
  m <- rbind(
    c(10, 1, 20, 5),
    c(11, 2, 21, 6),
    c(12, 1, 22, 5),
    c(1, 10, 2, 15),
    c(2, 11, 3, 16),
    c(1, 12, 4, 15)
  )
  dimnames(m) <- list(paste0("s", 1:6), paste0("f", 1:4))
  m
}

sa_heatmap_group <- function() rep(c("ctrl", "case"), each = 3L)

# The input the three reductions share, in the same wide layout: half the samples
# lifted on half the features, so PC1 has a clean split to find and no two rows are
# alike, which is what t-SNE's duplicate check requires. Kept small so that both
# stochastic engines run in well under a second.
sa_reduce_matrix <- function(n = 24L, p = 6L, seed = 20260809) {
  set.seed(seed)
  m <- matrix(stats::rnorm(n * p), nrow = n,
              dimnames = list(paste0("s", seq_len(n)), paste0("f", seq_len(p))))
  lifted <- seq_len(n / 2)
  m[lifted, seq_len(p / 2)] <- m[lifted, seq_len(p / 2)] + 4
  m
}

# The input the four clusterings share: three tight, well-separated blobs in the
# same wide layout. Separated far enough that all four methods have to agree on it
# — a fixture only one of them can solve would test the fixture rather than the
# contract — and tight enough that the two density methods find three groups and
# no noise at their derived defaults.
sa_cluster_matrix <- function(n_each = 10L, seed = 20260815) {
  set.seed(seed)
  centres <- rbind(c(0, 0, 0, 0), c(9, 9, 0, 0), c(0, 0, 9, 9))
  m <- do.call(rbind, lapply(seq_len(nrow(centres)), function(i) {
    matrix(stats::rnorm(n_each * 4, sd = 0.4), nrow = n_each, byrow = TRUE) +
      rep(centres[i, ], each = n_each)
  }))
  dimnames(m) <- list(paste0("s", seq_len(nrow(m))),
                      paste0("f", seq_len(ncol(m))))
  m
}

# Which blob each row of sa_cluster_matrix() came from. The labels a clustering
# returns are numbered by first appearance, so this is directly comparable.
sa_cluster_truth <- function(n_each = 10L) rep(1:3, each = n_each)

sa_heatmap_levels <- function() c("ctrl", "case")

# A 2 x 3 table with a real association, large enough that every expected count
# clears 5, so the approximation check is quiet and the tests do not have to
# wrap every call in suppressMessages().
sa_categorical_frame <- function() {
  data.frame(
    smoker = rep(c("y", "n"), each = 60),
    grade  = c(rep(c("high", "mid", "low"), c(10, 20, 30)),
               rep(c("high", "mid", "low"), c(30, 20, 10))),
    stringsAsFactors = FALSE
  )
}

sa_categorical_fixture <- function(...) {
  compare_categorical_groups(sa_categorical_frame(), diagnose = FALSE, ...)
}

# A matched 2 x 2 whose discordant cells are lopsided (b = 14, c = 2), so the
# symmetry null is clearly rejected and the residual identity has two distinct
# numbers to hold between rather than a pair that happens to be equal.
sa_categorical_matched_frame <- function() {
  data.frame(
    before = rep(c("pass", "fail"), c(20, 30)),
    after  = c(rep(c("pass", "fail"), c(18, 2)),
               rep(c("pass", "fail"), c(14, 16))),
    stringsAsFactors = FALSE
  )
}

sa_categorical_matched_fixture <- function(...) {
  compare_categorical_groups(sa_categorical_matched_frame(), paired = TRUE,
                             diagnose = FALSE, ...)
}

# Three repeated binary conditions with a rising response rate, which is the
# design Cochran's Q is about and the one whose table is condition by response.
sa_categorical_repeated_frame <- function() {
  data.frame(
    t1 = rep(c("n", "y"), c(30, 10)),
    t2 = rep(c("n", "y"), c(20, 20)),
    t3 = rep(c("n", "y"), c(8, 32)),
    stringsAsFactors = FALSE
  )
}

sa_categorical_repeated_fixture <- function(...) {
  compare_categorical_groups(sa_categorical_repeated_frame(), paired = TRUE,
                             diagnose = FALSE, ...)
}

# Held-out rows to score models on. Every third row rather than a random draw,
# so both halves keep the same class balance and the fixture needs no seed. The
# two species are interleaved in `iris` in blocks of fifty, which is why the
# stride is taken over the whole frame rather than within a block.
sa_perf_split <- function(data) {
  at <- seq(2L, nrow(data), by = 3L)
  list(train = data[-at, , drop = FALSE], test = data[at, , drop = FALSE])
}

# Three models of mpg, nested so that they are genuinely ordered: the baseline
# reads three columns, the second one, and the third a column that barely
# predicts it at all. Fitted without resampling, since what is being tested is
# the held-out scoring rather than the fit.
sa_perf_reg_parts <- function() {
  parts <- sa_perf_split(mtcars)
  train <- parts$train
  fit <- function(preds) {
    fit_linear_regression(train, outcome = "mpg", predictors = preds,
                          cv = FALSE)
  }
  list(
    train = train,
    test  = parts$test,
    models = list(
      full  = fit(c("wt", "hp", "disp")),
      wt    = fit("wt"),
      qsec  = fit("qsec")
    )
  )
}

sa_perf_reg_fixture <- function(...) {
  parts <- sa_perf_reg_parts()
  evaluate_regression_models(parts$models$full,
                             list(wt_only = parts$models$wt,
                                  qsec_only = parts$models$qsec),
                             newdata = parts$test, ...)
}

# The two species `iris` cannot quite separate on sepal measurements alone,
# which is the point: a separable fit warns and returns probabilities of 0 and
# 1, and an AUC of exactly 1 hides every difference the comparisons are for.
sa_perf_cls_parts <- function() {
  d <- iris[iris$Species != "setosa", ]
  d$Species <- factor(d$Species)
  parts <- sa_perf_split(d)
  train <- parts$train
  fit <- function(preds, lv = c("versicolor", "virginica")) {
    fit_logistic_regression(train, outcome = "Species", predictors = preds,
                            outcome_lv = lv, cv = FALSE)
  }
  list(
    train = train,
    test  = parts$test,
    levels = c("versicolor", "virginica"),
    models = list(
      both     = fit(c("Sepal.Length", "Sepal.Width")),
      length   = fit("Sepal.Length"),
      width    = fit("Sepal.Width"),
      # The same predictors pointed at the other class, which is the model that
      # has to be refused rather than quietly reversed.
      flipped  = fit(c("Sepal.Length", "Sepal.Width"),
                     lv = c("virginica", "versicolor"))
    )
  )
}

sa_perf_cls_fixture <- function(...) {
  parts <- sa_perf_cls_parts()
  evaluate_classification_models(parts$models$both,
                                 list(length_only = parts$models$length,
                                      width_only = parts$models$width),
                                 newdata = parts$test, ...)
}

# Every scenario the contract has to cover, so a contract test does not have to
# remember to add the newest one.
sa_all_scenarios <- function() {
  list(
    two_group  = sa_two_group_fixture(),
    multi      = sa_multi_group_fixture(),
    repeated   = sa_repeated_fixture(),
    one_sample = sa_one_sample_fixture(),
    factorial  = sa_factorial_fixture()
  )
}
