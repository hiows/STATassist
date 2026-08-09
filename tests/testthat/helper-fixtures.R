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

sa_heatmap_levels <- function() c("ctrl", "case")

# Every scenario the contract has to cover, so a contract test does not have to
# remember to add the newest one.
sa_all_scenarios <- function() {
  list(
    two_group  = sa_two_group_fixture(),
    multi      = sa_multi_group_fixture(),
    repeated   = sa_repeated_fixture(),
    one_sample = sa_one_sample_fixture()
  )
}
