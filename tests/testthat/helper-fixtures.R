# Fixtures shared by the test files. Built from the datasets that ship with R so
# that a failure can be reproduced without loading anything, and seeded where
# random so that a failure is the same failure twice.

sa_feats <- function() {
  c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
}

sa_two_group_fixture <- function() {
  d <- iris[iris$Species != "setosa", ]
  compare_two_groups(d, sa_feats(), d$Species, c("virginica", "versicolor"),
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
