# summarize_association_stats() is a screen rather than a contract, so what it
# has to get right is arithmetic and shape: every cell is the number the engine
# would give for that pair, the four matrices are symmetric and share their
# dimnames, and the diagonal says the same thing in all of them. The refusals it
# does not make are as important as the ones it does -- a pair cor.test() cannot
# test must not stop the rest of the screen.

sa_assoc_feats <- function() c("mpg", "disp", "hp", "wt")

sa_assoc_methods <- function() c("pearson", "spearman", "kendall")


test_that("the coefficient is the one stats::cor() gives", {
  feats <- sa_assoc_feats()
  res <- summarize_association_stats(mtcars, feats)

  for (m in sa_assoc_methods()) {
    engine <- stats::cor(as.matrix(mtcars[feats]), method = m)
    expect_equal(res[[m]]$corr, engine, info = m)
  }
})


test_that("the p-value is the one stats::cor.test() gives, for every pair", {
  feats <- sa_assoc_feats()
  res <- summarize_association_stats(mtcars, feats)

  for (m in sa_assoc_methods()) {
    for (j in seq_len(length(feats) - 1L)) {
      for (k in seq.int(j + 1L, length(feats))) {
        engine <- suppressWarnings(
          stats::cor.test(mtcars[[feats[j]]], mtcars[[feats[k]]], method = m)
        )
        expect_equal(res[[m]]$pvalue[j, k], engine$p.value,
                     info = paste(m, feats[j], feats[k]))
      }
    }
  }
})


test_that("the adjustment runs over the pairs rather than over the cells", {
  feats <- sa_assoc_feats()
  n_pairs <- as.integer(choose(length(feats), 2L))

  for (adj in c("BH", "bonferroni", "holm", "none")) {
    res <- summarize_association_stats(mtcars, feats, methods = "pearson",
                                       adj_type = adj)
    up <- upper.tri(res$pearson$pvalue)
    expect_identical(sum(up), n_pairs, info = adj)
    expect_equal(res$pearson$adj_pvalue[up],
                 stats::p.adjust(res$pearson$pvalue[up], method = adj),
                 info = adj)
  }
})


test_that("every matrix is symmetric and carries the features as its dimnames", {
  feats <- sa_assoc_feats()
  res <- summarize_association_stats(mtcars, feats)

  for (m in sa_assoc_methods()) {
    for (slot in c("corr", "pvalue", "adj_pvalue", "n")) {
      x <- res[[m]][[slot]]
      tag <- paste(m, slot)
      expect_identical(dim(x), c(4L, 4L), info = tag)
      expect_identical(dimnames(x), list(feats, feats), info = tag)
      expect_equal(x, t(x), info = tag)
    }
  }
})


test_that("the diagonal is set rather than tested", {
  feats <- sa_assoc_feats()
  res <- summarize_association_stats(mtcars, feats)

  for (m in sa_assoc_methods()) {
    expect_identical(diag(res[[m]]$corr), stats::setNames(rep(1, 4), feats),
                     info = m)
    # A feature is not tested against itself, so there is no p-value to adjust
    # either, and draw_corrplot() can mask by p-value without an exception.
    expect_true(all(is.na(diag(res[[m]]$pvalue))), info = m)
    expect_true(all(is.na(diag(res[[m]]$adj_pvalue))), info = m)
    expect_identical(diag(res[[m]]$n), stats::setNames(rep(32L, 4), feats),
                     info = m)
  }
})


test_that("`n` counts the observations a pair shares", {
  feats <- sa_assoc_feats()
  d <- mtcars[feats]
  d$mpg[1:3] <- NA
  d$wt[c(3, 30)] <- NA

  res <- summarize_association_stats(d, methods = "pearson")
  n <- res$pearson$n

  expect_identical(n["mpg", "disp"], 29L)
  expect_identical(n["mpg", "wt"], 28L)   # row 3 is missing in both
  expect_identical(n["disp", "hp"], 32L)
  expect_identical(diag(n), stats::setNames(c(29L, 32L, 32L, 30L), feats))
})


test_that("an Inf is missing for the same reason an NA is", {
  feats <- sa_assoc_feats()
  d <- mtcars[feats]
  d$mpg[1] <- Inf

  res <- summarize_association_stats(d, methods = "pearson")

  expect_identical(res$pearson$n["mpg", "disp"], 31L)
  expect_equal(res$pearson$corr["mpg", "disp"],
               stats::cor(mtcars$mpg[-1], mtcars$disp[-1]))
})


test_that("a feature with no variance comes back NA rather than aborting", {
  feats <- sa_assoc_feats()
  d <- mtcars[feats]
  d$flat <- 1

  expect_message(
    res <- summarize_association_stats(d, methods = "pearson"),
    "no variance"
  )

  expect_true(all(is.na(res$pearson$corr["flat", feats])))
  expect_true(all(is.na(res$pearson$pvalue["flat", feats])))
  expect_true(all(is.na(res$pearson$adj_pvalue["flat", feats])))
  # The diagonal is structural, so it survives a feature that has no test.
  expect_identical(res$pearson$corr["flat", "flat"], 1)
  # The pairs that could be tested are still there, and the family they were
  # adjusted over is those pairs rather than every cell of the triangle.
  up <- upper.tri(res$pearson$pvalue) & !is.na(res$pearson$pvalue)
  expect_identical(sum(up), 6L)
  expect_equal(res$pearson$adj_pvalue[up],
               stats::p.adjust(res$pearson$pvalue[up], method = "BH"))
})


test_that("`use = \"complete.obs\"` reads every pair on one set of rows", {
  feats <- sa_assoc_feats()
  d <- mtcars[feats]
  d$mpg[1:3] <- NA
  d$wt[30] <- NA
  kept <- stats::complete.cases(d)

  expect_message(
    res <- summarize_association_stats(d, methods = "pearson",
                                       use = "complete.obs"),
    "Dropped 4 row"
  )

  expect_identical(res$design$n_obs, sum(kept))
  expect_true(all(res$pearson$n == sum(kept)))
  expect_equal(res$pearson$corr,
               stats::cor(as.matrix(d[kept, feats])))
  # The p-value follows the same rows, which is the point of the argument.
  expect_equal(res$pearson$pvalue["mpg", "disp"],
               stats::cor.test(d$mpg[kept], d$disp[kept])$p.value)
})


test_that("only the methods that were asked for get a slot", {
  feats <- sa_assoc_feats()

  one <- summarize_association_stats(mtcars, feats, methods = "spearman")
  expect_identical(names(one), c("spearman", "design"))

  # The slots keep the order they were asked in.
  two <- summarize_association_stats(mtcars, feats,
                                     methods = c("kendall", "pearson"))
  expect_identical(names(two), c("kendall", "pearson", "design"))
})


test_that("`feats` picks the columns and their order, and NULL takes the numeric ones", {
  chosen <- summarize_association_stats(mtcars, c("wt", "mpg"),
                                        methods = "pearson")
  expect_identical(dimnames(chosen$pearson$corr), list(c("wt", "mpg"),
                                                       c("wt", "mpg")))

  # iris carries a factor column, which is not a candidate for a correlation.
  all_numeric <- summarize_association_stats(iris, methods = "pearson")
  expect_identical(all_numeric$design$feats,
                   c("Sepal.Length", "Sepal.Width", "Petal.Length",
                     "Petal.Width"))
})


test_that("`design` records the call the matrices were produced under", {
  res <- summarize_association_stats(mtcars, sa_assoc_feats(),
                                     methods = c("pearson", "kendall"),
                                     adj_type = "holm")

  expect_identical(res$design$feats, sa_assoc_feats())
  expect_identical(res$design$n_obs, 32L)
  expect_identical(res$design$methods, c("pearson", "kendall"))
  expect_identical(res$design$adj_type, "holm")
  expect_identical(res$design$use, "pairwise.complete.obs")
})


test_that("a matrix is read as well as a data.frame", {
  feats <- sa_assoc_feats()
  m <- as.matrix(mtcars[feats])

  from_matrix <- summarize_association_stats(m, methods = "pearson")
  from_frame <- summarize_association_stats(mtcars, feats, methods = "pearson")

  expect_equal(from_matrix$pearson$corr, from_frame$pearson$corr)
})


test_that("the arguments that cannot be honoured are refused by name", {
  feats <- sa_assoc_feats()

  expect_error(
    summarize_association_stats(mtcars, feats, methods = "kendal"),
    "Not recognised: kendal"
  )
  expect_error(
    summarize_association_stats(mtcars, feats,
                                methods = c("pearson", "pearson")),
    "duplicated"
  )
  expect_error(
    summarize_association_stats(mtcars, feats, methods = character(0)),
    "non-empty"
  )
  expect_error(
    summarize_association_stats(mtcars, feats, adj_type = "BH2"),
    "adj_type"
  )
  expect_error(
    summarize_association_stats(mtcars, "mpg"),
    "at least 2 features"
  )
  expect_error(
    summarize_association_stats(iris, c("Sepal.Length", "Species")),
    "must refer to numeric columns"
  )
  expect_error(
    summarize_association_stats(iris["Species"]),
    "no numeric column"
  )
  expect_error(
    summarize_association_stats(list(a = 1, b = 2)),
    "data.frame or a matrix"
  )
})
