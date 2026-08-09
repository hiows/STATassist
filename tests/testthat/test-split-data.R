# A split is judged on what it did not do: no row on both sides, no row on
# neither, and no sampling unit straddling the two. Those are the properties
# pinned here, together with the balance the stratifier is there to preserve.
# The whole object is never compared, because `metadata` carries a timestamp.

# Twenty subjects measured three times each, ten per arm, so every subject is
# entirely inside one arm and the stratifier survives being folded up.
sa_split_repeated <- function(n_subjects = 20L, per_subject = 3L) {
  data.frame(
    subject = rep(paste0("s", seq_len(n_subjects)), each = per_subject),
    arm     = rep(c("control", "treated"),
                  each = n_subjects * per_subject / 2L),
    value   = seq_len(n_subjects * per_subject),
    stringsAsFactors = FALSE
  )
}


test_that("the result has the same shape whatever `times` asks for", {
  one <- split_data(iris, stratified = "Species", seed = 1)
  three <- split_data(iris, stratified = "Species", times = 3, seed = 1)

  expect_s3_class(one, "sa_split")
  expect_named(one, c("full_data", "datasets", "train_idx", "design",
                      "parameters", "metadata"))
  expect_named(three, names(one))

  # One split is a list of one rather than an unwrapped pair, so a consumer
  # never has to ask how many were drawn before it can read the first.
  expect_length(one$datasets, 1L)
  expect_length(three$datasets, 3L)
  expect_named(one$datasets, "Resample1")
  expect_named(three$datasets, c("Resample1", "Resample2", "Resample3"))
  expect_named(one$train_idx, names(one$datasets))

  expect_named(one$datasets[[1]],
               c("train_data", "test_data", "train_rows", "test_rows"))
  expect_identical(one$full_data, iris)
})


test_that("the two halves partition the rows exactly", {
  sp <- split_data(iris, stratified = "Species", times = 3, seed = 2)

  for (d in sp$datasets) {
    expect_length(intersect(d$train_rows, d$test_rows), 0L)
    expect_identical(sort(c(d$train_rows, d$test_rows)), seq_len(nrow(iris)))
    expect_identical(nrow(d$train_data), length(d$train_rows))
    expect_identical(nrow(d$test_data), length(d$test_rows))
  }
})


test_that("the row numbers survive the row names being reset", {
  sp <- split_data(iris, stratified = "Species", seed = 3)
  d <- sp$datasets[[1]]

  expected <- iris[d$train_rows, , drop = FALSE]
  rownames(expected) <- NULL
  expect_identical(d$train_data, expected)
  # The reset is why `train_rows` is carried: it is the only remaining record
  # of where a row came from.
  expect_identical(rownames(d$train_data), as.character(seq_len(nrow(d$train_data))))
  expect_identical(d$train_rows, sp$train_idx[[1]])
})


test_that("each stratum is split in the requested proportion", {
  sp <- split_data(iris, stratified = "Species", p_train = 0.8, seed = 4)
  d <- sp$datasets[[1]]

  # 40 of the 50 in every species, which a simple random draw would only
  # approach on average.
  expect_equal(unname(c(table(d$train_data$Species))), c(40L, 40L, 40L))
  expect_equal(unname(c(table(d$test_data$Species))), c(10L, 10L, 10L))
  expect_identical(sp$design$strata_n,
                   c(setosa = 50L, versicolor = 50L, virginica = 50L))
})


test_that("an unstratified split draws from the rows as a whole", {
  sp <- split_data(iris, p_train = 0.6, seed = 5)

  expect_identical(sp$design$stratified, NA_character_)
  expect_null(sp$design$strata_n)
  expect_identical(nrow(sp$datasets[[1]]$train_data), 90L)
})


test_that("a numeric stratifier is binned rather than matched", {
  sp <- split_data(iris, stratified = "Sepal.Length", seed = 6)

  # The bins live inside createDataPartition() and are not reported back, so
  # there is no count to give even though the split was stratified.
  expect_identical(sp$design$stratified, "Sepal.Length")
  expect_null(sp$design$strata_n)
  expect_gt(nrow(sp$datasets[[1]]$train_data), 100L)

  expect_error(split_data(iris, stratified = rep(1, nrow(iris))),
               "defines no strata")
})


test_that("a column name and the column itself give the same split", {
  by_name <- split_data(iris, stratified = "Species", seed = 7)
  by_vector <- split_data(iris, stratified = iris$Species, seed = 7)

  expect_identical(by_name$train_idx, by_vector$train_idx)
  # Only the record of what was split on differs, since a resolved vector no
  # longer remembers which form it arrived in.
  expect_identical(by_vector$design$stratified, "<vector>")
})


test_that("every row of a sampling unit lands on the same side", {
  d <- sa_split_repeated()
  sp <- split_data(d, stratified = "arm", id = "subject", times = 5, seed = 8)

  expect_identical(sp$design$n_units, 20L)
  expect_identical(sp$design$n_rows, 60L)
  expect_identical(sp$design$id, "subject")

  for (s in sp$datasets) {
    # The whole point of the argument. A row-wise split of this frame would
    # put nearly every subject in both halves.
    expect_length(intersect(s$train_data$subject, s$test_data$subject), 0L)
    # And a unit that is in at all is in whole: three rows, never one or two.
    expect_true(all(table(s$train_data$subject) == 3L))
    expect_true(all(table(s$test_data$subject) == 3L))
  }
})


test_that("with `id` the proportion is over units and the rows are reported", {
  # Unequal unit sizes: one subject carries six rows and the rest carry two, so
  # the row proportion cannot equal the unit proportion.
  d <- data.frame(
    subject = c(rep("s1", 6L), rep(paste0("s", 2:11), each = 2L)),
    value   = seq_len(26L),
    stringsAsFactors = FALSE
  )
  sp <- split_data(d, id = "subject", p_train = 0.75, seed = 9)

  expect_identical(sp$design$n_units, 11L)
  expect_equal(length(sp$datasets[[1]]$train_rows) / 26,
               sp$parameters$achieved_p[["Resample1"]])
  expect_false(isTRUE(all.equal(sp$parameters$achieved_p[["Resample1"]], 0.75)))
})


test_that("a unit whose rows disagree about the stratum is an error", {
  d <- sa_split_repeated()
  d$arm[1] <- "treated"

  # Not a majority vote. A unit assigned as a whole can only carry one stratum,
  # and guessing which produces a split that looks stratified but is not.
  expect_error(split_data(d, stratified = "arm", id = "subject"),
               "must be constant within each `id`")
  expect_error(split_data(d, stratified = "arm", id = "subject"), "s1")
})


test_that("a seed makes the draw reproducible without stealing the stream", {
  expect_identical(split_data(iris, stratified = "Species", seed = 42)$train_idx,
                   split_data(iris, stratified = "Species", seed = 42)$train_idx)
  expect_false(identical(
    split_data(iris, stratified = "Species", seed = 1)$train_idx,
    split_data(iris, stratified = "Species", seed = 2)$train_idx
  ))

  set.seed(99)
  before <- stats::runif(3)
  set.seed(99)
  invisible(split_data(iris, stratified = "Species", seed = 7))
  expect_equal(stats::runif(3), before)
})


test_that("without a seed the draw follows the caller's stream", {
  set.seed(11)
  first <- split_data(iris, stratified = "Species")$train_idx
  set.seed(11)
  expect_identical(split_data(iris, stratified = "Species")$train_idx, first)
})


test_that("repeats are drawn independently", {
  sp <- split_data(iris, stratified = "Species", times = 3, seed = 12)

  # The draft this came from indexed an undefined `idx` inside the loop, so
  # every repeat would have been the same rows if it had run at all.
  expect_false(identical(sp$train_idx[[1]], sp$train_idx[[2]]))
  expect_false(identical(sp$train_idx[[2]], sp$train_idx[[3]]))
  expect_equal(unname(lengths(sp$train_idx)), rep(114L, 3L))
})


test_that("a matrix is accepted and comes back as a data.frame", {
  m <- matrix(seq_len(40L), ncol = 2L,
              dimnames = list(NULL, c("a", "b")))
  sp <- split_data(m, p_train = 0.5, seed = 13)

  expect_s3_class(sp$full_data, "data.frame")
  expect_identical(nrow(sp$datasets[[1]]$train_data), 10L)
})


test_that("the arguments are checked before anything is drawn", {
  expect_error(split_data(list(a = 1)), "data.frame or a matrix")
  expect_error(split_data(iris[0, ]), "zero rows")
  expect_error(split_data(iris, p_train = 1), "`p_train` must be in")
  expect_error(split_data(iris, p_train = 0), "`p_train` must be in")
  expect_error(split_data(iris, times = 0), "`times` must be in")
  expect_error(split_data(iris, times = 2.5), "whole number")
  expect_error(split_data(iris, stratified = rep("a", 3L)),
               "one entry per row")
  expect_error(split_data(iris, id = rep("a", 3L)), "one entry per row")

  na_strata <- as.character(iris$Species)
  na_strata[1] <- NA
  expect_error(split_data(iris, stratified = na_strata), "must not contain NA")

  expect_error(split_data(iris[1, , drop = FALSE]), "at least 2 sampling units")
  # ceiling() takes every unit into training once p_train is high enough, and
  # an empty test set is not a split.
  expect_error(split_data(data.frame(a = 1:4), p_train = 0.99),
               "leaves the test set empty")
})


test_that("print reports what the split was made on", {
  sp <- split_data(sa_split_repeated(), stratified = "arm", id = "subject",
                   times = 2, seed = 14)

  expect_output(print(sp), "<sa_split>")
  expect_output(print(sp), "20 unit\\(s\\) of `subject`")
  expect_output(print(sp), "control 10, treated 10")
  expect_output(print(sp), "Resample2")
  expect_invisible(print(sp))

  expect_output(print(split_data(iris, seed = 15)), "stratify : none")
})
