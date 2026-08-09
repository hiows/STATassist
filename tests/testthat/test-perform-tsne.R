# What an embedding can be held to is not where a point landed — that is the
# method's business and a different perplexity is a different answer — but that the
# engine was handed the matrix this function says it was handed. Every test that
# matters here reproduces the `Rtsne` call directly: if the scaling drifted, if the
# derived perplexity changed, or if either overridden default came back on, the
# numbers would stop matching. `perform_pca()`'s file is where the shared input
# reading is tested, since it does it without a random number in sight.
#
# The whole object is never compared, because `metadata` carries a timestamp.

test_that("the embedding is one row per point and carries no rotation", {
  m <- sa_reduce_matrix()
  res <- perform_tsne(m, seed = 1)

  expect_s3_class(res, "sa_reduction")
  expect_s3_class(res, "sa_result")
  expect_named(res, c("analysis", "points", "design", "parameters", "scores",
                      "engine", "fit", "metadata"))
  expect_identical(res$analysis, "tsne")
  expect_identical(res$points, rownames(m))

  # There is no second answer to read off this fit, so the two tables a rotation
  # carries are absent rather than present and empty.
  expect_null(res$variance)
  expect_null(res$loadings)

  expect_named(res$scores, c("points", "tSNE1", "tSNE2"))
  expect_identical(nrow(res$scores), nrow(m))
  expect_identical(res$scores$points, res$points)
  expect_equal(unname(as.matrix(res$scores[-1])), res$fit$Y)

  expect_identical(res$design$point_type, "sample")
  expect_identical(res$parameters$embedding_scale, "samples")
  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_feats, ncol(m))
  expect_identical(res$parameters$n_dim, 2L)
  expect_identical(res$parameters$perplexity, 7)
  expect_identical(res$parameters$theta, 0.5)
  expect_identical(res$parameters$seed, 1)

  # The samples are the default scale, so naming it changes nothing at all.
  expect_equal(res$scores,
               perform_tsne(m, embedding_scale = "samples", seed = 1)$scores)
})


test_that("the engine is handed the standardised matrix and its own defaults off", {
  m <- sa_reduce_matrix()
  xs <- base::scale(m)

  # Reproducing the engine call pins three things at once: the scaling, the
  # perplexity that was derived, and the two Rtsne defaults that are turned off.
  # Leaving `normalize` or `pca` on would change these numbers.
  res <- perform_tsne(m, seed = 11)
  set.seed(11)
  ref <- Rtsne::Rtsne(xs, dims = 2, perplexity = 7, theta = 0.5,
                      normalize = FALSE, pca = FALSE, verbose = FALSE)
  expect_equal(unname(as.matrix(res$scores[-1])), ref$Y)
  expect_equal(res$fit$Y, ref$Y)

  # Scaling is on by default so that this and `perform_pca()` see one matrix, and
  # turning it off reaches the engine rather than only the record of it.
  centred <- perform_tsne(m, scale = FALSE, seed = 11)
  set.seed(11)
  ref <- Rtsne::Rtsne(base::scale(m, center = TRUE, scale = FALSE), dims = 2,
                      perplexity = 7, theta = 0.5, normalize = FALSE,
                      pca = FALSE, verbose = FALSE)
  expect_equal(unname(as.matrix(centred$scores[-1])), ref$Y)
  expect_false(centred$parameters$scale)

  # And the two arguments that are the method's own.
  fine <- perform_tsne(m, theta = 0, perplexity = 4, seed = 11)
  set.seed(11)
  ref <- Rtsne::Rtsne(xs, dims = 2, perplexity = 4, theta = 0,
                      normalize = FALSE, pca = FALSE, verbose = FALSE)
  expect_equal(unname(as.matrix(fine$scores[-1])), ref$Y)

  expect_identical(res$engine$package, "Rtsne")
  expect_identical(res$engine$method, "Rtsne")
  expect_identical(res$engine$overridden,
                   c("normalize = FALSE", "pca = FALSE"))
})


test_that("on the feature scale the engine is handed the transpose", {
  m <- sa_reduce_matrix(n = 30L, p = 20L)
  res <- perform_tsne(m, embedding_scale = "features", seed = 11)

  expect_identical(res$points, colnames(m))
  expect_identical(res$design$point_type, "feature")
  expect_identical(nrow(res$scores), ncol(m))
  expect_identical(res$scores$points, res$points)

  # The features are standardised and the transpose is then embedded as it stands.
  # Standardising after the transpose would standardise samples, and these numbers
  # are what would change if it ever did.
  set.seed(11)
  ref <- Rtsne::Rtsne(t(base::scale(m)), dims = 2, perplexity = 6, theta = 0.5,
                      normalize = FALSE, pca = FALSE, verbose = FALSE)
  expect_equal(res$fit$Y, ref$Y)
  expect_false(isTRUE(all.equal(
    res$fit$Y,
    suppressMessages(perform_tsne(t(m), seed = 11))$fit$Y
  )))

  # `design` describes the input rather than the picture, so none of its counts
  # moves with the scale.
  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_used, nrow(m))
  expect_identical(res$design$n_feats, ncol(m))
})


test_that("the neighbourhood size is counted in points, not in samples", {
  m <- sa_reduce_matrix(n = 30L, p = 20L)

  # 30 samples would derive a perplexity of 9. The features are what has
  # neighbours on the feature scale, and there are 20 of them.
  expect_identical(perform_tsne(m, seed = 1)$parameters$perplexity, 9)
  expect_identical(
    perform_tsne(m, embedding_scale = "features", seed = 1)$parameters$perplexity,
    6
  )

  # A rejected value quotes the limit in the margin's own words rather than in
  # samples, since the caller is not looking at samples.
  expect_error(perform_tsne(m, embedding_scale = "features", perplexity = 8),
               "usable feature\\(s\\)")
  expect_error(perform_tsne(m, perplexity = 11), "usable sample\\(s\\)")

  # Few points is neither an error nor a warning: it runs, and the derived
  # neighbourhood is small enough that it is worth saying out loud.
  expect_message(
    perform_tsne(sa_reduce_matrix(), embedding_scale = "features", seed = 1),
    "6 feature\\(s\\) to embed \\(perplexity = 1\\)"
  )
  expect_silent(perform_tsne(m, embedding_scale = "features", seed = 1))
  expect_silent(perform_tsne(m, seed = 1))
})


test_that("the derived perplexity follows the engine's own limit", {
  # Rtsne refuses `n - 1 < 3 * perplexity`, so this is the widest neighbourhood
  # each point count admits, capped at the conventional 30.
  expect_identical(sa_tsne_perplexity(NULL, 24L), 7)
  expect_identical(sa_tsne_perplexity(NULL, 31L), 10)
  expect_identical(sa_tsne_perplexity(NULL, 200L), 30)
  expect_identical(sa_tsne_perplexity(5, 24L), 5)

  # A value that was asked for and cannot be honoured names the limit. One this
  # function derived and cannot honour says instead that the points admit no
  # embedding at all: with one method per call there is no other method left for a
  # skip to protect, so it is an error too.
  expect_error(sa_tsne_perplexity(8, 24L), "must not exceed")
  expect_error(sa_tsne_perplexity(0, 24L), "perplexity")
  expect_error(sa_tsne_perplexity(NULL, 3L), "admit no perplexity")
  expect_error(sa_tsne_perplexity(NULL, 3L), "perform_pca", fixed = TRUE)

  # Four points admit one and three do not, which is why the limit is read rather
  # than guessed at.
  expect_error(
    suppressMessages(perform_tsne(sa_reduce_matrix(n = 3L, p = 4L))),
    "admit no perplexity"
  )
  expect_s3_class(
    suppressMessages(perform_tsne(sa_reduce_matrix(n = 4L, p = 4L), seed = 1)),
    "sa_reduction"
  )
})


test_that("an engine that refuses says so rather than coming back empty", {
  # Two identical rows are what Rtsne refuses, and `iris` is such a data set. The
  # refusal is the engine's own check and is not disabled here, so it arrives as
  # the engine wrote it.
  expect_error(perform_tsne(iris[1:4], seed = 1), "duplicates")
})


test_that("a seed fixes the embedding and puts the stream back", {
  m <- sa_reduce_matrix()

  a <- perform_tsne(m, seed = 42)
  expect_equal(a$scores, perform_tsne(m, seed = 42)$scores)
  expect_false(isTRUE(all.equal(a$scores, perform_tsne(m, seed = 43)$scores)))

  set.seed(99)
  before <- .Random.seed
  invisible(perform_tsne(m, seed = 7))
  expect_identical(.Random.seed, before)

  # Without a seed each run consumes the stream, so two runs are two pictures.
  # That is the method rather than a fault, and it is what `seed` buys.
  set.seed(5)
  expect_false(isTRUE(all.equal(perform_tsne(m)$scores,
                                perform_tsne(m)$scores)))
  expect_null(perform_tsne(m)$parameters$seed)

  # An argument is rejected before the stream is touched, so a refused call and an
  # accepted one do not draw different numbers from the same seed.
  set.seed(3)
  before <- .Random.seed
  expect_error(perform_tsne(m, seed = 1, n_dim = 4))
  expect_identical(.Random.seed, before)
})


test_that("what cannot be embedded is refused by name", {
  m <- sa_reduce_matrix()

  expect_error(perform_tsne(m, embedding_scale = "genes"), "should be one of")
  expect_error(perform_tsne(list(a = 1)), "data.frame or a matrix")
  expect_error(perform_tsne(m, feats = "nope"), "not found in")
  expect_error(perform_tsne(m[, 1, drop = FALSE]),
               "at least 2 samples and 2 features")
  expect_error(perform_tsne(m[, 1, drop = FALSE]), "`perform_tsne()`",
               fixed = TRUE)

  # Rtsne embeds into at most three dimensions, and the message says whose limit
  # it is and which function does not have it.
  expect_error(perform_tsne(m, n_dim = 4), "at most 3")
  expect_error(perform_tsne(m, n_dim = 4), "perform_umap", fixed = TRUE)
  expect_identical(ncol(perform_tsne(m, n_dim = 3, seed = 1)$scores), 4L)

  expect_error(perform_tsne(m, n_dim = 0), "n_dim")
  expect_error(perform_tsne(m, theta = 2), "theta")
  expect_error(perform_tsne(m, center = NA), "center")
  expect_error(perform_tsne(m, scale = "yes"), "scale")
})


test_that("print describes the neighbourhood and not a rotation", {
  m <- sa_reduce_matrix()
  out <- capture.output(print(perform_tsne(m, seed = 1)))

  expect_match(out[1], "<sa_reduction> tsne", fixed = TRUE)
  expect_true(any(grepl("24 sample(s) x 6 feature(s)", out, fixed = TRUE)))
  expect_true(any(grepl("points   : 24 sample(s)", out, fixed = TRUE)))
  expect_true(any(grepl("centred and scaled", out, fixed = TRUE)))
  expect_true(any(grepl("perplexity = 7", out, fixed = TRUE)))
  expect_true(any(grepl("theta = 0.5", out, fixed = TRUE)))
  expect_true(any(grepl("(seed = 1)", out, fixed = TRUE)))

  # No components, so no line about the variance they would have carried.
  expect_false(any(grepl("variance", out)))
  expect_false(any(grepl("component", out)))

  seedless <- capture.output(print(perform_tsne(m)))
  expect_false(any(grepl("seed", seedless)))
})


test_that("dropping the engine slot leaves an object that goes out as JSON", {
  res <- perform_tsne(sa_reduce_matrix(), seed = 1)

  portable <- res[setdiff(names(res), "fit")]
  unportable <- rapply(portable, function(v) is.function(v) || is.environment(v),
                       how = "unlist")
  expect_false(any(unportable))

  skip_if_not_installed("jsonlite")
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$analysis, "tsne")
  expect_identical(round_trip$points, res$points)
  expect_identical(round_trip$design$point_type, "sample")
  expect_equal(round_trip$scores$tSNE1, res$scores$tSNE1)
  expect_equal(round_trip$parameters$perplexity, res$parameters$perplexity)
})
