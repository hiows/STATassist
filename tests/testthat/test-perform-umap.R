# The one reduction that standardises nothing by default, so the first thing to pin
# is that the engine really is handed the values as they arrived and that turning the
# two flags on reaches it. Everything else is the same standard the other embedding
# is held to: the `umap` call is reproduced directly, so the derived neighbourhood,
# the metric and the one overridden default are all fixed by the numbers matching.
#
# The whole object is never compared, because `metadata` carries a timestamp.

test_that("the embedding is one row per point and carries no rotation", {
  m <- sa_reduce_matrix()
  res <- perform_umap(m, seed = 1)

  expect_s3_class(res, "sa_reduction")
  expect_s3_class(res, "sa_result")
  expect_named(res, c("analysis", "points", "design", "parameters", "scores",
                      "engine", "fit", "metadata"))
  expect_identical(res$analysis, "umap")
  expect_identical(res$points, rownames(m))

  expect_null(res$variance)
  expect_null(res$loadings)

  expect_named(res$scores, c("points", "UMAP1", "UMAP2"))
  expect_identical(nrow(res$scores), nrow(m))
  expect_identical(res$scores$points, res$points)
  expect_equal(unname(as.matrix(res$scores[-1])), unname(res$fit$layout))
  expect_identical(rownames(res$fit$layout), res$points)

  expect_identical(res$design$point_type, "sample")
  expect_identical(res$parameters$embedding_scale, "samples")
  expect_identical(res$parameters$n_dim, 2L)
  expect_identical(res$parameters$n_neighbors, 15L)
  expect_identical(res$parameters$min_dist, 0.1)
  expect_identical(res$parameters$metric, "euclidean")

  expect_identical(res$engine$package, "umap")
  expect_identical(res$engine$method, "umap")
  expect_identical(res$engine$overridden, character(0))
  expect_identical(res$parameters$method, "naive")

  # The samples are the default scale, so naming it changes nothing at all.
  expect_equal(res$scores,
               perform_umap(m, embedding_scale = "samples", seed = 1)$scores)
})


test_that("the values arrive unscaled unless both flags are turned on", {
  m <- sa_reduce_matrix()

  # This is the difference from `perform_pca()` and `perform_tsne()`: the engine is
  # handed `data` itself, which is the `umap` package's own picture of it.
  res <- perform_umap(m, seed = 11)
  set.seed(11)
  ref <- umap::umap(m, method = "naive", n_neighbors = 15, n_components = 2,
                    min_dist = 0.1, metric = "euclidean")
  expect_equal(unname(as.matrix(res$scores[-1])), unname(ref$layout))
  expect_false(res$parameters$center)
  expect_false(res$parameters$scale)

  # Turned on, they are the same two flags the other two reductions default to, so
  # this is the call that makes the three embeddings comparable.
  std <- perform_umap(m, center = TRUE, scale = TRUE, seed = 11)
  set.seed(11)
  ref <- umap::umap(base::scale(m), method = "naive", n_neighbors = 15,
                    n_components = 2, min_dist = 0.1, metric = "euclidean")
  expect_equal(unname(as.matrix(std$scores[-1])), unname(ref$layout))
  expect_false(isTRUE(all.equal(res$scores$UMAP1, std$scores$UMAP1)))

  # And the arguments that are the method's own reach the engine as passed.
  cos <- perform_umap(m, n_neighbors = 5, min_dist = 0.5, metric = "cosine",
                      n_dim = 3, seed = 11)
  set.seed(11)
  ref <- umap::umap(m, method = "naive", n_neighbors = 5, n_components = 3,
                    min_dist = 0.5, metric = "cosine")
  expect_equal(unname(as.matrix(cos$scores[-1])), unname(ref$layout))
  expect_named(cos$scores, c("points", "UMAP1", "UMAP2", "UMAP3"))

  # Nothing here is limited to three dimensions, which is what the t-SNE error
  # about `n_dim` points at.
  expect_identical(ncol(perform_umap(m, n_dim = 4, seed = 1)$scores), 5L)
})


test_that("a feature of no variance is kept where nothing divides by it", {
  m <- sa_reduce_matrix()
  m[, 2] <- 3

  # The default divides by nothing, so there is nothing a constant feature can
  # break and it stays, contributing the same zero to every distance.
  expect_silent(perform_umap(m, seed = 1))
  kept <- perform_umap(m, seed = 1)
  expect_identical(kept$design$dropped_feats, character(0))
  expect_identical(kept$design$n_feats, ncol(m))

  expect_message(perform_umap(m, center = TRUE, scale = TRUE, seed = 1),
                 "no variance")
  dropped <- suppressMessages(
    perform_umap(m, center = TRUE, scale = TRUE, seed = 1)
  )
  expect_identical(dropped$design$dropped_feats, "f2")
  expect_identical(dropped$design$n_feats, ncol(m) - 1L)
})


test_that("on the feature scale the engine is handed the transpose", {
  m <- sa_reduce_matrix(n = 30L, p = 20L)
  res <- perform_umap(m, embedding_scale = "features", center = TRUE,
                      scale = TRUE, seed = 11)

  expect_identical(res$points, colnames(m))
  expect_identical(res$design$point_type, "feature")
  expect_identical(nrow(res$scores), ncol(m))
  expect_identical(res$scores$points, res$points)

  # The features are standardised and the transpose is then embedded as it stands.
  # Standardising after the transpose would standardise samples.
  set.seed(11)
  ref <- umap::umap(t(base::scale(m)), method = "naive", n_neighbors = 15,
                    n_components = 2, min_dist = 0.1, metric = "euclidean")
  expect_equal(unname(res$fit$layout), unname(ref$layout))

  # Unscaled, the transpose is all the argument does.
  raw <- perform_umap(m, embedding_scale = "features", seed = 11)
  set.seed(11)
  ref <- umap::umap(t(m), method = "naive", n_neighbors = 15, n_components = 2,
                    min_dist = 0.1, metric = "euclidean")
  expect_equal(unname(raw$fit$layout), unname(ref$layout))

  # `design` describes the input rather than the picture, so none of its counts
  # moves with the scale.
  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_used, nrow(m))
  expect_identical(res$design$n_feats, ncol(m))
})


test_that("the neighbourhood size is counted in points, not in samples", {
  m <- sa_reduce_matrix(n = 30L, p = 20L)

  expect_identical(perform_umap(m, seed = 1)$parameters$n_neighbors, 15L)
  expect_identical(
    suppressMessages(
      perform_umap(sa_reduce_matrix(), embedding_scale = "features", seed = 1)
    )$parameters$n_neighbors,
    6L
  )

  # A rejected value quotes the limit in the margin's own words rather than in
  # samples, since the caller is not looking at samples.
  expect_error(perform_umap(m, embedding_scale = "features", n_neighbors = 25),
               "usable feature\\(s\\)")
  expect_error(perform_umap(m, n_neighbors = 40), "usable sample\\(s\\)")

  # Few points is neither an error nor a warning: it runs, and the derived
  # neighbourhood is small enough that it is worth saying out loud.
  expect_message(
    perform_umap(sa_reduce_matrix(), embedding_scale = "features", seed = 1),
    "6 feature\\(s\\) to embed \\(n_neighbors = 6\\)"
  )
  expect_silent(perform_umap(m, seed = 1))
})


test_that("the derived neighbourhood follows the engine's own limits", {
  expect_identical(sa_umap_neighbors(NULL, 24L), 15L)
  expect_identical(sa_umap_neighbors(NULL, 9L), 9L)
  expect_identical(sa_umap_neighbors(4, 9L), 4L)

  # A value that was asked for and cannot be honoured names the limit. One this
  # function derived and cannot honour says instead that the points admit no
  # embedding at all, and is an error for the same reason t-SNE's is.
  expect_error(sa_umap_neighbors(10, 9L), "must not exceed")
  expect_error(sa_umap_neighbors(1, 9L), "n_neighbors")
  expect_error(sa_umap_neighbors(NULL, 1L), "no neighbourhood of 2 or more")
  expect_error(sa_umap_neighbors(NULL, 1L), "perform_pca", fixed = TRUE)
})


test_that("a seed is worth less here, because the engine restores the stream", {
  m <- sa_reduce_matrix()

  a <- perform_umap(m, seed = 42)
  expect_equal(a$scores, perform_umap(m, seed = 42)$scores)
  expect_false(isTRUE(all.equal(a$scores, perform_umap(m, seed = 43)$scores)))

  set.seed(99)
  before <- .Random.seed
  invisible(perform_umap(m, seed = 7))
  expect_identical(.Random.seed, before)

  # `umap` puts the stream back itself, so two seedless runs from the same state
  # agree. Two seedless t-SNE runs do not, which is the difference the two
  # documentations describe.
  set.seed(5)
  expect_equal(perform_umap(m)$scores, perform_umap(m)$scores)
  expect_null(perform_umap(m)$parameters$seed)
})


test_that("what cannot be embedded is refused by name", {
  m <- sa_reduce_matrix()

  expect_error(perform_umap(m, embedding_scale = "genes"), "should be one of")
  expect_error(perform_umap(m, method = "python"), "should be one of")
  expect_error(perform_umap(m, metric = "mahalanobis"), "should be one of")
  expect_error(perform_umap(list(a = 1)), "data.frame or a matrix")
  expect_error(perform_umap(m, feats = "nope"), "not found in")
  expect_error(perform_umap(m[, 1, drop = FALSE]),
               "at least 2 samples and 2 features")
  expect_error(perform_umap(m[, 1, drop = FALSE]), "`perform_umap()`",
               fixed = TRUE)

  expect_error(perform_umap(m, n_dim = 0), "n_dim")
  expect_error(perform_umap(m, min_dist = 0), "min_dist")
  expect_error(perform_umap(m, n_neighbors = 1), "n_neighbors")
  expect_error(perform_umap(m, center = NA), "center")
  expect_error(perform_umap(m, scale = "yes"), "scale")

  # An argument is rejected before the stream is touched, so a refused call and an
  # accepted one do not draw different numbers from the same seed.
  set.seed(3)
  before <- .Random.seed
  expect_error(perform_umap(m, seed = 1, min_dist = 0))
  expect_identical(.Random.seed, before)
})


test_that("print describes the neighbourhood and not a rotation", {
  m <- sa_reduce_matrix()
  out <- capture.output(print(perform_umap(m, seed = 1)))

  expect_match(out[1], "<sa_reduction> umap", fixed = TRUE)
  expect_true(any(grepl("24 sample(s) x 6 feature(s)", out, fixed = TRUE)))
  expect_true(any(grepl("points   : 24 sample(s)", out, fixed = TRUE)))
  expect_true(any(grepl("values as they arrived", out, fixed = TRUE)))
  expect_true(any(grepl("n_neighbors = 15", out, fixed = TRUE)))
  expect_true(any(grepl("min_dist = 0.1, euclidean", out, fixed = TRUE)))
  expect_true(any(grepl("(seed = 1)", out, fixed = TRUE)))
  expect_false(any(grepl("variance", out)))

  std <- capture.output(print(perform_umap(m, center = TRUE, scale = TRUE,
                                           seed = 1)))
  expect_true(any(grepl("centred and scaled", std, fixed = TRUE)))
})


test_that("dropping the engine slot leaves an object that goes out as JSON", {
  res <- perform_umap(sa_reduce_matrix(), seed = 1)

  # This is why `$fit` is declared as the contract's exception: the `umap` object
  # carries a function, so the result as a whole is not portable.
  expect_true(is.function(res$fit$config$metric.function))

  portable <- res[setdiff(names(res), "fit")]
  unportable <- rapply(portable, function(v) is.function(v) || is.environment(v),
                       how = "unlist")
  expect_false(any(unportable))

  skip_if_not_installed("jsonlite")
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$analysis, "umap")
  expect_identical(round_trip$points, res$points)
  expect_identical(round_trip$design$point_type, "sample")
  expect_equal(round_trip$scores$UMAP2, res$scores$UMAP2)
  expect_identical(round_trip$parameters$metric, "euclidean")
})
