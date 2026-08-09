# A rotation is judged on two things a picture cannot be judged on. First, that the
# rows of `scores` are the margin `embedding_scale` asked for and not the other one,
# since transposing the input by hand produces something that looks right and answers
# a third question. Second, that the arithmetic is `prcomp()`'s own, which is pinned
# by reproducing the call directly: if the scaling drifted, the numbers would stop
# matching.
#
# The feature scale is held to a stronger standard than a shape check, because it is
# the one margin with a known right answer: its scores have to be the loadings of the
# sample scale, rescaled, and `variance` has to be untouched.
#
# This is also where the shared input reading is tested, since all three reductions
# read a frame through `sa_reduce_input()` and this is the one that does it without a
# random number in sight.
#
# The whole object is never compared, because `metadata` carries a timestamp.

test_that("the coordinates are one row per point, aligned with `points`", {
  m <- sa_reduce_matrix()
  res <- perform_pca(m)

  expect_s3_class(res, "sa_reduction")
  expect_s3_class(res, "sa_result")
  expect_named(res, c("analysis", "points", "design", "parameters", "variance",
                      "loadings", "scores", "engine", "fit", "metadata"))
  expect_identical(res$analysis, "pca")
  expect_identical(res$points, rownames(m))

  expect_s3_class(res$scores, "data.frame")
  expect_identical(nrow(res$scores), nrow(m))
  expect_identical(res$scores$points, res$points)
  expect_named(res$scores, c("points", paste0("PC", seq_len(ncol(m)))))

  # The engine object is readable on its own terms, which is why it is kept: the
  # sample coordinates are `$fit$x` and they carry the sample names.
  expect_s3_class(res$fit, "prcomp")
  expect_identical(rownames(res$fit$x), res$points)
  expect_identical(rownames(res$fit$rotation), colnames(m))
  expect_equal(unname(as.matrix(res$scores[-1])), unname(res$fit$x))

  # The features are the other axis, and they are the rows of exactly one table.
  expect_identical(res$design$point_type, "sample")
  expect_identical(res$parameters$embedding_scale, "samples")
  expect_identical(nrow(res$loadings), ncol(m))
  expect_identical(res$loadings$variables, colnames(m))
  expect_identical(res$design$feats, colnames(m))
  expect_identical(res$design$n_feats, ncol(m))
  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_used, nrow(m))
  expect_identical(res$design$n_dropped, 0L)

  # The samples are the default scale, so naming it changes nothing at all.
  expect_equal(res$scores, perform_pca(m, embedding_scale = "samples")$scores)

  expect_identical(res$engine$package, "stats")
  expect_identical(res$engine$method, "prcomp")
  expect_identical(res$engine$overridden, character(0))
  expect_true(res$parameters$center)
  expect_true(res$parameters$scale)
})


test_that("the scores are prcomp of the scaled matrix and the variance adds up", {
  m <- sa_reduce_matrix()
  res <- perform_pca(m)

  # prcomp() builds this matrix out of the same call to scale(), so passing the
  # flags and passing the scaled columns have to agree exactly.
  ref <- stats::prcomp(base::scale(m), center = FALSE, scale. = FALSE)
  expect_equal(unname(as.matrix(res$scores[-1])), unname(ref$x))
  expect_equal(unname(as.matrix(res$loadings[-1])), unname(ref$rotation))
  expect_equal(res$fit$x, stats::prcomp(m, center = TRUE, scale. = TRUE)$x)

  v <- ref$sdev^2
  expect_identical(res$variance$component, paste0("PC", seq_len(ncol(m))))
  expect_equal(res$variance$sdev, ref$sdev)
  expect_equal(res$variance$prop_var, v / sum(v) * 100)
  expect_equal(res$variance$cum_var, cumsum(v) / sum(v) * 100)
  expect_equal(res$variance$cum_var[ncol(m)], 100)

  # Every component is reported, not just the two that get drawn, since a share of
  # the variance is only a share if the rest of it is there to be a share of.
  expect_identical(nrow(res$variance), ncol(m))

  # Both flags reach the engine, and off is off.
  raw <- perform_pca(m, center = FALSE, scale = FALSE)
  expect_equal(unname(as.matrix(raw$scores[-1])),
               unname(stats::prcomp(m, center = FALSE, scale. = FALSE)$x))
  expect_false(raw$parameters$center)
})


test_that("the default rotation is of the samples, not of the transpose", {
  m <- sa_reduce_matrix()
  res <- perform_pca(m)
  flipped <- stats::prcomp(t(m), center = TRUE, scale. = TRUE)

  # The transpose is the trap because it comes out the same shape: `rotation` of a
  # feature-by-sample PCA also has one row per sample.
  expect_identical(nrow(flipped$rotation), nrow(m))

  # Same shape, different analysis. scale() standardises columns, so one of the two
  # standardised features and the other standardised samples.
  expect_false(isTRUE(all.equal(res$scores$PC1, unname(flipped$rotation[, 1]))))
  flipped_v <- flipped$sdev^2
  expect_false(isTRUE(all.equal(res$variance$prop_var[1],
                                flipped_v[1] / sum(flipped_v) * 100)))

  # And in the orientation this function uses, PC1 is the planted sample split.
  lifted <- seq_len(nrow(m) / 2)
  expect_true(all(res$scores$PC1[lifted] > 0) ||
                all(res$scores$PC1[lifted] < 0))
  expect_identical(sign(res$scores$PC1[1]),
                   -sign(res$scores$PC1[nrow(m)]))
})


test_that("the feature scale is the loadings of the sample scale, rescaled", {
  m <- sa_reduce_matrix()
  by_samp <- perform_pca(m)
  by_feat <- perform_pca(m, embedding_scale = "features")

  n <- nrow(m)
  p <- ncol(m)
  expect <- by_samp$fit$rotation %*% diag(by_samp$fit$sdev * sqrt(n - 1))

  # The identity that makes this the right definition of a feature scale: the
  # unit-length loadings rescaled to variance-weighted length. It is one fit read
  # from the other end rather than a second fit, so the two agree exactly and the
  # sign of a component is the same sign on both sides.
  expect_equal(unname(as.matrix(by_feat$scores[-1])), unname(expect))
  expect_equal(abs(diag(stats::cor(as.matrix(by_feat$scores[-1]), expect))),
               rep(1, p), tolerance = 1e-10)

  # The matrix never turned, so the fit is the same fit and every number the
  # components carry is the same number. Turning the scale does not move a single
  # axis label.
  expect_equal(by_feat$fit$x, by_samp$fit$x)
  expect_equal(by_feat$variance, by_samp$variance)

  # And the two tables swap: the margin that was not embedded is the other one.
  expect_identical(by_feat$loadings$variables, rownames(m))
  expect_equal(unname(as.matrix(by_feat$loadings[-1])), unname(by_samp$fit$x))
  expect_identical(by_samp$loadings$variables, colnames(m))
})


test_that("`embedding_scale = \"features\"` is not the same as transposing", {
  m <- sa_reduce_matrix()
  by_feat <- perform_pca(m, embedding_scale = "features")
  hand <- perform_pca(t(m))

  # Both come out one row per feature, which is what lets the hand transpose pass
  # for the argument. The numbers are a third analysis: `scale()` standardises
  # columns, so the transpose standardised the 24 samples rather than the 6
  # features. This difference is the whole reason the argument exists.
  expect_identical(hand$points, by_feat$points)
  expect_false(isTRUE(all.equal(hand$scores$PC1, by_feat$scores$PC1)))
  expect_false(isTRUE(all.equal(hand$variance$prop_var,
                                by_feat$variance$prop_var)))

  # What the hand transpose actually did is standardise the samples, which is the
  # third analysis, and this is the call it is equal to.
  expect_equal(unname(as.matrix(hand$scores[-1])),
               unname(stats::prcomp(base::scale(t(m)), center = FALSE,
                                    scale. = FALSE)$x))
})


test_that("the feature scale turns the picture without turning the design", {
  m <- sa_reduce_matrix(n = 30L, p = 20L)
  res <- perform_pca(m, embedding_scale = "features")

  expect_identical(res$points, colnames(m))
  expect_identical(res$design$point_type, "feature")
  expect_identical(res$parameters$embedding_scale, "features")
  expect_identical(nrow(res$scores), ncol(m))
  expect_identical(res$scores$points, res$points)

  # `loadings` is whichever margin was not embedded, so here it is one row per
  # sample, and it is still the components' other side.
  expect_identical(nrow(res$loadings), nrow(m))
  expect_identical(res$loadings$variables, rownames(m))

  # `$fit` is the sample-oriented fit whichever margin was asked for, which is why
  # the points are read off `$scores` and not off `$fit$x`.
  expect_identical(rownames(res$fit$x), rownames(m))
  expect_identical(rownames(res$fit$rotation), colnames(m))

  # `design` describes the input rather than the picture, so none of its counts
  # moves: the samples are still the rows of `data` and `feats` still its columns.
  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_used, nrow(m))
  expect_identical(res$design$n_feats, ncol(m))
  expect_identical(res$design$feats, colnames(m))

  out <- capture.output(print(res))
  expect_true(any(grepl("30 sample(s) x 20 feature(s)", out, fixed = TRUE)))
  expect_true(any(grepl("points   : 20 feature(s)", out, fixed = TRUE)))
})


test_that("rows that are not complete and finite go before the engine runs", {
  m <- sa_reduce_matrix()
  m[2, 1] <- NA
  m[5, 3] <- Inf

  expect_message(perform_pca(m), "not complete and finite")
  res <- suppressMessages(perform_pca(m))

  expect_identical(res$design$n_samples, nrow(m))
  expect_identical(res$design$n_dropped, 2L)
  expect_identical(res$design$n_used, nrow(m) - 2L)
  expect_identical(res$points, rownames(m)[-c(2, 5)])
  expect_identical(nrow(res$scores), nrow(m) - 2L)

  # The deletion is listwise and happens once, so the coordinates are those of the
  # complete rows and of nothing else.
  expect_equal(unname(as.matrix(res$scores[-1])),
               unname(stats::prcomp(m[-c(2, 5), ], center = TRUE,
                                    scale. = TRUE)$x))
})


test_that("a feature of no variance is left out only where it cannot be scaled", {
  m <- sa_reduce_matrix()
  m[, 2] <- 3

  expect_message(perform_pca(m), "no variance")
  res <- suppressMessages(perform_pca(m))
  expect_identical(res$design$dropped_feats, "f2")
  expect_identical(res$design$feats, colnames(m)[-2])
  expect_identical(res$design$n_feats, ncol(m) - 1L)
  expect_identical(nrow(res$loadings), ncol(m) - 1L)
  expect_false(anyNA(res$scores$PC1))

  # Nothing is divided by zero when nothing is divided, so with `scale = FALSE` the
  # feature stays and contributes a component of no variance instead.
  kept <- perform_pca(m, scale = FALSE)
  expect_identical(kept$design$dropped_feats, character(0))
  expect_identical(kept$design$n_feats, ncol(m))
  expect_equal(kept$variance$prop_var[ncol(m)], 0)
})


test_that("a non-numeric column is left out rather than refused", {
  d <- as.data.frame(sa_reduce_matrix())
  d$label <- rep(c("a", "b"), length.out = nrow(d))

  expect_message(perform_pca(d), "non-numeric")
  res <- suppressMessages(perform_pca(d))
  expect_identical(res$design$feats, colnames(sa_reduce_matrix()))

  # `feats` both selects and orders, as it does everywhere else in the package.
  named <- perform_pca(d, feats = c("f3", "f1"))
  expect_identical(named$design$feats, c("f3", "f1"))
  expect_identical(named$loadings$variables, c("f3", "f1"))
  expect_identical(named$variance$component, c("PC1", "PC2"))
})


test_that("sample labels come from the row names, or from the position", {
  m <- sa_reduce_matrix()
  expect_identical(perform_pca(m)$points, rownames(m))

  # A matrix with no dimnames at all is still reducible: the features are named by
  # column and the samples by position.
  bare <- m
  dimnames(bare) <- NULL
  res <- perform_pca(bare)
  expect_identical(res$points, as.character(seq_len(nrow(m))))
  expect_identical(res$design$feats, paste0("V", seq_len(ncol(m))))

  # A repeated sample name is a naming choice rather than an error, so it is kept
  # as it arrived; the tables are aligned by position either way.
  dup <- m
  rownames(dup) <- rep("s1", nrow(m))
  expect_identical(perform_pca(dup)$points, rep("s1", nrow(m)))
})


test_that("what cannot be rotated is refused by name", {
  m <- sa_reduce_matrix()

  expect_error(perform_pca(m, embedding_scale = "genes"), "should be one of")
  expect_error(perform_pca(list(a = 1)), "data.frame or a matrix")
  expect_error(perform_pca(iris["Species"]), "no numeric column")
  expect_error(perform_pca(m, feats = "nope"), "not found in")
  expect_error(perform_pca(m[, 1, drop = FALSE]),
               "at least 2 samples and 2 features")
  expect_error(perform_pca(m[1, , drop = FALSE]),
               "at least 2 samples and 2 features")
  expect_error(perform_pca(m[0, , drop = FALSE]), "zero rows")
  expect_error(perform_pca(m, center = NA), "center")
  expect_error(perform_pca(m, scale = "yes"), "scale")

  # The function that could not run is the one named in the message, since that is
  # the one the caller called.
  expect_error(perform_pca(m[1, , drop = FALSE]), "`perform_pca()`", fixed = TRUE)
})


test_that("print says what was rotated and reports the variance it has", {
  m <- sa_reduce_matrix()
  out <- capture.output(print(perform_pca(m)))

  expect_match(out[1], "<sa_reduction> pca", fixed = TRUE)
  expect_true(any(grepl("24 sample(s) x 6 feature(s)", out, fixed = TRUE)))
  expect_true(any(grepl("points   : 24 sample(s)", out, fixed = TRUE)))
  expect_true(any(grepl("centred and scaled", out, fixed = TRUE)))
  expect_true(any(grepl("3 of 6 component(s)", out, fixed = TRUE)))

  # A rotation has no neighbourhood, so it has no line describing one.
  expect_false(any(grepl("perplexity", out)))
  expect_false(any(grepl("n_neighbors", out)))

  only <- capture.output(print(perform_pca(m, scale = FALSE)))
  expect_true(any(grepl("centred", only, fixed = TRUE)))
  expect_false(any(grepl("centred and scaled", only, fixed = TRUE)))
  bare <- capture.output(print(perform_pca(m, center = FALSE, scale = FALSE)))
  expect_true(any(grepl("values as they arrived", bare, fixed = TRUE)))

  # `n` counts components rather than truncating the object.
  one <- capture.output(print(perform_pca(m), n = 1))
  expect_true(any(grepl("1 of 6 component(s)", one, fixed = TRUE)))
  expect_false(any(grepl("PC2", one, fixed = TRUE)))

  flat <- m
  flat[, 2] <- 3
  dropped <- capture.output(print(suppressMessages(perform_pca(flat))))
  expect_true(any(grepl("f2 (no variance)", dropped, fixed = TRUE)))
})


test_that("dropping the engine slot leaves an object that goes out as JSON", {
  res <- perform_pca(sa_reduce_matrix())

  portable <- res[setdiff(names(res), "fit")]
  unportable <- rapply(portable, function(v) is.function(v) || is.environment(v),
                       how = "unlist")
  expect_false(any(unportable))

  skip_if_not_installed("jsonlite")
  round_trip <- jsonlite::fromJSON(
    jsonlite::toJSON(portable, na = "string", digits = NA)
  )
  expect_identical(round_trip$analysis, "pca")
  expect_identical(round_trip$points, res$points)
  expect_identical(round_trip$design$point_type, "sample")
  expect_identical(round_trip$design$feats, res$design$feats)
  expect_equal(round_trip$scores$PC1, res$scores$PC1)
  expect_equal(round_trip$variance$prop_var, res$variance$prop_var)
  expect_identical(round_trip$parameters$embedding_scale, "samples")
})


test_that("the contract refuses a result that disagrees with itself", {
  # These fire on a mistake inside the package rather than on a call, so they are
  # tested through the assembler directly.
  design <- list(point_type = "sample", n_samples = 3L, n_used = 3L,
                 n_dropped = 0L, n_feats = 2L, feats = c("a", "b"),
                 dropped_feats = character(0))
  scores <- data.frame(points = c("s1", "s2", "s3"), PC1 = 1:3)
  engine <- list(package = "stats", method = "prcomp", label = "l",
                 overridden = character(0))
  variance <- data.frame(component = "PC1", sdev = 1, prop_var = 100,
                         cum_var = 100)
  loadings <- data.frame(variables = c("a", "b"), PC1 = c(1, 0))
  ok <- function(...) {
    args <- list(analysis = "pca", points = c("s1", "s2", "s3"), design = design,
                 parameters = list(), scores = scores, variance = variance,
                 loadings = loadings, engine = engine, fit = NULL)
    replaced <- list(...)
    # Assigned one at a time rather than through modifyList(), which recurses into
    # a data.frame and would rebuild `scores` instead of replacing it.
    for (nm in names(replaced)) args[[nm]] <- replaced[[nm]]
    do.call(sa_new_reduction, args)
  }

  expect_s3_class(ok(), "sa_reduction")
  expect_error(ok(analysis = "pls"), "must be one of")
  expect_error(ok(points = character(0)), "non-empty character")
  expect_error(ok(scores = scores[-1, ]), "aligned with `points`")
  expect_error(ok(engine = engine[-1]), "missing `package`")
  expect_error(ok(variance = NULL), "present exactly when")
  expect_error(ok(loadings = loadings[1, ]), "row\\(s\\) for 2 variable")

  # An embedding is the other side of the same check: it has neither table.
  expect_s3_class(
    sa_new_reduction("tsne", c("s1", "s2", "s3"), design, list(), scores,
                     engine = engine, fit = NULL),
    "sa_reduction"
  )
  expect_error(
    sa_new_reduction("tsne", c("s1", "s2", "s3"), design, list(), scores,
                     variance = variance, engine = engine, fit = NULL),
    "present exactly when"
  )
})
