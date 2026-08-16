# center_by_control() exists to be handed to a comparison, so what it has to
# get right is the set of things the comparison reports afterwards: the fold
# change survives it, every p-value survives it, and the reference centre lands
# on 1. The rest of the file pins the structure it must not disturb, because a
# result whose rows no longer line up with `group` cannot be passed on at all.

sa_norm_frame <- function(scale = c("raw", "log2"),
                          species = levels(iris$Species)) {
  scale <- match.arg(scale)
  d <- iris[iris$Species %in% species, ]
  if (scale == "log2") {
    d[sa_feats()] <- log2(d[sa_feats()])
  }
  d
}

sa_norm_combinations <- function() {
  expand.grid(input_scale = c("raw", "log2"), fc_mean = c("arith", "geom"),
              stringsAsFactors = FALSE)
}


test_that("the two-group effect table is the one found before normalising", {
  combos <- sa_norm_combinations()
  two_lv <- c("versicolor", "virginica")

  for (i in seq_len(nrow(combos))) {
    scale <- combos$input_scale[i]
    mean_type <- combos$fc_mean[i]
    tag <- paste(scale, mean_type)
    d <- sa_norm_frame(scale, c("versicolor", "virginica"))

    norm <- center_by_control(d, sa_feats(), d$Species, two_lv,
                                 fc_mean = mean_type, input_scale = scale)
    before <- compare_two_groups(d, sa_feats(), d$Species, two_lv,
                                 fc_mean = mean_type, input_scale = scale,
                                 diagnose = FALSE)
    after <- compare_two_groups(norm, sa_feats(), d$Species, two_lv,
                                fc_mean = mean_type, input_scale = scale,
                                diagnose = FALSE)

    # Both centres are divided by the same baseline, so the ratio between them
    # cannot move.
    expect_equal(after$effect$fold_change, before$effect$fold_change,
                 info = tag)
    expect_equal(after$effect$log2fc, before$effect$log2fc, info = tag)
    # The reference centre is the one that was divided out, so it is now 1 and
    # `x_center` reads as the fold change on its own.
    expect_equal(after$effect$y_center, rep(1, length(sa_feats())), info = tag)
    expect_equal(after$effect$x_center, after$effect$fold_change, info = tag)
  }
})


test_that("no two-group p-value moves, in any test family", {
  combos <- sa_norm_combinations()
  two_lv <- c("versicolor", "virginica")

  for (i in seq_len(nrow(combos))) {
    scale <- combos$input_scale[i]
    mean_type <- combos$fc_mean[i]
    tag <- paste(scale, mean_type)
    d <- sa_norm_frame(scale, c("versicolor", "virginica"))

    norm <- center_by_control(d, sa_feats(), d$Species, two_lv,
                                 fc_mean = mean_type, input_scale = scale)
    before <- compare_two_groups(d, sa_feats(), d$Species, two_lv,
                                 fc_mean = mean_type, input_scale = scale,
                                 diagnose = FALSE)
    after <- compare_two_groups(norm, sa_feats(), d$Species, two_lv,
                                fc_mean = mean_type, input_scale = scale,
                                diagnose = FALSE)

    # A shift on the log2 scale, a positive rescaling on the raw one: the t
    # statistic is invariant to both and every rank is left where it was.
    for (nm in names(before$tests)) {
      expect_equal(after$tests[[nm]]$pval, before$tests[[nm]]$pval,
                   info = paste(tag, nm))
    }
  }
})


test_that("the multi-group effect table survives it the same way", {
  combos <- sa_norm_combinations()
  group_lv <- levels(iris$Species)

  for (i in seq_len(nrow(combos))) {
    scale <- combos$input_scale[i]
    mean_type <- combos$fc_mean[i]
    tag <- paste(scale, mean_type)
    d <- sa_norm_frame(scale)

    norm <- center_by_control(d, sa_feats(), d$Species, group_lv,
                                 fc_mean = mean_type, input_scale = scale)
    args <- list(feats = sa_feats(), group = d$Species, group_lv = group_lv,
                 fc_mean = mean_type, input_scale = scale, posthoc = FALSE,
                 diagnose = FALSE)
    before <- do.call(compare_multiple_groups, c(list(d), args))
    after <- do.call(compare_multiple_groups, c(list(norm), args))

    expect_equal(after$effect$fold_change, before$effect$fold_change,
                 info = tag)
    expect_equal(after$effect$log2fc, before$effect$log2fc, info = tag)
    expect_equal(after$effect$ref_center, rep(1, length(sa_feats())),
                 info = tag)
    # With the reference at 1 the extreme level's centre is the fold change.
    expect_equal(after$effect$extreme_center, after$effect$fold_change,
                 info = tag)
    for (nm in names(before$tests)) {
      expect_equal(after$tests[[nm]]$pval, before$tests[[nm]]$pval,
                   info = paste(tag, nm))
    }
  }
})


test_that("the quantity removed is the centre the comparison would divide", {
  ctrl <- iris$Species == "setosa"
  group_lv <- levels(iris$Species)
  raw <- iris
  logged <- sa_norm_frame("log2")

  expected <- list(
    list(scale = "raw", mean_type = "arith", data = raw,
         baseline = mean(raw$Petal.Length[ctrl])),
    list(scale = "raw", mean_type = "geom", data = raw,
         baseline = exp(mean(log(raw$Petal.Length[ctrl])))),
    # sa_fc_center() takes the centre on the raw scale, so on log2 input the
    # geometric mean collapses to the arithmetic mean of the logged values.
    list(scale = "log2", mean_type = "geom", data = logged,
         baseline = mean(logged$Petal.Length[ctrl])),
    list(scale = "log2", mean_type = "arith", data = logged,
         baseline = log2(mean(2^logged$Petal.Length[ctrl])))
  )

  for (case in expected) {
    tag <- paste(case$scale, case$mean_type)
    norm <- center_by_control(case$data, "Petal.Length", iris$Species,
                                 group_lv, fc_mean = case$mean_type,
                                 input_scale = case$scale)
    want <- if (case$scale == "log2") {
      case$data$Petal.Length - case$baseline
    } else {
      case$data$Petal.Length / case$baseline
    }
    expect_equal(norm$Petal.Length, want, info = tag)
  }
})


test_that("the control group lands on 1 raw and 0 on the log2 scale", {
  group_lv <- levels(iris$Species)
  ctrl <- iris$Species == "setosa"

  raw <- center_by_control(iris, sa_feats(), iris$Species, group_lv)
  expect_equal(colMeans(raw[ctrl, sa_feats()]),
               rep(1, length(sa_feats())), ignore_attr = TRUE)

  logged <- sa_norm_frame("log2")
  lg <- center_by_control(logged, sa_feats(), logged$Species, group_lv,
                             input_scale = "log2")
  expect_equal(colMeans(lg[ctrl, sa_feats()]),
               rep(0, length(sa_feats())), ignore_attr = TRUE)
})


test_that("fc_mean follows input_scale unless it is named", {
  group_lv <- levels(iris$Species)
  logged <- sa_norm_frame("log2")
  ctrl <- iris$Species == "setosa"
  geom_baseline <- mean(logged$Petal.Length[ctrl])
  arith_baseline <- log2(mean(2^logged$Petal.Length[ctrl]))

  # Different enough that the default cannot be mistaken for the other one.
  expect_false(isTRUE(all.equal(geom_baseline, arith_baseline)))

  by_default <- center_by_control(logged, "Petal.Length", logged$Species,
                                     group_lv, input_scale = "log2")
  expect_equal(by_default$Petal.Length,
               logged$Petal.Length - geom_baseline)

  named <- center_by_control(logged, "Petal.Length", logged$Species,
                                group_lv, fc_mean = "arith",
                                input_scale = "log2")
  expect_equal(named$Petal.Length, logged$Petal.Length - arith_baseline)

  # On the raw scale the default is the arithmetic mean instead.
  on_raw <- center_by_control(iris, "Petal.Length", iris$Species, group_lv)
  expect_equal(on_raw$Petal.Length,
               iris$Petal.Length / mean(iris$Petal.Length[ctrl]))
})


test_that("control_label re-points the baseline without rewriting group_lv", {
  group_lv <- levels(iris$Species)
  virg <- iris$Species == "virginica"

  out <- center_by_control(iris, sa_feats(), iris$Species, group_lv,
                              control_label = "virginica")
  expect_equal(colMeans(out[virg, sa_feats()]),
               rep(1, length(sa_feats())), ignore_attr = TRUE)

  # The error is the comparison's own, since the same helper resolves it.
  expect_error(
    center_by_control(iris, sa_feats(), iris$Species, group_lv,
                         control_label = "nope"),
    "`control_label` names a level `group_lv` does not hold"
  )
})


test_that("everything that is not a named feature comes back untouched", {
  group_lv <- levels(iris$Species)
  out <- center_by_control(iris, sa_feats(), iris$Species, group_lv)

  expect_s3_class(out, "data.frame")
  expect_identical(dim(out), dim(iris))
  expect_identical(names(out), names(iris))
  expect_identical(rownames(out), rownames(iris))
  expect_identical(out$Species, iris$Species)

  # A column left out of `feats` is not normalised at all, which is what lets
  # the result carry an id or a covariate through to the comparison.
  partial <- center_by_control(iris, sa_feats()[1:2], iris$Species, group_lv)
  expect_identical(partial$Petal.Length, iris$Petal.Length)
  expect_identical(partial$Petal.Width, iris$Petal.Width)
  expect_false(isTRUE(all.equal(partial$Sepal.Length, iris$Sepal.Length)))

  # A matrix is accepted and comes back as the data.frame the comparisons take.
  as_matrix <- center_by_control(as.matrix(iris[sa_feats()]), sa_feats(),
                                    iris$Species, group_lv)
  expect_s3_class(as_matrix, "data.frame")
  expect_equal(as_matrix, out[sa_feats()], ignore_attr = TRUE)
})


test_that("rows outside group_lv are kept but take no part in the baseline", {
  sub_lv <- c("setosa", "versicolor")
  ctrl <- iris$Species == "setosa"

  expect_message(
    out <- center_by_control(iris, sa_feats(), iris$Species, sub_lv),
    "Kept 50 row\\(s\\) belonging to a level outside `group_lv`"
  )

  # Dropping them here would leave the result shorter than the `group` vector
  # the comparison still has to be given.
  expect_identical(nrow(out), nrow(iris))
  expect_true(all(is.finite(as.matrix(out[!ctrl, sa_feats()]))))
  # The baseline is the setosa centre alone, so the excluded level neither
  # shifts it nor turns it NA through an NA index.
  expect_equal(out$Petal.Length,
               iris$Petal.Length / mean(iris$Petal.Length[ctrl]))
  expect_equal(colMeans(out[ctrl, sa_feats()]),
               rep(1, length(sa_feats())), ignore_attr = TRUE)

  # The comparison drops them itself, and the arguments do not change. Setosa
  # and versicolor petals do not overlap, so Brunner-Munzel has nothing to
  # estimate and says so; that is the comparison's business, not this one's.
  res <- suppressMessages(suppressWarnings(
    compare_two_groups(out, sa_feats(), iris$Species, sub_lv, diagnose = FALSE)
  ))
  expect_equal(res$effect$y_center, rep(1, length(sa_feats())))
  expect_identical(res$design$n_dropped, 50L)
})


test_that("a baseline that cannot be divided fails its feature and no other", {
  group_lv <- levels(iris$Species)
  ctrl <- iris$Species == "setosa"
  d <- iris[sa_feats()]
  d$zero <- ifelse(ctrl, 0, 5)
  d$negative <- ifelse(ctrl, -3, 5)
  feats <- c(sa_feats(), "zero", "negative")

  expect_warning(
    out <- center_by_control(d, feats, iris$Species, group_lv),
    "could not be taken for 2 of 6 feature\\(s\\)"
  )
  # Dividing by zero would send the column to infinity and dividing by a
  # negative centre would reverse every rank, so neither is carried forward.
  expect_true(all(is.na(out$zero)))
  expect_true(all(is.na(out$negative)))
  expect_equal(colMeans(out[ctrl, sa_feats()]),
               rep(1, length(sa_feats())), ignore_attr = TRUE)

  w <- tryCatch(center_by_control(d, feats, iris$Species, group_lv),
                warning = conditionMessage)
  expect_match(w, "zero: the setosa centre is 0")
  expect_match(w, "negative: the setosa centre is -3")
})


test_that("geom refuses a control group that is not strictly positive", {
  group_lv <- levels(iris$Species)
  d <- iris[sa_feats()]
  d$negative <- ifelse(iris$Species == "setosa", -3, 5)

  # The refusal is sa_fc_center()'s, the same one the comparison would give.
  expect_warning(
    out <- center_by_control(d, "negative", iris$Species, group_lv,
                                fc_mean = "geom"),
    "geometric mean is undefined"
  )
  expect_true(all(is.na(out$negative)))
})


test_that("missing values leave the baseline and the other rows alone", {
  group_lv <- levels(iris$Species)
  d <- iris
  d$Petal.Length[c(1L, 2L)] <- NA
  ctrl <- iris$Species == "setosa"
  kept <- ctrl & !is.na(d$Petal.Length)

  out <- center_by_control(d, sa_feats(), d$Species, group_lv)

  expect_equal(out$Petal.Length,
               d$Petal.Length / mean(d$Petal.Length[kept]))
  # An NA stays exactly where it was rather than spreading to the column.
  expect_identical(which(is.na(out$Petal.Length)), c(1L, 2L))
})


test_that("the argument checks are the ones the comparisons make", {
  group_lv <- levels(iris$Species)

  expect_error(center_by_control(iris, "Species", iris$Species, group_lv),
               "`feats` must refer to numeric columns")
  expect_error(center_by_control(iris, "nope", iris$Species, group_lv),
               "`feats` not found in `data`")
  expect_error(center_by_control(iris, sa_feats(), iris$Species[-1],
                                    group_lv),
               "`group` must have one entry per row of `data`")
  expect_error(center_by_control(iris, sa_feats(), iris$Species, "setosa"),
               "`group_lv` must contain at least 2 levels")
  expect_error(center_by_control(iris, sa_feats(), iris$Species,
                                    c("setosa", "nope")),
               "`group_lv` level\\(s\\) absent from `group`")
})
