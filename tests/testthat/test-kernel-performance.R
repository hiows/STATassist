# Heavy for CRAN check time; full suite still runs under `devtools::test()`
# (NOT_CRAN=true).
skip_on_cran()
# The performance kernels, checked two ways. Where `pROC` computes the same
# quantity it is the oracle, since agreeing with the reference implementation is
# the whole reason for writing these out rather than depending on it. Where it
# does not — IDI and NRI have no counterpart there — the same arithmetic is
# written out a second way beside the test, the way sa_svm_importance() is
# checked.

# Two models of one outcome with a real difference between them, and no ties,
# so the rank arithmetic is exercised on the easy case first.
sa_kernel_fixture <- function(n = 200L, seed = 20260815) {
  set.seed(seed)
  response <- stats::rbinom(n, 1, 0.35)
  list(
    response = response,
    old = stats::plogis(stats::rnorm(n, mean = response * 0.7)),
    new = stats::plogis(stats::rnorm(n, mean = response * 1.4))
  )
}

# The same predictions rounded to one decimal, so most rows share a value and
# every tie rule in the file has to be right.
sa_kernel_tied <- function() {
  d <- sa_kernel_fixture()
  d$old <- round(d$old, 1)
  d$new <- round(d$new, 1)
  d
}

sa_proc_roc <- function(response, predictor) {
  pROC::roc(response, predictor, levels = c(0, 1), direction = "<",
            quiet = TRUE)
}

test_that("the AUC is pROC's, with and without ties", {
  skip_if_not_installed("pROC")
  for (d in list(sa_kernel_fixture(), sa_kernel_tied())) {
    expect_equal(sa_auc(d$response, d$old),
                 as.numeric(pROC::auc(sa_proc_roc(d$response, d$old))))
    expect_equal(sa_auc(d$response, d$new),
                 as.numeric(pROC::auc(sa_proc_roc(d$response, d$new))))
  }
})

test_that("the AUC is the Mann-Whitney statistic it claims to be", {
  d <- sa_kernel_tied()
  event <- d$old[d$response == 1]
  other <- d$old[d$response == 0]
  # Written out as the definition rather than as ranks: the chance a drawn event
  # outranks a drawn non-event, a tie counting as half.
  pairs <- outer(event, other, function(a, b) (a > b) + 0.5 * (a == b))
  expect_equal(sa_auc(d$response, d$old), mean(pairs))
})

test_that("the ROC points are pROC's operating points", {
  skip_if_not_installed("pROC")
  for (d in list(sa_kernel_fixture(), sa_kernel_tied())) {
    mine <- sa_roc_points(d$response, d$old)
    theirs <- sa_proc_roc(d$response, d$old)
    # Compared as sets, since pROC holds them in the other order.
    expect_equal(sort(mine$sensitivity), sort(theirs$sensitivities))
    expect_equal(sort(mine$specificity), sort(theirs$specificities))
    # Ties are one point rather than several: a run of equal predictions cannot
    # be separated by any threshold.
    expect_identical(nrow(mine), length(unique(d$old)) + 1L)
  }
})

test_that("the area under the drawn curve is the AUC", {
  d <- sa_kernel_tied()
  points <- sa_roc_points(d$response, d$old)
  fpr <- 1 - points$specificity
  tpr <- points$sensitivity
  expect_equal(sum(diff(fpr) * (utils::head(tpr, -1) + utils::tail(tpr, -1)) / 2),
               sa_auc(d$response, d$old))
})

test_that("DeLong's standard error is the one pROC's interval is built on", {
  skip_if_not_installed("pROC")
  d <- sa_kernel_fixture()
  mine <- sa_auc_delong(d$response, d$old)
  theirs <- pROC::ci.auc(sa_proc_roc(d$response, d$old), method = "delong")
  expect_equal(mine[["auc"]], as.numeric(theirs[2]))
  expect_equal(mine[["auc"]] - stats::qnorm(0.975) * mine[["se"]],
               as.numeric(theirs[1]))
  expect_equal(mine[["auc"]] + stats::qnorm(0.975) * mine[["se"]],
               as.numeric(theirs[3]))
})

test_that("the placement values average to the AUC on both classes", {
  d <- sa_kernel_tied()
  placement <- sa_placement_values(d$response, d$old)
  auc <- sa_auc(d$response, d$old)
  expect_equal(mean(placement$event), auc)
  expect_equal(mean(placement$other), auc)
  expect_identical(length(placement$event), sum(d$response == 1))
  expect_identical(length(placement$other), sum(d$response == 0))
})

test_that("DeLong's paired test is pROC's roc.test", {
  skip_if_not_installed("pROC")
  for (d in list(sa_kernel_fixture(), sa_kernel_tied())) {
    mine <- sa_delong_test(d$response, d$new, d$old)
    theirs <- pROC::roc.test(sa_proc_roc(d$response, d$new),
                             sa_proc_roc(d$response, d$old),
                             method = "delong", paired = TRUE)
    expect_equal(mine[["statistic"]], as.numeric(theirs$statistic))
    expect_equal(mine[["pval"]], theirs$p.value)
    expect_equal(mine[["delta"]],
                 sa_auc(d$response, d$new) - sa_auc(d$response, d$old))
  }
})

test_that("the paired test is not the unpaired one", {
  d <- sa_kernel_fixture()
  # Two models that ranked the same rows covary, and ignoring that would make
  # every difference between them look more surprising than it is. If the
  # covariance were being dropped the two standard errors would agree.
  paired <- sa_delong_test(d$response, d$new, d$old)
  first <- sa_auc_delong(d$response, d$new)
  second <- sa_auc_delong(d$response, d$old)
  unpaired <- sqrt(first[["se"]]^2 + second[["se"]]^2)
  expect_lt(paired[["se"]], unpaired)
})

test_that("a model held against itself has nothing to report", {
  d <- sa_kernel_fixture()
  same <- sa_delong_test(d$response, d$old, d$old)
  expect_equal(same[["delta"]], 0)
  expect_equal(same[["se"]], 0)
  # Reporting a p-value of 1 would claim a test against a distribution with no
  # spread. There is no number here, and the column says so.
  expect_true(is.na(same[["statistic"]]))
  expect_true(is.na(same[["pval"]]))
  expect_equal(sa_idi(d$response, d$old, d$old)[["idi"]], 0)
  expect_equal(sa_nri(d$response, d$old, d$old)[["nri"]], 0)
})

test_that("the IDI is the two class-wise mean movements, subtracted", {
  d <- sa_kernel_fixture()
  moved <- d$new - d$old
  event <- moved[d$response == 1]
  other <- moved[d$response == 0]
  mine <- sa_idi(d$response, d$old, d$new)

  expect_equal(mine[["idi"]], mean(event) - mean(other))
  # The same thing said the other way round: how much further apart the two
  # classes' predictions ended up than they started.
  expect_equal(mine[["idi"]],
               (mean(d$new[d$response == 1]) - mean(d$new[d$response == 0])) -
                 (mean(d$old[d$response == 1]) - mean(d$old[d$response == 0])))
  # Two independent means, so their variances add.
  expect_equal(mine[["se"]],
               sqrt(stats::var(event) / length(event) +
                      stats::var(other) / length(other)))
  expect_equal(mine[["statistic"]], mine[["idi"]] / mine[["se"]])
  expect_equal(mine[["pval"]],
               2 * stats::pnorm(-abs(mine[["statistic"]])))
})

test_that("the IDI is signed by which way the classes moved", {
  # A new model that pushes every event up and every non-event down improved
  # the separation, whatever it did to the ranking.
  response <- rep(c(0, 1), each = 20)
  old <- rep(0.5, 40)
  new <- ifelse(response == 1, 0.7, 0.3)
  expect_equal(sa_idi(response, old, new)[["idi"]], 0.4)
  expect_equal(sa_idi(response, new, old)[["idi"]], -0.4)
})

test_that("the NRI counts directions and nothing else", {
  d <- sa_kernel_fixture()
  moved <- d$new - d$old
  event <- moved[d$response == 1]
  other <- moved[d$response == 0]
  up_event <- mean(event > 0)
  down_event <- mean(event < 0)
  up_other <- mean(other > 0)
  down_other <- mean(other < 0)
  mine <- sa_nri(d$response, d$old, d$new)

  expect_equal(mine[["nri_event"]], up_event - down_event)
  # A non-event whose probability rose was reclassified away from the truth, so
  # its component reads the other way round.
  expect_equal(mine[["nri_other"]], down_other - up_other)
  expect_equal(mine[["nri"]], mine[["nri_event"]] + mine[["nri_other"]])
  expect_equal(
    mine[["se"]],
    sqrt((up_event + down_event - mine[["nri_event"]]^2) / length(event) +
           (up_other + down_other - mine[["nri_other"]]^2) / length(other))
  )
  expect_equal(mine[["pval"]],
               2 * stats::pnorm(-abs(mine[["nri"]] / mine[["se"]])))
})

test_that("the NRI is blind to how far a probability moved", {
  # Which is what tells it from the IDI: a change of any size counts once.
  response <- rep(c(0, 1), each = 20)
  old <- rep(0.5, 40)
  small <- ifelse(response == 1, 0.51, 0.49)
  large <- ifelse(response == 1, 0.99, 0.01)
  expect_equal(sa_nri(response, old, small)[["nri"]],
               sa_nri(response, old, large)[["nri"]])
  expect_equal(sa_nri(response, old, small)[["nri"]], 2)
  expect_false(isTRUE(all.equal(sa_idi(response, old, small)[["idi"]],
                                sa_idi(response, old, large)[["idi"]])))
})

test_that("the NRI reaches its bounds when every row moves one way", {
  response <- rep(c(0, 1), each = 20)
  old <- rep(0.5, 40)
  # Every event up and every non-event down is the best it can do.
  expect_equal(sa_nri(response, old, ifelse(response == 1, 0.9, 0.1))[["nri"]],
               2)
  expect_equal(sa_nri(response, old, ifelse(response == 1, 0.1, 0.9))[["nri"]],
               -2)
})

test_that("the Brier score is the mean squared distance to the outcome", {
  d <- sa_kernel_fixture()
  expect_equal(sa_brier(d$response, d$old), mean((d$old - d$response)^2))
  # An AUC is blind to what this measures: a model that ranks perfectly and
  # hedges every prediction has an AUC of 1 and a poor Brier score.
  response <- rep(c(0, 1), each = 20)
  hedged <- ifelse(response == 1, 0.6, 0.4)
  expect_equal(sa_auc(response, hedged), 1)
  expect_equal(sa_brier(response, hedged), 0.16)
})

test_that("the threshold scores are counts at the stated cut", {
  d <- sa_kernel_fixture()
  mine <- sa_threshold_scores(d$response, d$old, 0.45)
  called <- d$old >= 0.45
  expect_equal(mine[["accuracy"]], mean(called == (d$response == 1)))
  expect_equal(mine[["sensitivity"]], mean(called[d$response == 1]))
  expect_equal(mine[["specificity"]], mean(!called[d$response == 0]))
  # A cut below everything calls every row an event.
  everything <- sa_threshold_scores(d$response, d$old, 0)
  expect_equal(everything[["sensitivity"]], 1)
  expect_equal(everything[["specificity"]], 0)
})

test_that("a perfect and a reversed ranking are the two ends", {
  response <- rep(c(0, 1), each = 15)
  expect_equal(sa_auc(response, response), 1)
  expect_equal(sa_auc(response, -response), 0)
  # A predictor carrying no information at all ranks nothing, which is the
  # all-ties case rather than a random one.
  expect_equal(sa_auc(response, rep(0.5, 30)), 0.5)
})

test_that("a class that is not there is a refusal rather than a NaN", {
  expect_error(sa_auc(rep(1, 10), stats::runif(10)), "single class")
  expect_error(sa_auc(rep(0, 10), stats::runif(10)), "single class")
  expect_error(sa_auc(c(0, 1, 2), c(0.1, 0.2, 0.3)), "internal error")
  expect_error(sa_auc(c(0, 1), c(0.1, 0.2, 0.3)), "internal error")
})

test_that("one row of a class leaves the spread unavailable, not zero", {
  response <- c(1, rep(0, 12))
  predictor <- seq(0.9, 0.1, length.out = 13)
  area <- sa_auc_delong(response, predictor)
  expect_equal(area[["auc"]], 1)
  expect_true(is.na(area[["se"]]))
  expect_true(is.na(sa_delong_test(response, predictor,
                                   rev(predictor))[["pval"]]))
})
