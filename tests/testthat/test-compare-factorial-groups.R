# compare_factorial_groups() has an external answer for almost everything it
# computes, so most of what follows is an equality against base R rather than a
# tolerance band. stats::aov() reports Type I sums of squares, which coincide
# with Type II and III on a balanced design, so a balanced fixture checks the
# default and an unbalanced one checks `ss_type = "I"`.


# A design whose cells can be given any sizes, so that balance is a property of
# the fixture rather than of the dataset that happened to be available.
sa_fact_frame <- function(sizes = rep(6L, 6L),
                          lv = list(a = c("a1", "a2", "a3"),
                                    b = c("b1", "b2")),
                          seed = 20260812) {
  set.seed(seed)
  cells <- sa_fact_grid(lv)
  blocks <- lapply(seq_len(nrow(cells)), function(k) {
    data.frame(
      a = lv$a[cells$a[k]],
      b = lv$b[cells$b[k]],
      y = stats::rnorm(sizes[k], mean = 0.8 * cells$a[k] + 1.5 * cells$b[k]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, blocks)
}

sa_aov_table <- function(dat) {
  summary(stats::aov(y ~ a * b, data = dat))[[1]]
}


# Cells whose centres are exactly what they are told to be. Two observations
# either side of the cell mean on the log2 scale, where the geometric centre is
# 2^mean and log2() of it is that mean back again, so the term components can be
# checked against arithmetic instead of against a tolerance. `mu` is in
# `sa_fact_grid()` order, the first factor varying fastest.
sa_fact_exact_frame <- function(mu, lv) {
  cells <- sa_fact_grid(lv)
  blocks <- lapply(seq_len(nrow(cells)), function(k) {
    data.frame(
      a = lv$a[cells$a[k]],
      b = lv$b[cells$b[k]],
      y = mu[k] + c(-0.5, 0.5),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, blocks)
}

sa_fact_exact_effect <- function(mu, lv) {
  dat <- sa_fact_exact_frame(mu, lv)
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b),
                                  input_scale = "log2", posthoc = FALSE,
                                  diagnose = FALSE)
  stats::setNames(res$terms$log2_effect, res$terms$terms)
}


test_that("a balanced design matches aov() term for term", {
  dat <- sa_fact_frame()
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b),
                                  posthoc = FALSE, diagnose = FALSE)
  tab <- sa_aov_table(dat)

  expect_identical(res$terms$terms, c("a", "b", "a:b"))
  expect_equal(res$terms$df, tab[1:3, "Df"])
  expect_equal(res$terms$ss, tab[1:3, "Sum Sq"])
  expect_equal(res$terms$ms, tab[1:3, "Mean Sq"])
  expect_equal(res$terms$f_stat, tab[1:3, "F value"])
  expect_equal(res$terms$pval, tab[1:3, "Pr(>F)"])
  expect_equal(unique(res$terms$df_error), tab[4, "Df"])
})

test_that("every sums of squares type agrees on a balanced design", {
  dat <- sa_fact_frame()
  args <- list(data = dat["y"], feats = "y",
               factors = list(a = dat$a, b = dat$b),
               posthoc = FALSE, diagnose = FALSE)
  three <- lapply(c("III", "II", "I"), function(tp) {
    do.call(compare_factorial_groups, c(args, list(ss_type = tp)))$terms$ss
  })
  expect_equal(three[[1]], three[[2]])
  expect_equal(three[[1]], three[[3]])
})

test_that("an unbalanced design matches aov() under ss_type = I", {
  dat <- sa_fact_frame(sizes = c(4L, 9L, 5L, 7L, 6L, 11L))
  tab <- sa_aov_table(dat)
  seq_res <- compare_factorial_groups(dat["y"], "y",
                                      list(a = dat$a, b = dat$b),
                                      ss_type = "I", posthoc = FALSE,
                                      diagnose = FALSE)
  expect_equal(seq_res$terms$ss, tab[1:3, "Sum Sq"])
  expect_equal(seq_res$terms$f_stat, tab[1:3, "F value"])

  # And the types part company there, which is the whole reason the argument
  # exists. The interaction is the one term the three always agree on, being
  # contained in nothing.
  adj_res <- compare_factorial_groups(dat["y"], "y",
                                      list(a = dat$a, b = dat$b),
                                      posthoc = FALSE, diagnose = FALSE)
  expect_false(isTRUE(all.equal(adj_res$terms$ss[1], seq_res$terms$ss[1])))
  expect_equal(adj_res$terms$ss[3], seq_res$terms$ss[3])
})

test_that("the sequential sums of squares add up to the total", {
  dat <- sa_fact_frame(sizes = c(4L, 9L, 5L, 7L, 6L, 11L))
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b),
                                  ss_type = "I", posthoc = FALSE,
                                  diagnose = FALSE)
  ss_total <- sum((dat$y - mean(dat$y))^2)
  ss_within <- sum(vapply(split(dat$y, paste(dat$a, dat$b)), function(v) {
    sum((v - mean(v))^2)
  }, numeric(1)))
  expect_equal(sum(res$terms$ss) + ss_within, ss_total)
  # eta squared is each term's share of the same total, so the shares of a
  # decomposition that adds up have to add up too.
  expect_equal(sum(res$terms$eta_sq), (ss_total - ss_within) / ss_total)
})

test_that("the whole-model test is the one-way ANOVA over the cells", {
  dat <- sa_fact_frame(sizes = c(4L, 9L, 5L, 7L, 6L, 11L))
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b),
                                  posthoc = FALSE, diagnose = FALSE)
  cell <- paste(dat$a, dat$b, sep = ".")
  samples <- split(dat$y, factor(cell, levels = res$design$group_lv))
  one_way <- sa_oneway_anova(samples)

  tbl <- res$tests$anova_test
  expect_equal(tbl$f_stat, unname(one_way[["f_stat"]]))
  expect_equal(tbl$df1, unname(one_way[["df1"]]))
  expect_equal(tbl$df2, unname(one_way[["df2"]]))
  expect_equal(tbl$pval, unname(one_way[["pval"]]))
  expect_equal(tbl$n_cells, 6)
})

test_that("marginal contrasts match TukeyHSD on the same fit", {
  dat <- sa_fact_frame()
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b),
                                  posthoc_alpha = 1, diagnose = FALSE)
  # TukeyHSD() on a factorial fit compares unweighted marginal means against the
  # whole-model error term, which is what the marginal contrasts here are.
  ref <- stats::TukeyHSD(stats::aov(y ~ a * b, data = dat), which = "a")$a
  ours <- subset(res$posthoc$anova_test, factor == "a" & is.na(stratum))

  expect_identical(ours$contrast, sub("-", " - ", rownames(ref)))
  expect_equal(ours$estimate, unname(ref[, "diff"]))
  expect_equal(ours$lower_conf, unname(ref[, "lwr"]))
  expect_equal(ours$upper_conf, unname(ref[, "upr"]))
  expect_equal(ours$pval, unname(ref[, "p adj"]))
  # Already family-wise, so nothing is adjusted a second time.
  expect_identical(ours$pval, ours$pval_adj)
})

test_that("a simple effect holds the other factor at one level", {
  dat <- sa_fact_frame()
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b),
                                  posthoc_alpha = 1, diagnose = FALSE)
  row <- subset(res$posthoc$anova_test,
                factor == "a" & stratum == "b2" & contrast == "a2 - a1")
  at <- function(av, bv) mean(dat$y[dat$a == av & dat$b == bv])
  expect_equal(row$estimate, at("a2", "b2") - at("a1", "b2"))
  expect_equal(row$n1, 6)
  expect_equal(row$n2, 6)
  # The error term is the model's, not the two cells' own, so the degrees of
  # freedom are those of the whole fit.
  expect_equal(row$df, res$tests$anova_test$df2)
})

test_that("a contrast runs only when the term that licenses it does", {
  dat <- sa_fact_frame()
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b),
                                  diagnose = FALSE)
  padj <- stats::setNames(res$terms$pval_adj, res$terms$terms)
  tbl <- res$posthoc$anova_test

  # This fixture has both main effects and no interaction, so the marginal
  # comparisons are there and the simple effects are not.
  expect_lt(padj[["a"]], 0.05)
  expect_gt(padj[["a:b"]], 0.05)
  expect_true(all(is.na(tbl$stratum)))
  expect_setequal(unique(tbl$factor), c("a", "b"))
})

test_that("posthoc_scope picks which contrasts are asked for", {
  dat <- sa_fact_frame()
  args <- list(data = dat["y"], feats = "y",
               factors = list(a = dat$a, b = dat$b),
               posthoc_alpha = 1, diagnose = FALSE)
  both <- do.call(compare_factorial_groups, args)$posthoc$anova_test
  marg <- do.call(compare_factorial_groups,
                  c(args, list(posthoc_scope = "marginal")))$posthoc$anova_test
  simple <- do.call(compare_factorial_groups,
                    c(args, list(posthoc_scope = "simple")))$posthoc$anova_test

  expect_true(all(is.na(marg$stratum)))
  expect_false(any(is.na(simple$stratum)))
  expect_identical(nrow(both), nrow(marg) + nrow(simple))
  expect_identical(nrow(marg), 4L)   # three pairs of a, one of b
})

test_that("the number of factors chooses the ANOVA and its name", {
  lv2 <- list(a = c("a1", "a2"), b = c("b1", "b2"))
  lv3 <- c(lv2, list(c = c("c1", "c2")))
  lv4 <- c(lv3, list(d = c("d1", "d2")))
  expected <- c("two_way", "three_way", "factorial")
  labels <- c("Two-way ANOVA", "Three-way ANOVA", "Factorial ANOVA")

  for (i in seq_along(list(lv2, lv3, lv4))) {
    lv <- list(lv2, lv3, lv4)[[i]]
    grid <- sa_fact_grid(lv)
    dat <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
      row <- lapply(names(lv), function(f) rep(lv[[f]][grid[[f]][k]], 4L))
      names(row) <- names(lv)
      set.seed(k)
      as.data.frame(c(row, list(y = stats::rnorm(4L))),
                    stringsAsFactors = FALSE)
    }))
    res <- compare_factorial_groups(dat["y"], "y", as.list(dat[names(lv)]),
                                    posthoc = FALSE, diagnose = FALSE)
    expect_identical(res$design$anova_type, expected[i])
    expect_identical(res$design$n_factors, length(lv))
    expect_true(startsWith(res$test_info$anova_test$label, labels[i]))
    # Every main effect and every interaction of every order, and nothing else.
    expect_identical(length(unique(res$terms$terms)),
                     as.integer(2^length(lv)) - 1L)
  }
})

test_that("the result keeps the contract on both axes", {
  res <- sa_factorial_fixture()
  expect_s3_class(res, c("sa_factorial", "sa_comparison", "sa_result"),
                  exact = TRUE)
  expect_identical(res$analysis, "factorial_comparison")
  expect_identical(res$tests$anova_test$features, res$features)
  expect_identical(res$effect$features, res$features)
  expect_true(all(sa_test_table_columns() %in% names(res$tests$anova_test)))
  expect_true(all(sa_term_table_columns() %in% names(res$terms)))
  expect_true(all(c(sa_posthoc_table_columns(), "factor", "stratum") %in%
                    names(res$posthoc$anova_test)))
  expect_identical(nrow(res$terms), length(res$features) * 3L)
  expect_identical(res$design$group_lv,
                   c("A.H", "B.H", "A.L", "B.L", "A.M", "B.M"))
  expect_false("pairwise" %in% names(res))
})

test_that("posthoc = FALSE leaves the slot out rather than empty", {
  res <- sa_factorial_fixture(posthoc = FALSE)
  expect_false("posthoc" %in% names(res))
  expect_identical(res$parameters$n_posthoc, 0L)
  expect_true("terms" %in% names(res))
})

test_that("the effect table names the most extreme cell", {
  res <- sa_factorial_fixture()
  cells <- split(warpbreaks$breaks,
                 paste(warpbreaks$wool, warpbreaks$tension, sep = "."))
  means <- vapply(cells[res$design$group_lv], mean, numeric(1))
  ratios <- means / means[1]
  furthest <- which.max(abs(log2(ratios[-1]))) + 1L

  expect_identical(res$effect$extreme_cell, names(means)[furthest])
  expect_equal(res$effect$ref_center, unname(means[1]))
  expect_equal(res$effect$fold_change, unname(ratios[furthest]))
  expect_equal(res$effect$log2fc, log2(unname(ratios[furthest])))
})

test_that("levels come from the data when factor_lv is not given", {
  res <- sa_factorial_fixture()
  expect_identical(res$design$factor_lv,
                   list(wool = c("A", "B"), tension = c("H", "L", "M")))
  # Naming them puts the reference where the caller wants it, and the cell
  # labels and the effect table follow.
  named <- compare_factorial_groups(
    warpbreaks, "breaks", list(wool = "wool", tension = "tension"),
    factor_lv = list(tension = c("L", "M", "H"), wool = c("A", "B")),
    posthoc = FALSE, diagnose = FALSE
  )
  expect_identical(named$design$group_lv[1], "L.A")
  expect_identical(unique(named$terms$terms), c("tension", "wool",
                                                "tension:wool"))
})

test_that("rows outside factor_lv are dropped and counted", {
  expect_message(
    res <- compare_factorial_groups(
      warpbreaks, "breaks", list(wool = "wool", tension = "tension"),
      factor_lv = list(wool = c("A", "B"), tension = c("L", "M")),
      posthoc = FALSE, diagnose = FALSE
    ),
    "Dropped 18 row"
  )
  expect_identical(res$design$n_dropped, 18L)
  expect_identical(res$tests$anova_test$n_used, 36)
  expect_identical(res$design$n_empty_cells, 0L)
})

test_that("a single factor is refused with the sibling function named", {
  expect_error(
    compare_factorial_groups(warpbreaks, "breaks", list(wool = "wool")),
    "compare_multiple_groups"
  )
  expect_error(
    compare_factorial_groups(warpbreaks, "breaks",
                             list(wool = "wool", tension = "tension"),
                             factor_lv = list(wool = c("A", "B"))),
    "factor_lv"
  )
  expect_error(
    compare_factorial_groups(warpbreaks, "breaks",
                             list(wool = "wool", tension = "tension"),
                             factor_lv = list(wool = c("A", "B"),
                                              tension = c("L", "X"))),
    "absent from"
  )
})

test_that("a factor named after a cell table column is refused", {
  # The cell table names a column after each factor, so a factor called `mean`
  # would silently overwrite the means. Better refused at the door than found
  # in the result.
  wb <- warpbreaks
  wb$mean <- wb$wool
  expect_error(
    compare_factorial_groups(wb, "breaks",
                             list(mean = "mean", tension = "tension")),
    "Reserved: features, cell, n, mean, sd, se"
  )
})

test_that("a within-subject factor points at what is missing", {
  expect_error(
    compare_factorial_groups(warpbreaks, "breaks",
                             list(wool = "wool", tension = "tension"),
                             within = "tension"),
    "not implemented yet"
  )
  expect_warning(
    compare_factorial_groups(warpbreaks, "breaks",
                             list(wool = "wool", tension = "tension"),
                             id = seq_len(nrow(warpbreaks)),
                             posthoc = FALSE, diagnose = FALSE),
    "ignored"
  )
})

test_that("an empty cell fails that feature rather than the run", {
  dat <- sa_fact_frame()
  dat <- dat[!(dat$a == "a3" & dat$b == "b2"), , drop = FALSE]
  dat$other <- dat$y
  # The empty cell is announced once, and then costs both feature-axis tables
  # their rows: the fold change has no reference to reach, the model has a cell
  # it cannot estimate.
  expect_warning(
    expect_warning(
      expect_message(
        res <- compare_factorial_groups(dat[c("y", "other")], c("y", "other"),
                                        list(a = dat$a, b = dat$b),
                                        factor_lv = list(a = c("a1", "a2", "a3"),
                                                         b = c("b1", "b2")),
                                        posthoc = FALSE, diagnose = FALSE),
        "hold no observation"
      ),
      "fold change could not be computed"
    ),
    "Two-way ANOVA .* could not be computed"
  )
  expect_identical(res$design$n_empty_cells, 1L)
  expect_true(all(is.na(res$tests$anova_test$pval)))
  # Named rather than read off `res$terms$pval` alone: a term table that lost its
  # statistics columns altogether answers `all(is.na(NULL))` with TRUE, which is
  # how a fit list shortened by a failed feature used to pass this unnoticed.
  expect_true(all(sa_term_table_columns() %in% names(res$terms)))
  expect_true(all(is.na(res$terms$pval)))
  expect_true(all(is.na(res$terms$f_stat)))
  # The rows are still there, since the questions were asked whether or not they
  # could be answered.
  expect_identical(nrow(res$terms), 6L)
  # And so are the cells, whose means came from the samples rather than the fit:
  # a cell that holds observations reports its mean even where the model could
  # not be estimated, and only the pooled `se` needs the fit.
  expect_identical(nrow(res$cells), 12L)
  expect_identical(sum(is.na(res$cells$mean)), 2L)
  expect_true(all(is.na(res$cells$se)))
})

test_that("a feature with missing values is fitted on what is left", {
  dat <- sa_fact_frame()
  dat$y[c(1L, 2L, 20L)] <- NA_real_
  # Dropping three observations leaves the cells unequal, so this is also the
  # case where the sequential type is the one aov() can be compared against.
  res <- compare_factorial_groups(dat["y"], "y",
                                  list(a = dat$a, b = dat$b),
                                  ss_type = "I", posthoc = FALSE,
                                  diagnose = FALSE)
  expect_identical(res$tests$anova_test$n_used, nrow(dat) - 3)
  expect_equal(res$terms$f_stat,
               sa_aov_table(dat[!is.na(dat$y), ])[1:3, "F value"])
})

test_that("the simulator's args run without translation", {
  sim <- simulate_factorial_groups(n_feats = 6, n_per_cell = 6, n_up = 2,
                                   n_down = 2, seed = 7)
  res <- suppressWarnings(suppressMessages(
    do.call(compare_factorial_groups, sim$args)
  ))
  expect_identical(res$features, sim$args$feats)
  expect_identical(res$parameters$input_scale, "log2")
  expect_identical(res$design$factor_lv, sim$args$factor_lv)

  # The two tables are keyed the same way, which is the point of the term axis.
  merged <- merge(res$terms, sim$truth_term, by = c("features", "terms"))
  expect_identical(nrow(merged), nrow(res$terms))
  expect_identical(sort(unique(merged$terms)), sort(unique(sim$truth_term$terms)))
})

test_that("the post-hoc table merges with truth_contrast", {
  sim <- simulate_factorial_groups(n_feats = 8, n_per_cell = 10, seed = 11)
  res <- suppressWarnings(suppressMessages(
    do.call(compare_factorial_groups, c(sim$args, list(posthoc_alpha = 1)))
  ))
  ph <- res$posthoc$anova_test
  merged <- merge(ph, sim$truth_contrast,
                  by = c("features", "factor", "stratum", "contrast"))
  # `stratum` is NA for a marginal contrast, and merge() matches NA to NA, so a
  # marginal row is not silently dropped.
  expect_identical(nrow(merged), nrow(ph))
  expect_true(all(merged$group1.x == merged$group1.y))
  expect_gt(cor(merged$estimate, merged$delta), 0.5)
})

test_that("log2_effect is the term's own share of the cell centres", {
  three <- list(a = c("a1", "a2", "a3"), b = c("b1", "b2"))
  # The same profile in both levels of `b`, so `a` carries everything. Its
  # largest component is the level furthest from the grand mean of 5/3.
  eff <- sa_fact_exact_effect(rep(c(0, 1, 4), 2), three)
  expect_equal(eff[["a"]], 4 - 5 / 3)
  expect_equal(eff[["b"]], 0)
  expect_equal(eff[["a:b"]], 0)

  # A crossover: `a` raises one level of `b` by as much as it lowers the other,
  # so both marginal profiles are flat and only the interaction is left.
  two <- list(a = c("a1", "a2"), b = c("b1", "b2"))
  eff <- sa_fact_exact_effect(c(1, -1, -1, 1), two)
  expect_equal(eff[["a"]], 0)
  expect_equal(eff[["b"]], 0)
  expect_equal(eff[["a:b"]], 1)

  # An additive design: the marginal difference is 2 in each factor, and a
  # component is half of it, being a deviation from the grand mean rather than a
  # difference between two levels.
  eff <- sa_fact_exact_effect(c(0, 2, 2, 4), two)
  expect_equal(abs(eff[["a"]]), 1)
  expect_equal(abs(eff[["b"]]), 1)
  expect_equal(eff[["a:b"]], 0)
})

test_that("the components of every term add back up to the cells", {
  lv <- list(a = c("a1", "a2", "a3"), b = c("b1", "b2"))
  cells <- sa_fact_grid(lv)
  eff <- c(0.3, 1.1, -0.4, 2, 0.7, -1.3)

  parts <- vapply(list("a", "b", c("a", "b")),
                  function(term) sa_fact_component(eff, cells, term),
                  numeric(nrow(cells)))
  # The decomposition is exhaustive: main effects, interaction and grand mean are
  # the cells again. Nothing is left over for a term that does not exist.
  expect_equal(rowSums(parts) + mean(eff), eff)
})

test_that("a crossover feature's effect lands on the interaction term", {
  # The spreads are tightened from the defaults so that the noise floor of a
  # component sits well below the smallest planted one. A component of a cell
  # mean is as noisy as the cell mean, so at the default spread the two overlap
  # and the test would be about the sample size rather than about the column.
  sim <- simulate_factorial_groups(
    n_feats = 20, n_per_cell = 20, n_up = 8, n_down = 8,
    factor_lv = list(a = c("a1", "a2"), b = c("b1", "b2")),
    ref_sd = c(0.6, 1.2), cell_sd = c(0.6, 1.2),
    term_mix = c(crossover = 1), seed = 20260813
  )
  res <- suppressWarnings(suppressMessages(do.call(
    compare_factorial_groups, c(sim$args, list(posthoc = FALSE,
                                               diagnose = FALSE))
  )))
  d <- merge(res$terms[c("features", "terms", "log2_effect")], sim$truth_term,
             by = c("features", "terms"))

  # `truth_term$max_abs_delta` is the same decomposition applied to the shifts
  # that were planted, so the two columns are one quantity measured twice and the
  # answer key is unsigned.
  expect_identical(unique(d$terms[d$is_effect]), "a:b")
  expect_gt(stats::cor(abs(d$log2_effect), d$max_abs_delta), 0.9)
  expect_gt(mean(abs(d$log2_effect[d$is_effect])),
            5 * mean(abs(d$log2_effect[!d$is_effect])))
})

test_that("by = 'term' returns one verdict table per term", {
  res <- sa_factorial_fixture(posthoc = FALSE)
  sig <- estimate_significance(res, by = "term", log2fc_cutoff = 0.1)

  expect_identical(sig$analysis_type, "factorial_comparison")
  expect_identical(names(sig$significance), unique(res$terms$terms))
  for (nm in names(sig$significance)) {
    one <- sig$significance[[nm]]
    expect_identical(names(one),
                     c("features", "log2fc", "pvalue", "adj_pvalue",
                       "is_signif"))
    expect_identical(one$features, res$features)
    at <- res$terms$terms == nm
    expect_equal(one$log2fc, res$terms$log2_effect[at])
    expect_equal(one$pvalue, res$terms$pval[at])
    expect_equal(one$adj_pvalue, res$terms$pval_adj[at])
    # Which term a table is for travels with it, which is what lets the plot
    # label its own panel.
    expect_identical(attr(one, "term"), nm)
    expect_identical(attr(one, "log2fc_cutoff"), 0.1)
  }
  expect_identical(attr(sig$significance[["wool:tension"]], "term_order"), 2L)
  expect_output(print(sig), "one table per term")
})

test_that("by = 'term' adjusts within the term either way", {
  res <- sa_factorial_fixture(posthoc = FALSE)
  sig <- estimate_significance(res, by = "term", adj_type = "bonferroni")
  raw <- res$terms$pval[res$terms$terms == "wool"]

  expect_equal(sig$significance$wool$adj_pvalue,
               stats::p.adjust(raw, "bonferroni"))
  expect_identical(attr(sig$significance$wool, "adj_type"), "bonferroni")
})

test_that("by = 'term' says so when there is no term axis", {
  expect_error(estimate_significance(sa_two_group_fixture(), test = "t_test",
                                     by = "term"),
               "needs a term axis")
  expect_error(estimate_significance(sa_multi_group_fixture(), by = "term"),
               "compare_factorial_groups")
})

test_that("factorial defaults to by = omnibus", {
  res <- sa_factorial_fixture(posthoc = FALSE)
  sig <- estimate_significance(res)
  expect_s3_class(sig$significance, "data.frame")
  expect_identical(sig$significance$features, res$features)
  expect_true("extreme_cell" %in% names(sig$significance))
  expect_identical(sig$significance$extreme_cell, res$effect$extreme_cell)

  by_term <- estimate_significance(res, by = "term")
  expect_identical(names(by_term$significance), unique(res$terms$terms))
  expect_identical(by_term$significance$wool$features, res$features)
})

test_that("the consumers of a comparison read a factorial result", {
  res <- sa_factorial_fixture()
  sig <- estimate_significance(res, test = "anova_test")
  expect_s3_class(sig, "sa_significance")
  expect_identical(sig$analysis_type, "factorial_comparison")
  expect_s3_class(sig$significance, "data.frame")
  expect_identical(sig$significance$features, res$features)

  by_term <- estimate_significance(res, test = "anova_test", by = "term")
  expect_type(by_term$significance, "list")
  expect_identical(names(by_term$significance), unique(res$terms$terms))

  expect_output(print(res), "factors  : wool \\(2\\) x tension \\(3\\)")
  expect_output(print(res), "anova    : two-way")
  expect_output(print(res), "wool:tension")

  skip_if_not_installed("ggplot2")
  expect_s3_class(draw_forest_plot(res, test = "anova_test", print = FALSE),
                  "data.frame")

  # `$pairwise` is the one slot a factorial result does not have, and the reader
  # that needs it says so rather than failing on a NULL.
  expect_error(estimate_significance(res, by = "contrast"),
               "needs a pairwise stage")
})

test_that("the arguments are checked before anything is fitted", {
  bad <- list(data = warpbreaks, feats = "breaks",
              factors = list(wool = "wool", tension = "tension"))
  expect_error(do.call(compare_factorial_groups,
                       c(bad, list(ss_type = "IV"))), "should be one of")
  expect_error(do.call(compare_factorial_groups,
                       c(bad, list(posthoc_scope = "everything"))),
               "should be one of")
  expect_error(do.call(compare_factorial_groups,
                       c(bad, list(conf_level = 1))), "conf_level")
  expect_error(do.call(compare_factorial_groups,
                       c(bad, list(p_adjust = "nope"))), "p_adjust")
  expect_error(do.call(compare_factorial_groups,
                       c(bad, list(posthoc = "yes"))), "posthoc")
  expect_error(
    compare_factorial_groups(warpbreaks, "wool",
                             list(wool = "wool", tension = "tension")),
    "numeric"
  )
})
