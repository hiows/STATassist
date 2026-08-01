# STATassist

**STATassist** runs every test that applies to a comparison scenario in one call and returns standardised, feature-wise tables: parametric, rank-based and robust tests side by side, with fold changes, intervals and multiplicity-adjusted p-values. Significance tables, forest plots, volcano plots, boxplots and back-to-back histograms all read from the same result object, so effect sizes and p-values always describe the same observations.

Two groups, three or more groups and a single sample all return the same object, so `plot()`, `estimate_significance()` and anything else that reads a result works across them without being told which scenario produced it.

Dependencies are limited to base R (`stats`, `graphics`, `grDevices`).

## Example

The volcano plot below is exactly what **Quick start §1–2** produce on `iris` (setosa removed): `group_lv = c("virginica", "versicolor")`, so positive `log2fc` means higher in virginica.

![Volcano plot from compare_two_groups and estimate_significance](man/figures/README-volcano.png)

`plot()` on the same result draws the estimates against their intervals, and on a multi-group result it drops through to the pairwise contrasts, since an omnibus test has no interval of its own:

| Two groups, `plot(res)` | Three groups, `plot(multi, type = "posthoc")` |
| --- | --- |
| ![Forest plot of mean differences](man/figures/README-forest.png) | ![Tukey HSD contrasts for petal length](man/figures/README-posthoc.png) |

Grouped boxplots and a back-to-back histogram use the same wide data interface:

| Boxplot (three species) | Butterfly histogram (two species) |
| --- | --- |
| ![Grouped boxplot](man/figures/README-boxplot.png) | ![Butterfly histogram](man/figures/README-butterfly.png) |

---

## Installation

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hiows/STATassist")          # latest
remotes::install_github("hiows/STATassist@v0.1.0")   # pinned to a release
```

Each release is tagged, so the pinned form keeps returning the same code no matter what lands on the default branch afterwards. `v0.0.1` is the previous release.

The package is not on CRAN yet. When it is submitted, this README will note the CRAN line as well.

---

## Quick start

```r
library(STATassist)
```

### 1. Compare two groups (all applicable tests)

Wide `data.frame`: one row per observation, numeric columns are features. Direction is fixed by `group_lv` (`group_lv[1]` vs `group_lv[2]` for differences and fold change).

```r
iris2 <- iris[iris$Species != "setosa", ]
feats <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")

res <- compare_two_groups(
  data     = iris2,
  feats    = feats,
  group    = iris2$Species,
  group_lv = c("virginica", "versicolor")
)

res
res$effect          # fold_change, log2fc per feature
res$tests$t_test    # Welch or paired t, depending on paired =
res$tests$wilcox_test
res$tests$robust_test
```

Paired example (`sleep`, same subjects under two drugs):

```r
paired_res <- compare_two_groups(
  data     = sleep["extra"],
  feats    = "extra",
  group    = sleep$group,
  group_lv = c("1", "2"),
  id       = sleep$ID,
  paired   = TRUE,
  alternative = "less"
)
paired_res$tests$t_test
```

### 2. Significance and volcano plot

`estimate_significance()` takes the comparison object and applies cutoffs to `log2fc` and p-values. By default it uses the adjusted p-values already stored in the result (`adj_type = NULL` avoids double adjustment).

```r
sig <- estimate_significance(
  res,
  test          = "t_test",
  log2fc_cutoff = 0.25,
  pval_cutoff   = 0.05
)

draw_volcano_plot(sig, main = "Virginica vs versicolor")
```

Use `test = "wilcox_test"` or `test = "robust_test"` to threshold on a different family; `log2fc` stays the same because it comes from `res$effect`.

### 3. Compare three or more groups, with the matching post-hoc stage

Four omnibus tests run side by side, each paired with the pairwise procedure that shares its assumptions: ANOVA with Tukey HSD, Welch's ANOVA with Games-Howell, Yuen's trimmed-mean ANOVA with pairwise Yuen, Kruskal-Wallis with Dunn's test. Pairing them in the result is what makes it impossible to follow a rank-based omnibus test with a parametric comparison by accident.

```r
multi <- compare_multiple_groups(
  data     = iris,
  feats    = feats,
  group    = iris$Species,
  group_lv = levels(iris$Species)   # the first level is the reference
)

multi
multi$tests$anova_test      # F, df, eta_sq, omega_sq, adjusted p
multi$posthoc$anova_test    # one row per feature and pair of levels
multi$posthoc$kruskal_test  # Dunn's, on the same features
```

```
      features n_used     f_stat    eta_sq     pval_adj
1 Sepal.Length    150  119.26450 0.6187057 2.226226e-31
2  Sepal.Width    150   49.16004 0.4007828 4.492017e-17
3 Petal.Length    150 1180.16118 0.9413717 1.142711e-90
4  Petal.Width    150  960.00715 0.9288829 8.338892e-85
```

An omnibus test reports that the levels are not all alike, not by how much, so its `lower_conf` and `upper_conf` are `NA` throughout and the intervals live in `$posthoc` instead. `estimate` there reads as `group1 - group2`:

```
      features               contrast estimate     pval_adj
1 Sepal.Length    setosa - versicolor   -0.930 3.386180e-14
2 Sepal.Length     setosa - virginica   -1.582 2.997602e-15
3 Sepal.Length versicolor - virginica   -0.652 8.287558e-09
```

The pairwise stage runs only for features whose omnibus test cleared `posthoc_alpha`. A feature that did not qualify is **absent** from the post-hoc table rather than present with `NA`, because "never asked" and "asked and unanswerable" are different facts; `multi$parameters$n_posthoc` records how many features entered.

Repeated conditions need `id` and a complete rectangle. Subjects missing any condition are dropped whole and listed in `design$unmatched_ids`:

```r
chicks <- ChickWeight[ChickWeight$Time %in% c(0, 6, 12, 18), ]

rm_res <- compare_multiple_groups(
  data     = data.frame(weight = chicks$weight),
  feats    = "weight",
  group    = paste0("day", chicks$Time),
  group_lv = c("day0", "day6", "day12", "day18"),
  id       = chicks$Chick,
  paired   = TRUE
)

# Mauchly's sphericity test and both epsilon corrections sit on the same row
rm_res$tests$anova_test[, c("f_stat", "pval", "mauchly_pval", "gg_eps", "pval_gg")]
```

```
    f_stat         pval mauchly_pval   gg_eps      pval_gg
1 267.5938 2.637852e-57 8.779229e-41 0.392112 5.018735e-24
```

Sphericity is badly violated here, which is exactly why the corrected p-value is reported next to the uncorrected one rather than instead of it. Repeated measures also swap in Friedman as `$tests$kruskal_test`, with Conover's pairwise comparisons behind it.

### 4. Compare one sample against a hypothesised value

```r
one <- compare_one_sample(iris, "Sepal.Length", mu = 5.8)
one$tests$t_test
```

```
      features n_used   center  mu       diff     stderr    t_stat  df
1 Sepal.Length    150 5.843333 5.8 0.04333333 0.06761132 0.6409184 149
    cohens_d      pval  pval_adj lower_conf upper_conf
1 0.05233076 0.5225603 0.5225603   5.709732   5.976934
```

`$tests$wilcox_test` adds the signed-rank test with a Hodges-Lehmann pseudo-median, and `$tests$prop_test` a score test with a Wilson interval for binary features. A feature that is not binary produces an `NA` row and one named warning rather than a silently coerced number:

```r
binary <- transform(mtcars, am = as.numeric(am))
compare_one_sample(binary, "am", mu = 0.5, p = 0.4)$tests$prop_test
```

### 5. Check the assumptions before choosing a test

Each assumption is checked twice, by tests that fail differently: Shapiro-Wilk against Kolmogorov-Smirnov for normality, median-centred Levene against Bartlett for homogeneity of variance.

```r
d <- diagnose_distribution(iris, feats, iris$Species)
d
d$normality   # one row per feature and level
d$variance    # one row per feature
d$summary     # normal_ok / variance_ok flags per feature
```

```
<sa_diagnosis> distribution_diagnosis
  features : 4
  groups   : setosa, versicolor, virginica
  settings : alpha = 0.05, outlier criterion = iqr

  checks
    normality  1 of 4 feature(s) have a group failing Shapiro-Wilk at 0.05
    variance   3 of 4 feature(s) fail Levene at 0.05
    outliers   13 observation(s) flagged across 4 feature(s)

  A failed check never changes which tests run. It changes which of
  them deserves the most weight, and that judgement stays with you.
```

A failed check never blocks an analysis and never swaps one test for another. It changes which member of the reported family deserves the most weight: skewed groups favour the rank-based and robust members, unequal variances favour Welch's and Brunner-Munzel's treatments of the same data.

`screen_outliers()` flags observations and **does not remove them**. `row` is the row number in the original `data`, so a flagged point can be looked up:

```r
screen_outliers(iris, feats, iris$Species)                    # 1.5 x IQR fences
screen_outliers(iris, feats, criterion = "robust_z")          # median and MAD
screen_outliers(iris, feats, criterion = "grubbs", alpha = 0.05)
```

```
      features     group row value    score
1 Sepal.Length virginica 107   4.9 1.962963
2  Sepal.Width    setosa  16   4.4 1.526316
3  Sepal.Width    setosa  42   2.3 1.894737
4  Sepal.Width virginica 118   3.8 1.666667
```

The same checks are attached to every comparison as `$diagnostics` unless you pass `diagnose = FALSE`.

### 6. `plot()`: one method for every scenario

`plot()` reads only the columns the result contract guarantees, which is why one method covers all three scenarios. `type = "auto"`, the default, picks the first view the chosen table can support.

```r
plot(res)                                        # estimates with intervals
plot(res, test = "wilcox_test", sort_by = "pvalue")
plot(multi, type = "posthoc", feature = "Petal.Length")
plot(multi, test = "kruskal_test", type = "pvalue")
plot(res, dark = TRUE)
```

The p-value view is the fallback for a table with no interval to draw, and marks the `alpha` threshold:

![Kruskal-Wallis p-values across the three species](man/figures/README-pvalue.png)

### 7. Grouped boxplot

```r
draw_grouped_boxplot(
  data     = iris,
  feats    = feats,
  group    = iris$Species,
  group_lv = levels(iris$Species),
  ylab     = "cm",
  main     = "Iris by species",
  dark     = TRUE   # or FALSE for a light theme
)
```

### 8. Back-to-back histogram (two groups only)

```r
draw_butterfly_hist(
  data     = iris2,
  feat     = "Petal.Length",
  group    = iris2$Species,
  group_lv = c("versicolor", "virginica"),
  breaks   = 12,
  scale    = "proportion",   # or "count"
  main     = "Petal length"
)
```

The call returns bin summaries and per-group `histogram` objects for further plotting.

### 9. Descriptive statistics

```r
summarize_descriptive_stats(iris, feats)

# By group (one row per feature × level)
summarize_descriptive_stats(
  iris, feats, iris$Species,
  group_lv = c("virginica", "versicolor", "setosa")
)
```

---

## Main functions

| Function | Purpose |
| --- | --- |
| `compare_two_groups()` | Welch / Wilcoxon / robust tests plus fold change for two groups |
| `compare_multiple_groups()` | Four omnibus tests for three or more groups, each with its matching post-hoc stage; independent or repeated |
| `compare_one_sample()` | One-sample t, signed-rank and proportion tests against a hypothesised value |
| `diagnose_distribution()` | Normality, homogeneity of variance and outliers for a set of features |
| `screen_outliers()` | Flag observations by IQR fences, robust z or Grubbs, without removing them |
| `estimate_significance()` | Filter features by log2FC and p-value from any comparison result |
| `plot()` (`sa_comparison`) | Forest plot of estimates, of pairwise contrasts, or of p-values |
| `draw_volcano_plot()` | Volcano plot from `estimate_significance()` output |
| `draw_grouped_boxplot()` | Boxplots for several features × group levels |
| `draw_butterfly_hist()` | Back-to-back histogram for exactly two groups |
| `summarize_descriptive_stats()` | Feature-wise (and optional group-wise) descriptive table |

---

## Author

**Wonseok Oh** ([ORCID: 0009-0002-0687-8466](https://orcid.org/0009-0002-0687-8466))

## License

MIT © 2026 Wonseok Oh. See [LICENSE.md](LICENSE.md) for details.
