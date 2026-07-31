# STATassist

**STATassist** runs every test that applies to a comparison scenario in one call and returns standardised, feature-wise tables: parametric, rank-based and robust tests side by side, with fold changes, intervals and multiplicity-adjusted p-values. Significance tables and volcano, boxplot and back-to-back histogram plots all read from the same result object, so effect sizes and p-values always describe the same observations.

Dependencies are limited to base R (`stats`, `graphics`, `grDevices`).

## Example

The volcano plot below is exactly what **Quick start §1–2** produce on `iris` (setosa removed): `group_lv = c("virginica", "versicolor")`, so positive `log2fc` means higher in virginica.

![Volcano plot from compare_two_groups and estimate_significance](man/figures/README-volcano.png)

Grouped boxplots and a back-to-back histogram use the same wide data interface:

| Boxplot (three species) | Butterfly histogram (two species) |
| --- | --- |
| ![Grouped boxplot](man/figures/README-boxplot.png) | ![Butterfly histogram](man/figures/README-butterfly.png) |

---

## Installation

Install the development version from GitHub (branch `STATassist_v0.0.1`):

```r
# install.packages("remotes")
remotes::install_github("hiows/STATassist", ref = "STATassist_v0.0.1")
```

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

### 3. Grouped boxplot

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

### 4. Back-to-back histogram (two groups only)

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

### 5. Descriptive statistics

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
| `estimate_significance()` | Filter features by log2FC and p-value from a comparison result |
| `draw_volcano_plot()` | Volcano plot from `estimate_significance()` output |
| `draw_grouped_boxplot()` | Boxplots for several features × group levels |
| `draw_butterfly_hist()` | Back-to-back histogram for exactly two groups |
| `summarize_descriptive_stats()` | Feature-wise (and optional group-wise) descriptive table |

---

## Author

**Wonseok Oh** ([ORCID: 0009-0002-0687-8466](https://orcid.org/0009-0002-0687-8466))

## License

MIT © 2026 Wonseok Oh. See [LICENSE.md](LICENSE.md) for details.
