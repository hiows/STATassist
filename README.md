# STATassist

**STATassist** runs every method that applies to a question in one call and returns standardised tables. A comparison reports parametric, rank-based and robust tests side by side, with fold changes, intervals and multiplicity-adjusted p-values. A model reports one row per term with its estimate and whatever inference that model honestly supports. A dimension reduction reports one row per point with its coordinates. A feature selection reports one row per candidate with what it was ranked by and whether it survived. A contingency table reports one row per cell under the null the test was read against. A performance evaluation reports one row per model on the same held-out rows. A clustering reports one row per point with the label it was assigned. Seven result contracts, and everything downstream reads them rather than the engine underneath.

Two groups, three or more groups and a single sample all return the same object, so `draw_forest_plot()`, `estimate_significance()` and anything else that reads a result works across them without being told which scenario produced it. Five models — linear, logistic, penalized, forest and kernel — return the same object too, so `coef()` and `predict(model, newdata = )` are one line each whichever of them was fitted, and `perform_rfe()` and `perform_stepwise()` hand the predictors they kept straight to any of the five.

Every example below runs on simulated data whose answer was planted on purpose, so a verdict can be **scored** rather than trusted: `simulate_two_groups()`, `simulate_multiple_groups()`, `simulate_factorial_groups()`, `simulate_categorical_groups()`, `simulate_regression()` and `simulate_classification()` hand back the effects and coefficients they put in.

The comparison, diagnostic and visualisation functions use base R only (`stats`, `graphics`, `grDevices`, `utils`). The modelling and dimension-reduction functions are built on `caret` and call `glmnet`, `randomForest`, `kernlab`, `Rtsne` and `umap` through it. The two density-based clustering functions call `dbscan` directly.

## Example

The volcano plot below is what **§1–§2** produce on 30 simulated genes, eight of them planted up and eight planted down in `case`. Since the answer is known, **§3** scores the plot against it: 13 of the 16 planted genes are called, and none of the 14 null ones.

![Volcano plot from compare_two_groups and estimate_significance](man/figures/README-volcano.png)

`draw_forest_plot()` on the same result draws either the estimates against their intervals or the p-values against the threshold, from the same table:

| `type = "estimate"` | `type = "pvalue"` |
| --- | --- |
| ![Forest plot of differences in log2 means](man/figures/README-forest-estimate.png) | ![Forest plot of adjusted p-values](man/figures/README-forest-pvalue.png) |

And the same wide input feeds the plots that look at the data instead of at a result:

| Grouped boxplot | Back-to-back histogram |
| --- | --- |
| ![Grouped boxplot of the first ten genes](man/figures/README-boxplot.png) | ![Butterfly histogram of gene_8](man/figures/README-butterfly.png) |

---

## Installation

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hiows/STATassist")          # latest
remotes::install_github("hiows/STATassist@v1.0.0")   # pinned to a release
```

Each release is tagged, so the pinned form keeps returning the same code no matter what lands on the default branch afterwards. `v0.6.0` is the previous release, and [NEWS.md](NEWS.md) says what changed between them: a fifth and sixth and seventh result contract for a table, a scored evaluation and a clustering; factorial and categorical comparisons with their own significance and plots; model evaluation on held-out rows; feature-pair association screens; and a grouped barplot beside the boxplot. Nothing an existing call reads has moved.

The package is not on CRAN yet. When it is submitted, this README will note the CRAN line as well.

---

# Part 1 — Comparison

```r
library(STATassist)
```

### 1. Compare two groups (all applicable tests)

Wide `data.frame`: one row per observation, numeric columns are features. Direction is fixed by `group_lv`, whose first level is the reference: differences read `group_lv[2] - group_lv[1]` and fold changes `group_lv[2] / group_lv[1]`. The same rule holds for three or more groups, so a control named first stays the reference whichever function reads it.

`simulate_two_groups()` returns data in exactly that shape along with the effects it planted. Thirty features, eight moved up and eight moved down in `case`, the other fourteen null. Its `args` element is named after the arguments of `compare_two_groups()`, so it can be spread out as below, or handed over whole with `do.call(compare_two_groups, sim$args)`.

```r
sim <- simulate_two_groups(n_feats = 30, n_up = 8, n_down = 8, seed = 2026)

comp_res <- compare_two_groups(
  data        = sim$args$data,
  feats       = sim$args$feats,
  group       = sim$args$group,
  group_lv    = sim$args$group_lv,
  input_scale = sim$args$input_scale
)

comp_res
comp_res$effect          # fold_change, log2fc per feature
comp_res$tests$t_test    # Welch or paired t, depending on paired =
comp_res$tests$wilcox_test
comp_res$tests$robust_test
```

```
<sa_two_group> two_group_comparison
  groups   : control vs case  (independent)
  features : 30
  settings : alternative = two.sided, conf_level = 0.95, p_adjust = BH

  tests
    $t_test       13 of 30 at pval_adj <= 0.05
                 Welch's t-test
    $wilcox_test  11 of 30 at pval_adj <= 0.05
                 Wilcoxon rank sum test (Mann-Whitney U test)
    $robust_test  12 of 30 at pval_adj <= 0.05
                 Brunner-Munzel test

  $diagnostics attached
```

`group_lv` is `c("control", "case")`, so `control` is the reference and a positive `log2fc` means higher in `case`, which is where the effects were planted.

```r
head(comp_res$effect, 4)
```

```
  features   x_center   y_center fold_change     log2fc
1   gene_1 116.645989 580.265253   0.2010218 -2.3145758
2   gene_2 273.339540 210.363380   1.2993685  0.3778106
3   gene_3   5.534795   9.386638   0.5896461 -0.7620787
4   gene_4 106.679403  37.811073   2.8213799  1.4964010
```

The centres are in the hundreds while the features themselves run from about 1 to 15, because the data is on the log2 scale, as gene expression usually is, which is what `input_scale = "log2"` says. Dividing two means of logged values is not a fold change and can even come out with the wrong sign: log2 centres of -1 and -2 are a two-fold increase, but their ratio reads as a two-fold decrease. Each observation is raised back through `2^x` before the centres are taken, and `fc_mean` then defaults to `"geom"`, which makes `log2fc` the difference of the two log2 means.

Only `comp_res$effect` is converted. The tests still run on the log2 values, which is the reason for logging them in the first place.

```r
# The same quantity reached from the raw side, since exp(mean(log(2^v)))
# is 2^mean(v).
sim_raw <- sim$args$data
sim_raw[] <- 2^sim_raw

res_geom <- compare_two_groups(
  data     = sim_raw,
  feats    = sim$args$feats,
  group    = sim$args$group,
  group_lv = sim$args$group_lv,
  fc_mean  = "geom"
)

all.equal(comp_res$effect, res_geom$effect)
#> [1] TRUE
```

Note that this is the geometric mean fold change. On raw data `"arith"`, the default there, is a different centre and gives a different number.

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

`estimate_significance()` takes the comparison object and applies cutoffs to `log2fc` and p-values. `adj_type = NULL`, the default, uses the adjusted p-values already stored in the result and so avoids double adjustment; naming a method re-adjusts from `pval`.

```r
sig <- estimate_significance(
  comp_res,
  test          = "t_test",
  log2fc_cutoff = 1,
  pval_cutoff   = 0.05,
  adj_type      = "BH"
)
sig

verdict <- sig$significance   # one row per feature
draw_volcano_plot(sig, xlim = c(-3, 3))
```

```
<sa_significance> two_group_comparison
  test     : t_test  (Welch's t-test)
  cutoffs  : abs(log2fc) >= 1, adj_pvalue <= 0.05  (BH)
  verdict  : 13 of 30 significant
```

The verdict comes back as `$significance`, a data.frame of `features`, `log2fc`, `pvalue`, `adj_pvalue` and `is_signif`, beside the `$analysis_type` it was read from. The scenario name travels with the table because `log2fc` does not mean the same thing in all three: with two groups it is the second level over the reference, with three or more it is the level furthest from the reference, which is why a multi-group volcano plot says so on its x axis.

Pass `test = "wilcox_test"` or `test = "robust_test"` to threshold on a different family; `log2fc` stays the same because it comes from `comp_res$effect`.

### 3. Score the verdict against the planted answer

The comparison above ran on data whose answer is known, so the verdict can be scored rather than trusted. Unplanted features have a true fold change of exactly zero, which makes anything called among them a false positive by definition.

```r
planted <- sim$truth$direction != "none"
table(planted = planted, called = verdict$is_signif %in% TRUE)
```

```
       called
planted FALSE TRUE
  FALSE    14    0
  TRUE      3   13
```

Thirteen of the sixteen planted features come back and none of the fourteen null ones is called. The three that were missed are worth looking up rather than guessing at, which is what the rest of `truth` is for:

```r
missed <- planted & !(verdict$is_signif %in% TRUE)
sim$truth[missed, ]
verdict[missed, ]
```

```
   features direction    log2fc baseline  sd_case sd_control
3    gene_3      down -1.091246 3.401400 2.823662   1.879291
23  gene_23      down -1.849667 4.325627 3.141577   2.138058
26  gene_26        up  1.212882 3.855122 2.802282   2.373063
   features     log2fc     pvalue adj_pvalue is_signif
3    gene_3 -0.7620787 0.09393639  0.2012923     FALSE
23  gene_23 -0.8510372 0.15229026  0.3045805     FALSE
26  gene_26  0.5736033 0.31638180  0.4745727     FALSE
```

`gene_3` and `gene_26` were planted at about 1.1 and 1.2, barely over the cutoff, and both were estimated under 0.8. Nothing went wrong: an estimate carries a sampling error of its own, so a feature planted near the cutoff lands below it a good share of the time, and the same noise keeps its p-value from clearing 0.05 either. That is the third reason a real volcano plot loses features, next to the p-value cutoff and the multiplicity adjustment, and it is the one a simulation that recovers everything would hide.

`gene_23` is the more interesting miss. It was planted at 1.85, well clear of the cutoff, and still came back at 0.85 — and its `sd_case` of 3.14 against a `sd_control` of 2.14 is why. A group whose spread was widened along with its centre is harder to distinguish, not easier, and `truth` records both so the two can be read together.

```r
## Recall differs between the three families on the same data and the same truth
vapply(names(comp_res$tests), function(nm) {
  hit <- estimate_significance(comp_res, test = nm)$significance$is_signif
  mean(hit[planted] %in% TRUE)
}, numeric(1))
```

```
     t_test wilcox_test robust_test 
     0.8125      0.6875      0.7500 
```

### 4. Grouped boxplot and back-to-back histogram

Both read the same wide input the comparison took, and both draw the levels in the order `group_lv` gives them, so the reference lands on the left.

```r
first_ten <- paste0("gene_", 1:10)

draw_grouped_boxplot(
  data     = sim$args$data,
  feats    = first_ten,
  group    = sim$args$group,
  group_lv = sim$args$group_lv,
  ylim     = c(-5, 20)
)
```

![Grouped boxplot of the first ten genes](man/figures/README-boxplot.png)

```r
draw_butterfly_hist(
  data     = sim$args$data,
  feat     = "gene_8",
  group    = sim$args$group,
  group_lv = sim$args$group_lv,
  breaks   = seq(-5, 20, by = 1),
  type     = "both"        # or "freq" for bars only, "dens" for the curve only
)
```

![Butterfly histogram of gene_8, bars and density](man/figures/README-butterfly.png)

The call returns bin summaries and per-group `histogram` objects for further plotting, plus per-group `density` objects when a density is drawn. A density and a bar can only be read against one axis when the bar is a density too: a count or a proportion per bin scales with the bin width, which the curve knows nothing about. So `type = "dens"` and `type = "both"` move `scale` to `"density"`, and reject a `scale` that was asked for explicitly and says otherwise rather than drawing two incomparable shapes.

```r
drawn <- draw_butterfly_hist(
  data        = sim$args$data,
  feat        = "gene_8",
  group       = sim$args$group,
  group_lv    = sim$args$group_lv,
  breaks      = seq(-5, 20, by = 1),
  type        = "both",
  dens_adjust = 1.8,                      # smooth the shape further
  dens_col    = c("#08306B", "#67000D"),  # one outline colour per level
  dens_alpha  = 0.45                      # fill opacity, so the bars show through
)
drawn$group_densities$case
```

### 5. `draw_forest_plot()`: one function for every scenario

`draw_forest_plot()` reads only the columns the result contract guarantees, which is why one function covers all three scenarios. `type = "auto"`, the default, picks the first view the chosen table can support. `plot()` on a `sa_comparison` is the same function under the name R users reach for first, so the first two lines below are interchangeable.

```r
draw_forest_plot(comp_res)                       # estimates with intervals
plot(comp_res)                                   # the same call

draw_forest_plot(comp_res, test = "wilcox_test", sort_by = "pvalue")
draw_forest_plot(comp_res, dark = TRUE)
```

`feats` picks the features to draw and the order to draw them in, from the top of the plot down, and `sort_by` reorders whatever `feats` selected. `xlim` fixes the axis instead of deriving it, so two plots can be read against each other.

```r
draw_forest_plot(
  comp_res, test = "t_test", type = "estimate",
  feats = first_ten, sort_by = "pvalue", xlim = c(-6, 6)
)
```

![Forest plot of mean differences for the first ten genes](man/figures/README-forest-estimate.png)

The p-value view is the fallback for a table with no interval to draw, and marks the `alpha` threshold. It is also worth asking for on purpose, since it puts the whole selection on one scale:

```r
draw_forest_plot(
  comp_res, test = "t_test", type = "pvalue",
  feats = first_ten, sort_by = "pvalue"
)
```

![Forest plot of adjusted p-values for the first ten genes](man/figures/README-forest-pvalue.png)

`use_adjusted = FALSE` reads the `pval` column instead of `pval_adj`. The colouring, the sorting, the p-value view and the legend all follow it, so the plot always names the p-value it actually used.

### 6. Clustered heatmap

The same wide input, transposed so that features run down the rows. `group` labels the samples, which become the columns, and the strip above them is drawn from it. `stats::heatmap()` draws the cells, the strip and the trees; the upright colour key and the group legend beside it are added afterwards, since it draws none:

```r
drawn <- draw_heatmap(
  data              = sim$args$data,   # a matrix works too
  group             = sim$args$group,
  group_lv          = sim$args$group_lv,
  scale             = "feature",       # or "sample" / "none"
  hclust_method     = "ward.D2",
  show_sample_names = FALSE            # 100 samples, no room for 100 labels
)
```

![Clustered heatmap of the simulated two-group data](man/figures/README-heatmap.png)

Nothing in the plot was told where the groups are, and it still puts the 17 leftmost columns in `control` and 29 of the 50 `case` samples together in one run. It is not two clean blocks, and it should not be: fourteen of these genes were planted with no effect at all, and the planted half carries enough noise that **§3** misses three of them.

Features are z-scored across the samples by default. One colour scale is shared by every cell and features are not measured on a common scale, so without it a single high-abundance feature takes the whole range and the rest of the plot is left white. The plot does not name which `scale` ran, but the numbers beside the colour key are whatever it produced, and `matrix` on the result is the scaled data.

`feats` picks the features to draw, and their order when they are not clustered:

```r
planted_feats <- sim$truth$features[sim$truth$direction != "none"]
draw_heatmap(
  data              = sim$args$data,
  group             = sim$args$group,
  group_lv          = sim$args$group_lv,
  feats             = planted_feats,
  cluster_feats     = FALSE,           # keep the order `feats` names
  dist_method       = "correlation",   # group samples by profile shape
  hclust_method     = "ward.D2",
  show_sample_names = FALSE
)
```

The clustering comes back on the result rather than staying inside the picture, so what the plot claims can be checked:

```r
rownames(drawn$matrix)          # features, top to bottom as drawn
drawn$feat_hclust               # the hclust object behind the row dendrogram
rle(as.character(sim$args$group)[drawn$sample_order])
```

Missing values are drawn as grey cells rather than dropped, and a feature with no variance is centred instead of divided by zero. When a pair of features shares no observed sample there is no distance between them at all, in which case that axis keeps its input order and says so rather than failing.

### 7. Compare three or more groups, with the matching post-hoc stage

Four omnibus tests run side by side, each paired with the pairwise procedure that shares its assumptions: ANOVA with Tukey HSD, Welch's ANOVA with Games-Howell, Yuen's trimmed-mean ANOVA with pairwise Yuen, Kruskal-Wallis with Dunn's test. Pairing them in the result is what makes it impossible to follow a rank-based omnibus test with a parametric comparison by accident.

`simulate_multiple_groups()` builds one control group and any number of treatment groups, and `n_treat` states how many by its length.

```r
sim_multi <- simulate_multiple_groups(
  n_feats   = 10,
  n_control = 50,
  n_treat   = c(50, 50, 50),
  n_up      = 3,
  n_down    = 3,
  seed      = 2026
)

multi <- compare_multiple_groups(
  data        = sim_multi$args$data,
  feats       = sim_multi$args$feats,
  group       = sim_multi$args$group,
  group_lv    = sim_multi$args$group_lv,
  input_scale = sim_multi$args$input_scale
)

multi
multi$tests$anova_test[, c("features", "n_used", "f_stat", "eta_sq", "pval_adj")]
```

```
   features n_used     f_stat      eta_sq     pval_adj
1    prot_1    200 15.5709484 0.192461364 4.007090e-08
2    prot_2    200  5.3653151 0.075889925 3.575055e-03
3    prot_3    200  4.3704284 0.062700036 1.054867e-02
4    prot_4    200  0.1977247 0.003017267 8.978535e-01
5    prot_5    200  0.4906963 0.007454669 7.657082e-01
6    prot_6    200  1.2876072 0.019327364 3.498056e-01
7    prot_7    200 13.3218691 0.169370477 2.988271e-07
8    prot_8    200  1.3413981 0.020118538 3.498056e-01
9    prot_9    200  2.5466629 0.037517133 9.533879e-02
10  prot_10    200 10.5506865 0.139037001 6.104458e-06
```

An omnibus test reports that the levels are not all alike, not by how much, so its `lower_conf` and `upper_conf` are `NA` throughout and the intervals live in `$posthoc` instead. `estimate` there reads as `group1 - group2`, and the reference is the level being subtracted, so a contrast against it points the same way the fold change does:

```r
ph <- multi$posthoc$anova_test
ph[ph$features == "prot_1", c("features", "contrast", "estimate", "pval_adj")]
```

```
  features          contrast   estimate     pval_adj
1   prot_1 treat_1 - control  0.8950192 1.524659e-01
2   prot_1 treat_2 - control  0.2647502 9.239116e-01
3   prot_1 treat_3 - control  2.6276154 1.918369e-08
4   prot_1 treat_2 - treat_1 -0.6302689 4.465873e-01
5   prot_1 treat_3 - treat_1  1.7325962 3.630417e-04
6   prot_1 treat_3 - treat_2  2.3628651 4.752954e-07
```

`draw_forest_plot()` reaches the same rows with `type = "posthoc"`, which is what `type = "auto"` falls through to on an omnibus table:

```r
draw_forest_plot(multi, test = "anova_test", type = "posthoc", feats = "prot_1",
                 sort_by = "pvalue")
```

![Tukey HSD contrasts for prot_1](man/figures/README-multi-posthoc.png)

Only `treat_3` moved this feature, which is one of the three shapes `simulate_multiple_groups()` plants: `"all"` moves every treatment group alike, `"gradient"` moves them in a ramp, and `"single"` moves one and leaves the rest at exactly zero. They are recovered at visibly different rates by the same omnibus test, which is the point of planting more than one.

The pairwise stage runs only for features whose omnibus test cleared `posthoc_alpha`. A feature that did not qualify is **absent** from the post-hoc table rather than present with `NA`, because "never asked" and "asked and unanswerable" are different facts; `multi$parameters$n_posthoc` records how many features entered.

`$pairwise` holds the same numbers one contrast at a time, keyed by test and then by contrast label:

```r
names(multi$pairwise$anova_test)
multi$pairwise$anova_test[["treat_3 - control"]][
  , c("features", "log2fc", "estimate", "pval_adj")
]
```

```
[1] "treat_1 - control" "treat_2 - control" "treat_3 - control"
[4] "treat_2 - treat_1" "treat_3 - treat_1" "treat_3 - treat_2"
   features       log2fc   estimate     pval_adj
1    prot_1  2.627615364  2.6276154 1.918369e-08
2    prot_2  1.562213521  1.5622135 4.409768e-03
3    prot_3 -1.459508771 -1.4595088 9.464962e-03
4    prot_4 -0.269499199         NA           NA
5    prot_5  0.007238073         NA           NA
6    prot_6 -0.383242425         NA           NA
7    prot_7 -1.743104602 -1.7431046 1.463146e-04
8    prot_8 -0.610917723         NA           NA
9    prot_9  1.147411959         NA           NA
10  prot_10  0.178635827  0.1786358 9.701071e-01
```

These tables are rectangular where `$posthoc` is ragged: each holds every feature, in the order the rest of the object uses, so a feature that did not qualify is present with its inference columns `NA`. They add `log2fc`, which no post-hoc procedure reports, being the ratio of the two group centres rather than anything a test produced. It divides `group1` by `group2`, so it agrees in sign with the `estimate` beside it, and it is filled even where the test was never run.

The omnibus verdict and its volcano plot work the same way they did with two groups, except that `log2fc` is now the level furthest from the reference, which the x axis says:

```r
sig_multi <- estimate_significance(multi, test = "anova_test",
                                   pval_cutoff = 0.05, adj_type = "BH")
draw_volcano_plot(sig_multi, xlim = c(-4, 4))

draw_grouped_boxplot(
  data     = sim_multi$args$data,
  feats    = sim_multi$args$feats,
  group    = sim_multi$args$group,
  group_lv = sim_multi$args$group_lv,
  ylim     = c(-10, 20)
)
```

| Volcano plot, four groups | Boxplot, four groups |
| --- | --- |
| ![Volcano plot of the multi-group verdict](man/figures/README-multi-volcano.png) | ![Grouped boxplot of ten proteins across four groups](man/figures/README-multi-boxplot.png) |

`estimate_significance(multi, by = "contrast")` reads the `$pairwise` tables instead, and its `$significance` is then one verdict table per contrast: `sig$significance[["treat_3 - control"]]` is what `draw_volcano_plot()` takes.

Repeated conditions need `id` and a complete rectangle; `simulate_multiple_groups(paired = TRUE)` builds one, and subjects missing any condition are dropped whole and listed in `design$unmatched_ids`.

```r
sim_paired <- simulate_multiple_groups(
  n_feats = 10, n_control = 50, n_treat = c(50, 50, 50),
  n_up = 3, n_down = 3, seed = 2026, paired = TRUE
)

rm_res <- compare_multiple_groups(
  data        = sim_paired$args$data,
  feats       = sim_paired$args$feats,
  group       = sim_paired$args$group,
  group_lv    = sim_paired$args$group_lv,
  id          = sim_paired$args$id,
  input_scale = sim_paired$args$input_scale,
  paired      = TRUE
)

# Mauchly's sphericity test and both epsilon corrections sit on the same row
rm_res$tests$anova_test[1:3, c("features", "f_stat", "pval", "mauchly_pval",
                               "gg_eps", "pval_gg")]
```

```
  features     f_stat         pval mauchly_pval    gg_eps      pval_gg
1   prot_1 21.1243904 1.958421e-11 0.0001389954 0.7346618 5.518760e-09
2   prot_2  0.9501076 4.181613e-01 0.0003986818 0.7597238 3.995928e-01
3   prot_3 11.6182657 7.062201e-07 0.0843851300 0.8906821 2.394895e-06
```

Sphericity is violated for the first two features here, which is exactly why the corrected p-value is reported next to the uncorrected one rather than instead of it. Repeated measures also swap in Friedman as `$tests$kruskal_test`, with Conover's pairwise comparisons behind it.

### 8. Factorial crossed design

Crossing two factors asks three questions at once — each main effect and their interaction — and the answer is planted **per model term**, not per cell. `simulate_factorial_groups()` returns `truth_term` beside the wide data so each row of the ANOVA table can be scored.

```r
sim_fact <- simulate_factorial_groups(seed = 2026)

fact_comp <- compare_factorial_groups(
  data          = sim_fact$args$data,
  feats         = sim_fact$args$feats,
  factors       = sim_fact$args$factors,
  factor_lv     = sim_fact$args$factor_lv,
  control_label = list(treatment = "control", sex = "male"),
  input_scale   = sim_fact$args$input_scale
)

fact_comp
head(fact_comp$effect, 3)
fact_comp$terms[fact_comp$terms$features == "prot_1", ]
```

```
<sa_factorial> factorial_comparison
  factors  : treatment (4) x sex (2)  (8 cells, independent)
  anova    : two-way, Type III sums of squares
  features : 100
  settings : alternative = two.sided, conf_level = 0.95, p_adjust = BH

  tests
    $anova_test  22 of 100 at pval_adj <= 0.05
                Two-way ANOVA (Type III sums of squares)
                post-hoc: 65 of 167 contrast(s) over 22 feature(s), Tukey HSD on marginal means and simple effects

  terms
    treatment      13 of 100 at pval_adj <= 0.05
    sex            9 of 100 at pval_adj <= 0.05
    treatment:sex  5 of 100 at pval_adj <= 0.05

  $diagnostics attached
```

```
  features n_used n_cells ref_center   extreme_cell extreme_center fold_change
1   prot_1    160       8 705.347026 treat_C.female      221.60839   0.3141835
2   prot_2    160       8 188.637301 treat_A.female       19.96096   0.1058166
3   prot_3    160       8   8.808372 treat_C.female       18.07013   2.0514720
     log2fc
1 -1.670321
2 -3.240362
3  1.036659
```

The whole-model F in `$tests$anova_test` is the same test `estimate_significance()` and `draw_forest_plot()` already know from multi-group comparisons. Term-wise inference lives in `$terms`, one row per feature and model term, with `pval_adj` corrected **across features within each term** rather than across terms.

`control_label` names the reference level of each factor without rewriting every level list. Here `treatment = "control"` and `sex = "male"` move those levels first, which is what `$effect` uses as the reference cell and what the volcano x-axis names as `most extreme cell vs control.male`.

```r
sig_fact <- estimate_significance(fact_comp)
sig_fact

sig_fact_term <- estimate_significance(fact_comp, by = "term")
draw_volcano_plot(sig_fact_term)

draw_forest_plot(
  fact_comp, type = "pvalue",
  feats = paste0("prot_", 1:20), sort_by = "pvalue"
)
```

```
<sa_significance> factorial_comparison
  test     : anova_test  (Two-way ANOVA (Type III sums of squares))
  cutoffs  : abs(log2fc) >= 1, adj_pvalue <= 0.05  (BH)
  verdict  : 22 of 100 significant
```

![Term-wise volcano plot of the factorial design](man/figures/README-factorial-volcano.png)

| Omnibus p-values, twenty features | Tukey contrasts for prot_14 |
| --- | --- |
| ![Forest plot of adjusted p-values across twenty proteins](man/figures/README-factorial-forest-pvalue.png) | ![Forest plot of Tukey contrasts for prot_14](man/figures/README-factorial-forest-estimate.png) |

`prot_14` was planted as a **crossover**: its treatment and sex main effects are exactly zero and only the interaction was moved. That is invisible in a one-factor read of treatment alone, and visible the moment the lines cross:

```r
subset(sim_fact$truth_term, features == "prot_14")

draw_interaction_plot(fact_comp, feats = "prot_14")

draw_grouped_boxplot(
  data          = sim_fact$args$data,
  feats         = "prot_14",
  factors       = sim_fact$args$factors,
  factor_lv     = sim_fact$args$factor_lv,
  control_label = list(treatment = "control", sex = "male"),
  ylim          = c(5, 25)
)
```

```
   features         terms term_order is_within max_abs_delta is_effect
40  prot_14     treatment          1     FALSE      1.682590      TRUE
41  prot_14           sex          1     FALSE      0.000000     FALSE
42  prot_14 treatment:sex          2     FALSE      1.346072      TRUE
```

| Interaction plot, prot_14 | Boxplot, prot_14 across eight cells |
| --- | --- |
| ![Interaction plot of cell means for prot_14](man/figures/README-factorial-interaction.png) | ![Grouped boxplot of prot_14 across the factorial cells](man/figures/README-factorial-boxplot.png) |

The interaction panel is the place to read a crossover: the treatment profile runs one way in `male` and the opposite way in `female`, inside a single feature panel rather than by comparing panels across the page.

### 9. Categorical contingency tables

A contingency table has no feature axis. The question is about the **table as a whole**, or about one cell at a time, and the result is `sa_categorical` rather than `sa_comparison`. `design$null` names the hypothesis the expected counts and residuals were read under — independence here — and every downstream function reads that same null.

```r
sim_cat <- simulate_categorical_groups(seed = 2026)

cat_comp <- compare_categorical_groups(
  data          = sim_cat$args$data,
  category_lv   = sim_cat$args$category_lv,
  control_label = list(cat_1 = "n", cat_2 = "mid"),
  paired        = sim_cat$args$paired
)

cat_comp
cat_comp$association
cat_comp$tests
```

```
<sa_categorical> categorical_comparison
  table    : cat_1 (2) x cat_2 (3)  (6 cells, independent)
  null     : independence -- a cell is expected at the product of its margins
  observed : 200 row(s)
  settings : conf_level = 0.95, correct = TRUE

  tests
    $chisq_test   pval = 5.52e-07  (null rejected at 0.05)
                 Chi-square test of independence
    $fisher_test  pval = 3.51e-07  (null rejected at 0.05)
                 Fisher's exact test

  association
    cramers_v                0.38
    contingency_coefficient  0.355

  $diagnostics attached, rule expected_count_min: met
```

```
                  measure  estimate lower_conf upper_conf
1               cramers_v 0.3796113         NA         NA
2 contingency_coefficient 0.3549002         NA         NA
```

```
$chisq_test
  n_used statistic df         pval lower_conf upper_conf
1    200  28.82095  2 5.515821e-07         NA         NA

$fisher_test
  n_used statistic df         pval lower_conf upper_conf odds_ratio_cond
1    200        NA NA 3.507415e-07         NA         NA              NA
```

`estimate_significance()` is refused here with a pointer to `estimate_categorical_significance()`, for the same reason `diagnose_distribution()` refuses a comparison it cannot read. The cell reading scores each `(row_level, col_level)` pair by how far `observed / expected` sits from one, with p-values from the Pearson residual:

```r
sig_cell <- estimate_categorical_significance(cat_comp, by = "cell")
head(sig_cell$significance, 4)

sig_table <- estimate_categorical_significance(
  cat_comp, by = "table", test = "chisq_test"
)
sig_table
```

```
  row_level col_level observed expected      lift  log2_lift std_residual
1         n       mid       42   28.025 1.4986619  0.5836750     4.339151
2         y       mid       17   30.975 0.5488297 -0.8655695    -4.339151
3         n      high       20   37.050 0.5398111 -0.8894735    -4.949778
4         y      high       58   40.950 1.4163614  0.5021894     4.949778
        pvalue   adj_pvalue is_signif
1 1.430340e-05 2.145510e-05     FALSE
2 1.430340e-05 2.145510e-05     FALSE
3 7.429824e-07 2.228947e-06     FALSE
4 7.429824e-07 2.228947e-06     FALSE
```

```
<sa_categorical_significance> categorical_comparison
  reading  : table  (2 x 3 table)
  null     : independence -- a cell is expected at the product of its margins
  test     : chisq_test  (Chi-square test of independence)
  cutoffs  : pvalue <= 0.05
  verdict  : cramers_v = 0.38  (significant)
```

At the default cutoffs every cell misses — the omnibus test already rejected independence, but no single cell clears both a fold-change and an adjusted p-value at once. That is a different verdict from the table reading, which is one row and one association measure.

`draw_mosaic_plot()` shades each tile by the residual under the same null, and draws the expected conditional proportion as a dashed line inside each strip so the eye reads distance from the null rather than distance from the neighbouring strip:

```r
draw_mosaic_plot(cat_comp)
```

![Mosaic plot shaded by Pearson residuals under independence](man/figures/README-mosaic.png)

`plot()` on an `sa_categorical` is the same function under the name R users reach for first.

### 10. Compare one sample against a hypothesised value

```r
one <- compare_one_sample(sim$args$data, "gene_8", mu = 8)
one$tests$t_test
```

```
  features n_used   center mu     diff    stderr   t_stat df cohens_d
1   gene_8    100 11.41152  8 3.411517 0.2390404 14.27172 99 1.427172
          pval     pval_adj lower_conf upper_conf
1 9.163982e-26 9.163982e-26   10.93721   11.88583
```

`$tests$wilcox_test` adds the signed-rank test with a Hodges-Lehmann pseudo-median, and `$tests$prop_test` a score test with a Wilson interval for binary features. `gene_8` is not binary, so the call above also emits one named warning and leaves that row `NA` rather than coercing a number out of it:

```r
flag <- data.frame(is_case = as.numeric(sim$args$group == "case"))
compare_one_sample(flag, "is_case", mu = 0.5, p = 0.5)$tests$prop_test
```

```
  features n_used n_success proportion   p diff chi_sq df cohens_h pval
1  is_case    100        50        0.5 0.5    0      0  1        0    1
  pval_adj lower_conf upper_conf
1        1  0.4038315  0.5961685
```

### 11. Check the assumptions before choosing a test

Each assumption is checked twice, by tests that fail differently: Shapiro-Wilk against Kolmogorov-Smirnov for normality, median-centred Levene against Bartlett for homogeneity of variance.

```r
d <- diagnose_distribution(sim$args$data, sim$args$feats, sim$args$group)
d
d$normality   # one row per feature and level
d$variance    # one row per feature
d$summary     # normal_ok / variance_ok flags per feature
```

```
<sa_diagnosis> distribution_diagnosis
  features : 30
  groups   : case, control
  settings : alpha = 0.05, outlier criterion = iqr

  checks
    normality  2 of 30 feature(s) have a group failing Shapiro-Wilk at 0.05
    variance   15 of 30 feature(s) fail Levene at 0.05
    outliers   47 observation(s) flagged across 17 feature(s)

  A failed check never changes which tests run. It changes which of
  them deserves the most weight, and that judgement stays with you.
```

Half of these features fail the variance check, which is not a flaw in the data: `simulate_two_groups()` widens the spread of a group along with its centre, and `gene_23` in **§3** is what that costs. A failed check never blocks an analysis and never swaps one test for another. It changes which member of the reported family deserves the most weight: skewed groups favour the rank-based and robust members, unequal variances favour Welch's and Brunner-Munzel's treatments of the same data.

`screen_outliers()` flags observations and **does not remove them**. `row` is the row number in the original `data`, so a flagged point can be looked up:

```r
screen_outliers(sim$args$data, first_ten, sim$args$group)      # 1.5 x IQR fences
screen_outliers(sim$args$data, first_ten, criterion = "robust_z")
screen_outliers(sim$args$data, first_ten, criterion = "grubbs", alpha = 0.05)
```

```
  features   group row     value    score
1   gene_2    case  69 12.962802 1.657094
2   gene_4    case  59  1.849220 1.977710
3   gene_4    case  74  2.335740 1.731842
4   gene_4 control  18  9.149657 1.613356
```

The same checks are attached to every comparison as `$diagnostics` unless you pass `diagnose = FALSE`.

### 12. Descriptive statistics

```r
summarize_descriptive_stats(sim$args$data, paste0("gene_", 1:3))

# By group (one row per feature x level)
summarize_descriptive_stats(sim$args$data, "gene_8", sim$args$group,
                            group_lv = sim$args$group_lv)
```

```
  features   group  n n_missing     mean       sd      var        se        cv
1   gene_8 control 50         0 10.33740 1.375964 1.893278 0.1945907 0.1331054
2   gene_8    case 50         0 12.48563 2.701269 7.296853 0.3820171 0.2163502
       min        q1   median       q3      max      iqr out_lower_bound
1 5.939601  9.461326 10.37261 11.42933 13.81878 1.968002        6.509322
2 6.528890 11.489627 12.72451 14.10401 19.53043 2.614387        7.568047
  out_upper_bound      mad   skewness excess_kurtosis
1        14.38133 1.514535 -0.2724662       1.1000757
2        18.02560 2.036055 -0.1622755       0.3640198
```

### 13. Grouped barplot

`draw_grouped_boxplot()` shows the spread a group's observations have. `draw_grouped_barplot()` shows one number standing for them — a mean, a median, or any column `summarize_descriptive_stats()` already computed — with an error bar whose meaning follows the height.

```r
sim_bar <- simulate_two_groups(
  n_feats = 10, n_up = 3, n_down = 3, seed = 2026
)

draw_grouped_barplot(
  data          = sim_bar$args$data,
  feats         = sim_bar$args$feats,
  group         = sim_bar$args$group,
  group_lv      = sim_bar$args$group_lv,
  control_label = "control",
  errorbar      = "se"
)
```

![Grouped barplot of ten features with standard-error bars](man/figures/README-grouped-barplot.png)

The bar and the table row are the same number because the heights are read from `summarize_descriptive_stats()` rather than recomputed. A median takes only a notch interval; a count or a spread takes none:

```r
draw_grouped_barplot(
  data     = sim_bar$args$data,
  feats    = sim_bar$args$feats,
  group    = sim_bar$args$group,
  group_lv = sim_bar$args$group_lv,
  mainbar  = "median",
  errorbar = "ci"
)
```

### 14. Feature-pair association

Nothing in Part 1 yet asked how **two** features move together. `summarize_association_stats()` is a screen, not a contract: Pearson, Spearman and Kendall come back side by side on the same pairs, each as four square matrices — coefficient, p-value, adjusted p-value and the observations the pair shared.

```r
assoc_cor <- make_block_cor(
  n_features = 10,
  blocks = list(
    list(features = 1:3, cor = 0.9),
    list(features = 4:5, cor = 0.5, against = 6:7)
  )
)

assoc_sim <- simulate_regression(
  n_pred = 10, n_factor_pred = 0, cor_mat = assoc_cor, seed = 2026
)

assoc <- summarize_association_stats(
  data = assoc_sim$args$data[, -1, drop = FALSE],
  feats = colnames(assoc_sim$args$data)[-1]
)

assoc$design
assoc$pearson$corr[1:3, 1:3]
```

```
$feats
 [1] "x_1"  "x_2"  "x_3"  "x_4"  "x_5"  "x_6"  "x_7"  "x_8"  "x_9"  "x_10"

$n_obs
[1] 200

$methods
[1] "pearson"  "spearman" "kendall" 

$adj_type
[1] "BH"

$use
[1] "pairwise.complete.obs"
```

```
          x_1       x_2       x_3
x_1 1.0000000 0.9169472 0.8981318
x_2 0.9169472 1.0000000 0.9275508
x_3 0.8981318 0.9275508 1.0000000
```

The `against` block planted a positive correlation inside `x_1`–`x_3` and `x_4`–`x_5`, and a negative one between those two sides, which is why the upper-left block reads near 0.9 and the cross-block cells read near -0.5.

`draw_corrplot()` is three decisions on top of `draw_heatmap()`: nothing is standardised, the colours are fixed at -1 to 1, and both axes share one clustering order so the diagonal stays diagonal.

```r
draw_corrplot(assoc$pearson$corr)

draw_corrplot(
  assoc$pearson$corr,
  pvalue = assoc$pearson$adj_pvalue
)
```

| All pairs | Pairs that cleared BH at 0.05 |
| --- | --- |
| ![Correlation heatmap of ten predictors](man/figures/README-corrplot.png) | ![Correlation heatmap with non-significant cells blanked](man/figures/README-corrplot-masked.png) |

Blanking happens **after** clustering, so the tree is built on the full matrix the reader is being shown. The distance is `1 - cor()`, the same rule `cluster_hclust()` and `draw_heatmap(dist_method = "correlation")` use.

---

# Part 2 — Modelling

A model has no feature axis. Every table in Part 1 repeats `features` in the same order; a model has one outcome and a set of **terms**, and the terms are not the columns that were handed in, since one factor predictor becomes several. `terms` takes the place of `features`, `coefficients$terms` repeats that order, and the eleven slots of an `sa_model` are the same eleven whichever of the five models produced it. §18 and §19 are the two sections here that are searches rather than fits, and they have an axis of their own again: `candidates`, the columns they were asked to choose between. §23 and §24 are the first sections that **score** fitted models on held-out rows rather than fit or search on training ones.

### 15. Data whose coefficients are known, and a split that does not leak

`make_block_cor()` builds the correlation the predictors are drawn with. Blocks may not overlap, and a matrix that is symmetric with a unit diagonal but describes no data that could exist is refused here rather than inside an engine. A block whose predictors do not all move the same way names the other side as `against`, which is how a negative correlation below `-1/(k - 1)` gets written down at all: that is the floor on one value shared by `k` predictors, so three of them cannot disagree past -0.5 while a split block has no such limit.

```r
cor_mat <- make_block_cor(
  n_features = 8,
  blocks = list(
    list(features = 1:2, cor = 0.8),
    list(features = 3:5, cor = 0.5)
  )
)

sim_reg <- simulate_regression(cor_mat = cor_mat, seed = 2026)
subset(sim_reg$truth, role == "signal")
```

```
  predictors   role       beta direction value_mean value_sd max_cor_signal
1        x_1 signal  1.1993458        up          0        1              0
5        x_5 signal  0.5376967        up          0        1              0
7        x_7 signal -1.7915160      down          0        1              0
8        x_8 signal -0.8787518      down          0        1              0
```

Four predictors carry a coefficient and the other four are exactly zero, so a false positive is a count rather than an estimate. `max_cor_signal` is why the correlation blocks are there at all: `x_2` is null but correlates with the planted `x_1` at 0.8, and a null predictor that correlates with a planted one is pulled off zero by data alone. No number of rows fixes that, and every section below runs into it.

`truth` has one row per predictor; `truth_term` has one row per term, aligned with `coefficients` by position, since a three-level factor is two terms and a constant predictor is none.

`split_data()` defines what "the training half" means, and closes the two ways a training set learns what it must not. `stratified` keeps the balance of the whole data on both sides, and `id` sends every row of one sampling unit to the same side.

```r
dataset <- split_data(
  data       = sim_reg$args$data,
  stratified = sim_reg$args$data$x_cat_1,
  p_train    = 0.75,
  times      = 1,
  seed       = 2026
)
dataset

train_data <- dataset$datasets[[1]]$train_data
test_data  <- dataset$datasets[[1]]$test_data
```

```
<sa_split> train/test partition
  rows     : 200
  stratify : <vector>
             high 66, low 67, mid 67
  settings : p_train = 0.75, times = 1, seed = 2026

  splits
    $Resample1  train 152 / test 48  (p = 0.760)
```

`p_train` is a proportion of rows, or of units when `id` is given, and the row proportion actually reached is reported as `p` above and stored in `parameters$achieved_p`. The shape does not depend on `times`: `datasets` is a list of one when one split was asked for.

`sim_reg$split_args` is named after the arguments of `split_data()` for the same reason `args` is named after the arguments of the model, so either can be handed over with `do.call()`.

### 16. Linear regression

Every model takes `data`, `outcome` and `predictors`, and every one resamples the same way. Cross-validation here scores the fit and does not choose it: the final model is fitted on all usable rows either way, so `cv = TRUE` and `cv = FALSE` give identical coefficients and differ only in `performance` and `resampling`.

```r
lin <- fit_linear_regression(
  data       = train_data,
  outcome    = sim_reg$args$outcome,
  predictors = sim_reg$args$predictors,
  cv         = TRUE,
  cv_method  = "repeated_kfold",
  n_fold     = 10,
  n_repeat   = 3,
  seed       = 2026
)

lin
```

```
<sa_model> linear_regression
  outcome  : y  (continuous)
  rows     : 152 used
  terms    : 11 over 9 predictor(s)
  settings : repeated_kfold, 10 fold(s) x 3 repeat(s), conf_level = 0.95

  coefficients
    (Intercept)      -1.036  [-1.85, -0.226]  p = 0.0125
    x_1              0.4689  [-0.376, 1.31]  p = 0.274
    x_2              0.7246  [-0.127, 1.58]  p = 0.0948
    x_3              0.2052  [-0.382, 0.792]  p = 0.491
    x_4             -0.3861  [-0.924, 0.152]  p = 0.158
    x_5              0.9453  [0.425, 1.47]  p = 0.000451
    x_6             0.08831  [-0.369, 0.545]  p = 0.703
    x_7              -1.442  [-1.88, -1]  p = 1.48e-09
    x_8              -1.197  [-1.72, -0.671]  p = 1.43e-05
    x_cat_1mid        3.135  [1.97, 4.3]  p = 3.71e-07
    ... and 1 more term(s) in $coefficients

  fit      : r_squared = 0.495, adj_r_squared = 0.459, sigma = 2.75, f_stat =
             13.8, df1 = 10, df2 = 141, pval = 9.29e-17, aic = 751, bic = 788
  resample : RMSE = 2.81 (SD 0.5), Rsquared = 0.451 (SD 0.15), MAE = 2.22 (SD
             0.42) over 30 resample(s)
```

`coef()` on the result is the whole table. `coef()` on `$fit` is the named vector `lm` would have given, for indexing and multiplying:

```r
coef(lin)[, c("terms", "estimate", "pval")]
```

```
         terms    estimate         pval
1  (Intercept) -1.03596473 1.250653e-02
2          x_1  0.46894277 2.743564e-01
3          x_2  0.72456532 9.476519e-02
4          x_3  0.20519712 4.907838e-01
5          x_4 -0.38610825 1.582901e-01
6          x_5  0.94533612 4.510358e-04
7          x_6  0.08831122 7.030701e-01
8          x_7 -1.44197889 1.479616e-09
9          x_8 -1.19679727 1.431263e-05
10  x_cat_1mid  3.13513207 3.711674e-07
11 x_cat_1high -0.18221871 7.500329e-01
```

Three of the four planted coefficients are picked out at any usual threshold, and `x_1` is not: it is the one that correlates with `x_2` at 0.8, and the pair splits the effect between them at p = 0.27 and p = 0.09. This is `max_cor_signal` doing exactly what it records, and no amount of `n_fold` changes it.

Refit on the terms that survived a threshold, predict the held-out half, and read the two against each other. A factor has to be named as its column again, since `x_cat_1mid` is a term and `x_cat_1` is a predictor:

```r
terms_kept <- coef(lin)$terms[-1][coef(lin)$pval[-1] < 0.01]
kept <- unique(sub("high$", "", sub("mid$", "", terms_kept)))
kept
#> [1] "x_5"     "x_7"     "x_8"     "x_cat_1"

lin_kept <- fit_linear_regression(
  data       = train_data,
  outcome    = sim_reg$args$outcome,
  predictors = kept,
  cv         = TRUE, cv_method = "repeated_kfold",
  n_fold     = 10, n_repeat = 3, seed = 2026
)

y_hat <- predict(lin_kept, newdata = test_data)
round(cor(test_data$y, y_hat), 3)
#> [1] 0.827
```

![Predicted against observed for the linear model](man/figures/README-linear-regression.png)

Dropping five predictors, one of them planted, cost 0.008 of correlation on the held-out half: 0.835 with all nine against 0.827 with four. The predictor that was lost was the one whose effect its correlated neighbour was already carrying.

`fit_stats` is a named list rather than a table, because these are quantities per model and not per term: `r_squared`, `adj_r_squared`, `sigma`, the F test, `aic` and `bic`.

### 17. Logistic regression

The same call with a two-class outcome. `outcome_lv` follows the `group_lv` rule — the first level is the reference — so the coefficients describe the odds of `outcome_lv[2]`, and a vector handed to both `compare_two_groups()` and this function points the same way in both.

```r
sim_cls <- simulate_classification(cor_mat = cor_mat, seed = 2026)

cls <- split_data(
  data       = sim_cls$args$data,
  stratified = sim_cls$args$data$y,   # about one row in four is an event
  p_train    = 0.75,
  times      = 1,
  seed       = 2026
)
cls_train <- cls$datasets[[1]]$train_data
cls_test  <- cls$datasets[[1]]$test_data

log_fit <- fit_logistic_regression(
  data       = cls_train,
  outcome    = sim_cls$args$outcome,
  predictors = sim_cls$args$predictors,
  outcome_lv = sim_cls$args$outcome_lv,
  cv         = TRUE, cv_method = "repeated_kfold",
  n_fold     = 10, n_repeat = 3, seed = 2026
)

log_fit
coef(log_fit)[, c("terms", "odds_ratio", "pval")]
```

```
<sa_model> logistic_regression
  outcome  : y  (two classes)
             modelling the odds of case against control, 39 of 150 row(s)
  rows     : 150 used
  terms    : 11 over 9 predictor(s)
  settings : repeated_kfold, 10 fold(s) x 3 repeat(s), conf_level = 0.95

  coefficients
    (Intercept)      -3.071  [-4.54, -1.6]  p = 4.28e-05
    x_1               2.081  [0.738, 3.42]  p = 0.00239
    x_2           -0.002994  [-1.12, 1.11]  p = 0.996
    x_3             -0.2176  [-0.972, 0.536]  p = 0.572
    x_4            -0.07222  [-0.829, 0.685]  p = 0.852
    x_5               0.972  [0.152, 1.79]  p = 0.0202
    x_6              0.3154  [-0.37, 1]  p = 0.367
    x_7              -2.377  [-3.39, -1.36]  p = 4.23e-06
    x_8               -0.74  [-1.37, -0.114]  p = 0.0205
    x_cat_1mid        2.756  [1.03, 4.48]  p = 0.00171
    ... and 1 more term(s) in $coefficients

  fit      : null_deviance = 172, residual_deviance = 76.1, df_null = 149,
             df_residual = 139, mcfadden_r2 = 0.558, lr_stat = 95.9, lr_df =
             10, lr_pval = 3.66e-16, aic = 98.1, bic = 131
  resample : Accuracy = 0.836 (SD 0.081), Kappa = 0.533 (SD 0.28) over 30
             resample(s)
```

```
         terms  odds_ratio         pval
1  (Intercept)  0.04637071 4.276342e-05
2          x_1  8.01038947 2.394441e-03
3          x_2  0.99701067 9.957929e-01
4          x_3  0.80444965 5.716053e-01
5          x_4  0.93032378 8.516444e-01
6          x_5  2.64318666 2.022411e-02
7          x_6  1.37075135 3.670408e-01
8          x_7  0.09286615 4.231261e-06
9          x_8  0.47711976 2.054874e-02
10  x_cat_1mid 15.73809479 1.706532e-03
11 x_cat_1high  0.33968065 1.743603e-01
```

All four planted predictors clear 0.05 here, and the two that were planted `down` come back with an odds ratio under one. The interval columns are `or_lower_conf` and `or_upper_conf`, exponentiated from the Wald interval on the log-odds scale rather than profiled, so the two numbers always agree with the standard error in the same row.

`predict(model, newdata = , type = "response")` is the probability of `outcome_lv[2]`, which is the second column of `type = "prob"` and the class in `type = "raw"`. Fitting the significant terms and the rest separately and drawing all three against the held-out half is what the figure below does:

```r
prob_all <- predict(log_fit, newdata = cls_test, type = "response")
roc_all <- pROC::roc(cls_test$y, prob_all, levels = c("control", "case"))
round(as.numeric(pROC::auc(roc_all)), 3)
#> [1] 0.917
```

![ROC curves for the logistic model](man/figures/README-logistic-regression.png)

The five predictors behind the significant terms reach 0.911 against 0.917 for all nine, and the predictors behind the terms that were **not** significant still reach 0.844 — because `x_cat_1` appears in both sets. One of its two levels cleared the threshold and the other did not, so naming the columns behind the terms puts the factor on both sides. A term is not a predictor, and this is where the difference shows.

### 18. Recursive feature elimination

The paragraph above selected predictors, and it did it the way most analyses do: fit once, read the p-values, keep what cleared 0.05. Two things are wrong with that even when the answer comes out right. The threshold is arbitrary, and the p-values that chose the predictors came from all 150 training rows, so the resampled score of the model that follows describes a fit whose predictors were already chosen — the choosing sits outside the resampling that reports on it. `perform_rfe()` asks the same question with the elimination **inside** the resampling: rank the candidates, drop the weakest, score what is left, and repeat until one predictor is standing. There is no `cv` argument, because an elimination with nothing held out has no score to choose a size by.

```r
rfe <- perform_rfe(
  data          = cls_train,
  outcome       = sim_cls$args$outcome,
  predictors    = sim_cls$args$predictors,
  outcome_lv    = sim_cls$args$outcome_lv,
  control_label = "control",
  model         = "logistic",
  seed          = 2026
)

rfe
```

```
<sa_selection> rfe
  outcome  : y  (two classes)
             modelling case against control, 39 of 150 row(s)
  rows     : 150 used
  search   : Binomial logistic regression over 9 candidate(s), size(s) 1, 2,
             3, 4, 5, 6, 7, 8, 9
  settings : repeated_kfold, 5 fold(s) x 5 repeat(s), Accuracy maximised
  selected : 5 of 9  (Accuracy = 0.849 (SD 0.059) over 25 resample(s))

  ranking  (absolute Wald z)
    x_7           4.048  selected
    x_cat_1       2.759  selected
    x_1           2.698  selected
    x_5           2.083  selected
    x_8           2.051  selected
    x_6          0.8837  dropped
    x_3          0.7253  dropped
    x_2          0.4517  dropped
    x_4          0.3858  dropped
```

This is `sa_selection`, the fourth result contract. `candidates` takes the place `features` holds in a comparison, `terms` in a model and `points` in a reduction, and two tables hang off it because "which predictors" and "how many" are two different answers: `ranking` has one row per candidate, `profile` one row per subset size.

Two details of the ranking are worth the space. It is the absolute Wald z rather than the coefficient, so a predictor measured in grams and the same predictor in kilograms are eliminated in the same order — a coefficient is an effect per unit, and ranking by its size ranks by the units. And `x_cat_1` is ranked as one candidate rather than as its two dummy terms, which is exactly the trap §17 ran into: one of its levels cleared 0.05 and the other did not, so naming the columns behind the terms put the factor on both sides of that comparison. Here a factor is kept or dropped as a column, which is the only thing a later `predictors =` could accept.

```r
rfe$selected
#> [1] "x_7"     "x_cat_1" "x_1"     "x_5"     "x_8"

subset(sim_cls$truth, role != "null")$predictors
#> [1] "x_1"     "x_5"     "x_7"     "x_8"     "x_cat_1"
```

The five it kept are the five that were planted, and no threshold was named anywhere in the call. `x_2` is the interesting one at the other end: it is null but correlates with the planted `x_1` at 0.8, which is what `max_cor_signal` in §15 warned about, and it still lands in the bottom four rather than being carried in by that correlation.

How many to keep is a resampled number like any other, and `profile` is where it is kept honest:

```r
rfe$profile
```

```
  n_vars  Accuracy     Kappa AccuracySD   KappaSD chosen
1      1 0.8386934 0.5266630 0.05018270 0.1408640  FALSE
2      2 0.8213482 0.5051976 0.06073051 0.1675061  FALSE
3      3 0.8481187 0.5814007 0.07990726 0.2200090  FALSE
4      4 0.8438087 0.5822900 0.06504897 0.1574801  FALSE
5      5 0.8494090 0.5992209 0.05894860 0.1566444   TRUE
6      6 0.8345124 0.5608401 0.07265709 0.1886907  FALSE
7      7 0.8438947 0.5835724 0.06967129 0.1781846  FALSE
8      8 0.8412310 0.5766993 0.07081240 0.1787119  FALSE
9      9 0.8426103 0.5809903 0.07156702 0.1811506  FALSE
```

Five won by 0.0013 accuracy over three and by 0.0068 over keeping everything, against standard deviations of 0.05 to 0.08 on those same rows. That is a table of near-ties, and reading it is the point of it being a table: the same call with `n_fold = 10, n_repeat = 3` keeps three predictors instead of five and reaches 0.861 on the held-out half rather than 0.911. The search is a better-behaved filter than a p-value threshold, not an oracle, and `profile` is what says which of those two it was on this data.

`$selected` is a set of column names and nothing else, so it goes straight back into a fit:

```r
rfe_fit <- fit_logistic_regression(
  data       = cls_train,
  outcome    = sim_cls$args$outcome,
  predictors = rfe$selected,
  outcome_lv = sim_cls$args$outcome_lv,
  cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)

prob_rfe <- predict(rfe_fit, newdata = cls_test, type = "response")
roc_rfe <- pROC::roc(cls_test$y, prob_rfe, levels = c("control", "case"))
round(as.numeric(pROC::auc(roc_rfe)), 3)
#> [1] 0.911
```

The same 0.911 the significant terms of §17 reached, against 0.917 for all nine. The gain is not a higher AUC — it is arriving there without a threshold, without the term-versus-predictor confusion, and with the size chosen on rows that did not score it. `control_label = "control"` fixes the direction of the search the way it does everywhere else, so the Wald z the ranking uses is about the odds of `case`, and `$fit` is the `caret::rfe()` object for anything the contract does not carry.

### 19. Stepwise selection by information criteria

§18 bought its answer with resampling, and paid the full price: nine subset sizes, twenty-five resamples each. `perform_stepwise()` asks the same question and settles the bill another way. It walks one term at a time — drop the one whose absence costs least, refit, stop when no move helps — and judges every move by an information criterion, which is the likelihood of the model with a flat charge levied against the number of parameters it spent. Nothing is held out, so there is no `cv` argument and no `seed` either: the path is a deterministic consequence of the data and the charge.

```r
step_sel <- perform_stepwise(
  data          = cls_train,
  outcome       = sim_cls$args$outcome,
  predictors    = sim_cls$args$predictors,
  outcome_lv    = sim_cls$args$outcome_lv,
  control_label = "control",
  model         = "logistic",
  criterion     = "AIC"
)

step_sel
```

```
<sa_selection> stepwise
  outcome  : y  (two classes)
             modelling case against control, 39 of 150 row(s)
  rows     : 150 used
  search   : Binomial logistic regression over 9 candidate(s), 4 step(s)
  settings : backward search, AIC minimised at 2 per parameter
  selected : 5 of 9  (AIC = 91.6951)

  ranking  (AIC increase when the predictor is left out)
    x_7           49.64  selected
    x_1           27.64  selected
    x_cat_1       23.25  selected
    x_5           5.164  selected
    x_8           3.841  selected
    x_6         -0.8116  dropped
    x_3          -1.256  dropped
    x_4          -1.785  dropped
    x_2          -1.996  dropped
```

Same contract, same `candidates` axis, and two slots that mean something different here. `ranking$estimate` is what leaving that one predictor out of the selected model would cost the criterion, so unlike §18's absolute Wald z it has a sign, and the sign is the verdict: positive for the five worth their parameters, negative for the four the model is better off without. `parameters$maximize` is `FALSE` for the same reason, since a criterion is a cost and not a score, and `resampling` is `NULL` because nothing was resampled.

```r
step_sel$selected
#> [1] "x_7"     "x_1"     "x_cat_1" "x_5"     "x_8"

subset(sim_cls$truth, role != "null")$predictors
#> [1] "x_1"     "x_5"     "x_7"     "x_8"     "x_cat_1"
```

The five that were planted, again, and `x_cat_1` is again one candidate rather than its two dummy terms. `x_2` is where the two searches read differently: the elimination of §18 left it in the bottom four, and here it is the first thing to go.

That first drop is visible because `profile` is a different table on this side. §18's is a ladder, one row per subset size with all nine of them scored; this one is a path, one row per step, and `step` names the move that reached it:

```r
step_sel$profile
```

```
  n_vars      AIC      BIC  step chosen
1      9 98.05835 131.1753        FALSE
2      8 96.05837 126.1647 - x_2  FALSE
3      7 94.09535 121.1911 - x_4  FALSE
4      6 92.50669 116.5918 - x_3  FALSE
5      5 91.69511 112.7696 - x_6   TRUE
```

The first row is the model the search started from, which is why its `step` is empty, and `chosen` is `TRUE` on the last row because a stepwise search stops where it chose. The four sizes below 5 are absent by construction: this is a record of where the walk went, not a survey of every size, and a size the path never reached has no criterion to report.

Both criteria sit on every row, so the same path is readable on either scale, and here they agree:

```r
perform_stepwise(
  data          = cls_train,
  outcome       = sim_cls$args$outcome,
  predictors    = sim_cls$args$predictors,
  outcome_lv    = sim_cls$args$outcome_lv,
  control_label = "control",
  model         = "logistic",
  criterion     = "BIC"
)$selected
#> [1] "x_7"     "x_1"     "x_cat_1" "x_5"     "x_8"
```

BIC charges `log(150)` = 5.01 per parameter against AIC's 2, and the `BIC` column above falls at every step of the path, so the heavier charge walks the same way and stops in the same place. It is closer than the identical answer suggests: `x_8` is worth 3.841 on the AIC scale and only 0.8303 on the BIC one, so the last predictor in is the first the charge would take. `direction` does not change the answer on this data either, since `"forward"` and `"both"` both arrive at the same five. That is what a well-separated signal looks like, and not something to count on — the AIC path and the BIC path are the same object only until one predictor sits near the charge.

`$selected` is a set of column names and nothing else, so it goes back into a fit the way §18's did:

```r
step_fit <- fit_logistic_regression(
  data       = cls_train,
  outcome    = sim_cls$args$outcome,
  predictors = step_sel$selected,
  outcome_lv = sim_cls$args$outcome_lv,
  cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)

prob_step <- predict(step_fit, newdata = cls_test, type = "response")
roc_step <- pROC::roc(cls_test$y, prob_step, levels = c("control", "case"))
round(as.numeric(pROC::auc(roc_step)), 3)
#> [1] 0.911
```

0.911, which is where §18 landed and where §17's significant terms landed, from a search that never fitted a model to anything but the full 150 rows. That cheapness is also the one thing to hold against the number above it: `AIC = 91.6951` was computed on exactly the rows the model was fitted to, so it ranks the models on this path against each other and claims nothing about a new one. §18's accuracy came from folds the elimination had not seen and is a claim of that kind; this criterion is not, which is why the AUC is measured on `cls_test`. `$fit` is the `stats::step()` result, the selected `glm` with the whole path attached as `$anova`, for anything the contract does not carry.

### 20. Elastic net

One function covers the three corners of one model: `"lasso"` is alpha 1, `"ridge"` is alpha 0, and `"elastic_net"` tunes alpha as well. The outcome type is read from the column, so the same call does regression and classification.

This is the first model where resampling **chooses** rather than scores. `parameters$lambda` and `parameters$alpha` are therefore the values that won, not the grid that was offered; the grid is the rows of `performance`.

```r
enet <- fit_elastic_net(
  data       = train_data,
  outcome    = sim_reg$args$outcome,
  predictors = sim_reg$args$predictors,
  penalty    = "lasso",
  lambda     = c(0.01, 0.1, 0.5, 1, 2),
  cv         = TRUE, cv_method = "repeated_kfold",
  n_fold     = 10, n_repeat = 3, seed = 2026
)

coef(enet)
```

```
         terms   estimate selected
1  (Intercept) -0.9499000     TRUE
2          x_1  0.4392503     TRUE
3          x_2  0.6072695     TRUE
4          x_3  0.0535525     TRUE
5          x_4 -0.1237904     TRUE
6          x_5  0.7752077     TRUE
7          x_6  0.0000000    FALSE
8          x_7 -1.3456974     TRUE
9          x_8 -1.0743007     TRUE
10  x_cat_1mid  2.8170913     TRUE
11 x_cat_1high -0.1891498     TRUE
```

There are no `stderr`, `pval` or interval columns here, and they are **absent rather than `NA`**. A penalized estimate is deliberately biased and the usual standard error assumes an unbiased one, so there is no honest number to put in them; a column of `NA` reads as a table with its values missing, which is a different claim. `selected` takes their place, and `is.null(coef(fit)$pval)` is how a consumer tells the two kinds of table apart. Every term keeps its row either way: a dropped term is `estimate = 0`, not a missing row.

At the winning `lambda = 0.1` the penalty drops only `x_6`, so it is a gentler filter than a p-value at 0.01 was — and on the held-out half the refit reaches 0.836, a shade above the linear model's 0.827.

![Predicted against observed for the penalized model](man/figures/README-elastic-net-regression.png)

The classification path is the same call with `outcome_lv`, and it separates the kept terms from the dropped ones as clearly as the p-value did:

![ROC curves for the penalized model](man/figures/README-elastic-net-roc.png)

The terms LASSO kept reach 0.906 on the held-out half, the same as the full fit, and the ones it dropped reach 0.682.

### 21. Random forest

The first model with no coefficients. A forest holds hundreds of trees and their splits, not one effect per predictor, so `estimate` is **permutation importance** — `%IncMSE` for regression, `MeanDecreaseAccuracy` for classification — and the table is sorted by it, since that is the order worth reading first. `impurity` carries the other measure the same fit reports, because the two disagree in a way worth seeing: permutation is measured on out-of-bag rows, impurity on the splits themselves.

```r
rf <- fit_rf(
  data       = train_data,
  outcome    = sim_reg$args$outcome,
  predictors = sim_reg$args$predictors,
  mtry       = c(2, 5, 8),
  ntree      = 500,
  cv         = TRUE, cv_method = "repeated_kfold",
  n_fold     = 10, n_repeat = 3, seed = 2026
)

rf
coef(rf)
```

```
<sa_model> random_forest
  outcome  : y  (continuous)
  rows     : 152 used
  terms    : 9 over 9 predictor(s)
  settings : repeated_kfold, 10 fold(s) x 3 repeat(s)
  forest   : 500 tree(s), mtry = 5, nodesize = 5  (mtry chosen from 3 candidate(s))

  importance  (permutation)
    x_cat_1       2.835
    x_7           2.707
    x_5           1.429
    x_8          0.9525
    x_2          0.8965
    x_1          0.6468
    x_4        -0.02487
    x_6        -0.09821
    x_3         -0.1449

  fit      : oob_r_squared = 0.327, oob_rmse = 3.06, oob_mae = 2.41, n_oob =
             152
  resample : RMSE = 2.99 (SD 0.49), Rsquared = 0.378 (SD 0.18), MAE = 2.39
             (SD 0.38) over 30 resample(s)
```

```
    terms    estimate  impurity
1 x_cat_1  2.83479166 268.89580
2     x_7  2.70674308 420.78325
3     x_5  1.42927220 323.38487
4     x_8  0.95250001 252.16450
5     x_2  0.89647721 208.29766
6     x_1  0.64680537 174.57364
7     x_4 -0.02486564  95.13807
8     x_6 -0.09820752 156.14747
9     x_3 -0.14494930 120.42547
```

The three predictors at the bottom are **negative**, and that is an answer rather than a missing value: a predictor that carries nothing can do worse than its own permutation. All three are null, and the null `x_2` outranks the planted `x_1` on both measures, which is the 0.8 correlation between the two showing up a third time. Above that the two measures do not agree: `x_cat_1` leads on permutation and sits third on impurity, behind `x_7` and `x_5`.

Importance is not scaled by the between-tree standard deviation, which `randomForest::importance()` does by default. That ratio is referred to no distribution, so the mean loss itself is what the table carries.

`fit_stats` is out-of-bag rather than in-sample, and says so in its names: `oob_r_squared`, `oob_rmse`, `oob_mae`, and for classification `oob_accuracy`, `oob_kappa`, `oob_sensitivity` and `oob_specificity` against `outcome_lv[2]`. A third of the rows are out of bag for each tree and the forest has already predicted them from trees that never saw them, which is an honest held-out score for free. Here it is 3.06 against an in-sample RMSE that would flatter the fit.

A forest splits factors by level directly, so there is no dummy coding and `x_cat_1` is one term rather than two — unlike every other model in this part.

| Regression | Classification |
| --- | --- |
| ![Predicted against observed for the forest](man/figures/README-random-forest-regression.png) | ![ROC curves for the forest](man/figures/README-random-forest-roc.png) |

The forest is the weakest of the five on this data, at 0.747 held-out correlation and 0.841 AUC, which is what a flexible model costs on 152 rows with a mostly linear truth. Its top five and low five separate cleanly all the same: 0.833 against 0.624.

### 22. Support vector machine

The second model with no coefficients, for the opposite reason. A forest has too many numbers per predictor to report one; a radial kernel machine has **none** — it holds support vectors and their weights, which are points in the data rather than directions in the predictor space. So `estimate` is permutation importance again, measured in the metric the resampling tuned on, so the table and `performance` read in the same unit.

```r
svm <- fit_svm(
  data       = train_data,
  outcome    = sim_reg$args$outcome,
  predictors = sim_reg$args$predictors,
  C          = 2^seq(-5, 10, by = 2),
  sigma      = NULL,           # read from the data by kernlab::sigest()
  cv         = TRUE, cv_method = "repeated_kfold",
  n_fold     = 10, n_repeat = 3, seed = 2026
)

svm
names(coef(svm))
#> [1] "terms"    "estimate"
```

```
<sa_model> svm
  outcome  : y  (continuous)
  rows     : 152 used
  terms    : 10 over 9 predictor(s)
  settings : repeated_kfold, 10 fold(s) x 3 repeat(s)
  kernel   : radial, C = 0.5, sigma = 0.0487  (chosen from 8 candidate(s))

  importance  (permutation)
    x_7              0.6094
    x_cat_1mid       0.5242
    x_8               0.292
    x_5              0.2586
    x_2              0.1404
    x_1              0.1066
    x_cat_1high     0.06732
    x_6              0.0615
    x_4             0.05414
    x_3             0.04743

  fit      : r_squared = 0.516, rmse = 2.59, mae = 1.9, n_support_vector =
             136, support_vector_rate = 0.895
  resample : RMSE = 2.95 (SD 0.47), Rsquared = 0.413 (SD 0.16), MAE = 2.33
             (SD 0.38) over 30 resample(s)
```

Unlike the forest, this importance is measured on the rows the machine was fitted to, because a machine sees every row at once and has no out-of-bag half to permute. A term fitted to noise therefore earns a little importance it could not have earned out of sample, and the numbers at the bottom of this table are small rather than negative for that reason. `sigma = NULL` reads the kernel width from the data as the median of `kernlab::sigest()`, and the predictors are centred and scaled before the kernel measures a distance, so `sigma` is a width on the standardised scale.

| Regression | Classification |
| --- | --- |
| ![Predicted against observed for the machine](man/figures/README-svm-regression.png) | ![ROC curves for the machine](man/figures/README-svm-roc.png)  |

### 23. Score regression models on held-out rows

§16–§22 each scored themselves by hand — a correlation or an AUC beside a scatter or an ROC snippet. `evaluate_regression_models()` is the layer above `predict.sa_model()`: the same held-out rows for every model, one table of metrics and one of deltas against a baseline. Rows are the **intersection** of what every model could predict, not the union, because a delta only means something when the two numbers came from the same rows.

The predictors are chosen once by RFE on the training half, then every model is refit on that same set so the comparison is about engines rather than about which columns each engine saw:

```r
eval_rfe <- perform_rfe(
  data       = train_data,
  outcome    = sim_reg$args$outcome,
  predictors = sim_reg$args$predictors,
  seed       = 2026
)
eval_sel <- eval_rfe$selected

eval_lin <- fit_linear_regression(
  data = train_data, outcome = sim_reg$args$outcome,
  predictors = eval_sel, cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)
eval_lasso <- fit_elastic_net(
  data = train_data, outcome = sim_reg$args$outcome,
  predictors = eval_sel, penalty = "lasso",
  cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)
eval_rf <- fit_rf(
  data = train_data, outcome = sim_reg$args$outcome,
  predictors = eval_sel, cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)
eval_svm <- fit_svm(
  data = train_data, outcome = sim_reg$args$outcome,
  predictors = eval_sel, C = 2^seq(-5, 10, by = 2), sigma = NULL,
  cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)

eval_reg <- evaluate_regression_models(
  baseline_model = eval_lin,
  new_models     = list(lasso = eval_lasso, rf = eval_rf, svm = eval_svm),
  newdata        = test_data,
  answer         = test_data$y,
  baseline_label = "linear"
)

eval_reg
eval_reg$metrics
eval_reg$comparisons
```

```
<sa_performance> regression_performance
  outcome  : <vector>  (continuous)
  rows     : 48 scored
  models   : 4, baseline = linear

  metrics
    linear  cor = 0.828, r_squared = 0.68, rmse = 2.14, mae = 1.67
    lasso   cor = 0.828, r_squared = 0.682, rmse = 2.13, mae = 1.65
    rf      cor = 0.761, r_squared = 0.553, rmse = 2.53, mae = 1.96
    svm     cor = 0.827, r_squared = 0.651, rmse = 2.23, mae = 1.75

  comparisons  (against linear)
    lasso  delta_cor = 5.72e-05, delta_r_squared = 0.00204, delta_rmse =
           -0.00681, delta_mae = -0.0121
    rf     delta_cor = -0.0669, delta_r_squared = -0.127, delta_rmse = 0.389,
           delta_mae = 0.297
    svm    delta_cor = -0.000417, delta_r_squared = -0.0293, delta_rmse =
           0.0958, delta_mae = 0.0832
```

```
   model n_used       cor r_squared     rmse      mae       bias calib_slope
1 linear     48 0.8276569 0.6801600 2.135962 1.665212 -0.1449732   0.7331519
2  lasso     48 0.8277140 0.6821956 2.129154 1.653141 -0.1599249   0.7128353
3     rf     48 0.7607476 0.5530111 2.525083 1.962046 -0.2248015   0.4654314
4    svm     48 0.8272403 0.6508357 2.231732 1.748396 -0.3199326   0.5501323
  calib_intercept
1    -0.001149596
2    -0.005151245
3     0.063315928
4    -0.077466584
```

```
  model     delta_cor delta_r_squared   delta_rmse   delta_mae
1 lasso  0.0000571794     0.002035625 -0.006808025 -0.01207140
2    rf -0.0669092668    -0.127148879  0.389120084  0.29683330
3   svm -0.0004165980    -0.029324210  0.095769772  0.08318366
```

There is no p-value beside the deltas: held-out error on rows the caller did not generate has no null the package can name. `cor` and `r_squared` sit together because their gap is what `calib_slope` and `calib_intercept` report — a slope under one is a model whose predictions are squeezed towards their own mean.

`draw_prediction_plot()` reads those calibration numbers rather than refitting, so the line and the table cannot drift apart. With more than one model, `type = "overlay"` and `points = FALSE` compares the calibration lines alone:

```r
draw_prediction_plot(
  eval_reg, type = "overlay",
  anno_corr = TRUE, anno_rsq = TRUE, anno_lm = TRUE, cex.anno = 0.85
)
```

![Predicted against observed for four regression models on the same held-out rows](man/figures/README-eval-regression.png)

`plot()` on an `sa_performance` with `analysis = "regression_performance"` is the same function.

### 24. Score classification models on held-out rows

The same intersection rule and the same RFE-first pipeline, with logistic regression as the baseline. Here the comparisons add three paired questions beside `delta_auc`: DeLong's test on the ranks, the IDI on the probabilities, and the NRI on how often each probability moved the right way.

```r
eval_rfe_cls <- perform_rfe(
  data = cls_train, outcome = sim_cls$args$outcome,
  predictors = sim_cls$args$predictors,
  outcome_lv = sim_cls$args$outcome_lv,
  control_label = "control", seed = 2026, model = "logistic"
)
eval_sel_cls <- eval_rfe_cls$selected

eval_log <- fit_logistic_regression(
  data = cls_train, outcome = sim_cls$args$outcome,
  predictors = eval_sel_cls, outcome_lv = sim_cls$args$outcome_lv,
  control_label = "control",
  cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)
eval_lasso_cls <- fit_elastic_net(
  data = cls_train, outcome = sim_cls$args$outcome,
  predictors = eval_sel_cls, outcome_lv = sim_cls$args$outcome_lv,
  penalty = "lasso", cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)
eval_rf_cls <- fit_rf(
  data = cls_train, outcome = sim_cls$args$outcome,
  predictors = eval_sel_cls, outcome_lv = sim_cls$args$outcome_lv,
  cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)
eval_svm_cls <- fit_svm(
  data = cls_train, outcome = sim_cls$args$outcome,
  predictors = eval_sel_cls, outcome_lv = sim_cls$args$outcome_lv,
  C = 2^seq(-5, 10, by = 2), sigma = NULL,
  cv = TRUE, cv_method = "repeated_kfold",
  n_fold = 10, n_repeat = 3, seed = 2026
)

eval_cls <- evaluate_classification_models(
  baseline_model = eval_log,
  new_models     = list(
    lasso = eval_lasso_cls, rf = eval_rf_cls, svm = eval_svm_cls
  ),
  newdata        = cls_test,
  answer         = cls_test$y,
  outcome_lv     = sim_cls$args$outcome_lv,
  control_label  = "control",
  baseline_label = "logistic"
)

eval_cls
eval_cls$comparisons
```

```
<sa_performance> classification_performance
  outcome  : <vector>  (two classes)
             scoring the probability of case against control, 13 of 50 row(s)
  rows     : 50 scored
  models   : 4, baseline = logistic
  threshold: 0.5  (accuracy, sensitivity and specificity only)

  metrics
    logistic  auc = 0.911 [0.824, 0.997], brier = 0.0965, accuracy = 0.88
    lasso     auc = 0.886 [0.778, 0.994], brier = 0.107, accuracy = 0.86
    rf        auc = 0.852 [0.717, 0.988], brier = 0.111, accuracy = 0.9
    svm       auc = 0.89 [0.769, 1.01], brier = 0.096, accuracy = 0.9

  comparisons  (against logistic)
    lasso
      delta_auc  = -0.0249 [-0.0663, 0.0164] p = 0.237
      IDI        = -0.167 [-0.261, -0.0723] p = 0.000541
      NRI        = -1.18 [-1.72, -0.638] p = 1.89e-05
    rf
      delta_auc  = -0.0582 [-0.151, 0.0342] p = 0.217
      IDI        = -0.136 [-0.244, -0.0274] p = 0.0141
      NRI        = -1.42 [-1.87, -0.972] p = 5.82e-10
    svm
      delta_auc  = -0.0208 [-0.109, 0.0676] p = 0.645
      IDI        = -0.0413 [-0.135, 0.0523] p = 0.387
      NRI        = -1.27 [-1.78, -0.76] p = 1.01e-06
```

```
  model   delta_auc delta_auc_lower_conf delta_auc_upper_conf delta_auc_pval
1 lasso -0.02494802          -0.06632472           0.01642867      0.2373018
2    rf -0.05821206          -0.15061219           0.03418807      0.2169136
3   svm -0.02079002          -0.10922284           0.06764280      0.6449596
          idi idi_lower_conf idi_upper_conf     idi_pval       nri  nri_event
1 -0.16683055     -0.2613504    -0.07231065 0.0005413859 -1.176715 -0.2307692
2 -0.13560688     -0.2438270    -0.02738675 0.0140507709 -1.422037 -0.6923077
3 -0.04127849     -0.1348295     0.05227251 0.3871400524 -1.268191 -0.5384615
  nri_nonevent nri_lower_conf nri_upper_conf     nri_pval
1   -0.9459459      -1.715864     -0.6375667 1.888282e-05
2   -0.7297297      -1.871932     -0.9721430 5.824685e-10
3   -0.7297297      -1.776480     -0.7599029 1.007561e-06
```

The logistic baseline's AUC of 0.911 is the same number §17 reached on the significant terms, now arrived at through one call. LASSO's AUC is lower but its IDI is significantly negative — probabilities moved the wrong way on average even when the ranks barely changed, which is why three statistics are reported rather than one.

```r
draw_roc_curve(eval_cls, anno_auc = TRUE, cex.anno = 1)
```

![ROC curves of four classifiers scored on the same held-out rows](man/figures/README-eval-classification.png)

Every number is about the odds of `case`, because that is what the fitted models predict; `outcome_lv` and `control_label` are read as statements about what was already fit, not as knobs to turn after the fact.

### 25. Predict on held-out data

`predict()` goes to the **result**, not to `$fit`. The engine object knows only the column names it was handed: `glmnet` and `kernlab` were given a design matrix and read it by position, so a frame whose numeric columns are in a different order is multiplied by the wrong coefficients without any error, and a factor predictor is dropped from it entirely. The result object is the only thing that knows which columns were predictors and what the levels of a factor were, so one line covers all five models:

```r
predict(lin,     newdata = test_data)                     # numeric
predict(log_fit, newdata = cls_test, type = "raw")        # factor, outcome_lv
predict(log_fit, newdata = cls_test, type = "prob")       # one column per class
predict(log_fit, newdata = cls_test, type = "response")   # P(outcome_lv[2])
```

Extra columns in `newdata` are ignored, a missing one is an error that names it, and a level the training data never saw is an error that names both the column and the level. A level that is missing from `newdata` is not an error: its dummy column is simply zero, since the levels come from `design$predictor_lv` rather than from the new rows. One prediction comes back per row, and a row with a missing value in a predictor is `NA` rather than dropped, so the answer stays aligned with `newdata`.

---

# Part 3 — Unsupervised learning

Everything above had an answer to score against. These sections have none for the coordinates, and §27 has one for the labels when a grouping was known before any algorithm ran. §26 asks where each point lands when many features are pressed into two dimensions; the three functions answer differently on purpose: PCA is a rotation, so it is reversible and says which feature moved a point, but it only finds straight structure; t-SNE and UMAP find curved structure but cannot say which feature made it. §27 asks which points belong together on the same coordinates, using the same `points` axis the reductions use.

### 26. `perform_pca()`, `perform_tsne()` and `perform_umap()`

Three functions rather than one call with a `methods` argument, because they answer in coordinates that share no scale — nothing but `points` could be joined between them — and because `perplexity`, `n_neighbors` and `metric` each belong to exactly one of them. What makes them comparable is the input: all three read `data` the same way, so the same rows drop for the same reason.

`embedding_scale` chooses which margin becomes the points. The input is one row per sample as everywhere else in the package, and `design$point_type` reports which axis was embedded.

```r
red_cor <- make_block_cor(
  n_features = 8,
  blocks = list(
    list(features = 1:2, cor = 0.8),
    list(features = 3:5, cor = 0.5),
    list(features = 7:8, cor = 0.9)
  )
)
red_data <- simulate_classification(cor_mat = red_cor, seed = 2026)$args$data

pca <- perform_pca(
  data            = red_data,
  feats           = paste0("x_", 1:8),
  embedding_scale = "features",
  center          = TRUE,
  scale           = TRUE
)

pca
pca$scores[, 1:3]
```

```
<sa_reduction> pca
  data     : 200 sample(s) x 8 feature(s)
  points   : 8 feature(s)
  scaling  : centred and scaled
  variance : PC1 25.88%, PC2 24.46%, PC3 22.3% (3 of 8 component(s), 72.64%
             cumulative)
```

```
  points        PC1      PC2
1    x_1 -1.1565840 9.030714
2    x_2  0.2042277 9.636624
3    x_3  9.9079045 4.897448
4    x_4  9.4190232 5.584413
5    x_5 10.1635916 3.756129
6    x_6 -3.0399404 1.583551
7    x_7 -7.8137355 8.413601
8    x_8 -7.0841229 8.507287
```

The three blocks come out as three groups and `x_6`, which belongs to none, sits on its own:

![PCA of the eight features](man/figures/README-pca.png)

PCA is a singular value decomposition, so one fit answers both margins at once and the matrix is never turned around. `embedding_scale = "features"` rescales the rotation from unit length to variance-weighted length and puts the sample scale in `$loadings`; `$variance` and `$fit` are not touched, so the axis labels are the same either way.

```r
by_sample <- perform_pca(data = red_data, feats = paste0("x_", 1:8))
all.equal(pca$variance, by_sample$variance)
#> [1] TRUE
```

Transposing the input by hand instead is a **third** analysis, not the same one: `prcomp()` always centres and scales the columns it is given, so `perform_pca(t(data))` standardises samples rather than features. The shapes match, the picture reads, and the answer is to a different question. This is the one mistake in these three functions that produces a plot instead of an error, which is why all three document it.

`perform_tsne()` sees literally the same matrix `perform_pca()` did, which needs two of `Rtsne`'s defaults turned off: `normalize = TRUE` would overwrite the `center` and `scale` that were asked for, and `pca = TRUE` would show t-SNE a rotation rather than the matrix. Both overrides are recorded in `engine$overridden`.

```r
tsne <- perform_tsne(
  data            = red_data,
  feats           = paste0("x_", 1:8),
  embedding_scale = "features",
  center          = TRUE,
  scale           = TRUE,
  seed            = 2026
)
tsne
```

```
<sa_reduction> tsne
  data     : 200 sample(s) x 8 feature(s)
  points   : 8 feature(s)
  scaling  : centred and scaled
  tsne     : 2 dimension(s), perplexity = 2, theta = 0.5  (seed = 2026)
```

`perform_umap()` is the one that standardises nothing by default, because `metric` is its own argument and `"cosine"` or `"pearson"` compares the shape of a row rather than its size, which already answers what standardising would. Both neighbourhood sizes are read from the engine's limits when they are `NULL`, and both are read from the number of **points** rather than the number of samples — with eight features that is small enough to be worth a message:

```r
umap_res <- perform_umap(
  data            = red_data,
  feats           = paste0("x_", 1:8),
  embedding_scale = "features",
  n_neighbors     = 3,
  center          = FALSE,
  scale           = FALSE,
  seed            = 2026
)
umap_res
#> Only 8 feature(s) to embed (n_neighbors = 3). This method describes a
#> neighbourhood, and below about 16 points there is not much of one to describe.
```

```
<sa_reduction> umap
  data     : 200 sample(s) x 8 feature(s)
  points   : 8 feature(s)
  scaling  : none, values as they arrived
  umap     : 2 dimension(s), method = naive, n_neighbors = 3, min_dist = 0.1, euclidean  (seed = 2026)
```

| t-SNE | UMAP |
| --- | --- |
| ![t-SNE of the eight features](man/figures/README-tsne.png) | ![UMAP of the eight features](man/figures/README-umap.png) |

The coordinates are `$scores` in all three and the engine object is `$fit`. `seed` buys different things in the two stochastic methods: `umap` restores the random stream itself, so two seedless calls agree, while two seedless `Rtsne` calls do not.

### 27. Cluster on an embedding

The four `cluster_*()` functions read their input through the same path as the reductions, so a clustering and an embedding of the same frame are about the same rows. `sa_cluster` assigns every point a label — `0` is noise for the density methods — and `draw_dim_reduction_plot()` is where those labels meet the coordinates.

Colour and shape are two channels on purpose. A clustering is what the data was found to say; a `group` is what was known before either algorithm ran. One colour per shape is a clustering that recovered the groups; one shape split across colours is a group the data does not see as one thing.

```r
clust_sim <- simulate_two_groups(
  n_feats = 50, deg_log2fc = c(5, 10), seed = 2026
)

clust_pca <- perform_pca(
  data            = clust_sim$args$data,
  feats           = clust_sim$args$feats,
  embedding_scale = "samples"
)

clust_km <- cluster_kmeans(
  data          = clust_pca$scores,
  feats         = c("PC1", "PC2"),
  cluster_scale = "samples",
  n_clust       = 2,
  seed          = 2026
)

clust_km
table(clust_sim$args$group, clust_km$assignments$cluster)
```

```
<sa_cluster> kmeans
  data     : 100 sample(s) x 2 feature(s)
  points   : 100 sample(s)
  scaling  : centred and scaled
  clusters : 2
  sizes    : #1 n = 50, s = 0.697; #2 n = 50, s = 0.385
  silhouette: mean 0.541 over the 100 assigned sample(s), on the euclidean
              distance
  kmeans   : k = 2, 25 start(s), 100.21 within-cluster ss  (seed = 2026)
```

```
         cluster
group      1  2
  case     0 50
  control 50  0
```

Fifty features were moved up or down on a log2 scale of five to ten, embedded in the first two principal components of the samples, and cut into two clusters. Every control lands in cluster 1 and every case in cluster 2, which is the grouping that was planted — read through shape on the first plot and through colour on the second:

```r
draw_dim_reduction_plot(
  clust_pca,
  group    = clust_sim$args$group,
  group_lv = clust_sim$args$group_lv,
  col      = c("black", "red3")
)

draw_dim_reduction_plot(
  clust_pca,
  cluster_result = clust_km,
  cluster_lv     = c("Cluster1", "Cluster2")
)
```

| PCA, shape = known group | PCA, colour = k-means |
| --- | --- |
| ![PCA of fifty features coloured by the planted group](man/figures/README-cluster-pca-group.png) | ![PCA of fifty features coloured by k-means clusters](man/figures/README-cluster-pca-cluster.png) |

The same read on a UMAP of the samples rather than a PCA:

```r
clust_umap <- perform_umap(
  data            = clust_sim$args$data,
  feats           = clust_sim$args$feats,
  embedding_scale = "samples",
  seed            = 2026
)
clust_km_umap <- cluster_kmeans(
  data          = clust_umap$scores,
  feats         = c("UMAP1", "UMAP2"),
  cluster_scale = "samples",
  n_clust       = 2,
  seed          = 2026
)

draw_dim_reduction_plot(
  clust_umap,
  group    = clust_sim$args$group,
  group_lv = clust_sim$args$group_lv,
  col      = c("black", "red3")
)

draw_dim_reduction_plot(
  clust_umap,
  cluster_result = clust_km_umap,
  cluster_lv     = c("Cluster1", "Cluster2")
)
```

| UMAP, shape = known group | UMAP, colour = k-means |
| --- | --- |
| ![UMAP of fifty features shaped by the planted group](man/figures/README-cluster-umap-group.png) | ![UMAP of fifty features coloured by k-means clusters](man/figures/README-cluster-umap-cluster.png) |

`plot()` on an `sa_reduction` is the same function. Noise from `cluster_dbscan()` or `cluster_snn()` is drawn grey rather than given a palette colour, because a point left out is the absence of a cluster rather than a cluster of its own.

---

## Main functions

| Function | Purpose |
| --- | --- |
| `compare_two_groups()` | Welch / Wilcoxon / robust tests plus fold change for two groups |
| `compare_multiple_groups()` | Four omnibus tests for three or more groups, each with its matching post-hoc stage; independent or repeated |
| `compare_factorial_groups()` | One two-way, three-way or factorial ANOVA for crossed factors, with an answer per model term and Tukey contrasts on the marginal means and inside each stratum |
| `compare_categorical_groups()` | Chi-square beside Fisher on two categorical variables, or McNemar / Cochran's Q on repeated binary conditions; `design$null` names the hypothesis the expected counts and residuals are read under, and the result is not an `sa_comparison` |
| `compare_one_sample()` | One-sample t, signed-rank and proportion tests against a hypothesised value |
| `diagnose_distribution()` | Normality, homogeneity of variance and outliers for a set of features |
| `screen_outliers()` | Flag observations by IQR fences, robust z or Grubbs, without removing them |
| `estimate_significance()` | Filter features by log2FC and p-value from any comparison result, over the omnibus test, one pairwise contrast at a time, or one model term at a time |
| `estimate_categorical_significance()` | The same verdict for a contingency table, read one cell at a time from `observed / expected` and the cell's standardized residual, or once for the table from an association measure and an omnibus p-value |
| `summarize_descriptive_stats()` | Feature-wise (and optional group-wise) descriptive table |
| `summarize_association_stats()` | Pearson, Spearman and Kendall on every pair of features, each as a square matrix of coefficients beside its p-values, the p-values adjusted across the pairs, and the observations the pair shared |
| `center_by_control()` | Remove the control group's centre from every feature, so each value reads as its distance from the control rather than as a measurement of its own |
| `draw_forest_plot()` | Forest plot of estimates, of pairwise contrasts, or of p-values; `plot()` on a `sa_comparison` calls it |
| `draw_volcano_plot()` | Volcano plot from `estimate_significance()` output, or one panel per model term from its `by = "term"` reading |
| `draw_grouped_boxplot()` | Boxplots for several features x group levels, or for a crossed design as one panel per feature with the remaining factors along the x axis and the primary factor the boxes, so an interaction is visible inside a panel (`panel_by = "factor"` transposes it) |
| `draw_grouped_barplot()` | One column of `summarize_descriptive_stats()` as clusters of bars, one cluster per feature and one bar per group level; `errorbar` is read under `mainbar`, so a mean takes a standard error, a standard deviation or Student's interval, a median takes the notch the boxplot notches with, and a count or a spread takes none |
| `draw_heatmap()` | Clustered heatmap of features x samples, with the sample groups annotated |
| `draw_corrplot()` | The correlation matrix as a heatmap: nothing standardised, the colours fixed at -1 to 1, one clustering shared by both axes so that the diagonal stays diagonal, and the pairs that did not clear their p-value drawn as blank cells |
| `draw_butterfly_hist()` | Back-to-back histogram, kernel density, or both, for exactly two groups |
| `draw_interaction_plot()` | Cell means of a crossed design joined across one factor, one line per level of another, in one pair of factors, every pair at once, or with a third factor kept in panels of its own |
| `draw_mosaic_plot()` | Mosaic of a contingency table, shaded by the residual of the null the result was tested against and marked where that null would have cut each strip; `plot()` on an `sa_categorical` calls it |
| `split_data()` | Train/test partition, stratified and leakage-aware through `id` |
| `fit_linear_regression()` | Linear model with coefficient inference and resampled performance |
| `fit_logistic_regression()` | Two-class logistic model with odds ratios and their intervals |
| `fit_elastic_net()` | LASSO, ridge or elastic net for either outcome type, with the selected terms |
| `fit_rf()` | Random forest with permutation and impurity importance and out-of-bag fit |
| `fit_svm()` | Radial-kernel support vector machine with permutation importance |
| `evaluate_regression_models()` | Score one or more fitted regressions on the same held-out rows, each against a baseline |
| `evaluate_classification_models()` | The same for a two-class outcome, with DeLong's test, the IDI and the NRI against the baseline |
| `draw_prediction_plot()` | Observed against predicted for a scored regression, with the identity line and the calibration line the metrics table already holds |
| `draw_roc_curve()` | ROC curves of the scored classifications, always overlaid; `plot()` on an `sa_performance` calls whichever of these two the result calls for |
| `perform_rfe()` | Recursive feature elimination inside the resampling, returning the predictors it kept, their ranking and the score at every subset size |
| `perform_stepwise()` | Stepwise search by AIC or BIC, returning the predictors it kept, what each one is worth to the criterion, and the path it walked |
| `perform_pca()` | Principal components of the samples or of the features, with loadings and variance |
| `perform_tsne()` | t-SNE embedding of either margin |
| `perform_umap()` | UMAP embedding of either margin |
| `cluster_hclust()` | Hierarchical clustering of the samples or of the features, returning the tree as well as the cut |
| `cluster_kmeans()` | k-means clustering of either margin, best of 25 starts |
| `cluster_dbscan()` | Density-based clustering, deriving the number of clusters and leaving sparse points as noise |
| `cluster_snn()` | Shared nearest neighbour clustering, which groups by how many neighbours two points have in common rather than by a radius |
| `draw_dim_reduction_plot()` | Two coordinates of a reduction as a scatter, coloured by a clustering of the same points and shaped by a grouping that was known already, so the two can be read against each other; `plot()` on an `sa_reduction` calls it |
| `simulate_two_groups()` | Two-group log2 expression data with the planted answer returned alongside it |
| `simulate_multiple_groups()` | One control and any number of treatment groups, scored per feature, per level and per contrast |
| `simulate_factorial_groups()` | Any number of crossed factors, each between subjects or within them, scored per model term as well as per cell and per contrast |
| `simulate_categorical_groups()` | A contingency table, or repeated binary conditions, with the planted association returned cell by cell, and the symmetric share as well for a matched pair |
| `simulate_regression()` | Continuous outcome from planted coefficients, with a correlation structure and the truth per predictor and per term |
| `simulate_classification()` | Two-class outcome at a chosen event rate from the same design |
| `make_block_cor()` | Block correlation matrix for the simulators, with `against` for the predictors of a block that move the other way, checked for positive definiteness |

---

## Author

**Wonseok Oh** ([ORCID: 0009-0002-0687-8466](https://orcid.org/0009-0002-0687-8466))

## License

MIT © 2026 Wonseok Oh. See [LICENSE.md](LICENSE.md) for details.
