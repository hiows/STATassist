# One-off: PNGs for README.md (run from repo root: Rscript tools/render_readme_figures.R)
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required to load the package for rendering.", call. = FALSE)
}
devtools::load_all(".", quiet = TRUE)

fig_dir <- file.path("man", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

feats4 <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
iris2 <- iris[iris$Species != "setosa", ]
lv2 <- c("versicolor", "virginica")

res <- compare_two_groups(
  data = iris2, feats = feats4, group = iris2$Species, group_lv = lv2
)
sig <- estimate_significance(res, test = "t_test", log2fc_cutoff = 0.25,
                             pval_cutoff = 0.05)

png(file.path(fig_dir, "README-boxplot.png"), width = 1100, height = 550, res = 110)
invisible(draw_grouped_boxplot(
  data = iris, feats = feats4, group = iris$Species,
  group_lv = levels(iris$Species),
  ylab = "cm", main = "Iris measurements by species", dark = TRUE
))
dev.off()

png(file.path(fig_dir, "README-butterfly.png"), width = 800, height = 650, res = 110)
invisible(draw_butterfly_hist(
  data = iris2, feat = "Petal.Length", group = iris2$Species, group_lv = lv2,
  breaks = 12, scale = "proportion",
  main = "Petal length: virginica vs versicolor"
))
dev.off()

png(file.path(fig_dir, "README-volcano.png"), width = 700, height = 650, res = 110)
invisible(draw_volcano_plot(
  sig, main = "Virginica vs versicolor (Welch t-test)"
))
dev.off()

cat("Wrote figures to", fig_dir, "\n")
