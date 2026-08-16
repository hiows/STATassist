# PNGs for README.md. Run from the package directory:
#   Rscript tools/render_readme_figures.R
#
# Every figure the README links to is drawn here, from the same calls the README
# quotes, so a figure cannot disagree with the text beside it. The scenarios and
# their arguments are the ones worked out in Test/readme_blueprint/: one seed
# (2026) for all of them, one correlation structure for every supervised model,
# and one train/test split per outcome type.
#
# The two plots the package does not draw itself, a predicted-against-observed
# scatter and an ROC curve, are base R here rather than new exports. They need
# pROC, and tools/ is in .Rbuildignore, so nothing in this file reaches a user.
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required to load the package for rendering.", call. = FALSE)
}
for (pkg in c("pROC", "colorspace")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " is required to draw the model figures.", call. = FALSE)
  }
}
devtools::load_all(".", quiet = TRUE)

fig_dir <- file.path("man", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# One entry point for every figure, so the size of a plot is stated where the
# plot is drawn and the device is always closed.
figure <- function(name, width, height, expr) {
  path <- file.path(fig_dir, paste0("README-", name, ".png"))
  grDevices::png(path, width = width, height = height, units = "px", res = 110)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  cat("  ", path, "\n", sep = "")
}

# Predicted against observed, with the identity line to read the fitted line
# against. The correlation in the title is the one number the README quotes.
draw_fit_scatter <- function(y, y_hat, main) {
  lim <- range(c(y, y_hat), finite = TRUE)
  plot(y, y_hat, type = "n", xlim = lim, ylim = lim,
       xlab = "y (test data)", ylab = "y_hat",
       main = paste0(main, "\nCorr = ", round(stats::cor(y, y_hat), 3)))
  graphics::abline(a = 0, b = 1, col = "gray60", lwd = 2, lty = 3)
  m <- stats::lm(y_hat ~ y)
  graphics::abline(m, col = "red", lwd = 2)
  graphics::points(y, y_hat, pch = 16)
  intercept <- round(stats::coef(m)[[1]], 3)
  slope <- round(stats::coef(m)[[2]], 3)
  graphics::legend(
    "bottomright",
    paste0("y = ", slope, "x ", if (intercept < 0) "- " else "+ ",
           abs(intercept)),
    col = "red", lwd = 2, bty = "n"
  )
  invisible(m)
}

# One ROC curve per set of predicted probabilities, all against the same
# response, with the AUC in the legend beside the name of the set.
draw_roc <- function(response, probs, main,
                     cols = c("black", "red", "gray70"),
                     ltys = c(3, 1, 1)) {
  aucs <- numeric(length(probs))
  for (i in seq_along(probs)) {
    roc <- pROC::roc(response = response, predictor = probs[[i]],
                     levels = c("control", "case"), quiet = TRUE)
    aucs[i] <- as.numeric(pROC::auc(roc))
    x <- 1 - roc$specificities
    y <- roc$sensitivities
    if (i == 1L) {
      plot(x, y, type = "l", lwd = 2, lty = ltys[i], col = cols[i],
           xlim = c(0, 1), ylim = c(0, 1),
           xlab = "1 - specificity", ylab = "sensitivity", main = main)
      graphics::abline(a = 0, b = 1, col = "gray60", lwd = 2, lty = 3)
    } else {
      graphics::lines(x, y, lwd = 2, lty = ltys[i], col = cols[i])
    }
  }
  graphics::legend(
    "bottomright",
    paste0(names(probs), " (", round(aucs, 3), ")"),
    col = cols[seq_along(probs)], lty = ltys[seq_along(probs)], lwd = 2,
    bty = "n"
  )
  invisible(stats::setNames(round(aucs, 3), names(probs)))
}

# Terms come back from a fitted model with the level pasted onto the column name,
# so a factor has to be named as its column again before it can be refitted.
as_predictors <- function(terms) {
  unique(sub("high$", "", sub("mid$", "", terms)))
}


cat("=== 1. two groups ===\n")

sim2 <- simulate_two_groups(n_feats = 30, n_up = 8, n_down = 8, seed = 2026)
args2 <- sim2$args

comp2 <- compare_two_groups(
  data        = args2$data,
  feats       = args2$feats,
  group       = args2$group,
  group_lv    = args2$group_lv,
  input_scale = args2$input_scale
)

sig2 <- estimate_significance(
  comparison_result = comp2,
  test              = "t_test",
  log2fc_cutoff     = 1,
  pval_cutoff       = 0.05,
  adj_type          = "BH"
)

first_ten <- paste0("gene_", 1:10)

figure("volcano", 700, 650, draw_volcano_plot(sig2, xlim = c(-3, 3)))

figure("boxplot", 1200, 620, draw_grouped_boxplot(
  data     = args2$data,
  feats    = first_ten,
  group    = args2$group,
  group_lv = args2$group_lv,
  ylim     = c(-5, 20)
))

figure("butterfly", 800, 650, draw_butterfly_hist(
  data     = args2$data,
  feat     = "gene_8",
  group    = args2$group,
  group_lv = args2$group_lv,
  breaks   = seq(-5, 20, by = 1),
  type     = "both"
))

figure("forest-pvalue", 900, 620, draw_forest_plot(
  comparison_result = comp2,
  test              = "t_test",
  type              = "pvalue",
  feats             = first_ten,
  sort_by           = "pvalue"
))

figure("forest-estimate", 900, 620, draw_forest_plot(
  comparison_result = comp2,
  test              = "t_test",
  type              = "estimate",
  feats             = first_ten,
  sort_by           = "pvalue",
  xlim              = c(-6, 6)
))

# 100 samples leave no room for 100 column labels, which is the one argument this
# call adds to the blueprint's.
figure("heatmap", 900, 770, draw_heatmap(
  data              = args2$data,
  group             = args2$group,
  group_lv          = args2$group_lv,
  hclust_method     = "ward.D2",
  show_sample_names = FALSE
))


cat("=== 2. three or more groups ===\n")

simN <- simulate_multiple_groups(
  n_feats   = 10,
  n_control = 50,
  n_treat   = c(50, 50, 50),
  n_up      = 3,
  n_down    = 3,
  seed      = 2026,
  paired    = FALSE
)
argsN <- simN$args

compN <- compare_multiple_groups(
  data        = argsN$data,
  feats       = argsN$feats,
  group       = argsN$group,
  group_lv    = argsN$group_lv,
  id          = argsN$id,
  input_scale = argsN$input_scale,
  paired      = FALSE
)

sigN <- estimate_significance(
  comparison_result = compN,
  test              = "anova_test",
  pval_cutoff       = 0.05,
  adj_type          = "BH"
)

figure("multi-volcano", 700, 650, draw_volcano_plot(sigN, xlim = c(-4, 4)))

figure("multi-boxplot", 1200, 620, draw_grouped_boxplot(
  data     = argsN$data,
  feats    = argsN$feats,
  group    = argsN$group,
  group_lv = argsN$group_lv,
  ylim     = c(-10, 20)
))

figure("multi-posthoc", 900, 520, draw_forest_plot(
  comparison_result = compN,
  test              = "anova_test",
  type              = "posthoc",
  feats             = "prot_1",
  sort_by           = "pvalue"
))


cat("=== 3. supervised learning: the shared data ===\n")

# Two predictors correlated at 0.8 and three at 0.5. A null predictor that
# correlates with a planted one is pulled off zero, which is what makes the
# feature selection in the sections below worth scoring.
cor_mat <- make_block_cor(
  n_features = 8,
  blocks = list(
    list(features = 1:2, cor = 0.8),
    list(features = 3:5, cor = 0.5)
  )
)

sim_reg <- simulate_regression(cor_mat = cor_mat, seed = 2026)
reg_args <- sim_reg$args
reg_split <- split_data(
  data       = reg_args$data,
  stratified = reg_args$data$x_cat_1,
  p_train    = 0.75,
  times      = 1,
  seed       = 2026
)
reg_train <- reg_split$datasets[[1]]$train_data
reg_test <- reg_split$datasets[[1]]$test_data

sim_cls <- simulate_classification(cor_mat = cor_mat, seed = 2026)
cls_args <- sim_cls$args
cls_split <- split_data(
  data       = cls_args$data,
  stratified = cls_args$data$y,
  p_train    = 0.75,
  times      = 1,
  seed       = 2026
)
cls_train <- cls_split$datasets[[1]]$train_data
cls_test <- cls_split$datasets[[1]]$test_data

cv <- list(cv = TRUE, cv_method = "repeated_kfold", n_fold = 10, n_repeat = 3,
           seed = 2026)


cat("=== 4. linear and logistic regression ===\n")

lin_all <- do.call(fit_linear_regression, c(list(
  data = reg_train, outcome = reg_args$outcome, predictors = reg_args$predictors
), cv))
lin_keep <- as_predictors(
  coef(lin_all)$terms[-1][coef(lin_all)$pval[-1] < 0.01]
)
lin <- do.call(fit_linear_regression, c(list(
  data = reg_train, outcome = reg_args$outcome, predictors = lin_keep
), cv))

figure("linear-regression", 700, 650, draw_fit_scatter(
  reg_test$y, predict(lin, newdata = reg_test), "Linear regression"
))

log_all <- do.call(fit_logistic_regression, c(list(
  data = cls_train, outcome = cls_args$outcome,
  predictors = cls_args$predictors, outcome_lv = cls_args$outcome_lv
), cv))
log_terms <- coef(log_all)$terms[-1]
log_signif <- as_predictors(log_terms[coef(log_all)$pval[-1] < 0.05])
log_null <- as_predictors(log_terms[coef(log_all)$pval[-1] >= 0.05])

log_signif_fit <- do.call(fit_logistic_regression, c(list(
  data = cls_train, outcome = cls_args$outcome, predictors = log_signif,
  outcome_lv = cls_args$outcome_lv
), cv))
log_null_fit <- do.call(fit_logistic_regression, c(list(
  data = cls_train, outcome = cls_args$outcome, predictors = log_null,
  outcome_lv = cls_args$outcome_lv
), cv))

figure("logistic-regression", 700, 650, draw_roc(
  cls_test$y,
  list(
    "All features"            = predict(log_all, cls_test, type = "response"),
    "Selected features"       = predict(log_signif_fit, cls_test,
                                        type = "response"),
    "Not significant features" = predict(log_null_fit, cls_test,
                                         type = "response")
  ),
  "Logistic regression"
))


cat("=== 5. elastic net ===\n")

lasso <- list(penalty = "lasso", lambda = c(0.01, 0.1, 0.5, 1, 2))

enet_reg_all <- do.call(fit_elastic_net, c(list(
  data = reg_train, outcome = reg_args$outcome,
  predictors = reg_args$predictors
), lasso, cv))
enet_reg_keep <- as_predictors(
  coef(enet_reg_all)$terms[-1][coef(enet_reg_all)$selected[-1]]
)
enet_reg <- do.call(fit_elastic_net, c(list(
  data = reg_train, outcome = reg_args$outcome, predictors = enet_reg_keep
), lasso, cv))

figure("elastic-net-regression", 700, 650, draw_fit_scatter(
  reg_test$y, predict(enet_reg, newdata = reg_test), "Elastic net (LASSO)"
))

enet_cls_all <- do.call(fit_elastic_net, c(list(
  data = cls_train, outcome = cls_args$outcome,
  predictors = cls_args$predictors, outcome_lv = cls_args$outcome_lv
), lasso, cv))
enet_cls_keep <- as_predictors(
  coef(enet_cls_all)$terms[-1][coef(enet_cls_all)$selected[-1]]
)
enet_cls <- do.call(fit_elastic_net, c(list(
  data = cls_train, outcome = cls_args$outcome, predictors = enet_cls_keep,
  outcome_lv = cls_args$outcome_lv
), lasso, cv))
enet_cls_dropped <- setdiff(cls_args$predictors, enet_cls_keep)
enet_cls_null <- do.call(fit_logistic_regression, c(list(
  data = cls_train, outcome = cls_args$outcome,
  predictors = enet_cls_dropped, outcome_lv = cls_args$outcome_lv
), cv))

figure("elastic-net-roc", 700, 650, draw_roc(
  cls_test$y,
  list(
    "All features"     = predict(enet_cls_all, cls_test, type = "response"),
    "Kept features"    = predict(enet_cls, cls_test, type = "response"),
    "Dropped features" = predict(enet_cls_null, cls_test, type = "response")
  ),
  "Elastic net (LASSO)"
))


cat("=== 6. random forest ===\n")

top_five <- function(model) {
  imp <- coef(model)
  as_predictors(imp$terms[order(imp$estimate, decreasing = TRUE)][1:5])
}

rf_reg_all <- do.call(fit_rf, c(list(
  data = reg_train, outcome = reg_args$outcome,
  predictors = reg_args$predictors, mtry = c(2, 5, 8), ntree = 500
), cv))
rf_reg <- do.call(fit_rf, c(list(
  data = reg_train, outcome = reg_args$outcome,
  predictors = top_five(rf_reg_all), ntree = 500
), cv))

figure("random-forest-regression", 700, 650, draw_fit_scatter(
  reg_test$y, predict(rf_reg, newdata = reg_test), "Random forest"
))

rf_cls_all <- do.call(fit_rf, c(list(
  data = cls_train, outcome = cls_args$outcome,
  predictors = cls_args$predictors, outcome_lv = cls_args$outcome_lv,
  mtry = c(2, 5, 8), ntree = 500
), cv))
rf_cls_keep <- top_five(rf_cls_all)
rf_cls <- do.call(fit_rf, c(list(
  data = cls_train, outcome = cls_args$outcome, predictors = rf_cls_keep,
  outcome_lv = cls_args$outcome_lv, ntree = 500
), cv))
rf_cls_rest <- do.call(fit_rf, c(list(
  data = cls_train, outcome = cls_args$outcome,
  predictors = setdiff(cls_args$predictors, rf_cls_keep),
  outcome_lv = cls_args$outcome_lv, ntree = 500
), cv))

figure("random-forest-roc", 700, 650, draw_roc(
  cls_test$y,
  list(
    "All features"            = predict(rf_cls_all, cls_test,
                                        type = "response"),
    "Top importance"          = predict(rf_cls, cls_test, type = "response"),
    "Low importance features" = predict(rf_cls_rest, cls_test,
                                        type = "response")
  ),
  "Random forest"
))


cat("=== 7. support vector machine ===\n")

svm_grid <- list(C = 2^seq(-5, 10, by = 2), sigma = NULL)

svm_reg_all <- do.call(fit_svm, c(list(
  data = reg_train, outcome = reg_args$outcome,
  predictors = reg_args$predictors
), svm_grid, cv))
svm_reg <- do.call(fit_svm, c(list(
  data = reg_train, outcome = reg_args$outcome,
  predictors = top_five(svm_reg_all)
), svm_grid, cv))

figure("svm-regression", 700, 650, draw_fit_scatter(
  reg_test$y, predict(svm_reg, newdata = reg_test), "Support vector machine"
))

svm_cls_all <- do.call(fit_svm, c(list(
  data = cls_train, outcome = cls_args$outcome,
  predictors = cls_args$predictors, outcome_lv = cls_args$outcome_lv
), svm_grid, cv))
svm_cls_keep <- top_five(svm_cls_all)
svm_cls <- do.call(fit_svm, c(list(
  data = cls_train, outcome = cls_args$outcome, predictors = svm_cls_keep,
  outcome_lv = cls_args$outcome_lv
), svm_grid, cv))
svm_cls_rest <- do.call(fit_svm, c(list(
  data = cls_train, outcome = cls_args$outcome,
  predictors = setdiff(cls_args$predictors, svm_cls_keep),
  outcome_lv = cls_args$outcome_lv
), svm_grid, cv))

figure("svm-roc", 700, 650, draw_roc(
  cls_test$y,
  list(
    "All features"            = predict(svm_cls_all, cls_test,
                                        type = "response"),
    "Top importance"          = predict(svm_cls, cls_test, type = "response"),
    "Low importance features" = predict(svm_cls_rest, cls_test,
                                        type = "response")
  ),
  "Support vector machine"
))


cat("=== 8. dimension reduction ===\n")

# A third block is added here so that the feature embedding has three groups of
# correlated predictors to separate and one predictor belonging to none.
red_cor <- make_block_cor(
  n_features = 8,
  blocks = list(
    list(features = 1:2, cor = 0.8),
    list(features = 3:5, cor = 0.5),
    list(features = 7:8, cor = 0.9)
  )
)
red_data <- simulate_classification(cor_mat = red_cor, seed = 2026)$args$data
red_feats <- paste0("x_", 1:8)
# One colour per block, in the order make_block_cor() laid them out: x_1-x_2,
# x_3-x_5, x_6 alone, x_7-x_8.
red_cols <- rep(colorspace::qualitative_hcl(4, "Dark 2"), times = c(2, 3, 1, 2))

label_plot <- function(x, y, labels, xlab, ylab, main) {
  plot(x, y, type = "n", xlab = xlab, ylab = ylab, main = main)
  graphics::text(x, y, labels = labels, col = red_cols, font = 2)
}

pca <- perform_pca(
  data = red_data, feats = red_feats, embedding_scale = "features",
  center = TRUE, scale = TRUE
)
figure("pca", 700, 650, label_plot(
  pca$scores$PC1, pca$scores$PC2, pca$scores$points,
  paste0("PC1 (", round(pca$variance$prop_var[1], 2), "%)"),
  paste0("PC2 (", round(pca$variance$prop_var[2], 2), "%)"),
  "PCA of the features"
))

tsne <- perform_tsne(
  data = red_data, feats = red_feats, embedding_scale = "features",
  center = TRUE, scale = TRUE, seed = 2026
)
figure("tsne", 700, 650, label_plot(
  tsne$scores$tSNE1, tsne$scores$tSNE2, tsne$scores$points,
  "t-SNE1", "t-SNE2", "t-SNE of the features"
))

umap_res <- perform_umap(
  data = red_data, feats = red_feats, embedding_scale = "features",
  n_neighbors = 3, center = FALSE, scale = FALSE, seed = 2026
)
figure("umap", 700, 650, label_plot(
  umap_res$scores$UMAP1, umap_res$scores$UMAP2, umap_res$scores$points,
  "UMAP1", "UMAP2", "UMAP of the features"
))


cat("=== 9. factorial crossed design ===\n")

sim_fact <- simulate_factorial_groups(seed = 2026)
fact_args <- sim_fact$args
fact_comp <- compare_factorial_groups(
  data          = fact_args$data,
  feats         = fact_args$feats,
  factors       = fact_args$factors,
  factor_lv     = fact_args$factor_lv,
  control_label = list(treatment = "control", sex = "male"),
  input_scale   = fact_args$input_scale
)
sig_fact <- estimate_significance(fact_comp)
sig_fact_term <- estimate_significance(fact_comp, by = "term")
interest_feats <- paste0("prot_", 1:20)
interest_feat <- "prot_14"

figure("factorial-forest-pvalue", 900, 520, draw_forest_plot(
  comparison_result = fact_comp,
  type              = "pvalue",
  feats             = interest_feats,
  sort_by           = "pvalue"
))

figure("factorial-volcano", 900, 650, draw_volcano_plot(sig_fact_term))

figure("factorial-forest-estimate", 900, 520, draw_forest_plot(
  comparison_result = fact_comp,
  feats             = interest_feat,
  sort_by           = "pvalue",
  xlim              = c(-6, 6)
))

figure("factorial-interaction", 800, 620, draw_interaction_plot(
  comparison_result = fact_comp,
  feats             = interest_feat
))

figure("factorial-boxplot", 800, 620, draw_grouped_boxplot(
  data          = fact_args$data,
  feats         = interest_feat,
  factors       = fact_args$factors,
  factor_lv     = fact_args$factor_lv,
  control_label = list(treatment = "control", sex = "male"),
  ylim          = c(5, 25)
))


cat("=== 10. categorical contingency table ===\n")

sim_cat <- simulate_categorical_groups(seed = 2026)
cat_args <- sim_cat$args
cat_comp <- compare_categorical_groups(
  data          = cat_args$data,
  category_lv   = cat_args$category_lv,
  control_label = list(cat_1 = "n", cat_2 = "mid"),
  paired        = cat_args$paired
)

figure("mosaic", 800, 620, draw_mosaic_plot(cat_comp))


cat("=== 11. grouped barplot ===\n")

sim_bar <- simulate_two_groups(
  n_feats = 10, n_up = 3, n_down = 3, seed = 2026
)
bar_args <- sim_bar$args

figure("grouped-barplot", 900, 620, draw_grouped_barplot(
  data          = bar_args$data,
  feats         = bar_args$feats,
  group         = bar_args$group,
  group_lv      = bar_args$group_lv,
  control_label = "control",
  errorbar      = "se"
))


cat("=== 12. feature-pair association ===\n")

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
assoc_data <- assoc_sim$args$data[, -1, drop = FALSE]
assoc_feats <- colnames(assoc_data)
assoc <- summarize_association_stats(data = assoc_data, feats = assoc_feats)

figure("corrplot", 700, 700, draw_corrplot(
  assoc$pearson$corr, cex.axis = 0.9
))

figure("corrplot-masked", 700, 700, draw_corrplot(
  assoc$pearson$corr,
  pvalue    = assoc$pearson$adj_pvalue,
  cex.axis  = 0.9
))


cat("=== 13. evaluate regression models ===\n")

eval_reg_sim <- simulate_regression(cor_mat = cor_mat, seed = 2026)
eval_reg_args <- eval_reg_sim$args
eval_reg_split <- split_data(
  data = eval_reg_args$data, p_train = 0.75, times = 1, seed = 2026
)
eval_reg_train <- eval_reg_split$datasets[[1]]$train_data
eval_reg_test <- eval_reg_split$datasets[[1]]$test_data

eval_rfe <- perform_rfe(
  data       = eval_reg_train,
  outcome    = eval_reg_args$outcome,
  predictors = eval_reg_args$predictors,
  seed       = 2026
)
eval_sel <- eval_rfe$selected

eval_lin <- do.call(fit_linear_regression, c(list(
  data = eval_reg_train, outcome = eval_reg_args$outcome,
  predictors = eval_sel
), cv))
eval_lasso <- do.call(fit_elastic_net, c(list(
  data = eval_reg_train, outcome = eval_reg_args$outcome,
  predictors = eval_sel, penalty = "lasso"
), cv))
eval_rf <- do.call(fit_rf, c(list(
  data = eval_reg_train, outcome = eval_reg_args$outcome,
  predictors = eval_sel
), cv))
eval_svm <- do.call(fit_svm, c(list(
  data = eval_reg_train, outcome = eval_reg_args$outcome,
  predictors = eval_sel
), svm_grid, cv))

eval_reg <- evaluate_regression_models(
  baseline_model = eval_lin,
  new_models     = list(lasso = eval_lasso, rf = eval_rf, svm = eval_svm),
  newdata        = eval_reg_test,
  answer         = eval_reg_test$y,
  baseline_label = "linear"
)

figure("eval-regression", 800, 650, draw_prediction_plot(
  performance_result = eval_reg,
  type               = "overlay",
  anno_corr          = TRUE,
  anno_rsq           = TRUE,
  anno_lm            = TRUE,
  cex.anno           = 0.85
))


cat("=== 14. evaluate classification models ===\n")

eval_cls_sim <- simulate_classification(cor_mat = cor_mat, seed = 2026)
eval_cls_args <- eval_cls_sim$args
eval_cls_split <- split_data(
  data = eval_cls_args$data, stratified = eval_cls_args$data$y,
  p_train = 0.75, times = 1, seed = 2026
)
eval_cls_train <- eval_cls_split$datasets[[1]]$train_data
eval_cls_test <- eval_cls_split$datasets[[1]]$test_data

eval_rfe_cls <- perform_rfe(
  data          = eval_cls_train,
  outcome       = eval_cls_args$outcome,
  predictors    = eval_cls_args$predictors,
  outcome_lv    = eval_cls_args$outcome_lv,
  control_label = "control",
  seed          = 2026,
  model         = "logistic"
)
eval_sel_cls <- eval_rfe_cls$selected

eval_log <- do.call(fit_logistic_regression, c(list(
  data = eval_cls_train, outcome = eval_cls_args$outcome,
  predictors = eval_sel_cls, outcome_lv = eval_cls_args$outcome_lv,
  control_label = "control"
), cv))
eval_lasso_cls <- do.call(fit_elastic_net, c(list(
  data = eval_cls_train, outcome = eval_cls_args$outcome,
  predictors = eval_sel_cls, outcome_lv = eval_cls_args$outcome_lv,
  penalty = "lasso"
), cv))
eval_rf_cls <- do.call(fit_rf, c(list(
  data = eval_cls_train, outcome = eval_cls_args$outcome,
  predictors = eval_sel_cls, outcome_lv = eval_cls_args$outcome_lv
), cv))
eval_svm_cls <- do.call(fit_svm, c(list(
  data = eval_cls_train, outcome = eval_cls_args$outcome,
  predictors = eval_sel_cls, outcome_lv = eval_cls_args$outcome_lv
), svm_grid, cv))

eval_cls <- evaluate_classification_models(
  baseline_model = eval_log,
  new_models     = list(
    lasso = eval_lasso_cls, rf = eval_rf_cls, svm = eval_svm_cls
  ),
  newdata        = eval_cls_test,
  answer         = eval_cls_test$y,
  outcome_lv     = eval_cls_args$outcome_lv,
  control_label  = "control",
  baseline_label = "logistic"
)

figure("eval-classification", 650, 650, draw_roc_curve(
  performance_result = eval_cls,
  anno_auc           = TRUE,
  cex.anno           = 1
))


cat("=== 15. cluster on an embedding ===\n")

clust_sim <- simulate_two_groups(
  n_feats = 50, deg_log2fc = c(5, 10), seed = 2026
)
clust_args <- clust_sim$args

clust_pca <- perform_pca(
  data = clust_args$data,
  feats = clust_args$feats,
  embedding_scale = "samples"
)
clust_km_pca <- cluster_kmeans(
  data = clust_pca$scores,
  feats = c("PC1", "PC2"),
  cluster_scale = "samples",
  n_clust = 2,
  seed = 2026
)

figure("cluster-pca-group", 700, 650, draw_dim_reduction_plot(
  reduction_result = clust_pca,
  group            = clust_args$group,
  group_lv         = clust_args$group_lv,
  col              = c("black", "red3")
))

figure("cluster-pca-cluster", 700, 650, draw_dim_reduction_plot(
  reduction_result = clust_pca,
  cluster_result   = clust_km_pca,
  cluster_lv       = c("Cluster1", "Cluster2")
))

clust_umap <- perform_umap(
  data = clust_args$data,
  feats = clust_args$feats,
  embedding_scale = "samples",
  seed = 2026
)
clust_km_umap <- cluster_kmeans(
  data = clust_umap$scores,
  feats = c("UMAP1", "UMAP2"),
  cluster_scale = "samples",
  n_clust = 2,
  seed = 2026
)

figure("cluster-umap-group", 700, 650, draw_dim_reduction_plot(
  reduction_result = clust_umap,
  group            = clust_args$group,
  group_lv         = clust_args$group_lv,
  col              = c("black", "red3")
))

figure("cluster-umap-cluster", 700, 650, draw_dim_reduction_plot(
  reduction_result = clust_umap,
  cluster_result   = clust_km_umap,
  cluster_lv       = c("Cluster1", "Cluster2")
))


cat("\nWrote", length(list.files(fig_dir, pattern = "\\.png$")),
    "figures to", fig_dir, "\n")
