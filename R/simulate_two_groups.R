#' Simulate a two-group experiment whose answer is known
#'
#' Generates log2-scale expression data for two independent groups with a fixed
#' number of features moved up and down on purpose, and returns the planted
#' answer alongside the data. Every quantity a comparison estimates can then be
#' checked against what was actually put there, which is what a real data set
#' can never offer.
#'
#' The point of the exercise is the gap. A comparison does not recover every
#' feature that was planted, and the reasons it misses them are the three things
#' worth understanding about a volcano plot: the p-value may not clear its
#' cutoff, the multiplicity adjustment may take it back, and the estimated
#' `log2fc` carries a sampling error of its own, so a feature planted just above
#' the magnitude cutoff lands below it about half the time. The defaults are set
#' so that a run recovers most but not all of what was planted, because a
#' simulation that recovers everything teaches none of this.
#'
#' @param n_feats Number of features to generate. Columns are named `gene_1`
#'   upwards.
#' @param n_case,n_control Observations in each group. They do not have to
#'   match.
#' @param n_up,n_down How many features are moved up and down in the case group.
#'   Their sum cannot exceed `n_feats`, and every other feature is left with a
#'   true fold change of exactly zero.
#' @param expr_range Range the baseline log2 expression of each feature is drawn
#'   from. The default spans what log2 CPM or RMA values usually cover. Both
#'   groups share the baseline, which is what makes an unplanted feature null.
#' @param case_sd,control_sd Ranges the per-feature standard deviation of each
#'   group is drawn from. They are drawn independently, so the groups end up
#'   with unequal variances, which is the situation Welch's t-test and the
#'   Brunner-Munzel test exist for. Pass the same range twice for a
#'   homoscedastic data set. The defaults leave roughly four planted features in
#'   five recoverable at the default cutoffs; narrowing them recovers nearly
#'   everything and widening them costs recall quickly.
#' @param deg_log2fc Range the magnitude of the planted effect is drawn from, on
#'   the log2 scale. The default of `c(1, 2.5)` is a two-fold to roughly
#'   six-fold change, which straddles the `log2fc_cutoff = 1` that
#'   [estimate_significance()] applies by default.
#' @param group_lv The two group labels, the first being the control and the
#'   second the one the effect is applied to. Passed straight through to the
#'   returned arguments, so it also fixes the direction
#'   [compare_two_groups()] reads: a planted increase comes back as a positive
#'   `log2fc` because the control is the reference.
#' @param seed Seed for the draw, or `NULL` to use the stream as it stands.
#'   Supplying one does not disturb the caller: the previous random number state
#'   is put back when the function returns.
#'
#' @return A list of two elements.
#'
#'   \describe{
#'     \item{`args`}{`data`, `feats`, `group`, `group_lv` and `input_scale`,
#'       named after the arguments of [compare_two_groups()] so that
#'       `do.call(compare_two_groups, sim$args)` runs the comparison.
#'       `input_scale` is `"log2"`, since that is the scale the data is on.}
#'     \item{`truth`}{One row per feature, aligned with `feats`, holding
#'       `features`, `direction` (`"up"`, `"down"` or `"none"`), `log2fc` (the
#'       effect that was planted, exactly `0` for `"none"`), `baseline` and the
#'       two group standard deviations. The last three are there so that a
#'       feature the comparison missed can be looked up rather than guessed at:
#'       a large `sd_case` explains a miss that the effect size alone does not.}
#'   }
#'
#' @section How the data is built:
#' Each feature gets a baseline `b` drawn from `expr_range`, two standard
#' deviations drawn from `case_sd` and `control_sd`, and a planted effect `d`
#' that is a positive draw from `deg_log2fc` when the feature is one of the
#' `n_up`, a negative draw when it is one of the `n_down`, and `0` otherwise.
#' Case observations are then normal around `b + d` and control observations
#' normal around `b`.
#'
#' Because the baseline is shared, the true log2 fold change of a feature is `d`
#' and nothing else. An unplanted feature is null in the strict sense, so a
#' feature called significant is a false positive by definition and the
#' multiplicity adjustment can be judged on it. A model that gave each group its
#' own random offset would look more lifelike and would make both the recall and
#' the false positive rate impossible to compute, since an unplanted feature
#' would then differ between the groups too.
#'
#' The data is on the log2 scale throughout, so `input_scale = "log2"` comes
#' back with it. The effect is added rather than multiplied, which is what makes
#' `deg_log2fc` a difference of log2 means rather than a ratio.
#'
#' @seealso [compare_two_groups()], which consumes `args` directly, and
#'   [estimate_significance()] for the verdict that `truth` is there to score.
#'
#' @examples
#' sim <- simulate_two_groups(seed = 1)
#' table(sim$truth$direction)
#'
#' ## The names in `args` are compare_two_groups()'s own, so the comparison is
#' ## one call away.
#' res <- do.call(compare_two_groups, sim$args)
#' sig <- estimate_significance(res, test = "t_test")$significance
#'
#' ## Scored against what was planted. The off-diagonal cells are the two kinds
#' ## of mistake: features that were planted and missed, and null features that
#' ## were called anyway.
#' planted <- sim$truth$direction != "none"
#' called <- sig$is_signif %in% TRUE
#' table(planted = planted, called = called)
#'
#' ## The direction is recovered too, not just the fact of a difference.
#' table(truth = sim$truth$direction[called], sign = sign(sig$log2fc[called]))
#'
#' ## Recall differs between the three families on the same data and the same
#' ## truth, which is the reason all three are reported.
#' vapply(names(res$tests), function(nm) {
#'   hit <- estimate_significance(res, test = nm)$significance$is_signif
#'   mean(hit[planted] %in% TRUE)
#' }, numeric(1))
#'
#' ## Turning the noise up costs recall without changing what was planted.
#' noisy <- simulate_two_groups(seed = 1, case_sd = c(4, 6),
#'                              control_sd = c(4, 6))
#' noisy_sig <- estimate_significance(
#'   do.call(compare_two_groups, noisy$args)
#' )$significance
#' mean(noisy_sig$is_signif[noisy$truth$direction != "none"] %in% TRUE)
#'
#' @export
simulate_two_groups <- function(n_feats = 100,
                                n_case = 50,
                                n_control = 50,
                                n_up = 15,
                                n_down = 15,
                                expr_range = c(2, 12),
                                case_sd = c(1.8, 3.2),
                                control_sd = c(1.2, 2.4),
                                deg_log2fc = c(1, 2.5),
                                group_lv = c("control", "case"),
                                seed = NULL) {

  n_feats <- sa_check_count(n_feats, "n_feats", 1)
  n_case <- sa_check_count(n_case, "n_case", 2)
  n_control <- sa_check_count(n_control, "n_control", 2)
  n_up <- sa_check_count(n_up, "n_up")
  n_down <- sa_check_count(n_down, "n_down")
  if (n_up + n_down > n_feats) {
    stop("`n_up` + `n_down` is ", n_up + n_down, ", which is more features ",
         "than the ", n_feats, " that `n_feats` asks for.", call. = FALSE)
  }
  sa_check_range(expr_range, "expr_range")
  sa_check_range(case_sd, "case_sd", 0)
  sa_check_range(control_sd, "control_sd", 0)
  sa_check_range(deg_log2fc, "deg_log2fc", 0)
  if (!is.character(group_lv) || length(group_lv) != 2L || anyNA(group_lv) ||
      group_lv[1] == group_lv[2]) {
    stop("`group_lv` must be two distinct non-missing group labels.",
         call. = FALSE)
  }

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  feats <- paste0("gene_", seq_len(n_feats))
  baseline <- stats::runif(n_feats, expr_range[1], expr_range[2])
  sd_case <- stats::runif(n_feats, case_sd[1], case_sd[2])
  sd_control <- stats::runif(n_feats, control_sd[1], control_sd[2])

  direction <- rep("none", n_feats)
  delta <- numeric(n_feats)
  if (n_up + n_down > 0L) {
    # Taken from the head and the tail of one shuffled draw. Selecting the down
    # set as the complement of the up set with `-idx` would return everything
    # rather than nothing when the up set is empty.
    picked <- sample.int(n_feats, n_up + n_down)
    up_idx <- utils::head(picked, n_up)
    down_idx <- utils::tail(picked, n_down)
    direction[up_idx] <- "up"
    direction[down_idx] <- "down"
    delta[up_idx] <- stats::runif(n_up, deg_log2fc[1], deg_log2fc[2])
    delta[down_idx] <- -stats::runif(n_down, deg_log2fc[1], deg_log2fc[2])
  }

  draw <- function(n, center, spread) {
    vapply(seq_len(n_feats), function(i) {
      stats::rnorm(n, mean = center[i], sd = spread[i])
    }, numeric(n))
  }
  # Drawn case first and stacked control first. `group_lv` names the control
  # ahead of the case group, so the rows have to follow, while the order the
  # draws consume the random stream is left alone: reversing that instead would
  # hand a seed that used to give one data set a different one.
  case_values <- draw(n_case, baseline + delta, sd_case)
  control_values <- draw(n_control, baseline, sd_control)
  values <- rbind(control_values, case_values)

  data <- as.data.frame(values)
  names(data) <- feats
  rownames(data) <- NULL

  list(
    args = list(
      data        = data,
      feats       = feats,
      group       = rep(group_lv, times = c(n_control, n_case)),
      group_lv    = group_lv,
      input_scale = "log2"
    ),
    truth = data.frame(
      features   = feats,
      direction  = direction,
      log2fc     = delta,
      baseline   = baseline,
      sd_case    = sd_case,
      sd_control = sd_control,
      stringsAsFactors = FALSE
    )
  )
}
