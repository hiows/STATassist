#' Simulate a control-versus-treatments experiment whose answer is known
#'
#' The multi-group counterpart of [simulate_two_groups()]. Generates log2-scale
#' abundance data for one control group and any number of treatment groups, and
#' returns the planted answer alongside the data so that a comparison can be
#' scored against what was actually put there.
#'
#' With two groups there is one thing to get right: whether a feature moved.
#' With three or more there are two, and they fail separately. The omnibus test
#' asks whether the levels are all alike, and the post-hoc stage asks which of
#' them differ. A feature can clear the first and be misread by the second, and
#' a shape of effect that the omnibus finds easy can be the one the post-hoc
#' stage finds hard. That is why the effect is planted in three shapes rather
#' than one, and why the answer comes back in three tables rather than one.
#'
#' @param n_feats Number of features to generate. Columns are named
#'   `prot_1` upwards, or whatever `feat_prefix` asks for.
#' @param n_control Number of observations in the control group.
#' @param n_treat Observations in each treatment group, one entry per group, so
#'   that its length is how many treatment groups there are. Pass
#'   `c(50, 40, 30)` for three treatment groups of different sizes and
#'   `rep(30, 5)` for five of the same size. There is no separate argument for
#'   the number of groups: the sizes already say it, and two arguments that
#'   could disagree would have to be settled somewhere.
#'
#'   [compare_multiple_groups()] needs at least three levels in all, so this
#'   needs at least two entries. A single number is allowed only when
#'   `group_lv` says how many groups to spread it over.
#' @param n_up,n_down How many features are moved up and down in the treatment
#'   groups. Their sum cannot exceed `n_feats`, and every other feature is left
#'   with a true effect of exactly zero in every group. The defaults are a
#'   fraction of `n_feats` rather than a fixed count, so that asking for fewer
#'   features plants fewer effects instead of failing; at the default
#'   `n_feats = 100` they are the 15 and 15 the defaults were tuned at.
#' @param pattern_mix Named vector of relative weights over the three shapes an
#'   effect can take, described under "The three shapes" below. Set a weight to
#'   zero to leave that shape out. The planted features are split between the
#'   shapes in these proportions by the largest remainder method rather than
#'   drawn at random, so the counts are exactly what the weights ask for and do
#'   not move with the seed.
#' @param expr_range Range the baseline log2 abundance of each feature is drawn
#'   from. Every level shares the baseline, which is what makes an unplanted
#'   feature null.
#' @param control_sd,treat_sd Ranges the per-feature standard deviation of the
#'   control group and of each treatment group are drawn from. Every group draws
#'   its own, so the design is heteroscedastic, which is the situation Welch's
#'   ANOVA and Games-Howell exist for. Pass the same range twice for equal
#'   variances. The defaults are the two ranges [simulate_two_groups()] uses and
#'   leave roughly four planted features in five recoverable by the omnibus test
#'   at the default cutoffs; narrowing them recovers nearly everything and
#'   widening them costs recall quickly.
#' @param deg_log2fc Range the magnitude of the planted effect is drawn from, on
#'   the log2 scale. This is the magnitude at the level that carries the full
#'   effect; the `"gradient"` shape places the intermediate levels below it.
#' @param paired Logical. If `TRUE`, the levels are treated as repeated
#'   conditions measured on the same subjects, and `args` gains `id` and
#'   `paired` so that [compare_multiple_groups()] runs the within-subject tests.
#' @param subject_sd Range the per-feature subject standard deviation is drawn
#'   from. A subject's offset is drawn once per feature and reused across every
#'   condition, which is what a within-subject test exists to remove. The
#'   default is deliberately of the same order as the residual spread, so that
#'   analysing the same table without `id` costs most of the recall. A small
#'   offset would let the unpaired analysis do nearly as well and would teach
#'   the opposite of why the design exists. Ignored when `paired = FALSE`.
#' @param group_lv Group labels, the first being the control that every effect
#'   is planted against. Defaults to `c("control", "treat_1", ...)`, one label
#'   per entry of `n_treat` plus the control. When supplied it says how many
#'   groups there are, so a `n_treat` of length one is spread over them and a
#'   `n_treat` left at its default is replaced by that many groups of the
#'   default size. Supplying labels and sizes that count differently is an
#'   error rather than a guess.
#' @param feat_prefix Prefix for the generated feature names. `"prot"` gives
#'   `prot_1`, `prot_2` and so on.
#' @param seed Seed for the draw, or `NULL` to use the stream as it stands.
#'   Supplying one does not disturb the caller: the previous random number state
#'   is put back when the function returns.
#'
#' @return A list of four elements.
#'
#'   \describe{
#'     \item{`args`}{`data`, `feats`, `group`, `group_lv` and `input_scale`,
#'       named after the arguments of [compare_multiple_groups()] so that
#'       `do.call(compare_multiple_groups, sim$args)` runs the comparison. Under
#'       `paired = TRUE` it also carries `id` and `paired`.}
#'     \item{`truth`}{One row per feature, aligned with `feats`, holding
#'       `features`, `pattern`, `direction`, `extreme_level`, `extreme_tied`,
#'       `log2fc`, `baseline` and `sd_subject`. This is the table that scores
#'       `$effect` and the omnibus tests.}
#'     \item{`truth_group`}{One row per feature and level, holding `features`,
#'       `group`, `is_ref`, `delta`, `center`, `sd` and `n`. A feature the
#'       comparison missed can be looked up here rather than guessed at: a large
#'       `sd` explains a miss that the effect size alone does not.}
#'     \item{`truth_contrast`}{One row per feature and pair of levels, in the
#'       row order and direction the post-hoc tables use, holding `features`,
#'       `contrast`, `group1`, `group2`, `delta` and `is_diff`. This is the
#'       table that scores `$posthoc`.}
#'   }
#'
#' @section The three shapes:
#' Each planted feature is given a magnitude `d` drawn from `deg_log2fc`,
#' positive for an up feature and negative for a down one, and one of three
#' shapes that decides what each treatment group does with it.
#'
#' \describe{
#'   \item{`"all"`}{Every treatment group is shifted by `d`. Only the control
#'     stands apart, so the omnibus test has the whole effect to work with and
#'     every contrast against the control should be found.}
#'   \item{`"gradient"`}{Treatment group `g` of `k` is shifted by `d * g / k`,
#'     so the last one carries the full effect and the ones before it carry a
#'     fraction. This is the dose-response shape, and its early contrasts are
#'     the ones a post-hoc stage loses first.}
#'   \item{`"single"`}{One treatment group, chosen at random, is shifted by `d`
#'     and the rest are left at exactly zero. The omnibus test is diluted here,
#'     since most of the levels it compares are alike, so this is the shape it
#'     misses most often. When it does clear the cutoff, exactly the contrasts
#'     involving that one level should come back.}
#' }
#'
#' A feature that was not planted has a delta of exactly zero in every group.
#' Both kinds of mistake are therefore defined: a contrast called significant on
#' a zero delta is a false positive, and a non-zero delta that was not called is
#' a miss.
#'
#' @section Directions:
#' `truth$log2fc` is the delta of whichever level sits furthest from the control,
#' which is the quantity `$effect$log2fc` estimates. A treatment group that went
#' up gives a positive value in both.
#'
#' `truth_contrast$delta` reads the same way, because the post-hoc tables do. A
#' post-hoc `estimate` is `group1 - group2` with `group1` the later level of
#' `group_lv`, and the control is the first level, so the contrast is
#' `treat_1 - control` and a feature whose treatment groups went up is positive
#' there too. `truth_contrast$delta` is computed as the same difference, so
#' every direction in the three tables points one way.
#'
#' Under the `"all"` shape every treatment group carries the same delta, so no
#' single level is furthest from the control. `extreme_level` then records the
#' first of the tied levels and `extreme_tied` is `TRUE`, which is the flag that
#' says to score the magnitude rather than the name of the level. It is also
#' `TRUE`, with `extreme_level` missing, for an unplanted feature.
#'
#' @section Repeated conditions:
#' Under `paired = TRUE` each subject is measured under every condition, so no
#' subject is dropped and `design$unmatched_ids` comes back empty. Each subject
#' gets an offset per feature, drawn once and added to all of its conditions,
#' which is the between-subject variation the within-subject tests remove. The
#' residual standard deviation still differs between conditions, so sphericity
#' does not hold and the Mauchly, Greenhouse-Geisser and Huynh-Feldt columns of
#' the repeated measures ANOVA have something to report.
#'
#' The same subjects appearing under every condition also means every group
#' holds the same number of them, so `n_control` and every entry of `n_treat`
#' have to agree. Unequal sizes are rejected rather than quietly levelled,
#' since the sizes are the clearest statement of which design was meant.
#'
#' @seealso [compare_multiple_groups()], which consumes `args` directly,
#'   [simulate_two_groups()] for the two-group case, and
#'   [estimate_significance()] for the verdict that `truth` is there to score.
#'
#' @examples
#' sim <- simulate_multiple_groups(n_feats = 30, n_up = 5, n_down = 5, seed = 1)
#' table(pattern = sim$truth$pattern, direction = sim$truth$direction)
#'
#' ## The names in `args` are compare_multiple_groups()'s own, so the comparison
#' ## is one call away.
#' res <- do.call(compare_multiple_groups, c(sim$args, diagnose = FALSE))
#' sig <- estimate_significance(res, test = "anova_test")$significance
#'
#' ## Scored against what was planted. The off-diagonal cells are the two kinds
#' ## of mistake: features that were planted and missed, and null features that
#' ## were called anyway.
#' planted <- sim$truth$direction != "none"
#' table(planted = planted, called = sig$is_signif %in% TRUE)
#'
#' ## The shape of the effect decides how hard the omnibus test finds it. A
#' ## "single" feature differs from the control in one level out of several, so
#' ## most of what the test compares is alike.
#' tapply(sig$is_signif %in% TRUE, sim$truth$pattern, mean)
#'
#' ## The pairwise stage is scored against `truth_contrast`, which is already in
#' ## the row order and the direction the post-hoc table uses.
#' ph <- merge(res$posthoc$anova_test, sim$truth_contrast,
#'             by = c("features", "contrast"))
#' table(differs = ph$is_diff, called = ph$pval_adj <= 0.05)
#'
#' ## A missed feature is looked up rather than guessed at. The row per level
#' ## carries the size the group was given, which is one of the reasons.
#' subset(sim$truth_group, features == sim$truth$features[1])
#'
#' ## `n_treat` is one size per treatment group, so its length is how many there
#' ## are. Five groups of unequal size, against a control of 40:
#' uneven <- simulate_multiple_groups(n_feats = 20, n_control = 40,
#'                                    n_treat = c(30, 25, 20, 15, 10), seed = 1)
#' table(uneven$args$group)[uneven$args$group_lv]
#'
#' ## Labels alone also say how many groups there are, so one size is spread
#' ## over them.
#' dose <- simulate_multiple_groups(n_feats = 20, n_treat = 25,
#'                                  group_lv = c("dmso", "low", "mid", "high"),
#'                                  seed = 1)
#' table(dose$args$group)[dose$args$group_lv]
#'
#' ## Repeated conditions: the same subjects under every treatment, so every
#' ## group holds the same number of them.
#' rep_sim <- simulate_multiple_groups(n_feats = 10, n_up = 2, n_down = 2,
#'                                     n_control = 12, n_treat = rep(12, 3),
#'                                     paired = TRUE, seed = 1)
#' rep_res <- do.call(compare_multiple_groups, c(rep_sim$args, diagnose = FALSE))
#' rep_res$tests$anova_test[1:3, c("features", "f_stat", "pval_adj")]
#'
#' @export
simulate_multiple_groups <- function(n_feats = 100,
                                     n_control = 50,
                                     n_treat = rep(50, 3),
                                     n_up = round(0.15 * n_feats),
                                     n_down = round(0.15 * n_feats),
                                     pattern_mix = c(all = 1, gradient = 1,
                                                     single = 1),
                                     expr_range = c(2, 12),
                                     control_sd = c(1.2, 2.4),
                                     treat_sd = c(1.8, 3.2),
                                     deg_log2fc = c(1, 2.5),
                                     paired = FALSE,
                                     subject_sd = c(2, 4),
                                     group_lv = NULL,
                                     feat_prefix = "prot",
                                     seed = NULL) {

  sa_check_flag(paired, "paired")
  design <- sa_sim_design(n_control, n_treat, group_lv, missing(n_treat),
                          paired)
  group_lv <- design$group_lv
  sizes <- design$sizes
  n_lv <- length(group_lv)
  n_treat_groups <- n_lv - 1L

  # `n_up` and `n_down` default to a fraction of `n_feats`, so their promises
  # are forced after this line and see the checked value rather than whatever
  # was passed in.
  n_feats <- sa_check_count(n_feats, "n_feats", 1)
  n_up <- sa_check_count(n_up, "n_up")
  n_down <- sa_check_count(n_down, "n_down")
  if (n_up + n_down > n_feats) {
    stop("`n_up` + `n_down` is ", n_up + n_down, ", which is more features ",
         "than the ", n_feats, " that `n_feats` asks for.", call. = FALSE)
  }
  mix <- sa_sim_pattern_mix(pattern_mix)
  sa_check_range(expr_range, "expr_range")
  sa_check_range(control_sd, "control_sd", 0)
  sa_check_range(treat_sd, "treat_sd", 0)
  sa_check_range(deg_log2fc, "deg_log2fc", 0)
  sa_check_range(subject_sd, "subject_sd", 0)
  if (!is.character(feat_prefix) || length(feat_prefix) != 1L ||
      is.na(feat_prefix) || !nzchar(feat_prefix)) {
    stop("`feat_prefix` must be a single non-empty string.", call. = FALSE)
  }

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  feats <- paste0(feat_prefix, "_", seq_len(n_feats))
  baseline <- stats::runif(n_feats, expr_range[1], expr_range[2])

  sd_mat <- matrix(0, nrow = n_feats, ncol = n_lv,
                   dimnames = list(feats, group_lv))
  sd_mat[, 1L] <- stats::runif(n_feats, control_sd[1], control_sd[2])
  for (g in seq_len(n_treat_groups)) {
    sd_mat[, 1L + g] <- stats::runif(n_feats, treat_sd[1], treat_sd[2])
  }
  sd_subject <- if (paired) {
    stats::runif(n_feats, subject_sd[1], subject_sd[2])
  } else {
    rep(NA_real_, n_feats)
  }

  delta <- matrix(0, nrow = n_feats, ncol = n_lv,
                  dimnames = list(feats, group_lv))
  direction <- rep("none", n_feats)
  pattern <- rep("none", n_feats)

  if (n_up + n_down > 0L) {
    # Taken from the head and the tail of one shuffled draw. Selecting the down
    # set as the complement of the up set with `-idx` would return everything
    # rather than nothing when the up set is empty.
    picked <- sample.int(n_feats, n_up + n_down)
    up_idx <- utils::head(picked, n_up)
    down_idx <- utils::tail(picked, n_down)
    direction[up_idx] <- "up"
    direction[down_idx] <- "down"

    plant_idx <- c(up_idx, down_idx)
    plant_mag <- c(stats::runif(n_up, deg_log2fc[1], deg_log2fc[2]),
                   -stats::runif(n_down, deg_log2fc[1], deg_log2fc[2]))
    # Each direction is split between the shapes on its own, so the mix holds
    # within the up set and within the down set rather than only in total. The
    # indices were drawn at random already, so handing the shapes out in blocks
    # still lands them on random features.
    plant_pat <- c(rep(names(mix), sa_sim_allocate(n_up, mix)),
                   rep(names(mix), sa_sim_allocate(n_down, mix)))

    for (k in seq_along(plant_idx)) {
      i <- plant_idx[k]
      delta[i, -1L] <- sa_sim_pattern_delta(plant_mag[k], plant_pat[k],
                                            n_treat_groups)
      pattern[i] <- plant_pat[k]
    }
  }

  center <- baseline + delta
  subjects <- if (paired) paste0("subject_", seq_len(sizes[1])) else NULL
  # One offset per subject and feature, drawn before the conditions and added to
  # all of them. Drawing it inside the condition loop would make it noise rather
  # than a subject effect, and the within-subject tests would have nothing to
  # gain over the independent ones.
  offsets <- if (paired) {
    vapply(seq_len(n_feats), function(i) {
      stats::rnorm(sizes[1], 0, sd_subject[i])
    }, numeric(sizes[1]))
  } else {
    NULL
  }

  blocks <- lapply(seq_len(n_lv), function(g) {
    values <- vapply(seq_len(n_feats), function(i) {
      stats::rnorm(sizes[g], mean = center[i, g], sd = sd_mat[i, g])
    }, numeric(sizes[g]))
    if (paired) values + offsets else values
  })

  data <- as.data.frame(do.call(rbind, blocks))
  names(data) <- feats
  rownames(data) <- NULL

  args <- list(
    data     = data,
    feats    = feats,
    group    = rep(group_lv, times = sizes),
    group_lv = group_lv
  )
  if (paired) {
    args$id <- rep(subjects, times = n_lv)
    args$paired <- TRUE
  }
  args$input_scale <- "log2"

  list(
    args           = args,
    truth          = sa_sim_truth(feats, delta, group_lv, pattern, direction,
                                  baseline, sd_subject),
    truth_group    = sa_sim_truth_group(feats, delta, center, sd_mat, group_lv,
                                        sizes),
    truth_contrast = sa_sim_truth_contrast(feats, delta, group_lv)
  )
}


#' Work out the group labels and the size of each one
#'
#' The count, the labels and the sizes all come from the same pair of arguments,
#' so they are settled together rather than in two passes that could disagree.
#' `n_treat` carries one size per treatment group, which makes its length the
#' number of them; `group_lv` carries the labels, which makes its length say the
#' same thing. When both are given they have to agree, and the one case where
#' they need not is a single size, which has an obvious number of groups to be
#' spread over as soon as the labels say how many there are.
#'
#' @param n_control The argument as received.
#' @param n_treat The argument as received.
#' @param group_lv The argument as received, possibly `NULL`.
#' @param use_default Result of `missing(n_treat)` in the calling function.
#' @param paired Whether the levels are repeated conditions.
#'
#' @return `list(group_lv = , sizes = )`, the sizes in level order with the
#'   control first.
#'
#' @keywords internal
#' @noRd
sa_sim_design <- function(n_control, n_treat, group_lv, use_default, paired) {
  if (!is.numeric(n_treat) || length(n_treat) == 0L) {
    stop("`n_treat` must be one or more group sizes, one per treatment group.",
         call. = FALSE)
  }

  if (is.null(group_lv)) {
    if (length(n_treat) < 2L) {
      stop("`n_treat` holds one size per treatment group, and there must be ",
           "at least two of them for a comparison of three or more levels. ",
           "Pass a size per group, such as `n_treat = rep(", n_treat[1],
           ", 3)`, or use simulate_two_groups() for two groups in all.",
           call. = FALSE)
    }
    group_lv <- c("control", paste0("treat_", seq_along(n_treat)))
  } else {
    if (!is.character(group_lv) || length(group_lv) < 3L || anyNA(group_lv) ||
        anyDuplicated(group_lv) > 0L) {
      stop("`group_lv` must be at least three distinct non-missing group ",
           "labels, the first being the control.", call. = FALSE)
    }
    n_wanted <- length(group_lv) - 1L
    # Labels say how many groups there are, so one size has somewhere to go.
    # The default is treated as one size for the same reason: it says how big a
    # group should be, not how many of them the caller wanted.
    if (use_default || length(n_treat) == 1L) {
      n_treat <- rep(n_treat[1], n_wanted)
    }
    if (length(n_treat) != n_wanted) {
      stop("`group_lv` names ", n_wanted, " treatment group(s) after the ",
           "control, but `n_treat` gives ", length(n_treat), " size(s).",
           call. = FALSE)
    }
  }

  sizes <- c(sa_check_count(n_control, "n_control", 2),
             vapply(seq_along(n_treat), function(k) {
               sa_check_count(n_treat[k], paste0("n_treat[", k, "]"), 2)
             }, integer(1)))

  if (paired && length(unique(sizes)) > 1L) {
    stop("`paired = TRUE` measures every condition on the same subjects, so ",
         "every group holds the same number of them, but the sizes given are ",
         paste(sizes, collapse = ", "), ".", call. = FALSE)
  }

  list(group_lv = group_lv, sizes = sizes)
}


#' Check the shape weights and drop the ones set to zero
#'
#' `simulate_factorial_groups()` weights two sets of shapes rather than one, and
#' the four things that make a set of weights wrong are the same for both, so the
#' catalogue and the argument name are parameters rather than a second copy of
#' this function.
#'
#' @param known The shape names the caller accepts.
#' @param arg Argument name to name in the error.
#'
#' @keywords internal
#' @noRd
sa_sim_pattern_mix <- function(pattern_mix,
                               known = c("all", "gradient", "single"),
                               arg = "pattern_mix") {
  if (!is.numeric(pattern_mix) || length(pattern_mix) == 0L ||
      is.null(names(pattern_mix)) || anyNA(pattern_mix) ||
      anyDuplicated(names(pattern_mix)) > 0L) {
    stop("`", arg, "` must be a named numeric vector of weights with one ",
         "entry per shape and no duplicates. Known shapes are: ",
         paste(known, collapse = ", "), ".", call. = FALSE)
  }
  unknown <- setdiff(names(pattern_mix), known)
  if (length(unknown) > 0L) {
    stop("`", arg, "` names unknown shape(s): ",
         paste(unknown, collapse = ", "), ". Known shapes are: ",
         paste(known, collapse = ", "), ".", call. = FALSE)
  }
  if (any(pattern_mix < 0)) {
    stop("`", arg, "` weights must not be negative.", call. = FALSE)
  }
  if (sum(pattern_mix) <= 0) {
    stop("`", arg, "` needs at least one positive weight, otherwise there ",
         "is no shape left to plant an effect in.", call. = FALSE)
  }
  pattern_mix[pattern_mix > 0]
}


#' Split `n` between weighted shapes without drawing lots
#'
#' The largest remainder method rather than [stats::rmultinom()], so that the
#' counts are exactly the proportions the weights ask for and are a function of
#' the arguments alone. A random split would make the number of features of each
#' shape move with the seed, which is the one thing about a simulation that
#' should not have to be looked up.
#'
#' @keywords internal
#' @noRd
sa_sim_allocate <- function(n, weights) {
  out <- stats::setNames(integer(length(weights)), names(weights))
  if (n == 0L) {
    return(out)
  }
  share <- n * weights / sum(weights)
  out[] <- as.integer(floor(share))
  short <- n - sum(out)
  if (short > 0L) {
    # `order()` leaves ties in the order the weights were given, so an even mix
    # hands the remainder to the earlier shapes rather than to an arbitrary one.
    take <- utils::head(order(share - out, decreasing = TRUE), short)
    out[take] <- out[take] + 1L
  }
  out
}


#' Spread one magnitude over the treatment groups according to its shape
#'
#' @param d Signed magnitude of the effect, on the log2 scale.
#' @param pattern `"all"`, `"gradient"` or `"single"`.
#' @param n_groups Number of treatment groups.
#'
#' @return Numeric vector of length `n_groups`, one delta per treatment group.
#'
#' @keywords internal
#' @noRd
sa_sim_pattern_delta <- function(d, pattern, n_groups) {
  switch(
    pattern,
    all      = rep(d, n_groups),
    gradient = d * seq_len(n_groups) / n_groups,
    single   = {
      out <- numeric(n_groups)
      out[sample.int(n_groups, 1L)] <- d
      out
    },
    stop("internal error: unknown effect shape `", pattern, "`.",
         call. = FALSE)
  )
}


#' Feature-level answer, aligned with the `effect` table
#'
#' @param delta Features by levels, the control column being zero throughout.
#'
#' @keywords internal
#' @noRd
sa_sim_truth <- function(feats, delta, group_lv, pattern, direction, baseline,
                         sd_subject) {
  treat_delta <- delta[, -1L, drop = FALSE]
  abs_delta <- abs(treat_delta)
  largest <- apply(abs_delta, 1, max)
  # A shape that moves every treatment group by the same amount leaves no single
  # level furthest from the control, so the tie is reported rather than broken
  # silently and scored on.
  tied <- rowSums(abs_delta == largest) > 1L
  which_max <- max.col(abs_delta, ties.method = "first")

  extreme_level <- group_lv[1L + which_max]
  extreme_level[largest == 0] <- NA_character_

  data.frame(
    features      = feats,
    pattern       = pattern,
    direction     = direction,
    extreme_level = extreme_level,
    extreme_tied  = tied,
    log2fc        = treat_delta[cbind(seq_along(feats), which_max)],
    baseline      = baseline,
    sd_subject    = sd_subject,
    stringsAsFactors = FALSE
  )
}


#' Per feature and level answer
#'
#' @keywords internal
#' @noRd
sa_sim_truth_group <- function(feats, delta, center, sd_mat, group_lv, sizes) {
  n_lv <- length(group_lv)
  data.frame(
    features = rep(feats, each = n_lv),
    group    = rep(group_lv, times = length(feats)),
    is_ref   = rep(seq_len(n_lv) == 1L, times = length(feats)),
    delta    = as.vector(t(delta)),
    center   = as.vector(t(center)),
    sd       = as.vector(t(sd_mat)),
    n        = rep(sizes, times = length(feats)),
    stringsAsFactors = FALSE
  )
}


#' Per feature and pair answer, in the post-hoc table's own direction
#'
#' The pairs come from `sa_level_pairs()` rather than from a second listing of
#' them here, so the row order and the `group1 - group2` direction cannot drift
#' apart from the tables this is meant to score.
#'
#' @keywords internal
#' @noRd
sa_sim_truth_contrast <- function(feats, delta, group_lv) {
  pairs <- sa_level_pairs(group_lv)
  diffs <- delta[, pairs$i, drop = FALSE] - delta[, pairs$j, drop = FALSE]
  flat <- as.vector(t(diffs))

  data.frame(
    features = rep(feats, each = nrow(pairs)),
    contrast = rep(pairs$contrast, times = length(feats)),
    group1   = rep(pairs$group1, times = length(feats)),
    group2   = rep(pairs$group2, times = length(feats)),
    delta    = flat,
    is_diff  = flat != 0,
    stringsAsFactors = FALSE
  )
}
