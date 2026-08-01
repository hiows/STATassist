# Internal input validation helpers shared by the exported functions. All checks
# run before any computation starts so that a bad call fails at the boundary
# instead of halfway through a feature loop.

#' Check that `feats` is a usable vector of feature names
#'
#' Shared by the functions that take `data` and by those that only take
#' per-feature vectors, so that a feature name is rejected for the same reasons
#' everywhere.
#'
#' @keywords internal
#' @noRd
sa_check_feat_names <- function(feats) {
  if (!is.character(feats) || length(feats) == 0L) {
    stop("`feats` must be a non-empty character vector of feature names.",
         call. = FALSE)
  }
  if (anyNA(feats)) {
    stop("`feats` must not contain NA.", call. = FALSE)
  }
  dup_feats <- unique(feats[duplicated(feats)])
  if (length(dup_feats) > 0L) {
    stop("`feats` contains duplicated names: ",
         paste(dup_feats, collapse = ", "), call. = FALSE)
  }
  invisible(feats)
}


#' Validate wide-format grouped input
#'
#' @param data Wide data.frame, one row per observation.
#' @param feats Character vector of numeric column names.
#' @param group Grouping vector, `length(group) == nrow(data)`.
#' @param group_lv Group levels to keep, in display order.
#' @param id Optional pairing key, filtered alongside `data`.
#' @param n_levels If not `NULL`, `group_lv` must have exactly this length.
#' @param min_levels Minimum acceptable number of levels.
#'
#' @return List with the row-filtered `data`, `feats`, `id`, a `group` factor
#'   whose levels are `group_lv`, and `n_dropped` rows removed for belonging to
#'   a level outside `group_lv`.
#'
#' @keywords internal
#' @noRd
sa_validate_wide_input <- function(data,
                                   feats,
                                   group,
                                   group_lv,
                                   id = NULL,
                                   n_levels = NULL,
                                   min_levels = 2L) {
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or a matrix.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("`data` has zero rows.", call. = FALSE)
  }

  sa_check_feat_names(feats)
  unknown_feats <- setdiff(feats, names(data))
  if (length(unknown_feats) > 0L) {
    stop("`feats` not found in `data`: ",
         paste(unknown_feats, collapse = ", "), call. = FALSE)
  }
  non_numeric <- feats[!vapply(data[feats], is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop("`feats` must refer to numeric columns. Not numeric: ",
         paste(non_numeric, collapse = ", "), call. = FALSE)
  }

  if (length(group) != nrow(data)) {
    stop("`group` must have one entry per row of `data`: got ",
         length(group), " for ", nrow(data), " rows.", call. = FALSE)
  }
  group_chr <- as.character(group)

  if (!is.null(id) && length(id) != nrow(data)) {
    stop("`id` must have one entry per row of `data`: got ",
         length(id), " for ", nrow(data), " rows.", call. = FALSE)
  }

  if (length(group_lv) == 0L) {
    stop("`group_lv` must be a non-empty vector of group levels.",
         call. = FALSE)
  }
  group_lv <- as.character(group_lv)
  if (anyNA(group_lv)) {
    stop("`group_lv` must not contain NA.", call. = FALSE)
  }
  dup_lv <- unique(group_lv[duplicated(group_lv)])
  if (length(dup_lv) > 0L) {
    stop("`group_lv` contains duplicated levels: ",
         paste(dup_lv, collapse = ", "), call. = FALSE)
  }
  if (!is.null(n_levels) && length(group_lv) != n_levels) {
    stop("`group_lv` must contain exactly ", n_levels, " levels, but ",
         length(group_lv), " were given: ",
         paste(group_lv, collapse = ", "), call. = FALSE)
  }
  if (length(group_lv) < min_levels) {
    stop("`group_lv` must contain at least ", min_levels, " levels.",
         call. = FALSE)
  }
  absent_lv <- group_lv[!group_lv %in% group_chr]
  if (length(absent_lv) > 0L) {
    stop("`group_lv` level(s) absent from `group`: ",
         paste(absent_lv, collapse = ", "), call. = FALSE)
  }

  # Rows outside group_lv must be dropped rather than coerced to NA, otherwise
  # logical indexing on the factor injects NA values into the test samples.
  keep <- group_chr %in% group_lv
  n_dropped <- sum(!keep)

  list(
    data      = data[keep, , drop = FALSE],
    feats     = feats,
    group     = factor(group_chr[keep], levels = group_lv),
    id        = if (is.null(id)) NULL else as.character(id)[keep],
    n_dropped = n_dropped
  )
}


#' Check a length-one logical argument
#'
#' @keywords internal
#' @noRd
sa_check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}


#' Check a length-one numeric argument against a range
#'
#' @param lower_open,upper_open Treat the corresponding bound as exclusive.
#'
#' @keywords internal
#' @noRd
sa_check_scalar_num <- function(x,
                                arg,
                                lower = -Inf,
                                upper = Inf,
                                lower_open = FALSE,
                                upper_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be a single non-missing number.", call. = FALSE)
  }
  too_low <- if (lower_open) x <= lower else x < lower
  too_high <- if (upper_open) x >= upper else x > upper
  if (too_low || too_high) {
    stop("`", arg, "` must be in ",
         if (lower_open) "(" else "[", lower, ", ", upper,
         if (upper_open) ")" else "]", ", but is ", x, ".", call. = FALSE)
  }
  invisible(x)
}


#' Check a multiplicity adjustment method
#'
#' Validated against [stats::p.adjust.methods] rather than a hand-written list of
#' choices. A `match.arg()` list once carried a misspelled entry that got through
#' here and was rejected by [stats::p.adjust()] much later.
#'
#' @keywords internal
#' @noRd
sa_check_p_adjust <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !x %in% stats::p.adjust.methods) {
    stop("`", arg, "` must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".", call. = FALSE)
  }
  invisible(x)
}


#' Check a plot margin argument
#'
#' @keywords internal
#' @noRd
sa_check_margin <- function(x, arg = "margin") {
  if (!is.numeric(x) || length(x) != 4L || anyNA(x) || any(x < 0)) {
    stop("`", arg, "` must be a numeric vector of 4 non-negative values.",
         call. = FALSE)
  }
  invisible(x)
}


#' Check an optional axis range argument
#'
#' `NULL` means the range is derived from the data, so it is always accepted.
#'
#' @keywords internal
#' @noRd
sa_check_lim <- function(x, arg) {
  if (!is.null(x) &&
      (!is.numeric(x) || length(x) != 2L || !all(is.finite(x)))) {
    stop("`", arg, "` must be NULL or a finite numeric vector of length 2.",
         call. = FALSE)
  }
  invisible(x)
}


#' Pair up two groups by row order
#'
#' Used when no `id` is supplied. Order is all the information available, so
#' unequal group sizes mean the pairs cannot be formed at all.
#'
#' @keywords internal
#' @noRd
sa_pair_by_order <- function(group, group_lv) {
  idx_x <- which(group == group_lv[1])
  idx_y <- which(group == group_lv[2])
  if (length(idx_x) != length(idx_y)) {
    stop("`paired = TRUE` without `id` pairs observations by row order, which ",
         "requires the same number of rows per group. Got ",
         group_lv[1], " = ", length(idx_x), ", ",
         group_lv[2], " = ", length(idx_y),
         ". Supply `id` to match on a pairing key instead.", call. = FALSE)
  }
  list(idx_x = idx_x, idx_y = idx_y, unmatched = character(0))
}


#' Pair up two groups by an explicit pairing key
#'
#' Row order carries no meaning here, so a reordered or partially incomplete
#' data set still yields the correct pairs. Pairs come out in the row order of
#' `group_lv[1]` so the result is deterministic.
#'
#' @return `idx_x` / `idx_y` are equal-length row indices forming the pairs, and
#'   `unmatched` lists ids that appear in only one group.
#'
#' @keywords internal
#' @noRd
sa_pair_by_id <- function(id, group, group_lv) {
  if (anyNA(id)) {
    stop("`id` must not contain NA when it is used to form pairs.",
         call. = FALSE)
  }
  idx_x <- which(group == group_lv[1])
  idx_y <- which(group == group_lv[2])
  id_x <- id[idx_x]
  id_y <- id[idx_y]

  repeated <- union(unique(id_x[duplicated(id_x)]),
                    unique(id_y[duplicated(id_y)]))
  if (length(repeated) > 0L) {
    stop("`id` must be unique within each group, otherwise the pairing is ",
         "ambiguous. Repeated id(s): ",
         paste(repeated, collapse = ", "), ".", call. = FALSE)
  }

  common <- intersect(id_x, id_y)
  if (length(common) < 2L) {
    stop("only ", length(common), " id(s) appear in both `", group_lv[1],
         "` and `", group_lv[2], "`; at least 2 pairs are needed.",
         call. = FALSE)
  }

  list(
    idx_x     = idx_x[match(common, id_x)],
    idx_y     = idx_y[match(common, id_y)],
    unmatched = setdiff(union(id_x, id_y), common)
  )
}


#' Line up three or more repeated conditions by subject
#'
#' The many-condition counterpart of `sa_pair_by_id()`. A within-subject omnibus
#' test needs a complete rectangle, so subjects missing any condition are dropped
#' rather than partially used, and the caller is told how many went.
#'
#' Row order pairing is deliberately not offered here. With two groups it is at
#' least well defined; with `k` conditions it would also have to assume that the
#' groups are stored in the same subject order, and there is no way to notice
#' when they are not.
#'
#' @param id Subject identifier, one entry per row.
#' @param group Condition factor whose levels are `group_lv`.
#' @param group_lv Condition levels, in display order.
#'
#' @return List with `idx`, a subjects-by-conditions matrix of row indices whose
#'   columns follow `group_lv`, `subjects` in row order, and `unmatched`.
#'
#' @keywords internal
#' @noRd
sa_align_by_subject <- function(id, group, group_lv) {
  if (anyNA(id)) {
    stop("`id` must not contain NA when it is used to align conditions.",
         call. = FALSE)
  }
  per_level <- lapply(group_lv, function(lv) {
    rows <- which(group == lv)
    ids <- id[rows]
    repeated <- unique(ids[duplicated(ids)])
    if (length(repeated) > 0L) {
      stop("`id` must be unique within each condition, otherwise the design ",
           "is ambiguous. Repeated id(s) in `", lv, "`: ",
           paste(repeated, collapse = ", "), ".", call. = FALSE)
    }
    stats::setNames(rows, ids)
  })
  names(per_level) <- group_lv

  all_ids <- unique(unlist(lapply(per_level, names), use.names = FALSE))
  complete <- Reduce(intersect, lapply(per_level, names))
  # Ordered by first appearance rather than sorted, so the subject order is the
  # data's own and a numeric id is not silently reordered as text.
  complete <- all_ids[all_ids %in% complete]

  if (length(complete) < 2L) {
    stop("only ", length(complete), " subject(s) have all ", length(group_lv),
         " condition(s); at least 2 complete subjects are needed.",
         call. = FALSE)
  }

  idx <- vapply(per_level, function(rows) unname(rows[complete]),
                integer(length(complete)))
  idx <- matrix(idx, nrow = length(complete), ncol = length(group_lv),
                dimnames = list(complete, group_lv))

  list(idx = idx, subjects = complete,
       unmatched = setdiff(all_ids, complete))
}


#' Check a vector of p-values
#'
#' @keywords internal
#' @noRd
sa_check_pvalues <- function(pvalue, arg = "pvalue") {
  bad <- which(!is.na(pvalue) & (!is.finite(pvalue) | pvalue < 0 | pvalue > 1))
  if (length(bad) > 0L) {
    shown <- bad[seq_len(min(5L, length(bad)))]
    stop("`", arg, "` must lie in [0, 1]. Offending position(s): ",
         paste(shown, collapse = ", "),
         if (length(bad) > 5L) ", ..." , ".", call. = FALSE)
  }
  invisible(pvalue)
}


#' Build an all-NA result row with the expected names
#'
#' @keywords internal
#' @noRd
sa_na_row <- function(nms) {
  stats::setNames(rep(NA_real_, length(nms)), nms)
}


#' Assemble one result row from named scalars
#'
#' Engines attach their own names to the values they return ([stats::t.test()]
#' names its statistic `t` and its parameter `df`, for instance), and
#' `c(a = x)` on a named `x` yields `a.t` rather than `a`. Forcing every value
#' through `as.numeric()` first keeps the row names exactly as written here.
#'
#' @keywords internal
#' @noRd
sa_row <- function(...) {
  vals <- list(...)
  stats::setNames(
    vapply(vals, function(v) as.numeric(v)[1], numeric(1)),
    names(vals)
  )
}


#' Run one test across all features and assemble a result table
#'
#' A single unusable feature must not abort a scan over hundreds of them, so
#' errors become an all-NA row and are reported together in one warning at the
#' end. Engine warnings (ties, inexact p-values) are collected into one grouped
#' `message()` instead of one warning per feature.
#'
#' @param feats Feature names, one row of the output per entry.
#' @param columns Numeric column names `fun` is expected to return.
#' @param label Human readable test name used in the warning text.
#' @param fun Function of the feature index returning a named numeric vector.
#' @param p_adjust Method passed to [stats::p.adjust()], or `NULL` for a table
#'   that holds no p-value at all, such as the effect estimates.
#'
#' @return data.frame with `features`, `columns` and, unless `p_adjust` is
#'   `NULL`, `pval_adj`.
#'
#' @keywords internal
#' @noRd
sa_feature_table <- function(feats, columns, label, fun, p_adjust = "none") {
  failures <- character(0)
  notes <- character(0)

  rows <- lapply(seq_along(feats), function(i) {
    caught <- character(0)
    row <- withCallingHandlers(
      tryCatch(
        fun(i),
        error = function(e) {
          failures[[feats[i]]] <<- conditionMessage(e)
          sa_na_row(columns)
        }
      ),
      warning = function(w) {
        caught <<- c(caught, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    if (length(caught) > 0L) {
      notes[[feats[i]]] <<- paste(unique(caught), collapse = "; ")
    }
    absent <- setdiff(columns, names(row))
    if (length(absent) > 0L) {
      stop("internal error: ", label, " row for `", feats[i],
           "` is missing column(s): ", paste(absent, collapse = ", "),
           call. = FALSE)
    }
    row[columns]
  })

  out <- as.data.frame(do.call(rbind, rows))
  out <- cbind(features = feats, out)
  rownames(out) <- NULL
  if (!is.null(p_adjust)) {
    out <- sa_add_padj(out, p_adjust)
  }

  if (length(notes) > 0L) {
    grouped <- table(notes)
    message(label, ": engine note(s) for ", length(notes), " of ",
            length(feats), " feature(s):\n",
            paste0("  [", grouped, " feature(s)] ", names(grouped),
                   collapse = "\n"))
  }
  if (length(failures) > 0L) {
    warning(label, " could not be computed for ", length(failures), " of ",
            length(feats), " feature(s); those rows are NA:\n",
            paste0("  ", names(failures), ": ", failures, collapse = "\n"),
            call. = FALSE)
  }

  out
}


#' Append a multiplicity adjusted p-value next to `pval`
#'
#' @keywords internal
#' @noRd
sa_add_padj <- function(df, method) {
  df$pval_adj <- stats::p.adjust(df$pval, method = method)
  nms <- setdiff(names(df), "pval_adj")
  df[append(nms, "pval_adj", after = match("pval", nms))]
}


#' Every unordered pair of group levels, in display order
#'
#' Pairs are generated so that the first member is always the earlier level of
#' `group_lv`. Every post-hoc estimate therefore reads as
#' `group_lv[i] - group_lv[j]` with `i < j`, which is the same direction rule
#' `group_lv` already fixes for the two-group case.
#'
#' @keywords internal
#' @noRd
sa_level_pairs <- function(group_lv) {
  idx <- utils::combn(seq_along(group_lv), 2L)
  data.frame(
    i        = idx[1, ],
    j        = idx[2, ],
    group1   = group_lv[idx[1, ]],
    group2   = group_lv[idx[2, ]],
    contrast = paste(group_lv[idx[1, ]], group_lv[idx[2, ]], sep = " - "),
    stringsAsFactors = FALSE
  )
}


#' Run one post-hoc procedure across features and assemble a table
#'
#' The counterpart of `sa_feature_table()` for the pairwise stage. Two things
#' differ. A post-hoc row is a feature and a pair rather than a feature, so `fun`
#' returns all pairs of one feature at once, in the row order of
#' `sa_level_pairs()`. And a feature whose omnibus test was not significant is
#' left out of the table entirely rather than filled with `NA`, because an absent
#' row means "not asked" while an `NA` row means "asked and failed".
#'
#' @param feats Feature names to run, already reduced to those that qualified.
#' @param group_lv Group levels, fixing the pair order and the direction.
#' @param columns Numeric column names `fun` is expected to return per pair.
#' @param label Human readable procedure name used in the warning text.
#' @param fun Function of the feature name returning a matrix or data.frame with
#'   one row per pair and `columns` as its columns.
#' @param p_adjust Method passed to [stats::p.adjust()], applied across the pairs
#'   within each feature. The pairwise family is the set of contrasts of one
#'   feature, so adjusting across features here as well would mix two families.
#'
#' @return data.frame with the columns of `sa_posthoc_table_columns()`, or a
#'   zero-row frame with the same columns when `feats` is empty.
#'
#' @keywords internal
#' @noRd
sa_posthoc_table <- function(feats, group_lv, columns, label, fun,
                             p_adjust = "holm") {
  pairs <- sa_level_pairs(group_lv)
  n_pairs <- nrow(pairs)
  empty <- data.frame(features = character(0), pairs[0, c("contrast", "group1",
                                                          "group2")],
                      stringsAsFactors = FALSE)

  if (length(feats) == 0L) {
    for (nm in setdiff(sa_posthoc_table_columns(), names(empty))) {
      empty[[nm]] <- numeric(0)
    }
    return(empty[sa_posthoc_table_columns()])
  }

  failures <- character(0)
  blocks <- lapply(feats, function(f) {
    mat <- tryCatch(
      fun(f),
      error = function(e) {
        failures[[f]] <<- conditionMessage(e)
        matrix(NA_real_, nrow = n_pairs, ncol = length(columns),
               dimnames = list(NULL, columns))
      }
    )
    mat <- as.data.frame(mat)
    if (nrow(mat) != n_pairs) {
      stop("internal error: ", label, " returned ", nrow(mat), " row(s) for `",
           f, "`, expected ", n_pairs, ".", call. = FALSE)
    }
    absent <- setdiff(columns, names(mat))
    if (length(absent) > 0L) {
      stop("internal error: ", label, " table for `", f,
           "` is missing column(s): ", paste(absent, collapse = ", "),
           call. = FALSE)
    }
    out <- data.frame(
      features = rep(f, n_pairs),
      pairs[c("contrast", "group1", "group2")],
      mat[columns],
      stringsAsFactors = FALSE
    )
    out$pval_adj <- stats::p.adjust(out$pval, method = p_adjust)
    out
  })

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL

  if (length(failures) > 0L) {
    warning(label, " could not be computed for ", length(failures), " of ",
            length(feats), " feature(s); those rows are NA:\n",
            paste0("  ", names(failures), ": ", failures, collapse = "\n"),
            call. = FALSE)
  }

  out[sa_posthoc_table_columns()]
}
