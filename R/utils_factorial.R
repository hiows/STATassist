# Facts about a crossed design rather than about what was planted in one. How
# many cells there are, the order they are counted in, which terms a fully
# crossed model has and which cells each pairwise contrast averages are all
# answers that `simulate_factorial_groups()` and `compare_factorial_groups()`
# have to agree on exactly, since the first writes the answer key the second is
# scored against. Two copies of the enumeration would line up on a balanced
# two-by-two design and drift apart on anything larger, and the drift would show
# as a recall figure rather than as an error.
#
# This is the same reason `sa_level_pairs()` is shared by every post-hoc kernel:
# the direction a contrast reads in is one decision, made once.


#' Point each factor at the reference level it was told to hold
#'
#' `sa_control_first()` once per factor named. A crossed design has one reference
#' per factor rather than one in total, and the cell where every one of them
#' lands is the reference cell that `effect` and the cell labels are read
#' against. Naming a level here moves it to the front of its own factor and
#' leaves the other factors in the order they arrived, so pointing one factor of
#' three at its control is a sentence rather than a rewrite of all three.
#'
#' Both shapes say the same thing, one level name per factor keyed by the
#' factor's name: a list is what `factors` and `factor_lv` already are, and a
#' character vector is what one name per factor fits in, so neither is made to
#' pretend to be the other.
#'
#' @param factor_lv Named list of factors and their levels, already validated.
#' @param control_label Named list or named character vector holding the level to
#'   hold as the reference for some or all of the factors, or `NULL` to leave
#'   every factor as it arrived.
#' @param lv_source What the levels came from, `"factor_lv"` when the caller
#'   named them and `"factors"` when they were sorted out of the data, so that an
#'   error names the argument the missing level is missing from.
#'
#' @return `factor_lv`, each named factor re-pointed and the rest untouched.
#'
#' @keywords internal
#' @noRd
sa_fact_control_first <- function(factor_lv, control_label,
                                  lv_source = "factor_lv") {
  if (is.null(control_label)) {
    return(factor_lv)
  }
  if (is.list(control_label)) {
    bad <- names(control_label)[lengths(control_label) != 1L]
    if (length(bad) > 0L || length(control_label) == 0L) {
      stop("`control_label` must hold one level name per factor. Entry/entries ",
           "holding something else: ", paste(bad, collapse = ", "), ".",
           call. = FALSE)
    }
    control_label <- unlist(control_label, use.names = TRUE)
  }
  if (!is.atomic(control_label) || length(control_label) == 0L ||
        is.null(names(control_label)) || anyNA(names(control_label)) ||
        !all(nzchar(names(control_label))) ||
        anyDuplicated(names(control_label)) > 0L) {
    stop("`control_label` must be a named list or named character vector, one ",
         "level name per factor it points, with the factor's name as the name. ",
         "Naming a factor twice, or none, is not a direction.", call. = FALSE)
  }
  unknown <- setdiff(names(control_label), names(factor_lv))
  if (length(unknown) > 0L) {
    stop("`control_label` names factor(s) the design does not hold: ",
         paste(unknown, collapse = ", "), ". Present: ",
         paste(names(factor_lv), collapse = ", "), ".", call. = FALSE)
  }
  for (nm in names(control_label)) {
    factor_lv[[nm]] <- sa_control_first(
      factor_lv[[nm]], control_label[[nm]],
      arg    = paste0("control_label$", nm),
      lv_arg = paste0(lv_source, "$", nm)
    )
  }
  factor_lv
}


#' Cross a set of factors, allowing the empty set
#'
#' [expand.grid()] of nothing is a frame of one row and no columns, which is the
#' answer wanted when no factor is between or none is within: one combination,
#' holding no constraints. Written out rather than relied on, because the empty
#' case is the one every index built on it has to survive.
#'
#' @param lv_list Named list of level vectors, possibly empty.
#'
#' @return data.frame of level indices, one row per combination, in the order
#'   the first factor varying fastest.
#'
#' @keywords internal
#' @noRd
sa_fact_grid <- function(lv_list) {
  if (length(lv_list) == 0L) {
    return(data.frame(row.names = 1L))
  }
  expand.grid(lapply(lv_list, seq_along), KEEP.OUT.ATTRS = FALSE)
}


#' The readable label of every cell of a grid
#'
#' The level names of each factor joined by a dot, in `factor_lv` order, which
#' is what `truth_cell` and a post-hoc stratum are keyed on.
#'
#' @param factor_lv Named list of factors and their levels.
#' @param cells Grid of level indices, as `sa_fact_grid()` returns.
#'
#' @keywords internal
#' @noRd
sa_fact_cell_labels <- function(factor_lv, cells) {
  n <- nrow(cells)
  apply(
    vapply(names(factor_lv), function(f) factor_lv[[f]][cells[[f]]],
           character(n)),
    1, paste, collapse = "."
  )
}


#' Which cell of the grid each row of level indices sits in
#'
#' The grid counts the first factor fastest, so the cell number is the level
#' indices read as a mixed-radix number. Doing it by arithmetic rather than by
#' matching label strings keeps the numbering identical to `sa_fact_grid()`'s row
#' order by construction, and holds for the empty set of factors, where every row
#' belongs to the single cell.
#'
#' @param level_idx Integer matrix, one row per observation and one column per
#'   factor, in `dims` order.
#' @param dims Number of levels of each factor.
#'
#' @keywords internal
#' @noRd
sa_fact_cell_index <- function(level_idx, dims) {
  if (length(dims) == 0L) {
    return(rep(1L, nrow(level_idx)))
  }
  strides <- as.integer(cumprod(c(1L, dims))[seq_along(dims)])
  1L + as.vector((level_idx - 1L) %*% strides)
}


#' Every term a fully crossed model of these factors has
#'
#' Main effects first, then the interactions in increasing order, and within an
#' order in the order the factors were declared. That is the order an ANOVA table
#' lists them in, so `truth_term` needs no reordering to sit beside one.
#'
#' @return List of character vectors, each the factors a term is over.
#'
#' @keywords internal
#' @noRd
sa_fact_terms <- function(fac_names) {
  out <- list()
  for (m in seq_along(fac_names)) {
    out <- c(out, utils::combn(fac_names, m, simplify = FALSE))
  }
  out
}


#' The readable label of every model term
#'
#' `a:b`, the way [stats::terms()] writes an interaction, so a row of `$terms`
#' and a row of an `aov()` summary can be matched on the name alone.
#'
#' @param terms List of character vectors, as `sa_fact_terms()` returns.
#'
#' @keywords internal
#' @noRd
sa_fact_term_labels <- function(terms) {
  vapply(terms, paste, character(1), collapse = ":")
}


#' Every subset of a term, the empty one included
#'
#' @keywords internal
#' @noRd
sa_fact_subsets <- function(term) {
  out <- list(character(0))
  for (m in seq_along(term)) {
    out <- c(out, utils::combn(term, m, simplify = FALSE))
  }
  out
}


#' Anything below this is the rounding left over from averaging, not an effect
#'
#' @keywords internal
#' @noRd
sa_fact_tol <- function() 1e-8


#' The ANOVA component of one term, cell by cell
#'
#' The inclusion-exclusion form of the decomposition: the component of a term is
#' the cell values averaged down to that term, minus everything already accounted
#' for by its sub-terms. Computing it rather than remembering what was planted is
#' the point, and it is the same computation whether the values are the shifts a
#' simulation put in the cells or the centres a comparison measured out of them,
#' which is why this lives here rather than beside either caller.
#'
#' @param eff One value per cell, in `cells` order.
#' @param cells Grid of level indices, as `sa_fact_grid()` returns.
#' @param term Character vector of the factors the term is over.
#'
#' @keywords internal
#' @noRd
sa_fact_component <- function(eff, cells, term) {
  total <- numeric(length(eff))
  for (sub in sa_fact_subsets(term)) {
    sign <- (-1)^(length(term) - length(sub))
    total <- total + sign * sa_fact_collapse(eff, cells, sub)
  }
  total[abs(total) < sa_fact_tol()] <- 0
  total
}


#' Average the cells down to the levels of a few factors
#'
#' Unweighted: every cell counts once, whatever it holds. That is what makes the
#' decomposition a statement about the levels rather than about how many
#' observations happened to land in each, which is the same choice `ss_type =
#' "III"` and the marginal means of the post-hoc stage make.
#'
#' @param keep Factors to keep, possibly none, in which case this is the grand
#'   mean repeated.
#'
#' @keywords internal
#' @noRd
sa_fact_collapse <- function(eff, cells, keep) {
  if (length(keep) == 0L) {
    return(rep(mean(eff), length(eff)))
  }
  do.call(stats::ave,
          c(list(eff), unname(as.list(cells[keep])), list(FUN = mean)))
}


#' Index of the largest `|component|`, with near-ties taking the earlier cell
#'
#' Values within `sa_fact_tol()` of the running maximum are treated as equal, so
#' a two-level factor's `-d/2` / `+d/2` pair keeps the earlier (reference) cell
#' even when floating point has made one side a hair larger.
#'
#' @keywords internal
#' @noRd
sa_fact_first_max_abs <- function(comp) {
  best <- -Inf
  chosen <- NA_integer_
  for (i in seq_along(comp)) {
    magnitude <- abs(comp[[i]])
    if (!is.finite(magnitude)) {
      next
    }
    if (is.na(chosen) || magnitude > best + sa_fact_tol()) {
      best <- magnitude
      chosen <- i
    }
  }
  if (is.na(chosen)) 1L else chosen
}


#' The largest effect each term accounts for, with its sign
#'
#' One number per term out of the whole component vector, so that a term has an
#' effect size that can be put on an axis beside its p-value. The component of
#' largest absolute value is the one taken, since it is the cell the term moved
#' furthest and a term that moved nothing is exactly zero there.
#'
#' A component is a deviation from what the other terms already predict, not the
#' difference between two levels: the components of a two-level factor whose
#' levels differ by `d` are `-d/2` and `+d/2`. That is the quantity
#' `simulate_factorial_groups()` records in `truth_term$max_abs_delta`, and the
#' reason to keep it rather than rescale to a pairwise difference is that the two
#' tables are then the same number and a result can be scored against the answer
#' key term by term.
#'
#' @param eff One value per cell, in `cells` order, on whatever scale the caller
#'   wants the answer on.
#' @param cells Grid of level indices, as `sa_fact_grid()` returns.
#' @param terms List of character vectors, as `sa_fact_terms()` returns.
#'
#' @return Numeric vector, one entry per term, in `terms` order. Absolute values
#'   within `sa_fact_tol()` of the maximum are treated as a tie (exact or near),
#'   and the earlier cell wins — the earlier level of the first factor, which is
#'   the reference after `control_label`.
#'
#' @keywords internal
#' @noRd
sa_fact_term_effect <- function(eff, cells, terms) {
  vapply(terms, function(term) {
    comp <- sa_fact_component(eff, cells, term)
    if (all(is.na(comp))) {
      return(NA_real_)
    }
    comp[[sa_fact_first_max_abs(comp)]]
  }, numeric(1), USE.NAMES = FALSE)
}


#' The pairwise contrasts a crossed design has, and which cells each is over
#'
#' Built once from the design and reused for every feature, since which cells a
#' contrast averages is a fact about the layout rather than about the data in it.
#' The simulator turns these selections into `truth_contrast$delta` and the
#' comparison turns them into estimates, so the two tables merge on
#' `factor` / `stratum` / `contrast` without either side sorting the other.
#'
#' A factorial design has two pairwise questions per factor and both are here.
#' The marginal contrast averages the other factors away, which is what an
#' estimated marginal mean does and what a main effect is a statement about. The
#' simple effect holds them at one combination, which is the only one of the two
#' that says anything when an interaction is real.
#'
#' @param design List holding at least `factor_lv` and the `cells` grid.
#'
#' @return List with `table`, the row labels, and `sel1` / `sel2`, the cells
#'   whose mean each side of the contrast is.
#'
#' @keywords internal
#' @noRd
sa_fact_contrast_skeleton <- function(design) {
  cells <- design$cells
  factor_lv <- design$factor_lv

  fac <- character(0)
  stratum <- character(0)
  contrast <- character(0)
  group1 <- character(0)
  group2 <- character(0)
  sel1 <- list()
  sel2 <- list()

  for (f in names(factor_lv)) {
    lv <- factor_lv[[f]]
    pairs <- sa_level_pairs(lv)
    others <- setdiff(names(factor_lv), f)
    strata <- sa_fact_grid(factor_lv[others])

    # The marginal contrast first, then one block per combination of the other
    # factors, so a table read top to bottom moves from the main effect to the
    # simple effects it may be hiding.
    add <- function(label, at) {
      fac <<- c(fac, rep(f, nrow(pairs)))
      stratum <<- c(stratum, rep(label, nrow(pairs)))
      contrast <<- c(contrast, pairs$contrast)
      group1 <<- c(group1, pairs$group1)
      group2 <<- c(group2, pairs$group2)
      sel1 <<- c(sel1, at[pairs$i])
      sel2 <<- c(sel2, at[pairs$j])
    }

    add(NA_character_, lapply(seq_along(lv), function(j) which(cells[[f]] == j)))
    for (s in seq_len(nrow(strata))) {
      held <- Reduce(`&`, lapply(others, function(o) {
        cells[[o]] == strata[[o]][s]
      }))
      label <- paste(vapply(others, function(o) {
        factor_lv[[o]][strata[[o]][s]]
      }, character(1)), collapse = ".")
      add(label, lapply(seq_along(lv), function(j) {
        which(held & cells[[f]] == j)
      }))
    }
  }

  list(
    table = data.frame(factor = fac, stratum = stratum, contrast = contrast,
                       group1 = group1, group2 = group2,
                       stringsAsFactors = FALSE),
    sel1 = sel1,
    sel2 = sel2
  )
}
