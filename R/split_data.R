# The train/test partition for the supervised learning family. A split is the
# first place a model can be handed something it is not supposed to know, and
# the two ways that happens are structural rather than accidental, so both are
# dealt with here instead of being left to the caller: `stratified` keeps the
# balance of the whole data set in both halves, and `id` keeps every row of one
# sampling unit on the same side of the split.

#' Fold rows into sampling units and carry the stratifier up with them
#'
#' Units come out in order of first appearance rather than sorted, so a numeric
#' id is not silently reordered as text, which is the same rule
#' `sa_align_by_subject()` follows.
#'
#' A unit is assigned to one side of the split as a whole, so it can only carry
#' one stratum. A unit whose rows disagree is an error rather than a majority
#' vote: the two are different designs, and guessing which one was meant would
#' produce a split that looks stratified but is not.
#'
#' @return List with `rows`, the row indices of each unit, and `stratum`, one
#'   entry per unit in the type it arrived as, or `NULL` when there is none.
#'
#' @keywords internal
#' @noRd
sa_fold_by_unit <- function(id, stratified, n) {
  if (is.null(id)) {
    return(list(rows = split(seq_len(n), factor(seq_len(n), levels = seq_len(n))),
                stratum = stratified))
  }

  units <- unique(id)
  key <- match(id, units)
  rows <- split(seq_len(n), factor(key, levels = seq_along(units)))
  names(rows) <- as.character(units)

  stratum <- NULL
  if (!is.null(stratified)) {
    mixed <- vapply(rows, function(r) length(unique(stratified[r])) > 1L,
                    logical(1))
    if (any(mixed)) {
      shown <- names(rows)[mixed]
      stop("`stratified` must be constant within each `id`, since a unit is ",
           "assigned to one side of the split as a whole. Offending id(s): ",
           paste(utils::head(shown, 5L), collapse = ", "),
           if (length(shown) > 5L) ", ..." , ".", call. = FALSE)
    }
    stratum <- stratified[vapply(rows, function(r) r[1L], integer(1))]
  }

  list(rows = rows, stratum = stratum)
}


#' Split data into training and test sets
#'
#' Partitions the rows of a data set into a training half and a test half,
#' optionally several times over. The partition is stratified, so the balance of
#' the whole data set is preserved in both halves rather than left to the draw,
#' and it can be taken over sampling units rather than over rows, so that
#' repeated measurements of one subject never end up on both sides.
#'
#' Nothing is fitted or transformed here. The point of splitting first is that
#' every later step — imputation, scaling, feature selection, hyperparameter
#' tuning — is fitted on the training half alone, and this function is what
#' makes that half well defined.
#'
#' @details
#' The partition itself is drawn by [caret::createDataPartition()], which takes
#' `ceiling(n * p_train)` observations from within each stratum. What it does
#' with `stratified` depends on its type, and the difference is worth knowing
#' because it is invisible from the outside:
#'
#' \describe{
#'   \item{Factor or character}{Each distinct value is a stratum. A value
#'     occurring only once is put into the training set and the test set has
#'     none of it, which `createDataPartition()` warns about.}
#'   \item{Numeric}{Cut into up to five quantile bins first, and the bins are
#'     the strata. This is how a continuous outcome is kept from landing
#'     entirely on one side, but it does mean a numeric stratifier is never
#'     matched exactly.}
#'   \item{`NULL`}{No strata. The split is a simple random draw of
#'     `ceiling(nrow(data) * p_train)` rows.}
#' }
#'
#' With `id` the whole thing moves up one level. Rows are folded into units,
#' each unit takes the stratum its rows agree on, the partition is drawn over
#' units, and the chosen units are expanded back into row indices. `p_train` is
#' then a proportion of units, not of rows, and the two differ whenever the
#' units have unequal sizes; `parameters$achieved_p` reports the row proportion
#' each repeat actually reached.
#'
#' @param data A data.frame (or matrix) in wide format, one row per
#'   observation.
#' @param stratified What to preserve the balance of, either the name of a
#'   column of `data` or a vector with one entry per row. Usually the outcome
#'   the model will predict. `NULL` draws a simple random split.
#' @param id Sampling unit, either a column name or a vector with one entry per
#'   row. Rows sharing a value are assigned to the same side of the split, which
#'   is what keeps repeated measurements or technical replicates of one subject
#'   from appearing in both halves. `NULL` treats every row as its own unit.
#' @param p_train Proportion allocated to the training set, strictly between 0
#'   and 1.
#' @param times Number of independent splits to draw. The shape of the result
#'   does not depend on it: one split still comes back as a list of one.
#' @param seed Seed for the draw, or `NULL` to use the stream as it stands.
#'   Supplying one does not disturb the caller: the previous random number state
#'   is put back when the function returns.
#'
#' @return An object of class `sa_split`, a plain list of five elements.
#'
#'   \describe{
#'     \item{`full_data`}{The input, exactly as it was passed in. Held once
#'       rather than once per repeat.}
#'     \item{`datasets`}{One element per repeat, named `Resample1` upwards, each
#'       a list of `train_data`, `test_data`, and the `train_rows` and
#'       `test_rows` they were taken from. The two frames have their row names
#'       reset, so the row numbers are the only record of where a row came from
#'       and they are kept beside it.}
#'     \item{`train_idx`}{The training row indices of every repeat, in the form
#'       `caret::createDataPartition(list = TRUE)` returns. A matrix is not
#'       offered because units of unequal size make the repeats different
#'       lengths.}
#'     \item{`design`}{What the split was made on: row and unit counts, the
#'       labels of `stratified` and `id`, and the number of units per stratum
#'       where the strata are discrete.}
#'     \item{`parameters`}{`p_train`, `times`, `seed`, and `achieved_p`, the row
#'       proportion each repeat actually reached.}
#'   }
#'
#' @seealso [caret::createDataPartition()], which draws the partition.
#'
#' @examples
#' ## Stratified on the outcome: both halves keep the species balance of the
#' ## whole data set rather than whatever the draw happened to give.
#' sp <- split_data(iris, stratified = "Species", seed = 1)
#' sp
#' table(sp$datasets[[1]]$train_data$Species)
#' table(sp$datasets[[1]]$test_data$Species)
#'
#' ## Three splits at once, for a repeated hold-out.
#' sp3 <- split_data(iris, stratified = "Species", times = 3, seed = 1)
#' vapply(sp3$datasets, function(d) nrow(d$train_data), numeric(1))
#'
#' ## Three measurements per subject. Splitting by row would put most subjects
#' ## in both halves; splitting by `id` cannot.
#' rep_data <- data.frame(
#'   subject = rep(paste0("s", 1:20), each = 3),
#'   arm     = rep(c("control", "treated"), each = 30),
#'   value   = seq_len(60)
#' )
#' sp_id <- split_data(rep_data, stratified = "arm", id = "subject", seed = 1)
#' intersect(sp_id$datasets[[1]]$train_data$subject,
#'           sp_id$datasets[[1]]$test_data$subject)
#'
#' ## `p_train` is a proportion of units once `id` is given, so the row
#' ## proportion it reaches is reported rather than assumed.
#' sp_id$parameters$achieved_p
#'
#' @export
split_data <- function(data,
                       stratified = NULL,
                       id = NULL,
                       p_train = 0.75,
                       times = 1,
                       seed = NULL) {

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or a matrix.", call. = FALSE)
  }
  n <- nrow(data)
  if (n == 0L) {
    stop("`data` has zero rows.", call. = FALSE)
  }

  sa_check_scalar_num(p_train, "p_train", 0, 1,
                      lower_open = TRUE, upper_open = TRUE)
  times <- sa_check_count(times, "times", 1)

  strat <- sa_resolve_row_vector(stratified, "stratified", data)
  unit <- sa_resolve_row_vector(id, "id", data)

  folded <- sa_fold_by_unit(unit$value, strat$value, n)
  unit_rows <- folded$rows
  n_units <- length(unit_rows)
  if (n_units < 2L) {
    stop("a split needs at least 2 sampling units, but `data` has ", n_units,
         ".", call. = FALSE)
  }

  # `createDataPartition()` needs something to stratify on either way. A single
  # level is the honest spelling of "no strata": every unit is in the one class
  # and the draw is a simple random sample of it.
  strata <- folded$stratum
  if (is.null(strata)) {
    y <- factor(rep("all", n_units))
  } else if (is.numeric(strata)) {
    if (length(unique(strata)) < 2L) {
      stop("`stratified` is numeric and constant, so it defines no strata. ",
           "Pass `stratified = NULL` for an unstratified split.", call. = FALSE)
    }
    y <- strata
  } else {
    # Dropped to the levels actually present. An unused factor level is a
    # property of how the column was built, not a stratum with no members.
    y <- factor(as.character(strata))
  }

  restore_seed <- sa_preserve_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  unit_idx <- caret::createDataPartition(y = y, times = times, p = p_train,
                                         list = FALSE)

  train_rows <- lapply(seq_len(times), function(i) {
    sort(unlist(unit_rows[unit_idx[, i]], use.names = FALSE))
  })
  names(train_rows) <- colnames(unit_idx)

  n_test <- vapply(train_rows, function(rows) n - length(rows), numeric(1))
  if (any(n_test == 0)) {
    stop("`p_train` = ", p_train, " leaves the test set empty: ",
         "the partition takes `ceiling(n * p_train)` units from each stratum, ",
         "which here is all ", n_units, " of them. Lower `p_train` or supply ",
         "more units.", call. = FALSE)
  }

  datasets <- lapply(train_rows, function(rows) {
    test <- setdiff(seq_len(n), rows)
    train_data <- data[rows, , drop = FALSE]
    test_data <- data[test, , drop = FALSE]
    rownames(train_data) <- NULL
    rownames(test_data) <- NULL

    list(
      train_data = train_data,
      test_data  = test_data,
      train_rows = rows,
      test_rows  = test
    )
  })

  structure(
    list(
      full_data = data,
      datasets  = datasets,
      train_idx = train_rows,
      design    = list(
        n_rows     = n,
        n_units    = n_units,
        stratified = strat$label,
        id         = unit$label,
        # A numeric stratifier is binned inside `createDataPartition()` and the
        # bins are not reported back, so there is no honest count to give.
        strata_n   = if (is.factor(y) && !is.null(strata)) c(table(y)) else NULL
      ),
      parameters = list(
        p_train    = p_train,
        times      = times,
        seed       = seed,
        achieved_p = vapply(train_rows, function(rows) length(rows) / n,
                            numeric(1))
      ),
      metadata = sa_metadata()
    ),
    class = c("sa_split", "sa_result")
  )
}


#' Print a train/test split
#'
#' Summarises what the split was made on and how large each half came out,
#' rather than printing the data. The frames themselves are in
#' `x$datasets[[i]]$train_data` and `$test_data`.
#'
#' @param x A split, as returned by [split_data()].
#' @param ... Ignored, present for consistency with [print()].
#'
#' @return `x` invisibly.
#'
#' @examples
#' split_data(iris, stratified = "Species", times = 2, seed = 1)
#'
#' @export
print.sa_split <- function(x, ...) {
  design <- x$design
  params <- x$parameters

  cat("<sa_split> train/test partition\n")
  cat("  rows     : ", design$n_rows,
      if (!is.na(design$id)) {
        paste0("  (", design$n_units, " unit(s) of `", design$id, "`)")
      },
      "\n", sep = "")
  cat("  stratify : ",
      if (is.na(design$stratified)) "none" else design$stratified, "\n",
      sep = "")
  if (!is.null(design$strata_n)) {
    cat("             ",
        paste0(names(design$strata_n), " ", design$strata_n, collapse = ", "),
        "\n", sep = "")
  }
  cat("  settings : p_train = ", params$p_train, ", times = ", params$times,
      if (!is.null(params$seed)) paste0(", seed = ", params$seed),
      "\n", sep = "")

  cat("\n  splits\n")
  width <- max(nchar(names(x$datasets)))
  for (nm in names(x$datasets)) {
    d <- x$datasets[[nm]]
    cat("    $", formatC(nm, width = -width), "  train ", nrow(d$train_data),
        " / test ", nrow(d$test_data),
        "  (p = ", format(round(params$achieved_p[[nm]], 3), nsmall = 3), ")\n",
        sep = "")
  }

  invisible(x)
}
