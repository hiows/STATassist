# The clustering kernel, written to the same rule as `kernel_performance.R`: plain
# input in, a numeric vector out, no fitted object kept anywhere and nothing said
# to the user. There is one of them, because there is one number all four methods
# can be read on.
#
# It is written out here rather than taken from the `cluster` package for the
# reason the robust and performance kernels are written out. The definition is
# four lines of arithmetic over a distance matrix this package already has in its
# hand, `cluster::silhouette()` would be a dependency bought for those four lines,
# and its object is a matrix with an attribute rather than the vector the result
# needs. `scikit-learn` spells the same definition `silhouette_samples`, so a
# Python transcription has a counterpart to check against.

#' Silhouette width of every point
#'
#' The silhouette of a point compares how far it sits from its own cluster with how
#' far it sits from the nearest cluster it is not in: `(b - a) / max(a, b)`, where
#' `a` is the mean distance to the other members of its own cluster and `b` is the
#' smallest mean distance to the members of another. It runs from 1, meaning the
#' point is far closer to its own group than to any other, through 0 at the border
#' between two, to -1 for a point that would be better off elsewhere.
#'
#' Three conventions, all Rousseeuw's:
#'
#' * A point alone in its cluster scores 0. It has no `a` to speak of, and the
#'   alternative — dividing by nothing and calling it 1 — would score a singleton
#'   as the best-placed point in the data.
#' * Noise scores `NA`, and takes no part in any other point's `a` or `b`. It is
#'   not a cluster, so a point cannot be near to it in the sense `b` measures.
#' * A single cluster scores `NA` throughout. There is no other cluster for `b` to
#'   be about, and the width is a comparison rather than a measurement.
#'
#' @param d A [stats::dist()] object over the points, in the order `cluster` is in.
#' @param cluster Integer cluster label per point, `0` for noise.
#'
#' @return Numeric vector, one silhouette width per point, `NA` where undefined.
#'
#' @keywords internal
#' @noRd
sa_silhouette <- function(d, cluster) {
  cluster <- as.integer(cluster)
  out <- rep(NA_real_, length(cluster))

  assigned <- cluster > 0L
  ids <- sort(unique(cluster[assigned]))
  if (length(ids) < 2L) {
    return(out)
  }

  m <- as.matrix(d)
  # One logical row per cluster, so the membership of each is worked out once
  # rather than once per point.
  members <- lapply(ids, function(g) assigned & cluster == g)

  for (i in which(assigned)) {
    own <- match(cluster[i], ids)
    same <- members[[own]]
    same[i] <- FALSE
    if (!any(same)) {
      out[i] <- 0
      next
    }
    a <- mean(m[i, same])
    b <- min(vapply(members[-own], function(g) mean(m[i, g]), numeric(1)))
    scale <- max(a, b)
    # Coincident points give a == b == 0, which is a tie rather than a division.
    out[i] <- if (scale > 0) (b - a) / scale else 0
  }

  out
}
