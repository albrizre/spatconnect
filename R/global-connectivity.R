#' Global spatial-connectivity analysis
#'
#' Computes the Betti-0 connectivity curve for areal data and assesses its
#' departure from spatial exchangeability using random relabeling.
#'
#' @param x Numeric vector with one value per areal unit.
#' @param graph Spatial adjacency information. Accepted formats are an
#'   `spdep` `nb` object, an `spdep` `listw` object, a two-column edge list,
#'   or an `n x n` adjacency matrix.
#' @param reference Reference level used to center `x`. Either a finite numeric
#'   value or one of `"mean"`, `"median"`, or `"zero"`.
#' @param scale Positive scale used for the deviations. Either a finite numeric
#'   value or one of `"sd"`, `"mad"`, or `"none"`.
#' @param side Analysis direction: `"superlevel"`, `"sublevel"`, or `"both"`.
#' @param n_perm Number of random relabelings used for inference.
#' @param seed Optional random seed.
#' @param conf Pointwise random-relabeling envelope level.
#' @param lambda_grid Optional non-negative threshold grid. When supplied, the
#'   integrated absolute discrepancy is evaluated by trapezoidal integration on
#'   this grid. This is useful for exact reproduction of analyses based on a
#'   fixed numerical grid.
#' @param n_grid Optional number of equally spaced thresholds from zero to the
#'   largest positive selected value. Ignored when `lambda_grid` is supplied.
#'   When both `lambda_grid` and `n_grid` are `NULL`, discrepancy integration is
#'   exact over the distinct positive data values.
#' @param keep_permutations Logical; if `TRUE`, retain the matrix of relabeled
#'   Betti-0 curves.
#'
#' @details
#' The superlevel analysis uses `z = (x - b) / s`. The sublevel analysis applies
#' the same superlevel construction to `-z`.
#'
#' The global Monte Carlo test treats the observed curve and the `n_perm`
#' relabeled curves as an exchangeable set of `n_perm + 1` curves. Every
#' integrated absolute discrepancy is computed from their common mean, and the
#' observed statistic is ranked among all `n_perm + 1` values. For graphical
#' presentation, the reported reference mean and pointwise envelope use only the
#' relabeled curves.
#'
#' @return For one side, an object of class `spatconnect_global`. For
#'   `side = "both"`, an object of class `spatconnect_global_both` containing
#'   `superlevel` and `sublevel` analyses.
#' @export
#'
#' @examples
#' x <- c(3, 2, 1, -1)
#' edges <- data.frame(i = c(1, 2, 3), j = c(2, 3, 4))
#' global_connectivity(x, edges, reference = 0, scale = "none",
#'                     side = "superlevel", n_perm = 99, seed = 1)
global_connectivity <- function(
    x,
    graph,
    reference = "mean",
    scale = "sd",
    side = c("both", "superlevel", "sublevel"),
    n_perm = 999,
    seed = NULL,
    conf = 0.95,
    lambda_grid = NULL,
    n_grid = NULL,
    keep_permutations = FALSE
) {
  side <- match.arg(side)

  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 1 || n_perm != floor(n_perm)) {
    stop("`n_perm` must be a positive integer.", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)

  if (length(conf) != 1L || !is.finite(conf) || conf <= 0 || conf >= 1) {
    stop("`conf` must lie strictly between 0 and 1.", call. = FALSE)
  }

  if (!is.null(lambda_grid)) {
    lambda_grid <- sort(unique(as.numeric(lambda_grid)))
    if (length(lambda_grid) < 2L || any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) {
      stop("`lambda_grid` must contain at least two finite non-negative thresholds.", call. = FALSE)
    }
  }

  if (!is.null(n_grid)) {
    if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 2 || n_grid != floor(n_grid)) {
      stop("`n_grid` must be NULL or an integer of at least 2.", call. = FALSE)
    }
    n_grid <- as.integer(n_grid)
  }

  prep <- .prepare_signal(x, reference = reference, scale = scale)
  edges <- .as_edge_list(graph, length(x))
  adj <- .adjacency_list(length(x), edges)

  run_one <- function(which_side) {
    y <- .side_vector(prep$z, which_side)
    .global_connectivity_one(
      y = y,
      x = prep$x,
      z = prep$z,
      edges = edges,
      adj = adj,
      side = which_side,
      reference = prep$reference,
      scale = prep$scale,
      n_perm = n_perm,
      seed = seed,
      conf = conf,
      lambda_grid = lambda_grid,
      n_grid = n_grid,
      keep_permutations = keep_permutations
    )
  }

  if (side == "both") {
    out <- list(
      superlevel = run_one("superlevel"),
      sublevel = run_one("sublevel"),
      reference = prep$reference,
      scale = prep$scale,
      n_areas = length(x),
      n_edges = nrow(edges)
    )
    class(out) <- "spatconnect_global_both"
    return(out)
  }

  run_one(side)
}

.global_connectivity_one <- function(
    y, x, z, edges, adj, side, reference, scale,
    n_perm, seed, conf, lambda_grid, n_grid, keep_permutations
) {
  if (!is.null(seed)) set.seed(seed)

  n <- length(y)
  positive_thresholds <- sort(unique(y[is.finite(y) & y > 0]))

  if (!is.null(lambda_grid)) {
    lambda_eval <- lambda_grid
    integration_method <- "trapezoidal_user_grid"
  } else if (!is.null(n_grid)) {
    ymax <- max(y, na.rm = TRUE)
    if (is.finite(ymax) && ymax > 0) {
      lambda_eval <- seq(0, ymax, length.out = n_grid)
    } else {
      lambda_eval <- seq(0, 1, length.out = n_grid)
    }
    integration_method <- paste0("trapezoidal_regular_grid_", n_grid)
  } else {
    lambda_eval <- sort(unique(c(0, positive_thresholds)))
    if (length(lambda_eval) == 1L) lambda_eval <- c(0, 1)
    integration_method <- "exact_positive_breakpoints"
  }

  beta_obs <- .beta0_values(y, edges, lambda_eval, adj = adj)
  beta_perm <- matrix(NA_real_, nrow = length(lambda_eval), ncol = n_perm)

  for (r in seq_len(n_perm)) {
    y_r <- sample(y, size = n, replace = FALSE)
    beta_perm[, r] <- .beta0_values(y_r, edges, lambda_eval, adj = adj)
  }

  alpha <- (1 - conf) / 2
  null_mean <- rowMeans(beta_perm)
  null_lower <- apply(beta_perm, 1, stats::quantile, probs = alpha, names = FALSE)
  null_upper <- apply(beta_perm, 1, stats::quantile, probs = 1 - alpha, names = FALSE)

  beta_all <- cbind(beta_obs, beta_perm)
  common_mean <- rowMeans(beta_all)

  if (integration_method == "exact_positive_breakpoints") {
    if (length(positive_thresholds) == 0L) {
      D_all <- rep(0, n_perm + 1L)
    } else {
      idx <- match(positive_thresholds, lambda_eval)
      widths <- diff(c(0, positive_thresholds))
      curves_event <- beta_all[idx, , drop = FALSE]
      mean_event <- rowMeans(curves_event)
      D_all <- colSums(abs(curves_event - mean_event) * widths)
    }
  } else {
    D_all <- apply(
      beta_all,
      2,
      function(curve) .trapz(abs(curve - common_mean), lambda_eval)
    )
  }

  D_obs <- D_all[1]
  D_perm <- D_all[-1]
  p_value <- mean(D_all >= D_obs)

  curve <- data.frame(
    lambda = lambda_eval,
    threshold_original = .original_threshold(
      lambda_eval, b = reference, s = scale, side = side
    ),
    observed = beta_obs,
    relabel_mean = null_mean,
    relabel_lower = null_lower,
    relabel_upper = null_upper,
    common_mean = common_mean
  )

  out <- list(
    side = side,
    n_areas = n,
    n_edges = nrow(edges),
    reference = reference,
    scale = scale,
    x = x,
    z = z,
    y = y,
    curve = curve,
    integrated_B = .integrated_beta0_exact(y, edges, adj = adj),
    discrepancy = D_obs,
    discrepancy_relabelings = D_perm,
    discrepancy_all = D_all,
    p_value = p_value,
    n_perm = n_perm,
    confidence = conf,
    integration_method = integration_method,
    global_calibration = "common_mean_R_plus_1"
  )

  if (isTRUE(keep_permutations)) out$beta0_relabelings <- beta_perm

  class(out) <- "spatconnect_global"
  out
}
