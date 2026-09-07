# Internal connectivity algorithms -----------------------------------------

.make_dsu <- function(n) {
  parent <- seq_len(n)
  rank <- integer(n)

  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }

  union <- function(a, b) {
    ra <- find(a)
    rb <- find(b)
    if (ra == rb) return(FALSE)

    if (rank[ra] < rank[rb]) {
      tmp <- ra
      ra <- rb
      rb <- tmp
    }

    parent[rb] <<- ra
    if (rank[ra] == rank[rb]) rank[ra] <<- rank[ra] + 1L
    TRUE
  }

  list(find = find, union = union)
}

.adjacency_list <- function(n, edges) {
  adj <- vector("list", n)
  if (nrow(edges) == 0L) return(adj)

  for (r in seq_len(nrow(edges))) {
    i <- as.integer(edges$i[r])
    j <- as.integer(edges$j[r])
    adj[[i]] <- c(adj[[i]], j)
    adj[[j]] <- c(adj[[j]], i)
  }

  adj
}

.check_distinct_positive <- function(y) {
  positive <- y[is.finite(y) & y > 0]
  if (anyDuplicated(positive) > 0L) {
    stop(
      "The local decomposition encountered exactly tied positive selected values. ",
      "The current implementation assumes distinct active-side values so that ",
      "areas enter the filtration one at a time.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.beta0_values <- function(y, edges, lambda, adj = NULL) {
  y <- as.numeric(y)
  lambda <- as.numeric(lambda)
  n <- length(y)

  if (length(lambda) == 0L) return(numeric(0))
  if (any(!is.finite(lambda))) {
    stop("All thresholds in `lambda` must be finite.", call. = FALSE)
  }

  if (is.null(adj)) adj <- .adjacency_list(n, edges)

  lambda_min <- min(lambda)
  eligible <- which(is.finite(y) & y >= lambda_min)
  beta0 <- numeric(length(lambda))

  if (length(eligible) == 0L) return(beta0)

  # Exact ties are processed deterministically. At a threshold equal to a tied
  # value, all vertices in the tied block are included, and the final component
  # count after that block is invariant to the within-block order.
  ord <- eligible[order(-y[eligible], eligible)]

  dsu <- .make_dsu(n)
  active <- logical(n)
  beta_after <- integer(length(ord))
  cumulative_merging <- 0L

  for (r in seq_along(ord)) {
    j <- ord[r]
    js <- adj[[j]]
    active_neighbours <- if (length(js) > 0L) js[active[js]] else integer(0)

    q_j <- 0L
    if (length(active_neighbours) > 0L) {
      roots <- vapply(active_neighbours, dsu$find, integer(1))
      q_j <- length(unique(roots))
    }

    cumulative_merging <- cumulative_merging + q_j
    active[j] <- TRUE

    if (length(active_neighbours) > 0L) {
      for (h in active_neighbours) dsu$union(j, h)
    }

    beta_after[r] <- r - cumulative_merging
  }

  y_ord <- y[ord]
  for (g in seq_along(lambda)) {
    k <- sum(y_ord >= lambda[g])
    if (k > 0L) beta0[g] <- beta_after[k]
  }

  as.numeric(beta0)
}

.integrated_beta0_exact <- function(y, edges, adj = NULL) {
  thresholds <- sort(unique(y[is.finite(y) & y > 0]))
  if (length(thresholds) == 0L) return(0)

  beta <- .beta0_values(y, edges, thresholds, adj = adj)
  widths <- diff(c(0, thresholds))
  sum(widths * beta)
}

.integrated_aml <- function(y, edges, adj = NULL) {
  # Exact component-specific activation--merging allocation used in the paper.
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(adj)) adj <- .adjacency_list(n, edges)

  .check_distinct_positive(y)

  A <- pmax(y, 0)
  M <- numeric(n)
  q_entry <- integer(n)
  positive <- which(is.finite(y) & y > 0)

  if (length(positive) == 0L) {
    return(list(
      A = A,
      M = M,
      L = A,
      B = sum(A),
      q_entry = q_entry,
      allocation_rule = "component_specific_symmetric"
    ))
  }

  ord <- positive[order(y[positive], decreasing = TRUE)]
  dsu <- .make_dsu(n)
  active <- logical(n)

  for (j in ord) {
    js <- adj[[j]]
    higher_neighbours <- if (length(js) > 0L) js[active[js]] else integer(0)

    if (length(higher_neighbours) > 0L) {
      roots <- vapply(higher_neighbours, dsu$find, integer(1))
      roots_unique <- unique(roots)
      q_entry[j] <- length(roots_unique)

      for (rr in roots_unique) {
        nbrs_rr <- higher_neighbours[roots == rr]
        k_rr <- length(nbrs_rr)
        w <- y[j]

        # One merging unit is created for each previously active component C
        # touched by j. Half is assigned to j and half is shared equally among
        # the neighbours of j belonging to C. Multiplication by y[j] performs
        # the exact integration over thresholds 0 < lambda <= y[j].
        M[j] <- M[j] + 0.5 * w
        M[nbrs_rr] <- M[nbrs_rr] + 0.5 * w / k_rr
      }
    }

    active[j] <- TRUE
    if (length(higher_neighbours) > 0L) {
      for (h in higher_neighbours) dsu$union(j, h)
    }
  }

  L <- A - M

  list(
    A = A,
    M = M,
    L = L,
    B = sum(L),
    q_entry = q_entry,
    allocation_rule = "component_specific_symmetric"
  )
}

.trapz <- function(y, x) {
  if (length(y) != length(x)) {
    stop("`y` and `x` must have the same length.", call. = FALSE)
  }
  if (length(y) < 2L) return(0)
  sum(0.5 * (y[-1] + y[-length(y)]) * diff(x))
}
