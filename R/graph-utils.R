# Internal graph utilities -------------------------------------------------

.as_edge_list <- function(graph, n) {
  if (inherits(graph, "listw")) {
    graph <- graph$neighbours
  }

  if (inherits(graph, "nb")) {
    edge_list <- vector("list", length(graph))
    counter <- 0L

    for (i in seq_along(graph)) {
      js <- as.integer(graph[[i]])
      js <- js[is.finite(js) & js > i]
      if (length(js) > 0L) {
        counter <- counter + 1L
        edge_list[[counter]] <- data.frame(i = i, j = js)
      }
    }

    if (counter == 0L) {
      return(data.frame(i = integer(0), j = integer(0)))
    }

    edges <- do.call(rbind, edge_list[seq_len(counter)])
    rownames(edges) <- NULL
    return(.validate_edges(edges, n))
  }

  if (is.matrix(graph) && nrow(graph) == n && ncol(graph) == n) {
    nz <- graph != 0
    if (!identical(nz, t(nz))) {
      stop("Adjacency matrices must represent an undirected (symmetric) graph.",
           call. = FALSE)
    }
    idx <- which(nz & upper.tri(nz), arr.ind = TRUE)
    if (nrow(idx) == 0L) {
      return(data.frame(i = integer(0), j = integer(0)))
    }
    edges <- data.frame(i = idx[, 1], j = idx[, 2])
    return(.validate_edges(edges, n))
  }

  if (is.data.frame(graph) || is.matrix(graph)) {
    graph <- as.data.frame(graph)
    if (ncol(graph) < 2L) {
      stop("An edge list must contain at least two columns.", call. = FALSE)
    }

    if (all(c("i", "j") %in% names(graph))) {
      edges <- graph[, c("i", "j"), drop = FALSE]
    } else if (all(c("from", "to") %in% names(graph))) {
      edges <- graph[, c("from", "to"), drop = FALSE]
      names(edges) <- c("i", "j")
    } else {
      edges <- graph[, 1:2, drop = FALSE]
      names(edges) <- c("i", "j")
    }

    return(.validate_edges(edges, n))
  }

  stop(
    "`graph` must be an spdep `nb`/`listw` object, a two-column edge list, ",
    "or an n x n adjacency matrix.",
    call. = FALSE
  )
}

.validate_edges <- function(edges, n) {
  if (nrow(edges) == 0L) {
    return(data.frame(i = integer(0), j = integer(0)))
  }

  i <- suppressWarnings(as.integer(edges[[1]]))
  j <- suppressWarnings(as.integer(edges[[2]]))

  if (anyNA(i) || anyNA(j)) {
    stop("Edge endpoints must be integer vertex indices.", call. = FALSE)
  }
  if (any(i < 1L | i > n | j < 1L | j > n)) {
    stop("Edge endpoints must lie between 1 and length(x).", call. = FALSE)
  }
  if (any(i == j)) {
    stop("Self-edges are not allowed.", call. = FALSE)
  }

  out <- data.frame(i = pmin(i, j), j = pmax(i, j))
  out <- unique(out)
  out <- out[order(out$i, out$j), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.area_ids <- function(x, graph) {
  if (!is.null(names(x)) && length(names(x)) == length(x) &&
      all(nzchar(names(x)))) {
    return(as.character(names(x)))
  }

  nb <- if (inherits(graph, "listw")) graph$neighbours else graph
  if (inherits(nb, "nb")) {
    ids <- attr(nb, "region.id")
    if (!is.null(ids) && length(ids) == length(x)) {
      return(as.character(ids))
    }
  }

  as.character(seq_along(x))
}
