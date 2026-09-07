# Print methods ------------------------------------------------------------

#' @noRd
#' @export
print.spatconnect_global <- function(x, ...) {
  cat("Global spatial-connectivity analysis\n")
  cat("  Analysis:    ", x$side, "\n", sep = "")
  cat("  Areas/edges: ", x$n_areas, "/", x$n_edges, "\n", sep = "")
  cat("  Reference:   ", format(x$reference, digits = 5), "\n", sep = "")
  cat("  Scale:       ", format(x$scale, digits = 5), "\n", sep = "")
  cat("  B:           ", format(x$integrated_B, digits = 5), "\n", sep = "")
  cat("  D:           ", format(x$discrepancy, digits = 5), "\n", sep = "")
  cat("  Monte Carlo p-value: ", format.pval(x$p_value, digits = 4), "\n", sep = "")
  cat("  Relabelings: ", x$n_perm, "\n", sep = "")
  cat("  Integration: ", x$integration_method, "\n", sep = "")
  invisible(x)
}

#' @noRd
#' @export
print.spatconnect_global_both <- function(x, ...) {
  cat("Global spatial-connectivity analysis\n\n")
  cat("Superlevel analysis:\n")
  print(x$superlevel)
  cat("\nSublevel analysis:\n")
  print(x$sublevel)
  invisible(x)
}

#' @noRd
#' @export
print.spatconnect_local <- function(x, ...) {
  n_unadj <- sum(x$results$p_value < x$alpha, na.rm = TRUE)
  n_adj <- sum(x$results$p_adjusted < x$alpha, na.rm = TRUE)

  cat("Local spatial-connectivity analysis\n")
  cat("  Analysis:    ", x$side, "\n", sep = "")
  cat("  Areas/edges: ", x$n_areas, "/", x$n_edges, "\n", sep = "")
  cat("  Allocation:  component-specific symmetric\n")
  cat("  Relabelings per area: ", x$n_perm, "\n", sep = "")
  cat("  Unadjusted p < ", x$alpha, ": ", n_unadj, "\n", sep = "")
  cat("  Adjusted p < ", x$alpha, ": ", n_adj, " (", x$adjust, ")\n", sep = "")
  invisible(x)
}

#' @noRd
#' @export
print.spatconnect_local_both <- function(x, ...) {
  cat("Local spatial-connectivity analysis\n\n")
  cat("Superlevel analysis:\n")
  print(x$superlevel)
  cat("\nSublevel analysis:\n")
  print(x$sublevel)
  invisible(x)
}
