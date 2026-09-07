# Internal signal utilities ------------------------------------------------

.prepare_signal <- function(x, reference = "mean", scale = "sd") {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (length(x) < 2L) stop("`x` must contain at least two areal observations.", call. = FALSE)
  if (any(!is.finite(x))) {
    stop("`x` must contain only finite values; handle missing values before analysis.", call. = FALSE)
  }

  if (is.numeric(reference) && length(reference) == 1L && is.finite(reference)) {
    b <- as.numeric(reference)
    reference_label <- "numeric"
  } else if (is.character(reference) && length(reference) == 1L) {
    reference <- match.arg(reference, c("mean", "median", "zero"))
    b <- switch(reference, mean = mean(x), median = stats::median(x), zero = 0)
    reference_label <- reference
  } else {
    stop("`reference` must be a finite number or one of 'mean', 'median', 'zero'.", call. = FALSE)
  }

  if (is.numeric(scale) && length(scale) == 1L && is.finite(scale) && scale > 0) {
    s <- as.numeric(scale)
    scale_label <- "numeric"
  } else if (is.character(scale) && length(scale) == 1L) {
    scale <- match.arg(scale, c("sd", "mad", "none"))
    s <- switch(scale, sd = stats::sd(x), mad = stats::mad(x), none = 1)
    scale_label <- scale
  } else {
    stop("`scale` must be a positive number or one of 'sd', 'mad', 'none'.", call. = FALSE)
  }

  if (!is.finite(s) || s <= 0) {
    stop(
      "The selected scale is zero or undefined. Use `scale = 'none'` or provide ",
      "a positive numeric scale.",
      call. = FALSE
    )
  }

  z <- (as.numeric(x) - b) / s

  list(
    x = as.numeric(x),
    z = z,
    reference = b,
    reference_method = reference_label,
    scale = s,
    scale_method = scale_label
  )
}

.side_vector <- function(z, side) {
  switch(
    side,
    superlevel = z,
    sublevel = -z,
    stop("Unknown analysis side.", call. = FALSE)
  )
}

.original_threshold <- function(lambda, b, s, side) {
  if (side == "superlevel") b + s * lambda else b - s * lambda
}
