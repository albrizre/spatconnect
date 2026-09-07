test_that("component-specific integrated decomposition is exact on a chain", {
  y <- c(3, 2, 1)
  edges <- data.frame(i = c(1, 2), j = c(2, 3))

  aml <- spatconnect::local_connectivity(
    y, edges,
    reference = 0,
    scale = "none",
    side = "superlevel",
    n_perm = 9,
    seed = 1,
    classify = FALSE,
    progress = FALSE
  )

  expect_equal(aml$results$A, c(3, 2, 1))
  expect_equal(aml$results$M, c(1, 1.5, 0.5))
  expect_equal(aml$results$L, c(2, 0.5, 0.5))
  expect_equal(sum(aml$results$L), 3)
  expect_equal(aml$integrated_B, 3)
})

test_that("allocation is normalized separately within an active component", {
  y <- c(4, 3, 2)
  edges <- data.frame(i = c(1, 1, 2), j = c(2, 3, 3))

  l <- local_connectivity(
    y, edges,
    reference = 0,
    scale = "none",
    side = "superlevel",
    n_perm = 9,
    seed = 1,
    classify = FALSE,
    progress = FALSE
  )

  # At y_3 = 2, vertices 1 and 2 already belong to the same active component,
  # so the receiving half is shared equally between them.
  expect_equal(l$results$M, c(2, 2, 1))
  expect_equal(l$results$q_entry, c(0, 1, 1))
  expect_equal(sum(l$results$L), l$integrated_B)
  expect_equal(l$integrated_B, 4)
})

test_that("an entering area can merge two previously disconnected components", {
  y <- c(4, 3, 2)
  edges <- data.frame(i = c(1, 2), j = c(3, 3))

  l <- local_connectivity(
    y, edges,
    reference = 0,
    scale = "none",
    side = "superlevel",
    n_perm = 9,
    seed = 1,
    classify = FALSE,
    progress = FALSE
  )

  expect_equal(l$results$q_entry, c(0, 0, 2))
  expect_equal(l$results$M, c(1, 1, 2))
  expect_equal(l$integrated_B, 5)
})

test_that("global analysis handles a selected side with no positive values", {
  x <- c(-3, -2, -1)
  edges <- data.frame(i = c(1, 2), j = c(2, 3))

  g <- global_connectivity(
    x, edges,
    reference = 0,
    scale = "none",
    side = "superlevel",
    n_perm = 9,
    seed = 1
  )

  expect_equal(g$integrated_B, 0)
  expect_equal(g$discrepancy, 0)
  expect_equal(g$p_value, 1)
})

test_that("both-sided API returns superlevel and sublevel analyses", {
  x <- c(-2, -1, 1, 3)
  edges <- data.frame(i = c(1, 2, 3), j = c(2, 3, 4))

  g <- global_connectivity(
    x, edges,
    reference = 0,
    scale = "none",
    side = "both",
    n_perm = 9,
    seed = 1
  )

  expect_s3_class(g, "spatconnect_global_both")
  expect_s3_class(g$superlevel, "spatconnect_global")
  expect_s3_class(g$sublevel, "spatconnect_global")
})

test_that("global common-mean calibration returns R plus one discrepancies", {
  x <- c(-2, -1, 1, 3)
  edges <- data.frame(i = c(1, 2, 3), j = c(2, 3, 4))

  g <- global_connectivity(
    x, edges,
    reference = 0,
    scale = "none",
    side = "superlevel",
    n_perm = 19,
    seed = 2
  )

  expect_equal(length(g$discrepancy_all), 20)
  expect_identical(g$global_calibration, "common_mean_R_plus_1")
  expect_equal(g$p_value, mean(g$discrepancy_all >= g$discrepancy))
})

test_that("fixed grid option is available for case-study reproduction", {
  x <- c(-2, -1, 1, 3)
  edges <- data.frame(i = c(1, 2, 3), j = c(2, 3, 4))

  g <- global_connectivity(
    x, edges,
    reference = 0,
    scale = "none",
    side = "superlevel",
    n_perm = 9,
    seed = 1,
    n_grid = 80
  )

  expect_equal(nrow(g$curve), 80)
  expect_identical(g$integration_method, "trapezoidal_regular_grid_80")
})

test_that("exact positive vertex ties are rejected by the local method", {
  x <- c(2, 2, 1)
  edges <- data.frame(i = c(1, 2), j = c(2, 3))

  expect_error(
    local_connectivity(
      x, edges,
      reference = 0,
      scale = "none",
      side = "superlevel",
      n_perm = 9,
      seed = 1,
      classify = FALSE,
      progress = FALSE
    ),
    "exactly tied positive selected values"
  )
})
