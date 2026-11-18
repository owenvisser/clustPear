test_that("optional sampling weights preserve unweighted behavior for unit weights", {
  d <- simulate_clustered(M = 16, nmax = 9, eta = 1, seed = 607)
  d$w <- 1

  a <- cluster_pearson(d, "cluster_id", "X", "Y", K = 5, er_max_pairs = 256)
  b <- cluster_pearson(
    d, "cluster_id", "X", "Y",
    K = 5,
    er_max_pairs = 256,
    survey_weights = "w"
  )

  expect_equal(a$pearson, b$pearson, tolerance = 1e-9)
  expect_equal(a$se, b$se, tolerance = 1e-9)
  expect_identical(a$er_exact, b$er_exact)
})

test_that("sampling-weight column and vector match; scale invariant", {
  d <- simulate_clustered(M = 12, nmax = 7, eta = 0, seed = 22)
  d$w <- 0.4 + 0.3 * (d$cluster_id %% 5)

  by_name <- cluster_pearson(
    d, "cluster_id", "X", "Y",
    K = 4,
    er_max_pairs = 100,
    survey_weights = "w"
  )

  by_vec <- cluster_pearson(
    d, "cluster_id", "X", "Y",
    K = 4,
    er_max_pairs = 100,
    survey_weights = d$w
  )

  by_scale <- cluster_pearson(
    d, "cluster_id", "X", "Y",
    K = 4,
    er_max_pairs = 100,
    survey_weights = 100 * d$w
  )

  expect_identical(by_name, by_vec)
  expect_equal(by_vec$pearson, by_scale$pearson, tolerance = 1e-9)
  expect_equal(by_vec$se, by_scale$se, tolerance = 1e-9)
  expect_true(any(abs(
    by_vec$pearson - cluster_pearson(
      d, "cluster_id", "X", "Y",
      K = 4,
      er_max_pairs = 100
    )$pearson
  ) > 1e-7))
})

test_that("weighted NW/CW agree with direct weighted Pearson", {
  d <- data.frame(
    id = c(1, 1, 2, 2, 2, 3, 3),
    X = c(0, 1, 2, 3, 0, 1, 4),
    Y = c(1, 3, 1, 4, 2, 5, 2),
    w = c(1, 1, 2, 2, 2, 4, 4)
  )

  z <- cluster_pearson(d, "id", "X", "Y", method = "NW", survey_weights = "w")

  weighted_cor <- function(x, y, wt) {
    mu_x <- sum(wt * x) / sum(wt);
    mu_y <- sum(wt * y) / sum(wt)
    cov_xy <- sum(wt * (x - mu_x) * (y - mu_y)) / sum(wt)
    var_x <- sum(wt * (x - mu_x)^2) / sum(wt)
    var_y <- sum(wt * (y - mu_y)^2) / sum(wt)
    cov_xy / sqrt(var_x * var_y)
  }

  expect_equal(z$pearson, weighted_cor(d$X, d$Y, d$w), tolerance = 1e-10)

  cw <- cluster_pearson(d, "id", "X", "Y", method = "CW", survey_weights = d$w)
  wt <- d$w / as.numeric(table(d$id)[as.character(d$id)])

  expect_equal(cw$pearson, weighted_cor(d$X, d$Y, wt), tolerance = 1e-10)
})

test_that("weights are cluster-constant, finite and positive", {
  d <- simulate_clustered(M = 10, nmax = 6, eta = 0, seed = 91)
  d$w <- rep(1, nrow(d))
  d$w[1] <- 2

  expect_error(
    cluster_pearson(d, "cluster_id", "X", "Y", survey_weights = "w"),
    "constant within each cluster"
  )

  d$w <- rep(1, nrow(d));
  d$w[1] <- -1

  expect_error(
    cluster_pearson(d, "cluster_id", "X", "Y", survey_weights = "w"),
    "positive and finite"
  )

  expect_error(
    cluster_pearson(d, "cluster_id", "X", "Y", survey_weights = rep(1, 2)),
    "one value per row"
  )

  expect_error(
    cluster_pearson(d, "cluster_id", "X", "Y", survey_weights = "absent"),
    "column not found"
  )
})

test_that("weighted EQ equals exhaustive ER for a sufficiently fine grid", {
  d <- simulate_clustered(M = 15, nmax = 6, eta = 0, seed = 113)

  d$X <- as.integer(cut(d$X, breaks = c(-Inf, -0.5, 0, 0.5, Inf)))
  d$Y <- as.integer(cut(d$Y, breaks = c(-Inf, -0.8, 0.4, Inf)))
  d$w <- 1 + (d$cluster_id %% 4) * 0.2

  z <- cluster_pearson(
    d, "cluster_id", "X", "Y",
    K = 20,
    er_max_pairs = 100,
    survey_weights = "w"
  )

  expect_true(z$er_exact[4])
  expect_equal(z$pearson[3], z$pearson[4], tolerance = 1e-10)
  expect_equal(z$se[3], z$se[4], tolerance = 1e-10)
})
