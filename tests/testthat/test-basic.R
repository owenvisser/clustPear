test_that("simulation returns a reproducible observed three-column dataset", {
  a <- simulate_clustered(M = 15, nmax = 8, eta = 2, seed = 42)
  b <- simulate_clustered(M = 15, nmax = 8, eta = 2, seed = 42)

  expect_identical(a, b)
  expect_s3_class(a, "tbl_df")
  expect_s3_class(a, "data.frame")
  expect_equal(ncol(a), 3L)
  expect_named(a, c("cluster_id", "X", "Y"))
  expect_true(nrow(a) <= 15 * 8)
  expect_false(anyNA(a))
})

test_that("four estimators return a simple data frame", {
  a <- simulate_clustered(M = 20, nmax = 9, eta = 0, seed = 19)
  z <- cluster_pearson(a, "cluster_id", "X", "Y", K = 5)

  expect_s3_class(z, "data.frame")
  expect_equal(z$method, c("NW", "CW", "EQ", "ER"))
  expect_true(all(is.finite(z$pearson)))
  expect_true(all(is.finite(z$se)))
  expect_true(all(z$lower <= z$upper))
  expect_true(z$er_exact[4])
  expect_equal(z$pearson[1], stats::cor(a$X, a$Y), tolerance = 1e-10)
})

test_that("other columns and cluster ID representation do not affect results", {
  a <- simulate_clustered(M = 22, nmax = 6, eta = 0, seed = 111)
  p <- cluster_pearson(a, "cluster_id", "X", "Y", K = 4)

  a$irrelevant <- rnorm(nrow(a))
  a$cluster_id <- paste0("participant-", a$cluster_id)

  p2 <- cluster_pearson(a, "cluster_id", "X", "Y", K = 4)

  expect_equal(p$pearson, p2$pearson, tolerance = 1e-10)
})

test_that("binary categorical works and EQ agrees with exhaustive ER", {
  a <- simulate_clustered(M = 25, nmax = 8, eta = 0, seed = 112)

  a$X <- factor(
    ifelse(a$X > 0, "high", "low"),
    levels = c("low", "high")
  )

  a$Y <- factor(
    ifelse(a$Y > 0, "high", "low"),
    levels = c("low", "high")
  )

  z <- cluster_pearson(
    a, "cluster_id", "X", "Y",
    x_type = "categorical",
    y_type = "categorical",
    K = 20
  )

  expect_equal(z$pearson[3], z$pearson[4], tolerance = 1e-10)
  expect_true(z$er_exact[4])
})

test_that("singletons are not silently discarded", {
  a <- data.frame(
    id = c(1, 2, 3, 3, 4, 4),
    X = c(1, 2, 3, 4, 5, 6),
    Y = c(4, 6, 5, 8, 7, 9)
  )

  z <- cluster_pearson(a, "id", "X", "Y", method = "CW")

  expect_equal(z$n_clusters, 4)
  expect_equal(z$n_obs, 6)
})

test_that("missing input is dropped only from named columns", {
  a <- simulate_clustered(M = 15, nmax = 8, eta = 0, seed = 222)

  a$ignored <- NA_real_
  a$X[1] <- NA_real_

  z <- cluster_pearson(a, "cluster_id", "X", "Y", method = "NW")

  expect_equal(z$n_obs, nrow(a) - 1)
})
