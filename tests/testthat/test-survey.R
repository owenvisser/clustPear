test_that("full person-level NHANES-style survey design gives Taylor SEs", {
  skip_if_not_installed("survey")
  person <- expand.grid(stratum = 1:3, psu = 1:2, within = 1:2)
  person$SEQN <- seq_len(nrow(person))
  person$weight <- 1 + person$SEQN / 9
  design <- survey::svydesign(
    ids = ~psu, strata = ~stratum, weights = ~weight,
    data = person, nest = TRUE
  )
  # Last person is NOT part of the dental domain. Keep in design.
  teeth <- person[rep(seq_len(nrow(person) - 1L), each = 3L), c("SEQN", "weight")]
  teeth$X <- rep(c(-1, 1, 3), nrow(teeth) / 3) +
    (teeth$SEQN %% 4) / 10
  teeth$Y <- rep(c(0, 3, 1), nrow(teeth) / 3) +
    (teeth$SEQN %% 5) / 5

  old <- cluster_pearson(teeth, "SEQN", "X", "Y", K = 4,
                         er_max_pairs = 100, survey_weights = "weight")
  new <- cluster_pearson(teeth, "SEQN", "X", "Y", K = 4,
                         er_max_pairs = 100, survey_design = design)
  expect_equal(new$pearson, old$pearson, tolerance = 1e-10)
  expect_identical(new$variance_method, rep("survey_taylor", 4L))
  expect_equal(new$design_df, rep(survey::degf(design), 4L))
  expect_equal(new$se_independent, old$se, tolerance = 1e-10)
  expect_true(all(is.finite(new$se)))
  expect_true(all(new$se >= 0))
  expect_equal(new$lower, new$pearson -
                 stats::qt(.975, survey::degf(design)) * new$se)
  expect_true(any(abs(new$se - old$se) > 1e-8))

  # Independent verification of NW's ratio/gradient linearization.
  d <- teeth
  w <- d$weight
  W <- sum(w)
  m <- colSums(cbind(d$X, d$Y, d$X^2, d$Y^2, d$X * d$Y) * w) / W
  vx <- m[3] - m[1]^2
  vy <- m[4] - m[2]^2
  den <- sqrt(vx * vy)
  rho <- (m[5] - m[1] * m[2]) / den
  grad <- c(-m[2] / den + rho * m[1] / vx,
            -m[1] / den + rho * m[2] / vy,
            -rho / (2 * vx), -rho / (2 * vy), 1 / den)
  d$if_nw <- drop(cbind(d$X - m[1], d$Y - m[2],
                         d$X^2 - m[3], d$Y^2 - m[4],
                         d$X * d$Y - m[5]) %*% grad) / W
  # Each person has 3 teeth. This is an UNWEIGHTED score, because survey
  # applies the person-level weights when calculating PSU totals.
  person$score <- 0
  tab <- tapply(d$if_nw, d$SEQN, sum)
  person$score[match(names(tab), person$SEQN)] <- tab
  psu_score <- aggregate(person$weight * person$score,
                         by = person[c("stratum", "psu")], FUN = sum)
  manual_var <- sum(vapply(split(psu_score$x, psu_score$stratum),
                           function(s) (s[1] - s[2])^2, numeric(1)))
  expect_equal(new$se[new$method == "NW"]^2, manual_var,
               tolerance = 1e-8)
})

test_that("survey mode validates all sampled persons and design object", {
  skip_if_not_installed("survey")
  p <- data.frame(SEQN = 1:8, psu = rep(1:2, 4),
                  strata = rep(1:2, each = 4), w = 1 + (1:8) / 10)
  design <- survey::svydesign(ids = ~psu, strata = ~strata,
                              weights = ~w, data = p, nest = TRUE)
  teeth <- data.frame(SEQN = rep(1:7, each = 2),
                      X = rep(1:7, each = 2) + rep(c(-1, 1), 7),
                      Y = rep(2:8, each = 2) + rep(c(1, -1), 7))
  expect_error(cluster_pearson(teeth, "SEQN", "X", "Y",
                               survey_weights = rep(1, nrow(teeth)),
                               survey_design = design), "either survey_design")
  expect_error(cluster_pearson(teeth, "SEQN", "X", "Y",
                               survey_design = data.frame()), "svydesign")
  expect_error(cluster_pearson(teeth, "SEQN", "X", "Y",
                               survey_design = design, survey_id = "missing"),
               "survey_id")
  wrong <- teeth
  wrong$SEQN[1] <- 999
  expect_error(cluster_pearson(wrong, "SEQN", "X", "Y",
                               survey_design = design), "missing from survey_design")
  expect_equal(nrow(design$variables), 8L) # should not be mutated by fit
})

test_that("patient linearizations are derivatives of weighted NW CW EQ ER", {
  d <- data.frame(id = rep(1:7, each = 3),
                  X = rep(c(1, 2, 4), 7) + rep(1:7, each = 3) %% 3,
                  Y = rep(c(3, 1, 5), 7) + rep(1:7, each = 3) %% 2,
                  w = rep(1 + (1:7) / 5, each = 3))
  # Work directly with the internal C++ backend to test the derivative used
  # by survey::svytotal without relying on survey for the numerical gradient.
  f <- function(weights) clustPear:::cluster_cor_weighted_cpp(
    d$id, d$X, d$Y, rep(weights, each = 3), 4L, 100L,
    c("NW", "CW", "EQ", "ER")
  )
  w <- unique(d$w)
  z <- f(w)
  expect_equal(dim(z$linearized), c(7L, 4L))
  expect_equal(drop(crossprod(w, z$linearized)), rep(0, 4),
               tolerance = 1e-9)
  eps <- 1e-5
  for (i in seq_along(w)) {
    w2 <- w
    w2[i] <- w2[i] + eps
    num <- (f(w2)$estimate - z$estimate) / eps
    expect_equal(as.numeric(num), as.numeric(z$linearized[i, ]),
                 tolerance = 1e-4)
  }
})
