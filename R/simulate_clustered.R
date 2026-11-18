#' Simulate an observed clustered data frame
#'
#' Generates one clustered bivariate sample using the V11 simulation model.
#' Only retained observations are returned; no missing/removed rows are
#' included. The output columns can be passed directly to cluster_pearson().
#'
#' @param M Number of clusters (default 100).
#' @param nmax Potential observations per cluster (default 20).
#' @param rho_uv Correlation between cluster effects (default 0.25).
#' @param eta Retention intensity (default 2); 0 retains every unit.
#' @param sigma2,tau2 Variance parameters from the simulation model.
#' @param seed Deterministic random seed.
#' @return A tibble with exactly three columns: cluster_id, X, Y.
#'   Printing the tibble shows only the first few rows; all observed rows remain available.
#' @export
#' @examples
#' d <- simulate_clustered()
#' head(d)
#' cluster_pearson(d, "cluster_id", "X", "Y")

simulate_clustered <- function(
  M = 100L,
  nmax = 20L,
  rho_uv = 0.25,
  eta = 2,
  sigma2 = 0.5,
  tau2 = 0.5,
  seed = 20261114L
) {
  if (
    length(M) != 1L || !is.numeric(M) || !is.finite(M) ||
    M != as.integer(M) || M < 2
  )
    stop("M must be an integer >= 2", call. = FALSE)

  if (
    length(nmax) != 1L || !is.numeric(nmax) || !is.finite(nmax) ||
    nmax != as.integer(nmax) || nmax < 1
  )
    stop("nmax must be an integer >= 1", call. = FALSE)

  tibble::as_tibble(simulate_clustered_cpp(
    as.integer(M),
    as.integer(nmax),
    rho_uv,
    eta,
    sigma2,
    tau2,
    seed
  ))
}
