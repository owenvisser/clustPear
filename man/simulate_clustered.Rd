\name{simulate_clustered}
\alias{simulate_clustered}
\title{Simulate observed clustered paired outcomes}
\usage{
simulate_clustered(M = 100L, nmax = 20L, rho_uv = 0.25,
  eta = 2, sigma2 = 0.5, tau2 = 0.5, seed = 20261114L)
}
\arguments{
\item{M}{Number of clusters.}
\item{nmax}{Potential units per cluster.}
\item{rho_uv}{Correlation of cluster effects.}
\item{eta}{Outcome-dependent retention parameter.}
\item{sigma2, tau2}{Variance parameters of the simulation model.}
\item{seed}{Nonnegative deterministic seed.}
}
\description{Generates one sample under the same clustered bivariate and outcome-dependent retention mechanism used in the accompanying methodology study. Only retained observations are returned.}
\value{A tibble containing cluster_id, X and Y, with no unobserved rows. Printing displays a limited preview without truncating the underlying data.}
\examples{
d <- simulate_clustered()
d                 # tibble preview; not all rows printed
head(d)
}
