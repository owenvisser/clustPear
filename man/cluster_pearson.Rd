\name{cluster_pearson}
\alias{cluster_pearson}
\title{Cluster-weighted Pearson correlations}
\usage{
cluster_pearson(data, cluster, x, y, x_type = "continuous",
  y_type = "continuous", method = "all", K = 20L,
  er_max_pairs = 100000L, conf_level = 0.95, survey_weights = NULL)
}
\arguments{
\item{data}{Data frame. Other columns are ignored.}
\item{cluster, x, y}{Column names supplied as strings.}
\item{x_type, y_type}{"continuous" (default) or "categorical".}
\item{method}{"all" (default), "NW", "CW", "EQ" or "ER".}
\item{K}{EQ quantile-grid resolution (default 20).}
\item{er_max_pairs}{ER pair budget; a larger exhaustive grid is replaced with a deterministic rank-grid approximation.}
\item{conf_level}{Nominal confidence level (default 0.95).}
\item{survey_weights}{Optional positive, cluster-constant sampling weights specified as a column name or numeric vector with one value per input row. Defaults to NULL (unweighted).}
}
\description{Computes marginal Pearson correlations for paired clustered outcomes using ordinary pooled weights (NW), inverse-cluster-size weights (CW), empirical quantile-grid averaged OPW weights (EQ), and empirical rank averaged OPW weights (ER).}
\details{Complete observed pairs are used; singleton clusters contribute. The returned standard errors are independent-cluster sandwich/delta estimates. EQ and ER use fixed observed cutpoints in the variance calculation, so their confidence bounds do not account for cutpoint-selection variability. Numerical ER is labeled by er_exact; if FALSE the computed result is an approximation. When \code{survey_weights} is supplied, it weights cluster estimating-equation contributions (not X/Y values) and must be constant within each cluster. Resulting SEs treat patients as independent clusters; they do not account for survey strata, PSU replication, or design-based degrees of freedom. In categorical mode, factors are assigned ordered integer scores from their level ordering and character vectors from alphabetical ordering. Pearson correlation depends on scores for outcomes with more than two levels. Binary/binary Pearson equals phi correlation.}
\value{A data frame with columns method, pearson, se, lower, upper, n_clusters, n_obs and er_exact.}
\examples{
d <- simulate_clustered()
cluster_pearson(d, "cluster_id", "X", "Y")
}
