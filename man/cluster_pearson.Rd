\name{cluster_pearson}
\alias{cluster_pearson}
\title{Cluster-weighted Pearson correlations}
\usage{
cluster_pearson(data, cluster, x, y, x_type = "continuous",
  y_type = "continuous", method = "all", K = 20L,
  er_max_pairs = 100000L, conf_level = 0.95,
  survey_weights = NULL, survey_design = NULL, survey_id = cluster)
}
\arguments{
\item{data}{Data frame. Other columns are ignored.}
\item{cluster, x, y}{Column names supplied as strings.}
\item{x_type, y_type}{"continuous" (default) or "categorical".}
\item{method}{"all" (default), "NW", "CW", "EQ" or "ER".}
\item{K}{EQ quantile-grid resolution (default 20).}
\item{er_max_pairs}{ER pair budget; a larger exhaustive grid is replaced with a deterministic rank-grid approximation.}
\item{conf_level}{Nominal confidence level (default 0.95).}
\item{survey_weights}{Optional positive, cluster-constant sampling weights specified as a column name or numeric vector with one value per input row. For point estimates and independent-patient sandwich SEs only. Not used together with \code{survey_design}.}
\item{survey_design}{Optional full person-level \code{survey::svydesign()} object built before dental and age restrictions, containing one record for each sampled person with a positive MEC sampling weight. Supplies survey weights and stratified-PSU Taylor-linearized variance. Requires the suggested \code{survey} package.}
\item{survey_id}{Name of the person identifier within \code{survey_design}, defaults to \code{cluster} (e.g. SEQN).}
}
\description{Computes marginal Pearson correlations for paired clustered outcomes using ordinary pooled weights (NW), inverse-cluster-size weights (CW), empirical quantile-grid averaged OPW weights (EQ), and empirical rank averaged OPW weights (ER).}
\details{Complete observed pairs are used; singleton clusters contribute. The returned standard errors are independent-cluster sandwich/delta estimates. EQ and ER use fixed observed cutpoints in the variance calculation, so their confidence bounds do not account for cutpoint-selection variability. Numerical ER is labeled by er_exact; if FALSE the computed result is an approximation. When \code{survey_weights} is supplied, it weights cluster estimating-equation contributions (not X/Y values) and must be constant within each cluster. These resulting SEs treat patients as independent clusters. When \code{survey_design} is supplied instead, the variance uses the full sampling design with patient-level estimating-equation linearizations, and confidence intervals use the design degrees of freedom. Persons outside the analysis domain contribute zero to the survey-total linearization. The empirical cutpoints in EQ/ER are held fixed in this Taylor approximation; this does not account for cutpoint estimation uncertainty. In categorical mode, factors are assigned ordered integer scores from their level ordering and character vectors from alphabetical ordering. Pearson correlation depends on scores for outcomes with more than two levels. Binary/binary Pearson equals phi correlation.}
\value{A data frame with columns method, pearson, se, lower, upper, n_clusters, n_obs and er_exact. With \code{survey_design}, also contains variance_method, design_df, and the old independent-patient sandwich SE as se_independent for comparison.}
\examples{
d <- simulate_clustered()
cluster_pearson(d, "cluster_id", "X", "Y")
}
