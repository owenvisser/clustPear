# clustPear (working name)

## Installation

The development version of `clustPear` is available
from GitHub:

```r
install.packages("remotes")
remotes::install_github("YOUR_USERNAME/clustPear")
```

## Example

```r
library(clustPear)

d <- simulate_clustered()

cluster_pearson(
  d,
  cluster = "cluster_id",
  x = "X",
  y = "Y",
  method = "all"
)
```

The package provides NW, CW, EQ, and ER estimators
of marginal Pearson association for clustered data.

## Further Description

Minimal clustered Pearson correlations, implemented with the V11 C++ estimator core.

```r
# Extract clustPear source into a local directory.
# Install Rcpp and remotes if needed:
# install.packages(c("Rcpp", "tibble", "remotes"))
# remotes::install_local("path/to/clustPear", dependencies = TRUE)
library(clustPear)
d <- simulate_clustered()   # one observed data frame, 100 clusters / 20 potential units
d                     # tibble prints a compact preview
head(d)
cluster_pearson(d, cluster="cluster_id", x="X", y="Y")

# Optional cluster-constant sampling weights: column or row-aligned vector
set.seed(10)
cluster_w <- setNames(runif(length(unique(d$cluster_id)), 0.5, 2),
  unique(d$cluster_id))
d$weight <- unname(cluster_w[as.character(d$cluster_id)])
cluster_pearson(d, "cluster_id", "X", "Y", survey_weights="weight")
cluster_pearson(d, "cluster_id", "X", "Y", survey_weights=d$weight)

# Binary or ordinal category examples:
d$X_binary <- factor(ifelse(d$X > 0, "above", "below"), levels=c("below","above"))
cluster_pearson(d, "cluster_id", "X_binary", "Y", x_type="categorical")
```

The simulated `d` is a **tibble**, so printing it shows a compact preview rather than every observed row. No observations are discarded; use `nrow(d)` to see the full size. The correlation result remains a plain data frame with one row per method and columns:
`method`, `pearson`, `se`, `lower`, `upper`, `n_clusters`, `n_obs`, and `er_exact`.

## Statistical interpretation

- `NW`: ordinary pooled Pearson correlation.
- `CW`: Pearson correlation with equal total contribution from each cluster.
- `EQ`: average of OPW Pearson correlations over a selected distinct-rank cutpoint grid (`K=20` default).
- `ER`: average over all distinct-rank cutpoint pairs if within `er_max_pairs`; otherwise deterministic capped rank-grid approximation.
- All SEs: independent-cluster delta-method sandwich. EQ and ER treat empirical cutpoints as fixed. CI: `estimate +/- qnorm((1+level)/2) * se`, not a survey-design interval.
- A categorical variable is converted to numeric scores. Binary categorical is valid (phi when both binary). For more than two unordered categories, Pearson depends on arbitrary ordering and should be interpreted only if scores/order are meaningful.
- `cluster`, `x`, and `y` are the only required input columns (plus `survey_weights` if supplied). All other columns are ignored. Rows with missing IDs or missing X/Y are dropped; singleton clusters are retained. At least two clusters and nonconstant X/Y are required.
- ER approximation status is in `er_exact`; `FALSE` means numerical ER reference rather than exhaustive averaging.
- Optional `survey_weights` adjusts the estimating equations using **one constant positive sampling weight per cluster**; it does not multiply weights into the outcomes. Supply a column name or row-aligned numeric vector.
- Weighted SEs still use a **patient/cluster-independent sandwich**; this is **not** a stratified-PSU NHANES design variance or a publication-ready NHANES confidence interval. Apply survey replication (e.g. a justified PSU jackknife) in a separate analysis script.
- The package does not apply dental-specific recoding or automatically implement complex-survey design inference.

## Before publication

This is a **starter package**, not a vetted CRAN release: R itself is not available in the creation environment. If you modify Rcpp-exported function signatures, run `Rcpp::compileAttributes()`, install, `testthat::test_local()`, and `R CMD check --as-cran` on a Windows R installation. Confirm numerical parity with published V11 reference results and validate singleton handling; V11 simulations excluded clusters with fewer than two observations, but this generic interface intentionally retains singleton clusters.

Replace the placeholder maintainer/License details before submission.
