# clustPear (working name)

## Installation

Install the local development source from the `clustPear/` package directory:

```r
install.packages(c("Rcpp", "tibble", "remotes"))
remotes::install_local("path/to/clustPear", dependencies = TRUE)
```

The repository URL and maintainer details are placeholders until the package
is published; do not use an unconfigured GitHub installation link.

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

## NHANES complex-survey inference (new in 0.1.3)

The `survey_weights` argument **only** provides unequal-weight point estimates
and independent-patient sandwich SEs. **Survey weights by themselves do not
identify sampling strata and PSUs.** For correct NHANES-style design-based
variance, supply a `survey_design` created on the *full MEC person sample* before
age, dental eligibility, or complete-pair restrictions.

```r
library(survey)
library(clustPear)

# person: one record per SEQN in both 2011–12 and 2013–14 DEMO files;
# include all positive-weight MEC examinees, not only dental participants.
person$WTMEC4YR <- person$WTMEC2YR / 2
person <- subset(person, is.finite(WTMEC4YR) & WTMEC4YR > 0 &
                        !is.na(SDMVSTRA) & !is.na(SDMVPSU))

design <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC4YR,
  data = person, nest = TRUE
)

# teeth: eligible tooth-level data, with any age restriction already applied.
fit <- cluster_pearson(
  teeth, cluster = "SEQN", x = "CAL", y = "DFS",
  method = "all", K = 20, survey_design = design
)
```

`survey_design` provides the weights automatically, so **do not also pass
`survey_weights`**. The result includes `variance_method = "survey_taylor"`
and `design_df`, plus `se_independent` for comparison. Confidence intervals use a Student-t reference distribution
with the survey design degrees of freedom. The design remains intact, with
zero contributions for persons outside the eligible dental/age domain.

NW uses the appropriate weighted ratio linearization for varying tooth counts.
CW, EQ, and ER use patient-normalized moment contributions; EQ/ER average
patient influence contributions *before* computing the design variance so the
covariance of dichotomization-specific correlations is included.

**Limitation:** Taylor linearization treats the empirically selected EQ/ER
cutpoints as fixed. The method accounts for sampling weights, strata, PSUs and
within-mouth contributions, but not for additional variation from learning
cutpoints on a new sample. Use a fully recalculated PSU replication analysis
as a sensitivity analysis if that uncertainty matters. It does not recover
the distribution of already lost teeth.


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
- Default or `survey_weights` SEs: independent-cluster delta-method sandwich and normal CIs. When `survey_design` is used: stratified-PSU Taylor-linearized SEs and survey-t CIs. Both modes treat empirical EQ/ER cutpoints as fixed.
- A categorical variable is converted to numeric scores. Binary categorical is valid (phi when both binary). For more than two unordered categories, Pearson depends on arbitrary ordering and should be interpreted only if scores/order are meaningful.
- `cluster`, `x`, and `y` are the only required input columns (plus `survey_weights` if supplied). All other columns are ignored. Rows with missing IDs or missing X/Y are dropped; singleton clusters are retained. At least two clusters and nonconstant X/Y are required.
- ER approximation status is in `er_exact`; `FALSE` means numerical ER reference rather than exhaustive averaging.
- Optional `survey_weights` adjusts the estimating equations using **one constant positive sampling weight per cluster**; it does not multiply weights into the outcomes. Supply a column name or row-aligned numeric vector.
- With `survey_weights` alone, SEs use the **independent-patient sandwich**. Pass `survey_design` for a stratified-PSU Taylor variance and survey t confidence interval. They cannot be passed together.
- The package does not apply dental-specific recoding. Complex-survey inference requires an explicitly supplied full person-level `survey_design` and the optional `survey` package.

## Before publication

This is a **starter package**, not a vetted CRAN release: R itself is not available in the creation environment. If you modify Rcpp-exported function signatures, run `Rcpp::compileAttributes()`, install, `testthat::test_local()`, and `R CMD check --as-cran` on a Windows R installation. Confirm numerical parity with published V11 reference results and validate singleton handling; V11 simulations excluded clusters with fewer than two observations, but this generic interface intentionally retains singleton clusters.

Replace the placeholder maintainer/License details before submission.
