# Example: use the generic clustPear package on the PREPARED NHANES tooth data.
# Run from the existing OPW project's root, not the clustPear package directory.
# This is an exploratory, unweighted patient-cluster calculation. It is NOT
# a replacement for the WTMEC4YR-weighted stratified-PSU jackknife NHANES
# manuscript analysis in scripts/17_nhanes_estimates_95CI.R.

stopifnot(requireNamespace("here", quietly = TRUE))
stopifnot(requireNamespace("clustPear", quietly = TRUE))

source(here::here("R", "nhanes_clinical_v3.R"))

# The existing preparation determines dental scoring, ages and eligibility.

d <- nhanes_prepare(nhanes_read_teeth(), require_design = FALSE)

# Change these four settings to test other endpoints/domains/methods.

age_group <- "All"  # one of "All", "30-49", "50-64", "65+"
x <- "CAL"        # or "PD"
y <- "DFS"        # or "DS", "DSI", "FS", "FSI", "DFSI"
K <- 20L

z <- nhanes_age_groups(d)[[age_group]]

if (is.null(z)) stop("Unknown NHANES age group: ", age_group)
if (!x %in% names(z) || !y %in% names(z)) stop("Unknown X/Y NHANES columns")

# Only the three named columns are read by clustPear. Complete pairs are
# filtered internally; patient clusters with one complete pair are retained.

fit <- clustPear::cluster_pearson(
  data = z,
  cluster = "SEQN",
  x = x,
  y = y,
  K = K,
  er_max_pairs = 8192L
)

fit$age_group <- age_group
fit$x <- x
fit$y <- y
fit$inference <- "unweighted independent-patient-cluster sandwich; exploratory"

print(fit, row.names = FALSE)

out <- here::here("results", "NHANES_clustPear_exploratory")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

file <- file.path(out, paste0("clustPear_", age_group, "_", x, "_", y, ".csv"))

utils::write.csv(fit, file, row.names = FALSE)

message("Saved exploratory clustPear results: ", file)
message("Do not substitute these confidence intervals for NHANES survey-design intervals.")
