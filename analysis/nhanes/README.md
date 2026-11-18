# NHANES example (not part of the installed package)

Copy `20_nhanes_clustPear_example.R` to the existing OPW research project's
`scripts/` directory and run `source("scripts/20_nhanes_clustPear_example.R")`
from that project root after installing clustPear. It loads the project's
existing prepared data and computes generic **unweighted** patient-cluster
correlations for one age group and outcome pair, without touching the
existing NHANES manuscript outputs. Change `age_group`, `x`, `y`, `K` in the
script to try another combination.

Do **not** use the package output as NHANES survey-design inference. The
existing manuscript workflow incorporates `WTMEC4YR`, `SDMVSTRA`, `SDMVPSU`
and stratified delete-one-PSU jackknife; clustPear does not. Additionally,
the generic engine retains singleton patient clusters and caps the example
ER calculation at 8,192 pairs, so even point estimates are not directly
comparable without controlling those differences.
