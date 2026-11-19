# NHANES examples (not part of the installed package)

`20_nhanes_clustPear_example.R` is an unweighted exploratory example using
prepared tooth data. It does **not** produce design-based NHANES standard errors.

For the journal analysis use `survey_design` as documented in the package's
root README. Build the R `survey::svydesign()` object from the complete merged
MEC person file, *before* filtering to adults, examined dental participants,
10+ surviving teeth, age groups, or complete paired records. Then use the
eligible tooth data in `cluster_pearson(..., survey_design = full_design)`.

The accompanying research-script bundle includes
`scripts/04_nhanes_clustPear_survey.R` and
`scripts/08_make_nhanes_tables.R`. These calculate the survey-Taylor standard
errors and write the corresponding manuscript tables. The optional
`04b_nhanes_design_JKn.R` reruns cutpoints inside PSU replicates for a
potential sensitivity analysis; it is slower and not the default estimator.
