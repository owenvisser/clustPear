# Internal helpers for NHANES-style design-based variance estimation.
# A full person-level survey design is required: the subset of patients with
# complete dental observations supplies estimating-equation contributions,
# while all other examined persons supply zero influence contributions.
.clustpear_validate_design <- function(design, id) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("The survey package is required for survey_design; install.packages('survey')", call. = FALSE)
  }
  if (!inherits(design, "survey.design")) {
    stop("survey_design must be a person-level survey::svydesign() object (not a replicate-weight design)", call. = FALSE)
  }
  if (!is.character(id) || length(id) != 1L || is.na(id) ||
      !id %in% names(design$variables)) {
    stop("survey_id must name a person identifier in survey_design$variables", call. = FALSE)
  }
  person_ids <- as.character(design$variables[[id]])
  if (anyNA(person_ids) || anyDuplicated(person_ids)) {
    stop("survey_design must have exactly one nonmissing row per person identifier", call. = FALSE)
  }
  w <- as.numeric(stats::weights(design, type = "sampling"))
  if (length(w) != length(person_ids) || any(!is.finite(w) | w <= 0)) {
    stop("survey_design requires a positive, finite examination weight for every person", call. = FALSE)
  }
  df <- survey::degf(design)
  if (length(df) != 1L || !is.finite(df) || df <= 0) {
    stop("survey_design has no positive finite design degrees of freedom", call. = FALSE)
  }
  list(design = design, ids = person_ids, weights = w, df = as.numeric(df))
}

.clustpear_design_se <- function(design, design_ids, cluster_ids, linearized) {
  if (!is.matrix(linearized) || nrow(linearized) != length(cluster_ids) ||
      ncol(linearized) < 1L || any(!is.finite(linearized))) {
    stop("Internal error: invalid patient-level linearization", call. = FALSE)
  }
  rows <- match(cluster_ids, design_ids)
  if (anyNA(rows)) stop("Internal error: participants missing from the survey design", call. = FALSE)
  ses <- numeric(ncol(linearized))
  # svytotal applies the design's sampling weights itself. Linearized values
  # are therefore UNWEIGHTED here; incorporating weights twice is incorrect.
  for (j in seq_len(ncol(linearized))) {
    score <- numeric(length(design_ids))
    score[rows] <- linearized[, j]
    work <- design
    variable <- ".clustpear_linearized_internal"
    while (variable %in% names(work$variables)) variable <- paste0(variable, "_")
    work$variables[[variable]] <- score
    f <- stats::as.formula(paste0("~", variable))
    total <- survey::svytotal(f, work)
    ses[j] <- unname(as.numeric(survey::SE(total))[1L])
    if (!is.finite(ses[j]) || ses[j] < 0) {
      stop("Survey Taylor-linearized standard error was not finite", call. = FALSE)
    }
  }
  ses
}
