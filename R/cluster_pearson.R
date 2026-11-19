
#' Cluster-weighted Pearson correlations
#'
#' Calculates NW, CW, EQ and ER Pearson correlations from three columns of an
#' arbitrary data frame. All columns other than the requested inputs are ignored. Rows with missing
#' cluster ID or missing X/Y are omitted pairwise.
#'
#' `categorical` values are converted to numerical scores. Factors use their
#' existing level order (1, 2, ...); character values use alphabetical order.
#' For unordered outcomes with more than two categories, Pearson's correlation
#' depends on the chosen scores/order; use ordered factors for meaningful results.
#' For binary outcomes, Pearson equals the phi correlation when both outcomes
#' are binary, and the point-biserial correlation when only one is binary.
#'
#' EQ and ER use the original numeric/ordinal scores for Pearson moments, but
#' use their ordered distinct values to construct OPW subgroup splits.
#' By default, SEs use the independent-cluster sandwich and delta method.
#' With `survey_design`, SEs use the full survey's strata/PSUs through a
#' Taylor-linearized survey total. EQ/ER SEs hold the empirical cutpoints fixed
#' in both modes; confidence intervals do not incorporate cutpoint-selection
#' uncertainty. ER reports whether its evaluated
#' rank grid was exhaustive or capped; a capped ER result is an approximation.
#'
#' @param data A data frame containing the cluster ID, X, and Y columns.
#' @param cluster,x,y Column names as strings.
#' @param x_type,y_type Either "continuous" (default) or "categorical".
#' @param method One of "all" (default), "NW", "CW", "EQ", "ER".
#' @param K EQ rank-grid resolution (default 20).
#' @param er_max_pairs Maximum evaluated split pairs for ER. If its exhaustive
#'   grid exceeds this budget, a deterministic rank-grid approximation is used.
#' @param conf_level Nominal confidence level; 0.95 by default.
#' @param survey_design Optional person-level `survey::svydesign()` object.
#'   When supplied, point estimates use the design weights and standard errors
#'   use stratified-PSU Taylor linearization via `survey::svytotal()`.
#'   The design must contain all sampled persons with positive examination
#'   weights, including those outside the dental or age analysis domain. It
#'   must contain exactly one row per person. Cannot be combined with
#'   `survey_weights`. Requires the suggested `survey` package.
#' @param survey_id Name of the person identifier in the design; defaults
#'   to the value of `cluster` (e.g., `SEQN` for NHANES).
#' @param survey_weights Optional positive cluster-constant sampling weights.
#'   Specify a column name as a string, or a numeric vector of length
#'   `nrow(data)` (one weight for each row). Default `NULL` uses equal
#'   cluster sampling weights, as in the original package. Weights modify
#'   the estimating equations; they are *not* multiplied into X or Y.
#'   The resulting SE remains an independent-cluster sandwich SE, not a
#'   complex-survey-design SE accounting for strata or PSUs. To obtain such
#'   SEs, pass a complete person-level `survey_design` instead.
#' @return A plain data.frame with correlation estimates, standard errors,
#'   confidence bounds, cluster counts, and ER exactness status. Survey mode
#'   additionally returns `variance_method`, `design_df`, and the original
#'   independent-patient sandwich SE as `se_independent` for comparison.
#' @export
#' @examples
#' d <- simulate_clustered()
#' cluster_pearson(d, "cluster_id", "X", "Y")

cluster_pearson <- function(
  data, cluster, x, y,
  x_type = "continuous", y_type = "continuous",
  method = "all", K = 20L,
  er_max_pairs = 100000L, conf_level = 0.95,
  survey_weights = NULL,
  survey_design = NULL, survey_id = cluster
) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }

  for (nm in c(cluster, x, y)) {
    if (!is.character(nm) || length(nm) != 1L || is.na(nm) || !nm %in% names(data)) {
      stop("cluster, x and y must name existing columns", call. = FALSE)
    }
  }

  x_type <- match.arg(x_type, c("continuous", "categorical"))
  y_type <- match.arg(y_type, c("continuous", "categorical"))
  method <- match.arg(method, c("all", "NW", "CW", "EQ", "ER"))

  if (
    length(K) != 1L ||
    !is.numeric(K) ||
    !is.finite(K) ||
    K != as.integer(K) ||
    K < 2 ||
    K > 200
  ) {
    stop("K must be an integer between 2 and 200", call. = FALSE)
  }

  if (
    length(er_max_pairs) != 1L ||
    !is.numeric(er_max_pairs) ||
    !is.finite(er_max_pairs) ||
    er_max_pairs != as.integer(er_max_pairs) ||
    er_max_pairs < 16
  ) {
    stop("er_max_pairs must be an integer >= 16", call. = FALSE)
  }

  if (
    length(conf_level) != 1L ||
    !is.numeric(conf_level) ||
    !is.finite(conf_level) ||
    conf_level <= 0 ||
    conf_level >= 1
  ) {
    stop("conf_level must lie strictly between 0 and 1", call. = FALSE)
  }

  encode <- function(z, type, label) {
    if (type == "continuous") {
      if (!is.numeric(z) || is.factor(z)) {
        stop(label, " must be numeric for type='continuous'", call. = FALSE)
      }

      return(as.numeric(z))
    }

    if (is.factor(z)) {
      if (!is.ordered(z) && nlevels(z) > 2L) {
        warning(
          label,
          ": unordered multi-level factor uses its stored level order as Pearson scores",
          call. = FALSE
        )
      }

      return(as.numeric(z))
    }

    if (is.character(z)) {
      if (length(unique(z[!is.na(z)])) > 2L) {
        warning(
          label,
          ": multi-level character categories use alphabetical scores for Pearson",
          call. = FALSE
        )
      }

      return(as.numeric(factor(z)))
    }

    if (is.logical(z)) {
      return(as.numeric(z))
    }

    if (is.numeric(z)) {
      return(as.numeric(z))
    }

    stop(label, " must be numeric, logical, factor, or character", call. = FALSE)
  }

  if (!is.null(survey_design) && !is.null(survey_weights)) {
    stop("Provide either survey_design (design-based SEs) or survey_weights (independent-cluster SEs), not both", call. = FALSE)
  }
  weighted <- !is.null(survey_weights) || !is.null(survey_design)
  survey_info <- NULL
  if (!is.null(survey_design)) {
    survey_info <- .clustpear_validate_design(survey_design, survey_id)
    target_id <- as.character(data[[cluster]])
    pos <- match(target_id, survey_info$ids)
    w <- rep.int(NA_real_, nrow(data))
    w[!is.na(pos)] <- survey_info$weights[pos[!is.na(pos)]]
  } else if (weighted) {
    if (
      is.character(survey_weights) &&
      length(survey_weights) == 1L &&
      !is.na(survey_weights)
    ) {
      if (!survey_weights %in% names(data)) {
        stop("survey_weights column not found: ", survey_weights, call. = FALSE)
      }

      w <- data[[survey_weights]]
    } else {
      w <- survey_weights
    }

    if (!is.numeric(w) || length(w) != nrow(data)) {
      stop(
        "survey_weights must be a column name or numeric vector with one value per row",
        call. = FALSE
      )
    }
  }

  xs <- encode(data[[x]], x_type, x)
  ys <- encode(data[[y]], y_type, y)
  ids <- data[[cluster]]

  keep <- !is.na(ids) & is.finite(xs) & is.finite(ys)
  cluster_ids <- unique(as.character(ids[keep]))
  if (!is.null(survey_info) && anyNA(match(cluster_ids, survey_info$ids))) {
    stop("Some analysis clusters are missing from survey_design; ensure its full person frame contains the dental participants", call. = FALSE)
  }

  ids <- match(ids[keep], unique(ids[keep]))
  xs <- xs[keep]
  ys <- ys[keep]

  if (weighted) {
    w <- as.numeric(w[keep])

    if (any(!is.finite(w) | w <= 0)) {
      stop("survey_weights must be positive and finite for included rows", call. = FALSE)
    }
  }

  if (length(xs) < 2L || length(unique(ids)) < 2L) {
    stop("At least two observed pairs in at least two clusters are required", call. = FALSE)
  }

  if (length(unique(xs)) < 2L || length(unique(ys)) < 2L) {
    stop("Both X and Y must have at least two distinct observed values", call. = FALSE)
  }

  requested <- if (method == "all") c("NW", "CW", "EQ", "ER") else method

  z <- if (weighted) {
    cluster_cor_weighted_cpp(
      ids, xs, ys, w,
      as.integer(K),
      as.integer(er_max_pairs),
      requested
    )
  } else {
    cluster_cor_cpp(
      ids, xs, ys,
      as.integer(K),
      as.integer(er_max_pairs),
      requested
    )
  }

  ses <- as.numeric(z$se)
  if (is.null(survey_info)) {
    quant <- stats::qnorm((1 + conf_level) / 2)
  } else {
    ses <- .clustpear_design_se(
      survey_info$design, survey_info$ids, cluster_ids, z$linearized
    )
    quant <- stats::qt((1 + conf_level) / 2, df = survey_info$df)
  }

  ans <- data.frame(
    method = as.character(z$method),
    pearson = as.numeric(z$estimate),
    se = ses,
    lower = as.numeric(z$estimate) - quant * ses,
    upper = as.numeric(z$estimate) + quant * ses,
    n_clusters = rep.int(z$clusters, length(requested)),
    n_obs = rep.int(z$observations, length(requested)),
    er_exact = as.logical(z$er_exact),
    row.names = NULL
  )
  if (!is.null(survey_info)) {
    ans$variance_method <- "survey_taylor"
    ans$design_df <- survey_info$df
    ans$se_independent <- as.numeric(z$se)
  }
  ans
}
