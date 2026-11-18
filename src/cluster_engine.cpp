// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <unordered_map>
#include "v11_core.hpp"
#include "weighted_core.hpp"

using namespace Rcpp;
using namespace opwsim;

// Internal engine: cluster labels are already encoded by R (arbitrary input ID types).
// All rows are complete observed X,Y pairs; singleton clusters are retained.
// [[Rcpp::export]]
List cluster_cor_cpp(IntegerVector id,
                     NumericVector x,
                     NumericVector y,
                     int K,
                     int max_pairs,
                     CharacterVector methods) {
  const int n = id.size();
  if (x.size() != n || y.size() != n || n < 2 || K < 2 || K > 200 || max_pairs < 16)
    stop("Invalid input size, K, or ER pair budget");
  std::unordered_map<int, int> index;
  std::vector<Patient> patients;
  for (int j = 0; j < n; j++) {
    if (IntegerVector::is_na(id[j]) || !std::isfinite(x[j]) || !std::isfinite(y[j]))
      stop("C++ input contains missing ID or non-finite outcomes");
    const int key = id[j];
    auto it = index.find(key);
    int pos;
    if (it == index.end()) {
      pos = int(patients.size());
      index.emplace(key, pos);
      patients.emplace_back();
    } else
      pos = it->second;
    Patient& p = patients[pos];
    p.teeth.push_back({x[j], y[j], x[j], y[j], 0, 0, 0, 0, 1, true});
    p.nobs++;
  }
  const int M = int(patients.size());
  if (M < 2)
    stop("At least two clusters with complete observed pairs are required");
  const std::vector<std::string> possible = {"NW", "CW", "EQ", "ER"};
  std::vector<std::string> requested;
  for (int j = 0; j < methods.size(); j++)
    requested.push_back(as<std::string>(methods[j]));
  const bool need_eq = std::find(requested.begin(), requested.end(), "EQ") != requested.end();
  const bool need_features =
      need_eq || std::find(requested.begin(), requested.end(), "NW") != requested.end() ||
      std::find(requested.begin(), requested.end(), "CW") != requested.end();
  std::vector<Feature> ff;
  if (need_features) {
    Cuts cuts;
    if (need_eq)
      cuts = empirical_cuts(patients, K, false, 1);
    ff.reserve(M);
    for (const Patient& p : patients)
      ff.push_back(features(p, cuts, false));
  }
  NumericVector estimates(requested.size()), ses(requested.size());
  LogicalVector er_exact(requested.size(), NA_LOGICAL);
  IntegerVector pairs_used(requested.size(), NA_INTEGER);
  for (size_t j = 0; j < requested.size(); j++) {
    const std::string& met = requested[j];
    Estimate z;
    if (met == "NW")
      z = evaluate(ff, NW, K);
    else if (met == "CW")
      z = evaluate(ff, CW, K);
    else if (met == "EQ")
      z = evaluate(ff, EQ, K);
    else if (met == "ER") {
      const ErResult er = estimate_er(patients, false, 1, max_pairs, max_pairs, 0, true);
      z = er.e;
      er_exact[j] = er.exact;
      pairs_used[j] = er.pairs_used;
    } else
      stop("Unknown estimation method");
    if (!std::isfinite(z.estimate) || !std::isfinite(z.se))
      stop("Correlation or sandwich SE undefined; check variation in X/Y and cluster count");
    estimates[j] = z.estimate;
    ses[j] = z.se;
  }
  return List::create(_["method"] = methods,
                      _["estimate"] = estimates,
                      _["se"] = ses,
                      _["er_exact"] = er_exact,
                      _["pairs_evaluated"] = pairs_used,
                      _["clusters"] = M,
                      _["observations"] = n);
}

// Optional survey-weighted point estimates and independent-cluster sandwich.
// Weights must be supplied for every observed pair and constant within cluster.
// [[Rcpp::export]]
List cluster_cor_weighted_cpp(IntegerVector id,
                              NumericVector x,
                              NumericVector y,
                              NumericVector w,
                              int K,
                              int max_pairs,
                              CharacterVector methods) {
  const int n = id.size();
  if (x.size() != n || y.size() != n || w.size() != n || n < 2 || K < 2 || K > 200 ||
      max_pairs < 16)
    stop("Invalid weighted input size, K, or ER pair budget");
  std::unordered_map<int, int> index;
  std::vector<Patient> patients;
  std::vector<double> weights;
  for (int j = 0; j < n; j++) {
    if (IntegerVector::is_na(id[j]) || !std::isfinite(x[j]) || !std::isfinite(y[j]) ||
        !std::isfinite(w[j]) || w[j] <= 0)
      stop("Weighted C++ input requires positive finite weights and observed outcomes");
    const int key = id[j];
    auto it = index.find(key);
    int pos;
    if (it == index.end()) {
      pos = int(patients.size());
      index.emplace(key, pos);
      patients.emplace_back();
      weights.push_back(w[j]);
    } else {
      pos = it->second;
      if (std::abs(weights[pos] - w[j]) >
          1e-8 * std::max({1.0, std::abs(w[j]), std::abs(weights[pos])}))
        stop("survey_weights must be constant within each cluster");
    }
    Patient& p = patients[pos];
    p.teeth.push_back({x[j], y[j], x[j], y[j], 0, 0, 0, 0, 1, true});
    p.nobs++;
  }
  const int M = int(patients.size());
  if (M < 2)
    stop("At least two clusters are required");
  std::vector<std::string> requested;
  for (int j = 0; j < methods.size(); j++)
    requested.push_back(as<std::string>(methods[j]));
  const bool need_splits = std::find(requested.begin(), requested.end(), "EQ") != requested.end() ||
                           std::find(requested.begin(), requested.end(), "ER") != requested.end();
  std::vector<RankIndex> indices;
  if (need_splits) {
    indices.reserve(M);
    for (const auto& p : patients)
      indices.emplace_back(p, false);
  }
  NumericVector estimates(requested.size()), ses(requested.size());
  LogicalVector er_exact(requested.size(), NA_LOGICAL);
  IntegerVector pairs_used(requested.size(), NA_INTEGER);
  for (size_t j = 0; j < requested.size(); j++) {
    const auto& met = requested[j];
    Estimate z;
    if (met == "NW" || met == "CW")
      z = clustpear_weighted::moment_estimator(patients, weights, met == "CW");
    else if (met == "EQ") {
      const auto pairs = clustpear_weighted::eq_pairs(patients, K);
      z = clustpear_weighted::pair_average(indices, weights, pairs);
    } else if (met == "ER") {
      bool exact = false;
      const auto pairs = clustpear_weighted::er_pairs(patients, max_pairs, exact);
      z = clustpear_weighted::pair_average(indices, weights, pairs);
      er_exact[j] = exact;
      pairs_used[j] = int(pairs.size());
    } else
      stop("Unknown estimation method");
    if (!std::isfinite(z.estimate) || !std::isfinite(z.se))
      stop("Weighted correlation or sandwich SE undefined; check outcome variation");
    estimates[j] = z.estimate;
    ses[j] = z.se;
  }
  return List::create(_["method"] = methods,
                      _["estimate"] = estimates,
                      _["se"] = ses,
                      _["er_exact"] = er_exact,
                      _["pairs_evaluated"] = pairs_used,
                      _["clusters"] = M,
                      _["observations"] = n);
}

// Generate one sample under the manuscript's clustered bivariate model,
// returning only observed records after outcome-dependent retention.
// [[Rcpp::export]]
DataFrame simulate_clustered_cpp(
    int M, int nmax, double rho_uv, double eta, double sigma2, double tau2, double seed) {
  if (M < 2 || nmax < 1 || std::abs(rho_uv) >= 1 || !std::isfinite(rho_uv) || eta < 0 ||
      !std::isfinite(eta) || sigma2 <= 0 || !std::isfinite(sigma2) || tau2 < 0 ||
      !std::isfinite(tau2) || !std::isfinite(seed) || seed < 0)
    stop("Invalid simulation parameters");
  Spec spec;
  spec.M = M;
  spec.nmax = nmax;
  spec.rho_uv = rho_uv;
  spec.eta = eta;
  spec.sigma2 = sigma2;
  spec.tau2 = tau2;
  spec.seed = uint64_t(seed);
  Rng rng(spec.seed);
  std::vector<int> ids;
  std::vector<double> xs, ys;
  for (int i = 1; i <= M; i++) {
    Patient p = make_patient(rng, spec);
    for (const Tooth& t : p.teeth)
      if (t.kept) {
        ids.push_back(i);
        xs.push_back(t.x);
        ys.push_back(t.y);
      }
  }
  return DataFrame::create(_["cluster_id"] = ids, _["X"] = xs, _["Y"] = ys);
}
