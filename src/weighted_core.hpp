#ifndef CLUSTPEAR_WEIGHTED_CORE_HPP
#define CLUSTPEAR_WEIGHTED_CORE_HPP
#include "v11_core.hpp"
#include <utility>

// Cluster-constant sampling weights multiply the estimating-equation
// contributions, not the observed outcome values. The mathematical formulas
// mirror OPW V11's nhanes_survey_core.hpp, without NHANES-specific assumptions.
namespace clustpear_weighted {
using opwsim::Mom;
using opwsim::Patient;
using opwsim::RankIndex;
using opwsim::RankPair;
using opwsim::Cuts;
using opwsim::Estimate;
using opwsim::moments;
using opwsim::pearson;
using opwsim::gradient;
using opwsim::dot;

inline Estimate moment_estimator(const std::vector<Patient>& patients,
                                 const std::vector<double>& weights,
                                 bool cw,
                                 std::vector<double>* linearized = nullptr) {
  const int M = int(patients.size());
  Mom numer;
  double denom = 0.0;
  std::vector<Mom> sums(M);
  for (int i = 0; i < M; i++) {
    for (const auto& t : patients[i].teeth)
      sums[i] += moments(t.x, t.y);
    const int n = patients[i].nobs;
    const Mom contrib = cw ? sums[i] / double(n) : sums[i];
    numer += contrib * weights[i];
    denom += weights[i] * (cw ? 1.0 : double(n));
  }
  const Mom mean = numer / denom;
  const double rho = pearson(mean);
  if (!std::isfinite(rho))
    throw std::runtime_error("Undefined weighted correlation");
  const Mom grad = gradient(mean);
  if (linearized) linearized->assign(M, 0.0);
  double sq = 0.0;
  for (int i = 0; i < M; i++) {
    const Mom contrib = cw ? sums[i] / double(patients[i].nobs) : sums[i];
    const double n = cw ? 1.0 : double(patients[i].nobs);
    const double unit_infl = dot(grad, contrib - mean * n) / denom;
    if (linearized) (*linearized)[i] = unit_infl;
    const double infl = unit_infl * weights[i];
    sq += infl * infl;
  }
  return {rho, std::sqrt(std::max(0.0, double(M) * sq / (M - 1.0)))};
}

inline Estimate pair_average(const std::vector<RankIndex>& indices,
                             const std::vector<double>& weights,
                             const std::vector<RankPair>& pairs,
                             std::vector<double>* linearized = nullptr) {
  const int M = int(indices.size()), Q = int(pairs.size());
  if (Q < 1)
    throw std::invalid_argument("No available split pairs");
  const double W = std::accumulate(weights.begin(), weights.end(), 0.0);
  std::vector<Mom> means(Q), grads(Q);
  for (int i = 0; i < M; i++) {
    for (int q = 0; q < Q; q++)
      means[q] += indices[i].query(pairs[q].x, pairs[q].y) * weights[i];
  }
  double rho = 0.0;
  for (int q = 0; q < Q; q++) {
    means[q] *= 1.0 / W;
    const double r = pearson(means[q]);
    if (!std::isfinite(r))
      throw std::runtime_error("Undefined weighted correlation at a split");
    rho += r;
    grads[q] = gradient(means[q]);
  }
  if (linearized) linearized->assign(M, 0.0);
  double sq = 0.0;
  for (int i = 0; i < M; i++) {
    double infl = 0.0;
    for (int q = 0; q < Q; q++)
      infl += dot(grads[q], indices[i].query(pairs[q].x, pairs[q].y) - means[q]);
    const double unit_infl = infl / (W * Q);
    if (linearized) (*linearized)[i] = unit_infl;
    infl = weights[i] * unit_infl;
    sq += infl * infl;
  }
  return {rho / Q, std::sqrt(std::max(0.0, double(M) * sq / (M - 1.0)))};
}

inline std::vector<RankPair> eq_pairs(const std::vector<Patient>& patients, int K) {
  const Cuts cuts = opwsim::empirical_cuts(patients, K, false, 1);
  std::vector<RankPair> out;
  out.reserve(cuts.x.size() * cuts.y.size());
  for (double x : cuts.x)
    for (double y : cuts.y)
      out.push_back({x, y, 0});
  return out;
}

// Same exhaustive/capped rank selection as V11 estimate_er(..., max_pairs,
// max_pairs, ...), but point estimates and influence functions use w_i.
inline std::vector<RankPair>
er_pairs(const std::vector<Patient>& patients, int max_pairs, bool& exact) {
  const auto vals = opwsim::rank_values(patients, false, 1);
  const int DX = int(vals.x.size()), DY = int(vals.y.size());
  const uint64_t Q = uint64_t(DX - 1) * uint64_t(DY - 1);
  exact = Q <= uint64_t(max_pairs);
  std::vector<double> xs, ys;
  if (exact) {
    xs.assign(vals.x.begin() + 1, vals.x.end());
    ys.assign(vals.y.begin() + 1, vals.y.end());
  } else {
    int lo = 2, hi = std::max(DX, DY);
    while (lo < hi) {
      const int mid = lo + (hi - lo + 1) / 2;
      const uint64_t count = uint64_t(std::min(mid, DX) - 1) * uint64_t(std::min(mid, DY) - 1);
      if (count <= uint64_t(max_pairs))
        lo = mid;
      else
        hi = mid - 1;
    }
    const Cuts cuts = opwsim::empirical_cuts(patients, lo, false, 1);
    xs = cuts.x;
    ys = cuts.y;
  }
  std::vector<RankPair> out;
  out.reserve(xs.size() * ys.size());
  for (double x : xs)
    for (double y : ys)
      out.push_back({x, y, 0});
  return out;
}
} // namespace clustpear_weighted
#endif
