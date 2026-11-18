#ifndef OPW_RANK_QUANTILE_SIM_HPP
#define OPW_RANK_QUANTILE_SIM_HPP
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace opwsim {
using std::size_t;
inline double expit(double x) {
  if (x >= 0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

struct Mom {
  std::array<double, 5> v{};
  Mom& operator+=(const Mom& b) {
    for (int i = 0; i < 5; i++)
      v[i] += b.v[i];
    return *this;
  }
  Mom& operator-=(const Mom& b) {
    for (int i = 0; i < 5; i++)
      v[i] -= b.v[i];
    return *this;
  }
  Mom& operator*=(double z) {
    for (double& a : v)
      a *= z;
    return *this;
  }
};

inline Mom operator+(Mom a, const Mom& b) {
  return a += b;
}

inline Mom operator-(Mom a, const Mom& b) {
  return a -= b;
}

inline Mom operator*(Mom a, double z) {
  return a *= z;
}

inline Mom operator/(Mom a, double z) {
  return a *= 1.0 / z;
}

inline double dot(const Mom& a, const Mom& b) {
  double z = 0;
  for (int i = 0; i < 5; i++)
    z += a.v[i] * b.v[i];
  return z;
}

inline Mom moments(double x, double y) {
  Mom h;
  h.v = {x, y, x * x, y * y, x * y};
  return h;
}

inline double pearson(const Mom& m) {
  const double vx = m.v[2] - m.v[0] * m.v[0], vy = m.v[3] - m.v[1] * m.v[1];
  if (vx <= 0 || vy <= 0)
    return std::numeric_limits<double>::quiet_NaN();
  return (m.v[4] - m.v[0] * m.v[1]) / std::sqrt(vx * vy);
}

inline Mom gradient(const Mom& m) {
  const double vx = m.v[2] - m.v[0] * m.v[0], vy = m.v[3] - m.v[1] * m.v[1];
  const double denom = std::sqrt(vx * vy), r = pearson(m);
  Mom g;
  g.v = {-m.v[1] / denom + r * m.v[0] / vx,
         -m.v[0] / denom + r * m.v[1] / vy,
         -r / (2 * vx),
         -r / (2 * vy),
         1 / denom};
  return g;
}
// Method indices and names stay aligned with Rcpp and all analysis scripts.
enum Method { NW = 0, CW = 1, OPW = 2, EQ = 3, ER = 4, NMETHODS = 5 };
inline constexpr const char* method_names[] = {"NW", "CW", "OPW", "EQ", "ER"};
struct Spec {
  int M = 100, nmax = 100, min_n = 2, reps = 1000, reference_clusters = 2000;
  int er_pairs = 100000, reference_er_pairs = 8192, er_exact_limit = 100000;
  int eq_boot_B = 0;
  double rho_uv = .25, sigma2 = .25, tau2 = .25, eta = 2;
  uint64_t seed = 20261103;
  std::vector<std::vector<double>> cuts;
};

struct Tooth {
  double xs, ys, x, y, u, v, z, r, pi;
  bool kept;
};

struct Patient {
  std::vector<Tooth> teeth;
  int nobs = 0;
};

struct Rng {
  std::mt19937_64 engine;
  std::normal_distribution<double> norm{0, 1};
  std::uniform_real_distribution<double> unif{0, 1};
  explicit Rng(uint64_t seed) : engine(seed) {
  }
  double n() {
    return norm(engine);
  }
  double u() {
    return unif(engine);
  }
};

inline void validate(const Spec& s) {
  if (s.M < 2 || s.nmax < 2 || s.min_n < 2 || s.min_n > s.nmax || s.reps < 1 ||
      s.reference_clusters < 8 || s.sigma2 <= 0 || s.tau2 < 0 || s.eta < 0 ||
      !std::isfinite(s.rho_uv) || std::abs(s.rho_uv) >= 1 || !std::isfinite(s.sigma2) ||
      !std::isfinite(s.tau2) || !std::isfinite(s.eta) || s.cuts.empty() || s.er_pairs < 16 ||
      s.reference_er_pairs < 16 || s.er_pairs % 4 != 0 || s.reference_er_pairs % 4 != 0 ||
      s.er_exact_limit < 1 || s.eq_boot_B != 0)
    throw std::invalid_argument(
        "EQ bootstrap is disabled: all simulation methods use fixed-cutpoint sandwich SEs");
  for (const auto& c : s.cuts) {
    if (c.empty() || c.size() > 199)
      throw std::invalid_argument("K must be 2..200");
    for (size_t i = 0; i < c.size(); i++) {
      if (!std::isfinite(c[i]) || (i > 0 && c[i] <= c[i - 1]))
        throw std::invalid_argument("Fixed cutpoints must be strictly increasing");
    }
  }
}

inline Patient make_patient(Rng& rng, const Spec& s) {
  const double u = rng.n(), v = s.rho_uv * u + std::sqrt(1 - s.rho_uv * s.rho_uv) * rng.n();
  const double sd = std::sqrt(s.sigma2), zsd = std::sqrt(s.tau2);
  const double marginal_sd = std::sqrt(1 + s.sigma2);
  Patient p;
  p.teeth.reserve(s.nmax);
  for (int j = 0; j < s.nmax; j++) {
    const double z = (u + v) * 0.5 + zsd * rng.n(), r = expit(-z);
    const double a = rng.n(), b = rng.n();
    const double x = u + sd * a, y = v + sd * (r * a + std::sqrt(std::max(0.0, 1 - r * r)) * b);
    const double xs = x / marginal_sd, ys = y / marginal_sd;
    // Matches manuscript Eq. (9): eta=0 retains every potential unit.
    const double pi = std::exp(-s.eta * expit(-xs) * expit(-ys));
    const bool kept = rng.u() < pi;
    p.teeth.push_back({xs, ys, x, y, u, v, z, r, pi, kept});
    p.nobs += int(kept);
  }
  return p;
}

inline double rho0(const Spec& s) {
  return (s.rho_uv + 0.5 * s.sigma2) / (1 + s.sigma2);
}

struct Cuts {
  std::vector<double> x, y;
};
// Select K_X-1 distinct cutpoint ranks from D_X-1 available ER cuts.
// ER cut r separates the first r distinct values from the remaining values.
// We store the first value in the UPPER group, rather than its midpoint:
// RankIndex::query uses x < upper_value, and Feature::features uses
// upper_bound(cuts,x), so both implement exactly the same dichotomization.
inline std::vector<double> distinct_rank_cutpoints(std::vector<double> values, int K) {
  if (K < 2)
    throw std::invalid_argument("Rank-grid EQ requires K >= 2");
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
  const int D = int(values.size());
  if (D < 2)
    throw std::runtime_error("EQ requires two distinct values per outcome");
  const int K_out = std::min(K, D);
  std::vector<double> cuts;
  cuts.reserve(K_out - 1);
  for (int h = 1; h < K_out; h++) {
    const int rank = int((int64_t(h) * (D - 1) + K_out - 1) / K_out);
    if (rank < 1 || rank >= D)
      throw std::logic_error("Invalid distinct-rank cut index");
    cuts.push_back(values[rank]); // first value of upper group (zero-indexed)
  }
  if (!std::is_sorted(cuts.begin(), cuts.end()) ||
      std::adjacent_find(cuts.begin(), cuts.end()) != cuts.end())
    throw std::logic_error("Distinct-rank grid repeated a cut");
  return cuts;
}
// The selected ranks are ceil[h(D_X-1)/K_X], h=1,...,K_X-1,
// K_X=min(K,D_X); Y is analogous. Ties do not receive extra weight.
inline Cuts
empirical_cuts(const std::vector<Patient>& patients, int K, bool use_complete, int min_n) {
  std::vector<double> xs, ys;
  for (const Patient& p : patients) {
    if (!use_complete && p.nobs < min_n)
      continue;
    for (const Tooth& t : p.teeth)
      if (use_complete || t.kept) {
        xs.push_back(t.xs);
        ys.push_back(t.ys);
      }
  }
  return {distinct_rank_cutpoints(std::move(xs), K), distinct_rank_cutpoints(std::move(ys), K)};
}

struct Feature {
  int n = 0, n_pairs = 0;
  Mom sum, cw, opw;
  std::vector<Mom> binary;
};
// One O(K^2) 2D prefix calculation yields every 2x2 binary split and OPW.
inline Feature features(const Patient& p, const Cuts& cuts, bool use_complete) {
  // Rectangular prefix grid accommodates K_X != K_Y with unique cuts.
  const int KX = int(cuts.x.size()) + 1, KY = int(cuts.y.size()) + 1;
  const int Q = (KX - 1) * (KY - 1), W = KY + 1;
  std::vector<int> cnt((KX + 1) * W, 0);
  std::vector<Mom> block((KX + 1) * W);
  Feature f;
  f.binary.resize(Q);
  for (const Tooth& t : p.teeth) {
    if (!use_complete && !t.kept)
      continue;
    const int k = int(std::upper_bound(cuts.x.begin(), cuts.x.end(), t.xs) - cuts.x.begin());
    const int l = int(std::upper_bound(cuts.y.begin(), cuts.y.end(), t.ys) - cuts.y.begin());
    const int id = (k + 1) * W + (l + 1);
    Mom h = moments(t.x, t.y);
    cnt[id]++;
    block[id] += h;
    f.sum += h;
    f.n++;
  }
  if (f.n == 0)
    return f;
  f.cw = f.sum / double(f.n);
  for (int k = 1; k <= KX; k++)
    for (int l = 1; l <= KY; l++) {
      int id = k * W + l;
      if (cnt[id]) {
        f.opw += block[id] / double(cnt[id]);
        f.n_pairs++;
      }
    }
  f.opw *= 1.0 / double(f.n_pairs);
  for (int k = 1; k <= KX; k++)
    for (int l = 1; l <= KY; l++) {
      int id = k * W + l;
      cnt[id] += cnt[id - W] + cnt[id - 1] - cnt[id - W - 1];
      block[id] += block[id - W] + block[id - 1] - block[id - W - 1];
    }
  const int N = cnt[KX * W + KY];
  const Mom H = block[KX * W + KY];
  for (int a = 1; a < KX; a++)
    for (int b = 1; b < KY; b++) {
      int i = a * W + b;
      int n00 = cnt[i], n10 = cnt[KX * W + b] - n00, n01 = cnt[a * W + KY] - n00;
      int n11 = N - n00 - n10 - n01;
      Mom h00 = block[i], h10 = block[KX * W + b] - h00;
      Mom h01 = block[a * W + KY] - h00, h11 = H - h00 - h10 - h01;
      Mom h;
      int parts = 0;
      if (n00) {
        h += h00 / double(n00);
        parts++;
      }
      if (n10) {
        h += h10 / double(n10);
        parts++;
      }
      if (n01) {
        h += h01 / double(n01);
        parts++;
      }
      if (n11) {
        h += h11 / double(n11);
        parts++;
      }
      if (!parts)
        throw std::runtime_error("Empty EQ grid cell");
      f.binary[(a - 1) * (KY - 1) + (b - 1)] = h / double(parts);
    }
  return f;
}

struct Estimate {
  double estimate = std::numeric_limits<double>::quiet_NaN();
  double se = std::numeric_limits<double>::quiet_NaN();
};
// Cumulative mean and sandwich IF across all selected binary split correlations.
inline Estimate evaluate(const std::vector<Feature>& patients, int method, int K) {
  int M = int(patients.size());
  if (M < 2)
    return {};
  if (method == EQ) {
    const int Q = int(patients[0].binary.size());
    if (Q < 1)
      throw std::runtime_error("EQ has no nondegenerate paired split");
    std::vector<Mom> means(Q), derivs(Q);
    int selected = 0;
    double corr_sum = 0;
    for (const Feature& f : patients)
      for (int q = 0; q < Q; q++)
        means[q] += f.binary[q];
    for (int q = 0; q < Q; q++) {
      means[q] *= 1.0 / double(M);
      corr_sum += pearson(means[q]);
      derivs[q] = gradient(means[q]);
      selected++;
    }
    double sq = 0;
    for (const Feature& f : patients) {
      double influence = 0;
      for (int q = 0; q < Q; q++)
        influence += dot(derivs[q], f.binary[q] - means[q]);
      influence /= double(selected);
      sq += influence * influence;
    }
    return {corr_sum / double(selected), std::sqrt(std::max(0.0, sq / (double(M) * (M - 1.0))))};
  }
  Mom mean;
  double total = 0;
  for (const Feature& f : patients) {
    if (method == NW) {
      mean += f.sum;
      total += f.n;
    } else
      mean += (method == CW ? f.cw : f.opw);
  }
  mean *= 1.0 / (method == NW ? total : double(M));
  Mom g = gradient(mean);
  double sq = 0;
  for (const Feature& f : patients) {
    Mom residual =
        (method == NW ? f.sum - mean * double(f.n) : (method == CW ? f.cw : f.opw) - mean);
    double influence = dot(g, residual);
    if (method == NW)
      influence *= double(M) / total;
    sq += influence * influence;
  }
  return {pearson(mean), std::sqrt(std::max(0.0, sq / (double(M) * (M - 1.0))))};
}
// Compute only the EQ point estimate for bootstrap replicates. Reusing
// evaluate(EQ) here would also compute a fixed-cutpoint delta SE and all
// cluster influence values, which are discarded by the percentile bootstrap.
inline double eq_point_only(const std::vector<Feature>& patients) {
  const int M = int(patients.size());
  if (M < 2)
    throw std::runtime_error("EQ requires two eligible clusters");
  const size_t Q = patients.front().binary.size();
  if (Q < 1)
    throw std::runtime_error("EQ has no nondegenerate paired split");
  std::vector<Mom> means(Q);
  for (const Feature& f : patients) {
    if (f.binary.size() != Q)
      throw std::logic_error("Inconsistent bootstrap EQ grid");
    for (size_t q = 0; q < Q; q++)
      means[q] += f.binary[q];
  }
  double sum = 0.0;
  for (Mom& m : means) {
    m *= 1.0 / double(M);
    const double r = pearson(m);
    if (!std::isfinite(r))
      throw std::runtime_error("Non-finite EQ bootstrap correlation");
    sum += r;
  }
  return sum / double(Q);
}
// Whole-cluster EQ bootstrap: re-draw clusters, rebuild distinct-rank cuts,
// and recompute every split-specific OPW correlation.
// B>=2; a missing SE is preferable to labelling the fixed-cutpoint delta SE
// as if it accounted for stochastic empirical cutpoints.
struct EqBootstrap {
  double se, lower, upper;
};

inline double percentile_sorted(const std::vector<double>& x, double p) {
  // Type-7 percentile, matching the default stats::quantile() in R.
  const double position = (x.size() - 1) * p;
  const size_t j = size_t(std::floor(position));
  if (j + 1 >= x.size())
    return x.back();
  const double fraction = position - j;
  return x[j] * (1 - fraction) + x[j + 1] * fraction;
}

inline EqBootstrap eq_bootstrap_interval(
    const std::vector<Patient>& patients, int K, int B, uint64_t seed, int min_n = 2) {
  if (B < 2)
    throw std::invalid_argument("EQ bootstrap requires B>=2");
  std::vector<const Patient*> eligible;
  for (const Patient& patient : patients)
    if (patient.nobs >= min_n)
      eligible.push_back(&patient);
  if (eligible.size() < 2)
    throw std::invalid_argument("EQ bootstrap needs two eligible clusters");
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<int> draw(0, int(eligible.size()) - 1);
  std::vector<double> estimates;
  estimates.reserve(B);
  for (int b = 0; b < B; b++) {
    std::vector<Patient> sampled;
    sampled.reserve(eligible.size());
    for (size_t i = 0; i < eligible.size(); i++)
      sampled.push_back(*eligible[draw(gen)]);
    const Cuts cuts = empirical_cuts(sampled, K, false, min_n);
    std::vector<Feature> features_b;
    features_b.reserve(sampled.size());
    for (const Patient& patient : sampled)
      features_b.push_back(features(patient, cuts, false));
    const double value = eq_point_only(features_b);
    if (!std::isfinite(value))
      throw std::runtime_error("Non-finite EQ bootstrap estimate");
    estimates.push_back(value);
  }
  const double avg = std::accumulate(estimates.begin(), estimates.end(), 0.0) / B;
  double s2 = 0;
  for (double r : estimates)
    s2 += (r - avg) * (r - avg);
  const double se = std::sqrt(s2 / (B - 1));
  std::sort(estimates.begin(), estimates.end());
  return {se, percentile_sorted(estimates, .025), percentile_sorted(estimates, .975)};
}
// For large independent reference samples, avoid storing features for all patients.
struct Accum {
  int n = 0;
  double total = 0, pairs = 0;
  Mom nw, cw, opw;
  std::vector<Mom> binary;
  explicit Accum(int Q = 0) : binary(Q) {
  }
  void add(const Feature& f) {
    if (f.n == 0)
      return;
    n++;
    total += f.n;
    pairs += f.n_pairs;
    nw += f.sum;
    cw += f.cw;
    opw += f.opw;
    for (size_t j = 0; j < binary.size(); j++)
      binary[j] += f.binary[j];
  }
  double get(int method, int /*K*/) const {
    if (n < 2)
      return std::numeric_limits<double>::quiet_NaN();
    if (method == NW)
      return pearson(nw / total);
    if (method == CW)
      return pearson(cw / double(n));
    if (method == OPW)
      return pearson(opw / double(n));
    double out = 0;
    int selected = 0;
    for (int j = 0; j < (int)binary.size(); j++) {
      out += pearson(binary[j] / double(n));
      selected++;
    }
    return out / selected;
  }
};

struct SimRow {
  int rep = 0, K = 0, method = 0, included = 0, excluded = 0;
  double estimate = 0, se = 0, se_fixed = 0, fullrep = 0, mean_n = 0, retained_fraction = 0,
         mean_pairs = 0;
  int eq_boot_B = 0;
  double ci_lower = std::numeric_limits<double>::quiet_NaN();
  double ci_upper = std::numeric_limits<double>::quiet_NaN();
  double numerical_se = 0;
  int er_pairs_used = 0, er_reference_K = 0;
  bool er_exact = false;
};

struct RefRow {
  int K = 0, method = 0, included = 0;
  double target_full = 0, se_full = 0, target_observed = 0, se_observed = 0;
  double retained_fraction = 0, mean_n = 0, excluded_fraction = 0, mean_pairs = 0;
  double numerical_se_full = 0, numerical_se_observed = 0;
  int er_pairs_used = 0;
};

struct SimOutput {
  std::vector<SimRow> rows;
};

struct RefOutput {
  std::vector<RefRow> rows;
};

inline bool is_empirical(int method) {
  return method == EQ || method == ER;
}

inline Cuts fixed_cuts(const std::vector<double>& c) {
  return {c, c};
}

// Fast four-quadrant queries: each patient has at most nmax teeth, so a local
// O(n_i^2) prefix grid gives O(log n_i) OPW moments for ANY global rank split.
struct Cell {
  int n = 0;
  Mom h;
};

struct RankIndex {
  int n = 0, w = 0;
  std::vector<double> x, y;
  std::vector<Cell> a;
  explicit RankIndex(const Patient& p, bool complete) {
    std::vector<const Tooth*> ts;
    for (const Tooth& t : p.teeth)
      if (complete || t.kept)
        ts.push_back(&t);
    n = int(ts.size());
    w = n + 1;
    x.reserve(n);
    y.reserve(n);
    for (const Tooth* t : ts) {
      x.push_back(t->xs);
      y.push_back(t->ys);
    }
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    a.resize(size_t(w) * w);
    for (const Tooth* t : ts) {
      int i = int(std::lower_bound(x.begin(), x.end(), t->xs) - x.begin()) + 1;
      int j = int(std::lower_bound(y.begin(), y.end(), t->ys) - y.begin()) + 1;
      Cell& c = a[size_t(i) * w + j];
      c.n++;
      c.h += moments(t->x, t->y);
    }
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++) {
        Cell& z = a[size_t(i) * w + j];
        const Cell& top = a[size_t(i - 1) * w + j];
        const Cell& left = a[size_t(i) * w + j - 1];
        const Cell& diag = a[size_t(i - 1) * w + j - 1];
        z.n += top.n + left.n - diag.n;
        z.h += top.h + left.h - diag.h;
      }
  }
  Mom query(double xcut, double ycut) const {
    int i = int(std::lower_bound(x.begin(), x.end(), xcut) - x.begin());
    int j = int(std::lower_bound(y.begin(), y.end(), ycut) - y.begin());
    const Cell& z00 = a[size_t(i) * w + j];
    const Cell& zx = a[size_t(n) * w + j];
    const Cell& zy = a[size_t(i) * w + n];
    const Cell& zall = a[size_t(n) * w + n];
    const int n00 = z00.n, n10 = zx.n - n00, n01 = zy.n - n00;
    const int n11 = n - n00 - n10 - n01;
    const Mom& h00 = z00.h;
    Mom h10 = zx.h - h00, h01 = zy.h - h00, h11 = zall.h - h00 - h10 - h01;
    Mom out;
    int parts = 0;
    if (n00) {
      out += h00 / double(n00);
      parts++;
    }
    if (n10) {
      out += h10 / double(n10);
      parts++;
    }
    if (n01) {
      out += h01 / double(n01);
      parts++;
    }
    if (n11) {
      out += h11 / double(n11);
      parts++;
    }
    if (!parts)
      throw std::runtime_error("Empty rank index");
    return out / double(parts);
  }
};

struct RankPair {
  double x, y;
  int block = 0;
};

struct ErResult {
  Estimate e;
  double numerical_se = 0;
  int pairs_used = 0, reference_K = 0;
  bool exact = false;
};

inline double radical_inverse(uint32_t n, uint32_t base) {
  double f = 1.0 / base, z = 0;
  while (n) {
    z += f * (n % base);
    n /= base;
    f /= base;
  }
  return z;
}
// Return the distinct attainable empirical splits. Cutpoint t is the first
// value in the upper group (so RankIndex::query tests x < t, equivalent to
// x <= the preceding distinct value). No stochastic integration is permitted.
struct RankValues {
  std::vector<double> x, y;
};

inline RankValues rank_values(const std::vector<Patient>& patients, bool complete, int min_n) {
  RankValues out;
  for (const Patient& p : patients) {
    if (!complete && p.nobs < min_n)
      continue;
    for (const Tooth& t : p.teeth)
      if (complete || t.kept) {
        out.x.push_back(t.xs);
        out.y.push_back(t.ys);
      }
  }
  for (auto* v : {&out.x, &out.y}) {
    std::sort(v->begin(), v->end());
    v->erase(std::unique(v->begin(), v->end()), v->end());
    if (v->size() < 2)
      throw std::runtime_error("ER requires two distinct values per outcome");
  }
  return out;
}

inline std::vector<RankPair> rank_pairs(const std::vector<Patient>& patients,
                                        bool complete,
                                        int min_n,
                                        int /*er_pairs*/,
                                        int /*exact_limit*/,
                                        uint64_t /*seed*/,
                                        bool& exact) {
  const auto values = rank_values(patients, complete, min_n);
  const size_t Q = (values.x.size() - 1) * (values.y.size() - 1);
  std::vector<RankPair> pairs;
  pairs.reserve(Q);
  for (size_t a = 1; a < values.x.size(); a++)
    for (size_t b = 1; b < values.y.size(); b++)
      pairs.push_back({values.x[a], values.y[b], 0});
  exact = true;
  return pairs;
}

// Adaptive ER reference for simulations. Exact below the pair threshold,
// otherwise a deterministic, high-resolution distinct-rank EQ grid whose
// number of evaluated pairs is <= er_pairs. No randomized cutpoint subsampling.
// The underlying rank-averaged ER estimand remains the full exhaustive sum.
// NHANES has a separate always-exact ER path in nhanes_survey_core.hpp.
// Fixed-split sandwich averages cluster influences BEFORE squaring.
inline ErResult estimate_er(const std::vector<Patient>& patients,
                            bool complete,
                            int min_n,
                            int er_pairs,
                            int exact_limit,
                            uint64_t /*seed*/,
                            bool compute_se) {
  const auto values = rank_values(patients, complete, min_n);
  const uint64_t full_Q = uint64_t(values.x.size() - 1) * uint64_t(values.y.size() - 1);
  const bool exact = full_Q <= uint64_t(exact_limit);
  std::vector<double> xc, yc;
  int reference_K = 0;
  if (exact) {
    xc.assign(values.x.begin() + 1, values.x.end());
    yc.assign(values.y.begin() + 1, values.y.end());
  } else {
    // Use the finest scalar rank grid whose actual rectangular pair count
    // fits the pair budget. This also uses the budget efficiently when
    // D_X and D_Y are unequal (so K_X and K_Y cap at different values).
    const int DX = int(values.x.size()), DY = int(values.y.size());
    int lo = 2, hi = std::max(DX, DY);
    while (lo < hi) {
      const int mid = lo + (hi - lo + 1) / 2;
      const uint64_t pairs = uint64_t(std::min(mid, DX) - 1) * uint64_t(std::min(mid, DY) - 1);
      if (pairs <= uint64_t(er_pairs))
        lo = mid;
      else
        hi = mid - 1;
    }
    reference_K = lo;
    const auto grid = empirical_cuts(patients, reference_K, complete, min_n);
    xc = grid.x;
    yc = grid.y;
  }
  const uint64_t Q = uint64_t(xc.size()) * uint64_t(yc.size());
  if (Q < 1 || Q > uint64_t(std::numeric_limits<int>::max()))
    throw std::length_error("ER reference split count invalid");
  if (!exact && Q > uint64_t(er_pairs))
    throw std::logic_error("ER reference exceeds specified pair budget");
  std::vector<RankIndex> indices;
  indices.reserve(patients.size());
  for (const Patient& p : patients) {
    if (!complete && p.nobs < min_n)
      continue;
    indices.emplace_back(p, complete);
  }
  const int M = int(indices.size());
  if (M < 2)
    throw std::runtime_error("Insufficient subjects for ER");
  double rho_sum = 0;
  std::vector<double> influences;
  if (compute_se)
    influences.assign(M, 0.0);
  for (double xx : xc)
    for (double yy : yc) {
      Mom mean;
      for (const auto& idx : indices)
        mean += idx.query(xx, yy);
      mean *= 1.0 / double(M);
      const double r = pearson(mean);
      if (!std::isfinite(r))
        throw std::runtime_error("ER reference undefined: zero variance");
      rho_sum += r;
      if (compute_se) {
        const Mom g = gradient(mean);
        for (int i = 0; i < M; i++)
          influences[i] += dot(g, indices[i].query(xx, yy) - mean);
      }
    }
  ErResult out;
  out.e.estimate = rho_sum / double(Q);
  out.pairs_used = int(Q);
  out.reference_K = reference_K;
  out.exact = exact;
  // Numerical approximation error is not estimated by the sandwich SE.
  out.numerical_se = exact ? 0.0 : std::numeric_limits<double>::quiet_NaN();
  if (compute_se) {
    double sq = 0;
    for (double u : influences) {
      u /= double(Q);
      sq += u * u;
    }
    out.e.se = std::sqrt(std::max(0.0, sq / (double(M) * (M - 1))));
  }
  return out;
}

inline SimOutput run_pilot(const Spec& s) {
  validate(s);
  Rng rng(s.seed);
  SimOutput out;
  out.rows.reserve(size_t(s.reps) * s.cuts.size() * NMETHODS);
  for (int rep = 1; rep <= s.reps; rep++) {
    std::vector<Patient> patients;
    patients.reserve(s.M);
    Mom all_full;
    int included = 0, retained = 0;
    for (int i = 0; i < s.M; i++) {
      patients.push_back(make_patient(rng, s));
      const Patient& p = patients.back();
      retained += p.nobs;
      included += int(p.nobs >= s.min_n);
      for (const Tooth& t : p.teeth)
        all_full += moments(t.x, t.y);
    }
    if (included < 2)
      throw std::runtime_error("Fewer than two retained patients");
    const double full_rho = pearson(all_full / double(s.M * s.nmax));
    const uint64_t erseed = s.seed + uint64_t(rep) * 0x9E3779B97F4A7C15ULL;
    const ErResult er =
        estimate_er(patients, false, s.min_n, s.er_pairs, s.er_exact_limit, erseed, true);
    for (const auto& cut : s.cuts) {
      const int K = int(cut.size()) + 1;
      const Cuts fixed = fixed_cuts(cut);
      const Cuts empirical = empirical_cuts(patients, K, false, s.min_n);
      std::vector<Feature> fixed_features, emp_features;
      fixed_features.reserve(included);
      emp_features.reserve(included);
      double mean_n = 0, mean_pairs = 0;
      for (const Patient& p : patients)
        if (p.nobs >= s.min_n) {
          Feature a = features(p, fixed, false), b = features(p, empirical, false);
          mean_n += a.n;
          mean_pairs += a.n_pairs;
          fixed_features.push_back(std::move(a));
          emp_features.push_back(std::move(b));
        }
      for (int method = 0; method < NMETHODS; method++) {
        Estimate e;
        if (method == OPW)
          continue; // OPW is internal to EQ and ER, not a comparator.
        if (method == ER)
          e = er.e;
        else
          e = evaluate(method == EQ ? emp_features : fixed_features, method, K);
        // All methods use the within-replicate fixed-cutpoint cluster-sandwich
        // delta SE. For EQ and ER, influences are averaged over their selected
        // splits *before* computing the sandwich variance.
        if (!std::isfinite(e.estimate) || !std::isfinite(e.se) || e.se < 0)
          throw std::runtime_error("Non-finite estimate or sandwich SE");
        SimRow row;
        row.rep = rep;
        row.K = K;
        row.method = method;
        row.included = included;
        row.excluded = s.M - included;
        row.estimate = e.estimate;
        row.se = e.se;
        row.se_fixed = e.se;
        row.ci_lower = e.estimate - 1.959963984540054 * e.se;
        row.ci_upper = e.estimate + 1.959963984540054 * e.se;
        row.eq_boot_B = 0;
        row.fullrep = full_rho;
        row.mean_n = mean_n / included;
        row.retained_fraction = double(retained) / (s.M * s.nmax);
        row.mean_pairs = mean_pairs / included;
        if (method == ER) {
          row.numerical_se = er.numerical_se;
          row.er_pairs_used = er.pairs_used;
          row.er_reference_K = er.reference_K;
          row.er_exact = er.exact;
        }
        out.rows.push_back(row);
      }
    }
  }
  return out;
}

inline RefOutput run_reference(const Spec& s) {
  validate(s);
  Rng rng(s.seed);
  std::vector<Patient> patients;
  patients.reserve(s.reference_clusters);
  int included = 0, retained = 0;
  for (int i = 0; i < s.reference_clusters; i++) {
    patients.push_back(make_patient(rng, s));
    included += int(patients.back().nobs >= s.min_n);
    retained += patients.back().nobs;
  }
  if (included < 8)
    throw std::runtime_error("Insufficient retained reference clusters");
  // Full-population ER reference would require roughly (reference_clusters*nmax)^2
  // distinct splits. It is NOT approximated: report NA for ER-specific
  // reference targets instead of silently substituting QMC. The simulation
  // itself always computes exact ER for every replicate.
  constexpr int B = 4;
  RefOutput out;
  out.rows.reserve(s.cuts.size() * NMETHODS);
  for (const auto& cut : s.cuts) {
    const int K = int(cut.size()) + 1, Q = (K - 1) * (K - 1);
    const Cuts fixed = fixed_cuts(cut), emp_full = empirical_cuts(patients, K, true, s.min_n),
               emp_observed = empirical_cuts(patients, K, false, s.min_n);
    std::array<std::array<Accum, 4>, B> batches;
    for (int b = 0; b < B; b++)
      for (int mode = 0; mode < 4; mode++)
        batches[b][mode] = Accum(Q);
    std::array<Accum, 4> all = {Accum(Q), Accum(Q), Accum(Q), Accum(Q)};
    for (int i = 0; i < s.reference_clusters; i++) {
      const Patient& p = patients[i];
      int b = i % B;
      Feature fixed_full = features(p, fixed, true);
      Feature rank_full = features(p, emp_full, true);
      all[0].add(fixed_full);
      batches[b][0].add(fixed_full);
      all[2].add(rank_full);
      batches[b][2].add(rank_full);
      if (p.nobs >= s.min_n) {
        Feature fixed_obs = features(p, fixed, false);
        Feature rank_obs = features(p, emp_observed, false);
        all[1].add(fixed_obs);
        batches[b][1].add(fixed_obs);
        all[3].add(rank_obs);
        batches[b][3].add(rank_obs);
      }
    }
    auto batch_se = [](const std::array<double, B>& a) {
      double mean = std::accumulate(a.begin(), a.end(), 0.0) / double(B), sq = 0;
      for (double x : a)
        sq += (x - mean) * (x - mean);
      return std::sqrt(sq / (double(B) * (B - 1)));
    };
    for (int method = 0; method < NMETHODS; method++) {
      if (method == OPW)
        continue;
      RefRow r;
      r.K = K;
      r.method = method;
      r.included = included;
      if (method == ER) {
        const double na = std::numeric_limits<double>::quiet_NaN();
        r.target_full = na;
        r.target_observed = na;
        r.se_full = na;
        r.se_observed = na;
        r.numerical_se_full = 0;
        r.numerical_se_observed = 0;
        r.er_pairs_used = 0;
      } else {
        const int idx = is_empirical(method) ? 2 : 0;
        std::array<double, B> full{}, obs{};
        for (int b = 0; b < B; b++) {
          full[b] = batches[b][idx].get(method, K);
          obs[b] = batches[b][idx + 1].get(method, K);
        }
        r.target_full = method <= CW ? rho0(s) : all[idx].get(method, K);
        r.se_full = method <= CW ? 0 : batch_se(full);
        r.target_observed = all[idx + 1].get(method, K);
        r.se_observed = batch_se(obs);
      }
      r.retained_fraction = double(retained) / (s.reference_clusters * s.nmax);
      r.mean_n = all[1].total / included;
      r.excluded_fraction = 1.0 - double(included) / s.reference_clusters;
      r.mean_pairs = all[1].pairs / included;
      if (method != ER && (!std::isfinite(r.target_full) || !std::isfinite(r.target_observed)))
        throw std::runtime_error("Non-finite reference target");
      out.rows.push_back(r);
    }
  }
  return out;
}
} // namespace opwsim
#endif
