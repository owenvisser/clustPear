#include "../src/core.h"
#include "../src/weighted_core.h"
#include <cassert>
#include <cmath>
#include <iostream>
using namespace opwsim;
namespace cp = clustpear_weighted;
int main() {
  Spec s;
  s.M = 12; s.nmax = 7; s.rho_uv = .4; s.eta = 1.2;
  s.seed = 415; s.sigma2 = .3; s.tau2 = .4;
  Rng rng(s.seed);
  std::vector<Patient> p;
  std::vector<double> w;
  for (int i = 0; i < s.M; i++) {
    auto person = make_patient(rng, s);
    if (person.nobs) {
      p.push_back(std::move(person));
      w.push_back(1.0 + 0.25 * i);
    }
  }
  std::vector<RankIndex> idx;
  for (const auto& person : p) idx.emplace_back(person, false);
  auto eq = cp::eq_pairs(p, 4);
  bool exact = false;
  auto er = cp::er_pairs(p, 100, exact);
  for (int type = 0; type < 4; type++) {
    std::vector<double> influence;
    auto fun = [&](const std::vector<double>& ww, std::vector<double>* inf) {
      if (type == 0 || type == 1)
        return cp::moment_estimator(p, ww, type == 1, inf);
      return cp::pair_average(idx, ww, type == 2 ? eq : er, inf);
    };
    const auto base = fun(w, &influence);
    assert(std::isfinite(base.estimate));
    assert(std::isfinite(base.se));
    assert(influence.size() == w.size());
    double score = 0.0;
    for (size_t i = 0; i < w.size(); i++) score += w[i] * influence[i];
    assert(std::abs(score) < 1e-9);
    for (size_t i = 0; i < w.size(); i++) {
      auto ww = w; ww[i] += 1e-5;
      const double numerical = (fun(ww, nullptr).estimate - base.estimate) / 1e-5;
      if (std::abs(numerical - influence[i]) > 1e-5) {
        std::cerr << "Derivative failure: method=" << type << ", i=" << i
                  << ", analytical=" << influence[i] << ", numerical=" << numerical << "\n";
        return 1;
      }
    }
    std::cout << "PASS method=" << type << " weighted Taylor derivative\n";
  }
}
