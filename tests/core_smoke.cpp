#include "../src/core.h"
#include <cassert>
#include <iostream>

using namespace opwsim;

int main() {
  Spec s;
  s.M = 25;
  s.nmax = 10;
  s.rho_uv = .25;
  s.eta = 2;
  s.seed = 12345;
  s.sigma2 = .5;
  s.tau2 = .5;
  Rng rng(s.seed);
  std::vector<Patient> p;
  for (int i = 0; i < s.M; i++) {
    auto person = make_patient(rng, s);
    if (person.nobs)
      p.push_back(std::move(person));
  }
  Cuts c = empirical_cuts(p, 20, false, 1);
  std::vector<Feature> ff;
  for (const auto& subj : p)
    ff.push_back(features(subj, c, false));
  for (int method : {NW, CW, EQ}) {
    Estimate z = evaluate(ff, method, 20);
    assert(std::isfinite(z.estimate) && std::isfinite(z.se));
  }
  auto er = estimate_er(p, false, 1, 4096, 4096, 0, true);
  assert(std::isfinite(er.e.estimate) && std::isfinite(er.e.se));
  std::cout << "PASS core: NW/CW/EQ/ER estimates + sandwich SE" << std::endl;
}
