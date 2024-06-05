// C/C++
#include <algorithm>

// canoe
#include <air_parcel.hpp>

// snap
#include "thermodynamics.hpp"

// Calculates phase equilibrium of
// Vapor <=> Precipitate (surface)
//
// Example phase equilibrium:
// H2O -> H2O(l)
//
RealArrayX Thermodynamics::TryEquilibriumTP_VaporCloud(
    AirParcel const& qfrac, int i, Real cv_hat, bool misty, Real layerSclFact,
    Real amd, Real btemp) const {
  Real xv = qfrac.w[i], xg = 1. - xv;
  Real t = btemp;
  std::vector<Real> rates(1 + cloud_index_set_[i].size(), 0.);
  Real xp = amd * ? ? ? ;  // precipitate mole fraction

#pragma omp simd reduction(+ : xg)
  for (int n = 0; n < NCLOUD; ++n) xg += -qfrac.c[n];

  for (int n = 0; n < cloud_index_set_[i].size(); ++n) {
    int j = cloud_index_set_[i][n];
    Real xs = svp_func1_[i][n](qfrac, i, j) / qfrac.w[IPR];
    ? ? ?  // need to pass btemp into svp_func1_ rather than qfrac
        Real xc = qfrac.c[j];

    if (misty) {  // in a cloudy ambient environment
      rates[0] += xs - xv / (xg + xv);

      // tried this as a potential fix to negative water when T very low, did
      // not seem to work
      // if (rates[0] < 0.) {
      //   rates[0] += -std::min(-rates[0], xv);
      //   rates[1 + n] = std::min(-rates[0], xv);
      // }
      continue;
    }

    // if saturation vapor pressure is larger than the total pressure
    // evaporate all condensates
    if (xs > 1.) {
      rates[0] += xc;
      rates[1 + n] = -xc;
      continue;
    }

    Real alpha = 0.;

    Real lv = beta_[1 + NVAPOR + j] / t - delta_[1 + NVAPOR + j];
    if (cv_hat > 0.) alpha = (lv - 1.) / cv_hat;

    Real s1 = xs / (1. - xs);
    Real rate = (s1 * xg - xv) / (1. + alpha * xg * lv * s1 / (1. - xs));

    // condensate at most xv vapor
    if (rate < 0.) {
      // only condensing from a ~1m layer of vapor near surf, so scale the
      // rate/total amount of vapor by layerSclFact
      rates[0] += -std::min(-rate / layerSclFact, xv / layerSclFact);
      rates[1 + n] = std::min(-rate / layerSclFact, xv / layerSclFact);
    }

    // evaporate at most xc precip
    if (rate > 0.) {
      rates[0] += std::min(rate, xp);
      rates[1 + n] = -std::min(rate, xp);
    }
  }

  // scale total rate
  if (rates[0] < 0. && std::abs(rates[0]) > xv) {
    Real r = xv / std::abs(rates[0]);
    for (auto& rate : rates) rate *= r;
  }

  return rates;
}