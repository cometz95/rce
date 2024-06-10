// C/C++
#include <algorithm>
#include <cmath>

// canoe
#include <air_parcel.hpp>

// snap
#include "thermodynamics.hpp"

// Calculates phase equilibrium of
// Vapor <=> Precipitate (surface)
//
// Example phase equilibrium:
// H2O -> H2O(l)
// this function follows the notation of Hartmann ch. 4 and 5 (see eq 5.12
// and 4.33) Hartmann, Global Physical Climatology
//
// i is defined by the position of the vapor in the airparcel (for (int i = 1; i
// <= NVAPOR; ++i) )
RealArrayX Thermodynamics::TryEquilibriumTP_VaporPrecip(
    AirParcel const& qfrac, int i, Real layerSclFact, Real& amd, Real btemp,
    Real dTs, Real cSurf, Real dt, Real Cde) const {
  Real xv = qfrac.w[i], xg = 1. - xv;

  // need to pass an airparcel object to svp_func1 to avoid making a new svp
  // function which takes T directly
  AirParcel btemp_container;
  btemp_container[IDN] = btemp;

  std::vector<Real> rates(1 + cloud_index_set_[i].size(), 0.);

#pragma omp simd reduction(+ : xg)
  for (int n = 0; n < NCLOUD; ++n) xg += -qfrac.c[n];

  for (int n = 0; n < cloud_index_set_[i].size(); ++n) {
    int j = cloud_index_set_[i][n];
    Real es_surf = svp_func1_[i][n](btemp_container, i, j);

    // getting des/dt locally around btemp with finite difference
    Real T2 = btemp + 0.1;
    btemp_container[IDN] = T2;
    Real es_2 = svp_func1_[i][n](btemp_container, i, j);
    Real T1 = btemp - 0.1;
    btemp_container[IDN] = T1;
    Real es_1 = svp_func1_[i][n](btemp_container, i, j);
    Real des_dt = (es_2 - es_1) / (T2 - T1);

    Real Mi = GetMu(i);
    Real Mbar = 0;
    for (int k = 1; k <= NVAPOR; ++k) Mbar += GetMu(k) * qfrac.w[k];

    // this comes from the defintion qstar = (es/p)*(Mi/Mbar), and a bit of
    // calculus evaluated at the surface
    Real dqstar_dt = (Mi / (qfrac.w[IPR] * Mbar)) * (des_dt - es_surf / btemp);

    Real L = GetLatentEnergyMass(i, btemp);
    Real Be = (GetCpMassRef(i) / L) / dqstar_dt;
    Real rhoatm = (qfrac.w[IPR] * Mbar) / (Constants::Rgas * qfrac.w[IDN]);
    // Real En = (cSurf/dt)*(dTs - dTs) / L;
    Real En = 0;
    Real Ubar =
        std::sqrt(std::pow(qfrac.w[IVX], 2) + std::pow(qfrac.w[IVY], 2) +
                  std::pow(qfrac.w[IVZ], 2));
    Real Eair = rhoatm * Cde * Ubar *
                ((svp_func1_[i][n](qfrac, i, j) / qfrac.w[IPR]) * (Mi / Mbar) -
                 qfrac.w[i]);
    // we calculate the evaporation rate, so we subtract the rate from amd and
    // add it to the vapor pool
    Real rate = En / (1 + Be) + Be * Eair / (1 + Be);  // kg/m^2/s

    if (amd > 0) {
      // if rate is negative, condensation occurs
      //  condense at most xv/layerSclFact vapor
      if (rate < 0.) {
        // only condensing from a ~1m layer of vapor near surf, so scale the
        // rate/total amount of vapor by layerSclFact

        // CONVERT xv/layerSclFact TO kg/m^2
        rates[0] += -std::min(-rate * dt, xv / layerSclFact);
        rates[1 + n] = std::min(-rate * dt, xv / layerSclFact);
      }
      // if rate is positive, evaporation occurs
      //  evaporate at most xc precip
      if (rate > 0.) {
        rates[0] += std::min(rate * dt, amd);
        rates[1 + n] = -std::min(rate * dt, amd);
      }
    } else {
      amd = 0;
      if (rate > 0.) {
        rates[0] = 0;
        rates[1 + n] = 0;
      } else {
        // CONVERT xv/layerSclFact TO kg/m^2
        rates[0] += -std::min(-rate * dt, xv / layerSclFact);
        rates[1 + n] = std::min(-rate * dt, xv / layerSclFact);
      }
    }
  }

  // scale total rate
  if (rates[0] < 0. && std::abs(rates[0]) > xv) {
    Real r = xv / std::abs(rates[0]);
    for (auto& rate : rates) rate *= r;
  }

  // subtract from amd
  return rates;
  //        //CONVERT rates to mol/mol
  // qfrac.w[i] += rates[0];
  // if (qfrac.w[i] < 0) {
  //  qfrac.w[i] = 0;
  //}
}
