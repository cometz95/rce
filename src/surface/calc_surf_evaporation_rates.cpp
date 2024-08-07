// C/C++
#include <algorithm>
#include <cmath>

// canoe
#include <air_parcel.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// surface
#include "surface.hpp"

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

// for now, amd is a scalar and associated with one phase. later, make amd[n]
// for each phase

// qfrac has T in IDN slot
RealArrayX Surface::CalcSurfEvapRates(AirParcel const& qfrac, int i,
                                      Real& amd_x, Real btemp, Real cSurf,
                                      Real dt, Real Cde, Real Mbar) const {
  auto pthermo = Thermodynamics::GetInstance();
  Real xv = qfrac.w[i];

  // need to pass an airparcel object to svp_func1 to avoid making a new svp
  // function which takes T directly
  AirParcel btemp_container(AirParcel::Type::MoleFrac);
  btemp_container.w[IDN] = btemp;

  std::vector<Real> rates(1 + pthermo->GetCloudIndexSet(i).size(), 0.);

  for (int n = 0; n < pthermo->GetCloudIndexSet(i).size(); ++n) {
    int j = pthermo->GetCloudIndex(i, n);
    Real es_surf = pthermo->GetSVPFunc1()[i][n](btemp_container, i, j);

    // getting des/dt locally around btemp with definition of the derivative
    Real T2 = btemp + 0.01;
    btemp_container.w[IDN] = T2;
    Real es_2 = pthermo->GetSVPFunc1()[i][n](btemp_container, i, j);
    Real T1 = btemp - 0.01;
    btemp_container.w[IDN] = T1;
    Real es_1 = pthermo->GetSVPFunc1()[i][n](btemp_container, i, j);
    Real des_dt = (es_2 - es_1) / (T2 - T1);
    // CHECKME (cmetz)
    int thermoIndex = 1 + n + NVAPOR;
    Real Mi = pthermo->GetMu(thermoIndex);

    // this comes from the defintion qstar = (es/p)*(Mi/Mbar), and a bit of
    // calculus evaluated at the surface
    Real dqstar_dt = (Mi / (qfrac.w[IPR] * Mbar)) * (des_dt - es_surf / btemp);

    // needs to be 2 or 3 (NOT i, which I set to zero in inp file). each phase
    // has a different latent heat of vaporization
    Real L = pthermo->GetLatentEnergyMass(thermoIndex, btemp);
    Real Be = (pthermo->GetCpMassRef(thermoIndex) / L) / dqstar_dt;
    Real rhoatm = (qfrac.w[IPR] * Mbar) / (qfrac.w[IDN] * Constants::Rgas);
    // Real En = (cSurf/dt)*(dTs - dTs) / L;
    Real En = 0;  // with a different parameterization of G, this wouldn't be 0

    //(cmetz) check this, should the vertical componant of the wind be included
    // or only X and Y?
    Real Ubar =
        std::sqrt(std::pow(qfrac.w[IVX], 2) + std::pow(qfrac.w[IVY], 2) +
                  std::pow(qfrac.w[IVZ], 2));
    Real Eair = rhoatm * Cde * Ubar *
                ((pthermo->GetSVPFunc1()[i][n](qfrac, i, j) / qfrac.w[IPR]) *
                     (Mi / Mbar) -
                 qfrac.w[i]);
    // we calculate the evaporation rate, so we subtract the rate from amd and
    // add it to the vapor pool
    Real rate = En / (1 + Be) + Be * Eair / (1 + Be);  // kg/m^2/s
    // std::cout << "thermoIndex: " << thermoIndex << std::endl;
    // std::cout << "rate: " << rate << std::endl;
    // std::cout << "Be: " << Be << std::endl;
    // std::cout << "Eair: " << Eair << std::endl;
    // std::cout << "dqstar_dt: " << dqstar_dt << std::endl;
    // std:: cout<< "L: " << L << std::endl;
    // std::cout << "cp mass ref: " << GetCpMassRef(2) << std::endl;
    // std::cout << "Mbar: " << Mbar << std::endl;
    // std::cout << "des_dt: " << des_dt << std::endl;
    // std::cout << "es_surf: " << es_surf << std::endl;
    // std::cout << "Mi: " << Mi << std::endl;
    // std::cout << "qfrac.w[IDN]: " << qfrac.w[IDN] << std::endl;
    // std::cout << "btemp: " << btemp << std::endl;

    if (amd_x > 0) {
      // if rate is negative, condensation occurs. condense at most
      // xv/layerThickness vapor
      if (rate < 0.) {
        // can only condense from a 1m layer of vapor near surf, so scale the
        // rate/total amount of vapor in the box by layerThickness in unit
        // conversion, the layerThickness cancels out
        rates[0] += -std::min(-rate * dt, xv * rhoatm * (Mi / Mbar));
        rates[1 + n] = std::min(-rate * dt, xv * rhoatm * (Mi / Mbar));
      }
      // if rate is positive, evaporation occurs. can evaporate at most amd
      // precip
      else if (rate > 0.) {
        rates[0] += std::min(rate * dt, amd_x);
        rates[1 + n] = -std::min(rate * dt, amd_x);
      }
    } else if (amd_x <= 0) {
      amd_x = 0;
      // if rate is positive and amd = 0, no more evaporation can occur
      if (rate > 0.) {
        rates[0] = 0;
        rates[1 + n] = 0;
      } else if (rate <
                 0.) {  // rate is negative and condensation can still occur
        rates[0] += -std::min(-rate * dt, xv * rhoatm * (Mi / Mbar));
        rates[1 + n] = std::min(-rate * dt, xv * rhoatm * (Mi / Mbar));
      }
    }
  }

  // scale total rate
  // if (rates[0] < 0. && std::abs(rates[0]) > xv) {
  //  Real r = xv / std::abs(rates[0]);
  //  for (auto& rate : rates) rate *= r;
  //}

  return rates;
}
