// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <impl.hpp>

// microphysics
#include <microphysics/microphysics.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation.hpp>

// surface
#include "surface.hpp"

Surface::Surface(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("surface");
  app->Log("Initialize Surface");

  is = pmb->is, js = pmb->js, ks = pmb->ks;
  ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  time = pmb->pmy_mesh->time;
  dzPBL = pmb->pcoord->x1f(is + 1) - pmb->pcoord->x1f(is);

  alpha_s = pin->GetReal("surface", "alpha_s");
  alpha_a = pin->GetReal("surface", "alpha_a");
  cSurf = pin->GetReal("surface", "csurf");
  rholH2O = pin->GetReal("surface", "rholH2O");
  rhosH2O = pin->GetReal("surface", "rhosH2O");
  omega = pin->GetReal("surface", "omega");
  s0 = pin->GetReal("surface", "s0");

  btempArray.NewAthenaArray(pmb->ncells2);
  amd.resize(NVAPOR);
  gel.resize(NVAPOR);

  for (int i = 0; i < NVAPOR; ++i) {
    // there are only 2 phases which can be stored on the surface
    amd[i].resize(2);
    gel[i].resize(2);
    for (int n = 0; n <= 1; ++n) {
      amd[i][n].NewAthenaArray(pmb->ncells2);
      gel[i][n].NewAthenaArray(pmb->ncells2);
    }
  }

  for (int j = js; j <= je; ++j) {
    btempArray(j) = pin->GetReal("surface", "init_btemp");
    amd[0][0](j) = pin->GetReal("surface", "init_amd_s_vapor1");
    amd[0][1](j) = pin->GetReal("surface", "init_amd_l_vapor1");
    gel[0][0](j) = amd[0][0](j) / rhosH2O;
    gel[0][1](j) = amd[0][1](j) / rholH2O;
  }
}

Surface::~Surface() {
  Application::Logger app("surface");
  app->Log("Destroy Surface");
}

Real Surface::ChangeTempFromForcing(MeshBlock *pmb, int j, Real dt,
                                    Real accumPrecipAmd, RealArrayX AmdEvap) {
  double swin = s0 * (1 + std::sin(omega * time));
  auto pthermo = Thermodynamics::GetInstance();

  // get the flux from each band and add it up
  double tot_fluxdn = 0;
  int numBands = pmb->pimpl->prad->GetNumBands();

  for (int n = 0; n < numBands; ++n)
    tot_fluxdn += pmb->pimpl->prad->GetBand(n)->bflxdn(ks, j, is);

  int sign = 0;
  // FIXME (cmetz) should loop over all precipitates, not just 0
  int thermoIndex = 1 + 0 + NVAPOR;
  Real cp = pthermo->GetCpMassRef(thermoIndex);
  Real L = pthermo->GetLatentEnergyMass(thermoIndex, btempArray(j));
  Real Mbar = pthermo->GetMu(pmb, ks, j, is);
  Real Tip = (pmb->phydro->w(IPR, ks, j, is) * Mbar) /
             (pmb->phydro->w(IDN, ks, j, is) * Constants::Rgas);
  // FIXME (cmetz) are these signs right?
  Real Tf = (cSurf * btempArray(j) - accumPrecipAmd * cp * Tip +
             sign * accumPrecipAmd * L) /
            (cSurf - accumPrecipAmd * cp);
  std::cout << "accumPrecipAmd: " << accumPrecipAmd << std::endl;
  std::cout << "(guess) Tf - Tis: " << Tf - btempArray(j) << std::endl;
  if (Tf < 273.)
    sign = 1;
  else if (Tf > 273.)
    sign = 0;
  Tf = (cSurf * btempArray(j) - accumPrecipAmd * cp * Tip +
        sign * accumPrecipAmd * L) /
       (cSurf - accumPrecipAmd * cp);
  std::cout << "Tf - Tis: " << Tf - btempArray(j) << std::endl;

  // FIXME (cmetz) should account for melting and vaporization of solids, not
  // just vaporization.
  //  right now don't have access to L for solid phase
  //  AmdEvap is the amount evaporated off of the surface, and its energy
  //  leaving the surface
  Real evapCooling = (AmdEvap[0] * L + AmdEvap[1] * L) / dt;

  double dTs =
      (swin * (1 - alpha_a) * (1 - alpha_s) + tot_fluxdn -
       Constants::stefanBoltzmann * pow(btempArray(j), 4) - evapCooling) *
      (dt / cSurf);

  btempArray(j) = Tf + dTs;
  return dTs;
}

Real Surface::AccumulatePrecipitates(MeshBlock *pmb, int iSkim, int j) {
  double precip = 0;
  Real accumPrecipAmd = 0.;

  if (btempArray(j) > 273)
    H2OisLiquid = true;
  else
    H2OisLiquid = false;

  for (int i = is; i <= iSkim; ++i) {
    precip = 0;
    for (int n = NCLOUD / 2; n < NCLOUD; ++n) {
      precip = pmb->pscalars->r(n, ks, j, i);
      pmb->pscalars->r(n, ks, j, i) = 0;
      accumPrecipAmd = precip * pmb->phydro->w(IDN, ks, j, is) * dzPBL;
      // FIXME (cmetz) should loop over all precipitates, not just H2O, and
      // return an array of accum precipitates
      if (n == 1) {
        if (H2OisLiquid)
          amd[0][1](j) += accumPrecipAmd;
        else
          amd[0][0](j) += accumPrecipAmd;
      }
    }
  }

  gel[0][1](j) = amd[0][1](j) / rholH2O;
  gel[0][0](j) = amd[0][0](j) / rhosH2O;

  return accumPrecipAmd;
}

RealArrayX Surface::EvapPrecip(MeshBlock *pmb, int j, Real dt) {
  auto pthermo = Thermodynamics::GetInstance();

  // get the current air parcel
  AirParcel air(AirParcel::Type::MoleFrac);
#pragma omp simd
  for (int n = 0; n < NHYDRO; ++n) air.w[n] = pmb->phydro->w(n, ks, j, is);
#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n)
    air.c[n] = pmb->pimpl->pmicro->u(n, ks, j, is);

  Real Mbar = pthermo->GetMu(pmb, ks, j, is);
  // make sure air has T in IDN slot
  air.w[IDN] = (air.w[IPR] * Mbar) / (air.w[IDN] * Constants::Rgas);

  std::vector<Real> AmdEvap(2, 0.);

  for (int n = 1; n <= NVAPOR; ++n) {
    std::vector<Real> rates(1 + pthermo->GetCloudIndexSet(n).size(), 0.);
    //(cmetz) check this value of CDE, see Hartmann pg 117
    rates = this->CalcSurfEvapRates(air, n, amd[0][1](j), btempArray(j), cSurf,
                                    dt, 3e-3, Mbar);
    // std::cout << "liquid evap rates: " << rates[0] << std::endl;

    // evaporate off of liquid amd
    amd[0][1](j) += -rates[0];
    AmdEvap[1] = rates[0];

    // add or subtract from vapor pool
    pmb->phydro->w(n, ks, j, is) += rates[0] * (Mbar / pthermo->GetMu(n)) /
                                    (pmb->phydro->w(IDN, ks, j, is) * dzPBL);
    if (pmb->phydro->w(n, ks, j, is) < 0) pmb->phydro->w(n, ks, j, is) = 0;

    rates = this->CalcSurfEvapRates(air, n, amd[0][0](j), btempArray(j), cSurf,
                                    dt, 3e-3, Mbar);
    // std::cout << "solid evap rates: " << rates[0] << std::endl;
    // evaporate off of solid amd
    amd[0][0](j) += -rates[0];
    AmdEvap[0] = rates[0];

    // add or subtract from vapor pool
    pmb->phydro->w(n, ks, j, is) += rates[0] * (Mbar / pthermo->GetMu(n)) /
                                    (pmb->phydro->w(IDN, ks, j, is) * dzPBL);
    if (pmb->phydro->w(n, ks, j, is) < 0) pmb->phydro->w(n, ks, j, is) = 0;
  }

  gel[0][1](j) = amd[0][1](j) / rholH2O;
  gel[0][0](j) = amd[0][0](j) / rhosH2O;

  return AmdEvap;
}
