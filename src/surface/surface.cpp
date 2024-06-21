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

  alpha_s = pin->GetOrAddReal("surface", "alpha_s", 0.3);
  alpha_a = pin->GetOrAddReal("surface", "alpha_a", 0.5);
  cSurf = pin->GetOrAddReal("surface", "csurf", 100000);
  rho_l_vapor1 = pin->GetOrAddReal("surface", "rho_l_vapor1", 1000);
  rho_s_vapor1 = pin->GetOrAddReal("surface", "rho_s_vapor1", 910);
  omega = pin->GetOrAddReal("surface", "omega", 0.00007272205);
  s0 = pin->GetOrAddReal("surface", "s0", 1360);
  // FIXME (cmetz) should be calculated from thermodynamics, not hard coded
  meltingPointVapor1 = 273.15;

  btempArray.NewAthenaArray(pmb->ncells3, pmb->ncells2);
  amd.NewAthenaArray(NVAPOR, 2, pmb->ncells3, pmb->ncells2);
  gel.NewAthenaArray(NVAPOR, 2, pmb->ncells3, pmb->ncells2);

  for (int j = js; j <= je; ++j) {
    for (int k = ks; k <= ke; ++k) {
      btempArray(k, j) = pin->GetOrAddReal("surface", "init_btemp", 0);
      amd(0, 0, k, j) = pin->GetOrAddReal("surface", "init_amd_s_vapor1", 0);
      amd(0, 1, k, j) = pin->GetOrAddReal("surface", "init_amd_l_vapor1", 0);
      gel(0, 0, k, j) = amd(0, 0, k, j) / rho_s_vapor1;
      gel(0, 1, k, j) = amd(0, 1, k, j) / rho_l_vapor1;
    }
  }
}

Surface::~Surface() {
  Application::Logger app("surface");
  app->Log("Destroy Surface");
}

void Surface::DoSurfaceProcesses(MeshBlock *pmb) {
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      Real dt = pmb->pmy_mesh->dt;
      Real accumPrecipAmd =
          pmb->pimpl->psurf->AccumulatePrecipitates(pmb, k, j, is);
      RealArrayX AmdEvap = pmb->pimpl->psurf->EvapPrecip(pmb, k, j, dt);
      Real dTs = pmb->pimpl->psurf->ChangeTempFromForcing(
          pmb, k, j, dt, accumPrecipAmd, AmdEvap);
    }
  }
}

Real Surface::ChangeTempFromForcing(MeshBlock *pmb, int k, int j, Real dt,
                                    Real accumPrecipAmd, RealArrayX AmdEvap) {
  Real swin = s0 * (1 + std::sin(omega * time));
  auto pthermo = Thermodynamics::GetInstance();

  // get the flux from each band and add it up
  Real tot_fluxdn = 0;
  int numBands = pmb->pimpl->prad->GetNumBands();

  for (int n = 0; n < numBands; ++n)
    tot_fluxdn += pmb->pimpl->prad->GetBand(n)->bflxdn(k, j, is);

  int sign = 0;
  // FIXME (cmetz) should loop over all precipitates, not just 0
  // e.g., int thermoIndex = 1 + nCloud + NVAPOR;
  // FIXME (cmetz) maybe later incorporate the q=mcat term (fix sign issue), for
  // now just include latent
  int thermoIndex = 1 + 0 + NVAPOR;
  Real cp = pthermo->GetCpMassRef(thermoIndex);
  Real L = pthermo->GetLatentEnergyMass(thermoIndex, btempArray(k, j));
  Real Mbar = pthermo->GetMu(pmb, k, j, is);
  // Real Tip = (pmb->phydro->w(IPR, k, j, is) * Mbar) /
  //            (pmb->phydro->w(IDN, k, j, is) * Constants::Rgas);
  //  FIXME (cmetz) are these signs right?
  // Real Tf = (cSurf * btempArray(k, j) - accumPrecipAmd * cp * Tip +
  //            sign * accumPrecipAmd * L) /
  //           (cSurf - accumPrecipAmd * cp);
  if (btempArray(k, j) < meltingPointVapor1)
    sign = 1;
  else if (btempArray(k, j) > meltingPointVapor1)
    sign = 0;
  Real dT_latent_precip = sign * (accumPrecipAmd * L) / cSurf;

  // Tf = (cSurf * btempArray(k, j) - accumPrecipAmd * cp * Tip +
  //       sign * accumPrecipAmd * L) /
  //      (cSurf - accumPrecipAmd * cp);

  // FIXME (cmetz) should account for melting and vaporization of solids, not
  // just vaporization.
  // right now don't have access to L for solid phase
  // AmdEvap is the amount evaporated off of the surface, and its energy
  // leaving the surface
  Real dT_evapCooling = (AmdEvap[0] * L + AmdEvap[1] * L) / cSurf;

  Real dTs = (swin * (1 - alpha_a) * (1 - alpha_s) + tot_fluxdn -
              Constants::stefanBoltzmann * pow(btempArray(k, j), 4)) *
                 (dt / cSurf) +
             dT_latent_precip - dT_evapCooling;

  btempArray(k, j) += dTs;
  return dTs;
}

Real Surface::AccumulatePrecipitates(MeshBlock *pmb, int k, int j, int iSkim) {
  Real precip = 0;
  Real accumPrecipAmd = 0.;

  for (int i = is; i <= iSkim; ++i) {
    precip = 0;
    for (int n = NCLOUD / 2; n < NCLOUD; ++n) {
      precip = pmb->pscalars->r(n, k, j, i);
      pmb->pscalars->r(n, k, j, i) = 0;
      accumPrecipAmd = precip * pmb->phydro->w(IDN, k, j, is) * dzPBL;
      // FIXME (cmetz) should loop over all precipitates, not just H2O, and
      // return an array of accum precipitates
      if (n == 1) {
        if (btempArray(k, j) > meltingPointVapor1)
          amd(0, 1, k, j) += accumPrecipAmd;
        else
          amd(0, 0, k, j) += accumPrecipAmd;
      }
    }
  }

  gel(0, 1, k, j) = amd(0, 1, k, j) / rho_l_vapor1;
  gel(0, 0, k, j) = amd(0, 0, k, j) / rho_s_vapor1;

  return accumPrecipAmd;
}

RealArrayX Surface::EvapPrecip(MeshBlock *pmb, int k, int j, Real dt) {
  auto pthermo = Thermodynamics::GetInstance();

  // get the current air parcel
  AirParcel air(AirParcel::Type::MoleFrac);
#pragma omp simd
  for (int n = 0; n < NHYDRO; ++n) air.w[n] = pmb->phydro->w(n, k, j, is);
#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n)
    air.c[n] = pmb->pimpl->pmicro->u(n, k, j, is);

  Real Mbar = pthermo->GetMu(pmb, k, j, is);
  // make sure air has T in IDN slot
  air.w[IDN] = (air.w[IPR] * Mbar) / (air.w[IDN] * Constants::Rgas);

  std::vector<Real> AmdEvap(2, 0.);

  for (int n = 1; n <= NVAPOR; ++n) {
    std::vector<Real> rates(1 + pthermo->GetCloudIndexSet(n).size(), 0.);
    //(cmetz) check this value of CDE, see Hartmann pg 117
    rates = this->CalcSurfEvapRates(air, n, amd(0, 1, k, j), btempArray(k, j),
                                    cSurf, dt, 3e-3, Mbar);
    // std::cout << "liquid evap rates: " << rates[0] << std::endl;

    // evaporate off of liquid amd
    amd(0, 1, k, j) += -rates[0];
    AmdEvap[1] = rates[0];

    // add or subtract from vapor pool
    pmb->phydro->w(n, k, j, is) += rates[0] * (Mbar / pthermo->GetMu(n)) /
                                   (pmb->phydro->w(IDN, k, j, is) * dzPBL);
    if (pmb->phydro->w(n, k, j, is) < 0) pmb->phydro->w(n, k, j, is) = 0;

    rates = this->CalcSurfEvapRates(air, n, amd(0, 0, k, j), btempArray(k, j),
                                    cSurf, dt, 3e-3, Mbar);
    // std::cout << "solid evap rates: " << rates[0] << std::endl;
    // evaporate off of solid amd
    amd(0, 0, k, j) += -rates[0];
    AmdEvap[0] = rates[0];

    // add or subtract from vapor pool
    pmb->phydro->w(n, k, j, is) += rates[0] * (Mbar / pthermo->GetMu(n)) /
                                   (pmb->phydro->w(IDN, k, j, is) * dzPBL);
    if (pmb->phydro->w(n, k, j, is) < 0) pmb->phydro->w(n, k, j, is) = 0;
  }

  gel(0, 1, k, j) = amd(0, 1, k, j) / rho_l_vapor1;
  gel(0, 0, k, j) = amd(0, 0, k, j) / rho_s_vapor1;

  return AmdEvap;
}
