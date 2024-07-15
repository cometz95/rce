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
  // C = cp*rho*dz, where dz is the depth of the surface layer that changes
  // temperature
  Real cLand = pin->GetOrAddReal("surface", "cLand", 100000);
  Real cOcean = pin->GetOrAddReal("surface", "cOcean", 300000);
  Real fOcean = pin->GetOrAddReal("surface", "fOcean", 0.5);
  // weight C by the fraction of ocean and land
  cSurf = cLand * (1 - fOcean) + cOcean * fOcean;
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
  fin = 0;
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      Real dt = pmb->pmy_mesh->dt;
      RealArrayX AmdEvap = pmb->pimpl->psurf->EvapPrecip(pmb, k, j, dt);
      AthenaArray<Real> accumPrecipAmd =
          pmb->pimpl->psurf->AccumulatePrecipitates(pmb, k, j, is);
      Real dTs = pmb->pimpl->psurf->ChangeTempFromForcing(
          pmb, k, j, dt, accumPrecipAmd, AmdEvap);
    }
  }
}

Real Surface::ChangeTempFromForcing(MeshBlock *pmb, int k, int j, Real dt,
                                    AthenaArray<Real> &accumPrecipAmd,
                                    RealArrayX AmdEvap) {
  Real swin = s0 * (1 + std::sin(omega * time));
  auto pthermo = Thermodynamics::GetInstance();

  // get the flux from each band and add it up
  Real tot_fluxdn = 0;
  int numBands = pmb->pimpl->prad->GetNumBands();

  for (int n = 0; n < numBands; ++n)
    tot_fluxdn += pmb->pimpl->prad->GetBand(n)->bflxdn(k, j, is);

  fin += swin * (1 - alpha_a) * (1 - alpha_s) + tot_fluxdn -
         Constants::stefanBoltzmann * pow(btempArray(k, j), 4);

  // FIXME (cmetz) should loop over all precipitates, not just 0
  // e.g., int thermoIndex = 1 + nCloud + NVAPOR;
  // FIXME (cmetz) maybe later incorporate the q=mcat term (fix sign issue), for
  // now just include latent
  //

  int thermoIndex = 1 + 0 + NVAPOR;
  Real L = pthermo->GetLatentEnergyMass(thermoIndex, btempArray(k, j));

  // FIXME (cmetz) when snow is added, remove this hardcoded value
  Real Lfusion = 334000;
  Real f_evapCooling = (AmdEvap[0] * (L + Lfusion) + AmdEvap[1] * L) / dt;

  fin -= f_evapCooling;

  int sign = 0;
  // new precip will freeze, so LH of precip will be positive into surf
  if (btempArray(k, j) < meltingPointVapor1) sign = 1;
  // new precip will remain liquid, no new LH
  else if (btempArray(k, j) >= meltingPointVapor1)
    sign = 0;

  // energy it would take to freeze all the liquid
  Real f_latent_precip = sign * (accumPrecipAmd(0, 1) * Lfusion) / dt;

  // FIXME (cmetz) later, insert lines to check if just fin is enough to get the
  // temp to cross the melting point
  Real dTs = (fin + f_latent_precip) * (dt / cSurf);
  Real AmdToMelt = 0;
  Real deltaFin = 0;

  // heating up, crossing the melting point
  if (btempArray(k, j) + dTs > meltingPointVapor1 &&
      btempArray(k, j) < meltingPointVapor1) {
    Real dTGap = meltingPointVapor1 - btempArray(k, j);
    Real fgap = dTGap * cSurf / dt;
    fin -= fgap;
    dTs = 0;
    btempArray(k, j) = meltingPointVapor1;
  }
  // cooling down, crossing the melting point
  else if (btempArray(k, j) + dTs < meltingPointVapor1 &&
           btempArray(k, j) > meltingPointVapor1) {
    Real dTGap = btempArray(k, j) - meltingPointVapor1;
    Real fgap = dTGap * cSurf / dt;
    fin -= fgap;
    dTs = 0;
    btempArray(k, j) = meltingPointVapor1;
  }

  // at the freezing point, all fin goes towards latent
  if (btempArray(k, j) == meltingPointVapor1) {
    AmdToMelt =
        fin * dt /
        Lfusion;  // will be negative if fin is negative, implying freezing
    amd(0, 1, k, j) += accumPrecipAmd(0, 1);
    amd(0, 0, k, j) += accumPrecipAmd(0, 0);

    if (amd(0, 1, k, j) + AmdToMelt > 0 && amd(0, 0, k, j) - AmdToMelt > 0) {
      amd(0, 1, k, j) += AmdToMelt;
      amd(0, 0, k, j) -= AmdToMelt;
      //std::cout << "AmdToMelt: " << AmdToMelt << std::endl;
      fin = 0;
      deltaFin = 0;
    }
    // during freezing, more liquid would be frozen than exists
    else if (amd(0, 1, k, j) + AmdToMelt < 0) {
      amd(0, 0, k, j) += amd(0, 1, k, j);
      deltaFin = amd(0, 1, k, j) * Lfusion / dt;
      amd(0, 1, k, j) = 0;
    }
    // during melting, more solid would melt than exists
    else if (amd(0, 0, k, j) - AmdToMelt < 0) {
      amd(0, 1, k, j) += amd(0, 0, k, j);
      deltaFin = -amd(0, 0, k, j) * Lfusion / dt;
      amd(0, 0, k, j) = 0;
    }
    fin += deltaFin;
    dTs = fin * (dt / cSurf);
  }
  // the normal case when we are not near the freezing point
  else if (btempArray(k, j) > meltingPointVapor1) {
    // new precip will remain liquid
    amd(0, 1, k, j) += accumPrecipAmd(0, 1);
    amd(0, 0, k, j) += accumPrecipAmd(0, 0);
  } else if (btempArray(k, j) < meltingPointVapor1) {
    // new precip will freeze
    amd(0, 0, k, j) += accumPrecipAmd(0, 0) + accumPrecipAmd(0, 1);
  }

  btempArray(k, j) += dTs;

  gel(0, 1, k, j) = amd(0, 1, k, j) / rho_l_vapor1;
  gel(0, 0, k, j) = amd(0, 0, k, j) / rho_s_vapor1;

  return dTs;
}

AthenaArray<Real> Surface::AccumulatePrecipitates(MeshBlock *pmb, int k, int j,
                                                  int iSkim) {
  Real precip = 0;
  AthenaArray<Real> accumPrecipAmd;
  accumPrecipAmd.NewAthenaArray(NVAPOR, 2);
  for (int v = 0; v < NVAPOR; ++v) {
    for (int n = 0; n < 2; ++n) accumPrecipAmd(v, n) = 0;
  }

  for (int v = 0; v < NVAPOR; ++v) {
    for (int i = is; i <= iSkim; ++i) {
      precip = 0;
      for (int n = 1; n < 2; ++n) {
        precip = pmb->pscalars->r(n, k, j, i);
        pmb->pscalars->r(n, k, j, i) = 0;
        accumPrecipAmd(v, 1) = precip * pmb->phydro->w(IDN, k, j, is) * dzPBL;
      }
    }
  }
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
  // set air to have T in IDN slot
  air.w[IDN] = (air.w[IPR] * Mbar) / (air.w[IDN] * Constants::Rgas);

  std::vector<Real> AmdEvap(2, 0.);

  for (int n = 1; n <= NVAPOR; ++n) {
    std::vector<Real> rates(1 + pthermo->GetCloudIndexSet(n).size(), 0.);
    //(cmetz) check this value of CDE, see Hartmann pg 117
    rates = this->CalcSurfEvapRates(air, n, amd(0, 1, k, j), btempArray(k, j),
                                    cSurf, dt, 3e-3, Mbar);

    // evaporate off of liquid amd
    if (amd(0, 1, k, j) - rates[0] > 0) {
      amd(0, 1, k, j) += -rates[0];
      AmdEvap[1] = rates[0];
    } else {
      AmdEvap[1] = amd(0, 1, k, j);
      rates[0] = amd(0, 1, k, j);
      amd(0, 1, k, j) = 0;
    }

    // add or subtract from vapor pool
    pmb->phydro->w(n, k, j, is) += rates[0] * (Mbar / pthermo->GetMu(n)) /
                                   (pmb->phydro->w(IDN, k, j, is) * dzPBL);
    if (pmb->phydro->w(n, k, j, is) < 0) {
      pmb->phydro->w(n, k, j, is) = 0;
      std::cout << "liquid evaporation induced negative water in vapor pool"
                << std::endl;
    }

    // solid evaporation
    //  FIXME (cmetz) later, when we have snow, add this back in
    // rates = this->CalcSurfEvapRates(air, n, amd(0, 0, k, j), btempArray(k,
    // j),
    //                                 cSurf, dt, 3e-3, Mbar);

    // if (amd(0, 0, k, j) - rates[0] > 0) {
    //   amd(0, 0, k, j) += -rates[0];
    //   AmdEvap[0] = rates[0];
    // }
    // else {
    //   AmdEvap[0] = amd(0, 0, k, j);
    //   rates[0] = amd(0, 0, k, j);
    //   amd(0, 0, k, j) = 0;
    // }

    // add or subtract from vapor pool
    // pmb->phydro->w(n, k, j, is) += rates[0] * (Mbar / pthermo->GetMu(n)) /
    //                               (pmb->phydro->w(IDN, k, j, is) * dzPBL);
    // if (pmb->phydro->w(n, k, j, is) < 0) {
    //  pmb->phydro->w(n, k, j, is) = 0;
    //  std::cout << "solid evaporation induced negative water in vapor pool" <<
    //  std::endl;
    //}
  }

  gel(0, 1, k, j) = amd(0, 1, k, j) / rho_l_vapor1;
  gel(0, 0, k, j) = amd(0, 0, k, j) / rho_s_vapor1;

  return AmdEvap;
}
