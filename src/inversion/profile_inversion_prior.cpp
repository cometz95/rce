/** @file profile_inversion_prior.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Nov 20, 2021 09:21:41 EST
 * @bug No known bugs.
 */

// C/C++ header
#include <iostream>

// Athena++ header
#include <mesh/mesh.hpp>
#include <hydro/hydro.hpp>

// harp2 headers
#include <configure.hpp>
#include "../debugger/debugger.hpp"
#include "gaussian_process.hpp"
#include "profile_inversion.hpp"

extern std::unique_ptr<Debugger> pdebug;

Real ProfileInversion::LogPriorProbability(Real **XpSample) const
{
  int nsample = plevel_.size();
  std::vector<Real> zlev(nsample);
  std::vector<Real> zstd(nsample);

  Real P0 = Constants::ReferencePressure;
  Real H0 = Constants::PressureScaleHeight;

  for (int i = 0; i < nsample; ++i)
    zlev[i] = -H0*log(plevel_[i]/P0);

  Real lnprior = 0.;
  for (auto m : idx_) {
    for (int i = 0; i < nsample; ++i)
      zstd[i] = Xstd_[m]*pow(exp(zlev[i]/H0), chi_);
    lnprior += gp_lnprior(SquaredExponential, XpSample[m],
      zlev.data(), zstd.data(), nsample, Xlen_[m]);
  }

  pdebug->Message("Log prior probability", lnprior);

  return lnprior;
}
