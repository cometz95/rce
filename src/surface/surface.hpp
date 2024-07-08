#ifndef SRC_SURFACE_SURFACE_HPP_
#define SRC_SURFACE_SURFACE_HPP_

// C/C++ headers
#include <memory>
#include <string>
#include <vector>

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <constants.hpp>
#include <virtual_groups.hpp>

// Forward declarations
class MeshBlock;
class ParameterInput;

using RealArrayX = std::vector<Real>;

class Surface : public ParameterGroup {
 public:
  Surface(MeshBlock *pmb, ParameterInput *pin);
  ~Surface();

  bool hasSurface = true;

  void DoSurfaceProcesses(MeshBlock *pmb);

  // calculates precipitates evaporation rates
  //  FIXME (cmetz) pass in the amd arrays, not just one amd(j) value
  RealArrayX CalcSurfEvapRates(AirParcel const &qfrac, int i, Real &amd_x,
                               Real btemp, Real cSurf, Real dt, Real Cde,
                               Real Mbar) const;

  Real ChangeTempFromForcing(MeshBlock *pmb, int k, int j, Real dt,
                             AthenaArray<Real> &accumPrecipAmd,
                             RealArrayX AmdEvap);
  AthenaArray<Real> AccumulatePrecipitates(MeshBlock *pmb, int k, int j,
                                           int iSkim);
  RealArrayX EvapPrecip(MeshBlock *pmb, int k, int j, Real dt);

  size_t RestartDataSizeInBytes() const { return 0; }
  size_t DumpRestartData(char *pdst) const { return 0; }
  size_t LoadRestartData(char *psrc) { return 0; }

  AthenaArray<Real> GetBTempArray() const { return btempArray; }
  AthenaArray<Real> GetAmd() const { return amd; }
  AthenaArray<Real> GetGel() const { return gel; }

 protected:
  AthenaArray<Real> btempArray;

  // dimensions (NVAPOR, 2, ncells3, ncells2)
  // the amd and gel containers hold NVAPOR arrays, each of which are 2
  // (numphases) long. 0 slot is sold, 1 slot is liquid.
  // and each of those hold an array ncells3/ncells2 long (number of
  // cells in the k/j direction on surface)
  AthenaArray<Real>
      amd;  // areal mass density of precipitates on surface #kg/m^2
  AthenaArray<Real> gel;  // global equivalent layer #m

  // init variables
  int is, js, ks;
  int ie, je, ke;
  Real dzPBL;
  Real meltingPointVapor1;
  Real rho_l_vapor1;
  Real rho_s_vapor1;
  Real cSurf;
  Real dt;
  Real time;
  Real alpha_s;
  Real alpha_a;
  Real omega;
  Real s0;

  // state variables
  // forcing (positive into surface), units of W/m^2
  Real fin;
};

using SurfacePtr = std::shared_ptr<Surface>;

#endif  // SRC_SURFACE_SURFACE_HPP_
