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

  // accumulates precipitates at one surface box (from only the 1st atm layer,
  // if iSkim = is)
  void AccumulatePrecipitates(int iSkim, int j);

  // calculates precipitates evaporation rates
  //  FIXME (cmetz) pass in the amd arrays, not just one amd(j) value
  RealArrayX CalcSurfEvapRates(AirParcel const &qfrac, int i, Real &amd_x,
                               Real btemp, Real dTs, Real cSurf, Real dt,
                               Real Cde, Real Mbar) const;

  Real ChangeTempFromForcing(MeshBlock *pmb, int j, Real dt);
  void AccumulatePrecipitates(MeshBlock *pmb, int iSkim, int j);
  void EvapPrecip(MeshBlock *pmb, int j, double dTs, Real dt);

  size_t RestartDataSizeInBytes() const { return 0; }
  size_t DumpRestartData(char *pdst) const { return 0; }
  size_t LoadRestartData(char *psrc) { return 0; }

  AthenaArray<Real> GetBTempArray() const { return btempArray; }
  std::vector<std::vector<AthenaArray<Real>>> GetAmd() const { return amd; }
  std::vector<std::vector<AthenaArray<Real>>> GetGel() const { return gel; }

 protected:
  AthenaArray<Real> btempArray;

  // the amd and gel containers hold NVAPOR arrays, each of which are 2
  // (NPHASE-SURF) long, and each of those hold an array ncells2 long (number of
  // cells in the j direction on surface) 0 slot is sold, 1 slot is liquid
  std::vector<std::vector<AthenaArray<Real>>> amd;
  std::vector<std::vector<AthenaArray<Real>>> gel;

  int is, js, ks;
  int ie, je, ke;
  bool H2OisLiquid;
  double dzPBL;
  double rholH2O;
  double rhosH2O;
  double cSurf;
  double dt;
  double time;
  double alpha_s;
  double alpha_a;
  double omega;
  double s0;
};

using SurfacePtr = std::shared_ptr<Surface>;

#endif  // SRC_SURFACE_SURFACE_HPP_
