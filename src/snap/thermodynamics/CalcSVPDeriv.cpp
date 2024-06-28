// C/C++
#include <algorithm>

// canoe
#include <air_parcel.hpp>

// snap
#include "thermodynamics.hpp"

// i is defined by the position of the vapor in the airparcel (for (int i = 1; i
// <= NVAPOR; ++i) )
// n is defined by   for (int n = 0; n < pthermo->GetCloudIndexSet(i).size();
// ++n) j is defined by int j = pthermo->GetCloudIndex(i, n);

Real Thermodynamics::CalcSVPDeriv(const Real temp, const int i, const int j,
                                  const int n) const {
  // airparcel which holds temp in the IDN slot
  AirParcel temp_container(AirParcel::Type::MoleFrac);
  temp_container.w[IDN] = temp;

  // getting des/dt locally around btemp with definition of the derivative
  Real T2 = temp + 0.01;
  temp_container.w[IDN] = T2;
  Real ln_es_2 = std::log(svp_func1_[i][n](temp_container, i, j));
  Real T1 = temp - 0.01;
  temp_container.w[IDN] = T1;
  Real ln_es_1 = std::log(svp_func1_[i][n](temp_container, i, j));
  Real dlnes_dt = (ln_es_2 - ln_es_1) / (T2 - T1);

  return dlnes_dt;
}
