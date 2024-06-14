#ifndef SRC_SURFACE_SURFACE_HPP_
#define SRC_SURFACE_SURFACE_HPP_

// C/C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include <athena/athena.hpp>

// canoe
#include <virtual_groups.hpp>

// Forward declarations
class MeshBlock;
class ParameterInput;

//! \brief manages all physics package data and functions
class Surface : public ParameterGroup {
 public:
  // functions
  Surface(MeshBlock *pmb, ParameterInput *pin){};
  ~Surface(){};

  size_t RestartDataSizeInBytes() const { return 0; }
  size_t DumpRestartData(char *pdst) const { return 0; }
  size_t LoadRestartData(char *psrc) { return 0; }
};

using SurfacePtr = std::shared_ptr<Surface>;

#endif  // SRC_SURFACE_SURFACE_HPP_
