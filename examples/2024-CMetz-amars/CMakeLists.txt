
# 1. Compile amars
setup_problem(amars)

# 2. Copy input files to run directory
file(GLOB inputs *.inp *.yaml *.py)
foreach(input ${inputs})
  execute_process(COMMAND ln -sf ${input} ${CMAKE_BINARY_DIR}/bin/${inp})
endforeach()
