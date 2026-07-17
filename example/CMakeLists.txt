# --------------------------
# CMakeLists.txt for HydroGrav/example
# --------------------------

# -------------------------
# Output directory for binaries
# -------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# --------------------------
# Link libraries
# --------------------------
link_libraries(
  HydroGravLib
)

# --------------------------
# Compiler warnings macro
# --------------------------
macro(add_example_executable name)
  add_executable(${name} ${ARGN})
  if(_hydrograv_warning_flags)
    target_compile_options(${name} PRIVATE ${_hydrograv_warning_flags})
  endif()
endmacro()

# --------------------------
# Example executables
# --------------------------
add_example_executable(run_fluid_profile run_fluid_profile.cpp)
add_example_executable(run_eos_gw_spectrum run_eos_gw_spectrum.cpp)
add_example_executable(run_kinetic_spectrum run_kinetic_spectrum.cpp)
add_example_executable(run_gw_spectrum run_gw_spectrum.cpp)
add_example_executable(run_bps run_bps.cpp)
add_example_executable(run_inv_pt run_inv_pt.cpp)