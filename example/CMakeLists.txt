# --------------------------
# CMakeLists.txt for DeepPhase/example
# --------------------------

# -------------------------
# Output directory for binaries
# -------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# --------------------------
# Link libraries
# --------------------------
link_libraries(
  DeepPhaseLib
)

# --------------------------
# Compiler warnings macro
# --------------------------
macro(add_example_executable name)
  add_executable(${name} ${ARGN})
  if(_deepphase_warning_flags)
    target_compile_options(${name} PRIVATE ${_deepphase_warning_flags})
  endif()
endmacro()

# --------------------------
# Example executables
# --------------------------
add_example_executable(run_fluid_profile run_fluid_profile.cpp)
add_example_executable(run_kinetic_spectrum run_kinetic_spectrum.cpp)
add_example_executable(run_gw_spectrum run_gw_spectrum.cpp)