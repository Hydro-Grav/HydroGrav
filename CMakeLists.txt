cmake_minimum_required(VERSION 3.11)

project(DeepPhase VERSION 1.0.0 LANGUAGES CXX)

# --------------------------
# Output directories
# --------------------------
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
endif()

if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
endif()

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
endif()

if(NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
endif()

# --------------------------
# Installation directories (use GNUInstallDirs for standard paths)
# --------------------------
include(GNUInstallDirs)

# --------------------------
# C++ standard and compile options
# --------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --------------------------
# Compiler warnings
# --------------------------
option(ENABLE_COMPILER_WARNINGS "Enable compiler warnings" ON)
if(ENABLE_COMPILER_WARNINGS)
  include(CompilerWarnings)
  set(_deepphase_warning_flags "${WARNING_FLAGS_TO_ADD}")
  # Propagate warning flags to subdirectories
  set(_deepphase_warning_flags "${_deepphase_warning_flags}" CACHE INTERNAL "DeepPhase compiler warning flags")
endif()

# --------------------------
# Other compile options
# --------------------------
option(BUILD_WITH_UNIT_TESTS "Enable unit tests" ON)
option(BUILD_WITH_EXAMPLES "Enable examples" ON)


# --------------------------
# Dependencies
# --------------------------
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(ALGLIB REQUIRED)

# --------------------------
# Source and header files
# --------------------------
file(GLOB SOURCES "src/*.cpp")
set(HEADERS 
  include/constants.hpp
  include/CubicSpline.tpp
  include/deepphase.hpp
  include/hydrodynamics.hpp
  include/maths_ops.hpp
  include/phasetransition.hpp
  include/physics.hpp
  include/profile.hpp
  include/spectrum.hpp
  include/tests.hpp
)

# --------------------------
# Check for matplotlibcpp header
# --------------------------
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/matplotlibcpp.h")
  set(ENABLE_MATPLOTLIB ON)
  list(APPEND HEADERS include/matplotlibcpp.h)
  message(STATUS "Matplotlibcpp header found. Plotting functionality will be enabled.")
else()
  set(ENABLE_MATPLOTLIB OFF)
  message(STATUS "Matplotlibcpp header not found. Plotting functionality will be disabled.")
endif()

# --------------------------
# Library target
# --------------------------

add_library(DeepPhaseLib STATIC ${SOURCES} ${HEADERS})
add_library(DeepPhase::DeepPhaseLib ALIAS DeepPhaseLib)

if(ENABLE_COMPILER_WARNINGS AND _deepphase_warning_flags)
  target_compile_options(DeepPhaseLib PRIVATE ${_deepphase_warning_flags})
endif()

target_include_directories(DeepPhaseLib
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/DeepPhase>
  ${GSL_INCLUDE_DIRS}
  ${ALGLIB_INCLUDES}
)

target_link_libraries(DeepPhaseLib
  PUBLIC
  OpenMP::OpenMP_CXX
  GSL::gsl
  GSL::gslcblas
  ${ALGLIB_LIBRARIES}
)

# --------------------------
# Optional Matplotlib (Python/Numpy) integration
# --------------------------
if(ENABLE_MATPLOTLIB)
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message(STATUS "Found NumPy include directory: ${NumPy_INCLUDE_DIR}")

  target_include_directories(DeepPhaseLib
    PUBLIC
      ${Python3_INCLUDE_DIRS}
      ${NumPy_INCLUDE_DIR}
  )

  target_link_libraries(DeepPhaseLib
    PUBLIC
      Python3::Python
  )

  target_compile_definitions(DeepPhaseLib PUBLIC ENABLE_MATPLOTLIB)
endif()

# --------------------------
# Optional profiler library
# --------------------------
find_library(PROFILER_LIB profiler)
if(PROFILER_LIB)
  target_link_libraries(DeepPhaseLib PUBLIC ${PROFILER_LIB})
else()
  message(WARNING "Profiler library not found; skipping.")
endif()

# --------------------------
# Installation rules
# --------------------------
install(TARGETS DeepPhaseLib
  EXPORT DeepPhaseTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/DeepPhase)

install(EXPORT DeepPhaseTargets
  FILE DeepPhaseTargets.cmake
  NAMESPACE DeepPhase::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DeepPhase
)

# --------------------------
# Config and version files
# --------------------------
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/DeepPhaseConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeepPhaseConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/DeepPhaseConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DeepPhase
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/DeepPhaseConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/DeepPhaseConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DeepPhase
)

# --------------------------
# Other subdirectories
# --------------------------
if(BUILD_WITH_EXAMPLES)
  add_subdirectory(example)
endif()
if(BUILD_WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif()