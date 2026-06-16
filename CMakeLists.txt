cmake_minimum_required(VERSION 3.11)

project(HydroGrav VERSION 1.0.0 LANGUAGES CXX)

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

# For gpu offloading
# set(CMAKE_C_COMPILER clang-18)
# set(CMAKE_CXX_COMPILER clang++-18)

# --------------------------
# Compiler warnings
# --------------------------
option(ENABLE_COMPILER_WARNINGS "Enable compiler warnings" ON)
if(ENABLE_COMPILER_WARNINGS)
  include(CompilerWarnings)
  set(_hydrograv_warning_flags "${WARNING_FLAGS_TO_ADD}")
  # Propagate warning flags to subdirectories
  set(_hydrograv_warning_flags "${_hydrograv_warning_flags}" CACHE INTERNAL "HydroGrav compiler warning flags")
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
find_package(ALGLIB REQUIRED)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED)

# --------------------------
# Source and header files
# --------------------------
file(GLOB SOURCES "src/*.cpp")
set(HEADERS 
  include/constants.hpp
  include/hydrograv.hpp
  include/maths.hpp
  include/phasetransition.hpp
  include/physics.hpp
  include/profile.hpp
  include/ssm.hpp
  include/sici.hpp
  include/sici_arrays.hpp
  include/solvers.tpp
  include/snr.hpp
  include/detectors.hpp
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

add_library(HydroGravLib STATIC ${SOURCES} ${HEADERS})
add_library(HydroGrav::HydroGravLib ALIAS HydroGravLib)

if(ENABLE_COMPILER_WARNINGS AND _hydrograv_warning_flags)
  target_compile_options(HydroGravLib PRIVATE ${_hydrograv_warning_flags})
endif()

target_include_directories(HydroGravLib
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/HydroGrav>
  ${ALGLIB_INCLUDES}
  ${GSL_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
)

target_link_libraries(HydroGravLib
  PUBLIC
  OpenMP::OpenMP_CXX
  ${ALGLIB_LIBRARIES}
  GSL::gsl GSL::gslcblas
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

  target_include_directories(HydroGravLib
    PUBLIC
      ${Python3_INCLUDE_DIRS}
      ${NumPy_INCLUDE_DIR}
  )

  target_link_libraries(HydroGravLib
    PUBLIC
      Python3::Python
  )

  target_compile_definitions(HydroGravLib PUBLIC ENABLE_MATPLOTLIB)
endif()

# --------------------------
# Optional profiler library
# --------------------------
find_library(PROFILER_LIB profiler)
if(PROFILER_LIB)
  target_link_libraries(HydroGravLib PUBLIC ${PROFILER_LIB})
else()
  message(WARNING "Profiler library not found; skipping.")
endif()

# --------------------------
# Installation rules
# --------------------------
install(TARGETS HydroGravLib
  EXPORT HydroGravTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/HydroGrav)

install(EXPORT HydroGravTargets
  FILE HydroGravTargets.cmake
  NAMESPACE HydroGrav::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/HydroGrav
)

# --------------------------
# Config and version files
# --------------------------
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/HydroGravConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/HydroGravConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/HydroGravConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/HydroGrav
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/HydroGravConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/HydroGravConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/HydroGrav
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