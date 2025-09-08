# Compiler Flags in DeepPhase

This document explains the main compiler flags and options used in the DeepPhase project, as set in the top-level `CMakeLists.txt` and related CMake modules. Generated using ChatGPT.

## C++ Standard
- **CMAKE_CXX_STANDARD 17**
  - The project is built using C++17 for modern language features and compatibility.

## Compiler Warnings
- **ENABLE_COMPILER_WARNINGS** (default: ON)
  - When enabled, includes extra warning flags to help catch potential issues during compilation.
  - Flags are set via the `cmake/CompilerWarnings.cmake` module and stored in the variable `${_deepphase_warning_flags}`.
  - These are applied to all targets. Typical flags may include:
    - `-Wall` : Enable most warning messages.
    - `-Wextra` : Enable additional warning messages.
    - `-Wpedantic` : Enforce strict ISO compliance.
    - (See `cmake/CompilerWarnings.cmake` for the full list.)

## Compiler Optimizations
- **ENABLE_COMPILER_OPTIMIZATIONS** (default: ON)
  - When enabled, adds optimization flags for performance.
  - Flags are set via the `cmake/CompilerOptimisations.cmake` module and stored in `${_deepphase_optimization_flags}`.
  - Typical flags may include:
    - `-O2` or `-O3` : Enable optimization levels 2 or 3.
    - `-march=native` : Optimize for the local machine architecture.
    - (See `cmake/CompilerOptimisations.cmake` for the full list.)

## OpenMP
- **find_package(OpenMP REQUIRED)**
  - Enables OpenMP support for parallelization.
  - Adds the necessary flags for your compiler (e.g., `-fopenmp` for GCC/Clang).

## Other Notable Options
- **BUILD_WITH_UNIT_TESTS** (default: ON)
  - Enables building unit tests.
- **BUILD_WITH_EXAMPLES** (default: ON)
  - Enables building example programs.

## Python/Numpy Integration (Optional)
- If `include/matplotlibcpp.h` is found, Python and NumPy headers are included for plotting support.
- The macro `ENABLE_MATPLOTLIB` is defined for conditional compilation.

## How to Change Flags
- To disable warnings or optimizations, configure with:
  ```sh
  cmake -DENABLE_COMPILER_WARNINGS=OFF -DENABLE_COMPILER_OPTIMIZATIONS=OFF ..
  ```
- To use a different compiler (e.g., Clang 18):
  Uncomment and set:
  ```cmake
  set(CMAKE_C_COMPILER clang-18)
  set(CMAKE_CXX_COMPILER clang++-18)
  ```

## References
- See `cmake/CompilerWarnings.cmake` and `cmake/CompilerOptimisations.cmake` for the full list of flags and logic.
- See the top-level `CMakeLists.txt` for all options and their default values.
