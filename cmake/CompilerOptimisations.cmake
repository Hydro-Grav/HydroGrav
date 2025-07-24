# CompilerOptimisations.cmake
# 
# This file defines compiler optimization flags for the DeepPhase project.
# It provides options for enabling various performance optimizations and
# sets up the appropriate flags based on the compiler being used.

# --------------------------
# Optimization Options
# --------------------------
option(ENABLE_FAST_MATH "Enable fast math optimizations" OFF)
option(ENABLE_NATIVE_ARCH "Enable native architecture optimizations" OFF)
option(ENABLE_AGGRESSIVE_OPTS "Enable aggressive optimizations for Release builds" ON)

# --------------------------
# Initialize optimization flags
# --------------------------
set(OPTIMIZATION_FLAGS_TO_ADD "")

# --------------------------
# Fast Math Optimizations
# --------------------------
if(ENABLE_FAST_MATH)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        list(APPEND OPTIMIZATION_FLAGS_TO_ADD "-ffast-math")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        list(APPEND OPTIMIZATION_FLAGS_TO_ADD "/fp:fast")
    endif()
    message(STATUS "Fast math optimizations enabled")
endif()

# --------------------------
# Native Architecture Optimizations
# --------------------------
if(ENABLE_NATIVE_ARCH)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        list(APPEND OPTIMIZATION_FLAGS_TO_ADD "-march=native" "-mtune=native")
    endif()
    message(STATUS "Native architecture optimizations enabled")
endif()

# --------------------------
# Build Type Specific Optimizations
# --------------------------
if(ENABLE_AGGRESSIVE_OPTS AND CMAKE_BUILD_TYPE STREQUAL "Release")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        list(APPEND OPTIMIZATION_FLAGS_TO_ADD "-O3" "-funroll-loops")
        # Additional aggressive optimizations
        list(APPEND OPTIMIZATION_FLAGS_TO_ADD "-fomit-frame-pointer")
        list(APPEND OPTIMIZATION_FLAGS_TO_ADD "-finline-functions")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        list(APPEND OPTIMIZATION_FLAGS_TO_ADD "/O2" "/Ob2")
    endif()
    message(STATUS "Aggressive Release optimizations enabled")
endif()

# --------------------------
# Vectorization Optimizations
# --------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Enable auto-vectorization reports (helpful for debugging performance)
    if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            list(APPEND OPTIMIZATION_FLAGS_TO_ADD "-ftree-vectorize")
        endif()
    endif()
endif()

# --------------------------
# Summary
# --------------------------
if(OPTIMIZATION_FLAGS_TO_ADD)
    message(STATUS "Optimization flags: ${OPTIMIZATION_FLAGS_TO_ADD}")
else()
    message(STATUS "No additional optimization flags enabled")
endif()
