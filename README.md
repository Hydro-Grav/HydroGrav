<h1 align="center">
HydroGrav
</h1>

<div align="center">
<i>Precise hydrodynamics for gravitational waves from phase transitions.</i>
</div>

**HydroGrav** is a C++17 software package for calculating fluid profiles and gravitational wave spectra from a first-order electroweak phase transition using either a simplified equation of state (bag or $\mu\nu$ models) or a via a direct calculation of thermodynamics from the effective potential.

## Dependencies

You need a C++17 compliant compiler, git, and the following dependencies:
* CMake, version 3.11 or higher
* OpenMP, version 4.5 or higher
* GSL, version 2.7.1 or higher
* ALGLIB, version 4.0 or higher
* BOOST, version 1.83 or higher

On *Ubuntu/Debian*-based distributions, `ALGLIB` and `GSL` can be installed by running:

    sudo apt install libalglib-dev libgsl-dev libboost-all-dev

On *Fedora*-based distributions, instead use:

    sudo dnf install alglib-devel gsl-devel boost-devel

Finally on *Mac*:

    brew install gsl alglib boost


**Note:**
- OpenMP is usually included with your C++ compiler (e.g., `g++` or `clang++`). If you encounter errors related to OpenMP, ensure your compiler supports it and is up to date.


## Building
To build the shared library and examples, use: ***URL subject to change***

    git clone https://github.com/William-Searle/HydroGrav
    cd HydroGrav
    mkdir build
    cd build
    cmake ..
    make


### CMake Options

You can customize the build using the following CMake options:

| Option                    | Default | Description                                      |
|---------------------------|---------|--------------------------------------------------|
| `BUILD_WITH_UNIT_TESTS`   | ON      | Build and enable unit tests                      |
| `BUILD_WITH_EXAMPLES`     | ON      | Build example executables                        |
| `ENABLE_COMPILER_WARNINGS`| ON      | Enable additional compiler warnings              |

To set an option, pass it to CMake with `-D`, for example:

    cmake -D BUILD_WITH_UNIT_TESTS=OFF ..


## MatPlotLib-CPP
<details>
<summary>Click me</summary>

`HydroGrav` includes optional plotting functionality by utilizing the `matplotlib-cpp` library. To enable these features:

1. Download the [`matplotlibcpp.h`](https://github.com/lava/matplotlib-cpp/blob/master/matplotlibcpp.h) header file from the [matplotlib-cpp repository](https://github.com/lava/matplotlib-cpp).
2. Place it in the `HydroGrav/include/` directory, so the file structure is:

   ```
   HydroGrav/include/matplotlibcpp.h
   ```

3. Run `cmake ..` in your `build` directory. During the build process, you should see the message:

   ```
   -- Matplotlib header found. Plotting functionality will be enabled.
   ```

This indicates that all optional plotting methods have been successfully compiled.


**Note:**
The current version of `matplotlib-cpp` requires `numpy` version **less than 2.0.0** (`v<2.x.x`) and `matplotlib` to function.
These can be installed on Ubuntu/Debian with (not recommended):

    sudo apt install python3 python3-numpy python3-matplotlib

Or by using (recommended):

    python3 -m venv venv
    source venv/bin/activate
    pip install 'numpy<2.0.0' matplotlib

`HydroGrav` has been tested with `numpy` `1.26.4`.

</details>

## Running
If the library and examples were successfully built, the examples and tests are available to run using the following executables:

    ./bin/run_fluid_profile
    ./bin/run_kinetic_spectrum
    ./bin/run_gw_spectrum
    ./bin/run_eos_gw_spectrum
    ./bin/unit_tests
