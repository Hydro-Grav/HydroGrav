<h1 align="center">
HydroGrav
</h1>

<div align="center">
<i>Precision gravitational waves from phase transitions.</i>
</div>

**HydroGrav** is a C++17 software package for calculating gravitational wave spectra from a first-order phase transition using either the Bag model or a via a direct calculation of thermodynamics from the effective potential.

## Dependencies

You need a C++17 compliant compiler, git, and the following dependencies:
* CMake, version 3.11 or higher
* OpenMP, version 3.0 or higher
* GSL, version 2.0 or higher
* ALGLIB, version 3.17 or higher

On *Ubuntu/Debian*-based distributions, `ALGLIB` and `GSL` can be installed by running:

    sudo apt install libalglib-dev libgsl-dev

On *Fedora*-based distributions, instead use:

    sudo dnf install alglib-devel gsl-devel

Finally on *Mac*:

    brew install gsl alglib


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

</details>

## Running
If the library and examples were successfully built, the examples and tests are available to run using the following executables: ***subject to change***

    ./bin/run_fluid_profile
    ./bin/run_kinetic_spectrum
    ./bin/run_gw_spectrum
    ./bin/unit_tests

## Bugs
* Bugs still present in GW calculation (not ready for use)
* Fluid profile & kinetic power spectrum working well!
