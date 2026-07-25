# BioSpring 

![GitHub Release](https://img.shields.io/github/v/release/LBT-CNRS/biospring)


BioSpring is a molecular simulation software based on a spring network, especially designed for
interactive simulations. The main goal of BioSpring is to allow user to interactively explore the biomechanical properties biomolecule using spring network model at different scales (all atom, coarsed grain, Carbon-Alpha, ...).

## Description

### Input

Biospring is a simulation program which takes two files as input. The first one describes the structure and the biophysical properties of the system you want to study using the NetCDF binary format to describe this data (`.nc`), and the second one a text file which describes the parameters and settings for the spring network simulation (`.msp`).

See the [User manual](doc/User_Manual.md) in the doc directory for more information.

See the [MSP Options](doc/MSP_Options.md) file for the full list of parameters.


### Features

Here is a list of the most important features added in this branch:

- Add Martini3 data for IMPALA coarse grained simulations
- Accept Martini3 structure input with bonds
- Update FreeSASA configuration in command line
- Add interactors & manage indepedent threads for MDDriver and FreeSASA
- Add CustomData class to InteractorMDDriver
- Extend IMPALA equations to switch between single flat to curved double membrane
- Add Rigidbody dynamics
- Add roll angle computation to InsertionVector
- Update and reformat other scripts in the project ...


## Compilation

BioSpring uses CMake and requires CMake 3.21 or newer, a C++20 compiler,
NetCDF C and NetCDF C++4. OpenMP, OpenCL, the OpenGL viewer, MDDriver and
FreeSASA are optional and disabled unless explicitly enabled.

```sh
cmake -S . -B build
cmake --build build
```

NetCDF is searched through `FindNetCDF.cmake` and `FindNetCDFCXX.cmake`.
The modules first try the official config packages and then accept explicit
include and library paths. See [BUILDING.md](BUILDING.md) for the available
variables and optional features.

Tests are disabled by default. They can be enabled with `-DBUILD_TESTING=ON` and run with:

```sh
ctest --test-dir build --output-on-failure
```

### Usage

A detailed explanation is available in the [User Manual](doc/User_Manual.md).

#### Examples

A simple 2 particles system is available in `example/2particles`.

## References


## License

BioSpring is licenced under the [CeCILL-C license](LICENSE.txt).

## Contributors :

- Marc Baaden
- Nicolas Ferey
- Olivier Delalande
- Benoist Laurent
- André Lanrezac
- Hubert Santuz


## Building

See [BUILDING.md](BUILDING.md).
