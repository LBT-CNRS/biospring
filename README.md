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


## Compilation & Installation

### Requirements

- CMake (http://www.cmake.org)
- NetCDF library (http://www.unidata.ucar.edu/software/netcdf/)
- NetCDF library headers for the C++ library (https://github.com/Unidata/netcdf-cxx4)

### Optional tools

- FreeSASA library (https://github.com/mittinatten/freesasa/) version 2.1.2
- MDDriver library (https://github.com/LBT-CNRS/MDDriver) version 1.0


- UnityMol (http://unitymol.sourceforge.net/) or VMD visualization software (http://www.ks.uiuc.edu/Research/vmd/)


### Conda

A `conda-env.yml` is provided to easily set up a conda environnement with BioSpring (with MDDriver and FreeSASA support) installed:

    conda env create -f conda-env.yml -n biospring
    conda activate biospring

This is the recommended way for MacOS and Windows OS.


### Docker

A docker file `Dockerfile` is also provided to set up a docker container with BioSpring (with MDDriver and FreeSASA support).

Also, docker images are built and provided through the Container registry of GitHub. They are available trough the [Packages link](https://github.com/LBT-CNRS/biospring/pkgs/container/biospring).

Docker images are built only when tags or Release are created.

Pull the Docker Image:

    docker pull ghcr.io/lbt-cnrs/biospring:latest

Please read [User Manual](doc/User_Manual.md) for detailed usage of the docker image.

### Compilation

#### Linux

Cmake & NetCDF library installation is straightforward on Linux:

    apt install cmake libnetcdf-c++4-dev

Then, for BioSpring:

```
git clone https://github.com/LBT-CNRS/biospring
cd biospring/

# On Mac/Linux
cmake -S . -B build -DCMAKE_INSTALL_PREFIX:PATH=/path/to/biospring/installation

# Compilation
cmake --build build -j4
cmake --install
```

### BioSpring with MDDriver and FreeSASA

BioSpring can also be compiled with MDDriver and/or FreeSASA support with those cmake options:

    cmake -S . -B build -DMDDRIVER_SUPPORT=ON -DMDDriver_DIR=/path/to/mddriver/installation/MDDriver/cmake/ -DFREESASA_SUPPORT=ON -DFREESASA_PREFIX=/path/to/freesasa/installation/freesasa-2.1.2/


### Test

BioSpring integrates unit & regressions tests. You can test the software with those commands:

    cmake -S . -B build -DBUILD_TESTS=ON
    cmake --build build -j4
    ctest --test-dir build


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
