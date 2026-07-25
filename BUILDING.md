# Building BioSpring

BioSpring requires CMake 3.21 or newer, a C++20 compiler, NetCDF C and
NetCDF C++4. OpenMP and the external integrations are optional and disabled unless explicitly enabled.

```sh
cmake -S . -B build
cmake --build build
```

## NetCDF discovery

The two APIs are discovered separately:

```cmake
find_package(NetCDF REQUIRED)
find_package(NetCDFCXX REQUIRED)
```

The find modules use package configuration files when they are found in normal
CMake search paths, then fall back to regular CMake header/library discovery.
Most NetCDF C++4 installations are found through the ``netcdf`` header and the
``netcdf-cxx4`` library; no C++ package-config directory is expected.

If automatic discovery fails, provide the headers and libraries directly:

```sh
cmake -S . -B build \
  -DNetCDF_INCLUDE_DIR=/path/to/netcdf-c/include \
  -DNetCDF_LIBRARY=/path/to/libnetcdf \
  -DNetCDFCXX_INCLUDE_DIR=/path/to/netcdf-cxx4/include \
  -DNetCDFCXX_LIBRARY=/path/to/libnetcdf-cxx4
```

The imported targets are `NetCDF::NetCDF` and
`NetCDFCXX::NetCDFCXX`. The C++ target links NetCDF C transitively.

## Optional dependencies

OpenMP is disabled by default so CMake does not populate `OpenMP_*` cache
entries in CMake GUI unless requested. Enable it explicitly with:

```sh
cmake -S . -B build -DOPENMP_SUPPORT=ON
```

When enabled, CMake requires `OpenMP::OpenMP_CXX`; when disabled, BioSpring
builds sequentially and does not call `find_package(OpenMP)`.

Other integrations are enabled explicitly:

```sh
cmake -S . -B build \
  -DMDDRIVER_SUPPORT=ON \
  -DFREESASA_SUPPORT=ON \
  -DOPENCL_SUPPORT=ON \
  -DOPENGL_SUPPORT=ON
```

Only the packages required by enabled features are searched. For MDDriver,
`MDDriver_DIR` may point to the directory containing `MDDriverConfig.cmake`.

## Tests

Tests are disabled by default. Enable them explicitly with:

```sh
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```
