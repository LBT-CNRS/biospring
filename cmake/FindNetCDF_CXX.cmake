#
# Find NetCDF CXX include directories and libraries
#
# NetCDF_CXX_FOUND             - Do not attempt to use NetCDF if "no", "0", or undefined.
# NetCDF_CXX_INCLUDE_DIR       - where to find NetCDF C headers
# NetCDF_CXX_LIBRARY           - list of libraries to link against when using NetCDF

# Exports also the target:
#     - netCDF::netCDF_CXX


set(NetCDF_CXX_PREFIX ""
  CACHE
  PATH
  "Path to search for NetCDF CXX header and library files"
)


# Locates netcdf-cxx headers.
find_path(NetCDF_CXX_INCLUDE_DIR
  NAMES netcdf
  HINTS
    ENV CPATH
    ENV C_INCLUDE_PATH
    ENV CPLUS_INCLUDE_PATH
    /usr/local/include
    /usr/include
    /sw/include
    /opt/local/include
)

# Locates netcdf-cxx library.
find_library(NetCDF_CXX_LIBRARY
  NAMES netcdf_c++4 netcdf-cxx4
  HINTS
    ENV LIBRARY_PATH
    ENV LD_LIBRARY_PATH
    ENV DYLD_LIBRARY_PATH
    ${NetCDF_CXX_PREFIX}
    ${NetCDF_CXX_PREFIX}/lib64
    ${NetCDF_CXX_PREFIX}/lib
    /usr/local/lib64
    /usr/lib64
    /usr/local/lib
    /usr/lib
    /sw/lib
    /opt/local/lib
)


find_package(NetCDF_C REQUIRED)


# handle the QUIETLY and REQUIRED arguments and set NetCDF_CXX_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NetCDF_CXX
    DEFAULT_MSG
    NetCDF_CXX_LIBRARY
    NetCDF_CXX_INCLUDE_DIR
)

mark_as_advanced(
    NetCDF_CXX_PREFIX
    NetCDF_CXX_INCLUDE_DIR
    NetCDF_CXX_LIBRARY
)

if (NetCDF_CXX_FOUND AND NOT TARGET netCDF::netCDF_CXX)
  add_library(netCDF::netCDF_CXX SHARED IMPORTED)
  set_property(TARGET netCDF::netCDF_CXX PROPERTY IMPORTED_LOCATION ${NetCDF_CXX_LIBRARY})
  target_include_directories(netCDF::netCDF_CXX INTERFACE ${NetCDF_CXX_INCLUDE_DIR})
  target_include_directories(netCDF::netCDF_CXX INTERFACE ${NetCDF_C_INCLUDE_DIR})
endif()
