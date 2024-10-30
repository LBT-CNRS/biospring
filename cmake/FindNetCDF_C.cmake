#
# Find NetCDF C include directories and libraries
#
# NetCDF_C_FOUND             - Do not attempt to use NetCDF if "no", "0", or undefined.
# NetCDF_C_INCLUDE_DIR       - where to find NetCDF C headers
# NetCDF_C_LIBRARY           - list of libraries to link against when using NetCDF

# Exports also the target:
#     - netCDF::netCDF_C


set(NetCDF_C_PREFIX ""
  CACHE
  PATH
  "Path to search for NetCDF C header and library files"
)


# Locates netcdf-c headers.
find_path(NetCDF_C_INCLUDE_DIR
  NAMES netcdf.h
  HINTS
    ENV CPATH
    ENV C_INCLUDE_PATH
    /usr/local/include
    /usr/include
    /sw/include
    /opt/local/include
)


# Locates netcdf-c library.
find_library(NetCDF_C_LIBRARY
  NAMES netcdf
  HINTS
    ENV LIBRARY_PATH
    ENV LD_LIBRARY_PATH
    ENV DYLD_LIBRARY_PATH
    ${NetCDF_C_PREFIX}
    ${NetCDF_C_PREFIX}/lib64
    ${NetCDF_C_PREFIX}/lib
    /usr/local/lib64
    /usr/lib64
    /usr/local/lib
    /usr/lib
    /sw/lib
    /opt/local/lib
)


# handle the QUIETLY and REQUIRED arguments and set NetCDF_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NetCDF_C
    DEFAULT_MSG
    NetCDF_C_LIBRARY
    NetCDF_C_INCLUDE_DIR
)

mark_as_advanced(
    NetCDF_C_PREFIX
    NetCDF_C_INCLUDE_DIR
    NetCDF_C_LIBRARY
)


if (NetCDF_C_FOUND AND NOT TARGET netCDF::netCDF_C)
  add_library(netCDF::netCDF_C SHARED IMPORTED)
  set_property(TARGET netCDF::netCDF_C PROPERTY IMPORTED_LOCATION ${NetCDF_C_LIBRARY})
  target_include_directories(netCDF::netCDF_C INTERFACE ${NetCDF_C_INCLUDE_DIR})
endif()
