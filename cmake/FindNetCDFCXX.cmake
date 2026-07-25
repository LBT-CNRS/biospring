#[=======================================================================[.rst:
FindNetCDFCXX
--------------

Find the NetCDF C++4 library.

This module deliberately avoids CONFIG-package lookup.  Many NetCDF C++4
installations provide only the C++ header named ``netcdf`` and the
``netcdf-cxx4`` library.  Discovery is therefore based only on headers and
libraries, which are the values users can provide manually when automatic
search fails.

NetCDF C is a required dependency and must be found first by ``FindNetCDF``.

Imported target
^^^^^^^^^^^^^^^

``NetCDFCXX::NetCDFCXX``

Hints accepted from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``NetCDFCXX_INCLUDE_DIR``
  Directory containing the C++4 header named ``netcdf``.
``NetCDFCXX_LIBRARY``
  Full path to the NetCDF C++4 library.
#]=======================================================================]

include(FindPackageHandleStandardArgs)

# Remove stale cache entries produced by earlier versions of this module that
# tried CONFIG mode.  NetCDF C++4 usually has no CMake config file, so this
# variable should not be exposed to users.
unset(netCDFCxx_DIR CACHE)
unset(netCDFCxx_DIR)
unset(netCDFCxx_FOUND)

find_package(NetCDF REQUIRED)

set(_NetCDFCXX_PREFIXES
    /usr
    /usr/local
    /opt/local
    /opt/netcdf
    /opt/netcdf-cxx
    /opt/netcdf-cxx4
    /opt/NetCDF
)

# Use the NetCDF C installation as the strongest hint.  NetCDF C++4 is often
# installed in the same prefix or a sibling prefix.
if(NetCDF_INCLUDE_DIR)
    list(APPEND _NetCDFCXX_PREFIXES
        "${NetCDF_INCLUDE_DIR}"
        "${NetCDF_INCLUDE_DIR}/.."
        "${NetCDF_INCLUDE_DIR}/../netcdf-cxx"
        "${NetCDF_INCLUDE_DIR}/../netcdf-cxx4"
    )
endif()

if(NetCDF_LIBRARY)
    get_filename_component(_NetCDF_LIBRARY_DIR "${NetCDF_LIBRARY}" DIRECTORY)
    if(_NetCDF_LIBRARY_DIR)
        list(APPEND _NetCDFCXX_PREFIXES
            "${_NetCDF_LIBRARY_DIR}"
            "${_NetCDF_LIBRARY_DIR}/.."
            "${_NetCDF_LIBRARY_DIR}/../netcdf-cxx"
            "${_NetCDF_LIBRARY_DIR}/../netcdf-cxx4"
        )
    endif()
endif()

if(APPLE)
    list(APPEND _NetCDFCXX_PREFIXES
        /opt/homebrew
        /opt/homebrew/opt/netcdf
        /opt/homebrew/opt/netcdf-cxx
        /opt/homebrew/opt/netcdf-cxx4
        /usr/local
        /usr/local/opt/netcdf
        /usr/local/opt/netcdf-cxx
        /usr/local/opt/netcdf-cxx4
        /opt/local
    )
    file(GLOB _NetCDFCXX_APPLE_PREFIXES LIST_DIRECTORIES TRUE
        /opt/homebrew/Cellar/netcdf/*
        /opt/homebrew/Cellar/netcdf-cxx/*
        /opt/homebrew/Cellar/netcdf-cxx4/*
        /usr/local/Cellar/netcdf/*
        /usr/local/Cellar/netcdf-cxx/*
        /usr/local/Cellar/netcdf-cxx4/*
        /opt/local/libexec/netcdf*
        /opt/local/lib/netcdf*
    )
    list(APPEND _NetCDFCXX_PREFIXES ${_NetCDFCXX_APPLE_PREFIXES})
endif()

if(WIN32)
    foreach(_root "$ENV{ProgramW6432}" "$ENV{ProgramFiles}" "$ENV{ProgramFiles(x86)}")
        if(_root)
            list(APPEND _NetCDFCXX_PREFIXES
                "${_root}/netCDF"
                "${_root}/NetCDF"
                "${_root}/netCDF-CXX"
                "${_root}/netCDF-CXX4"
                "${_root}/NetCDF-CXX"
                "${_root}/NetCDF-CXX4"
            )
            file(GLOB _NetCDFCXX_WINDOWS_PREFIXES LIST_DIRECTORIES TRUE
                "${_root}/netCDF*"
                "${_root}/NetCDF*"
            )
            list(APPEND _NetCDFCXX_PREFIXES ${_NetCDFCXX_WINDOWS_PREFIXES})
        endif()
    endforeach()
    list(APPEND _NetCDFCXX_PREFIXES
        C:/netCDF
        C:/NetCDF
        C:/netCDF-CXX
        C:/netCDF-CXX4
        C:/NetCDF-CXX
        C:/NetCDF-CXX4
    )
endif()

list(REMOVE_DUPLICATES _NetCDFCXX_PREFIXES)

find_path(NetCDFCXX_INCLUDE_DIR
    NAMES netcdf
    HINTS ${_NetCDFCXX_PREFIXES}
    PATH_SUFFIXES
        include
        Include
        Library/include
        include/netcdf
        include/netcdf-cxx
        include/netcdf-cxx4
    DOC "Directory containing the NetCDF C++4 header named netcdf"
)

find_library(NetCDFCXX_LIBRARY
    NAMES
        netcdf_c++4
        netcdf-cxx4
        netcdf_c++4_static
        netcdf-cxx4_static
        libnetcdf_c++4
        libnetcdf-cxx4
        netcdf_c++4_d
        netcdf-cxx4_d
        netcdf_c++4d
        netcdf-cxx4d
        netcdf_c++4_debug
        netcdf-cxx4_debug
        libnetcdf_c++4_d
        libnetcdf-cxx4_d
    HINTS ${_NetCDFCXX_PREFIXES}
    PATH_SUFFIXES
        lib
        lib64
        Library/lib
        bin
        Library/bin
    DOC "NetCDF C++4 library"
)

find_package_handle_standard_args(NetCDFCXX
    REQUIRED_VARS NetCDFCXX_LIBRARY NetCDFCXX_INCLUDE_DIR
    REASON_FAILURE_MESSAGE
        "Could not find NetCDF C++4. Provide NetCDFCXX_INCLUDE_DIR, the directory containing the C++ header <netcdf>, and NetCDFCXX_LIBRARY, the full path to the netcdf-cxx4 library."
)

if(NetCDFCXX_FOUND)
    set(NetCDFCXX_INCLUDE_DIRS "${NetCDFCXX_INCLUDE_DIR}")
    set(NetCDFCXX_LIBRARIES NetCDFCXX::NetCDFCXX)

    if(NOT TARGET NetCDFCXX::NetCDFCXX)
        add_library(NetCDFCXX::NetCDFCXX UNKNOWN IMPORTED)
        set_target_properties(NetCDFCXX::NetCDFCXX PROPERTIES
            IMPORTED_LOCATION "${NetCDFCXX_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${NetCDFCXX_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES NetCDF::NetCDF
        )
    endif()
endif()

mark_as_advanced(
    NetCDFCXX_INCLUDE_DIR
    NetCDFCXX_LIBRARY
)

unset(_NetCDFCXX_PREFIXES)
unset(_NetCDFCXX_APPLE_PREFIXES)
unset(_NetCDFCXX_WINDOWS_PREFIXES)
unset(_NetCDF_LIBRARY_DIR)
