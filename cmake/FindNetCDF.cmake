#[=======================================================================[.rst:
FindNetCDF
----------

Find the NetCDF C library.

Discovery is based on the C header and library that are present in nearly all
installations.  This module deliberately avoids CONFIG-package lookup because
a failed CONFIG lookup leaves a useless cache entry visible in CMake GUI.

No package manager executable is invoked.

Imported target
^^^^^^^^^^^^^^^

``NetCDF::NetCDF``

Hints accepted from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``NetCDF_INCLUDE_DIR``
  Directory containing ``netcdf.h``.
``NetCDF_LIBRARY``
  Full path to the NetCDF C library.
#]=======================================================================]

include(FindPackageHandleStandardArgs)

# Remove stale cache entries produced by older versions of this module that
# tried CONFIG mode.  NetCDF is intentionally resolved from its public header
# and library here, so this variable should not be exposed to users.
unset(netCDF_DIR CACHE)
unset(netCDF_DIR)
unset(netCDF_FOUND)

# Remove stale split Release/Debug cache entries from older versions.  The
# module now exposes only the two values a user may need to set manually:
# NetCDF_INCLUDE_DIR and NetCDF_LIBRARY.
unset(NetCDF_LIBRARY_RELEASE CACHE)
unset(NetCDF_LIBRARY_DEBUG CACHE)
unset(NetCDF_LIBRARY_RELEASE)
unset(NetCDF_LIBRARY_DEBUG)

set(_NetCDF_PREFIXES
    /usr
    /usr/local
    /opt/local
    /opt/netcdf
    /opt/NetCDF
)

if(APPLE)
    list(APPEND _NetCDF_PREFIXES
        /opt/homebrew
        /opt/homebrew/opt/netcdf
        /usr/local/opt/netcdf
    )
    file(GLOB _NetCDF_APPLE_PREFIXES LIST_DIRECTORIES TRUE
        /opt/homebrew/Cellar/netcdf/*
        /usr/local/Cellar/netcdf/*
        /opt/local/libexec/netcdf*
        /opt/local/lib/netcdf*
    )
    list(APPEND _NetCDF_PREFIXES ${_NetCDF_APPLE_PREFIXES})
endif()

if(WIN32)
    foreach(_root "$ENV{ProgramW6432}" "$ENV{ProgramFiles}" "$ENV{ProgramFiles\(x86\)}")
        if(_root)
            list(APPEND _NetCDF_PREFIXES "${_root}/netCDF" "${_root}/NetCDF")
            file(GLOB _NetCDF_WINDOWS_PREFIXES LIST_DIRECTORIES TRUE
                "${_root}/netCDF*"
                "${_root}/NetCDF*"
            )
            list(APPEND _NetCDF_PREFIXES ${_NetCDF_WINDOWS_PREFIXES})
        endif()
    endforeach()
    list(APPEND _NetCDF_PREFIXES C:/netCDF C:/NetCDF)
endif()

list(REMOVE_DUPLICATES _NetCDF_PREFIXES)

find_path(NetCDF_INCLUDE_DIR
    NAMES netcdf.h
    HINTS ${_NetCDF_PREFIXES}
    PATH_SUFFIXES include Include Library/include
    DOC "Directory containing netcdf.h"
)

find_library(NetCDF_LIBRARY
    NAMES netcdf netcdf_static libnetcdf netcdf_d netcdfd netcdf_debug netcdf_static_d libnetcdf_d
    HINTS ${_NetCDF_PREFIXES}
    PATH_SUFFIXES lib lib64 Library/lib
    DOC "NetCDF C library"
)

if(NetCDF_INCLUDE_DIR AND EXISTS "${NetCDF_INCLUDE_DIR}/netcdf_meta.h")
    file(STRINGS "${NetCDF_INCLUDE_DIR}/netcdf_meta.h" _NetCDF_version_lines
        REGEX "^#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH)[ \t]+")
    foreach(_component MAJOR MINOR PATCH)
        string(REGEX MATCH
            "NC_VERSION_${_component}[ \t]+\\(?([0-9]+)\\)?"
            _match "${_NetCDF_version_lines}")
        set(_NetCDF_version_${_component} "${CMAKE_MATCH_1}")
    endforeach()
    if(_NetCDF_version_MAJOR MATCHES "^[0-9]+$"
       AND _NetCDF_version_MINOR MATCHES "^[0-9]+$"
       AND _NetCDF_version_PATCH MATCHES "^[0-9]+$")
        set(NetCDF_VERSION
            "${_NetCDF_version_MAJOR}.${_NetCDF_version_MINOR}.${_NetCDF_version_PATCH}")
    endif()
endif()

find_package_handle_standard_args(NetCDF
    REQUIRED_VARS NetCDF_LIBRARY NetCDF_INCLUDE_DIR
    VERSION_VAR NetCDF_VERSION
    REASON_FAILURE_MESSAGE
        "Could not find NetCDF C. Provide NetCDF_INCLUDE_DIR, the directory containing netcdf.h, and NetCDF_LIBRARY, the full path to the NetCDF C library."
)

if(NetCDF_FOUND)
    set(NetCDF_INCLUDE_DIRS "${NetCDF_INCLUDE_DIR}")
    set(NetCDF_LIBRARIES NetCDF::NetCDF)

    if(NOT TARGET NetCDF::NetCDF)
        add_library(NetCDF::NetCDF UNKNOWN IMPORTED)
        set_target_properties(NetCDF::NetCDF PROPERTIES
            IMPORTED_LOCATION "${NetCDF_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(
    NetCDF_INCLUDE_DIR
    NetCDF_LIBRARY
)

unset(_NetCDF_PREFIXES)
unset(_NetCDF_APPLE_PREFIXES)
unset(_NetCDF_WINDOWS_PREFIXES)
