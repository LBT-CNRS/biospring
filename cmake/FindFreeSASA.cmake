# This cmake module aims to find the FreeSASA depencies in the context
# of the BioSpring project.
#
# Sets:
#     - FREESASA_FOUND
#     - FREESASA_LIBRARY
#     - FREESASA_INCLUDE_DIR

SET(FREESASA_PREFIX ""
    CACHE
    PATH
    "Path to search for FreeSASA header and library files")

FIND_PATH(FREESASA_INCLUDE_DIR
    NAMES     freesasa.h
    HINTS ENV C_INCLUDE_PATH CPLUS_INCLUDE_PATH
          ${FREESASA_PREFIX}/include
          /usr/local/include
          /usr/include
          /sw/include
          /opt/local/include
    )

FIND_LIBRARY(FREESASA_LIBRARY
    NAMES freesasa
    HINTS ENV LD_LIBRARY_PATH
          ENV DYLD_LIBRARY_PATH
          ${FREESASA_PREFIX}
          ${FREESASA_PREFIX}/lib64
          ${FREESASA_PREFIX}/lib
          /usr/local/lib64
          /usr/lib64
          /usr/lib64/netcdf-3
          /usr/local/lib
          /usr/lib
          /usr/lib/netcdf-3
          /sw/lib
          /opt/local/lib
)

MARK_AS_ADVANCED(
  FREESASA_LIBRARY
  FREESASA_INCLUDE_DIR
)

# handle the QUIETLY and REQUIRED arguments and set FREESASA_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FreeSASA
    DEFAULT_MSG
    FREESASA_LIBRARY
    FREESASA_INCLUDE_DIR
)

# Export target FreeSASA::FreeSASA
if (FreeSASA_FOUND AND NOT TARGET FreeSASA::FreeSASA)
  add_library(FreeSASA::FreeSASA STATIC IMPORTED)
  set_property(TARGET FreeSASA::FreeSASA PROPERTY IMPORTED_LOCATION ${FREESASA_LIBRARY})
  target_include_directories(FreeSASA::FreeSASA INTERFACE ${FREESASA_INCLUDE_DIR})
endif()
