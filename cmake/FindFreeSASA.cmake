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

set(_FreeSASA_INCLUDE_HINTS
    ${FREESASA_PREFIX}/include
    /usr/local/include
    /usr/include
    /sw/include
    /opt/local/include
)
set(_FreeSASA_LIBRARY_HINTS
    ${FREESASA_PREFIX}
    ${FREESASA_PREFIX}/lib64
    ${FREESASA_PREFIX}/lib
    /usr/local/lib64
    /usr/lib64
    /usr/local/lib
    /usr/lib
    /sw/lib
    /opt/local/lib
)

# MacPorts installs under /opt/local (already listed above). Homebrew's
# prefix differs between Apple Silicon and Intel Macs and is keg-only, so it
# is not on the compiler's default search path; add both explicitly.
if(APPLE)
    list(APPEND _FreeSASA_INCLUDE_HINTS /opt/homebrew/include /usr/local/opt/freesasa/include)
    list(APPEND _FreeSASA_LIBRARY_HINTS /opt/homebrew/lib /usr/local/opt/freesasa/lib)
    file(GLOB _FreeSASA_APPLE_PREFIXES LIST_DIRECTORIES TRUE
        /opt/homebrew/Cellar/freesasa/*
        /usr/local/Cellar/freesasa/*
    )
    foreach(_prefix ${_FreeSASA_APPLE_PREFIXES})
        list(APPEND _FreeSASA_INCLUDE_HINTS "${_prefix}/include")
        list(APPEND _FreeSASA_LIBRARY_HINTS "${_prefix}/lib")
    endforeach()
    unset(_FreeSASA_APPLE_PREFIXES)
    unset(_prefix)
endif()

if(WIN32)
    foreach(_root "$ENV{ProgramW6432}" "$ENV{ProgramFiles}" "$ENV{ProgramFiles\(x86\)}")
        if(_root)
            list(APPEND _FreeSASA_INCLUDE_HINTS "${_root}/freesasa/include" "${_root}/FreeSASA/include")
            list(APPEND _FreeSASA_LIBRARY_HINTS "${_root}/freesasa/lib" "${_root}/FreeSASA/lib")
        endif()
    endforeach()
endif()

FIND_PATH(FREESASA_INCLUDE_DIR
    NAMES     freesasa.h
    HINTS ENV C_INCLUDE_PATH CPLUS_INCLUDE_PATH
          ${_FreeSASA_INCLUDE_HINTS}
    )

FIND_LIBRARY(FREESASA_LIBRARY
    NAMES freesasa
    HINTS ENV LD_LIBRARY_PATH
          ENV DYLD_LIBRARY_PATH
          ${_FreeSASA_LIBRARY_HINTS}
)

unset(_FreeSASA_INCLUDE_HINTS)
unset(_FreeSASA_LIBRARY_HINTS)

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
  add_library(FreeSASA::FreeSASA UNKNOWN IMPORTED)
  set_property(TARGET FreeSASA::FreeSASA PROPERTY IMPORTED_LOCATION ${FREESASA_LIBRARY})
  target_include_directories(FreeSASA::FreeSASA INTERFACE ${FREESASA_INCLUDE_DIR})

  # Package-manager builds of FreeSASA (e.g. Homebrew) are commonly built
  # with XML and JSON output support, which pulls libxml2 and json-c symbols
  # into libfreesasa's static archive without recording either as a
  # dependency anywhere CMake can see. Link them whenever available so the
  # final executables do not fail with unresolved xmlNewDoc/xmlFreeDoc/... or
  # json_object_.../... symbols; this is harmless when the FreeSASA build in
  # use does not actually need them.
  find_package(LibXml2 QUIET)
  if(LibXml2_FOUND)
      target_link_libraries(FreeSASA::FreeSASA INTERFACE LibXml2::LibXml2)
  else()
      message(STATUS "FreeSASA found, but LibXml2 was not: if linking fails with "
          "undefined xmlNewDoc/xmlFreeDoc/... symbols, install libxml2 "
          "development files (the FreeSASA build in use was compiled with "
          "XML output support).")
  endif()

  # json-c has no CMake module bundled with CMake itself; prefer pkg-config
  # (how Homebrew, MacPorts and most Linux distributions expose it), with a
  # manual header/library search as a fallback for setups without
  # pkg-config (e.g. some Windows toolchains).
  find_package(PkgConfig QUIET)
  set(_FreeSASA_JSONC_FOUND FALSE)
  if(PkgConfig_FOUND)
      pkg_check_modules(JSONC QUIET IMPORTED_TARGET json-c)
      if(JSONC_FOUND)
          target_link_libraries(FreeSASA::FreeSASA INTERFACE PkgConfig::JSONC)
          set(_FreeSASA_JSONC_FOUND TRUE)
      endif()
  endif()
  if(NOT _FreeSASA_JSONC_FOUND)
      find_path(JSONC_INCLUDE_DIR NAMES json.h PATH_SUFFIXES json-c)
      find_library(JSONC_LIBRARY NAMES json-c)
      mark_as_advanced(JSONC_INCLUDE_DIR JSONC_LIBRARY)
      if(JSONC_INCLUDE_DIR AND JSONC_LIBRARY)
          if(NOT TARGET JSONC::JSONC)
              add_library(JSONC::JSONC UNKNOWN IMPORTED)
              set_target_properties(JSONC::JSONC PROPERTIES
                  IMPORTED_LOCATION "${JSONC_LIBRARY}"
                  INTERFACE_INCLUDE_DIRECTORIES "${JSONC_INCLUDE_DIR}"
              )
          endif()
          target_link_libraries(FreeSASA::FreeSASA INTERFACE JSONC::JSONC)
          set(_FreeSASA_JSONC_FOUND TRUE)
      endif()
  endif()
  if(NOT _FreeSASA_JSONC_FOUND)
      message(STATUS "FreeSASA found, but json-c was not: if linking fails with "
          "undefined json_object_.../... symbols, install json-c development "
          "files (the FreeSASA build in use was compiled with JSON output "
          "support).")
  endif()
  unset(_FreeSASA_JSONC_FOUND)
endif()
