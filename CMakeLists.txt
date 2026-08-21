cmake_minimum_required(VERSION 3.21)

project(BioSpring
    VERSION 1.1.1
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Without this, a plain `cmake -S . -B build` leaves CMAKE_BUILD_TYPE empty and
# single-config generators then pass no optimisation flags at all: the whole
# force loop, which is header-only and relies on inlining, builds at -O0. The
# guard on CMAKE_CONFIGURATION_TYPES leaves multi-config generators (MSVC,
# Xcode) alone, where the build type is chosen at build time instead.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(MDDRIVER_SUPPORT "Enable MDDriver coupling" OFF)
option(FREESASA_SUPPORT "Enable FreeSASA coupling" OFF)
option(OPENCL_SUPPORT "Enable the OpenCL backend" OFF)
option(OPENGL_SUPPORT "Enable the OpenGL viewer" OFF)
option(OPENMP_SUPPORT "Enable OpenMP parallelism" OFF)

# Keep the original, simple install-prefix behavior from the initial CMakeLists.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/BioSpring
        CACHE
        PATH
        "Install path prefix, prepended onto install directories."
        FORCE
    )
endif()

option(BUILD_TESTING "Build BioSpring tests" OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(NetCDF REQUIRED)
find_package(NetCDFCXX REQUIRED)
find_package(Threads REQUIRED)

if(MDDRIVER_SUPPORT)
    find_package(MDDriver REQUIRED)
endif()

if(FREESASA_SUPPORT)
    find_package(FreeSASA REQUIRED)
endif()

if(OPENCL_SUPPORT)
    find_package(OpenCL REQUIRED)
    # The bundled legacy OpenCL C++ header exposes OpenGL interop types.
    find_package(OpenGL REQUIRED)
endif()

if(OPENGL_SUPPORT)
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)
    find_package(GLEW REQUIRED)
endif()

if(OPENMP_SUPPORT)
    # Apple's system Clang (AppleClang) ships without OpenMP support at all,
    # so find_package(OpenMP) fails out of the box on macOS even when a
    # working OpenMP runtime is installed. Homebrew's `libomp` is the common
    # way to get one, but it is keg-only (not linked into a compiler search
    # path), so CMake cannot discover it without help. Auto-detect it here,
    # unless the user already told CMake how to find OpenMP explicitly.
    if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND NOT OpenMP_CXX_FLAGS)
        find_program(_biospring_brew_executable brew NO_CACHE)
        if(_biospring_brew_executable)
            execute_process(
                COMMAND "${_biospring_brew_executable}" --prefix libomp
                OUTPUT_VARIABLE _biospring_libomp_prefix
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
        endif()

        if(_biospring_libomp_prefix AND EXISTS "${_biospring_libomp_prefix}/include/omp.h")
            message(STATUS "Using Homebrew libomp for OpenMP: ${_biospring_libomp_prefix}")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${_biospring_libomp_prefix}/include")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY "${_biospring_libomp_prefix}/lib/libomp.dylib")
        endif()
        unset(_biospring_libomp_prefix)
    endif()

    find_package(OpenMP REQUIRED COMPONENTS CXX)
else()
    # Keep CMake GUI clean when OpenMP support is disabled.
    # These cache entries are created by CMake's FindOpenMP module after a
    # previous configure with OPENMP_SUPPORT=ON. They are not BioSpring options.
    get_cmake_property(_biospring_cache_variables CACHE_VARIABLES)
    foreach(_biospring_cache_variable IN LISTS _biospring_cache_variables)
        if(_biospring_cache_variable MATCHES "^OpenMP(_|$)")
            unset(${_biospring_cache_variable} CACHE)
        endif()
    endforeach()
    unset(_biospring_cache_variables)
    unset(_biospring_cache_variable)
endif()

set(BIOSPRING_VERSION "${PROJECT_VERSION}")
find_package(Git QUIET)
if(Git_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --tags --always --dirty
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_DESCRIPTION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(GIT_DESCRIPTION)
        set(BIOSPRING_VERSION "${PROJECT_VERSION}-${GIT_DESCRIPTION}")
    endif()
endif()

configure_file(
    "${PROJECT_SOURCE_DIR}/src/version.h.in"
    "${PROJECT_BINARY_DIR}/generated/version.h"
    @ONLY
)

add_subdirectory(src/IO/xdrfile)

set(BIOSPRING_CORE_SOURCES
    src/cli/argparse.cpp
    src/configuration/SafeConfigurationReader.cpp
    src/forcefield/ForceField.cpp
    src/grid/GridCoordinatesSystem.cpp
    src/grid/PotentialGrid.cpp
    src/reduce/ParticleGroup.cpp
    src/reduce/Reduce.cpp
    src/reduce/Reducer.cpp
    src/reduce/ReduceRule.cpp
    src/IO/io.cpp
    src/IO/CSVSampleWriter.cpp
    src/IO/ForceFieldReader.cpp
    src/IO/NetCDFReader.cpp
    src/IO/NetCDFWriter.cpp
    src/IO/OpenDXReader.cpp
    src/IO/PDBReader.cpp
    src/IO/PDBTrajectoryWriter.cpp
    src/IO/PDBWriter.cpp
    src/IO/PQRReader.cpp
    src/IO/PQRWriter.cpp
    src/IO/ReduceRuleReader.cpp
    src/IO/RigidBodyRuleReader.cpp
    src/rigidbodygroup/RigidBodyRule.cpp
    src/rigidbodygroup/RigidBodyBuilder.cpp
    src/IO/XTCTrajWriter.cpp
    src/IO/modern/CSVTrajectoryWriter.cpp
    src/IO/modern/PDBTrajectoryWriter.cpp
    src/IO/modern/XTCTrajectoryWriter.cpp
    src/IO/modern/TrajectoryManager.cpp
    src/Constraint.cpp
    src/InsertionVector.cpp
    src/Vector3f.cpp
    src/interactor/Interactor.cpp
    src/logging.cpp
    src/measure.cpp
    src/spn/Particle.cpp
    src/spn/ParticleProperty.cpp
    src/spn/Spring.cpp
    src/spn/SpringNetwork.cpp
    src/topology/Particle.cpp
    src/topology/ParticleCollection.cpp
    src/topology/ParticleProperties.cpp
    src/rigidbody/RigidBody.cpp
    src/rigidbody/matrix.cpp
    src/rigidbody/Quaternion.cpp
    src/rigidbody/RigidBodiesManager.cpp
)

if(MDDRIVER_SUPPORT)
    list(APPEND BIOSPRING_CORE_SOURCES
        src/interactor/mddriver/InteractorMDDriver.cpp
        src/interactor/mddriver/CustomData.cpp
    )
endif()

if(FREESASA_SUPPORT)
    list(APPEND BIOSPRING_CORE_SOURCES
        src/interactor/freesasa/InteractorFreeSASA.cpp
    )
endif()

add_library(biospring-core STATIC ${BIOSPRING_CORE_SOURCES})
add_library(BioSpring::core ALIAS biospring-core)

target_include_directories(biospring-core PUBLIC
    "${PROJECT_BINARY_DIR}/generated"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src/spn"
)

target_link_libraries(biospring-core PUBLIC
    NetCDFCXX::NetCDFCXX
    BioSpring::xdrfile
    Threads::Threads
)

if(MSVC)
    target_compile_options(biospring-core PRIVATE /W4)
    target_compile_definitions(biospring-core PRIVATE
        _CRT_SECURE_NO_WARNINGS NOMINMAX _USE_MATH_DEFINES)
else()
    target_compile_options(biospring-core PRIVATE -Wall -Wextra)
endif()

if(OPENMP_SUPPORT)
    target_compile_definitions(biospring-core PUBLIC OPENMP_SUPPORT)
    target_link_libraries(biospring-core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MDDRIVER_SUPPORT)
    target_compile_definitions(biospring-core PUBLIC MDDRIVER_SUPPORT)
    if(TARGET MDDriver::MDDriver)
        target_link_libraries(biospring-core PUBLIC MDDriver::MDDriver)
    elseif(MDDRIVER_INCLUDE_DIR AND MDDRIVER_LIBRARY)
        target_include_directories(biospring-core PUBLIC "${MDDRIVER_INCLUDE_DIR}")
        target_link_libraries(biospring-core PUBLIC "${MDDRIVER_LIBRARY}")
    else()
        message(FATAL_ERROR "MDDriverConfig.cmake did not expose its include directory and library")
    endif()
endif()

if(FREESASA_SUPPORT)
    target_compile_definitions(biospring-core PUBLIC FREESASA_SUPPORT)
    target_link_libraries(biospring-core PUBLIC FreeSASA::FreeSASA)
endif()

if(OPENCL_SUPPORT)
    add_subdirectory(src/opencl)
endif()

if(OPENGL_SUPPORT)
    add_subdirectory(src/viewer)
endif()

function(add_biospring_cli_tool tool)
    add_executable(${tool} ${ARGN})
    target_link_libraries(${tool} PRIVATE BioSpring::core)
    if(MSVC)
        target_compile_options(${tool} PRIVATE /W4)
        target_compile_definitions(${tool} PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
    else()
        target_compile_options(${tool} PRIVATE -Wall -Wextra)
    endif()
endfunction()

foreach(tool IN ITEMS biospring editspn mergespn pdb2spn)
    add_biospring_cli_tool(${tool}
        "src/cli/${tool}-cli.cpp"
        "src/cli/${tool}.cpp"
    )
endforeach()

foreach(tool IN ITEMS pdb2cdl pdb2pqr pdb2pdbconect)
    add_biospring_cli_tool(${tool}
        "src/cli/pdb2spn-cli.cpp"
        "src/cli/pdb2spn-legacy.cpp"
    )
endforeach()

if(BUILD_TESTING)
    # Pulled in only when tests are actually requested: `include(CTest)` adds a
    # batch of CDash/Dart cache variables (DART_TESTING_TIMEOUT,
    # MEMORYCHECK_*, COVERAGE_*, SITE, BUILDNAME, ...) that are irrelevant
    # clutter in cmake-gui/ccmake when BUILD_TESTING is off.
    include(CTest)
    add_subdirectory(src/tests)

    # GoogleTest (fetched above, in src/tests/CMakeLists.txt) creates a batch
    # of cache variables that are internal implementation details: nothing a
    # BioSpring user configuring "do I want tests" needs to see or change.
    # They cannot be suppressed at the source since they come from CMake's
    # own modules or fetched, third-party CMake code; mark them advanced
    # instead so they stay out of the default cmake-gui/ccmake view, matching
    # the OpenMP/CTest cleanup above:
    #  - CMAKE_INSTALL_* (~16 vars): created by GoogleTest's unconditional
    #    include(GNUInstallDirs), regardless of INSTALL_GTEST. BioSpring's
    #    own install() rules do not use them.
    #  - FETCHCONTENT_*: generic FetchContent module housekeeping (download
    #    location, offline mode, per-dependency overrides); only relevant to
    #    someone mirroring or vendoring the googletest checkout by hand.
    #  - GTEST_HAS_ABSL: GoogleTest's own Abseil/RE2 integration switch;
    #    BioSpring does not use Abseil.
    #  - GTest_DIR: left over from the find_package(GTest CONFIG QUIET)
    #    probe above, always GTest_DIR-NOTFOUND since BioSpring fetches its
    #    own copy instead of using a system install.
    #  - INSTALL_GTEST, BUILD_GMOCK: forced to OFF by src/tests/CMakeLists.txt
    #    already; showing them as editable checkboxes would be misleading
    #    since any change would be overwritten on the next configure.
    mark_as_advanced(
        CMAKE_INSTALL_BINDIR
        CMAKE_INSTALL_DATADIR
        CMAKE_INSTALL_DATAROOTDIR
        CMAKE_INSTALL_DOCDIR
        CMAKE_INSTALL_INCLUDEDIR
        CMAKE_INSTALL_INFODIR
        CMAKE_INSTALL_LIBDIR
        CMAKE_INSTALL_LIBEXECDIR
        CMAKE_INSTALL_LOCALEDIR
        CMAKE_INSTALL_LOCALSTATEDIR
        CMAKE_INSTALL_MANDIR
        CMAKE_INSTALL_OLDINCLUDEDIR
        CMAKE_INSTALL_RUNSTATEDIR
        CMAKE_INSTALL_SBINDIR
        CMAKE_INSTALL_SHAREDSTATEDIR
        CMAKE_INSTALL_SYSCONFDIR
        FETCHCONTENT_BASE_DIR
        FETCHCONTENT_FULLY_DISCONNECTED
        FETCHCONTENT_QUIET
        FETCHCONTENT_SOURCE_DIR_GOOGLETEST
        FETCHCONTENT_UPDATES_DISCONNECTED
        FETCHCONTENT_UPDATES_DISCONNECTED_GOOGLETEST
        GTEST_HAS_ABSL
        GTest_DIR
        INSTALL_GTEST
        BUILD_GMOCK
    )
endif()


# ############################################################################
#
# Installation rules...
#
# ############################################################################

# Set the installation location of biospring components.  This deliberately
# keeps the simple layout from the original project instead of exporting a
# CMake package.
set(BIN_INSTALL_DIR bin)
set(LIB_INSTALL_DIR lib)
set(INCLUDE_INSTALL_DIR include)

set(BIOSPRING_INSTALL_SHARE_DIR share/biospring)

set(DATA_INSTALL_DIR ${BIOSPRING_INSTALL_SHARE_DIR}/data)
set(SCRIPTS_INSTALL_DIR ${BIOSPRING_INSTALL_SHARE_DIR}/scripts)
set(EXAMPLES_INSTALL_DIR ${BIOSPRING_INSTALL_SHARE_DIR}/examples)
set(DOC_INSTALL_DIR ${BIOSPRING_INSTALL_SHARE_DIR}/doc)

set(CMAKEMOD_INSTALL_DIR ${BIOSPRING_INSTALL_SHARE_DIR}/cmake)

# install headers
install(DIRECTORY "${PROJECT_SOURCE_DIR}/src/"
    DESTINATION ${INCLUDE_INSTALL_DIR}
    FILES_MATCHING PATTERN "*.h"
)

# install targets
set(TARGETS biospring editspn mergespn pdb2spn pdb2cdl pdb2pqr pdb2pdbconect biospring-core)
install(TARGETS ${TARGETS}
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
)

# install configure files
install(FILES
    ${PROJECT_BINARY_DIR}/BioSpring-config.sh
    DESTINATION ${BIN_INSTALL_DIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
)

# install data
install(DIRECTORY ${CMAKE_SOURCE_DIR}/data/
    DESTINATION ${DATA_INSTALL_DIR}
    USE_SOURCE_PERMISSIONS
)

# install examples
install(DIRECTORY ${CMAKE_SOURCE_DIR}/example/
    DESTINATION ${EXAMPLES_INSTALL_DIR}
    USE_SOURCE_PERMISSIONS
)

# install scripts
install(DIRECTORY ${CMAKE_SOURCE_DIR}/scripts/
    DESTINATION ${SCRIPTS_INSTALL_DIR}
    USE_SOURCE_PERMISSIONS
)

# install doc
install(DIRECTORY ${CMAKE_SOURCE_DIR}/doc/
    DESTINATION ${DOC_INSTALL_DIR}
)

# ############################################################################
#
# Configure files...
#
# ############################################################################
configure_file(${PROJECT_SOURCE_DIR}/config/BioSpring-config.sh.in
    ${PROJECT_BINARY_DIR}/BioSpring-config.sh @ONLY)

#
# CPack
#
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Biospring - A spring network simulation engine.")
set(CPACK_PACKAGE_VENDOR "Marc Baaden")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")

if(NOT DEFINED CPACK_SYSTEM_NAME)
    set(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif()

set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "biospring ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "biospring-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}.${SVN_WC_REVISION}")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
    /\\.svn/;
    /\\.swp$/;
    /.~*/;
    /build.*/;
)

include(CPack)
