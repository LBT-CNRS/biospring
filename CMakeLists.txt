
cmake_minimum_required(VERSION 3.12)

# Project name
project(BioSpring)

# Program version.
set(BIOSPRING_VERSION_MAJOR 1)
set(BIOSPRING_VERSION_MINOR 1)
set(BIOSPRING_VERSION_PATCH 1)

# ############################################################################
#
# General settings.
#
# ############################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Where to find additionnal cmake modules.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Enable/disable MDDriver support.
option(MDDRIVER_SUPPORT "MDDriver support for interactive simulations" OFF)

# Enable/disable OPENMP support.
option(OPENMP_SUPPORT "OpenMP support for parallelism" ON)

# Enable/disable OpenGL support.
option(OPENGL_SUPPORT "OpenGL support" OFF)

# Enable/disable FreeSASA support.
option(FREESASA_SUPPORT "FreeSASA support" OFF)

# Unit tests.
option(BUILD_TESTS "build unit tests" OFF)

# By default, we'll be installing BioSpring in /some/standard/path/BioSpring.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/BioSpring
        CACHE
        PATH
        "Install path prefix, prepended onto install directories."
        FORCE
    )
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# ############################################################################
#
# Search for external packages.
#
# ############################################################################

# Search for NetCDF.
find_package(NetCDF_CXX REQUIRED)

# Search for MDDriver if required.
if(MDDRIVER_SUPPORT)
    find_package(MDDriver REQUIRED)
endif(MDDRIVER_SUPPORT)

# Search for OpenMP if required.
if(OPENMP_SUPPORT)
    find_package(OpenMP REQUIRED)
endif(OPENMP_SUPPORT)

# Search for OpenGL if required.
if(OPENGL_SUPPORT)
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)
    find_package(GLEW REQUIRED)
endif(OPENGL_SUPPORT)

# Search for Freesasa if required.
if(FREESASA_SUPPORT)
    find_package(FreeSASA REQUIRED)
endif(FREESASA_SUPPORT)

# ############################################################################
#
# Build the projet.
#
# ############################################################################

# ############################################################################
#
# If this is a git clone, uses the commit hash / tag instead of the revision
# number.
#
# ############################################################################
set(BIOSPRING_VERSION "${BIOSPRING_VERSION_MAJOR}.${BIOSPRING_VERSION_MINOR}.${BIOSPRING_VERSION_PATCH}")

if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
    # Searches git, dies if cannot find it.
    find_package(Git REQUIRED)

    # Does the current commit as a tag ?
    execute_process(COMMAND git tag --points-at HEAD OUTPUT_VARIABLE CURRENT_GIT_VERSION
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    if(CURRENT_GIT_VERSION)
    # Else uses the short git hash
    else()
        execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE CURRENT_GIT_VERSION
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
    endif()

    # Remove end-of-line carriage return.
    string(REGEX REPLACE "\n$" "" CURRENT_GIT_VERSION "${CURRENT_GIT_VERSION}")

    # Replaces classical version string with git version.
    set(BIOSPRING_VERSION "${BIOSPRING_VERSION}-${CURRENT_GIT_VERSION}")
endif()

# Creates version.h which contains the version string.
configure_file(${PROJECT_SOURCE_DIR}/src/version.h.in ${PROJECT_SOURCE_DIR}/src/version.h)

add_subdirectory(src/IO/xdrfile)

add_library(biospring-core STATIC
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
    src/IO/XTCTrajWriter.cpp
    src/IO/modern/CSVTrajectoryWriter.cpp
    src/IO/modern/PDBTrajectoryWriter.cpp
    src/IO/modern/XTCTrajectoryWriter.cpp
    src/IO/modern/TrajectoryManager.cpp

    src/Constraint.cpp
    src/InsertionVector.cpp
    src/Vector3f.cpp

    src/interactor/Interactor.cpp
    src/interactor/mddriver/InteractorMDDriver.cpp
    src/interactor/mddriver/CustomData.cpp
    src/interactor/freesasa/InteractorFreeSASA.cpp

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

target_include_directories(biospring-core PUBLIC
    ${PROJECT_SOURCE_DIR}/src/
    ${PROJECT_SOURCE_DIR}/src/spn/
)

target_compile_options(biospring-core PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
)

target_link_libraries(biospring-core
    PUBLIC
    netCDF::netCDF_CXX
    xdrfile
)

if(OpenMP_CXX_FOUND)
    target_compile_definitions(biospring-core PUBLIC -DOPENMP_SUPPORT)
    target_link_libraries(biospring-core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MDDRIVER_FOUND)
    target_compile_definitions(biospring-core PUBLIC -DMDDRIVER_SUPPORT)
    target_link_libraries(biospring-core PUBLIC MDDriver::MDDriver)
endif()

if(FREESASA_SUPPORT)
    target_compile_definitions(biospring-core PUBLIC -DFREESASA_SUPPORT)
    target_link_libraries(biospring-core PUBLIC FreeSASA::FreeSASA)
endif()

# if(OPENGL_FOUND)
# target_compile_definitions(biospring-core PUBLIC -DOPENGL_SUPPORT)
# target_link_libraries(biospring-core PUBLIC OpenGL::GL GLUT::GLUT GLEW::glew)
# endif()

# Executable compilation.
set(TOOLS biospring editspn mergespn pdb2spn)

foreach(TARGET ${TOOLS})
    add_executable(${TARGET} src/cli/${TARGET}-cli.cpp src/cli/${TARGET}.cpp)
    target_compile_options(${TARGET} PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W4>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
    )
    target_link_libraries(${TARGET} PRIVATE biospring-core)
endforeach()

# Unit tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/tests)
endif()



# ############################################################################
#
# Installation rules...
#
# ############################################################################

# Set the installation location of biospring components.
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
install(DIRECTORY ${BIOSPRING_INCLUDE_DIRS_SOURCE_TREE}
    DESTINATION ${INCLUDE_INSTALL_DIR}
    FILES_MATCHING PATTERN "*.h"
)

# install targets
set(TARGETS ${TOOLS} biospring-core)
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
    USE_SOURCE_PERMISSIONS
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
endif(NOT DEFINED CPACK_SYSTEM_NAME)

set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION_MAJOR "${BIOSPRING_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${BIOSPRING_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${BIOSPRING_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "biospring ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "biospring-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.${SVN_WC_REVISION}")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
    /\\\\.svn/;
    /\\\\.swp$/;
    /.~*/;
    /build.*/;
)

include(CPack)
