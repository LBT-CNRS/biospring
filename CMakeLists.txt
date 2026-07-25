cmake_minimum_required(VERSION 3.21)

project(BioSpring
    VERSION 1.1.1
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
include(CTest)

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

foreach(tool IN ITEMS biospring editspn mergespn pdb2spn)
    add_executable(${tool}
        "src/cli/${tool}-cli.cpp"
        "src/cli/${tool}.cpp"
    )
    target_link_libraries(${tool} PRIVATE BioSpring::core)
    if(MSVC)
        target_compile_options(${tool} PRIVATE /W4)
        target_compile_definitions(${tool} PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
    else()
        target_compile_options(${tool} PRIVATE -Wall -Wextra)
    endif()
endforeach()

if(BUILD_TESTING)
    add_subdirectory(src/tests)
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
install(DIRECTORY ${BIOSPRING_INCLUDE_DIRS_SOURCE_TREE}
    DESTINATION ${INCLUDE_INSTALL_DIR}
    FILES_MATCHING PATTERN "*.h"
)

# install targets
set(TARGETS biospring editspn mergespn pdb2spn biospring-core)
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
