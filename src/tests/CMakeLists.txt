# Setup Google Test framework. Prefer a GoogleTest installation already
# known to CMake. If none is available, fetch a version whose CMake project is
# compatible with recent CMake releases.
find_package(GTest CONFIG QUIET)
if(NOT GTest_FOUND)
    include(FetchContent)

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    # GoogleTest installs itself by default when embedded this way. BioSpring
    # only needs it to build the test binaries; it must not end up in
    # CMAKE_INSTALL_PREFIX alongside BioSpring's own files.
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    # BioSpring's tests only use gtest/gtest_main, never gmock: skip building
    # it entirely instead of just hiding the option, since it also saves real
    # compile time.
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
    )
    FetchContent_MakeAvailable(googletest)

    # GoogleTest is vendored, third-party code built as its own CMake
    # targets: it does not inherit BioSpring's own compile options, but the
    # compiler's baseline warnings (e.g. char8_t/char32_t conversion in
    # gtest-printers.h) still fire while building it. Silence warnings on its
    # targets rather than patching fetched source.
    foreach(_gtest_target IN ITEMS gtest gtest_main)
        if(TARGET ${_gtest_target})
            if(MSVC)
                target_compile_options(${_gtest_target} PRIVATE /w)
            else()
                target_compile_options(${_gtest_target} PRIVATE -w)
            endif()
        endif()
    endforeach()
    unset(_gtest_target)
endif()

include(GoogleTest)

# Copy test data for unit tests.
set(TEST_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)
file(COPY ${TEST_DATA_DIR}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

function(add_test_module MODULE)
    add_executable(test-${MODULE} ${CMAKE_CURRENT_SOURCE_DIR}/test-${MODULE}.cpp)
    target_link_libraries(test-${MODULE} PRIVATE biospring-core GTest::gtest_main)
    gtest_discover_tests(test-${MODULE})
endfunction()

function(add_test_module_cli MODULE)
    add_executable(test-${MODULE} ${CMAKE_CURRENT_SOURCE_DIR}/test-${MODULE}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../cli/${MODULE}-cli.cpp)
    target_link_libraries(test-${MODULE} PRIVATE biospring-core GTest::gtest_main)
    gtest_discover_tests(test-${MODULE})
endfunction()

# List of modules to be tested.
set(TEST_MODULES
    logging
    measure
    nsearch
    utils

    Box
    Configuration
    ForceFieldReader
    NetCDFRoundTrip
    OpenDXReader
    PDBReader
    Reducer
    RigidBodiesManager
    ReduceRuleReader
    Vector3f
)

# Creates the test executables for each module.
foreach(MODULE ${TEST_MODULES})
    add_test_module(${MODULE})
endforeach()

add_subdirectory(grid)
add_subdirectory(regression)
add_subdirectory(topology)