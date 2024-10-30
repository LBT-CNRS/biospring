include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.11.0
)
FetchContent_MakeAvailable(googletest)

# Setup Google Test framework.
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
    measure
    nsearch
    utils

    Box
    Configuration
    ForceFieldReader
    OpenDXReader
    Reducer
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