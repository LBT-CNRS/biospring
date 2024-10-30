include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../cli)

# List of modules to be tested.
set(TEST_CLI_MODULES
    pdb2spn
    biospring
)

set(TEST_MODULES
    electrostatic-energy
    steric-energy
    spring-energy
)

foreach(MODULE ${TEST_CLI_MODULES})
    add_test_module_cli(${MODULE})
endforeach()

foreach(MODULE ${TEST_MODULES})
    add_test_module(${MODULE})
endforeach()
