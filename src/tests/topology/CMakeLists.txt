set(TEST_MODULES
    Particle
    ParticleCollection
    ParticleProperties
    Spring
    SpringCollection
    Topology
)

foreach(MODULE ${TEST_MODULES})
    add_test_module(${MODULE})
endforeach()
