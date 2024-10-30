set(TEST_MODULES
    DenseGrid
    GridCoordinatesSystem
    InfiniteGrid
    PotentialGrid
    SparseGrid
)

foreach(MODULE ${TEST_MODULES})
    add_test_module(${MODULE})
endforeach()
