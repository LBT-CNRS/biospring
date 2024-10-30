#include <gtest/gtest.h>

#include "grid/SparseGrid.hpp"

using namespace biospring::grid;

TEST(SparseGridTest, InitializeGrid)
{
    SparseGrid<int> grid;
    grid.reshape({0.0, 0.0, 0.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0});

    EXPECT_FLOAT_EQ(grid.cell_size()[0], 1.0);
    EXPECT_FLOAT_EQ(grid.cell_size()[1], 1.0);
    EXPECT_FLOAT_EQ(grid.cell_size()[2], 1.0);

    EXPECT_EQ(grid.shape()[0], 10);
    EXPECT_EQ(grid.shape()[1], 10);
    EXPECT_EQ(grid.shape()[2], 10);

    EXPECT_FLOAT_EQ(grid.boundaries().min_x(), 0.0);
    EXPECT_FLOAT_EQ(grid.boundaries().min_y(), 0.0);
    EXPECT_FLOAT_EQ(grid.boundaries().min_z(), 0.0);
    EXPECT_FLOAT_EQ(grid.boundaries().max_x(), 10.0);
    EXPECT_FLOAT_EQ(grid.boundaries().max_y(), 10.0);
    EXPECT_FLOAT_EQ(grid.boundaries().max_z(), 10.0);

    EXPECT_EQ(grid.size(), 0);
    EXPECT_EQ(grid.max_size(), 10 * 10 * 10); // shape[0] * shape[1] * shape[2]
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
