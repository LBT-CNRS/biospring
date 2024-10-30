#include <gtest/gtest.h>

#include "grid/InfiniteGrid.hpp"

using namespace biospring::grid;

TEST(InfiniteGridTest, InitializeGrid)
{
    InfiniteGrid<int> grid;
    grid.initialize({1.0, 1.0, 1.0});

    EXPECT_FLOAT_EQ(grid.cell_size()[0], 1.0);
    EXPECT_FLOAT_EQ(grid.cell_size()[1], 1.0);
    EXPECT_FLOAT_EQ(grid.cell_size()[2], 1.0);

    EXPECT_EQ(grid.size(), 0);
}

TEST(InfiniteGridTest, AddAndRetrieveData)
{
    InfiniteGrid<int> grid;
    grid.initialize({1.0, 1.0, 1.0});

    discrete_coordinates cell{1, 2, 3};
    int value = 42;
    grid.add(cell, value);

    EXPECT_TRUE(grid.has_cell(cell));
    EXPECT_EQ(grid.at(cell), value);

    real_coordinates position{1.5, 2.5, 3.5};
    EXPECT_TRUE(grid.has_cell(grid.cell_coordinates(position)));
    EXPECT_EQ(grid.at(position), value);
}

TEST(InfiniteGridTest, ClearGrid)
{
    InfiniteGrid<int> grid;
    grid.initialize({1.0, 1.0, 1.0});

    discrete_coordinates cell{1, 2, 3};
    int value = 42;
    grid.add(cell, value);

    EXPECT_EQ(grid.size(), 1);
    grid.clear();
    EXPECT_EQ(grid.size(), 0);
    EXPECT_FALSE(grid.has_cell(cell));
}

TEST(InfiniteGridOfContainersTest, AddToContainer)
{
    InfiniteGridOfContainers<int> grid;
    grid.initialize({1.0, 1.0, 1.0});

    discrete_coordinates cell{1, 2, 3};
    int element = 42;
    grid.add(cell, element);

    EXPECT_EQ(grid.size(), 1);
    EXPECT_EQ(grid.at(cell).size(), 1);
    EXPECT_EQ(grid.at(cell)[0], element);
}

// TODO: test cells_within_radius

int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
