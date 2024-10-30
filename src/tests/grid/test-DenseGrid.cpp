#include <gtest/gtest.h>

#include "grid/DenseGrid.hpp"

using namespace biospring::grid;

TEST(DenseGridTest, InitializeGrid)
{
    DenseGrid<int> grid;
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

    EXPECT_EQ(grid.size(), 10 * 10 * 10); // shape[0] * shape[1] * shape[2]
}

TEST(DenseGridTest, AddAndRetrieveData)
{
    DenseGrid<int> grid;
    grid.reshape({0.0, 0.0, 0.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0});

    discrete_coordinates cell{1, 2, 3};
    int value = 42;
    grid.at(cell) = value;

    EXPECT_TRUE(grid.has_cell(cell));
    EXPECT_EQ(grid.at(cell), value);

    real_coordinates position{1.5, 2.5, 3.5};
    EXPECT_TRUE(grid.has_cell(grid.cell_coordinates(position)));
    EXPECT_EQ(grid.at(position), value);
}

TEST(DenseGridTest, ClearGrid)
{
    DenseGrid<int> grid;
    grid.reshape({0.0, 0.0, 0.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0});
    ASSERT_EQ(grid.size(), 1000);

    discrete_coordinates cell{1, 2, 3};
    int value = 42;
    grid.at(cell) = value;
    ASSERT_EQ(grid.at(cell), value);

    grid.clear();
    EXPECT_EQ(grid.at(cell), 0);
    EXPECT_EQ(grid.size(), 1000); // grid is not resized when cleared
}

TEST(DenseGridTest, is_out_of_grid_discrete_coordinates)
{
    DenseGrid<int> grid;
    grid.reshape({0.0, 0.0, 0.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0});

    EXPECT_FALSE(grid.is_out_of_grid(discrete_coordinates({0, 0, 0})));
    EXPECT_FALSE(grid.is_out_of_grid(discrete_coordinates({9, 9, 9})));
    EXPECT_TRUE(grid.is_out_of_grid(discrete_coordinates({-1, 0, 0})));
    EXPECT_TRUE(grid.is_out_of_grid(discrete_coordinates({0, -1, 0})));
    EXPECT_TRUE(grid.is_out_of_grid(discrete_coordinates({0, 0, -1})));
    EXPECT_TRUE(grid.is_out_of_grid(discrete_coordinates({10, 0, 0})));
    EXPECT_TRUE(grid.is_out_of_grid(discrete_coordinates({0, 10, 0})));
    EXPECT_TRUE(grid.is_out_of_grid(discrete_coordinates({0, 0, 10})));
}

TEST(DenseGridTest, is_out_of_grid_real_coordinates)
{
    DenseGrid<int> grid;
    grid.reshape({0.0, 0.0, 0.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0});

    EXPECT_FALSE(grid.is_out_of_grid(real_coordinates({0.0, 0.0, 0.0})));
    EXPECT_FALSE(grid.is_out_of_grid(real_coordinates({9.0, 9.0, 9.0})));
    EXPECT_TRUE(grid.is_out_of_grid(real_coordinates({-1.0, 0.0, 0.0})));
    EXPECT_TRUE(grid.is_out_of_grid(real_coordinates({0.0, -1.0, 0.0})));
    EXPECT_TRUE(grid.is_out_of_grid(real_coordinates({0.0, 0.0, -1.0})));
    EXPECT_TRUE(grid.is_out_of_grid(real_coordinates({10.0, 0.0, 0.0})));
    EXPECT_TRUE(grid.is_out_of_grid(real_coordinates({0.0, 10.0, 0.0})));
    EXPECT_TRUE(grid.is_out_of_grid(real_coordinates({0.0, 0.0, 10.0})));
}

TEST(DenseGridTest, at_fails_when_out_of_grid)
{
    DenseGrid<int> grid;
    grid.reshape({0.0, 0.0, 0.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0});

    EXPECT_THROW(grid.at(discrete_coordinates({-1, 0, 0})), std::out_of_range);
    EXPECT_THROW(grid.at(discrete_coordinates({0, -1, 0})), std::out_of_range);
    EXPECT_THROW(grid.at(discrete_coordinates({0, 0, -1})), std::out_of_range);
    EXPECT_THROW(grid.at(discrete_coordinates({10, 0, 0})), std::out_of_range);
    EXPECT_THROW(grid.at(discrete_coordinates({0, 10, 0})), std::out_of_range);
    EXPECT_THROW(grid.at(discrete_coordinates({0, 0, 10})), std::out_of_range);
}

TEST(DenseGridOfContainersTest, AddToContainer)
{
    DenseGridOfContainers<int> grid;
    grid.reshape({0.0, 0.0, 0.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0});

    discrete_coordinates cell{1, 2, 3};
    int element = 42;
    grid.add(cell, element);

    EXPECT_EQ(grid.size(), 1000); // shape[0] * shape[1] * shape[2]
    EXPECT_EQ(grid.number_of_elements(), 1);

    EXPECT_EQ(grid.at(cell).size(), 1);
    EXPECT_EQ(grid.at(cell)[0], element);
}

// TODO: test cells_within_radius

// -- Main function  ----------------------------------------------------------
int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
