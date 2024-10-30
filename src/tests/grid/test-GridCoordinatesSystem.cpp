#include <gtest/gtest.h>

#include <array>
#include <memory>

#include "grid/GridCoordinatesSystem.hpp"

using namespace biospring::grid;

struct TestGridCoordinatesSystem : public ::testing::Test
{
    std::unique_ptr<GridCoordinatesSystem> grid;

    void SetUp() override
    {
        ::testing::Test::SetUp();

        std::array<double, 6> box{-5, -5, -5, 6, 6, 6};
        std::array<double, 3> spacing{1, 1, 1};
        grid = std::make_unique<GridCoordinatesSystem>(box, spacing);

        ASSERT_EQ(11, grid->shape()[0]);
        ASSERT_EQ(11, grid->shape()[1]);
        ASSERT_EQ(11, grid->shape()[2]);

        ASSERT_FLOAT_EQ(1, grid->cell_size()[0]);
        ASSERT_FLOAT_EQ(1, grid->cell_size()[1]);
        ASSERT_FLOAT_EQ(1, grid->cell_size()[2]);

        ASSERT_FLOAT_EQ(-5.0, grid->origin()[0]);
        ASSERT_FLOAT_EQ(-5.0, grid->origin()[1]);
        ASSERT_FLOAT_EQ(-5.0, grid->origin()[2]);
    }
};

template <typename T> static bool contains(const T & container, const typename T::value_type & element)
{
    return std::find(container.begin(), container.end(), element) != container.end();
}

TEST_F(TestGridCoordinatesSystem, cells_within_offset)
{
    discrete_coordinates cell{0, 0, 0};
    std::vector<discrete_coordinates> neighbors = grid->cells_within_offset(cell, {1, 1, 1});

    ASSERT_EQ(neighbors.size(), 8);

    EXPECT_TRUE(contains(neighbors, discrete_coordinates(0, 0, 0)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(0, 0, 1)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(0, 1, 0)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(0, 1, 1)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(1, 0, 0)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(1, 0, 1)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(1, 1, 0)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(1, 1, 1)));

    neighbors = grid->cells_within_offset(cell, {2, 2, 2});

    ASSERT_EQ(neighbors.size(), 27);
}

TEST_F(TestGridCoordinatesSystem, cells_within_radius)
{
    discrete_coordinates cell{0, 0, 0};
    std::vector<discrete_coordinates> neighbors = grid->cells_within_radius(cell, 1);

    ASSERT_EQ(neighbors.size(), 8);

    EXPECT_TRUE(contains(neighbors, discrete_coordinates(0, 0, 0)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(0, 0, 1)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(0, 1, 0)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(0, 1, 1)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(1, 0, 0)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(1, 0, 1)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(1, 1, 0)));
    EXPECT_TRUE(contains(neighbors, discrete_coordinates(1, 1, 1)));
}

TEST_F(TestGridCoordinatesSystem, at_real_coordinates)
{
    // Success: retrieving a cell within the grid boundaries.
    EXPECT_EQ(discrete_coordinates(0, 0, 0), grid->at(real_coordinates(-5, -5, -5)));

    // Failure: retrieving a cell outside the grid boundaries.
    ASSERT_TRUE(grid->is_out_of_grid(real_coordinates(-5.1, -5.1, -5.1)));
    EXPECT_THROW(grid->at(real_coordinates(-5.1, -5.1, -5.1)), std::out_of_range);
}

TEST_F(TestGridCoordinatesSystem, is_out_of_grid_discrete_coordinates)
{
    EXPECT_FALSE(grid->is_out_of_grid(discrete_coordinates(0, 0, 0)));
    EXPECT_FALSE(grid->is_out_of_grid(discrete_coordinates(5, 5, 5)));
    EXPECT_FALSE(grid->is_out_of_grid(discrete_coordinates(10, 10, 10)));
    EXPECT_TRUE(grid->is_out_of_grid(discrete_coordinates(-1, -1, -1)));
    EXPECT_TRUE(grid->is_out_of_grid(discrete_coordinates(11, 11, 11)));
}

TEST_F(TestGridCoordinatesSystem, is_out_of_grid_real_coordinates)
{
    EXPECT_FALSE(grid->is_out_of_grid(real_coordinates(-5, -5, -5)));
    EXPECT_FALSE(grid->is_out_of_grid(real_coordinates(5.9, 5.9, 5.9)));
    EXPECT_TRUE(grid->is_out_of_grid(real_coordinates(-5.1, -5.1, -5.1)));

    auto b = grid->boundaries();

    EXPECT_TRUE(grid->is_out_of_grid(real_coordinates(6.0, 6.0, 6.0)));
}

TEST(GridCoordinatesSystem, initialization_from_boundaries)
{
    std::array<double, 6> box{0, 0, 0, 10, 10, 10};
    std::array<double, 3> spacing{1, 1, 1};

    GridCoordinatesSystem grid(box, spacing);

    EXPECT_EQ(grid.shape()[0], 10);
    EXPECT_EQ(grid.shape()[1], 10);
    EXPECT_EQ(grid.shape()[2], 10);

    EXPECT_EQ(grid.cell_size()[0], 1);
    EXPECT_EQ(grid.cell_size()[1], 1);
    EXPECT_EQ(grid.cell_size()[2], 1);

    EXPECT_EQ(grid.boundaries().min_x(), 0);
    EXPECT_EQ(grid.boundaries().min_y(), 0);
    EXPECT_EQ(grid.boundaries().min_z(), 0);

    EXPECT_EQ(grid.boundaries().max_x(), 10);
    EXPECT_EQ(grid.boundaries().max_y(), 10);
    EXPECT_EQ(grid.boundaries().max_z(), 10);
}

TEST(GridCoordinatesSystem, initialization_from_origin_and_shape)
{
    std::array<double, 3> origin{0, 0, 0};
    std::array<size_t, 3> shape{10, 10, 10};
    std::array<double, 3> spacing{1, 1, 1};

    GridCoordinatesSystem grid(origin, shape, spacing);

    EXPECT_EQ(grid.shape()[0], 10);
    EXPECT_EQ(grid.shape()[1], 10);
    EXPECT_EQ(grid.shape()[2], 10);

    EXPECT_EQ(grid.cell_size()[0], 1);
    EXPECT_EQ(grid.cell_size()[1], 1);
    EXPECT_EQ(grid.cell_size()[2], 1);

    EXPECT_EQ(grid.boundaries().min_x(), 0);
    EXPECT_EQ(grid.boundaries().min_y(), 0);
    EXPECT_EQ(grid.boundaries().min_z(), 0);

    EXPECT_EQ(grid.boundaries().max_x(), 10);
    EXPECT_EQ(grid.boundaries().max_y(), 10);
    EXPECT_EQ(grid.boundaries().max_z(), 10);
}

TEST(GridCoordinatesSystem, max_size)
{
    GridCoordinatesSystem grid;
    EXPECT_EQ(grid.max_size(), 0); // should be 0 when not initialized.

    std::array<double, 6> box{0, 0, 0, 10, 10, 10};
    std::array<double, 3> spacing{1, 1, 1};
    grid.initialize(box, spacing);

    EXPECT_EQ(grid.max_size(), 1000); // 10 x 10 x 10
}

// ============================================================================
//
//   Tests for coordinates.
//
// ============================================================================

TEST(DiscreteCoordinates, initialization)
{
    discrete_coordinates a{1, 2, 3};
    discrete_coordinates b{std::array<int, 3>{1, 2, 3}};
    discrete_coordinates c{std::vector<int>{1, 2, 3}};
    discrete_coordinates d{std::initializer_list<int>{1, 2, 3}};

    EXPECT_EQ(a.x, 1);
    EXPECT_EQ(a.y, 2);
    EXPECT_EQ(a.z, 3);

    EXPECT_EQ(b.x, 1);
    EXPECT_EQ(b.y, 2);
    EXPECT_EQ(b.z, 3);

    EXPECT_EQ(c.x, 1);
    EXPECT_EQ(c.y, 2);
    EXPECT_EQ(c.z, 3);

    EXPECT_EQ(d.x, 1);
    EXPECT_EQ(d.y, 2);
    EXPECT_EQ(d.z, 3);
}

TEST(DiscreteCoordinates, equality)
{
    discrete_coordinates a{1, 2, 3};
    discrete_coordinates b{1, 2, 3};
    discrete_coordinates c{1, 2, 4};

    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
}

TEST(RealCoordinates, equality)
{
    real_coordinates a{1.f, 2.f, 3.f};
    real_coordinates b{1.f, 2.f, 3.f};
    real_coordinates c{1.f, 2.f, 3.f};
    real_coordinates d{1.f, 2.f, 3.00001f};

    EXPECT_EQ(a, b);
    EXPECT_EQ(a, c);
    EXPECT_NE(a, d);
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
