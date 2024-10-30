#include <gtest/gtest.h>

#include "grid/PotentialGrid.hpp"

#include "timeit.hpp"

using namespace biospring::grid;
struct TestPotentialGrid : public ::testing::Test
{
    PotentialGrid grid;

    void SetUp() override
    {
        ::testing::Test::SetUp();

        grid.reshape({0.0, 0.0, 0.0, 10, 10, 10}, {0.1, 0.1, 0.1});

        ASSERT_EQ(100, grid.shape()[0]);
        ASSERT_EQ(100, grid.shape()[1]);
        ASSERT_EQ(100, grid.shape()[2]);

        ASSERT_EQ(0.1, grid.cell_size()[0]);
        ASSERT_EQ(0.1, grid.cell_size()[1]);
        ASSERT_EQ(0.1, grid.cell_size()[2]);

        ASSERT_EQ(0.0, grid.origin()[0]);
        ASSERT_EQ(0.0, grid.origin()[1]);
        ASSERT_EQ(0.0, grid.origin()[2]);
    }
};

TEST_F(TestPotentialGrid, compute_gradient)
{
    grid.at(discrete_coordinates(0, 0, 0)).scalar = 0.0;
    grid.at(discrete_coordinates(1, 0, 0)).scalar = 5.0;
    grid.at(discrete_coordinates(2, 0, 0)).scalar = 10.0;

    grid.compute_gradient();

    EXPECT_FLOAT_EQ(-50.0 * grid.GRADIENT_SCALE, grid.at(discrete_coordinates(1, 0, 0)).vector.getX());
}

TEST_F(TestPotentialGrid, accessToElements)
{
    biospring::timeit::Timer t1;

    for (size_t i = 0; i < grid.shape()[0]; i++)
        for (size_t j = 0; j < grid.shape()[1]; j++)
            for (size_t k = 0; k < grid.shape()[2]; k++)
            {
                discrete_coordinates cell(i, j, k);

                float scalar = static_cast<float>(i);
                Vector3f vector(static_cast<double>(i), static_cast<double>(i + 1), static_cast<double>(i + 2));

                grid.at(cell).scalar = scalar;
                grid.at(cell).vector = vector;
            }

    for (size_t i = 0; i < grid.shape()[0]; i++)
        for (size_t j = 0; j < grid.shape()[1]; j++)
            for (size_t k = 0; k < grid.shape()[2]; k++)
            {
                discrete_coordinates cell(i, j, k);

                EXPECT_FLOAT_EQ(static_cast<float>(i), grid.at(cell).scalar);

                const Vector3f expected(static_cast<float>(i), static_cast<float>(i + 1), static_cast<float>(i + 2));
                const Vector3f & actual = grid.at(cell).vector;
                EXPECT_FLOAT_EQ(expected.getX(), actual.getX());
                EXPECT_FLOAT_EQ(expected.getY(), actual.getY());
                EXPECT_FLOAT_EQ(expected.getZ(), actual.getZ());
            }
}

TEST(PotentialGrid, compute_gradient_dies_when_grid_not_initialized)
{
    PotentialGrid grid;
    ASSERT_EQ(grid.size(), 0);
    EXPECT_DEATH(grid.compute_gradient(), "cannot compute gradient of grid which size is 0");
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
