#include <gtest/gtest.h>

#include "box.hpp"

TEST(Box, initialization_from_box)
{
    std::array<double, 6> box = {3.0, 3.0, 3.0, 10.0, 10.0, 10.0};
    biospring::Box b(box);

    EXPECT_FLOAT_EQ(3.0, b.boundaries()[0]);
    EXPECT_FLOAT_EQ(3.0, b.boundaries()[1]);
    EXPECT_FLOAT_EQ(3.0, b.boundaries()[2]);
    EXPECT_FLOAT_EQ(10.0, b.boundaries()[3]);
    EXPECT_FLOAT_EQ(10.0, b.boundaries()[4]);
    EXPECT_FLOAT_EQ(10.0, b.boundaries()[5]);

    EXPECT_EQ(7.0, b.length_x());
    EXPECT_EQ(7.0, b.length_y());
    EXPECT_EQ(7.0, b.length_z());
}

TEST(Box, initialization_from_origin_and_length)
{
    std::array<double, 3> origin = {3.0, 3.0, 3.0};
    std::array<double, 3> length = {7.0, 7.0, 7.0};
    biospring::Box b(origin, length);

    EXPECT_FLOAT_EQ(3.0, b.boundaries()[0]);
    EXPECT_FLOAT_EQ(3.0, b.boundaries()[1]);
    EXPECT_FLOAT_EQ(3.0, b.boundaries()[2]);
    EXPECT_FLOAT_EQ(10.0, b.boundaries()[3]);
    EXPECT_FLOAT_EQ(10.0, b.boundaries()[4]);
    EXPECT_FLOAT_EQ(10.0, b.boundaries()[5]);

    EXPECT_EQ(7.0, b.length_x());
    EXPECT_EQ(7.0, b.length_y());
    EXPECT_EQ(7.0, b.length_z());
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
