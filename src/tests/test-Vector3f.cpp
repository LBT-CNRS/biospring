#include <gtest/gtest.h>
#include "Vector3f.h"

// =====================================================================================
// Vector3f:: vector addition.
TEST(TestVector3f, vector_addition)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(4.0, 5.0, 6.0);

    Vector3f result = lhs + rhs;

    EXPECT_FLOAT_EQ(5.0, result.getX());
    EXPECT_FLOAT_EQ(7.0, result.getY());
    EXPECT_FLOAT_EQ(9.0, result.getZ());
}

// =====================================================================================
// Vector3f:: vector subtraction.
TEST(TestVector3f, vector_subtraction)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(4.0, 5.0, 6.0);

    Vector3f result = lhs - rhs;

    EXPECT_FLOAT_EQ(-3.0, result.getX());
    EXPECT_FLOAT_EQ(-3.0, result.getY());
    EXPECT_FLOAT_EQ(-3.0, result.getZ());
}

// =====================================================================================
// Vector3f:: vector scaling.
TEST(TestVector3f, vector_scaling)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    float scale = 2.0;

    Vector3f result = lhs * scale;

    EXPECT_FLOAT_EQ(2.0, result.getX());
    EXPECT_FLOAT_EQ(4.0, result.getY());
    EXPECT_FLOAT_EQ(6.0, result.getZ());
}

// =====================================================================================
// Vector3f:: vector cross product.
TEST(TestVector3f, vector_cross_product)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(4.0, 5.0, 6.0);

    Vector3f result = lhs ^ rhs;

    EXPECT_FLOAT_EQ(-3.0, result.getX());
    EXPECT_FLOAT_EQ(6.0, result.getY());
    EXPECT_FLOAT_EQ(-3.0, result.getZ());
}

// =====================================================================================
// Vector3f:: vector term-by-term product.
TEST(TestVector3f, vector_term_by_term_product)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(4.0, 5.0, 6.0);

    Vector3f result = lhs * rhs;

    EXPECT_FLOAT_EQ(4.0, result.getX());
    EXPECT_FLOAT_EQ(10.0, result.getY());
    EXPECT_FLOAT_EQ(18.0, result.getZ());
}

// =====================================================================================
// Vector3f:: vector norm.
TEST(TestVector3f, vector_norm)
{
    Vector3f lhs(1.0, 2.0, 3.0);

    float result = lhs.norm();

    EXPECT_FLOAT_EQ(std::sqrt(14.0), result);
}

// =====================================================================================
// Vector3f:: vector normalization.
TEST(TestVector3f, vector_normalization)
{
    Vector3f lhs(1.0, 2.0, 3.0);

    lhs.normalize();

    EXPECT_FLOAT_EQ(1.0 / std::sqrt(14.0), lhs.getX());
    EXPECT_FLOAT_EQ(2.0 / std::sqrt(14.0), lhs.getY());
    EXPECT_FLOAT_EQ(3.0 / std::sqrt(14.0), lhs.getZ());
}

// =====================================================================================
// Vector3f:: vector opposite.
TEST(TestVector3f, vector_opposite)
{
    Vector3f lhs(1.0, 2.0, 3.0);

    Vector3f result = -lhs;

    EXPECT_FLOAT_EQ(-1.0, result.getX());
    EXPECT_FLOAT_EQ(-2.0, result.getY());
    EXPECT_FLOAT_EQ(-3.0, result.getZ());
}

// =====================================================================================
// Vector3f:: vector division.
TEST(TestVector3f, vector_division)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    float scale = 2.0;

    Vector3f result = lhs / scale;

    EXPECT_FLOAT_EQ(0.5, result.getX());
    EXPECT_FLOAT_EQ(1.0, result.getY());
    EXPECT_FLOAT_EQ(1.5, result.getZ());
}

// =====================================================================================
// Vector3f:: vector addition assignment.
TEST(TestVector3f, vector_addition_assignment)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(4.0, 5.0, 6.0);

    lhs += rhs;

    EXPECT_FLOAT_EQ(5.0, lhs.getX());
    EXPECT_FLOAT_EQ(7.0, lhs.getY());
    EXPECT_FLOAT_EQ(9.0, lhs.getZ());
}

// =====================================================================================
// Vector3f:: vector subtraction assignment.
TEST(TestVector3f, vector_subtraction_assignment)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(4.0, 5.0, 6.0);

    lhs -= rhs;

    EXPECT_FLOAT_EQ(-3.0, lhs.getX());
    EXPECT_FLOAT_EQ(-3.0, lhs.getY());
    EXPECT_FLOAT_EQ(-3.0, lhs.getZ());
}

// =====================================================================================
// Vector3f:: vector scaling assignment.
TEST(TestVector3f, vector_scaling_assignment)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    float scale = 2.0;

    lhs *= scale;

    EXPECT_FLOAT_EQ(2.0, lhs.getX());
    EXPECT_FLOAT_EQ(4.0, lhs.getY());
    EXPECT_FLOAT_EQ(6.0, lhs.getZ());
}

// =====================================================================================
// Vector3f:: term-by-term product assignment.
TEST(TestVector3f, vector_term_by_term_product_assignment)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(4.0, 5.0, 6.0);

    lhs *= rhs;

    EXPECT_FLOAT_EQ(4.0, lhs.getX());
    EXPECT_FLOAT_EQ(10.0, lhs.getY());
    EXPECT_FLOAT_EQ(18.0, lhs.getZ());
}

// =====================================================================================
// Vector3f:: dot product
TEST(TestVector3f, vector_dot_product)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(4.0, 5.0, 6.0);

    float result = lhs.dot(rhs);

    EXPECT_FLOAT_EQ(32.0, result);
}

// =====================================================================================
// Vector3f:: vector division assignment.
TEST(TestVector3f, vector_division_assignment)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    float scale = 2.0;

    lhs /= scale;

    EXPECT_FLOAT_EQ(0.5, lhs.getX());
    EXPECT_FLOAT_EQ(1.0, lhs.getY());
    EXPECT_FLOAT_EQ(1.5, lhs.getZ());
}

// =====================================================================================
// Vector3f:: vector equality/inequality.
TEST(TestVector3f, vector_equality)
{
    Vector3f lhs(1.0, 2.0, 3.0);
    Vector3f rhs(1.0, 2.0, 3.0);

    EXPECT_TRUE(lhs == rhs);

    rhs.setX(2.0);
    EXPECT_FALSE(lhs == rhs);
}


// =====================================================================================
// Vector3f:: vector conversion to array.
TEST(TestVector3f, conversion_to_array)
{
    Vector3f lhs(1.0, 2.0, 3.0);

    std::array<double, 3> result = lhs.to_array();

    EXPECT_FLOAT_EQ(1.0, result[0]);
    EXPECT_FLOAT_EQ(2.0, result[1]);
    EXPECT_FLOAT_EQ(3.0, result[2]);
}


// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
