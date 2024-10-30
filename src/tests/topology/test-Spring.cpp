#include <gtest/gtest.h>
#include "topology/Spring.hpp"

using namespace biospring;


// ======================================================================================
//
// Tests for the uniqueness of particle ids.
//
// ======================================================================================

// Ensures that the unique id of a particle is not changed when it used in a Spring
TEST(Spring, spring_does_not_copy_particles)
{
    // Create two particles.
    topology::Particle p1;
    topology::Particle p2;

    // Create a spring between the two particles.
    topology::Spring spring(p1, p2, 42.0, 17.0);

    // Check that the particles have the same unique id as before.
    EXPECT_EQ(p1.unique_id(), spring.first().unique_id());
    EXPECT_EQ(p2.unique_id(), spring.second().unique_id());
}


// ======================================================================================
//
// Tests for the `topology::Spring` constructors.
//
// ======================================================================================

TEST(Spring, constructor)
{
    // Create two particles.
    topology::Particle p1;
    topology::Particle p2;

    // Create a spring between the two particles.
    topology::Spring spring(p1, p2, 42.0, 17.0);

    // Check that the spring is initialized correctly.
    EXPECT_FLOAT_EQ(spring.equilibrium(), 42.0);
    EXPECT_FLOAT_EQ(spring.stiffness(), 17.0);
}

TEST(Spring, constructor_with_no_equilibrium_distance_and_no_stiffness)
{
    // Create two particles.
    topology::Particle p1;
    topology::Particle p2;

    // Sets the particle positions (equilibrium distance will be != 0).
    p1.properties().set_position({0.0, 0.0, 0.0});
    p2.properties().set_position({1.0, 0.0, 0.0});

    // Create a spring between the two particles.
    topology::Spring spring(p1, p2);

    // Check that the spring is initialized correctly.
    EXPECT_FLOAT_EQ(spring.equilibrium(), measure::distance(p1, p2));
    EXPECT_FLOAT_EQ(spring.stiffness(), 1.0);
}

TEST(Spring, copy_constructor)
{
    // Create two particles.
    topology::Particle p1;
    topology::Particle p2;

    // Create a spring between the two particles.
    topology::Spring spring(p1, p2, 42.0, 17.0);

    // Copy the spring.
    topology::Spring spring_copy(spring);

    // Check that the spring is initialized correctly.
    EXPECT_EQ(spring.first().unique_id(), p1.unique_id());
    EXPECT_EQ(spring.second().unique_id(), p2.unique_id());
    EXPECT_FLOAT_EQ(spring_copy.equilibrium(), 42.0);
    EXPECT_FLOAT_EQ(spring_copy.stiffness(), 17.0);
}



// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
