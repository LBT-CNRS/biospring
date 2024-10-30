//
// This file contains tests for the `topology::Particle` class.
//

#include <vector>

#include <gtest/gtest.h>
#include "topology/Particle.hpp"
#include "measure.hpp"

// ======================================================================================
// Tests for Particle::unique_id().
// ======================================================================================

TEST(TestParticle, unique_id)
{
    biospring::topology::Particle p1;
    biospring::topology::Particle p2;
    biospring::topology::Particle p3;

    EXPECT_NE(p1.unique_id(), p2.unique_id());
    EXPECT_NE(p1.unique_id(), p3.unique_id());
    EXPECT_NE(p2.unique_id(), p3.unique_id());
}

TEST(TestParticle, unique_id_on_copy)
{
    biospring::topology::Particle p1;
    biospring::topology::Particle p2 = p1.copy();
    biospring::topology::Particle p3 = p1.copy();

    EXPECT_NE(p1.unique_id(), p2.unique_id());
    EXPECT_NE(p1.unique_id(), p3.unique_id());
    EXPECT_NE(p2.unique_id(), p3.unique_id());
}

TEST(TestParticle, unique_id_unchanged_on_assignement)
{
    biospring::topology::Particle p1;
    biospring::topology::Particle p2;
    biospring::topology::Particle p3;

    pid_t p1_uid = p1.unique_id();
    pid_t p2_uid = p2.unique_id();
    pid_t p3_uid = p3.unique_id();

    p1 = p2;
    p2 = p3;
    p3 = p1;

    EXPECT_EQ(p1_uid, p1.unique_id());
    EXPECT_EQ(p2_uid, p2.unique_id());
    EXPECT_EQ(p3_uid, p3.unique_id());
}


TEST(TestParticle, unique_id_in_container)
{
    biospring::topology::Particle p;
    std::vector<biospring::topology::Particle> container;

    container.push_back(p.copy());
    ASSERT_EQ(container.capacity(), 1);

    pid_t expected = container.front().unique_id();

    container.push_back(p.copy());
    ASSERT_EQ(container.capacity(), 2);

    EXPECT_EQ(container.front().unique_id(), expected);
}


// ======================================================================================


TEST(TestParticle, copy)
{
    biospring::topology::Particle p1;
    biospring::topology::Particle p2 = p1;

    EXPECT_EQ(p1.properties(), p2.properties());

    p1.properties().set_position(Vector3f(1.0, 2.0, 3.0));
    EXPECT_NE(p1.properties(), p2.properties());
}


// Ensures that `Particle` meet the requirements to be passed as argument to
// `measure::distance`.
TEST(TestParticle, distance)
{
    biospring::topology::Particle p1;
    biospring::topology::Particle p2;

    ASSERT_NO_THROW(biospring::measure::distance(p1, p2));
    EXPECT_FLOAT_EQ(biospring::measure::distance(p1, p2), 0.0);
}
