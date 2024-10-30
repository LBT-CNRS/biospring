#include "topology/Particle.hpp"
#include "topology/ParticleCollection.hpp"
#include <gtest/gtest.h>

using namespace biospring;

TEST(ParticleCollection, find)
{
    topology::ParticleCollection container;
    topology::Particle p1, p2;

    container.push_back(p1);
    container.push_back(p2);

    EXPECT_EQ(container.find(container[0]), 0);
    EXPECT_EQ(container.find(container[1]), 1);
    EXPECT_EQ(container.find(p1), topology::ParticleCollection::PARTICLE_NOT_FOUND);
}

TEST(ParticleCollection, contains)
{
    topology::ParticleCollection container;
    topology::Particle p1;

    container.push_back(p1);

    EXPECT_TRUE(container.contains(container[0]));
    EXPECT_FALSE(container.contains(topology::Particle()));
}

// Tests that we can retrieve a particle from a topology using its unique id.
TEST(ParticleCollection, access_to_particles_using_uid)
{
    topology::ParticleCollection container;
    topology::Particle p1, p2;

    container.push_back(p1);
    container.push_back(p2);

    topology::pid_t expected_1 = container[0].unique_id();
    topology::pid_t expected_2 = container[1].unique_id();

    EXPECT_EQ(expected_1, container.at_uid(expected_1).unique_id());
    EXPECT_EQ(expected_2, container.at_uid(expected_2).unique_id());
}

TEST(ParticleCollection, addParticleMakesCopies)
{
    topology::ParticleCollection container;
    topology::Particle p1, p2;

    container.push_back(p1);
    container.push_back(p2);

    EXPECT_EQ(container.size(), 2);
    EXPECT_NE(container[0].unique_id(), p1.unique_id());
    EXPECT_NE(container[1].unique_id(), p2.unique_id());
}

// ====================================================================================
//
// Tests for append methods.
//
// ====================================================================================

TEST(ParticleCollection, plus_operator)
{
    topology::ParticleCollection container;
    topology::Particle p1, p2, p3;

    container += p1;
    container += p2;
    container += p3;

    EXPECT_EQ(container.size(), 3);

    topology::ParticleCollection container2 = container + p1;
    EXPECT_EQ(container2.size(), 4);
}

TEST(ParticleCollection, append_with_initializer_list)
{
    topology::ParticleCollection container;
    topology::Particle p1, p2, p3;

    container.push_back({p1, p2, p3});
    EXPECT_EQ(container.size(), 3);
}

TEST(ParticleCollection, append_with_vector)
{
    topology::ParticleCollection container;
    topology::Particle p1, p2, p3;
    std::vector<topology::Particle> particles = {p1, p2, p3};

    container.push_back(particles);
    EXPECT_EQ(container.size(), 3);
}

TEST(ParticleCollection, append_with_ParticleCollection)
{
    topology::ParticleCollection container;
    topology::Particle p1, p2, p3;
    topology::ParticleCollection particles = {p1, p2, p3};

    container.push_back(particles);
    EXPECT_EQ(container.size(), 3);
}

// ====================================================================================
//
// Tests for construction methods.
//
// ====================================================================================

TEST(ParticleCollection, construct_with_no_arguments)
{
    topology::ParticleCollection container;
    EXPECT_EQ(container.size(), 0);
}

TEST(ParticleCollection, construct_with_vector)
{
    topology::Particle p1, p2, p3;
    std::vector<topology::Particle> particles = {p1, p2, p3};
    topology::ParticleCollection container(particles);

    EXPECT_EQ(container.size(), 3);
}

TEST(ParticleCollection, construct_with_initializer_list)
{
    topology::Particle p1, p2, p3;
    topology::ParticleCollection container = {p1, p2, p3};

    EXPECT_EQ(container.size(), 3);
}

// ====================================================================================
//
// Tests for particle properties
//
// ====================================================================================

TEST(ParticleCollection, set_particle_property)
{
    topology::ParticleCollection particles;
    topology::Particle p1, p2;

    p1.properties().set_mass(42);
    p2.properties().set_mass(17);

    particles.push_back(p1);
    particles.push_back(p2);

    auto masses = particles.masses();

    EXPECT_FLOAT_EQ(masses[0], 42);
    EXPECT_FLOAT_EQ(masses[1], 17);

    masses = {12, 13};
    particles.set_masses(masses);

    EXPECT_FLOAT_EQ(particles[0].properties().mass(), 12.0);
    EXPECT_FLOAT_EQ(particles[1].properties().mass(), 13.0);

    // Should work with single value.
    particles.set_masses(42);
    EXPECT_FLOAT_EQ(particles[0].properties().mass(), 42.0);
    EXPECT_FLOAT_EQ(particles[1].properties().mass(), 42.0);
}

TEST(ParticleCollection, set_particle_property_Vector3f)
{
    topology::ParticleCollection particles;
    particles.push_back({topology::Particle(), topology::Particle()});
    ASSERT_EQ(particles.size(), 2);

    auto positions = particles.positions();
    EXPECT_FLOAT_EQ(positions[0].getX(), 0.0);
    EXPECT_FLOAT_EQ(positions[0].getY(), 0.0);
    EXPECT_FLOAT_EQ(positions[0].getZ(), 0.0);
    EXPECT_FLOAT_EQ(positions[1].getX(), 0.0);
    EXPECT_FLOAT_EQ(positions[1].getY(), 0.0);
    EXPECT_FLOAT_EQ(positions[1].getZ(), 0.0);

    positions = {{1, 2, 3}, {4, 5, 6}};

    particles.set_positions(positions);
    EXPECT_FLOAT_EQ(particles[0].getX(), 1.0);
    EXPECT_FLOAT_EQ(particles[0].getY(), 2.0);
    EXPECT_FLOAT_EQ(particles[0].getZ(), 3.0);
    EXPECT_FLOAT_EQ(particles[1].getX(), 4.0);
    EXPECT_FLOAT_EQ(particles[1].getY(), 5.0);
    EXPECT_FLOAT_EQ(particles[1].getZ(), 6.0);

    Vector3f position = {0.0, 0.0, 0.0};
    particles.set_positions(position);
    positions = particles.positions();
    EXPECT_FLOAT_EQ(positions[0].getX(), 0.0);
    EXPECT_FLOAT_EQ(positions[0].getY(), 0.0);
    EXPECT_FLOAT_EQ(positions[0].getZ(), 0.0);
    EXPECT_FLOAT_EQ(positions[1].getX(), 0.0);
    EXPECT_FLOAT_EQ(positions[1].getY(), 0.0);
    EXPECT_FLOAT_EQ(positions[1].getZ(), 0.0);
}

TEST(ParticleCollection, set_particle_string_property)
{
    topology::ParticleCollection particles;
    particles.push_back({topology::Particle(), topology::Particle()});
    ASSERT_EQ(particles.size(), 2);

    std::vector<std::string> names = {"foo", "bar"};
    particles.set_names(names);
    EXPECT_EQ(particles[0].properties().name(), names[0]);
    EXPECT_EQ(particles[1].properties().name(), names[1]);

    particles.set_residue_names("baz");
    EXPECT_EQ(particles[0].properties().residue_name(), "baz");
    EXPECT_EQ(particles[1].properties().residue_name(), "baz");
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
