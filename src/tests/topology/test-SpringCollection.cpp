#include <gtest/gtest.h>

#include "topology/Particle.hpp"
#include "topology/ParticleCollection.hpp"
#include "topology/Spring.hpp"
#include "topology/SpringCollection.hpp"

using namespace biospring;


// Returns a collection of particles of size `n`.
static topology::ParticleCollection generate_particles(size_t n)
{
    topology::ParticleCollection particles;
    for (size_t i = 0; i < n; ++i)
        particles.push_back(topology::Particle());
    return particles;
}

TEST(SpringCollection, plus_operator)
{
    topology::ParticleCollection particles = generate_particles(3);
    topology::SpringCollection container(particles);
    topology::SpringCollection container2(particles);

    container.add_spring(particles[0], particles[1], 1.0, 1.0);
    container2.add_spring(particles[1], particles[2], 1.0, 1.0);

    topology::SpringCollection container3 = container + container2;

    EXPECT_EQ(container3.size(), 2);
    EXPECT_TRUE(container3.exists(particles[0], particles[1]));
    EXPECT_TRUE(container3.exists(particles[1], particles[2]));
    EXPECT_FALSE(container3.exists(particles[0], particles[2]));
}

TEST(SpringCollection, plus_assignment_operator)
{
    topology::ParticleCollection particles = generate_particles(3);
    topology::SpringCollection container(particles);
    topology::SpringCollection container2(particles);

    container.add_spring(particles[0], particles[1], 1.0, 1.0);
    container2.add_spring(particles[1], particles[2], 1.0, 1.0);

    container += container2;

    EXPECT_EQ(container.size(), 2);
    EXPECT_TRUE(container.exists(particles[0], particles[1]));
    EXPECT_TRUE(container.exists(particles[1], particles[2]));
    EXPECT_FALSE(container.exists(particles[0], particles[2]));
}

TEST(SpringCollection, plus_assignment_operator_fails_if_particles_are_not_the_same)
{
    topology::ParticleCollection particles = generate_particles(3);
    topology::ParticleCollection particles2 = generate_particles(3);
    topology::SpringCollection container(particles);
    topology::SpringCollection container2(particles2);

    container.add_spring(particles[0], particles[1], 1.0, 1.0);
    container2.add_spring(particles2[1], particles2[2], 1.0, 1.0);

    EXPECT_THROW(container += container2, std::invalid_argument);
}

TEST(SpringCollection, add_spring)
{
    topology::ParticleCollection particles = generate_particles(3);
    topology::SpringCollection container(particles);

    container.add_spring(particles[0], particles[1], 1.0, 1.0);

    EXPECT_EQ(container.size(), 1);
    EXPECT_EQ(container[0].unique_id().first, particles[0].unique_id());
    EXPECT_EQ(container[0].unique_id().second, particles[1].unique_id());
}


TEST(SpringCollection, add_spring_fails_if_particles_are_not_in_collection)
{
    topology::ParticleCollection particles = generate_particles(3);
    topology::SpringCollection container(particles);

    topology::Particle p1, p2;

    EXPECT_THROW(container.add_spring(p1, p2, 1.0, 1.0), topology::InvalidParticleException);

    particles.push_back(p1);
    particles.push_back(p2);

    // Still fails: the particles are not in the topology: copies are.
    EXPECT_THROW(container.add_spring(p1, p2, 1.0, 1.0), topology::InvalidParticleException);
}


TEST(SpringCollection, add_spring_fails_if_particles_are_the_same)
{
    topology::ParticleCollection particles;
    topology::SpringCollection container(particles);

    particles.push_back(topology::Particle());

    EXPECT_THROW(container.add_spring(particles[0], particles[0], 1.0, 1.0), topology::SelfSpringException);
}


TEST(SpringCollection, add_spring_fails_if_spring_already_exists)
{
    topology::ParticleCollection particles = generate_particles(3);
    topology::SpringCollection container(particles);

    container.add_spring(particles[0], particles[1], 1.0, 1.0);

    EXPECT_THROW(container.add_spring(particles[0], particles[1], 1.0, 1.0), topology::SpringAlreadyExistsException);
}

TEST(SpringCollection, exists)
{
    topology::ParticleCollection particles = generate_particles(3);
    topology::SpringCollection container(particles);

    container.add_spring(particles[0], particles[1], 1.0, 1.0);

    EXPECT_TRUE(container.exists(particles[0], particles[1]));
    EXPECT_TRUE(container.exists(particles[1], particles[0]));
    EXPECT_FALSE(container.exists(particles[0], particles[2]));
    EXPECT_FALSE(container.exists(particles[2], particles[0]));
    EXPECT_FALSE(container.exists(particles[1], particles[2]));
    EXPECT_FALSE(container.exists(particles[2], particles[1]));
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
