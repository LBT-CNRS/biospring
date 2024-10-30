#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "Particle.h"
#include "measure.hpp"
#include "nsearch.hpp"

using Particle = biospring::spn::Particle;

bool is_position_duplicate(const Particle & candidate, const std::vector<Particle> & system);
bool has_position_duplicate(const std::vector<Particle> & system);
std::array<double, 3> random_position();
std::vector<Particle> generate_random_particles(size_t n);
std::vector<Particle> generate_particles(size_t n);
std::vector<Particle> generate_particle_group(size_t n, const std::array<double, 3> & position);
std::vector<Particle> generate_particle_groups(size_t n);

// =====================================================================================
//
// Test for `NeighborSearch2` class.
//
// =====================================================================================

TEST(TestNeighborSearch2, NeighbourSearchingOneNeighbor)
{
    // 10 groups of 10 particles spaced 10 units apart.
    const auto particles = generate_particle_groups(10);
    double cutoff = 1.0;
    biospring::nsearch::NeighborSearch2 ns(particles, cutoff);

    for (size_t i = 0; i < particles.size() - 1; i++)
    {
        const auto & neighbors = ns.get_neighbors(particles[i]);
        // Groups have 10 particles so each particle should have 9 neighbors.
        ASSERT_EQ(neighbors.size(), 9);
    }
}

// =====================================================================================
//
// Test for `NeighborSearchDynamic` class.
//
// =====================================================================================
TEST(TestNeighborSearchDynamic, is_dynamic)
{
    std::vector<Particle> particles = generate_particles(10);
    double cutoff = 2.0;
    biospring::nsearch::NeighborSearchDynamic ns(particles, cutoff);

    // Particle 0 has only one neighbor, particle 1.
    ASSERT_LE(biospring::measure::distance(particles[0], particles[1]), cutoff);

    std::vector<size_t> neighbors = ns.get_neighbors(0);
    ASSERT_EQ(neighbors.size(), 1);
    EXPECT_EQ(neighbors[0], 1);

    // Changes Particle 0's coordinates to be neighbor with particle 9 (and only
    // particle 9).
    particles[0].setPosition(Vector3f(10, 10, 10));
    ns.update();

    neighbors = ns.get_neighbors(0);
    ASSERT_EQ(neighbors.size(), 1);
    EXPECT_EQ(neighbors[0], 9);
}

// =====================================================================================
//
// Test for `NeighborSearch` class.
//
// =====================================================================================

TEST(TestNeighborSearch, NeighbourSearchingOneNeighbor)
{
    // 10 groups of 10 particles spaced 10 units apart.
    const auto particles = generate_particle_groups(10);
    double cutoff = 1.0;
    biospring::nsearch::NeighborSearch ns(particles, cutoff);

    for (size_t i = 0; i < particles.size() - 1; i++)
    {
        const auto & neighbors = ns.get_neighbors(particles[i]);
        // Groups have 10 particles so each particle should have 9 neighbors.
        EXPECT_EQ(neighbors.size(), 9);
    }
}

TEST(TestNeighborSearch, NeighbourSearchingNoNeighbor)
{
    const auto particles = generate_particles(10);
    biospring::nsearch::NeighborSearch ns(particles, 0.5);

    for (size_t i = 0; i < particles.size(); i++)
    {
        const auto & neighbors = ns.get_neighbors(particles[i]);
        EXPECT_TRUE(neighbors.empty());
    }
}

TEST(TestNeighborSearch, NeighbourSearchingRandom)
{
    const auto particles = generate_random_particles(1000);

    // makes sure there are no duplicates in the system.
    ASSERT_FALSE(has_position_duplicate(particles));

    double cutoff = 10.0;
    biospring::nsearch::NeighborSearch ns(particles, cutoff);

    for (size_t i = 0; i < particles.size() - 1; i++)
    {
        const auto & neighbors = ns.get_neighbors(particles[i]);
        for (size_t j = i + 1; j < particles.size(); j++)
        {
            const auto & p1 = particles[i];
            const auto & p2 = particles[j];
            double distance = biospring::measure::distance(p1, p2);

            bool neighbors_contains_j = std::find(neighbors.begin(), neighbors.end(), j) != neighbors.end();
            if (distance < cutoff)
                EXPECT_TRUE(neighbors_contains_j);
            else
                EXPECT_FALSE(neighbors_contains_j);
        }
    }
}

// =====================================================================================
//
// Test for `NeighborSearchBase` class.
//
// =====================================================================================

// Excepts a failure when the particle list is empty.
TEST(TestNeighborSearch, InitFailsEmptyParticleList)
{
    const auto particles = generate_particles(0);
    EXPECT_THROW(biospring::nsearch::NeighborSearch(particles, 1.0), std::invalid_argument);
}

// Excepts a failure when the cutoff distance is negative.
TEST(TestNeighborSearch, InitFailsNegativeCutoff)
{
    const auto particles = generate_particles(10);
    EXPECT_THROW(biospring::nsearch::NeighborSearch(particles, -1.0);, std::invalid_argument);
}

// Excepts a failure when the cutoff distance is too close to zero.
TEST(TestNeighborSearch, InitFailsCutoffIsZero)
{
    const auto particles = generate_particles(10);
    EXPECT_THROW(biospring::nsearch::NeighborSearch(particles, 1e-8), std::invalid_argument);
}

// =====================================================================================
//
// Helper functions.
//
// =====================================================================================

bool is_position_duplicate(const Particle & candidate, const std::vector<Particle> & system)
{
    for (size_t i = 0; i < system.size(); i++)
    {
        const auto & existing = system.at(i);
        const auto & candidate_position = candidate.getPosition();
        const auto & existing_position = existing.getPosition();

        // Checks that candidate and existing are not the same particle.
        if (&candidate == &existing)
            continue;

        if (std::abs(existing_position.getX() - candidate_position.getX()) < 1e-6 &&
            std::abs(existing_position.getY() - candidate_position.getY()) < 1e-6 &&
            std::abs(existing_position.getZ() - candidate_position.getZ()) < 1e-6)
        {
            return true; // The position is a duplicate
        }
    }
    return false; // The position is not a duplicate
}

bool has_position_duplicate(const std::vector<Particle> & system)
{
    for (size_t i = 0; i < system.size(); i++)
    {
        if (is_position_duplicate(system.at(i), system))
            return true;
    }
    return false;
}

// Returns a random 3D-position.
// Each coordinate is drawn from a uniform distribution on the interval [-100, 100].
std::array<double, 3> random_position()
{
    std::array<double, 3> pos;
    for (double & x : pos)
        x = 200.0 * (double)rand() / RAND_MAX - 100.0;
    return pos;
}

// Generates a collection of particles with random positions.
std::vector<Particle> generate_random_particles(size_t n)
{
    std::vector<Particle> particles(n);
    for (Particle & p : particles)
        p.setPosition(random_position());
    return particles;
}

// Generates a list of particles of size `n`.
// Particle coordinates are set to (0, 0, 0), (1, 1, 1), ..., (n - 1, n - 1, n - 1).
std::vector<Particle> generate_particles(size_t n)
{
    std::vector<Particle> particles(n);
    for (size_t i = 0; i < n; i++)
        particles[i].setPosition(Vector3f(i, i, i));
    return particles;
}

// Generate a list of particles of size `n` with the same coordinates.
std::vector<Particle> generate_particle_group(size_t n, const std::array<double, 3> & position)
{
    std::vector<Particle> particles(n);
    for (size_t i = 0; i < n; i++)
        particles[i].setPosition(position);
    return particles;
}

// Generates a list of n particle groups, spaced 10 units apart.
std::vector<Particle> generate_particle_groups(size_t n)
{
    std::vector<Particle> particles;
    for (size_t i = 0; i < n; i++)
    {
        std::array<double, 3> position = {10 * double(i), 0.0, 0.0};
        const auto group = generate_particle_group(10, position);
        particles.insert(particles.end(), group.begin(), group.end());
    }
    return particles;
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
