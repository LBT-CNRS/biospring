#include <gtest/gtest.h>

#include <vector>

#include "IO/io.h"
#include "SpringNetwork.h"
#include "measure.hpp"
#include "topology/Spring.hpp"
#include "topology/Topology.hpp"

using namespace biospring;

// ====================================================================================
//
// Helper functions.
//
// ====================================================================================

static auto generate_random_particles(size_t n);
// Returns a random 3D-position.
// Each coordinate is drawn from a uniform distribution on the interval [-100, 100].
static std::array<double, 3> random_position()
{
    std::array<double, 3> pos;
    for (double & x : pos)
        x = 200.0 * (double)rand() / RAND_MAX - 100.0;
    return pos;
}

// Generates a collection of particles with random positions.
static auto generate_random_particles(size_t n)
{
    std::vector<topology::Particle> particles(n);
    for (topology::Particle & p : particles)
        p.set_position(random_position());
    return particles;
}

// ====================================================================================
//
// Test for I/O functions.
//
// ====================================================================================
TEST(Topology, readTopology)
{
    std::string path = "../data/model.pdb";
    topology::Topology top = io::readTopology(path);
    EXPECT_EQ(top.number_of_particles(), 2814);
    EXPECT_EQ(top.number_of_springs(), 0);
    EXPECT_EQ(top.particles().size(), top.springs().particles().size());
}

// ====================================================================================
//
// Merge topologies.
//
// ====================================================================================

TEST(Topology, merge_static)
{
    topology::Topology top1, top2;
    top1.add_particles(generate_random_particles(100));
    top2.add_particles(generate_random_particles(100));

    top1.add_springs_from_cutoff(10.0);
    top2.add_springs_from_cutoff(10.0);

    ASSERT_GE(top1.number_of_springs(), 1);
    ASSERT_GE(top2.number_of_springs(), 1);

    size_t expected_number_of_particles = top1.number_of_particles() + top2.number_of_particles();
    size_t expected_number_of_springs = top1.number_of_springs() + top2.number_of_springs();

    topology::Topology merged = topology::Topology::merge(top1, top2);
    EXPECT_EQ(merged.number_of_particles(), expected_number_of_particles);
    EXPECT_EQ(merged.number_of_springs(), expected_number_of_springs);
}

TEST(Topology, merge)
{
    topology::Topology top1, top2;
    top1.add_particles(generate_random_particles(100));
    top2.add_particles(generate_random_particles(100));

    top1.add_springs_from_cutoff(10.0);
    top2.add_springs_from_cutoff(10.0);

    ASSERT_GE(top1.number_of_springs(), 1);
    ASSERT_GE(top2.number_of_springs(), 1);

    size_t expected_number_of_particles = top1.number_of_particles() + top2.number_of_particles();
    size_t expected_number_of_springs = top1.number_of_springs() + top2.number_of_springs();

    topology::Topology merged = top1.merge(top2);
    EXPECT_EQ(merged.number_of_particles(), expected_number_of_particles);
    EXPECT_EQ(merged.number_of_springs(), expected_number_of_springs);
}

// ====================================================================================
//
// Add particles
//
// ====================================================================================

TEST(Topology, add_particle_one_by_one)
{
    topology::Topology top;
    topology::Particle p1, p2;

    top.particles().push_back(p1);
    EXPECT_EQ(top.number_of_particles(), 1);

    top.add_particle(p2);
    EXPECT_EQ(top.number_of_particles(), 2);
}

TEST(Topology, add_particle_container)
{
    topology::Topology top;
    std::vector<topology::Particle> particles(2);

    top.add_particles(particles);
    EXPECT_EQ(top.number_of_particles(), 2);
}

TEST(Topology, add_particle_initializer_list)
{
    topology::Topology top;
    topology::Particle p1, p2;

    top.add_particles({p1, p2});
    EXPECT_EQ(top.number_of_particles(), 2);
}

TEST(Topology, add_spring)
{
    topology::Topology top;
    topology::Particle p1, p2;
    top.add_particles({p1, p2});

    top.add_spring(top.get_particle(0), top.get_particle(1), 1.0, 1.0);

    EXPECT_EQ(top.number_of_springs(), 1);
    EXPECT_FLOAT_EQ(top.get_spring(0).equilibrium(), 1.0);
    EXPECT_FLOAT_EQ(top.get_spring(0).stiffness(), 1.0);
}

// Checks that springs are added between all particles within a cutoff distance.
TEST(Topology, add_springs_from_cutoff)
{
    topology::Topology top;
    top.add_particles(generate_random_particles(1000));
    top.add_springs_from_cutoff(10.0);
    ASSERT_GE(top.number_of_springs(), 1);

    for (size_t i = 0; i < top.number_of_particles(); ++i)
    {
        for (size_t j = i + 1; j < top.number_of_particles(); ++j)
        {
            const topology::Particle & p1 = top.get_particle(i);
            const topology::Particle & p2 = top.get_particle(j);
            double distance = measure::distance(p1, p2);

            if (distance < 10.0)
                EXPECT_TRUE(top.has_spring_between(p1, p2));

            // else
            //     EXPECT_FALSE(top.has_spring_between(p1, p2));
        }
    }
}

TEST(Topology, to_spring_network)
{
    topology::Topology top;
    top.add_particles(generate_random_particles(1000));
    top.add_springs_from_cutoff(10.0);
    ASSERT_GE(top.number_of_springs(), 1);

    spn::SpringNetwork spn;
    top.to_spring_network(spn);

    EXPECT_EQ(spn.getNumberOfParticles(), top.number_of_particles());
    EXPECT_EQ(spn.getNumberOfSprings(), top.number_of_springs());
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
