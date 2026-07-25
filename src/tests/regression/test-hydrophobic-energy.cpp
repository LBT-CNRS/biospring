
#include <gtest/gtest.h>

#include <cmath>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "forcefield/constants.hpp"
#include "measure.hpp"

using namespace biospring;

// Reference function that independently calculates the hydrophobic energy
// between two particles (mirrors biospring::forcefield::hydrophobic_energy
// without reusing it, so a bug in the production formula would not be masked
// here).
double expected_hydrophobic_energy(const spn::Particle & p1, const spn::Particle & p2)
{
    double distance = measure::distance(p1, p2);
    double energy = -(p1.getHydrophobicity() * p2.getHydrophobicity()) * exp(-distance);
    energy = energy * forcefield::AVOGADRO_NUMBER; // J/mol
    energy = energy * 1.0E-3;                      // kJ/mol
    return energy;
}

struct TestHydrophobicEnergy : public ::testing::Test
{
    configuration::Configuration config;
    spn::SpringNetwork spn;
    spn::Particle p1, p2;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;
        config.hydrophobicity.enable = true;
        config.hydrophobicity.cutoff = 16.0;
    }

    void SetUpSpn()
    {
        p1.setPosition(Vector3f(0.0, 0.0, 0.0));
        p2.setPosition(Vector3f(2.0, 0.0, 0.0));

        p1.setHydrophobicity(0.5);
        p2.setHydrophobicity(0.3);

        spn.addParticle(p1);
        spn.addParticle(p2);

        spn.setup(config);
    }
};

// ============================================================================
// Two dynamic particles: both independently visit the pair (subject to the
// unique-pair dedup rule), so the total energy must equal the full pairwise
// energy exactly once.
TEST_F(TestHydrophobicEnergy, two_dynamic_particles)
{
    SetUpSpn();

    spn.idleRun();
    spn.computeParticleForces();

    const auto & a = spn.getParticle(0);
    const auto & b = spn.getParticle(1);

    const float expected = expected_hydrophobic_energy(a, b);
    EXPECT_FLOAT_EQ(spn.getHydrophobicEnergy(), expected);
}

// A static particle never re-visits its pairs on its own, so it never
// contributes its own share of the pair energy: the dynamic side must credit
// the full pairwise energy, not half of it.
TEST_F(TestHydrophobicEnergy, static_neighbor_contributes_full_pair_energy)
{
    p2.setStatic(true);
    SetUpSpn();

    spn.idleRun();
    spn.computeParticleForces();

    const auto & a = spn.getParticle(0);
    const auto & b = spn.getParticle(1);

    const float expected = expected_hydrophobic_energy(a, b);
    EXPECT_FLOAT_EQ(spn.getHydrophobicEnergy(), expected);

    EXPECT_FLOAT_EQ(b.getForce().getX(), 0.0f);
    EXPECT_FLOAT_EQ(b.getForce().getY(), 0.0f);
    EXPECT_FLOAT_EQ(b.getForce().getZ(), 0.0f);
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
