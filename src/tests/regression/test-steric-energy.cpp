
#include <gtest/gtest.h>
#include <iostream>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "forcefield/constants.hpp"

// Reference function that calculates the steric energy between two particles.
// Linear.
float steric_energy_linear(const biospring::spn::Particle & lhs, const biospring::spn::Particle & rhs)
{
    float stiffness = 1.0;
    float equilibrium = lhs.getRadius() + rhs.getRadius();
    float distance = lhs.distance(rhs);
    float distancevar = distance - equilibrium;
    float energy = 0.0;

    if (distancevar < 0.0)
        energy = 0.5 * stiffness * distancevar * distancevar;

    return energy;
}

// Reference function that calculates the steric energy between two particles.
// Amber 12-6.
float steric_energy_amber(const biospring::spn::Particle & lhs, const biospring::spn::Particle & rhs)
{
    float distance = lhs.distance(rhs);

    if (distance < biospring::forcefield::MINIMAL_DISTANCE_VDW_CUTOFF)
        return 0.0;

    float epsilon = sqrt(lhs.getEpsilon() * rhs.getEpsilon());
    float sigma = sqrt(lhs.getRadius() * rhs.getRadius());

    float a = pow(sigma / distance, 12);
    float b = 2 * pow(sigma / distance, 6);
    float V = epsilon * (a - b);

    return V;
}

struct TestStericEnergy : public ::testing::Test
{
    biospring::configuration::Configuration config;
    biospring::spn::SpringNetwork spn;
    biospring::spn::Particle p1, p2;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;
        config.steric.enable = true;
        config.steric.cutoff = 16.0;
    }

    // Sets up the SpringNetwork with two particles.
    void SetUpSpn()
    {
        p1.setPosition(Vector3f(0.0, 0.0, 0.0));
        p2.setPosition(Vector3f(0.0, 0.0, 0.0));

        // Parameters corresponding to the CC in amber.ff.
        p1.setCharge(0.5973);
        p2.setCharge(0.5973);

        p1.setRadius(1.908);
        p2.setRadius(1.908);

        p1.setEpsilon(0.086);
        p2.setEpsilon(0.086);

        p1.setMass(12.01);
        p2.setMass(12.01);

        spn.addParticle(p1);
        spn.addParticle(p2);

        spn.setup(config);
    }
};

struct TestStericEnergyLinear : public TestStericEnergy
{
    void SetUp() override
    {
        TestStericEnergy::SetUp();
        config.steric.mode = "linear";
        SetUpSpn();
    }
};

struct TestStericEnergyAmber : public TestStericEnergy
{
    void SetUp() override
    {
        TestStericEnergy::SetUp();
        config.steric.mode = "lennard-jones-8-6Amber";
        SetUpSpn();
    }
};

// ============================================================================

TEST_F(TestStericEnergyLinear, linear)
{
    for (float x = 0.01f; x < 5; x += 0.01)
    {
        biospring::spn::Particle & lhs = spn.getParticle(0);
        biospring::spn::Particle & rhs = spn.getParticle(1);

        rhs.setPosition(Vector3f(x, 0.0, 0.0));
        spn.idleRun();
        spn.computeParticleForces();
        lhs.resetForce();
        rhs.resetForce();

        float distance = lhs.distance(rhs);
        float actual = spn.getStericEnergy();
        float expected = steric_energy_linear(lhs, rhs);

        EXPECT_FLOAT_EQ(actual, expected);
    }
}

TEST_F(TestStericEnergyAmber, amber)
{
    for (float x = 0.01f; x < 5; x += 0.01)
    {
        biospring::spn::Particle & lhs = spn.getParticle(0);
        biospring::spn::Particle & rhs = spn.getParticle(1);

        rhs.setPosition(Vector3f(x, 0.0, 0.0));
        spn.idleRun();
        spn.computeParticleForces();
        lhs.resetForce();
        rhs.resetForce();

        float distance = lhs.distance(rhs);
        float actual = spn.getStericEnergy();
        float expected = steric_energy_amber(lhs, rhs);

        // A certain degree of incertainty is to be expected.
        if (actual < 1e4)
            EXPECT_LE(abs(actual - expected), 1e-3);

        else if (actual < 1e5)
            EXPECT_LE(abs(actual - expected), 1e-2);

        // Very large energies are expected to be less precise.
        else if (actual < 1e8)
            EXPECT_LE(abs(actual - expected), 1e1);

        else if (actual < 1e15)
            EXPECT_LE(abs(actual - expected), 1e7);
    }
}

// ============================================================================
// Regression tests for unique-pair evaluation (each pair's force/energy is
// computed once, and Newton's third law is applied explicitly instead of
// letting each side recompute the pair independently).
// ============================================================================

struct TestStericEnergyDedup : public ::testing::Test
{
    biospring::configuration::Configuration config;
    biospring::spn::SpringNetwork spn;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;
        config.steric.enable = true;
        config.steric.cutoff = 16.0;
        config.steric.mode = "linear";
    }

    biospring::spn::Particle makeParticle(float x, bool isStatic = false)
    {
        biospring::spn::Particle p;
        p.setPosition(Vector3f(x, 0.0, 0.0));
        p.setCharge(0.5973);
        p.setRadius(1.908);
        p.setEpsilon(0.086);
        p.setMass(12.01);
        p.setStatic(isStatic);
        return p;
    }
};

// Three mutually-visible dynamic particles: every unique pair must contribute
// its full energy exactly once, whichever side (lower or higher id) triggers
// the computation.
TEST_F(TestStericEnergyDedup, three_dynamic_particles_sum_all_unique_pairs)
{
    spn.addParticle(makeParticle(0.0));
    spn.addParticle(makeParticle(1.0));
    spn.addParticle(makeParticle(2.0));
    spn.setup(config);

    spn.idleRun();
    spn.computeParticleForces();

    const auto & a = spn.getParticle(0);
    const auto & b = spn.getParticle(1);
    const auto & c = spn.getParticle(2);

    const float expected = steric_energy_linear(a, b) + steric_energy_linear(a, c) + steric_energy_linear(b, c);
    EXPECT_FLOAT_EQ(spn.getStericEnergy(), expected);
}

// A static particle never re-visits its pairs on its own, so it never
// contributes its own share of the pair energy: the dynamic side must credit
// the full pairwise energy, not half of it (a static neighbor is not double
// counted the way a dynamic one is, since only one side ever computes it).
TEST_F(TestStericEnergyDedup, static_neighbor_contributes_full_pair_energy)
{
    spn.addParticle(makeParticle(0.0, /*isStatic=*/false));
    spn.addParticle(makeParticle(1.0, /*isStatic=*/true));
    spn.setup(config);

    spn.idleRun();
    spn.computeParticleForces();

    const auto & a = spn.getParticle(0);
    const auto & b = spn.getParticle(1);

    const float expected = steric_energy_linear(a, b);
    EXPECT_FLOAT_EQ(spn.getStericEnergy(), expected);

    // The static particle must never receive a force: it is never integrated
    // or reset, so any stray write would silently accumulate across steps.
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
