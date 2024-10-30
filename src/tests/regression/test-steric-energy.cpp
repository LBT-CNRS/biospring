
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
    float stiffness = 100.0;
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

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
