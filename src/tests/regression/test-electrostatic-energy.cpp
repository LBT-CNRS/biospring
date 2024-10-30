
#include <gtest/gtest.h>
#include <iostream>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "forcefield/constants.hpp"
#include "measure.hpp"

using namespace biospring;

struct TestElectrostaticEnergy : public ::testing::Test
{
    configuration::Configuration config;
    spn::SpringNetwork spn;
    spn::Particle p1, p2;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;
        config.electrostatic.enable = true;
        config.electrostatic.cutoff = 16.0;
        config.electrostatic.dielectric = 1.0;
        SetUpSpn();
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

double expected_electrostatic_energy(const spn::Particle & p1, const spn::Particle & p2, float dielectric)
{
    // Coulomb's constant.
    double k = 8.9875517873681764E9; // N.m2.C-2

    // Charge of the particles (in Coulomb).
    double q1 = p1.getCharge() * forcefield::ELECTRONCHARGE_TO_COULOMB;
    double q2 = p2.getCharge() * forcefield::ELECTRONCHARGE_TO_COULOMB;

    // Distance between the particles (in meter).
    double distance = measure::distance(p1, p2) * forcefield::ANGSTROM_TO_METER;

    // Energy (in Joule).
    double energy = k * q1 * q2 / (dielectric * distance);

    // Conversion to kJ/mol.
    energy = energy * forcefield::AVOGADRO_NUMBER * forcefield::JOULE_TO_KJOULE;

    return energy;
}

// ============================================================================
TEST_F(TestElectrostaticEnergy, energy)
{
    for (float x = 0.01f; x < 5; x += 0.01)
    {
        spn::Particle & lhs = spn.getParticle(0);
        spn::Particle & rhs = spn.getParticle(1);

        rhs.setPosition(Vector3f(x, 0.0, 0.0));
        spn.idleRun();
        spn.computeParticleForces();
        lhs.resetForce();
        rhs.resetForce();

        float distance = lhs.distance(rhs);
        float actual = spn.getElectrostaticEnergy();
        float expected = expected_electrostatic_energy(lhs, rhs, config.electrostatic.dielectric);

        EXPECT_FLOAT_EQ(actual, expected);
    }
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
