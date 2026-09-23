
#include <gtest/gtest.h>
#include <iostream>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLinear.h"
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

// A static particle never re-visits its pairs on its own, so it never
// contributes its own share of the pair energy: the dynamic side must credit
// the full pairwise energy, not half of it.
TEST(TestElectrostaticEnergyStatic, static_neighbor_contributes_full_pair_energy)
{
    configuration::Configuration config;
    config.sim.nbsteps = 1;
    config.sim.timestep = 0.01;
    config.electrostatic.enable = true;
    config.electrostatic.cutoff = 16.0;
    config.electrostatic.dielectric = 1.0;

    spn::Particle p1, p2;
    p1.setPosition(Vector3f(0.0, 0.0, 0.0));
    p2.setPosition(Vector3f(2.0, 0.0, 0.0));
    p1.setCharge(0.5973);
    p2.setCharge(-0.5973);
    p2.setStatic(true);

    spn::SpringNetwork spn;
    spn.addParticle(p1);
    spn.addParticle(p2);
    spn.setup(config);

    spn.idleRun();
    spn.computeParticleForces();

    const auto & a = spn.getParticle(0);
    const auto & b = spn.getParticle(1);

    const float expected = expected_electrostatic_energy(a, b, config.electrostatic.dielectric);
    EXPECT_FLOAT_EQ(spn.getElectrostaticEnergy(), expected);

    EXPECT_FLOAT_EQ(b.getForce().getX(), 0.0f);
    EXPECT_FLOAT_EQ(b.getForce().getY(), 0.0f);
    EXPECT_FLOAT_EQ(b.getForce().getZ(), 0.0f);
}

// -- Main function  ----------------------------------------------------------

// A distance-dependent dielectric changes the exponent, not just a constant.
// With epsilon(r) = e0*r the energy goes as 1/r^2, so -dE/dr goes as 2/r^3 --
// twice what substituting e0*r into the ordinary 1/r^2 force expression
// gives. Getting that wrong halves every electrostatic force while leaving
// every reported energy correct, which no energy check would ever reveal.
//
// Energy and force module are returned in different unit systems (kJ/mol
// against Da.A.fs-2), so the check is expressed as a RATIO between the two
// modes: whatever the conversion factor is, it is the same on both sides and
// cancels. The constant dielectric is the reference, its gradient being
// long-established.
TEST(TestElectrostaticEnergyStatic, distance_dependent_dielectric_force_matches_its_own_gradient)
{
    biospring::forcefield::ForceFieldElectrostaticCoulombAndStericLinear ff;
    const float q1 = 0.5f, q2 = -0.7f, h = 1e-3f;

    auto ratio = [&](float r) {
        const float fd = -(ff.computeElectrostaticEnergy(q1, q2, r + h) -
                           ff.computeElectrostaticEnergy(q1, q2, r - h)) / (2.0f * h);
        return fd / ff.computeElectrostaticForceModule(q1, q2, r);
    };

    ff.setDielectric(4.0f);
    ff.setDistanceDependentDielectric(false);
    const float reference = ratio(5.0f);

    ff.setDistanceDependentDielectric(true);
    for (float r : {3.0f, 5.0f, 8.0f})
        EXPECT_NEAR(ratio(r) / reference, 1.0f, 0.02f) << "at r = " << r;

    // ... and the energy really does fall off faster: doubling the distance
    // divides it by four, not by two.
    EXPECT_NEAR(ff.computeElectrostaticEnergy(q1, q2, 3.0f) /
                ff.computeElectrostaticEnergy(q1, q2, 6.0f), 4.0f, 0.05f);
}

int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
