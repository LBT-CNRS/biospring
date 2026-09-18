
#include <gtest/gtest.h>

#include <cmath>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "forcefield/ForceField.h"
#include "forcefield/constants.hpp"
#include "measure.hpp"

using namespace biospring;

// Not 10, which is the default: a hard-coded decay length would still pass
// every assertion below if the fixture used the value it is hard-coded to.
static const double DECAY_LENGTH = 7.0;

// Reference function that independently calculates the hydrophobic energy
// between two particles (mirrors biospring::forcefield::hydrophobic_energy
// without reusing it, so a bug in the production formula would not be masked
// here).
//
// The factor DECAY_LENGTH is the integral's, not decoration: the force module
// is A exp(-r/L), so the potential it is the gradient of is -A L exp(-r/L).
// Dropping L would leave a quantity in kJ.mol-1.A-1, which is not an energy.
// This function used to multiply by Avogadro's number instead, which made the
// reported energy 6.022e20 times the potential the integrator was actually
// using -- and this test confirmed the arithmetic of that while it was wrong,
// because it reproduced the same mistake independently. Hence
// energy_is_the_integral_of_the_force below: agreeing with a second copy of a
// formula is not the same as the formula being right.
double expected_hydrophobic_energy(const spn::Particle & p1, const spn::Particle & p2)
{
    double distance = measure::distance(p1, p2);
    return -(p1.getHydrophobicity() * p2.getHydrophobicity()) * DECAY_LENGTH *
           exp(-distance / DECAY_LENGTH);
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
        config.hydrophobicity.decaylength = DECAY_LENGTH;
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

// ============================================================================
// The one that has teeth against the class of bug the two tests above missed
// for as long as they existed: an energy and a force that are not the same
// function. Neither of them compares the two, so both passed while the energy
// was off by N_A * 1e-3 = 6.022e20.
//
// BioSpring's sign convention for a "force module" is the one documented in
// electrostatic_shared.h: it is projected on the axis from a particle to its
// NEIGHBOUR, so a positive module means attraction. That makes the module
// +dE/dr rather than the physical radial force -dE/dr -- check it against
// Coulomb if in doubt: like charges give a positive energy falling with r, so
// dE/dr < 0, and the code indeed returns a negative module for them.
//
// The energy is in kJ.mol-1 and the module in Da.A.fs-2, so the comparison
// carries GLOBAL_SPRING_FORCE_CONVERT.
TEST(TestHydrophobicLaw, energy_is_the_integral_of_the_force)
{
    forcefield::ForceField ff;
    ff.setHydrophobicityScale(1.0f);
    ff.setHydrophobicityDecayLength(static_cast<float>(DECAY_LENGTH));

    const float h1 = 0.5f, h2 = 0.3f;
    const double convert = forcefield::GLOBAL_SPRING_FORCE_CONVERT;
    // computeHydrophobicityEnergy returns a float, so a central difference has
    // to balance truncation against cancellation: 1e-3 A leaves only three
    // significant digits in the difference. 0.05 A puts both errors near 1e-5
    // relative, which is why the comparison below is relative too -- and still
    // three orders tighter than anything a missing decay length, a sign flip or
    // an Avogadro factor would survive.
    const double step = 0.05;

    // Across the range the term actually acts over, including well inside the
    // decay length and well outside it.
    for (double r : {1.0, 2.5, 5.0, 7.0, 10.0, 15.0})
    {
        const double eplus = ff.computeHydrophobicityEnergy(h1, h2, r + step);
        const double eminus = ff.computeHydrophobicityEnergy(h1, h2, r - step);
        const double gradient = (eplus - eminus) / (2.0 * step); // kJ.mol-1.A-1
        const double module = ff.computeHydrophobicityForceModule(h1, h2, r);

        EXPECT_NEAR(module, gradient * convert, 1.0e-3 * std::fabs(module))
            << "at r = " << r << " A the force is not the gradient of the reported energy: module "
            << module << " vs dE/dr " << gradient * convert;
    }

    // And the energy has to be an energy: negative (the pair attracts), and
    // decaying by exactly e over one decay length.
    const double at_zero = ff.computeHydrophobicityEnergy(h1, h2, 0.0);
    const double at_one_length = ff.computeHydrophobicityEnergy(h1, h2, DECAY_LENGTH);
    EXPECT_LT(at_zero, 0.0);
    EXPECT_NEAR(at_zero / at_one_length, M_E, 1.0e-5);

    // The decay length must come from the force field, not from a constant in
    // the law: halving it must halve the range, not merely rescale. With L/2 in
    // force, going from L/2 to L is two decay lengths rather than one, so the
    // energy has fallen by e once more -- the ratio is 1/e, the reciprocal of
    // the one just above.
    ff.setHydrophobicityDecayLength(static_cast<float>(DECAY_LENGTH / 2.0));
    const double at_far = ff.computeHydrophobicityEnergy(h1, h2, DECAY_LENGTH);
    const double at_near = ff.computeHydrophobicityEnergy(h1, h2, DECAY_LENGTH / 2.0);
    EXPECT_NEAR(at_far / at_near, 1.0 / M_E, 1.0e-5);
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
