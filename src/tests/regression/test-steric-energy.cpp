
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
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
    // rmin of the pair. The radii in a .ff file are AMBER's R* = rmin/2, so two
    // of them ADD; see biospring_lorentz_berthelot_radius in steric_shared.h.
    float sigma = lhs.getRadius() + rhs.getRadius();

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
        config.steric.mode = "lennard-jones-12-6Amber";
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

        // RELATIVE, and that is not a weakening of this test but the only bound
        // it can carry. `expected` above is computed in double through pow();
        // the network computes in float. Over 0.01 to 5 A a 12-6 potential
        // spans fourteen orders of magnitude, so the difference between the two
        // is a float ULP, which is a FRACTION of the value and not a number of
        // kJ/mol.
        //
        // The ladder of absolute bounds this replaces -- 1e-3 below 1e4 kJ/mol,
        // then 1e-2, 1e1, 1e7 -- was one ULP at the top of its range: at
        // x = 0.21 A the energy is 1.11e+14 kJ/mol, where a float's ULP is
        // 8.4e+06, against a tolerance of 1e+07. It therefore passed or failed
        // on which way the last bit of pow() happened to round, and it did
        // round differently on Linux and on macOS: green here and red in CI for
        // three branches, on a change that moved no steric physics at all.
        //
        // 1e-6 is about ten times the worst relative gap measured over the
        // whole sweep (1.0e-7, at x = 1.22 A) and still far inside float's own
        // seven digits. The absolute floor keeps it meaningful where the
        // potential has decayed to nothing and a relative bound would not be.
        EXPECT_NEAR(actual, expected, std::max(1.0e-6f * std::abs(expected), 1.0e-6f));
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

// ============================================================================
// THE RADIUS COMBINING RULE, which belongs to the force-field FILE and not to
// the law. data/forcefield/amber*.ff stores AMBER's R* = rmin/2, whose radii
// ADD; CAonlyLewitt.ff stores a CA-CA contact distance already, whose radii are
// MEANED. Getting this backwards does not shift the well, it destroys it: two
// AMBER carbons meaned sit at 1.91 A instead of 3.82 A, so a stack of
// nucleobases 3.4 A apart has no minimum at all, and two of Levitt's residues
// summed sit at 10.4 A where a real protein's closest non-bonded CA pair is
// 3.7 A. Each law therefore defaults to the rule of its own parameter set, and
// steric.radiusrule overrides it.
// ============================================================================

struct TestRadiusRule : public ::testing::Test
{
    // A pair of UNEQUAL radii, because every rule agrees on equal ones: the
    // geometric mean of r and r is r, which is exactly why the defect survived.
    static constexpr float R1 = 1.0f;
    static constexpr float R2 = 4.0f;
    static constexpr float E1 = 0.5f;
    static constexpr float E2 = 0.5f;

    // Each rule's own pair minimum for (R1, R2).
    static constexpr float SUM = 5.0f;              // 1 + 4
    static constexpr float GEOMETRIC = 2.0f;        // sqrt(1 * 4)
    static constexpr float ARITHMETIC = 2.5f;       // (1 + 4) / 2

    biospring::configuration::Configuration config;
    biospring::spn::SpringNetwork spn;

    void SetUp() override
    {
        config = biospring::configuration::defaultConfiguration();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;
        config.steric.enable = true;
        config.steric.cutoff = 16.0;
    }

    // The steric energy the network reports for two particles `distance` apart.
    float energyAt(const std::string & mode, const std::string & rule, float distance)
    {
        biospring::spn::SpringNetwork net;
        for (int i = 0; i < 2; ++i)
        {
            biospring::spn::Particle p;
            p.setPosition(Vector3f(i * distance, 0.0, 0.0));
            p.setCharge(0.0);
            p.setRadius(i == 0 ? R1 : R2);
            p.setEpsilon(i == 0 ? E1 : E2);
            p.setMass(12.0);
            net.addParticle(p);
        }
        biospring::configuration::Configuration c = config;
        c.steric.mode = mode;
        c.steric.radiusrule = rule;
        net.setup(c);
        net.idleRun();
        net.computeParticleForces();
        return net.getStericEnergy();
    }

    // Where the reported energy is lowest, scanned finely. This is the pair
    // minimum as the RUNNING CODE sees it -- not as a reimplementation of the
    // formula says it should be, which is how the defect hid in the first place.
    float minimumOf(const std::string & mode, const std::string & rule)
    {
        float best = 1e30f, at = 0.0f;
        for (int k = 50; k <= 1200; ++k)
        {
            const float d = k * 0.01f;
            const float e = energyAt(mode, rule, d);
            if (e < best) { best = e; at = d; }
        }
        return at;
    }
};

// AMBER's 12-6 with no radiusrule given: its radii must ADD, because amber*.ff
// stores R* = rmin/2.
TEST_F(TestRadiusRule, amber_defaults_to_the_sum_of_the_radii)
{
    EXPECT_NEAR(minimumOf("lennard-jones-12-6Amber", "default"), SUM, 0.02f);
}

// Levitt's 8-6 with no radiusrule given: its radii must be MEANED, because
// CAonlyLewitt.ff already stores the contact distance. Summing them would put
// the minimum at 5 A instead of 2 A.
TEST_F(TestRadiusRule, lewitt_defaults_to_the_geometric_mean)
{
    EXPECT_NEAR(minimumOf("lennard-jones-8-6Lewitt", "default"), GEOMETRIC, 0.02f);
}

// And the default is overridable in both directions, which is what lets a .ff
// be paired with a law it did not come from.
TEST_F(TestRadiusRule, the_msp_overrides_each_law_default)
{
    EXPECT_NEAR(minimumOf("lennard-jones-12-6Amber", "geometric-mean"), GEOMETRIC, 0.02f);
    EXPECT_NEAR(minimumOf("lennard-jones-12-6Amber", "arithmetic-mean"), ARITHMETIC, 0.02f);
    EXPECT_NEAR(minimumOf("lennard-jones-8-6Lewitt", "sum"), SUM, 0.02f);
    EXPECT_NEAR(minimumOf("lennard-jones-8-6Lewitt", "arithmetic-mean"), ARITHMETIC, 0.02f);
}

// The linear law always summed its radii and must keep doing so: two beads
// overlap when they touch, and the overlap is zero beyond that.
TEST_F(TestRadiusRule, the_linear_law_ignores_the_rule_and_keeps_summing)
{
    EXPECT_FLOAT_EQ(energyAt("linear", "default", SUM + 0.1f), 0.0f);
    EXPECT_GT(energyAt("linear", "default", SUM - 1.0f), 0.0f);
    // Even asked for something else: the law's equilibrium is a contact between
    // two spheres, and no combining rule applies to it.
    EXPECT_FLOAT_EQ(energyAt("linear", "geometric-mean", SUM + 0.1f), 0.0f);
    EXPECT_GT(energyAt("linear", "geometric-mean", SUM - 1.0f), 0.0f);
}

// Zacharias' law pairs a product on the radius with a product on epsilon, so
// its rule is not separable from the law and the setting must not reach it.
TEST_F(TestRadiusRule, zacharias_keeps_its_own_product_whatever_is_asked)
{
    const float own = minimumOf("lennard-jones-8-6Zacharias", "default");
    EXPECT_NEAR(own, sqrt(8.0f / 6.0f) * R1 * R2, 0.02f);
    EXPECT_FLOAT_EQ(minimumOf("lennard-jones-8-6Zacharias", "sum"), own);
    EXPECT_FLOAT_EQ(minimumOf("lennard-jones-8-6Zacharias", "geometric-mean"), own);
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
