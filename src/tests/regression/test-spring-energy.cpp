
#include <gtest/gtest.h>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "forcefield/constants.hpp"
#include "topology.hpp"

// Reference function that calculates the spring energy for a given spring.
float spring_energy(const biospring::spn::Spring & spring)
{
    const biospring::spn::Particle & lhs = spring.getParticle1();
    const biospring::spn::Particle & rhs = spring.getParticle2();

    float stiffness = biospring::spn::Spring::DEFAULT_STIFFNESS; // as designed in TestSpringEnergy.
    float equilibrium = 0.0; // initial distance between p1 and p2 (as designed in TestSpringEnergy).
    float distance = biospring::spn::Particle::distance(lhs, rhs);
    float r2 = (distance - equilibrium) * (distance - equilibrium);

    return 0.5 * stiffness * r2;
}

struct TestSpringEnergy : public ::testing::Test
{
    biospring::configuration::Configuration config;
    biospring::spn::SpringNetwork spn;
    biospring::spn::Particle p1, p2;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;
        config.spring.enable = true;
        config.spring.cutoff = 16.0;

        SetUpSpn();

        ASSERT_EQ(spn.getNumberOfParticles(), 2);
        ASSERT_EQ(spn.getNumberOfSprings(), 1);

        const biospring::spn::Spring & spring = spn.getSpring(0);
        ASSERT_FLOAT_EQ(spring.getStiffness(), biospring::spn::Spring::DEFAULT_STIFFNESS);
        ASSERT_FLOAT_EQ(spring.getEquilibrium(), 0.0); // distance between p1 and p2.
    }

    // Sets up the SpringNetwork with two particles.
    void SetUpSpn()
    {
        biospring::topology::Topology top;
        biospring::topology::Particle p1, p2;

        p1.properties().set_position(Vector3f(0.0, 0.0, 0.0));
        p2.properties().set_position(Vector3f(0.0, 0.0, 0.0));

        // Parameters corresponding to the CC in amber.ff.
        p1.properties().set_charge(0.5973);
        p2.properties().set_charge(0.5973);

        p1.properties().set_radius(1.908);
        p2.properties().set_radius(1.908);

        p1.properties().set_epsilon(0.086);
        p2.properties().set_epsilon(0.086);

        p1.properties().set_mass(12.01);
        p2.properties().set_mass(12.01);

        top.add_particle(p1);
        top.add_particle(p2);

        top.add_spring(top.get_particle(0), top.get_particle(1));

        top.to_spring_network(spn);
        spn.setup(config);
    }
};

// ============================================================================

TEST_F(TestSpringEnergy, spring)
{
    for (float x = 0.01f; x < 5; x += 0.01)
    {
        const biospring::spn::Particle & lhs = spn.getParticle(0);
        biospring::spn::Particle & rhs = spn.getParticle(1);

        rhs.setPosition(Vector3f(x, 0.0, 0.0));
        spn.idleRun();
        spn.computeForces();

        float distance = lhs.distance(rhs);
        float actual = spn.getSpringEnergy();
        float expected = spring_energy(spn.getSpring(0));

        EXPECT_FLOAT_EQ(actual, expected);
    }
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
