#include <gtest/gtest.h>

#include <cmath>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"

using namespace biospring;

namespace
{

constexpr float C_MASS = 12.011f;
constexpr float H_MASS = 1.008f;

// Adds a particle at (x, y, 0) with the given mass, and returns its index.
unsigned addAtom(spn::SpringNetwork & network, float x, float y, float mass)
{
    spn::Particle p;
    p.setPosition(Vector3f(x, y, 0.0f));
    p.setMass(mass);
    p.setRadius(1.9f);
    p.setEpsilon(0.1f);
    p.setCharge(0.0f);
    network.addParticle(p);
    return network.getNumberOfParticles() - 1;
}

configuration::Configuration bareConfig(bool scale14)
{
    configuration::Configuration config = configuration::defaultConfiguration();
    config.sim.nbsteps = 1;
    config.sim.timestep = 0.0;
    config.spring.enable = true;
    config.scale14.enable = scale14;
    return config;
}

// A straight chain of four carbons, meshed the way --rigidbody meshes one: a
// spring on each bond AND on each valence angle, which is exactly the spring
// set the 1-4 derivation has to read. 1.5 A bonds at 120 degrees put the 1-3
// pairs at 2.60 A and the 1-4 pair at 3.90.
void buildChainOfFour(spn::SpringNetwork & network)
{
    const unsigned a = addAtom(network, 0.0f, 0.0f, C_MASS);
    const unsigned b = addAtom(network, 1.5f, 0.0f, C_MASS);
    const unsigned c = addAtom(network, 2.25f, 1.299f, C_MASS);
    const unsigned d = addAtom(network, 3.75f, 1.299f, C_MASS);
    network.addSpring(a, b, -1.0f, 100.0f);
    network.addSpring(b, c, -1.0f, 100.0f);
    network.addSpring(c, d, -1.0f, 100.0f);
    network.addSpring(a, c, -1.0f, 100.0f); // the a-b-c angle
    network.addSpring(b, d, -1.0f, 100.0f); // the b-c-d angle
}

} // namespace

// The pair three bonds apart is the one that gets the factor, and only it.
TEST(Scale14, FindsTheOnePairThreeBondsApart)
{
    spn::SpringNetwork network;
    buildChainOfFour(network);
    configuration::Configuration config = bareConfig(true);
    network.setup(config);

    EXPECT_TRUE(network.getParticle(0).isInOneFourNeighbors(3));
    EXPECT_TRUE(network.getParticle(3).isInOneFourNeighbors(0));
    // 1-2 and 1-3 are excluded by their spring, never scaled.
    EXPECT_FALSE(network.getParticle(0).isInOneFourNeighbors(1));
    EXPECT_FALSE(network.getParticle(0).isInOneFourNeighbors(2));
    EXPECT_FALSE(network.getParticle(1).isInOneFourNeighbors(3));
    EXPECT_EQ(network.getParticle(0).getNumberOfOneFourPairs(), 1u);
    EXPECT_EQ(network.getParticle(1).getNumberOfOneFourPairs(), 0u);
}

// Off unless asked for: enabling it changes the energy of every network
// already built, so a run that does not set it must see no 1-4 pair at all.
TEST(Scale14, RecordsNothingWhenDisabled)
{
    spn::SpringNetwork network;
    buildChainOfFour(network);
    configuration::Configuration config = bareConfig(false);
    network.setup(config);

    for (unsigned i = 0; i < network.getNumberOfParticles(); i++)
        EXPECT_EQ(network.getParticle(i).getNumberOfOneFourPairs(), 0u) << "particle " << i;
}

// TWO HYDROGENS ARE NEVER BONDED TO EACH OTHER. The mesh springs the H-C-H
// angle, and in a methyl that pair sits at 1.78 A -- SHORTER than the 1.5 A
// C-C bond plus its tolerance. A perception that went by distance alone would
// call it a bond, and then every atom on the far side of each hydrogen would
// become a spurious 1-4 partner.
TEST(Scale14, NeverPerceivesABondBetweenTwoHydrogens)
{
    spn::SpringNetwork network;
    const unsigned c1 = addAtom(network, 0.0f, 0.0f, C_MASS);
    const unsigned h1 = addAtom(network, 0.63f, 0.89f, H_MASS);
    const unsigned h2 = addAtom(network, 0.63f, -0.89f, H_MASS);
    const unsigned c2 = addAtom(network, -1.5f, 0.0f, C_MASS);
    network.addSpring(c1, h1, -1.0f, 100.0f);
    network.addSpring(c1, h2, -1.0f, 100.0f);
    network.addSpring(c1, c2, -1.0f, 100.0f);
    network.addSpring(h1, h2, -1.0f, 100.0f); // the H-C-H angle, 1.78 A
    network.addSpring(h1, c2, -1.0f, 100.0f); // the H-C-C angle
    network.addSpring(h2, c2, -1.0f, 100.0f);
    // The precondition this test rests on: the H...H 1-3 distance falls inside
    // the band a carbon-carbon bond occupies (0.76 + 0.76 + 0.45 = 1.97 A), so
    // a rule that went by distance alone would take it for a bond.
    const float hh = (network.getParticle(h1).getPosition() - network.getParticle(h2).getPosition()).norm();
    ASSERT_GT(hh, 1.5f) << "precondition: the H...H pair must be longer than a C-C bond";
    ASSERT_LT(hh, 1.97f) << "precondition: the H...H pair must be inside the C-C bond threshold";

    configuration::Configuration config = bareConfig(true);
    network.setup(config);

    // Were H1-H2 taken for a bond, C2 would be three bonds from H2 through it.
    for (unsigned i = 0; i < network.getNumberOfParticles(); i++)
        EXPECT_EQ(network.getParticle(i).getNumberOfOneFourPairs(), 0u)
            << "particle " << i << " gained a 1-4 partner through an H...H pair";
}

// A pair that is BOTH 1-4 and 1-3 -- which every ring of five produces -- is
// excluded by AMBER, not scaled. Going round a cyclopentane, atom 0 reaches
// atom 2 in two bonds one way and three the other.
TEST(Scale14, ARingPairThatIsAlsoOneThreeIsNotScaled)
{
    spn::SpringNetwork network;
    const unsigned n = 5;
    for (unsigned i = 0; i < n; i++)
    {
        const float a = 2.0f * static_cast<float>(M_PI) * static_cast<float>(i) / static_cast<float>(n);
        addAtom(network, 1.3f * std::cos(a), 1.3f * std::sin(a), C_MASS);
    }
    for (unsigned i = 0; i < n; i++)
        network.addSpring(i, (i + 1) % n, -1.0f, 100.0f); // the five bonds
    for (unsigned i = 0; i < n; i++)
        network.addSpring(i, (i + 2) % n, -1.0f, 100.0f); // the five angles

    configuration::Configuration config = bareConfig(true);
    network.setup(config);

    for (unsigned i = 0; i < n; i++)
        EXPECT_EQ(network.getParticle(i).getNumberOfOneFourPairs(), 0u)
            << "ring atom " << i << ": every pair of a 5-ring is 1-2 or 1-3, so none may be scaled";
}

// The factor reaches the energy, and it is exactly the factor asked for.
TEST(Scale14, StericEnergyOfTheOneFourPairIsScaled)
{
    const float factor = 0.5f;
    float energies[2] = {0.0f, 0.0f};
    for (int pass = 0; pass < 2; pass++)
    {
        spn::SpringNetwork network;
        buildChainOfFour(network);
        configuration::Configuration config = bareConfig(pass == 1);
        config.scale14.steric = factor;
        config.steric.enable = true;
        config.steric.mode = "lennard-jones-12-6Amber";
        config.steric.cutoff = 8.0;
        network.setup(config);
        network.run();
        energies[pass] = network.getStericEnergy();
    }
    // The chain's ONLY non-excluded pair is the 1-4 one, so the whole steric
    // energy is that pair's and the ratio is the factor itself.
    ASSERT_NE(energies[0], 0.0f) << "no steric energy at all: the test measures nothing";
    EXPECT_NEAR(energies[1] / energies[0], factor, 1.0e-5f)
        << "unscaled " << energies[0] << ", scaled " << energies[1];
}
