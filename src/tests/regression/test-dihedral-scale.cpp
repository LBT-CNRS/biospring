#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"

// dihedralphi.scale and its seven siblings weigh the torsion term. They existed
// in the parser long before anything read them -- EnergySetting carries a scale
// and every dihedral family is one -- so setting one parsed, printed in the
// configuration dump, and changed nothing. This is what makes that not happen
// again.
//
// The scale matters because there is no other way to weigh torsions against the
// rest of the model: spring.scale reaches them, but the rigid-body mesh that
// holds every bond and angle reads the same number, so turning the torsions
// down turns the bonds down with them.
struct DihedralScale : public ::testing::Test
{
    biospring::configuration::Configuration config;

    // A table that is FLAT: every bin the same energy and the same torque. The
    // point is the scaling, not the interpolation, and a flat table makes the
    // expected value exact whatever angle the four atoms happen to make.
    static constexpr float ENERGY = 8.0f;
    static constexpr float TORQUE = 3.0f;

    void build(biospring::spn::SpringNetwork & spn, unsigned family)
    {
        // Four atoms, deliberately not coplanar, so the dihedral is defined and
        // the gradient is non-zero on all four.
        const float xyz[4][3] = {{0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 0.0f},
                                 {1.5f, 0.0f, 0.0f}, {1.5f, 0.7f, 0.7f}};
        for (int i = 0; i < 4; ++i)
        {
            biospring::spn::Particle p;
            p.setPosition(Vector3f(xyz[i][0], xyz[i][1], xyz[i][2]));
            p.setMass(12.0f);
            spn.addParticle(p);
        }
        std::vector<std::array<std::vector<float>, 2>> tables(1);
        tables[0][0].assign(64 + 1, ENERGY);
        tables[0][1].assign(64 + 1, TORQUE);
        spn.setTorsionTables(tables);
        spn.addTorsion(family, 0, 1, 2, 3, 0);
        spn.setup(config);
    }

    float energyAtScale(unsigned family, double scale)
    {
        config = biospring::configuration::defaultConfiguration();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;
        config.dihedralphi.enable = true;
        config.dihedralpsi.enable = true;
        config.dihedralphi.scale = scale;
        config.dihedralpsi.scale = scale;
        biospring::spn::SpringNetwork spn;
        build(spn, family);
        spn.computeDihedralForces();
        return spn.getDihedralEnergy();
    }
};

TEST_F(DihedralScale, weighs_the_energy_in_proportion)
{
    const float full = energyAtScale(0, 1.0);
    ASSERT_NEAR(full, ENERGY, 1e-4f) << "a flat table should give exactly its own value at scale 1";

    EXPECT_NEAR(energyAtScale(0, 0.5), 0.5f * ENERGY, 1e-4f);
    EXPECT_NEAR(energyAtScale(0, 2.0), 2.0f * ENERGY, 1e-4f);
    EXPECT_NEAR(energyAtScale(0, 0.0), 0.0f, 1e-4f)
        << "scale 0 has to silence the family completely, not merely weaken it";
}

TEST_F(DihedralScale, weighs_the_force_in_proportion)
{
    // The force, not just the number the energy reports: a scale applied to the
    // energy and not to the torque would pass the test above and change nothing
    // about the trajectory, which is the whole point of having the parameter.
    const auto forceNorm = [this](double scale) {
        config = biospring::configuration::defaultConfiguration();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;
        config.dihedralphi.enable = true;
        config.dihedralphi.scale = scale;
        biospring::spn::SpringNetwork spn;
        build(spn, 0);
        spn.computeDihedralForces();
        float worst = 0.0f;
        for (unsigned i = 0; i < 4; ++i)
            worst = std::max(worst, spn.getParticle(i).getForce().norm());
        return worst;
    };
    const float full = forceNorm(1.0);
    ASSERT_GT(full, 0.0f) << "the flat table has a torque, so there must be a force";
    EXPECT_NEAR(forceNorm(0.5), 0.5f * full, 1e-3f * full);
    EXPECT_NEAR(forceNorm(0.0), 0.0f, 1e-5f * full);
}

TEST_F(DihedralScale, weighs_each_family_on_its_own)
{
    // Eight families, eight weights. phi/psi decide a backbone's secondary
    // structure, chi a side chain's rotamer, the planarity impropers hold a
    // ring flat -- scaling them together would make the parameter useless for
    // the one question it exists to answer.
    config = biospring::configuration::defaultConfiguration();
    config.sim.nbsteps = 1;
    config.sim.timestep = 0.01;
    config.dihedralphi.enable = true;
    config.dihedralpsi.enable = true;
    config.dihedralphi.scale = 0.25;
    config.dihedralpsi.scale = 1.0;

    biospring::spn::SpringNetwork phi, psi;
    build(phi, 0);   // family 0 is phi
    config = biospring::configuration::defaultConfiguration();
    config.sim.nbsteps = 1;
    config.sim.timestep = 0.01;
    config.dihedralphi.enable = true;
    config.dihedralpsi.enable = true;
    config.dihedralphi.scale = 0.25;
    config.dihedralpsi.scale = 1.0;
    build(psi, 1);   // family 1 is psi

    phi.computeDihedralForces();
    psi.computeDihedralForces();
    EXPECT_NEAR(phi.getDihedralEnergy(), 0.25f * ENERGY, 1e-4f)
        << "phi should have taken its own weight";
    EXPECT_NEAR(psi.getDihedralEnergy(), 1.0f * ENERGY, 1e-4f)
        << "psi should have been left alone by phi's weight";
}

int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
