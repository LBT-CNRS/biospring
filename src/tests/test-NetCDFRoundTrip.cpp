
#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include "IO/NetCDFReader.h"
#include "IO/NetCDFWriter.h"
#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"

using namespace biospring;

namespace
{

// Three separate per-particle floats that have historically been confused with
// one another. They are given distinct values so that writing one and reading
// back another is a visible failure rather than a coincidence:
//   hydrophobicity  -> the pairwise hydrophobic force term
//   transfer energy -> IMPALA, stored as "hydrophobicityscale"
//   accessible surface -> stored as "surfaceaccessibility"
constexpr float HYDROPHOBICITY = 0.25f;
constexpr float TRANSFER_ENERGY = 7.5f;
constexpr float ACCESSIBLE_SURFACE = 123.0f;

struct TestNetCDFRoundTrip : public ::testing::Test
{
    configuration::Configuration config;
    spn::SpringNetwork spn;
    std::string path;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;

        path = (std::filesystem::temp_directory_path() / "roundtrip.nc").string();

        for (int i = 0; i < 2; ++i)
        {
            spn::Particle p;
            p.setPosition(Vector3f(static_cast<float>(i), 0.0, 0.0));
            p.setMass(12.0);
            p.setRadius(1.5);
            p.setHydrophobicity(HYDROPHOBICITY);
            p.setTransferEnergyByAccessibleSurface(TRANSFER_ENERGY);
            p.setSolventAccessibilitySurface(ACCESSIBLE_SURFACE);
            spn.addParticle(p);
        }
        spn.setup(config);

        NetCDFWriter writer(path, &spn);
        writer.writeBinary();
    }
};

} // namespace

// The regression: `hydrophobicity` must survive a write/read cycle. It was lost
// when a refactor reused its slot in the particle buffer for the accessible
// surface, which left the pairwise hydrophobic force unable to be fed from a
// .nc file at all.
TEST_F(TestNetCDFRoundTrip, HydrophobicitySurvivesTheRoundTrip)
{
    NetCDFReader reader(path);
    reader.read();

    ASSERT_EQ(reader.getTopology().number_of_particles(), 2u);
    for (size_t i = 0; i < 2; ++i)
        EXPECT_FLOAT_EQ(reader.getTopology().get_particle(i).properties().hydrophobicity(), HYDROPHOBICITY);
}

// ... and it must not be confused with either of its neighbours. All three are
// written, so a test that only checked two of them could pass while the other
// two were crossed.
TEST_F(TestNetCDFRoundTrip, ThePerParticleFloatsAreNotSwapped)
{
    NetCDFReader reader(path);
    reader.read();

    ASSERT_EQ(reader.getTopology().number_of_particles(), 2u);
    const auto & properties = reader.getTopology().get_particle(0).properties();

    EXPECT_FLOAT_EQ(properties.hydrophobicity(), HYDROPHOBICITY);
    EXPECT_FLOAT_EQ(properties.imp().transfert_energy_by_accessible_surface(), TRANSFER_ENERGY);
    EXPECT_FLOAT_EQ(properties.imp().solvent_accessible_surface(), ACCESSIBLE_SURFACE);
}

// The axis substituent lists have to survive the .nc, because that is the only
// way they reach biospring at all: they come from DIHEDRALAXIS records read by
// pdb2spn, and biospring never sees a .bi.ff. Lose them here and
// dihedral.distributetorque silently has nothing to act on.
//
// Written as a flat atom list sliced by per-axis counts, which is where this
// can go wrong quietly: the two sides must come back with the right lengths,
// in the right order, on the right axis. So the sides are given DIFFERENT
// lengths -- a symmetric case would pass even with the two swapped.
TEST(NetCDFRoundTrip, DihedralAxisSubstituentsSurviveTheRoundTrip)
{
    const std::string path = (std::filesystem::temp_directory_path() / "roundtrip-axes.nc").string();

    spn::SpringNetwork spn;
    for (int i = 0; i < 5; ++i)
    {
        spn::Particle p;
        p.setPosition(Vector3f(static_cast<float>(i), 0.0, 0.0));
        p.setMass(12.0);
        spn.addParticle(p);
    }
    // A ghost on the (0, 1) axis: without one there is no axis to attach to.
    spn.addGhostParticle(static_cast<unsigned>(spn::GhostPlacement::AxisRotation), 0, 1, 2, 1.0f, 75.0f, 0.0f);
    spn.setAxisSubstituents(0, 1, {2, 3}, {4});
    ASSERT_EQ(spn.getNumberOfAxesWithSubstituents(), 1u);

    configuration::Configuration config;
    config.sim.nbsteps = 1;
    spn.setup(config);
    NetCDFWriter(path, &spn).writeBinary();

    NetCDFReader reader(path);
    reader.read();
    spn::SpringNetwork reloaded;
    reader.getTopology().to_spring_network(reloaded);

    std::vector<std::array<unsigned, 2>> anchors, counts;
    std::vector<unsigned> atoms;
    reloaded.getAxisSubstituents(anchors, counts, atoms);

    ASSERT_EQ(anchors.size(), 1u);
    EXPECT_EQ(anchors[0][0], 0u);
    EXPECT_EQ(anchors[0][1], 1u);
    EXPECT_EQ(counts[0][0], 2u);                      // B side, two atoms
    EXPECT_EQ(counts[0][1], 1u);                      // C side, one -- deliberately different
    ASSERT_EQ(atoms.size(), 3u);
    EXPECT_EQ(atoms[0], 2u);
    EXPECT_EQ(atoms[1], 3u);
    EXPECT_EQ(atoms[2], 4u);
}
