#include <cmath>
#include <gtest/gtest.h>

#include "staticbond/StaticBondBuilder.h"
#include "topology/Topology.hpp"

using namespace biospring;

namespace
{

// Adds one particle and returns its INDEX, never a reference: the particle
// vector reallocates as it grows, so a reference taken before the next
// add_particle dangles. Fetch the reference at the point of use instead.
size_t place(topology::Topology & top, const std::string & name, const std::string & resname,
             const std::string & chain, int resid, float x, float y, float z)
{
    topology::ParticleProperties properties;
    properties.set_name(name);
    properties.set_residue_name(resname);
    properties.set_chain_name(chain);
    properties.set_residue_id(resid);
    properties.set_position({x, y, z});
    top.add_particle(topology::Particle(properties));
    return top.number_of_particles() - 1;
}

staticbond::DonorAcceptorTable backboneTable()
{
    return {{{"ALA", "N"}, {1, 0}}, {{"ALA", "O"}, {0, 2}}};
}

} // namespace

// The whole point of the calibration: a hydrogen-bond spring must have the
// curvature of the Morse well the dynamic model used, at its minimum. If
// someone retunes De or a without retuning the spring, this catches it --
// nothing else would, because the energy at equilibrium is zero either way.
TEST(StaticBondBuilder, hydrogen_bond_stiffness_matches_the_morse_curvature)
{
    const float de = 16.7f;   // kJ.mol-1
    const float a = 1.34f;    // A-1
    EXPECT_NEAR(staticbond::HYDROGEN_BOND_STIFFNESS, 2.0f * de * a * a, 0.05f);
}

// AMBER ff14SB's S-S force constant, in BioSpring's own 0.5.k.dr^2 convention.
TEST(StaticBondBuilder, disulfide_stiffness_is_the_amber_value)
{
    EXPECT_NEAR(staticbond::DISULFIDE_STIFFNESS, 1389.1f, 0.1f);
}

TEST(StaticBondBuilder, springs_a_donor_acceptor_pair_within_the_cutoff)
{
    topology::Topology top;
    place(top, "N", "ALA", "A", 1, 0.0f, 0.0f, 0.0f);
    place(top, "O", "ALA", "A", 10, 2.9f, 0.0f, 0.0f);

    const size_t added = staticbond::addHydrogenBondSprings(top, backboneTable(), 3.5f,
                                                            staticbond::HYDROGEN_BOND_STIFFNESS);
    ASSERT_EQ(added, 1u);
    ASSERT_EQ(top.number_of_springs(), 1u);

    // Equilibrium is the OBSERVED distance, so the input structure starts at
    // rest: turning the option on must not deform what it was given.
    const auto & spring = top.get_spring(0);
    EXPECT_NEAR(spring.equilibrium(), 2.9, 1e-4);
    EXPECT_NEAR(spring.stiffness(), staticbond::HYDROGEN_BOND_STIFFNESS, 1e-3);
}

TEST(StaticBondBuilder, ignores_a_pair_beyond_the_cutoff)
{
    topology::Topology top;
    place(top, "N", "ALA", "A", 1, 0.0f, 0.0f, 0.0f);
    place(top, "O", "ALA", "A", 10, 4.2f, 0.0f, 0.0f);

    EXPECT_EQ(staticbond::addHydrogenBondSprings(top, backboneTable(), 3.5f, 60.0f), 0u);
    EXPECT_EQ(top.number_of_springs(), 0u);
}

// An i/i+1 donor-acceptor contact is the backbone's own geometry. Springing it
// would quietly stiffen the chain in a way no one asked for.
TEST(StaticBondBuilder, ignores_neighbours_in_sequence)
{
    topology::Topology top;
    place(top, "N", "ALA", "A", 5, 0.0f, 0.0f, 0.0f);
    place(top, "O", "ALA", "A", 6, 2.9f, 0.0f, 0.0f);

    EXPECT_EQ(staticbond::addHydrogenBondSprings(top, backboneTable(), 3.5f, 60.0f), 0u);
}

// The same residue id in ANOTHER chain is not a sequence neighbour.
TEST(StaticBondBuilder, spans_two_chains_at_the_same_residue_id)
{
    topology::Topology top;
    place(top, "N", "ALA", "A", 5, 0.0f, 0.0f, 0.0f);
    place(top, "O", "ALA", "B", 5, 2.9f, 0.0f, 0.0f);

    EXPECT_EQ(staticbond::addHydrogenBondSprings(top, backboneTable(), 3.5f, 60.0f), 1u);
}

// Additivity is the contract: -static-hbond must compose with a network built
// by -cutoff or -rigidbody, and must never overwrite one of its springs.
TEST(StaticBondBuilder, leaves_an_existing_spring_alone)
{
    topology::Topology top;
    const size_t donor = place(top, "N", "ALA", "A", 1, 0.0f, 0.0f, 0.0f);
    const size_t acceptor = place(top, "O", "ALA", "A", 10, 2.9f, 0.0f, 0.0f);
    top.add_spring(top.get_particle(donor), top.get_particle(acceptor), -1.0, 8000.0);

    EXPECT_EQ(staticbond::addHydrogenBondSprings(top, backboneTable(), 3.5f, 60.0f), 0u);
    ASSERT_EQ(top.number_of_springs(), 1u);
    EXPECT_NEAR(top.get_spring(0).stiffness(), 8000.0, 1e-3);
}

TEST(StaticBondBuilder, springs_a_disulfide_between_two_cysteine_sulfurs)
{
    topology::Topology top;
    place(top, "SG", "CYS", "A", 22, 0.0f, 0.0f, 0.0f);
    place(top, "SG", "CYS", "A", 157, 2.04f, 0.0f, 0.0f);

    EXPECT_EQ(staticbond::addDisulfideSprings(top, staticbond::DISULFIDE_CUTOFF, staticbond::DISULFIDE_STIFFNESS), 1u);
    ASSERT_EQ(top.number_of_springs(), 1u);
    EXPECT_NEAR(top.get_spring(0).stiffness(), staticbond::DISULFIDE_STIFFNESS, 1e-2);
}

// amber.grp renames CYS's sulfur to CSG. A test on the raw PDB name alone
// would pass while the option silently found nothing on every reduced system.
TEST(StaticBondBuilder, finds_the_sulfur_under_its_reduced_name)
{
    topology::Topology top;
    place(top, "CSG", "CYS", "A", 22, 0.0f, 0.0f, 0.0f);
    place(top, "CSG", "CYS", "A", 157, 2.04f, 0.0f, 0.0f);

    EXPECT_EQ(staticbond::addDisulfideSprings(top, staticbond::DISULFIDE_CUTOFF, staticbond::DISULFIDE_STIFFNESS), 1u);
}

// A free thiol pair is not a bridge.
TEST(StaticBondBuilder, ignores_sulfurs_too_far_apart_to_be_bonded)
{
    topology::Topology top;
    place(top, "SG", "CYS", "A", 22, 0.0f, 0.0f, 0.0f);
    place(top, "SG", "CYS", "A", 157, 6.0f, 0.0f, 0.0f);

    EXPECT_EQ(staticbond::addDisulfideSprings(top, staticbond::DISULFIDE_CUTOFF, staticbond::DISULFIDE_STIFFNESS), 0u);
    EXPECT_EQ(top.number_of_springs(), 0u);
}

// The opposite choice from a hydrogen bond, and the one that matters in
// practice: a PDB CONECT record turns the bridge into a spring carrying the
// global -stiffness (8000 in the rigid-body examples). That is the same bond,
// so it must end up with the bridge's force constant, not the mesh's.
TEST(StaticBondBuilder, recalibrates_a_disulfide_already_sprung_by_a_conect_record)
{
    topology::Topology top;
    const size_t first = place(top, "SG", "CYS", "A", 21, 0.0f, 0.0f, 0.0f);
    const size_t second = place(top, "SG", "CYS", "A", 174, 2.08f, 0.0f, 0.0f);
    top.add_spring(top.get_particle(first), top.get_particle(second), -1.0, 8000.0);

    EXPECT_EQ(staticbond::addDisulfideSprings(top, staticbond::DISULFIDE_CUTOFF, staticbond::DISULFIDE_STIFFNESS), 1u);
    ASSERT_EQ(top.number_of_springs(), 1u);
    EXPECT_NEAR(top.get_spring(0).stiffness(), staticbond::DISULFIDE_STIFFNESS, 1e-2);
}
