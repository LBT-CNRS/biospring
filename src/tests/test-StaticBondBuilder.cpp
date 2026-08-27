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

// Stands in for a CONECT record: pdb2spn turns each one into a spring carrying
// the global -stiffness, and that is what these functions retune.
void declare(topology::Topology & top, size_t a, size_t b, double stiffness = 8000.0)
{
    top.add_spring(top.get_particle(a), top.get_particle(b), -1.0, stiffness);
}

staticbond::DonorAcceptorTable backboneTable()
{
    staticbond::DonorAcceptorTable t;
    t[{"ALA", "N"}] = {1, 0};
    t[{"ALA", "O"}] = {0, 2};
    t[{"ALA", "CB"}] = {0, 0};
    return t;
}

} // namespace

// The calibration: a hydrogen-bond spring must carry the curvature of a Morse
// well at its minimum. If someone retunes De or a without retuning the spring,
// this catches it -- nothing else would, because the energy at equilibrium is
// zero either way. Away from equilibrium a harmonic has no claim to be right
// and none is made; the constant is a stiffness AROUND the equilibrium length.
TEST(StaticBondBuilder, hydrogen_bond_stiffness_matches_the_morse_curvature)
{
    const float de = 16.7f; // kJ.mol-1
    const float a = 1.34f;  // A-1
    EXPECT_NEAR(staticbond::HYDROGEN_BOND_STIFFNESS, 2.0f * de * a * a, 0.05f);
}

// AMBER ff14SB's S-S force constant, in BioSpring's own 0.5.k.dr^2 convention.
TEST(StaticBondBuilder, disulfide_stiffness_is_the_amber_value)
{
    EXPECT_NEAR(staticbond::DISULFIDE_STIFFNESS, 1389.1f, 0.1f);
}

TEST(StaticBondBuilder, retunes_a_declared_donor_acceptor_bond)
{
    topology::Topology top;
    const size_t n = place(top, "N", "ALA", "A", 1, 0.0f, 0.0f, 0.0f);
    const size_t o = place(top, "O", "ALA", "A", 10, 2.9f, 0.0f, 0.0f);
    declare(top, n, o);

    EXPECT_EQ(staticbond::retuneHydrogenBondSprings(top, backboneTable(),
                                                    staticbond::HYDROGEN_BOND_STIFFNESS),
              1u);
    ASSERT_EQ(top.number_of_springs(), 1u);
    EXPECT_NEAR(top.get_spring(0).stiffness(), staticbond::HYDROGEN_BOND_STIFFNESS, 1e-3);
    // Equilibrium is the declared bond's own length: retuning a force constant
    // must not move the structure it was given.
    EXPECT_NEAR(top.get_spring(0).equilibrium(), 2.9, 1e-4);
}

// The orientation of a CONECT record carries no meaning -- it does not say
// which end donates -- so both must be recognised.
TEST(StaticBondBuilder, recognises_a_declared_bond_in_either_orientation)
{
    topology::Topology top;
    const size_t o = place(top, "O", "ALA", "A", 1, 0.0f, 0.0f, 0.0f);
    const size_t n = place(top, "N", "ALA", "A", 10, 2.9f, 0.0f, 0.0f);
    declare(top, o, n);

    EXPECT_EQ(staticbond::retuneHydrogenBondSprings(top, backboneTable(), 60.0f), 1u);
}

// THE POINT OF THE WHOLE DESIGN: nothing is searched for. Two atoms at a
// perfect hydrogen-bond distance that the structure does not declare bonded
// stay unbonded.
TEST(StaticBondBuilder, invents_no_bond_the_structure_did_not_declare)
{
    topology::Topology top;
    place(top, "N", "ALA", "A", 1, 0.0f, 0.0f, 0.0f);
    place(top, "O", "ALA", "A", 10, 2.9f, 0.0f, 0.0f);

    EXPECT_EQ(staticbond::retuneHydrogenBondSprings(top, backboneTable(), 60.0f), 0u);
    EXPECT_EQ(top.number_of_springs(), 0u);
}

// A CONECT record also carries ordinary covalent bonds. Only a donor/acceptor
// pair is a hydrogen bond; everything else keeps the stiffness it had.
TEST(StaticBondBuilder, leaves_a_declared_bond_that_is_not_a_hydrogen_bond_alone)
{
    topology::Topology top;
    const size_t cb = place(top, "CB", "ALA", "A", 1, 0.0f, 0.0f, 0.0f);
    const size_t o = place(top, "O", "ALA", "A", 10, 2.9f, 0.0f, 0.0f);
    declare(top, cb, o);

    EXPECT_EQ(staticbond::retuneHydrogenBondSprings(top, backboneTable(), 60.0f), 0u);
    EXPECT_NEAR(top.get_spring(0).stiffness(), 8000.0, 1e-3);
}

TEST(StaticBondBuilder, retunes_a_declared_disulfide)
{
    topology::Topology top;
    const size_t a = place(top, "SG", "CYS", "A", 22, 0.0f, 0.0f, 0.0f);
    const size_t b = place(top, "SG", "CYS", "A", 157, 2.04f, 0.0f, 0.0f);
    declare(top, a, b);

    EXPECT_EQ(staticbond::retuneDisulfideSprings(top, staticbond::DISULFIDE_STIFFNESS), 1u);
    EXPECT_NEAR(top.get_spring(0).stiffness(), staticbond::DISULFIDE_STIFFNESS, 1e-2);
    EXPECT_NEAR(top.get_spring(0).equilibrium(), 2.04, 1e-4);
}

// amber.grp renames CYS's sulfur to CSG. A test on the raw PDB name alone
// would pass while the option silently found nothing on every reduced system.
TEST(StaticBondBuilder, finds_the_sulfur_under_its_reduced_name)
{
    topology::Topology top;
    const size_t a = place(top, "CSG", "CYS", "A", 22, 0.0f, 0.0f, 0.0f);
    const size_t b = place(top, "CSG", "CYS", "A", 157, 2.04f, 0.0f, 0.0f);
    declare(top, a, b);

    EXPECT_EQ(staticbond::retuneDisulfideSprings(top, staticbond::DISULFIDE_STIFFNESS), 1u);
}

// Two cysteine sulfurs the structure does not declare bonded are two free
// thiols, however close they sit.
TEST(StaticBondBuilder, invents_no_disulfide_the_structure_did_not_declare)
{
    topology::Topology top;
    place(top, "SG", "CYS", "A", 22, 0.0f, 0.0f, 0.0f);
    place(top, "SG", "CYS", "A", 157, 2.04f, 0.0f, 0.0f);

    EXPECT_EQ(staticbond::retuneDisulfideSprings(top, staticbond::DISULFIDE_STIFFNESS), 0u);
    EXPECT_EQ(top.number_of_springs(), 0u);
}

// Both options read the same declared bonds, so a structure that declares one
// of each must come out with one spring of each stiffness -- neither pass may
// claim the other's bond.
TEST(StaticBondBuilder, the_two_passes_do_not_claim_each_others_bonds)
{
    topology::Topology top;
    const size_t n = place(top, "N", "ALA", "A", 1, 0.0f, 0.0f, 0.0f);
    const size_t o = place(top, "O", "ALA", "A", 10, 2.9f, 0.0f, 0.0f);
    const size_t s1 = place(top, "SG", "CYS", "A", 22, 10.0f, 0.0f, 0.0f);
    const size_t s2 = place(top, "SG", "CYS", "A", 157, 12.04f, 0.0f, 0.0f);
    declare(top, n, o);
    declare(top, s1, s2);

    EXPECT_EQ(staticbond::retuneHydrogenBondSprings(top, backboneTable(),
                                                    staticbond::HYDROGEN_BOND_STIFFNESS),
              1u);
    EXPECT_EQ(staticbond::retuneDisulfideSprings(top, staticbond::DISULFIDE_STIFFNESS), 1u);
    EXPECT_NEAR(top.get_spring(0).stiffness(), staticbond::HYDROGEN_BOND_STIFFNESS, 1e-2);
    EXPECT_NEAR(top.get_spring(1).stiffness(), staticbond::DISULFIDE_STIFFNESS, 1e-2);
}
