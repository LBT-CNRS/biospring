#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include "IO/io.h"
#include "reduce/Reducer.h"
#include "measure.hpp"
#include "topology.hpp"

static const std::string PATH_FORCEFIELD = "data/model.ff";
static const std::string PATH_REDUCE_RULES = "data/model.grp";
static const std::string PATH_TOPOLOGY = "data/model.pdb";

struct TestReducer : public ::testing::Test
{
    std::string path_topology;
    std::string path_forcefield;
    std::string path_reduce_rules;
    biospring::topology::Topology topology;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        path_topology = PATH_TOPOLOGY;
        path_forcefield = PATH_FORCEFIELD;
        path_reduce_rules = PATH_REDUCE_RULES;

        ASSERT_TRUE(std::filesystem::exists(path_topology));
        ASSERT_TRUE(std::filesystem::exists(path_forcefield));
        ASSERT_TRUE(std::filesystem::exists(path_reduce_rules));

        topology = biospring::io::readTopology(path_topology);
        ASSERT_EQ(topology.number_of_particles(), 2814);
    }
};

// ===========================================================================
// Reduce tests
// ===========================================================================

TEST_F(TestReducer, reduce)
{
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(path_reduce_rules);
    reducer.reduce();

    // model.grp's SER rule used to require a nonexistent "CO" atom, so every SER
    // grain was always incomplete and dropped by default; now that it correctly
    // requires "OG", the 9 SER residues form complete grains and are kept too.
    EXPECT_EQ(reducer.target_topology().number_of_particles(), 183);
}

TEST_F(TestReducer, reduce_ignore_missing)
{
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(path_reduce_rules);
    reducer.set_ignore_missing_particles(true);
    reducer.reduce();

    EXPECT_EQ(reducer.target_topology().number_of_particles(), 183);
}

// ===========================================================================
// Basic tests
// ===========================================================================

TEST_F(TestReducer, initialize_forcefield)
{
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    EXPECT_EQ(reducer.forcefield().getNumberOfProperties(), 20);
}

TEST_F(TestReducer, initialize_rules)
{
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_rules(path_reduce_rules);
    EXPECT_EQ(reducer.rules().size(), 20);
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// A declared bond -- a PDB CONECT record -- must survive reduction AND the
// copy pdb2spn makes immediately after it.
//
// This is the test that was missing, and its absence let a real bug through.
// A Spring holds a Particle&, so one created inside the reducer points into
// storage the copy leaves behind. The failure is SILENT: the spring keeps its
// equilibrium length and simply joins two different atoms. Measured on 074's
// duplex, 48 base-pair springs with a correct r0 of 2.7-3.3 A ended up on
// atoms 5 to 12 A apart, which at mesh stiffness threw the structure to 1e8 A.
//
// The invariant that catches it: a carried spring's equilibrium is the length
// it was declared with, so it must still equal the distance between the two
// particles it joins. Every other check -- spring count, equilibrium value,
// even the endpoints' residue ids -- passes while the bond is wrong.
TEST_F(TestReducer, declared_bonds_survive_the_copy_out_of_the_reducer)
{
    // The fixture's PDB declares no bond, so make one: two backbone nitrogens
    // far enough apart in sequence to land in different grains.
    size_t first = topology.number_of_particles(), second = topology.number_of_particles();
    int first_residue = -1;
    for (size_t i = 0; i < topology.number_of_particles(); ++i)
    {
        const auto & p = topology.get_particle(i).properties();
        if (p.name() != "N")
            continue;
        if (first == topology.number_of_particles())
        {
            first = i;
            first_residue = p.residue_id();
        }
        else if (p.residue_id() > first_residue + 4)
        {
            second = i;
            break;
        }
    }
    ASSERT_LT(first, topology.number_of_particles());
    ASSERT_LT(second, topology.number_of_particles());
    topology.add_spring(topology.get_particle(first), topology.get_particle(second));
    ASSERT_EQ(topology.number_of_springs(), 1u);

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(path_reduce_rules);
    reducer.reduce();

    // Exactly what pdb2spn does, copy included.
    biospring::topology::Topology reduced = reducer.target_topology();
    reducer.carry_declared_bonds(reduced);

    ASSERT_EQ(reduced.number_of_springs(), 1u) << "the declared bond was lost by reduction";
    const auto & spring = reduced.get_spring(0);
    EXPECT_NEAR(spring.equilibrium(), biospring::measure::distance(spring.first(), spring.second()), 1e-3)
        << "the spring kept its length but moved to different atoms";
}
