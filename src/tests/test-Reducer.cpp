#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <utility>

#include "IO/io.h"
#include "IO/ReduceRuleReader.h"
#include "reduce/Reducer.h"
#include "TempPath.h"
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

// A declared bond -- a PDB CONECT record -- must survive reduction AND the
// assignment pdb2spn makes immediately after it.
//
// This is the test that was missing, and its absence let a real bug through.
// A Spring holds a Particle&, not an index, so any spring that outlives a
// replacement of the particle list points at whatever landed at its old
// position. The failure is SILENT: the spring keeps its equilibrium length and
// simply joins two different atoms. Measured on 074's duplex, 48 base-pair
// springs with a correct r0 of 2.7-3.3 A ended up on atoms 5 to 12 A apart,
// which at mesh stiffness threw the structure to 1e8 A.
//
// The invariant that catches it: a carried spring's equilibrium is the length
// it was declared with, so it must still equal the distance between the two
// particles it joins. Every other check -- spring count, equilibrium value,
// even the endpoints' residue ids -- passes while the bond is wrong.
TEST_F(TestReducer, declared_bonds_survive_reduction_and_the_assignment)
{
    // Learn first which grains this model actually keeps: most residues here
    // lose one, and a bond whose end lands in a dropped grain is correctly
    // discarded -- which would make this test pass for the wrong reason.
    std::set<std::pair<int, std::string>> kept;
    {
        biospring::reduce::Reducer probe(topology);
        probe.initialize_forcefield(path_forcefield);
        probe.initialize_rules(path_reduce_rules);
        probe.reduce();
        const auto & reduced = probe.target_topology();
        for (size_t i = 0; i < reduced.number_of_particles(); ++i)
        {
            const auto & p = reduced.get_particle(i).properties();
            kept.insert({p.residue_id(), p.name()});
        }
    }
    ASSERT_FALSE(kept.empty());

    // Two atoms whose grains both survive, in residues far enough apart that
    // they cannot merge into one grain.
    biospring::reduce::ReduceRuleReader rules_reader(path_reduce_rules);
    rules_reader.read();
    const auto all_rules = rules_reader.rules();
    auto grain_of = [&](const biospring::topology::Particle & particle) -> std::string
    {
        const auto & p = particle.properties();
        const auto rules = all_rules.get_rules_for_residue(p.residue_name());
        for (size_t r = 0; r < rules.size(); ++r)
            if (rules[r].hasAtomNamed(p.name()))
                return rules[r].getName();
        return "";
    };

    size_t first = topology.number_of_particles(), second = topology.number_of_particles();
    int first_residue = -1;
    for (size_t i = 0; i < topology.number_of_particles(); ++i)
    {
        const auto & p = topology.get_particle(i).properties();
        const std::string grain = grain_of(topology.get_particle(i));
        if (grain.empty() || kept.find({p.residue_id(), grain}) == kept.end())
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

    // Exactly what pdb2spn does, assignment included -- and the assignment is
    // where this used to break. Topology::operator= kept its own springs
    // whenever the two particle counts happened to match, while
    // _copy_particles replaced, in place, the very list those springs hold
    // references into. Equal counts say nothing about equal ordering, and the
    // reducer emits grains in rule order.
    topology = reducer.target_topology();

    {
        std::set<biospring::topology::pid_t> uids;
        size_t dup = 0;
        for (size_t i = 0; i < topology.number_of_particles(); ++i)
            if (!uids.insert(topology.get_particle(i).unique_id()).second)
                ++dup;
        EXPECT_EQ(dup, 0u) << "the reduced topology reuses particle unique ids, which is what the .nc "
                              "writer resolves spring endpoints through";
    }
    ASSERT_EQ(topology.number_of_springs(), 1u) << "the declared bond was lost by reduction";
    const auto & spring = topology.get_spring(0);
    EXPECT_NEAR(spring.equilibrium(), biospring::measure::distance(spring.first(), spring.second()), 1e-3)
        << "the spring kept its length but moved to different atoms";
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// ===========================================================================
// Rules that match nothing
// ===========================================================================
//
// amber.grp declares the N-terminal ammonium hydrogens and the C-terminal OXT
// for all 20 residues, so most of those rules necessarily match nothing in any
// given model. Reported per residue, as it was, that came to 42% of the whole
// example suite's build output and 10334 of 076's 10434 lines, which buried the
// four warnings that example actually had. These two tests pin down the shape of
// the replacement: a rule matching nowhere is COUNTED, costs no grain, and a
// rule matching somewhere is not counted at all.

static std::string grp_plus(const std::string & path, const std::string & extra)
{
    std::ifstream in(path);
    std::string out((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    out += extra;
    const std::string tmp = biospring::tests::temp_path(".grp");
    std::ofstream(tmp) << out;
    return tmp;
}

TEST_F(TestReducer, a_rule_matching_no_atom_anywhere_is_counted_and_costs_no_grain)
{
    // "ZZZ" is in no residue of model.pdb, so this rule can never build.
    const std::string grp = grp_plus(path_reduce_rules, "AOXT ALA ZZZ\n");

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(grp);
    reducer.reduce();

    EXPECT_EQ(reducer.grains_never_built(), 1u);
    // The same model as without the rule: an unmatched rule removes nothing.
    EXPECT_EQ(reducer.target_topology().number_of_particles(), 183);
}

TEST_F(TestReducer, a_rule_matching_somewhere_is_not_counted_however_many_residues_skip_it)
{
    // CB exists in every residue of model.pdb but glycine, so this rule builds
    // for one residue name and is skipped by none of ITS OWN -- what must not
    // happen is for it to be reported as never built.
    const std::string grp = grp_plus(path_reduce_rules, "AXTRA ALA CB\n");

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(grp);
    reducer.reduce();

    EXPECT_EQ(reducer.grains_never_built(), 0u);
}

// A grain name is keyed by residue: one residue type satisfying its own
// terminal rule must not silence another type's. ALA's extra rule matches, ARG's
// cannot, and the count has to be exactly one.
TEST_F(TestReducer, the_count_is_per_residue_not_per_grain_name)
{
    const std::string grp = grp_plus(path_reduce_rules, "XX ALA CB\nXX ARG ZZZ\n");

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(grp);
    reducer.reduce();

    EXPECT_EQ(reducer.grains_never_built(), 1u);
}
